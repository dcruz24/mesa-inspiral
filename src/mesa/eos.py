import os
from typing import Iterable, Optional, Sequence

import numpy as np
import pyMesa as pym


class PyMesaEOS:
    """Checked wrapper for MESA EOS/chem/const/math calls.

    This follows the standalone pyMesaMain interface style, but keeps the old
    attributes used by the inspiral energy-deposition code.
    """

    def __init__(self, cacheDir: str, eos_inlist: str = 'eos.inlist', workDir: Optional[str] = None):
        self.cache_dir = cacheDir
        self.work_dir = os.path.abspath(workDir or os.getcwd())
        self.eos_inlist = self._resolve_path(eos_inlist)

        self.eos_lib, self.eos_def = self._load_module('eos')
        self.const_lib, self.const_def = self._load_module('const')
        self.math_lib, _ = self._load_module('math')
        self.chem_lib, self.chem_def = self._load_module('chem')

        missing = [
            name
            for name, module in (
                ('eos', self.eos_lib),
                ('const', self.const_lib),
                ('math', self.math_lib),
                ('chem', self.chem_lib),
            )
            if module is None
        ]
        if missing:
            raise RuntimeError(f"pyMesa could not load required modules: {', '.join(missing)}")

        self.math_lib.math_init()
        self._init_const()
        self._init_chem()
        self._init_eos()
        self.eos_handle = self._alloc_eos_handle()
        print('Loaded Modules for eos, const, math, chem...')

    def _resolve_path(self, path: str) -> str:
        if os.path.isabs(path):
            return path
        return os.path.join(self.work_dir, path)

    def _load_module(self, name: str):
        try:
            return pym.loadMod(name, self.cache_dir)
        except TypeError:
            return pym.loadMod(name)

    def _check_ierr(self, res, message: str):
        ierr = getattr(res, 'args', {}).get('ierr')
        if ierr is not None and ierr != 0:
            raise pym.MesaError(f'{message} (ierr={ierr})')

    def _init_const(self):
        ierr = 0
        res = self.const_lib.const_init(pym.MESA_DIR, ierr)
        self._check_ierr(res, 'const_init failed')

    def _init_chem(self):
        ierr = 0
        res = self.chem_lib.chem_init('isotopes.data', ierr)
        self._check_ierr(res, 'chem_init failed')

    def _init_eos(self):
        ierr = 0
        res = self.eos_lib.eos_init(pym.EOSDT_CACHE, True, ierr)
        self._check_ierr(res, 'eos_init failed')

    def _alloc_eos_handle(self):
        if not os.path.exists(self.eos_inlist):
            raise FileNotFoundError(f'EOS inlist not found: {self.eos_inlist}')
        ierr = 0
        res = self.eos_lib.alloc_eos_handle_using_inlist(self.eos_inlist, ierr)
        self._check_ierr(res, 'alloc_eos_handle_using_inlist failed')
        self.res = res
        return res.result

    def lookup_chem_ids(self, isotopes: Iterable[str]) -> np.ndarray:
        chem_ids = []
        for isotope in isotopes:
            chem_ids.append(self.chem_lib.chem_get_iso_id(isotope).result)
        return np.asarray(chem_ids)

    def get_basic_eos_result(
        self,
        rho: float,
        temperature: float,
        chem_ids: Sequence[int],
        xa: Sequence[float],
        net_iso: Optional[Sequence[int]] = None,
    ):
        chem_ids_arr = np.asarray(chem_ids)
        xa_arr = np.asarray(xa, dtype=float)
        species = xa_arr.size
        if chem_ids_arr.size != species:
            raise ValueError('chem_ids and xa must contain the same number of species')

        net_iso_arr = np.asarray(net_iso if net_iso is not None else chem_ids_arr)
        nres = self.eos_def.num_eos_basic_results
        res = np.zeros(nres)
        d_dlnRho = np.zeros(nres)
        d_dlnT = np.zeros(nres)
        d_dxa = np.zeros((nres, species))
        ierr = 0

        eos_result = self.eos_lib.eosdt_get(
            self.eos_handle,
            species,
            chem_ids_arr,
            net_iso_arr,
            xa_arr,
            rho,
            np.log10(rho),
            temperature,
            np.log10(temperature),
            res,
            d_dlnRho,
            d_dlnT,
            d_dxa,
            ierr,
        )
        self._check_ierr(eos_result, 'eosdt_get failed')
        return eos_result

    def get_internal_energy(
        self,
        rho: float,
        temperature: float,
        chem_ids: Sequence[int],
        xa: Sequence[float],
        net_iso: Optional[Sequence[int]] = None,
    ) -> float:
        eos_result = self.get_basic_eos_result(rho, temperature, chem_ids, xa, net_iso=net_iso)
        i_lne = self.eos_def.i_lne - 1
        return float(np.exp(eos_result.args['res'][i_lne]))

    def solve_temperature_from_internal_energy(
        self,
        rho: float,
        temperature_guess: float,
        chem_ids: Sequence[int],
        xa: Sequence[float],
        internal_energy: float,
        net_iso: Optional[Sequence[int]] = None,
        logT_tol: float = 1e-8,
        other_tol: float = 1e-8,
        max_iter: int = 20,
    ) -> float:
        chem_ids_arr = np.asarray(chem_ids)
        xa_arr = np.asarray(xa, dtype=float)
        species = xa_arr.size
        if chem_ids_arr.size != species:
            raise ValueError('chem_ids and xa must contain the same number of species')

        net_iso_arr = np.asarray(net_iso if net_iso is not None else chem_ids_arr)
        nres = self.eos_def.num_eos_basic_results
        res = np.zeros(nres)
        d_dlnRho_const_T = np.zeros(nres)
        d_dlnT_const_Rho = np.zeros(nres)
        d_dxa_const_TRho = np.zeros((nres, species))
        eos_calls = 0
        ierr = 0

        arg_not_provided = -9e99
        logRho = np.log10(rho)
        logT_guess = np.log10(temperature_guess)

        eos_result = self.eos_lib.eosdt_get_t(
            self.eos_handle,
            species,
            chem_ids_arr,
            net_iso_arr,
            xa_arr,
            logRho,
            self.eos_def.i_lne,
            np.log(internal_energy),
            logT_tol,
            other_tol,
            max_iter,
            logT_guess,
            arg_not_provided,
            arg_not_provided,
            arg_not_provided,
            arg_not_provided,
            logT_guess,
            res,
            d_dlnRho_const_T,
            d_dlnT_const_Rho,
            d_dxa_const_TRho,
            eos_calls,
            ierr,
        )
        self._check_ierr(eos_result, 'eosdt_get_t failed')
        return float(10 ** eos_result.args['logt_result'])

    def free_eos_handle(self):
        free_handle = getattr(self.eos_lib, 'free_eos_handle', None)
        if free_handle is None:
            return None
        ierr = 0
        res = free_handle(self.eos_handle, ierr)
        self._check_ierr(res, 'free_eos_handle failed')
        return res


class pyMesa_context(PyMesaEOS):
    """Backward-compatible name used by the inspiral driver."""


def getEosResult_IE(context, rho, T, species, chem_ids, net_iso, xa):
    return context.get_basic_eos_result(rho, T, chem_ids, xa, net_iso=net_iso)


def getEosResult_Temperature(context, rho, T, species, chem_ids, net_iso, xa, oldIE, IE, i_lnE):
    return context.solve_temperature_from_internal_energy(rho, T, chem_ids, xa, IE, net_iso=net_iso)
