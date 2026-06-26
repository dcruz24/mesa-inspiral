# Fit Check Plots

These scripts regenerate visual checks for drag-fit and neutron-star EOS helper functions.
They are kept under `tests/` because their outputs are verification artifacts rather than examples for the main workflow.

Run the static checks directly with:

```bash
MPLBACKEND=Agg python tests/fit_checks/macleod_fit_curves.py
MPLBACKEND=Agg python tests/fit_checks/macleod_fit_surfaces.py
MPLBACKEND=Agg python tests/fit_checks/neutron_star_eos_curves.py
```

Generated figures are written to `tests/fit_checks/plots/`.
The interactive Plotly script is optional and exits cleanly if Plotly is not installed.
