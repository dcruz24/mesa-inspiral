# ----------------------------------------------
# Physical Constants in CGS Units
# ----------------------------------------------
G_cgs        = 6.6743e-8          # Gravitational constant [cm^3 g^-1 s^-2]
c_cgs        = 2.99792458e10      # Speed of light [cm/s]
Msun_cgs     = 1.9891e33          # Solar mass [g]
Rsun_cgs     = 6.957e10           # Solar radius [cm]
pc_cgs       = 3.08567758128e18   # Parsec [cm]
kpc_cgs      = 1.0e3 * pc_cgs     # Kiloparsec [cm]
Mpc_cgs      = 1.0e6 * pc_cgs     # Megaparsec [cm]
au_cgs       = 1.49598073e13      # Astronomical Unit [cm]
yr_cgs       = 3.1556952e7        # Julian year [s]
kyr_cgs      = 1.0e3 * yr_cgs     # Kiloyear [s]
Myr_cgs      = 1.0e6 * yr_cgs     # Megayear [s]
Gyr_cgs      = 1.0e9 * yr_cgs     # Gigayear [s]
m_p_cgs      = 1.67262192e-24     # Proton mass [g]
sigmaT_cgs   = 6.652e-25          # Thomson cross section [cm^2]
sigma_SB_cgs = 5.670374419e-5     # Stefan–Boltzmann constant [erg cm^-2 s^-1 K^-4]

# ----------------------------------------------
# Conversion Factors to Geometrized Units (G = c = Msun = 1)
# Source: Hector O. Silva
# ----------------------------------------------

# -------------------------
# Fundamental constants
# -------------------------
c_geo = 1.0     # Speed of light (normalized)
G_geo = 1.0     # Gravitational constant (normalized)

# -------------------------
# Length Units
# -------------------------

cm_geo     = 6.7706e-6         # 1 cm in code units
fm_geo     = 1.0e-13 * cm_geo
cm3_geo    = cm_geo**3
cmm3_geo   = 1.0 / cm3_geo
R_km_geo   = 0.67706           # 1 km in code units

# -------------------------
# Mass Units
# -------------------------
g_geo      = 5.0279e-34        # 1 g in code units
mb_geo     = 1.66e-24 * g_geo
Msun_geo   = 1.0               # Msun = 1 by definition

# -------------------------
# Time Units
# -------------------------

s_geo      = 2.0296e5          # 1 s in code units

# -------------------------
# Derived Units
# -------------------------
gcm3_geo   = 1.6199e-18        # Density conversion
dyncm2_geo = 1.806e-39         # Pressure conversion
Inertia_geo= g_geo * cm_geo**2
int_e_geo  = 1.1126e-21        # Energy unit (1 erg)


# ------------------------
commonIsotopes = ['h1', 'h2', 'he3','he4', 'li7', 'be7', 'b8', 'c12', 'c13', 'n13', 'n14', 
                     'n15', 'o14', 'o15','o16', 'o17', 'o18', 'f17', 'f18', 'f19', 'ne18', 'ne19', 
                     'ne20', 'ne22', 'mg22','mg24', 'fe56', 'si28']