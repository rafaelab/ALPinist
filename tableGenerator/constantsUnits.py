"""
Define constants and units.
Reexport the ones already defined in CRPropa.
"""

from crpropa import mass_electron
from crpropa import sigma_thomson
from crpropa import alpha_finestructure
from crpropa import c_light, c_squared
from crpropa import eV, keV, MeV, GeV, TeV, PeV, EeV
from crpropa import kpc, Mpc, Gpc


EPl = 1.9561e9 # J
# me2 = 6.72236e-27 # J^2
me2 = (mass_electron * c_squared) ** 2  # squared electron mass [J^2/c^4]
me4 = (mass_electron * c_squared) ** 4  # squared electron mass [J^2/c^4]
me6 = (mass_electron * c_squared) ** 6  # squared electron mass [J^2/c^4]
sigmaThomson = sigma_thomson  # Thomson cross section [m^2]
alpha = alpha_finestructure  # fine structure constant

# define a dictionary with particle masses
particleMassesDictionary = {}
particleMassesDictionary[ 11] = mass_electron
particleMassesDictionary[-11] = mass_electron
particleMassesDictionary[ 22] = 0.