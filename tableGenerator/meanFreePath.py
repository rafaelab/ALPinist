import sys
import numpy as np

# SciPy's integration routines
from scipy.integrate import cumulative_trapezoid, romb

# import from CRPropa3-data folder
sys.path.append('CRPropa3-data/')
import gitHelp as gh
from interactionRate import mean_log_spacing, calculateDensityIntegral

# import units and constants from CRPropa
from crpropa import eV, Mpc
from constantsUnits import *



def computeInteractionRate(sKin, xs, E, field, z = 0, cdf = False):
	"""
	Calculate the interaction rate for given tabulated cross sections against an isotropic photon background.
	The tabulated cross sections need to be of length n = 2^i + 1 and the tabulation points log-linearly spaced.
	This is equivalent to the CRPropa3-data function found in interactionRate.py called `calc_rate_s` (with different arguments).

	# Input
	. sKin    : tabulated (s - m**2) for cross sections [J^2]
	. xs      : tabulated cross sections [m^2]
	. E       : (array of) cosmic ray energies [J]
	. field   : photon background, see photonField.py
	. z       : redshift
	. cdf     : calculate cumulative differential rate

	# Output
	. interaction rate 1/lambda(gamma) [1/Mpc] or
	. cumulative differential rate d(1/lambda)/d(sKin) [1/Mpc/J^2]
	"""
	# size of table
	nE, nS = len(E), len(sKin)

	if cdf:
		# precalculate the field integral if it not exists and load it afterwards
		calculateDensityIntegral(field)
		file = 'temp/fieldDensity/' + field.name + '.txt'
		densityIntegral = np.loadtxt(file)

		# interpolate
		I = np.zeros((nE, nS))

		for i in range(nE):
			I[i, :] = np.interp(sKin / 4 / E[i], densityIntegral[:, 0], densityIntegral[:, 1])
		y = np.array([xs * sKin for i in range(nE)]) * I

		return cumulative_trapezoid(y = y, x = sKin, initial = 0) / 8. / np.expand_dims(E, -1) ** 2 * Mpc    
		
	else:
		F = cumulative_trapezoid(x = sKin, y = sKin * xs, initial = 0)
		n = field.getDensity(np.outer(1. / (4. * E), sKin), z)
		y = n * F / sKin

		return romb(y, dx = mean_log_spacing(sKin)) / 2. / E * Mpc

