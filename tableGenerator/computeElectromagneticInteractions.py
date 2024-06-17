import os
import sys
import numpy as np

from warnings import filterwarnings

# import from CRPropa3-data folder
sys.path.append('/Users/rab/Dropbox/softwares/CRPropa/CRPropa3-data-current')
import gitHelp as gh
import photonField

from interactionsElectromagnetic import *
from meanFreePath import computeInteractionRate
from constantsUnits import *

filterwarnings('ignore')


###############################################################################
###############################################################################
def process(interaction, field, folder = '../data', logTabSMin = 4., logTabSMax = 23.):
	""" 
	Calculate the interaction rates for a given process on a given photon field.

	# Input
	. interaction: EM process
	. field   : photon field as defined in photonField.py
	. name    : name of the process which will be calculated. Necessary for the naming of the data folder
	"""
	name = interaction.label 
	nIntegral = 21
	epsilon = 0.
	if interaction.label == 'InverseComptonScattering':
		nIntegral = 22
		if 'CMB' in field.name:
			epsilon = 1e-4
		elif 'IRB' in field.name:
			epsilon = 1e-4
			logTabSMin = 8. # not accurate but prevents divergences
	

	if not folder.endswith('/'):
		folder += '/'
	folder = folder + name
	if not os.path.exists(folder):
		os.makedirs(folder)

	# tabulated energies, limit to energies where the interaction is possible
	E = np.logspace(9, 23, 281) * eV
	Erange = (E[0], E[-1])

	Emin = interaction.minimumEnergyLab(field, Erange)
	E = E[E > Emin]
	
	# -------------------------------------------
	# calculate interaction rates
	# -------------------------------------------
	# tabulated values of s_kin = s - mc^2
	# Note: integration method (Romberg) requires 2^n + 1 log-spaced tabulation points
	sKin = np.logspace(logTabSMin, logTabSMax, 2 ** nIntegral + 1) * eV ** 2
	xs = interaction.computeCrossSections(sKin, epsilon = epsilon)
	rate = computeInteractionRate(sKin, xs, E, field)

	# save
	if interaction.helicities is None:
		fname = folder + '/rate_%s.txt' % (field.name)
	else:
		h = interaction.polarisationIdentifier()
		fname = folder + '/rate_pol%i_%s.txt' % (h, field.name)

	data = np.c_[np.log10(E / eV), rate]
	fmt = '%.2f\t%8.7e'
	try:
		git_hash = gh.get_git_revision_hash()
		headerStr = '%s interaction rates\nphoton field: %s\n' % (name, field.info)
		headerStr += 'Produced with crpropa-data version: ' + git_hash + '\n'
		headerStr += 'log10(E/eV), 1/lambda [1/Mpc]'
		header = (headerStr)
	except:
		headerStr = '%s interaction rates\n' % (name)
		headerStr += 'photon field: %s\n' % (field.info)
		headerStr += 'log10(E/eV), 1/lambda [1/Mpc]'
		header = (headerStr)
	np.savetxt(fname, data, fmt = fmt, header = header)


	# -------------------------------------------
	# calculate cumulative differential interaction rates for sampling s values
	# -------------------------------------------

	sKinMin = max(interaction.thresholdEnergy2(), 4. * field.getEmin() * E[0])

	# tabulated values of s_kin = s - mc^2, limit to relevant range
	# Note: use higher resolution and then downsample
	sKin = np.logspace(logTabSMin, logTabSMax, 380000 + 1) * eV ** 2
	sKin = sKin[sKin > sKinMin]

	xs = interaction.computeCrossSections(sKin, epsilon = epsilon)
	rate = computeInteractionRate(sKin, xs, E, field, cdf = True)

	# downsample
	sKin_save = np.logspace(logTabSMin, logTabSMax, 190 + 1) * eV ** 2
	sKin_save = sKin_save[sKin_save > sKinMin]
	rate_save = np.array([np.interp(sKin_save, sKin, r) for r in rate])
	
	# save
	data = np.c_[np.log10(E / eV), rate_save]  # prepend log10(E/eV) as first column
	row0 = np.r_[0, np.log10(sKin_save / eV ** 2)][np.newaxis]
	data = np.r_[row0, data]  # prepend log10(s_kin/eV^2) as first row

	# save
	if interaction.helicities is None:
		fname = folder + '/cdf_%s.txt' % (field.name)
	else:
		h = interaction.polarisationIdentifier()
		fname = folder + '/cdf_pol%i_%s.txt' % (h, field.name)

	fmt = '%.2f' + '\t%6.5e' * np.shape(rate_save)[1]
	try:
		git_hash = gh.get_git_revision_hash()
		headerStr = '%s cumulative differential rate\n' % name
		headerStr += 'photon field: %s\n' % field.info
		headerStr += 'Produced with crpropa-data version: ' + git_hash + '\n'
		headerStr += 'log10(E/eV), d(1/lambda)/ds_kin [1/Mpc/eV^2] for log10(s_kin/eV^2) as given in first row'
		header = (headerStr)
	except:
		headerStr = '%s cumulative differential rate' % name
		headerStr += 'photon field: %s\n' % field.info
		headerStr +='log10(E/eV), d(1/lambda)/ds_kin [1/Mpc/eV^2] for log10(s_kin/eV^2) as given in first row' 
		header = (header)
	np.savetxt(fname, data, fmt = fmt, header = header)

	del data, rate, sKin, sKin_save, rate_save




###########################################################################################
###########################################################################################
if __name__ == "__main__":

	fields = [
		photonField.CMB(),
		# photonField.EBL_Saldana21('best')
		# photonField.EBL_Saldana21('upper'),
		# photonField.EBL_Saldana21('lower')
	]
	

	interactionsPP = []
	interactionsICS = []

	interactionsPP.append(PairProduction())
	for h in  (1, 2, 3, 4, 5):
		interactionsPP.append(PairProduction(identifier = h))

	interactionsICS.append(InverseComptonScattering())
	for h in  (1, 2, 3, 4, 5, 6):
		interactionsICS.append(InverseComptonScattering(identifier = h))



	for field in fields:
		print('##########################################')
		print('==> photon field = ', field.name)

		# print('  - pair production ')
		# for interaction in interactionsPP:
		# 	print('    . ', interaction.helicities)
		# 	process(interaction, field)

		print('  - inverse Compton scattering ')
		for interaction in interactionsICS:
			print('    . ', interaction.helicities)
			process(interaction, field)
