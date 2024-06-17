import numpy as np
from abc import ABC, abstractmethod
from warnings import filterwarnings

from constantsUnits import *


filterwarnings('ignore')


###############################################################################
###############################################################################
class ElectromagneticInteraction(ABC):
	"""
	Abstract base class for electromagnetic processes of the type:
	  X + gamma -> ...,
	wherein X is a photon or charged lepton.
	"""
	@property
	def name(self):
		return self._name
	
	@name.setter
	def name(self, n):
		self._name = n
	
	@property
	def label(self):
		return self._label
	
	@label.setter
	def label(self, l):
		self._label = l

	@property
	def sMin(self):
		return self._sMin
	
	@sMin.setter
	def sMin(self, s):
		self._sMin = s

	@property
	def incidentParticle(self):
		return self._incidentParticle
	
	@incidentParticle.setter
	def incidentParticle(self, p):
		self._incidentParticle = p

	@property
	def massInitial(self):
		return self._massInitial
	
	@massInitial.setter
	def massInitial(self, m):
		self._massInitial = m

	@property
	def helicities(self):
		return self._helicities
	
	@helicities.setter
	def helicities(self, h):
		if h not in (None, '++++', '+++-', '++-+', '++--', '+-++', '+-+-', '+--+', '+---', '-+++', '-++-', '-+-+', '-+--', '--++', '--+-', '---+', '----'):
			raise ValueError('Helicities must be either None (unpolarised) or a combination of \'+\' and \'-\', containing 4 signs.')
		self._helicities = h

	def thresholdEnergy2(self):
		return self._sMin
	
	@abstractmethod
	def polarisationIdentifier(self):
		pass

	@abstractmethod
	def crossSection(self, s):
		pass

	@abstractmethod
	def computeCrossSections(self, sKin):
		""" 
		Get cross section for tabulated s_kin 
		"""
		pass

	def minimumEnergyLab(self, field, Erange):
		""" 
		Return minimum required for interaction *sigma* with *field* 
		"""
		return self.thresholdEnergy2() / 4. / field.getEmax()

	def getHelicitiesLabel(self):
		"""
		Return a string corresponding to the helicity combination.
		"""
		if self._helicities is None:
			return ''
		else:
			return self._helicities.replace('+', 'p').replace('-', 'm')


###############################################################################
###############################################################################
class PairProduction(ElectromagneticInteraction):
	"""
	Breit-Wheeler pair production:
	  gamma + gamma -> e+ + e-
	"""
	def __init__(self, helicities = None, identifier = None):
		self.name = 'pair production'
		self.label = 'PairProduction'
		self.sMin = 4. * me2
		self.incidentParticle = 22
		self.massInitial = 0.
		self.helicities = helicities

		# initialise based on case; one helicity combination for each individual cross section
		if identifier is not None:
			if identifier == 1:
				self.helicities = '++++'
			elif identifier == 2:
				self.helicities = '++--'
			elif identifier == 3:
				self.helicities = '-+--'
			elif identifier == 4:
				self.helicities = '-+-+'
			elif identifier == 5:
				self.helicities = '---+'
			else:
				raise(ValueError('Cannot initialise pair production with helicity combination identifier %i.' % identifier))


	def polarisationIdentifier(self):
		"""
		"""
		if self.helicities is None or self.helicities in ('0', '0000'):
			return 0
		elif self.helicities in ('++++', '----'):
			return 1
		elif self.helicities in ('++--', '--++'):
			return 2
		elif self.helicities in ('-+--', '+---', '-+++', '+-++'):
			return 3
		elif self.helicities in ('-+-+', '-++-', '+--+', '+-+-'):
			return 4
		elif self.helicities in ('---+', '--+-', '++-+', '+++-'):
			return 5

	def crossSection(self, s, epsilon = 0.):
		""" 
		Pair production cross section (Breit-Wheeler), see Lee 1996 
		"""
		if (s < self.sMin):
			return 0.
		
		b = np.sqrt(1. - self.sMin / s)
		k = 0.5 * (np.log1p(b) - np.log1p(-b))

		if self.helicities is None:
			return sigmaThomson * 3 / 16 * (1 - b * b) * ((3 - b ** 4) * (np.log1p(b) - np.log1p(-b)) - 2. * b * (2 - b * b))
		
		elif self.helicities in ('++++', '----'):
			return sigmaThomson * 3 / 2. * me2 * (s - 2 * me2 + s * b) * (s + 4. * me2 / b * k) / (s ** 3)
		
		elif self.helicities in ('++--', '--++'):
			return sigmaThomson * 3. / 2. * me2 * b * (s - 2 * me2 - s * b) * (s + 4 * me2 / b * k) / (s ** 3)
		
		elif self.helicities in ('-+--', '+---', '-+++', '+-++'):
			return sigmaThomson * 3. * me4 * (2 * b * s * (s + 2 * me2) + 16 * me2 * (me2 - s) * k) / (s ** 4 * b ** 2)
		
		elif self.helicities in ('-+-+', '-++-', '+--+', '+-+-'):
			return sigmaThomson * 3. * me2 * ((me2 - s) * s * b + (4 * me4 - 2 * me2 * s + s ** 2) * k) / (s ** 3 * b ** 2)
		
		elif self.helicities in ('---+', '--+-', '++-+', '+++-'):
			return 0.
		
	def computeCrossSections(self, sKin, epsilon = 0.):
		""" 
		Get cross section for tabulated s_kin 
		"""
		return np.array([self.crossSection(s, epsilon = epsilon) for s in sKin])

	def computeBeta(self, s):
		"""
		"""
		return np.sqrt(1. - s / self.sMin)



###############################################################################
###############################################################################
class InverseComptonScattering(ElectromagneticInteraction):
	"""
	Inverse Compton Scattering:
	  e + gamma -> e + gamma
	"""
	def __init__(self, helicities = None, identifier = None):
		self.name = 'inverse Compton scattering'
		self.label = 'InverseComptonScattering'
		self.sMin = me2
		self.incidentParticle = 11
		self.massInitial = mass_electron
		self.helicities = helicities

		# initialise based on case; one helicity combination for each individual cross section
		if identifier is not None:
			if identifier == 1:
				self.helicities = '++++'
			elif identifier == 2:
				self.helicities = '++--'
			elif identifier == 3:
				self.helicities = '-+--'
			elif identifier == 4:
				self.helicities = '-+-+'
			elif identifier == 5:
				self.helicities = '--+-'
			elif identifier == 6:
				self.helicities = '-++-'
			else:
				raise(ValueError('Cannot initialise inverse Compton scattering with helicity combination identifier %i.' % identifier))


	def polarisationIdentifier(self):
		"""
		"""
		if self.helicities is None or self.helicities in ('0', '0000'):
			return 0
		elif self.helicities in ('++++', '----'):
			return 1		
		elif self.helicities in ('++--', '--++'):
			return 2
		elif self.helicities in ('-+--', '---+', '+++-', '+-++'):
			return 3
		elif self.helicities in ('-+-+', '+-+-'):
			return 4
		elif self.helicities in ('--+-', '-+++', '+---', '++-+'):
			return 	5
		elif self.helicities in ('-++-', '+--+'):
			return 6

	def crossSection(self, s, epsilon = 0.):
		"""
		Inverse Compton scattering cross sections, see Lee 1996
		"""
		sMin = self.thresholdEnergy2()
		if s < sMin:  # numerically unstable close to sMin
			return 0.

		# note: formula unstable for (s - smin) / smin < 1e-5
		a = 3. / 8. * sigmaThomson * (me2 / (me2 - s) ** 4)
		k = np.log(s / me2)

		# the cross section flattens at low energies
		# to prevent numerical problems, we return the ratio
		sigmas = {'1': 0.9673117210572482, '2': 0.021846706047965028, '3': 0.01704508535119666, '4': 0.37837252326629683, '5': 0.5097246271735112, '6': 0.04423928290101834}

		if self.helicities is None:
			b = (s - sMin) / (s + sMin)
			f1 = 2. / (b * (1 + b)) * (2 + 2 * b - b ** 2 - 2 * b ** 3)
			f2 = 1. / b ** 2 * (2. - 3. * b ** 2 - b ** 3) * np.log((1. + b) / (1. - b))
			return 3. * sigmaThomson / 8. * me2 / s / b * (f1 - f2)
		
		elif self.helicities in ('++++', '----'):
			if s < (1. + epsilon) * sMin:
				return sigmas['1'] * sigmaThomson
			else:
				f1 = 2. * (me2 - s) * (me6 + 15. * me4 * s - 12. * me2 * s ** 2 + 2. * s ** 3)
				f2 = 4. * s * (4. * me6 + 2. * me4 * s - 4. * me2 * s ** 2 + s ** 3) * k
				return a / s * (f1 + f2)
		
		elif self.helicities in ('++--', '--++'):
			if s < (1. + epsilon) * sMin:
				return sigmas['2'] * sigmaThomson
			else:
				return 2. * a * me4 / s / s * (me6 - 6. * me4 * s + 3. * me2 * s * s + 2. * s ** 3 - 6. * me2 * s ** 2 * k) 
		
		elif self.helicities in ('-+--', '---+', '+++-', '+-++'):
			if s < (1. + epsilon) * sMin:
				return sigmas['3'] * sigmaThomson
			else:
				return 2. * a * me4 / s * (me4 + 4. * me2 * s - 5. * s ** 2 + 2. * s * (2. * me2 + s) * k)

		elif self.helicities in ('-+-+', '+-+-'):
			if s < (1. + epsilon) * sMin:
				return sigmas['4'] * sigmaThomson
			else:
				return 2. * a * (2. * me6 + 3. * me4 * s - 6. * me2 * s * s + s ** 3 + 6. * me4 * s * k)

		elif self.helicities in ('--+-', '-+++', '+---', '++-+'):
			if s < (1. + epsilon) * sMin:
				return sigmas['5'] * sigmaThomson
			else:
				return - a * (2. * me2 * (me2 - s) * (5. * me2 + s) + 4 * (me6 + 2. * me4 * s) * k)

		elif self.helicities in ('-++-', '+--+'):
			if s < (1. + epsilon) * sMin:
				return sigmas['6'] * sigmaThomson
			else:
				return 2. * a * (me6 - 6. * me4 * s + 3. * me2 * s ** 2 + 2. * s * s * s - 6. * me2 * s * s * k)
		

	def computeCrossSections(self, sKin, epsilon = 0.):
		""" 
		Get cross section for tabulated s_kin 
		"""
		return np.array([self.crossSection(sK + me2, epsilon = epsilon) for sK in (sKin)])
	
	def computeBeta(self, s):
		"""
		"""
		return (s - self.sMin) / (s - self.sMin)

###############################################################################
###############################################################################


