#include "alpinist/ALPPhotonMixing.h"


namespace alpinist {


ALPPhotonMixing::ALPPhotonMixing(double mass, double coupling, ref_ptr<MagneticField> field, ref_ptr<PlasmaDensity> density, double limit) : Module() {
	setAxionMass(mass);
	setCouplingConstant(coupling);
	setMagneticField(field);
	setPlasmaDensity(density);
	setInteractionTag("ALPPh");
	setLimit(limit);
	setDescription("ALPPhotonMixing::ALPPhotonMixing");	
}

void ALPPhotonMixing::setAxionMass(double mass) {
	axionMass = mass;
}

void ALPPhotonMixing::setCouplingConstant(double g) {
	couplingConstant = g;
}

void ALPPhotonMixing::setMagneticField(ref_ptr<MagneticField> bField) {
	magneticField = bField;
}

void ALPPhotonMixing::setPlasmaDensity(ref_ptr<PlasmaDensity> density) {
	plasmaDensity = density;
}

void ALPPhotonMixing::setInteractionTag(std::string tag) {
	interactionTag = tag;
}

void ALPPhotonMixing::setLimit(double l) {
	limit = l;
}

double ALPPhotonMixing::getAxionMass() const {
	return axionMass;
}

double ALPPhotonMixing::getCouplingConstant() const {
	return couplingConstant;
}

ref_ptr<MagneticField> ALPPhotonMixing::getMagneticField() const {
	return magneticField;
}

ref_ptr<PlasmaDensity> ALPPhotonMixing::getPlasmaDensity() const {
	return plasmaDensity;
}

std::string ALPPhotonMixing::getInteractionTag() const {
	return interactionTag;
}

double ALPPhotonMixing::getLimit() const {
	return limit;
}

WaveFunction3c ALPPhotonMixing::getWaveFunction(const Candidate& candidate) const {
	Eigen::RowVector3cd fields;
	std::complex<double> em2 = (candidate.getProperty(varFieldEM2)).toComplexDouble();
	std::complex<double> em1 = (candidate.getProperty(varFieldEM1)).toComplexDouble();
	std::complex<double> alp = (candidate.getProperty(varFieldALP)).toComplexDouble();
	fields << em2, em1, alp;
	WaveFunction3c psi(fields);
	psi.normalise();

	return psi;
}

void ALPPhotonMixing::process(Candidate* candidate) const {
	if (axionMass == 0. || couplingConstant == 0.)
		return;

	int id = candidate->current.getId();
	if (id != 22 && id != 51)
		return;

	// retrieve candidate properties
	double redshift = candidate->getRedshift();
	double energy = candidate->current.getEnergy() * (1 + redshift);
	Vector3d position2 = candidate->current.getPosition();
	Vector3d position1 = candidate->previous.getPosition();
	// Vector3d position0 = candidate->created.getPosition();
	Vector3d position = (position1 + position2) / 2.;
	double step = candidate->getCurrentStep();
	double w0 = candidate->getWeight();

	// get fields at the beginning of the step
	WaveFunction3c initialState = getWaveFunction(*candidate);

	// local magnetic field and plasma density
	// note that propagation is along the x-direction, so transverse field is along y and z
	double density = plasmaDensity->getDensity(position, redshift);
	Vector3d bField = (magneticField->getField(position, redshift)).getPerpendicularTo(candidate->current.getDirection());

	// compute mixing matrix
	MixingParameters mixing(couplingConstant, axionMass, energy, density, bField);
	double oscillationLength = 2. / mixing.deltaOscillation;


	// allow oscillation and adjust step lengths  accordingly
	// this is important when the module is called with other interaction modules (e.g., pair production)
	if (limit > 0.) {
		if (step > oscillationLength * limit) {
			evolve(candidate, mixing, initialState, step);	
			candidate->limitNextStep(oscillationLength * limit);
		} else {
			evolve(candidate, mixing, initialState, step);	
		}
	} else {
		evolve(candidate, mixing, initialState, step);
	}
}

void ALPPhotonMixing::evolve(Candidate* candidate, MixingParameters& mixing, const WaveFunction3c& initialState, const double& distance) const {
	int id = candidate->current.getId();

	// solve equations of motion for the three states
	Eigen::Matrix3d mixingMatrix = mixing.getMixingMatrix();
	mixing.solveEigenSystem(mixingMatrix);
	Eigen::Matrix3cd evolutionMatrix = mixing.getSimilarityTransformingMatrix(distance);

	// compute (spatial) evolution of the system
	Eigen::Vector3cd finalStateVector =  evolutionMatrix * initialState.getCoefficients().transpose();
	WaveFunction3c finalState(finalStateVector);
	finalState.normalise();

	// estimate oscillation probability
	double probPA = std::norm(std::conj(initialState[0]) * finalState[2] + std::conj(initialState[1]) * finalState[2]); 
	double probAP = std::norm(std::conj(initialState[2]) * finalState[0] + std::conj(initialState[2]) * finalState[1]);
	double probOsc = (id == 22) ? probPA : probAP;

	// decide whether oscillation occurs (only particle id changes)
	Random& random = Random::instance();
	if (random.rand() < probOsc) {
		candidate->current.setId((id == 22) ? 51 : 22);
	} 

	// return to original basis
	finalState.operate(evolutionMatrix);

	// assign weights to ALPs/photons (separately from the usual weights)
	Eigen::Matrix3cd densityMatrix = finalState.getDensityMatrix();
	candidate->setProperty(varProbabilityPhoton, Variant::fromDouble(std::abs(densityMatrix(0, 0)) + std::abs(densityMatrix(1, 1))));
	candidate->setProperty(varProbabilityALP, Variant::fromDouble(std::abs(densityMatrix(2, 2))));

	// store the new fields at the end of the step
	candidate->setProperty(varFieldEM1, Variant::fromComplexDouble(finalState[1]));
	candidate->setProperty(varFieldEM2, Variant::fromComplexDouble(finalState[0]));
	candidate->setProperty(varFieldALP, Variant::fromComplexDouble(finalState[2]));	
}


/*********************************************************************************************************************/

MixingParameters::MixingParameters() {
}

MixingParameters::MixingParameters(const double& couplingConstant, const double& axionMass, const double& energy, const double& densityMedium, const Vector3d& magneticField) {
	computeParameters(couplingConstant, axionMass, energy, densityMedium, magneticField);
}

void MixingParameters::computeParameters(const double& couplingConstant, const double& axionMass, const double& energy, const double& densityMedium, const Vector3d& bField) {
	double magneticField = bField.getR();
	double 	deltaQED = fineStructure * energy / (90. * M_PI) * pow_integer<2>(magneticField / B_schwinger) / (h_dirac * c_light);
	double plasmaFrequency = computePlasmaFrequency(densityMedium);
	deltaPlasma = - pow_integer<2>(plasmaFrequency) / (2. * energy) * (h_dirac / c_light);
	deltaAxion = - pow_integer<2>(axionMass) / (2. * energy) / (h_dirac * c_light);
	deltaParallel = deltaPlasma + 4. * deltaQED;	
	deltaPerpendicular = deltaPlasma + 7. * deltaQED;
	deltaMixing = 0.5 * couplingConstant * (magneticField / pow(h_dirac * c_light, - 1.5) / sqrt(mu0)) / (h_dirac * c_light);
	deltaOscillation = sqrt(pow_integer<2>(deltaAxion - deltaPerpendicular) + pow_integer<2>(2 * deltaMixing));
	mixingAngle = 0.5 * asin(2. * deltaMixing / deltaOscillation);
	angleB = acos(bField.getZ() / bField.getR()); 
}

void MixingParameters::printParameters() {
	std::cout << "plasma: " << deltaPlasma << std::endl;
	std::cout << "axion: " << deltaAxion << std::endl;
	// std::cout << "QED: " << deltaQED << std::endl;
	std::cout << "mixing: " << deltaMixing << std::endl;
	std::cout << "parallel: " << deltaParallel << std::endl;
	std::cout << "perpendicular: " << deltaPerpendicular << std::endl;
	std::cout << "oscillation: " << deltaOscillation << std::endl;
	std::cout << "mixing angle: " << mixingAngle << std::endl;
	// std::cout << "mixing angle: " <<  << std::endl;
}

Eigen::Matrix3d MixingParameters::getMixingMatrix() const {
	Eigen::Matrix3d M;
	double sQ = sin(angleB);
	double cQ = cos(angleB);
	double deltaXX = deltaParallel * sQ * sQ + deltaPerpendicular * cQ * cQ;
	double deltaYY = deltaParallel * cQ * cQ + deltaPerpendicular * sQ * sQ;
	double deltaXY = (deltaParallel - deltaPerpendicular) * sQ * cQ;
	M << deltaXX, deltaXY, deltaMixing * sQ, deltaXY, deltaYY, deltaMixing * cQ, deltaMixing * sQ, deltaMixing * cQ, deltaAxion;
	M << deltaYY, deltaXY, deltaMixing * sQ, deltaXY, deltaXX, deltaMixing * cQ, deltaMixing * sQ, deltaMixing * cQ, deltaAxion;

	return M;
}

void MixingParameters::solveEigenSystem(const Eigen::Matrix3d& M) {
	eigenSolver.compute(M);
}

Eigen::Vector3d MixingParameters::getEigenValues() const {
	return eigenSolver.eigenvalues();
}

Eigen::Matrix3cd MixingParameters::getEigenVectors() const {
	return eigenSolver.eigenvectors().cast<std::complex<double>>();
} 

Eigen::Matrix3cd MixingParameters::getExponentialDiagonalEigenvalueMatrix(const double& x) const {
	Eigen::Vector3d eigenValues = getEigenValues();
	Eigen::Vector3cd v;
	for (size_t i = 0; i < 3; i++) {
		std::complex<double> l = {0., - 0.5 * x * eigenValues(i) - energy / (h_dirac * c_light) * 0.5 * x};
		v(i) = exp(l);
	} 

	Eigen::Matrix3cd D = v.asDiagonal();
	return D;
}

Eigen::Matrix3cd MixingParameters::getSimilarityTransformingMatrix(const double& x) const {
	return getEigenVectors().inverse() * getExponentialDiagonalEigenvalueMatrix(x) * getEigenVectors();
}



} // namespace