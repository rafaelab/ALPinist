#include "alpinist/PolarisationCorrections.h"

namespace alpinist {


PolarisationCorrections::PolarisationCorrections() {
}

PolarisationCorrections::PolarisationCorrections(int pId, std::string interaction) {
	setParticleId(pId);
	setInteraction(interaction);
}

void PolarisationCorrections::setParticleId(int pId) {
	particleId = pId;
}

void PolarisationCorrections::setInteraction(std::string interactionProcess) {
	interaction = interactionProcess;
}

int PolarisationCorrections::getParticleId() const {
	return particleId;
}

std::string PolarisationCorrections::getInteraction() const {
	return interaction;
}

void PolarisationCorrections::process(Candidate* candidate) const {
	if (abs(candidate->current.getId()) != abs(particleId)) 
		return;

	if (candidate->getTagOrigin() != interaction)
		return;

	// probabilities are weights and are thus preserved
	double pALP = (candidate->parent)->getProperty(varProbabilityALP);
	double pPhoton = (candidate->parent)->getProperty(varProbabilityPhoton);
	candidate->setProperty(varProbabilityALP, Variant::fromDouble(pALP));
	candidate->setProperty(varProbabilityPhoton, Variant::fromDouble(pPhoton));


	// the initial particle interacts with another one with random polarisation, so final state is random
	// even if pairs are produced, statistically, the correlation between each one can be neglected.
	Random& random = Random::instance();
	double theta = random.randUniform(0, 2 * M_PI);
	double phi1 = random.randUniform(0, 2 * M_PI);
	double phi2 = random.randUniform(0, 2 * M_PI);
	std::complex<double> em1 = {cos(theta) * cos(phi1), cos(theta) * sin(- phi1)};
	std::complex<double> em2 = {sin(theta) * cos(phi2), sin(theta) * sin(- phi2)};
	std::complex<double> alp = {0., 0.};
	candidate->setProperty(varFieldEM1, Variant::fromComplexDouble(em1));
	candidate->setProperty(varFieldEM2, Variant::fromComplexDouble(em2));
	candidate->setProperty(varFieldALP, Variant::fromComplexDouble(alp));
}

} // namespace alpinist