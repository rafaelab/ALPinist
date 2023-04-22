#include "alpinist/SourceFeatures.h"


namespace alpinist {


/*********************************************************************************************************************/

SourceNoPolarisation::SourceNoPolarisation() : SourceFeature() { 
}

void SourceNoPolarisation::prepareCandidate(Candidate& candidate) const {
	Random &random = Random::instance();

	double theta = random.randUniform(0, 2 * M_PI);
	double phi1 = random.randUniform(0, 2 * M_PI);
	double phi2 = random.randUniform(0, 2 * M_PI);

	setPolarisedCandidate(candidate, theta, phi1, phi2);
}


/*********************************************************************************************************************/

SourceLinearPolarisation::SourceLinearPolarisation(double angle, bool coherent, double phase) : SourceFeature() { 
	setCoherent(coherent);
	setPolarisation(angle);
	setPhase(phase);
}

void SourceLinearPolarisation::prepareCandidate(Candidate& candidate) const {
	Random &random = Random::instance();

	double theta = polarisationAngle;
	double phi1 = phaseAngle;
	if (not coherent) {
		phi1 = random.randUniform(0, 2 * M_PI);
	}
	double phi2 = phi1;

	setPolarisedCandidate(candidate, theta, phi1, phi2);
}

void SourceLinearPolarisation::setCoherent(bool b) {
	coherent = b;
}

void SourceLinearPolarisation::setPhase(double phase) {
	phaseAngle = phase;
}

void SourceLinearPolarisation::setPolarisation(double pol) {
	polarisationAngle = pol;
}



/*********************************************************************************************************************/

SourceCircularPolarisation::SourceCircularPolarisation(double angle, bool coherent, double phase) : SourceFeature() { 
	setPhase(phase);
	setCoherent(coherent);
	setPolarisation(angle);
}

void SourceCircularPolarisation::prepareCandidate(Candidate& candidate) const {
	Random &random = Random::instance();

	double theta = polarisationAngle;
	double phi1 = phaseAngle;
	if (not coherent) {
		phi1 = random.randUniform(0, 2 * M_PI);
	}
	double phi2 = phi1 + M_PI / 2.;

	setPolarisedCandidate(candidate, theta, phi1, phi2);

}

void SourceCircularPolarisation::setCoherent(bool b) {
	coherent = b;
}

void SourceCircularPolarisation::setPhase(double phase) {
	phaseAngle = phase;
}

void SourceCircularPolarisation::setPolarisation(double pol) {
	polarisationAngle = pol;
}


/*********************************************************************************************************************/

SourceALPState::SourceALPState(double a0r, double a0i) : SourceFeature() { 
	std::complex<double> a0 = {a0r, a0i};
	setInitialState(a0);
}

SourceALPState::SourceALPState(std::complex<double> a0) : SourceFeature() { 
	setInitialState(a0);
}

void SourceALPState::setInitialState(std::complex<double> a0) {
	alpField = a0;
}

void SourceALPState::prepareCandidate(Candidate& candidate) const {
	candidate.setProperty(varFieldALP, Variant::fromComplexDouble(alpField));	
	candidate.setProperty(varProbabilityALP, Variant::fromDouble(std::norm(alpField)));
}


} // namespace alpinist