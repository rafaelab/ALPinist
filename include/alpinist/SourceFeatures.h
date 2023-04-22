#ifndef ALPINIST_SOURCEFEATURES_H
#define ALPINIST_SOURCEFEATURES_H

#include <cmath>
#include <complex>
#include <initializer_list>

#include <crpropa/Candidate.h>
#include <crpropa/Common.h>
#include <crpropa/Random.h>
#include <crpropa/Variant.h>
#include <crpropa/Source.h>


#include "alpinist/CandidateProperties.h"
#include "alpinist/Constants.h"


using namespace crpropa;

namespace alpinist {

/**
 @class SourceNoPolarisation
 @brief Assigns random polarisations to the candidates.
 The EM wave is assumed to be propagating along the x direction.
 The angle wrt to the z axis is theta.
 */
class SourceNoPolarisation: public SourceFeature {
public:
	SourceNoPolarisation();
	void prepareCandidate(Candidate& candidate) const;
};


/**
 @class SourceLinearPolarisation
 @brief Defines the polarisation of the candidate.
  The polarisation angle theta is assumed to be 0 when it is aligned with the z-axis.
 */
class SourceLinearPolarisation: public SourceFeature {
	double polarisationAngle; 
	double phaseAngle;
	bool coherent;
public:
	SourceLinearPolarisation(double angle, bool coherent = true, double phase = 0);
	void prepareCandidate(Candidate& candidate) const;
	void setPhase(double phase);
	void setCoherent(bool b);
	void setPolarisation(double angle);
};


/**
 @class SourceCircularPolarisation
 @brief Defines the (circular) polarisation of the wave
  The polarisation angle theta is assumed to be 0 when it is aligned with the z-axis.
  The phase angle can be constant, but this is redundant with the polarisation angle.
 */
class SourceCircularPolarisation: public SourceFeature {
	double polarisationAngle;
	double phaseAngle;
	bool coherent;
public:
	SourceCircularPolarisation(double angle, bool coherent = true, double phase = 0);
	void prepareCandidate(Candidate& candidate) const;
	void setPhase(double phase);
	void setCoherent(bool b);
	void setPolarisation(double angle);
};


/**
 @class SourceALPState
 @brief Sets the initial conditions for the ALP field.
  By default, a(0)=0.
 */
class SourceALPState: public SourceFeature {
private:
	std::complex<double> alpField;
public:
	SourceALPState(double realALPField = 0, double imagALPField = 0);
	SourceALPState(std::complex<double> alp);
	void setInitialState(std::complex<double> a0);
	void prepareCandidate(Candidate& candidate) const;
};


/**
 Convenient method to initialise the polarisation of candidates.
*/
inline void setPolarisedCandidate(Candidate& candidate, const double& theta, const double& phi1, const double& phi2) {
	std::complex<double> em1 = {cos(theta) * cos(phi1), cos(theta) * sin(- phi1)};
	std::complex<double> em2 = {sin(theta) * cos(phi2), sin(theta) * sin(- phi2)};
	candidate.setProperty(varFieldEM1, Variant::fromComplexDouble(em1));
	candidate.setProperty(varFieldEM2, Variant::fromComplexDouble(em2));
	candidate.setProperty(varProbabilityPhoton, Variant::fromDouble(1.));
}


} // namespace alpinist

#endif