#ifndef ALPINIST_POLARISATIONCORRECTIONS_H
#define ALPINIST_POLARISATIONCORRECTIONS_H


#include <crpropa/Candidate.h>
#include <crpropa/Common.h>
#include <crpropa/Module.h>
#include <crpropa/Random.h>
#include <crpropa/Referenced.h>
#include <crpropa/Units.h>

#include "alpinist/CandidateProperties.h"

namespace alpinist {

using namespace crpropa;

/**
 @class PolarisationCorrections
 @brief Correct the polarisation of secondary particles.
 This must be called if, for instance, the electrons from pair production are allowed to inverse-Compton-scatter background photons. In fact, the pair-produced electrons have specific helicities that are yet to be calculated. 
 For now, only the unpolarised case can be treated due to the changes in cross sections.
 Pair production and inverse Compton scattering have to be rewritten to include helicities.
 */
class PolarisationCorrections : public Module {
	protected:
		int particleId;
		std::string interaction;

	public:
		PolarisationCorrections();
		PolarisationCorrections(int pId, std::string interaction);
		void setParticleId(int pId);
		void setInteraction(std::string interaction);
		int getParticleId() const;
		std::string getInteraction() const;
		void process(Candidate* candidate) const;
};



} // namespace alpinist


#endif 