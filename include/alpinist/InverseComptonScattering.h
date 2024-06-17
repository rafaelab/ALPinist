#ifndef ALPINIST_INVERSECOMPTONSCATTERING_H
#define ALPINIST_INVERSECOMPTONSCATTERING_H

#include <cmath>
#include <fstream>
#include <sstream>

#include "crpropa/Module.h"
#include "crpropa/PhotonBackground.h"
#include "crpropa/Random.h"

#include "alpinist/Data.h"
#include "alpinist/Constants.h"



using crpropa::Candidate;
using crpropa::Module;
using crpropa::PhotonField;
using crpropa::Random;
using crpropa::Referenced;
using crpropa::Vector3d;
using crpropa::ref_ptr;
using crpropa::pow_integer;
using crpropa::interpolate;
using crpropa::closestIndex;



namespace alpinist {
/**
 * \addtogroup EnergyLosses
 * @{
 */

/**
 @class InverseComptonScattering
 @brief Inverse Compton scattering of electrons with background photons.

 This module simulates inverse Compton scattering of electrons with background photons for several photon fields.
 The upscattered photons are optionally created as secondary particles (default = false).
 The module limits the propagation step size to a fraction of the mean free path (default = 0.1).
 Thinning is available. A thinning of 0 means that all particles are tracked. 
 For the maximum thinning of 1, only a few representative particles are added to the list of secondaries.
 Note that for thinning>0 the output must contain the column "weights", which should be included in the post-processing.
*/
class InverseComptonScattering: public Module {
	private:
		ref_ptr<PhotonField> photonField;
		bool havePhotons;
		double limit;
		double thinning;
		std::string interactionTag = "EMIC";
		std::array<std::vector<double>, 7> tabEnergy;  //!< electron energy in [J]
		std::array<std::vector<double>, 7> tabRate;  //!< interaction rate in [1/m]
		std::array<std::vector<double>, 7> tabE;  //!< electron energy in [J]
		std::array<std::vector<double>, 7> tabs;  //!< s_kin = s - m^2 in [J**2]
		std::array<std::vector<std::vector<double>>, 7> tabCDF;  //!< cumulative interaction rate

	public:
		InverseComptonScattering(ref_ptr<PhotonField> photonField, bool havePhotons = false, double thinning = 0, double limit = 0.1);
		void setPhotonField(ref_ptr<PhotonField> photonField);
		void setHavePhotons(bool havePhotons);
		void setLimit(double limit);
		void setThinning(double thinning);
		void setInteractionTag(std::string tag);
		std::string getInteractionTag() const;
		void initRate(std::string filename, uint8_t polarisation);
		void initCumulativeRate(std::string filename, uint8_t polarisation);
		void process(Candidate* candidate) const;
		void performInteraction(Candidate* candidate) const;
};
/** @}*/



// Class to calculate the energy distribution of the ICS photon and to sample from it
class ICSSecondariesEnergyDistribution {
	private:
		std::vector< std::vector<double> > data;
		std::vector<double> s_values;
		size_t Ns = 1000;
		size_t Nrer = 1000;
		double s_min = mec2 * mec2;
		double s_max = 2e23 * eV * eV;
		double dls;
		uint8_t polarisation;

	public:
		// create the cumulative energy distribution of the up-scattered photon
		ICSSecondariesEnergyDistribution(uint8_t pol) {
			polarisation = pol;

			dls = (log(s_max) - log(s_min)) / Ns;
			data = std::vector< std::vector<double> >(1000, std::vector<double>(1000));
			std::vector<double> data_i(1000);

			// tabulate s bin borders
			s_values = std::vector<double>(1001);
			for (size_t i = 0; i < Ns + 1; ++i)
				s_values[i] = s_min * exp(i*dls);


			// for each s tabulate cumulative differential cross section
			for (size_t i = 0; i < Ns; i++) {
				double s = s_min * exp((i + 0.5) * dls);
				double beta = (s - s_min) / (s + s_min);
				double x0 = (1 - beta) / (1 + beta);
				double dlx = -log(x0) / Nrer;

				// cumulative midpoint integration
				data_i[0] = dSigmadE(x0, beta) * expm1(dlx);
				for (size_t j = 1; j < Nrer; j++) {
					double x = x0 * exp((j + 0.5) * dlx);
					double dx = exp((j + 1) * dlx) - exp(j * dlx);
					data_i[j] = dSigmadE(x, beta) * dx;
					data_i[j] += data_i[j - 1];
				}
				data[i] = data_i;
			}
		}

		double dSigmadE(double x, double beta) {
			switch (polarisation) {
				case 0:
					return dSigmadE0(x, beta);
				case 1:
					return dSigmadE0(x, beta);
				case 2:
					return dSigmadE0(x, beta);
				case 3:
					return dSigmadE0(x, beta);
				case 4:
					return dSigmadE0(x, beta);
				case 5:
					return dSigmadE0(x, beta);
				default:
					return dSigmadE0(x, beta);
			}
		}

		// differential cross-section, see Lee '96 (arXiv:9604098), eq. 23 for x = Ee'/Ee
		double dSigmadE0(double x, double beta) {
			double q = ((1 - beta) / beta) * (1 - 1. / x);
			return ((1 + beta) / beta) * (x + 1. / x + 2 * q + q * q);
		}

		// draw random energy for the up-scattered photon Ep(Ee, s)
		double sample(double Ee, double s) {
			size_t idx = std::lower_bound(s_values.begin(), s_values.end(), s) - s_values.begin();
			std::vector<double> s0 = data[idx];
			Random& random = Random::instance();
			size_t j = random.randBin(s0) + 1; // draw random bin (upper bin boundary returned)
			double beta = (s - s_min) / (s + s_min);
			double x0 = (1 - beta) / (1 + beta);
			double dlx = -log(x0) / Nrer;
			double binWidth = x0 * (exp(j * dlx) - exp((j - 1) * dlx));
			double Ep = (x0 * exp((j - 1) * dlx) + binWidth) * Ee;
			return std::min(Ee, Ep); // prevent Ep > Ee from numerical inaccuracies
		}
};



} // namespace alpinist

#endif // ALPINIST_EMINVERSECOMPTONSCATTERING_H
