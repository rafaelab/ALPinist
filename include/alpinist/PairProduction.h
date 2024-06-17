#ifndef ALPINIST_PAIRPRODUCTION_H
#define ALPINIST_PAIRPRODUCTION_H

#include <array>
#include <cmath>
#include <fstream>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>


#include <crpropa/Candidate.h>
#include <crpropa/Common.h>
#include <crpropa/Module.h>
#include <crpropa/PhotonBackground.h>
#include <crpropa/Random.h>
#include <crpropa/Referenced.h>
#include <crpropa/Units.h>
#include <crpropa/Vector3.h>

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
 @class PairProduction
 @brief Breit-Wheeler pair production of photons with background photons considering polarisation effects.

 This module simulates electron-pair production of photons with background photons:
 	gamma + gamma_b --> e+ + e- 
 The resulting electron-positron pair is optionally created (default = false).
 The module limits the propagation step size to a fraction of the mean free path (default = 0.1).
 Thinning is available. A thinning of 0 means that all particles are tracked. 
 For the maximum thinning of 1, only a few representative particles are added to the list of secondaries.
 Note that for thinning>0 the output must contain the column "weights", which should be included in the post-processing.
 */
class PairProduction: public Module {
	protected:
		ref_ptr<PhotonField> photonField; 	// target photon field
		bool haveElectrons; // add secondary electrons to simulation
		double limit; // limit the step to a fraction of the mean free path
		double thinning; // factor of the thinning (0: no thinning, 1: maximum thinning)
		std::string interactionTag;  //!< tag corresponding to this interaction
		// std::vector<double> tabEnergy;  //!< electron energy in [J]
		// std::vector<double> tabRate;  //!< interaction rate in [1/m]
		// std::vector<double> tabE;  //!< electron energy in [J]
		// std::vector<double> tabs;  //!< s_kin = s - m^2 in [J**2]
		// std::vector<std::vector<double>> tabCDF;  //!< cumulative interaction rate
		std::array<std::vector<double>, 5> tabEnergy;  //!< electron energy in [J]
		std::array<std::vector<double>, 5> tabRate;  //!< interaction rate in [1/m]
		std::array<std::vector<double>, 5> tabE;  //!< electron energy in [J]
		std::array<std::vector<double>, 5> tabs;  //!< s_kin = s - m^2 in [J**2]
		std::array<std::vector<std::vector<double>>, 5> tabCDF;  //!< cumulative interaction rate
		
	public:
		PairProduction(ref_ptr<PhotonField> photonField,  bool haveElectrons = false, double thinning = 0, double limit = 0.1);
		void setPhotonField(ref_ptr<PhotonField> photonField);
		void setHaveElectrons(bool haveElectrons);
		void setLimit(double limit);
		void setThinning(double thinning);
		void setInteractionTag(std::string tag);
		std::string getInteractionTag() const;
		void initRate(std::string filename, uint8_t polarisation);
		void initCumulativeRate(std::string filename, uint8_t polarisation);
		void performInteraction(Candidate* candidate) const;
		void process(Candidate* candidate) const;
};


// Hold an data array to interpolate the energy distribution 
class PPSecondariesEnergyDistribution {
	private:
		std::vector<double> tab_s;
		std::vector<std::vector<double>> data;
		uint8_t polarisation;
		size_t N = 1000;

	public:
		PPSecondariesEnergyDistribution(int pol) {
			polarisation = pol;

			size_t Ns = 1000;
			double s_min = 4 * mec2 * mec2;
			double s_max = 1e23 * eV * eV;
			double dls = log(s_max / s_min) / Ns;
			data = std::vector< std::vector<double> >(Ns, std::vector<double>(N));
			tab_s = std::vector<double>(Ns + 1);

			for (size_t i = 0; i < Ns + 1; ++i)
				tab_s[i] = s_min * exp(i * dls); // tabulate s bin borders

			for (size_t i = 0; i < Ns; i++) {
				double s = s_min * exp(i * dls + 0.5 * dls);
				double beta = sqrt(1 - s_min/s);
				double x0 = (1. - beta) / 2;
				double dx = log((1 + beta) / (1 - beta)) / N;

				// cumulative midpoint integration
				std::vector<double> data_i(1000);
				data_i[0] = dSigmadE(x0, beta) * expm1(dx);
				for (size_t j = 1; j < N; j++) {
					double x = x0 * exp(j * dx + 0.5 * dx);
					double binWidth = exp((j + 1) * dx) - exp(j * dx);
					data_i[j] = dSigmadE(x, beta) * binWidth + data_i[j - 1];
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

		// differential cross section for pair production for x = Epositron/Egamma, compare Lee 96 arXiv:9604098
		double dSigmadE0(double x, double beta) {
			double A = (x / (1. - x) + (1. - x) / x);
			double B =  (1. / x + 1. / (1. - x));
			double y = (1 - beta * beta);
			return A + y * B - y * y / 4 * B * B;
		}

		// sample positron energy from cdf(E, s_kin)
		double sample(double E0, double s) {
			// get distribution for given s
			size_t idx = std::lower_bound(tab_s.begin(), tab_s.end(), s) - tab_s.begin();
			if (idx > data.size())
				return NAN;
				
			std::vector<double> s0 = data[idx];

			// draw random bin
			Random& random = Random::instance();
			size_t j = random.randBin(s0) + 1;

			double s_min = 4. * mec2 * mec2;
			double beta = sqrtl(1. - s_min / s);
			double x0 = (1. - beta) / 2.;
			double dx = log((1 + beta) / (1 - beta)) / N;
			double binWidth = x0 * (exp(j * dx) - exp((j - 1) * dx));
			if (random.rand() < 0.5)
				return E0 * (x0 * exp((j - 1) * dx) + binWidth);
			else
				return E0 * (1 - (x0 * exp((j - 1) * dx) + binWidth));
		}
};



} // namespace alpinist

#endif // ALPINIST_PAIRPRODUCTION_H
