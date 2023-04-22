#ifndef ALPINIST_ALPPHOTONMIXING_H
#define ALPINIST_ALPPHOTONMIXING_H

#include <algorithm>
#include <cmath>
#include <complex>
#include <numeric>

#include <crpropa/Candidate.h>
#include <crpropa/Common.h>
#include <crpropa/Module.h>
#include <crpropa/Random.h>
#include <crpropa/Referenced.h>
#include <crpropa/Units.h>
#include <crpropa/magneticField/MagneticField.h>
#include <Eigen/Core>
#include <Eigen/Eigenvalues>

#include "alpinist/CandidateProperties.h"
#include "alpinist/Constants.h"
#include "alpinist/PlasmaDensity.h"
#include "alpinist/WaveFunction.h"



namespace alpinist {


// forward declarion for latter use
class MixingParameters;


/**
 @class ALPPhotonMixing
 @brief Perform the conversion of an ALP into a photon (and vice-versa).
 */
class ALPPhotonMixing : public Module {
	protected:
		double axionMass; //!< ALP mass [kg]
		double couplingConstant; //!< ALP coupling constant [J^-1]
		bool havePhotons; //!<  whether to add photons to the simulations 
		bool haveALPs; //!<  whether to add ALPs to the simulations 
		double limit; //!< adjust limit of step size in terms of oscillation length
		ref_ptr<PlasmaDensity> plasmaDensity; //!<  density of the plasma in the medium (mostly electrons)
		ref_ptr<MagneticField> magneticField; //!< CRPropa object of type MagneticField
		std::string interactionTag; //!< tag indicating this interaction
	public:
		ALPPhotonMixing(double axionMass, double couplingConstant, ref_ptr<MagneticField> magneticField, ref_ptr<PlasmaDensity> density, double limit = 0.);
		void setAxionMass(double axionMass);
		void setCouplingConstant(double couplingConstant);
		void setMagneticField(ref_ptr<MagneticField> bField);
		void setPlasmaDensity(ref_ptr<PlasmaDensity> density);
		void setInteractionTag(std::string tag);
		void setLimit(double tol);
		double getAxionMass() const;
		double getCouplingConstant() const;
		ref_ptr<PlasmaDensity> getPlasmaDensity() const;
		ref_ptr<MagneticField> getMagneticField() const;
		std::string getInteractionTag() const;
		double getLimit() const;
		WaveFunction3c getWaveFunction(const Candidate& candidate) const;
		void process(Candidate *candidate) const;
		void evolve(Candidate* candidate, MixingParameters& mixing, const WaveFunction3c& initialState, const double& distance) const;
};


/**
 @class MixingParameters
 @brief Class holding information about the mixing.
 */
class MixingParameters {
	public:
		double deltaPlasma;
		double deltaAxion;
		double deltaParallel;
		double deltaPerpendicular;
		double deltaMixing;
		double deltaOscillation;
		double mixingAngle;
		double angleB;
		double energy;
		Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> eigenSolver;

	public:
		MixingParameters();
		MixingParameters(const double& couplingConstant, const double& axionMass, const double& energy, const double& densityMedium, const Vector3d& magneticField);
		void computeParameters(const double& couplingConstant, const double& axionMass, const double& energy, const double& densityMedium, const Vector3d& magneticField);	
		Eigen::Matrix3d getMixingMatrix() const;
		void printParameters();
		void solveEigenSystem(const Eigen::Matrix3d& M);
		Eigen::Vector3d getEigenValues() const;
		Eigen::Matrix3cd getEigenVectors() const;
		Eigen::Matrix3cd getExponentialDiagonalEigenvalueMatrix(const double& x) const;
		Eigen::Matrix3cd getSimilarityTransformingMatrix(const double& x) const;
};




} // namespace alpinist 

#endif