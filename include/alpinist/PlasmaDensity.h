#ifndef ALPINIST_PLASMADENSITY_H
#define ALPINIST_PLASMADENSITY_H

#include <algorithm>
#include <iostream>
#include <stdexcept>

#include <crpropa/Cosmology.h>
#include <crpropa/Common.h>
#include <crpropa/Grid.h>
#include <crpropa/GridTools.h>
#include <crpropa/Referenced.h>

#include "alpinist/Constants.h"

using namespace crpropa;


namespace alpinist {

/**
 @class PlasmaDensity
 @brief Generic object describing the plasma density of the medium.
 It inherits from `Referenced` to be called as a `ref_ptr`.
 */
class PlasmaDensity : public Referenced {
	public:
		virtual ~PlasmaDensity() = default;
		virtual double getDensity(const Vector3d& position, const double& z = 0) const = 0;
};


/**
 @class PlasmaDensityUniform
 @brief Constant density of the medium.
 */
class PlasmaDensityUniform : public PlasmaDensity {
	protected:
		double density;
	public:
		PlasmaDensityUniform(double density);
		~PlasmaDensityUniform();
		void setDensityValue(double n);
		double getDensityValue() const;
		double getDensity(const Vector3d& position, const double& z = 0.) const;
};


/**
 @class PlasmaDensityGrid
 @brief Plasma density of medium position-dependent, given by the grid.
 */
class PlasmaDensityGrid : public PlasmaDensity {
	protected:
		ref_ptr<Grid1f> densityGrid;
		double normalisation;
	public:
		PlasmaDensityGrid();
		PlasmaDensityGrid(const ref_ptr<Grid1f> &grid, const double& norm = 1.);
		~PlasmaDensityGrid();
		void setNormalisation(double n);
		void setGrid(const ref_ptr<Grid1f>& grid);
		ref_ptr<Grid1f> getGrid() const;
		double getNormalisation() const;
		double getDensity(const Vector3d& position, const double &z = 0.) const;
};

} // namespace alpinist



/**
 Estimate the plasma frequency (in s^-1).
 */
inline double computePlasmaFrequency(double densityMedium) {
	return sqrt(densityMedium * eplus * eplus / mass_electron / epsilon0);
}


#endif // ALPINIST_PLASMADENSITY_H