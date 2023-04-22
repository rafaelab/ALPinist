#include "alpinist/PlasmaDensity.h"

namespace alpinist {

/*********************************************************************************************************************/

PlasmaDensityUniform::PlasmaDensityUniform(double n) {
	setDensityValue(n);
}

PlasmaDensityUniform::~PlasmaDensityUniform() {
}

void PlasmaDensityUniform::setDensityValue(double n) {
	density = n;
}

double PlasmaDensityUniform::getDensityValue() const {
	return density;
}

double PlasmaDensityUniform::getDensity(const Vector3d& position, const double& z) const {
	return density * pow_integer<3>(1 + z);
}
	

/*********************************************************************************************************************/

PlasmaDensityGrid::PlasmaDensityGrid() {
}

PlasmaDensityGrid::PlasmaDensityGrid(const ref_ptr<Grid1f>& grid, const double& norm) {
	setGrid(grid);
	setNormalisation(norm);
}

PlasmaDensityGrid::~PlasmaDensityGrid() {
}

void PlasmaDensityGrid::setGrid(const ref_ptr<Grid1f>& grid) {
	densityGrid = grid;
}

void PlasmaDensityGrid::setNormalisation(double norm) {
	normalisation = norm;
}

ref_ptr<Grid1f> PlasmaDensityGrid::getGrid() const {
	return densityGrid;
}

double PlasmaDensityGrid::getNormalisation() const {
	return normalisation;
}

double PlasmaDensityGrid::getDensity(const Vector3d& position, const double& z) const {
	return normalisation * densityGrid->interpolate(position) * pow_integer<3>(1 + z);
}





} // namespace alpinist