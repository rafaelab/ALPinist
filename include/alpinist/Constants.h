#ifndef ALPINIST_CONSTANTS_H
#define ALPINIST_CONSTANTS_H

#include <cmath>

#include <crpropa/Common.h>
#include <crpropa/Units.h>


namespace alpinist {

using namespace crpropa;

static const double h_dirac = h_planck / 2. / M_PI;
static const double fineStructure = eplus * eplus / 4. / M_PI / epsilon0 / h_dirac / c_light;
static const double radius_bohr = h_dirac * h_dirac / mass_electron / eplus;
static const double B_schwinger = pow_integer<2>(mass_electron * c_light) / eplus / h_dirac;
static const double E_schwinger = B_schwinger * c_light;


} // namespace


#endif // ALPINIST_CONSTANTS_H