%module(directors = "1", threads = "1", allprotected = "1") alpinist


%include "attribute.i"
%include "exception.i"
%include "stdint.i"
%include "std_array.i"
%include "std_complex.i"
%include "std_container.i"
%include "std_iostream.i"
%include "std_list.i"
%include "std_map.i"
%include "std_set.i"
%include "std_shared_ptr.i"
%include "std_string.i"
%include "std_vector.i"
%include "stl.i"
%include "typemaps.i"


/* Ignore list */
%ignore operator alpinist::PlasmaDensity*;


/* Headers */
%{
	#include "CRPropa.h"
	#include "alpinist/Constants.h"
	#include "alpinist/Data.h"
	#include "alpinist/CandidateProperties.h"
	#include "alpinist/SourceFeatures.h"
	#include "alpinist/PlasmaDensity.h"
	#include "alpinist/WaveFunction.h"
	#include "alpinist/PolarisationCorrections.h"
	#include "alpinist/PairProduction.h"
	#include "alpinist/InverseComptonScattering.h"
	#include "alpinist/ALPPhotonMixing.h"
	
	using namespace alpinist;
%}


/* Import CRPropa in wrapper */
%import (module = "crpropa") "crpropa.i"


/* Templates for ref_ptr */
%implicitconv crpropa::ref_ptr<alpinist::PlasmaDensity>;
%template(PlasmaDensityRefPtr) crpropa::ref_ptr<alpinist::PlasmaDensity>;
%feature("director") alpinist::PlasmaDensity;


/* Include plugin parts to generate wrappers for */
%include "alpinist/Constants.h"
%include "alpinist/Data.h"
%include "alpinist/CandidateProperties.h"
%include "alpinist/SourceFeatures.h"
%include "alpinist/PlasmaDensity.h"
%include "alpinist/WaveFunction.h"
%include "alpinist/PolarisationCorrections.h"
%include "alpinist/PairProduction.h"
%include "alpinist/InverseComptonScattering.h"
%include "alpinist/ALPPhotonMixing.h"



/* Hide warnings */
#pragma SWIG nowarn=312,325,361,389,401,508,509
