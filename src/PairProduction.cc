#include "alpinist/PairProduction.h"


namespace alpinist {


PairProduction::PairProduction(ref_ptr<PhotonField> photonField, bool haveElectrons, double thinning, double limit) {
	setPhotonField(photonField);
	setThinning(thinning);
	setLimit(limit);
	setHaveElectrons(haveElectrons);
}

void PairProduction::setPhotonField(ref_ptr<PhotonField> photonField) {
	photonField = photonField;

	std::string label = photonField->getFieldName();
	setDescription("PairProduction: " + label);

	initRate(getDataPath("PairProduction/rate_" + label + ".txt"), 0);
	initCumulativeRate(getDataPath("PairProduction/cdf_" + label + ".txt"), 0);

	// there are 5 polarisation modes, one of which has 0 cross section
	for (size_t i = 1; i < 5; i++) {
		std::ostringstream ssDif;
		std::ostringstream ssCum;
     	ssDif << "PairProduction/rate_pol" << i << "_" << label << ".txt";
		ssCum << "PairProduction/cdf_pol" << i << "_" << label << ".txt";
		std::string fnDif = getDataPath(ssDif.str());
		std::string fnCum = getDataPath(ssCum.str());
		initRate(fnDif, i);
		initCumulativeRate(fnCum, i);
	}
}

void PairProduction::setHaveElectrons(bool electrons) {
	haveElectrons = electrons;
}

void PairProduction::setLimit(double l) {
	limit = l;
}

void PairProduction::setThinning(double t) {
	thinning = t;
}

void PairProduction::setInteractionTag(std::string tag) {
	interactionTag = tag;
}

std::string PairProduction::getInteractionTag() const {
	return interactionTag;
}

void PairProduction::initRate(std::string filename, uint8_t polarisation) {
	std::ifstream infile(filename.c_str());

	if (! infile.good())
		throw std::runtime_error("PairProduction: could not open file " + filename);

	// clear previously loaded interaction rates
	tabEnergy[polarisation].clear();
	tabRate[polarisation].clear();

	while (infile.good()) {
		if (infile.peek() != '#') {
			double a, b;
			infile >> a >> b;
			if (infile) {
				tabEnergy[polarisation].push_back(pow(10, a) * eV);
				tabRate[polarisation].push_back(b / Mpc);
			}
		}
		infile.ignore(std::numeric_limits <std::streamsize> ::max(), '\n');
	}
	infile.close();
}

void PairProduction::initCumulativeRate(std::string filename, uint8_t polarisation) {
	std::ifstream infile(filename.c_str());

	if (! infile.good())
		throw std::runtime_error("PairProduction: could not open file " + filename);

	// clear previously loaded tables
	tabE[polarisation].clear();
	tabs[polarisation].clear();
	tabCDF[polarisation].clear();
	
	// skip header
	while (infile.peek() == '#')
		infile.ignore(std::numeric_limits <std::streamsize> ::max(), '\n');

	// read s values in first line
	double a;
	infile >> a; // skip first value
	while (infile.good() and (infile.peek() != '\n')) {
		infile >> a;
		tabs[polarisation].push_back(pow(10, a) * eV * eV);
	}

	// read all following lines: E, cdf values
	while (infile.good()) {
		infile >> a;
		if (! infile)
			break;  // end of file
		tabE[polarisation].push_back(pow(10, a) * eV);
		std::vector<double> cdf;
		for (int i = 0; i < tabs.size(); i++) {
			infile >> a;
			cdf.push_back(a / Mpc);
		}
		tabCDF[polarisation].push_back(cdf);
	}
	infile.close();
}

void PairProduction::performInteraction(Candidate *candidate) const {
	
	// photon is lost after interacting
	candidate->setActive(false);


	uint8_t polarisation = 0;

	// check if secondary electron pair needs to be produced
	if (not haveElectrons)
		return;
		
	// scale particle energy instead of background photon energy
	double z = candidate->getRedshift();
	double E = candidate->current.getEnergy() * (1 + z);

	// check if in tabulated energy range
	if (E < tabE[polarisation].front() or (E > tabE[polarisation].back()))
		return;

	// sample the value of s
	Random& random = Random::instance();
	size_t i = closestIndex(E, tabE[polarisation]);  // find closest tabulation point
	size_t j = random.randBin(tabCDF[polarisation][i]);
	double lo = std::max(4 * mec2 * mec2, tabs[polarisation][j - 1]);  // first s-tabulation point below min(s_kin) = (2 me c^2)^2; ensure physical value
	double hi = tabs[polarisation][j];
	double s = lo + random.rand() * (hi - lo);

	// sample electron / positron energy
	static PPSecondariesEnergyDistribution interpolation(0);
	double Ee = interpolation.sample(E, s);
	double Ep = E - Ee;
	double f = Ep / E;

	// for some backgrounds Ee=nan due to precision limitations.
	if (not std::isfinite(Ee) || not std::isfinite(Ep))
		return;

	// sample random position along current step
	Vector3d pos = random.randomInterpolatedPosition(candidate->previous.getPosition(), candidate->current.getPosition());
	// apply sampling
	if (random.rand() < pow(f, thinning)) {
		double w = 1. / pow(f, thinning);
		candidate->addSecondary(11, Ep / (1 + z), pos, w, interactionTag);
	}
	if (random.rand() < pow(1 - f, thinning)){
		double w = 1. / pow(1 - f, thinning);
		candidate->addSecondary(-11, Ee / (1 + z), pos, w, interactionTag);	
	}
}

void PairProduction::process(Candidate *candidate) const {
	// check if photon
	if (candidate->current.getId() != 22)
		return;

	uint8_t polarisation = 0;

	// scale particle energy instead of background photon energy
	double z = candidate->getRedshift();
	double E = candidate->current.getEnergy() * (1 + z);

	// check if in tabulated energy range
	if ((E < tabEnergy[polarisation].front()) or (E > tabEnergy[polarisation].back()))
		return;

	// interaction rate
	double rate = interpolate(E, tabEnergy[0], tabRate[0]);
	rate *= pow_integer<2>(1 + z) * photonField->getRedshiftScaling(z);

	// run this loop at least once to limit the step size 
	double step = candidate->getCurrentStep();
	Random& random = Random::instance();
	do {
		double randDistance = -log(random.rand()) / rate;
		// check for interaction; if it doesn't occur, limit next step
		if (step < randDistance) { 
			candidate->limitNextStep(limit / rate);
		} else {
			performInteraction(candidate);
			return;
		}
		step -= randDistance; 
	} while (step > 0.);

}



} // namespace alpinist
