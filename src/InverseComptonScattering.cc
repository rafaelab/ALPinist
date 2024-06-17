#include "alpinist/InverseComptonScattering.h"


namespace alpinist {


InverseComptonScattering::InverseComptonScattering(ref_ptr<PhotonField> photonField, bool havePhotons, double thinning, double limit) {
	setPhotonField(photonField);
	setHavePhotons(havePhotons);
	setLimit(limit);
	setThinning(thinning);
}

void InverseComptonScattering::setPhotonField(ref_ptr<PhotonField> field) {
	photonField = field;

	std::string label = photonField->getFieldName();
	setDescription("InverseComptonScattering: " + label);

	initRate(getDataPath("InverseComptonScattering/rate_" + label + ".txt"), 0);
	initCumulativeRate(getDataPath("InverseComptonScattering/cdf_" + label + ".txt"), 0);

	// there are 6 polarisation modes
	for (size_t i = 1; i <= 6; i++) {
		std::ostringstream ssDif;
		std::ostringstream ssCum;
     	ssDif << "InverseComptonScattering/rate_pol" << i << "_" << label << ".txt";
		ssCum << "InverseComptonScattering/cdf_pol" << i << "_" << label << ".txt";
		std::string fnDif = getDataPath(ssDif.str());
		std::string fnCum = getDataPath(ssCum.str());
		initRate(fnDif, i);
		initCumulativeRate(fnCum, i);
	}
}

void InverseComptonScattering::setHavePhotons(bool photons) {
	havePhotons = photons;
}

void InverseComptonScattering::setLimit(double l) {
	limit = l;
}

void InverseComptonScattering::setThinning(double t) {
	thinning = t;
}

void InverseComptonScattering::setInteractionTag(std::string tag) {
	interactionTag = tag;
}

std::string InverseComptonScattering::getInteractionTag() const {
	return interactionTag;
}

void InverseComptonScattering::initRate(std::string filename, uint8_t polarisation) {
	std::ifstream infile(filename.c_str());

	if (! infile.good())
		throw std::runtime_error("InverseComptonScattering: could not open file " + filename);

	// clear previously loaded tables
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

void InverseComptonScattering::initCumulativeRate(std::string filename, uint8_t polarisation) {
	std::ifstream infile(filename.c_str());

	if (! infile.good())
		throw std::runtime_error("InverseComptonScattering: could not open file " + filename);

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

void InverseComptonScattering::performInteraction(Candidate* candidate) const {
	uint8_t polarisation = 0;

	// scale the particle energy instead of background photons
	double z = candidate->getRedshift();
	double E = candidate->current.getEnergy() * (1 + z);

	if (E < tabE[polarisation].front() or E > tabE[polarisation].back())
		return;

	// sample the value of s
	Random& random = Random::instance();
	size_t i = closestIndex(E, tabE[polarisation]);
	size_t j = random.randBin(tabCDF[polarisation][i]);
	double s_kin = pow(10, log10(tabs[polarisation][j]) + (random.rand() - 0.5) * 0.1);
	double s = s_kin + mec2 * mec2;

	// sample electron energy after scattering
	static ICSSecondariesEnergyDistribution distribution(polarisation);
	double Enew = distribution.sample(E, s);

	// add up-scattered photon
	if (havePhotons) {
		double Esecondary = E - Enew;
		double f = Enew / E;
		if (random.rand() < pow(1 - f, thinning)) {
			double w = 1. / pow(1 - f, thinning);
			Vector3d pos = random.randomInterpolatedPosition(candidate->previous.getPosition(), candidate->current.getPosition());
			candidate->addSecondary(22, Esecondary / (1 + z), pos, w, interactionTag);
		}
	}

	// update the primary particle energy; do this after adding the secondary to correctly set the secondary's parent
	candidate->current.setEnergy(Enew / (1 + z));
}

void InverseComptonScattering::process(Candidate* candidate) const {
	// check if electron / positron
	int id = candidate->current.getId();
	if (abs(id) != 11)
		return;

	// scale the particle energy instead of background photons
	double z = candidate->getRedshift();
	double E = candidate->current.getEnergy() * (1 + z);

	uint8_t polarisation = 0;

	if (E < tabEnergy[polarisation].front() or (E > tabEnergy[polarisation].back()))
		return;

	// interaction rate
	double rate = interpolate(E, tabEnergy[polarisation], tabRate[polarisation]);
	rate *= pow_integer<2>(1 + z) * photonField->getRedshiftScaling(z);

	// run this loop at least once to limit the step size
	double step = candidate->getCurrentStep();
	Random& random = Random::instance();
	do {
		double randDistance = -log(random.rand()) / rate;

		// check for interaction; if it doesn't ocurr, limit next step
		if (step < randDistance) {
			candidate->limitNextStep(limit / rate);
			return;
		}
		performInteraction(candidate);

		// repeat with remaining step
		step -= randDistance;
	} while (step > 0);
}


} // namespace crpropa
