//read all kinds of TrajReader file
#include"trajectory.h"

#define DATA_LENGTH		24
#define LINE_LENGTH		62
#define MAX_ATOM_TYPE	256

inline int approx(double d) {
	if (d - (int)d > 0.5) {
		return (int)d + 1;
	}
	else {
		return (int)d;
	}
}

template<class Name, class Val>
inline const Val& get(const map<Name, Val>& m, const Name& n, const Val& d) {
	auto k = m.find(n);
	if (k == m.end()) {
		return d;
	}
	else {
		return k->second;
	}
}


int TrajReader::Open(const string& filename) {
	trjfile.open(filename, ios_base::in);
	if (!trjfile.is_open()) {
		return 1;
	}
	else {
		return 0;
	}
}

int TrajReader::ReadTrjHead(Simulation* simulation) {

	//read file head
	char buffer[LINE_LENGTH];
	int traj_length;
	trjfile >> buffer >> traj_length;
	trjfile >> buffer >> simulation->name;
	trjfile >> buffer >> simulation->atomNumber;
	trjfile.seekg(2 * LINE_LENGTH, ios_base::cur);

	if (simulation->timeStep == 0.0) {
		trjfile >> buffer >> simulation->timeStep;
		simulation->timeStep /= 1000.0;
	}
	else {
		trjfile >> buffer >> buffer;
	}

	printf("Name: %s\n", simulation->name);
	printf("Atom numbers: %d\n", simulation->atomNumber);
	printf("Time step: %.2e ps\n", simulation->timeStep);

	trjfile.seekg(traj_length, ios_base::beg);
	//read atom type
	trjfile.seekg(2 * LINE_LENGTH, ios_base::cur);

	atomTypes.resize(simulation->atomNumber + 1);
	atomWeights.resize(simulation->atomNumber + 1);
	for (int i = 0; i < simulation->atomNumber; i++) {
		int atomIndex, atomType; double atomWeight;
		trjfile >> atomIndex >> atomType >> atomWeight;
		atomTypes[atomIndex] = atomType;
		atomWeights[atomIndex] = atomWeight;
	}

	//record atom types
	for (auto w = atomWeights.begin() + 1; w != atomWeights.end(); w++) {
		if (find(simulation->atomWeights.begin(), simulation->atomWeights.end(), *w) == simulation->atomWeights.end())
			simulation->atomWeights.push_back(*w);
	}
	return 0;
}

int TrajReader::ReadTrjFrame(Frame& frameOut) {
	if (trjfile.peek() == EOF)return 0;
	
	int frameHeadLength;
	char buffer[LINE_LENGTH];
	memset(buffer, 0, LINE_LENGTH);

	trjfile >> buffer >> frameHeadLength;
	if (buffer[0] == '\0')return 0;
	/*	infile.seekg(frameHeadLength, ios_base::cur);	*/

	for (int i = 0; i < frameHeadLength / LINE_LENGTH + 2; i++)
		trjfile.getline(buffer, LINE_LENGTH, '\n');

	if (config.read_atompos) {
		frameOut.atoms.clear();
		for (int i = 0; i < atomTypes.size() - 1; i++) {
			atom crtAtom; double q;
			trjfile >> crtAtom.id >> crtAtom.x >> crtAtom.y >> crtAtom.z >> q;
			frameOut.atoms.push_back(crtAtom);
		}
	}
	else {
		for (int i = 0; i < atomTypes.size() - 1; i++) {
			trjfile.getline(buffer, LINE_LENGTH, '\n');
		}
	}

	int bondNumber;
	trjfile.getline(buffer, LINE_LENGTH, ',');
	trjfile >> bondNumber;
	frameOut.bonds.clear();
	for (int i = 0; i < bondNumber; i++) {
		bond crtBond; double length, raw_order;
		trjfile >> crtBond.id_1 >> crtBond.id_2 >> length >> raw_order;

		int type1 = int(atomWeights[crtBond.id_1] + 0.1), type2 = int(atomWeights[crtBond.id_2] + 0.1);
		int bond_identifier = type1 > type2 ? (type1 * MAX_ATOM_TYPE + type2) : (type2 * MAX_ATOM_TYPE + type1);

		crtBond.order = floor_vec(get(config.bondorder_cutoff, bond_identifier, config.bondorder_cutoff_default), raw_order);

		if (crtBond.order >= 0) {
			frameOut.bonds.push_back(crtBond);
		}
		if (config.count_bondorder) {
			bondorders[bond_identifier].push_back(raw_order);
		}
	}
	return 1;
}
