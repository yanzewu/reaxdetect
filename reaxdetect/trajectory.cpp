//read all kinds of TrajReader file
#include"trajectory.h"

#define DATA_LENGTH	24
#define LINE_LENGTH	62

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
	trjfile >> buffer >> simulation->timeStep;
	simulation->timeStep /= 1000.0;

	printf("Name: %s\n", simulation->name);
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

	trjfile >> buffer >> frameHeadLength;
	/*	infile.seekg(frameHeadLength, ios_base::cur);	*/

	for (int i = 0; i < frameHeadLength / LINE_LENGTH + 2; i++)
		trjfile.getline(buffer, LINE_LENGTH, '\n');

	if (config.read_atompos) {
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
	for (int i = 0; i < bondNumber; i++) {
		bond crtBond; double length, order;
		trjfile >> crtBond.id_1 >> crtBond.id_2 >> length >> order;
		if (order > config.bondorder_cutoff)frameOut.bonds.push_back(crtBond);
	}
	return 1;
}


/*
int TrajReader::ReadTrjFile(const string path, ReaxDetect::GlobalData* simulation, ReaxDetect::Config* cfg) {
	// updated in 2.2.2, separate details in two other functions.

	// open file
	trjfile.open(path, ios_base::in);
	if (!trjfile.is_open())return 1;

	// read head
	if (ReadTrjHead(path, simulation))return 1;

	// read frames
	int cnt = 0;
	int frmNumber = (int)(cfg->valueOf("FrameNumber"));
	double bondOrderCutoff = cfg->valueOf("BondOrderCutoff");
	int readAtomPos = (int)cfg->valueOf("ReadAtomPos");

	while (trjfile.peek() != EOF && (frmNumber == 0 || cnt < frmNumber)) {
		printf("Reading frame %d\r", cnt + 1);
		frame f;
		ReadTrjFrame(f, bondOrderCutoff, readAtomPos);
		if(cnt > 0)frames.push_back(f);	//updated in 2.2.0, to remove uncorrect molecule pack in first frame.
		cnt++;
	}

	trjfile.close();
	frames.pop_back(); //last element will cause uncertainty
	return 0;
}*/