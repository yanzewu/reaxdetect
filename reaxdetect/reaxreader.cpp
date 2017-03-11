//contains algorithm about molecule and Reaction generate
#include "reaxreader.h"
#include "smiles.h"
#include "filenamehandler.h"
#include "listhandler.h"

#define FLAG_REAGANT 0
#define FLAG_PRODUCT 1

//---------------top api----------------------

int ReaxReader::HandleData(TrajReader& reader, const Simulation& simulation)
{
	_prev_buffer_idx = 1;
	_crt_buffer_idx = 0;
	_mol_of_atom[0].resize(simulation.atomNumber + 1);
	_mol_of_atom[1].resize(simulation.atomNumber + 1);
	_bond_matrix.resize(simulation.atomNumber + 1);
	_atom_score.resize(simulation.atomNumber + 1);
	
	TrajReader::Frame _frame;
	int i = 0;
	while (reader.ReadTrjFrame(_frame)) {
		if (i % 1000 == 0) {
			printf("Reading frame %d\r", i);
		}
		FrameStat fstat;
		RecognizeMolecule(_frame, fstat, _crt_buffer_idx, simulation.atomNumber, reader.atomWeights);
		if (i > 0)RecognizeReaction(_frame, fstat, _crt_buffer_idx, _prev_buffer_idx);
		fss.push_back(fstat);
		SwapBuffer();
		i++;
	}

	for (auto& mol : molecules) {
		mol.compile();
		species.push_back(mol.toString());
	}

	//alignment
	for (auto& fs : fss) {
		fs.mol_freq.resize(molecules.size());
		fs.reaction_freq.resize(reactions.size());
	}
	return 0;
}


//----------------middle api------------------

void ReaxReader::SwapBuffer() {
	swap(_crt_buffer_idx, _prev_buffer_idx);
	_mol_buffer[_crt_buffer_idx].clear();
	_mol_index_buffer[_crt_buffer_idx].clear();
	_frame_buffer.atoms.clear();
	_frame_buffer.bonds.clear();

	for (auto& m : _mol_of_atom[_crt_buffer_idx])
		m = -1;
	for (auto& a : _atom_score)
		a = 0;
	for (auto& b : _bond_matrix)
		b.clear();
}
void ReaxReader::RecognizeMolecule(const TrajReader::Frame& frm, FrameStat& fs, tsize_t crtBufferIndex, int atomNumber, const vector<double>& atomWeights)
{
	for (auto bnd = frm.bonds.begin(); bnd != frm.bonds.end(); bnd++) {
		_bond_matrix[bnd->id_1].push_back(bnd->id_2); _bond_matrix[bnd->id_2].push_back(bnd->id_1);
	}
	//1.first scan(DFS): get score, root, and mol of atoms: O(N)
	Array molRoots;//root of molecules
	bool* marks = new bool[atomNumber + 1]{ false };
	for (auto i = 1; i <= atomNumber; i++) {
		if (marks[i])continue;
		molRoots.push_back(i);
		scan_score(i, _bond_matrix, atomWeights, marks, molRoots, _atom_score, _mol_of_atom[crtBufferIndex]);
	}
	//2.sort bondmatrix: O(N*k^2)
	for (auto bonds = _bond_matrix.begin() + 1; bonds != _bond_matrix.end(); bonds++) {
		sortby(*bonds, _atom_score, greater<int>());
	}
	memset(marks, 0, sizeof(bool)*(atomNumber + 1));

	//3.second scan(DFS): get molecule list: O(N)
	for (const auto& root : molRoots) {
		MatMolecule molecule(&_bond_matrix, &_atom_score);
		scan_molecule(root, _bond_matrix, marks, molecule);
		_mol_buffer[crtBufferIndex].push_back(molecule);
	}
	delete[] marks;
	
	//4.compress trees: O(N*m)
	vector<tsize_t> uroots;		//index in molbuffer
	Array uFreq;		//freq corresponding to molecules.
	_mol_index_buffer[crtBufferIndex].resize(molRoots.size());

	for (tsize_t i = 0; i < _mol_buffer[crtBufferIndex].size(); i++) {
		_mol_index_buffer[crtBufferIndex][i] = pushnew2index(i, _mol_buffer[crtBufferIndex], uroots, uFreq);
	}

	//5.final: get smiles: O(m^2*log(m))
	Array uindex(uroots.size());
	fs.mol_freq.resize(molecules.size());
	for (tsize_t i = 0; i < uroots.size(); i++) {
		smiles newsms = _mol_buffer[crtBufferIndex][uroots[i]].to_smiles();
			
		auto findpos = find(molecules.begin(), molecules.end(), newsms);
		if(findpos != molecules.end()){
			fs.mol_freq[findpos - molecules.begin()] += uFreq[i];
			uindex[i] = findpos - molecules.begin();
		}
		else {
			molecules.push_back(newsms);
			fs.mol_freq.push_back(uFreq[i]);
			uindex[i] = molecules.size() - 1;
		}
	}
	//6.restore uindex to index: O(N)
	for (auto& molindex : _mol_index_buffer[crtBufferIndex]) {
		molindex = uindex[molindex];
	}
}
void ReaxReader::RecognizeReaction(const TrajReader::Frame& frm, FrameStat& fs, tsize_t crtBufferIndex, tsize_t prevBufferIndex)
{
	vector<bool> marks[2]{ vector<bool>(_mol_buffer[prevBufferIndex].size(), false),
		vector<bool>(_mol_buffer[crtBufferIndex].size(), false) }; //marks of product and reagant
	fs.reaction_freq.resize(reactions.size());

	for (tsize_t i = 0; i < _mol_buffer[crtBufferIndex].size(); i++) {
		if (marks[1][i])continue;
		Reaction reaction;
		create_reaction(i, reaction, crtBufferIndex, prevBufferIndex, marks, FLAG_PRODUCT);
		if (reaction.check_valid()) {
			vector<Reaction>::const_iterator iter;
			if ((iter = find(reactions.begin(), reactions.end(), reaction)) != reactions.end()) {
				fs.reaction_freq[iter - reactions.begin()].first++;
			}
			else if ((iter = find(reactions.begin(), reactions.end(), -reaction)) != reactions.end()){
				fs.reaction_freq[iter - reactions.begin()].second++;
			}
			else {
				reactions.push_back(reaction);
				fs.reaction_freq.push_back(diint(1, 0));
			}
		}
	}
}


void ReaxReader::Check()
{
	ofstream log("debug.log");
	cerr << "Starting ReaxReader debug session...\n";	log << "Starting ReaxReader debug session...\n";
	cerr << "Checking molecule names...\n";	log << "Checking molecule names...\n";

	for (auto i = 0; i < molecules.size(); i++) {
		for (auto j = 0; j < i; j++) {
			if (molecules[i].toString() == molecules[j].toString()) {
				cerr << "Warning: Molecules " << i << " and " << j << " have the same name!\n";
				log << "Warning: Molecules " << i << " and " << j << " have the same name!\n";
				cerr << "name is: " << molecules[i].toString() << "\n";
				log << "name is: " << molecules[i].toString() << "\n";
			}
		}
	}

	cerr << "Checking molecule numbers...\n"; log << "Checking molecule numbers...\n";

	for (tsize_t i = 1; i < fss.size(); i++) {
		Array dc(molecules.size(), 0);
		for (tsize_t j = 0; j < fss[i].reaction_freq.size(); j++) {
			for (const auto& s : reactions[j].reagants)dc[s] -= (fss[i].reaction_freq[j].first - fss[i].reaction_freq[j].second);
			for (const auto& s : reactions[j].products)dc[s] += (fss[i].reaction_freq[j].first - fss[i].reaction_freq[j].second);
		}
		for (auto k = 0; k < fss[i - 1].mol_freq.size(); k++) {
			if (dc[k] != fss[i].mol_freq[k] - fss[i - 1].mol_freq[k]) {
				cerr << "Warning: At frame " << i << ":\n";
				log << "Warning: At frame " << i << ":\n";
				cerr << "Number of molecule " << k << " cannot match. (" << dc[k]
					<< "," << fss[i - 1].mol_freq[k] << "," << fss[i].mol_freq[k] << ")\n";
				log << "Number of molecule " << k << " cannot match. (" << dc[k]
					<< "," << fss[i - 1].mol_freq[k] << "," << fss[i].mol_freq[k] << ")\n";
				cerr << "name is: " << molecules[k].toString() << "\n";
				log << "name is: " << molecules[k].toString() << "\n";
			}
		}
	}
	log.close();
}

//----------------bottom--------------------------------

void ReaxReader::scan_score(int index, const Matrix& bondMatrix, const Arrayd& atomweight, Mark marks,
	Array& molRoots, Array& scores, Array& molofAtoms) {

	marks[index] = true;
	molofAtoms[index] = molRoots.size() - 1;
	int weight = (int)(atomweight[index] + 0.1);
	int nhcon = 0;
	for (const auto& child : bondMatrix[index]) {
		if (atomweight[child] >= 2 || bondMatrix[child].size() > 1)nhcon++;//not terminal H
		if (!marks[child])scan_score(child, bondMatrix, atomweight, marks, molRoots, scores, molofAtoms);
	}
	scores[index] = totmpscore(weight, nhcon, bondMatrix[index].size() - nhcon);
	if (scores[index] > scores[molRoots.back()])molRoots.back() = index;
}
void ReaxReader::scan_molecule(int index, const Matrix& bondMatrix, Mark marks, MatMolecule& molecule) {
	marks[index] = true;
	molecule.push_back(index);
	for (const auto& child : bondMatrix[index]) {
		if (!marks[child])scan_molecule(child, bondMatrix, marks, molecule);
	}
}

void ReaxReader::create_reaction(int index, Reaction& reac, tsize_t bufferIndex_1, tsize_t bufferIndex_2, vector<bool>* marks, unsigned char flag) {
	//flag = 0/1, means reagant/product
	//it is a quick step, so don't take it too seriously.
	//DFS
	marks[flag][index] = true;
	if (flag == FLAG_PRODUCT)reac.products.push_back(_mol_index_buffer[bufferIndex_1][index]);
	else reac.reagants.push_back(_mol_index_buffer[bufferIndex_1][index]);

	for (const auto& atom : _mol_buffer[bufferIndex_1][index].atoms) {
		//for each atom, find the molecule who contains it, and find other atoms together with it
		int index_2 = _mol_of_atom[bufferIndex_2][atom];		//here atom->data is the absolute index of atom
															//index_2 is the molecule relative index in the new buffer
		if (!marks[!flag][index_2])create_reaction(index_2, reac, bufferIndex_2, bufferIndex_1, marks, !flag);
	}
	return;
}


