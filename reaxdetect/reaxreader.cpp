//contains algorithm about molecule and Reaction generate

#include <functional>
#include <iostream>
#include <string.h>

#include "reaxreader.h"
#include "smiles/smiles.h"
#include "util/algorithmutil.h"
#include "util/path.h"

#define FLAG_REAGANT 0
#define FLAG_PRODUCT 1

//---------------top api----------------------

int ReaxReader::HandleData(TrajReader& reader, const Simulation& simulation)
{

    /* Main logic of ReaxDetect */

	printf("Buffer size=%zd\n", config.buffer_size);
	printf("Read interval=%d\n", config.recognize_interval);
	printf("Read begin=%d\n", config.recognize_begin);
	printf("Read end=%d\n", config.recognize_limit);

    TrajFrame _frame;
	size_t i = 0;           // how many frames readed
	size_t frame_i = 0;     // how many frames analyzed
    _crt_buffer = nullptr;
    _prev_buffer.clear();

    // Main loop
	while (reader.ReadTrjFrame(_frame, simulation)) {   // read a frame from simulation trajectory

		if (i % config.recognize_interval == 0 && i >= config.recognize_begin) {    // if should parse
			FrameStat fstat;
			fstat.t = simulation.timeStep * i;  // time

            RenewBuffer(config.buffer_size, simulation.atomNumber + 1);
			RecognizeMolecule(_frame, simulation.atomWeights, simulation.atomNumber, fstat);    // analyze connectivity; write molecule info, connectivity to buffer
			RecognizeReaction(_frame);  // detect reactions; write reaction info to buffer; check reverse reactions;
			fss.push_back(fstat);
			if (frame_i >= config.buffer_size - 1) {
				CommitReaction(fss[frame_i - config.buffer_size + 1]); // the reverse reactions have been checked, write the final reaction frequency in.
			}
			frame_i++;
		}
		i++;
        if (i % config.verbose_interval == 0) {
            printf("Frame read: %zd\n", i);
        }
		if (i == config.recognize_limit) {
			break;
		}
        _frame.clear();
	}
    
    RenewBuffer(config.buffer_size, simulation.atomNumber + 1);
	for (size_t j = frame_i - config.buffer_size + 1; j < frame_i; j++) {
		CommitReaction(fss[j]);
	}
	printf("Total %zd frames read, with %zd molecules and %zd reactions.\n", fss.size(), molecules.size(), reactions.size());

    if (_crt_buffer) {
        delete _crt_buffer;
    }

	printf("Encoding smiles...\n");
	for (auto& mol : molecules) {
		species.push_back(mol.to_smiles());
	}

	//alignment
	for (auto& fs : fss) {
		fs.mol_freq.resize(molecules.size());
		fs.reaction_freq.resize(reactions.size());
	}

#ifdef _DEBUG
	Check();
#endif // _DEBUG
	return 0;
}

//----------------middle api------------------

void ReaxReader::RenewBuffer(unsigned max_buffer_num, size_t max_atom_idx) {

    if (_prev_buffer.size() >= config.buffer_size - 1) {
        delete _prev_buffer.back();
        _prev_buffer.pop_back();
    }
    if (_crt_buffer) {
        _prev_buffer.push_front(_crt_buffer);
    }

    _crt_buffer = new BufferPage;
    _crt_buffer->mol_of_atom.assign(max_atom_idx, -1);
    _crt_buffer->bond_matrix.resize(max_atom_idx);
    _crt_buffer->atom_score.resize(max_atom_idx);
}

void ReaxReader::RecognizeMolecule(const TrajFrame& frm, const Arrayd& atomWeights, int atomNumber, FrameStat& fs)
{
    // 1. build bonding matrix
	for (const auto& bnd : frm.bonds) {
		_crt_buffer->bond_matrix[bnd.id_1].push_back({ bnd.id_2, bnd.order });
		_crt_buffer->bond_matrix[bnd.id_2].push_back({ bnd.id_1, bnd.order });
	}

	// 2. first scan(DFS): Scan atom groups. Assign to molRoots, mol_of_atom
	Array molRoots;//root of molecules
	Mark marks = new bool[atomNumber + 1]{ false };

	for (auto i = 1; i <= atomNumber; i++) {
		if (marks[i])continue;
		molRoots.push_back(i);
		scan_score(i, atomWeights, _crt_buffer, marks, molRoots);
	}

	// 3. sort bondmatrix: O(N*k^2) for MatMolecule comparision;
	for (auto bonds = _crt_buffer->bond_matrix.begin() + 1; bonds != _crt_buffer->bond_matrix.end(); bonds++) {
		sortby_pair(*bonds, _crt_buffer->atom_score, greater<int>());
	}
	memset(marks, 0, sizeof(bool)*(atomNumber + 1));

	// 4. second scan(DFS): assign MatMolecule (intermediate representation) list
	for (const auto& root : molRoots) {
		MatMolecule molecule(&_crt_buffer->bond_matrix, &_crt_buffer->atom_score);
		scan_molecule(root, _crt_buffer, marks, molecule);
		_crt_buffer->molecule.push_back(move(molecule));
	}
	delete[] marks;
	
	// 5. count MatMolecule frequency O(N*m). This step is to reduce the frequency of SMILES generation (only performed at unique MatMolecule) rather than all MatMolecules.
	vector<tsize_t> uroots;     // Unique Matmolecule index list. Other molecules in buffer->molecule must be same as one of them.
	Array uFreq;		        // Frequency of each unique MatMolecule in frame.
	_crt_buffer->mol_id.resize(molRoots.size());

	for (tsize_t i = 0; i < _crt_buffer->molecule.size(); i++) {
		_crt_buffer->mol_id[i] = (int)pushnew2index(i, _crt_buffer->molecule, uroots, uFreq);  // here mol_id stores frequency of MatMolecule
	}

	// 6. generating smiles: O(m^2*log(m)). Only performed on each MatMolecule.
	Array uindex(uroots.size());
	fs.mol_freq.resize(molecules.size());
	for (tsize_t i = 0; i < uroots.size(); i++) {
		smiles newsms = _crt_buffer->molecule[uroots[i]].to_smiles();
			
		auto findpos = find(molecules.begin(), molecules.end(), newsms);
		if(findpos != molecules.end()){
			fs.mol_freq[findpos - molecules.begin()] += uFreq[i];
			uindex[i] = findpos - molecules.begin();
		}
		else {
			molecules.push_back(move(newsms));
			fs.mol_freq.push_back(uFreq[i]);
			uindex[i] = (int)molecules.size() - 1;
		}
	}

	// 7. Replacing MatMolecule frequency to SMILES frequency.
	for (auto& molindex : _crt_buffer->mol_id) {
		molindex = uindex[molindex];
	}
}
void ReaxReader::RecognizeReaction(const TrajFrame& frm)
{
	if (_prev_buffer.empty())return;

	vector<bool> marks[2]{ 
		vector<bool>(_prev_buffer.front()->molecule.size(), false),
		vector<bool>(_crt_buffer->molecule.size(), false) 
	}; //marks of product and reagant

	for (tsize_t i = 0; i < _crt_buffer->molecule.size(); i++) {
		if (marks[1][i])continue;
		Reaction raw_reaction, reaction;
		create_reaction(i, _crt_buffer, _prev_buffer.front(), marks, reaction, raw_reaction, FLAG_PRODUCT);

		//stage reaction
		if (reaction.check_valid()) {
			Reaction raw_reaction_rev = -raw_reaction;
			Reaction reaction_rev = -reaction;

			auto next_buffer = _prev_buffer.begin();
			auto scan_buffer = next_buffer++;

			// look up previous buffer to remove reverse reactions
			for (; next_buffer != _prev_buffer.end(); scan_buffer = next_buffer++) {
				auto reaction_scanned = (*scan_buffer)->reaction.begin();
				auto raw_reaction_scanned = (*scan_buffer)->raw_reaction.begin();

				for (; reaction_scanned != (*scan_buffer)->reaction.end(); reaction_scanned++, raw_reaction_scanned++)
				{
					if (reaction_rev != *reaction_scanned)continue; // two reactions does not share same formula ==> must be different

                    // share same formula ==> chceck if atoms are identical.
					if (check_reaction(raw_reaction_rev, _prev_buffer.front(), _crt_buffer, *raw_reaction_scanned, *scan_buffer, *next_buffer)) {
						(*scan_buffer)->raw_reaction.erase(raw_reaction_scanned);
						(*scan_buffer)->reaction.erase(reaction_scanned);
						return;
					}
				}
			}
			_crt_buffer->raw_reaction.push_back(raw_reaction);
			_crt_buffer->reaction.push_back(reaction);
		}
	}
}
void ReaxReader::CommitReaction(FrameStat & fs_commit)
{
	fs_commit.reaction_freq.resize(reactions.size());

	// commit reaction into global list
	for(const auto& reaction : _prev_buffer.back()->reaction){
		vector<Reaction>::const_iterator iter;
		if ((iter = find(reactions.begin(), reactions.end(), reaction)) != reactions.end()) {
			fs_commit.reaction_freq[iter - reactions.begin()].first++;
		}
		else if ((iter = find(reactions.begin(), reactions.end(), -reaction)) != reactions.end()){
			fs_commit.reaction_freq[iter - reactions.begin()].second++;
		}
		else {
			reactions.push_back(move(reaction));
			fs_commit.reaction_freq.push_back(diint(1, 0));
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
			if (molecules[i].to_string() == molecules[j].to_string()) {
				cerr << "Warning: Molecules " << i << " and " << j << " have the same name!\n";
				log << "Warning: Molecules " << i << " and " << j << " have the same name!\n";
				cerr << "name is: " << molecules[i].to_string() << "\n";
				log << "name is: " << molecules[i].to_string() << "\n";
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
				cerr << "name is: " << molecules[k].to_string() << "\n";
				log << "name is: " << molecules[k].to_string() << "\n";
			}
		}
	}
	log.close();
	cerr << "Check finished." << endl;
}

//----------------bottom--------------------------------

void ReaxReader::scan_score(int index, const Arrayd& atomweight, BufferPage* buffer, Mark marks, Array& mol_roots) {

	marks[index] = true;
	buffer->mol_of_atom[index] = (int)mol_roots.size() - 1;
	int weight = (int)(atomweight[index] + 0.1);
	int nhcon = 0;
	int bondcon = 0;
	for (const auto& child : buffer->bond_matrix[index]) {
		bondcon += child.second;
		if (child.second == 0 || atomweight[child.first] >= 2 || buffer->bond_matrix[child.first].size() > 1)nhcon++;//not terminal H
		if (!marks[child.first])scan_score(child.first, atomweight, buffer, marks, mol_roots);
	}
	buffer->atom_score[index] = totmpscore(weight, nhcon, bondcon, (int)buffer->bond_matrix[index].size() - nhcon);
	if (buffer->atom_score[index] > buffer->atom_score[mol_roots.back()])mol_roots.back() = index;
}

void ReaxReader::scan_molecule(int index, BufferPage* buffer, Mark marks, MatMolecule& molecule) {
	marks[index] = true;
	molecule.push_back(index);
	for (const auto& child : buffer->bond_matrix[index]) {
		if (!marks[child.first])scan_molecule(child.first, buffer, marks, molecule);
	}
}

void ReaxReader::create_reaction(int index, const BufferPage* buffer_host, const BufferPage* buffer_guest, vector<bool>* marks, Reaction& reac, Reaction& raw_reac, unsigned char flag) {
	//flag = 0/1, means reagant/product
	//it is a quick step, so don't take it too seriously.
	//DFS
	marks[flag][index] = true;
	if (flag == FLAG_PRODUCT) {
		reac.products.push_back(buffer_host->mol_id[index]);
		raw_reac.products.push_back(index);
	}
	else {
		reac.reagants.push_back(buffer_host->mol_id[index]);
		raw_reac.reagants.push_back(index);
	}

	for (const auto& atom : buffer_host->molecule[index].atoms) {
		//for each atom, find the molecule who contains it, and find other atoms together with it
		int index_2 = buffer_guest->mol_of_atom[atom];		//here atom->data is the absolute index of atom
															//index_2 is the molecule relative index in the new buffer
		if (!marks[!flag][index_2])create_reaction(index_2, buffer_guest, buffer_host, marks, reac, raw_reac, !flag);
	}
	return;
}

bool ReaxReader::check_reaction(Reaction & reaction_1, BufferPage* buffer_1_prod, BufferPage* buffer_1_reac,
	Reaction & reaction_2, BufferPage* buffer_2_prod, BufferPage* buffer_2_reac)
{
	sortby(reaction_1.reagants, buffer_1_reac->mol_id, less<int>());
	sortby(reaction_2.reagants, buffer_2_reac->mol_id, less<int>());

	auto mol_1 = reaction_1.reagants.begin(), mol_2 = reaction_2.reagants.begin();
	for (; mol_1 != reaction_1.reagants.end(); mol_1++, mol_2++) {
		if (!buffer_1_reac->molecule[*mol_1].has_same_atoms(buffer_2_reac->molecule[*mol_2])) {  // if *EXACTLY* identical ==> remove.
			return false;
		}
	}

	sortby(reaction_1.products, buffer_1_prod->mol_id, less<int>());
	sortby(reaction_2.products, buffer_2_prod->mol_id, less<int>());

	mol_1 = reaction_1.products.begin(); 
	mol_2 = reaction_2.products.begin();
	for (; mol_1 != reaction_1.products.end(); mol_1++, mol_2++) {
		if (!buffer_1_prod->molecule[*mol_1].has_same_atoms(buffer_2_prod->molecule[*mol_2])) {
			return false;
		}
	}
	return true;
}
