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
	printf("Buffer size=%zd\n", config.buffer_size);
	printf("Read interval=%d\n", config.recognize_interval);
	printf("Read begin=%d\n", config.recognize_begin);
	printf("Read end=%d\n", config.recognize_limit);

	_buffer_pages = new BufferPage[config.buffer_size];

	InitBuffer(simulation.atomNumber + 1);
	TrajReader::Frame _frame;
	size_t i = 0;
	size_t frame_i = 0;
	while (reader.ReadTrjFrame(_frame)) {
		if (i % config.recognize_interval == 0 && i >= config.recognize_begin) {
			FrameStat fstat;
			fstat.t = simulation.timeStep * i;
			RecognizeMolecule(_frame, reader.atomWeights, simulation.atomNumber, fstat);

			// fast return: First frame
			if (frame_i == 0) {
				SwapBuffer();
				fss.push_back(fstat);
				frame_i++;
				break;
			}
			RecognizeReaction(_frame);
			i++;

			// mark last buffer as clean; At this time prevbuffer should have only one item
			_prev_buffer.front()->dirty = false;
			_prev_buffer.clear();

			for (size_t j = 0; j < config.buffer_size - 2; j++) {
				for (int k = 0; k < config.buffer_interval; k++) {
					if (!reader.ReadTrjFrame(_frame))goto recend;
					i++;
				}

				SwapBuffer();

				FrameStat fstat_tmp;
				RecognizeMolecule(_frame, reader.atomWeights, simulation.atomNumber, fstat_tmp);
				RecognizeReaction(_frame);
			}

			CommitReaction(fstat);
			fss.push_back(fstat);

			// mark those bs-2 buffers as clean
			while (_prev_buffer.size() > 1) {
				_prev_buffer.front()->dirty = false;
				_prev_buffer.pop_front();
			}

			// load a new buffer; push last into prevbuffer
			SwapBuffer();
			// however last is useless
			_prev_buffer.front()->dirty = false;
			_prev_buffer.pop_front();
		}

		if (i == config.recognize_limit) {
			break;
		}
	}

recend:

	printf("Total %zd frames read, with %zd molecules and %zd reactions.\n", fss.size(), molecules.size(), reactions.size());

	printf("Encoding smiles...\n");
	for (auto& mol : molecules) {
		species.push_back(mol.to_smiles());
	}

	//alignment
	for (auto& fs : fss) {
		fs.mol_freq.resize(molecules.size());
		fs.reaction_freq.resize(reactions.size());
	}

	delete[] _buffer_pages;

#ifdef _DEBUG
	Check();
#endif // _DEBUG
	return 0;
}

//----------------middle api------------------

void ReaxReader::InitBuffer(size_t max_atom_idx)
{
	_prev_buffer.clear();
	_crt_buffer = &_buffer_pages[0];
	for (BufferPage* page = _buffer_pages; page != _buffer_pages + config.buffer_size; page++) {
		page->mol_of_atom.assign(max_atom_idx, -1);
		page->bond_matrix.resize(max_atom_idx);
		page->atom_score.resize(max_atom_idx);
	}
}

void ReaxReader::SwapBuffer() {
	if (_prev_buffer.size() >= config.buffer_size - 1) {
		_prev_buffer.pop_back();
	}
	_prev_buffer.push_front(_crt_buffer);
	
	auto b = _buffer_pages;
	for (; b < _buffer_pages + config.buffer_size; b++) {
		if (!b->dirty) {
			_crt_buffer = b;
			break;
		}
	}
	if (b == _buffer_pages + config.buffer_size) {
		throw runtime_error("No clean buffer");
	}

	_crt_buffer->molecule.clear();
	_crt_buffer->mol_idx.clear();
	_crt_buffer->raw_reaction.clear();
	_crt_buffer->reaction.clear();
	_crt_buffer->dirty = true;

	for (auto& m : _crt_buffer->mol_of_atom)
		m = -1;
	for (auto& a : _crt_buffer->atom_score)
		a = 0;
	for (auto& b : _crt_buffer->bond_matrix)
		b.clear();
}
void ReaxReader::RecognizeMolecule(const TrajReader::Frame& frm, const Arrayd& atomWeights, int atomNumber, FrameStat& fs)
{
	for (const auto& bnd : frm.bonds) {
		_crt_buffer->bond_matrix[bnd.id_1].push_back({ bnd.id_2, bnd.order });
		_crt_buffer->bond_matrix[bnd.id_2].push_back({ bnd.id_1, bnd.order });
	}
	//1.first scan(DFS): get score, root, and mol of atoms: O(N)
	Array molRoots;//root of molecules
	Mark marks = new bool[atomNumber + 1]{ false };

	for (auto i = 1; i <= atomNumber; i++) {
		if (marks[i])continue;
		molRoots.push_back(i);
		scan_score(i, atomWeights, _crt_buffer, marks, molRoots);
	}
	//2.sort bondmatrix: O(N*k^2)
	for (auto bonds = _crt_buffer->bond_matrix.begin() + 1; bonds != _crt_buffer->bond_matrix.end(); bonds++) {
		sortby_pair(*bonds, _crt_buffer->atom_score, greater<int>());
	}
	memset(marks, 0, sizeof(bool)*(atomNumber + 1));

	//3.second scan(DFS): get molecule list: O(N)
	for (const auto& root : molRoots) {
		MatMolecule molecule(&_crt_buffer->bond_matrix, &_crt_buffer->atom_score);
		scan_molecule(root, _crt_buffer, marks, molecule);
		_crt_buffer->molecule.push_back(move(molecule));
	}
	delete[] marks;
	
	//4.compress trees: O(N*m)
	vector<tsize_t> uroots;		//index in molbuffer
	Array uFreq;		//freq corresponding to molecules.
	_crt_buffer->mol_idx.resize(molRoots.size());

	for (tsize_t i = 0; i < _crt_buffer->molecule.size(); i++) {
		_crt_buffer->mol_idx[i] = (int)pushnew2index(i, _crt_buffer->molecule, uroots, uFreq);
	}

	//5.final: get smiles: O(m^2*log(m))
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
	//6.restore uindex to index: O(N)
	for (auto& molindex : _crt_buffer->mol_idx) {
		molindex = uindex[molindex];
	}
}
void ReaxReader::RecognizeReaction(const TrajReader::Frame& frm)
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

			// look up previous buffer for quick-reverse reactions
			for (; next_buffer != _prev_buffer.end(); scan_buffer = next_buffer++) {
				auto reaction_scanned = (*scan_buffer)->reaction.begin();
				auto raw_reaction_scanned = (*scan_buffer)->raw_reaction.begin();

				for (; reaction_scanned != (*scan_buffer)->reaction.end(); reaction_scanned++, raw_reaction_scanned++)
				{
					if (reaction_rev != *reaction_scanned)continue;
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

void ReaxReader::create_reaction(int index, BufferPage* buffer_host, BufferPage* buffer_guest, vector<bool>* marks, Reaction& reac, Reaction& raw_reac, unsigned char flag) {
	//flag = 0/1, means reagant/product
	//it is a quick step, so don't take it too seriously.
	//DFS
	marks[flag][index] = true;
	if (flag == FLAG_PRODUCT) {
		reac.products.push_back(buffer_host->mol_idx[index]);
		raw_reac.products.push_back(index);
	}
	else {
		reac.reagants.push_back(buffer_host->mol_idx[index]);
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
	sortby(reaction_1.reagants, buffer_1_reac->mol_idx, less<int>());
	sortby(reaction_2.reagants, buffer_2_reac->mol_idx, less<int>());

	auto mol_1 = reaction_1.reagants.begin(), mol_2 = reaction_2.reagants.begin();
	for (; mol_1 != reaction_1.reagants.end(); mol_1++, mol_2++) {
		if (!buffer_1_reac->molecule[*mol_1].equals_to(buffer_2_reac->molecule[*mol_2])) {
			return false;
		}
	}

	sortby(reaction_1.products, buffer_1_prod->mol_idx, less<int>());
	sortby(reaction_2.products, buffer_2_prod->mol_idx, less<int>());

	mol_1 = reaction_1.products.begin(); 
	mol_2 = reaction_2.products.begin();
	for (; mol_1 != reaction_1.products.end(); mol_1++, mol_2++) {
		if (!buffer_1_prod->molecule[*mol_1].equals_to(buffer_2_prod->molecule[*mol_2])) {
			return false;
		}
	}
	return true;
}
