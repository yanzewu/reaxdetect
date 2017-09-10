#pragma once

#ifndef REACTIONSTAT_H
#define REACTIONSTAT_H

#include <vector>
#include <list>

#include "typedef.h"
#include "trajectory.h"
#include "smiles/smiles.h"

using namespace std;

#define MAX_ATOM_TYPE	256

/* Append in 9/5/2016. Version 2.1.3 */

//recording indices of reagants and products in a reaction
//basic statistics about the reaction and molecules
class ReaxReader {
public:

	typedef vector<list<pair<int, int> > > Matrix;
	typedef bool* Mark;

	struct Config {
		size_t buffer_size;
		int recognize_interval;
		int recognize_limit;
		int recognize_begin;
		Config() : buffer_size(2), recognize_begin(0) {
		}
	};

	//Frequecy of molecules and reactions in a frame
	struct FrameStat {
		double t;	// time stamp
		Array mol_freq;
		vector<diint> reaction_freq;
	};

	class Reaction {
	public:
		Array reagants;
		Array products;
		Reaction() {}
		Reaction(const Array& reagants, const Array& products) :
			reagants(reagants), products(products)
		{
		}
		bool check_valid();
		bool operator==(const Reaction&)const;
		bool operator!=(const Reaction&)const;
		Reaction operator-()const;
		string to_string(const vector<string>&)const;
	};

	vector<FrameStat> fss;			//Raw data. Do not use directly.
	vector<smiles> molecules;		//Smiles of all molecules in the trajectory. Do not use directly.
	vector<string> species;			// name of species
	vector<Reaction> reactions;		//List of all different reactions

	ReaxReader() : config() { }
	ReaxReader(const Config& c) : config(c) {
	}

	//handles data from trajectory
	int HandleData(TrajReader&, const Simulation&);//do before output
	
	//check (use only in debug)
	void Check();

protected:

	//Data structure storing a raw molecule
	struct MatMolecule {
		list<int> atoms;
		Matrix* pbondmatrix;
		Array* pscores;

		MatMolecule(Matrix* const pbm, Array* const pscr) :
		    pbondmatrix(pbm), pscores(pscr)
		{
		}
		void push_back(int index);
		bool operator==(const MatMolecule&)const;
		bool equals_to(MatMolecule&);
		smiles to_smiles()const;
	};


	struct BufferPage {
		Array mol_idx;	//map of relative index to absolute index	
		vector<MatMolecule> molecule;		//molecule buffer
		vector<Reaction> raw_reaction;	// reaction buffer for raw reactions
		vector<Reaction> reaction;		// real reaction buffer
		Array mol_of_atom;		//map of atom to the molecule who contains it	
		Array atom_score;	//raw data from the bond in one frame	
		Matrix bond_matrix;			//atom score buffer in one frame
	};

	Config config;
	BufferPage* _buffer_pages;

	BufferPage* _crt_buffer;
	list<BufferPage*> _prev_buffer;	//buffer index

	//detect molecule and write crtmol, crtfmi, fss.molfreq
	void RecognizeMolecule(const TrajReader::Frame& frm, const Arrayd& atomWeights, int atomNumber, FrameStat& fs);

	//detect reaction and write fss.reactions
	void RecognizeReaction(const TrajReader::Frame& frm);
	
	//commit reaction into lists
	void CommitReaction(FrameStat& fs_commit);

	//count bond order for different types
	void CountBondOrder(const TrajReader::Frame& frm, const Arrayd& atomWeights);

	//initialize buffer pages
	void InitBuffer(size_t max_atom_idx);

	//swap current buffer and previous buffer. Also do initialize work at the same time.
	void SwapBuffer();

	/*molecule generating subfunctions*/

	//[DFS] Scan raw scores of each atom. [Do not use directly]
	void scan_score(int index, const Arrayd& atom_weights, BufferPage* buffer, Mark marks, Array& molRoots);
	
	//[DFS] Scan molecule chain. [Do not use directly]
	void scan_molecule(int index, BufferPage* buffer, Mark marks, MatMolecule& molecule);
	
	//[DFS] Scan reaction between two frames. [Do not use directly]
	void create_reaction(int index, BufferPage* page_1, BufferPage* page_2,
		vector<bool>* marks, Reaction& reac, Reaction& raw_reac, unsigned char flag);

	//checking if two raw reactions contains same molecules
	bool check_reaction(Reaction & reaction_1, BufferPage* buffer_1_prod, BufferPage* buffer_1_reac,
		Reaction & reaction_2, BufferPage* buffer_2_prod, BufferPage* buffer_2_reac);

	//generate score to raw score
	static int totmpscore(int weight, int nhcon, int bond, int terminalh);

	//make tmp score to smiles score
	static int tmpscore2score(int tmpscore);

	//test if an atom is terminal hydrogen
	static bool isterminalhydrogen(int tmpscore);

};

#endif // !REACTIONSTAT_H
