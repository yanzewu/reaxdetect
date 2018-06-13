#pragma once

#ifndef REACTIONSTAT_H
#define REACTIONSTAT_H


#include "typedef.h"
#include "trajectory.h"
#include "smiles/smiles.h"

#include <deque>
#include <vector>
#include <list>

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
        unsigned verbose_interval;
		Config() : buffer_size(2), recognize_begin(0), verbose_interval(100) {
		}
	};

	//Frequecy of molecules and reactions in a frame
	struct FrameStat {
		double t;	// time stamp
		Array mol_freq;
		vector<diint> reaction_freq;
	};

    // Inner representation of reaction, using molecule ID (int) of reactants and products.
	struct Reaction {
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
		Reaction operator-()const; // reverse reaction
		string to_string(const vector<string>&)const;
	};

	vector<FrameStat> fss;			// reaction and species frequency in each frame
	vector<smiles> molecules;		// global species list
	vector<string> species;			// global species name list
	vector<Reaction> reactions;		// global reaction list

	ReaxReader() : config() { }
	ReaxReader(const Config& c) : config(c) {
	}

	//handles data from trajectory, including the file reading loop.
	int HandleData(TrajReader&, const Simulation&);
	
	//check (use only in debug)
	void Check();

private:

	// Intermediate representation of molecule. Faster than SMILES, but cannot guarantee uniqueness -- 
    // Same MatMolecule ==> maybe same molecule; Different MatMolecule ==> must be different molecule.
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
		bool has_same_atoms(MatMolecule&);   // If two MatMolecules are Identical --- with exactly same atom IDs
		smiles to_smiles()const;    // generate canonical SMILES representation
	};

    // Working environment of ReaxReader in one frame
	struct BufferPage {
		Array mol_id;	                    // global id of each molecule
		vector<MatMolecule> molecule;		// MatMolecule list
		vector<Reaction> raw_reaction;	    // reaction buffer for raw reactions
		vector<Reaction> reaction;		    // real reaction buffer
		Array mol_of_atom;		            // which molecule the atom belongs to
		Array atom_score;	                // atom score, used by MatMolecule
		Matrix bond_matrix;			        // connectivity adjacent list, used by MatMolecule
	};

	Config config;

	BufferPage* _crt_buffer;                // pointer to current buffer
    std::deque<BufferPage*> _prev_buffer;   // Previous buffer queue

	//detect molecule and write into buffer page.
	void RecognizeMolecule(const TrajFrame& frm, const Arrayd& atomWeights, int atomNumber, FrameStat& fs);

	//detect reaction and remove repeatance.
	void RecognizeReaction(const TrajFrame& frm);
	
	//commit reaction into lists
	void CommitReaction(FrameStat& fs_commit);

    //initiating a new buffer, push the old one into _prev_buffer
    void RenewBuffer(unsigned max_buffer_num, size_t max_atom_idx);

	/*molecule generating subfunctions*/

	//[DFS] Scan raw scores of each atom. [Do not use directly]
	static void scan_score(int index, const Arrayd& atom_weights, BufferPage* buffer, Mark marks, Array& molRoots);
	
	//[DFS] Scan molecule chain. [Do not use directly]
	static void scan_molecule(int index, BufferPage* buffer, Mark marks, MatMolecule& molecule);
	
	//[DFS] Scan reaction between two frames. [Do not use directly]
	static void create_reaction(int index, const BufferPage* page_1, const BufferPage* page_2,
		vector<bool>* marks, Reaction& reac, Reaction& raw_reac, unsigned char flag);

	//checking if two raw reactions contains same molecules
	static bool check_reaction(Reaction & reaction_1, BufferPage* buffer_1_prod, BufferPage* buffer_1_reac,
		Reaction & reaction_2, BufferPage* buffer_2_prod, BufferPage* buffer_2_reac);

	//generate score to raw score
	static int totmpscore(int weight, int nhcon, int bond, int terminalh);

	//make tmp score to smiles score
	static int tmpscore2score(int tmpscore);

	//test if an atom is terminal hydrogen
	static bool isterminalhydrogen(int tmpscore);

};

#endif // !REACTIONSTAT_H
