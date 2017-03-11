#pragma once

#ifndef REACTIONSTAT_H
#define REACTIONSTAT_H

#include "includes.h"
#include "typedef.h"
#include "trajectory.h"
#include "smiles.h"

/* Append in 9/5/2016. Version 2.1.3 */

//recording indices of reagants and products in a reaction
//basic statistics about the reaction and molecules
class ReaxReader {
public:

	typedef vector<list<int> > Matrix;
	typedef bool* Mark;

	//Frequecy of molecules and reactions in a frame
	struct FrameStat {
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
		bool operator==(const Reaction& reac)const;
		Reaction operator-()const;
		string to_string(const vector<string>&)const;
	};

	vector<FrameStat> fss;			//Raw data. Do not use directly.
	vector<smiles> molecules;		//Smiles of all molecules in the trajectory. Do not use directly.
	vector<string> species;			// name of species
	vector<Reaction> reactions;		//List of all different reactions

	ReaxReader() { }

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
		smiles to_smiles()const;
	};

	Array _mol_index_buffer[2];	//map of relative index to absolute index	
	vector<MatMolecule> _mol_buffer[2];		//molecule buffer	
	Array _mol_of_atom[2];		//map of atom to the molecule who contains it	
	Matrix _bond_matrix;	//raw data from the bond in one frame	
	Array _atom_score;			//atom score buffer in one frame
	TrajReader::Frame _frame_buffer;	//raw frame data

	tsize_t _crt_buffer_idx, _prev_buffer_idx;	//buffer index

	//detect molecule and write crtmol, crtfmi, fss.molfreq
	void RecognizeMolecule(const TrajReader::Frame& frm, FrameStat& fs, tsize_t crtBufferIndex, int atomNumber, const vector<double>& atomWeights);

	//detect reaction and write fss.reactions, calculate inflow
	void RecognizeReaction(const TrajReader::Frame& frm, FrameStat& fs, tsize_t crtBufferIndex, tsize_t prevBufferIndex);
	
	//swap current buffer and previous buffer. Also do initialize work at the same time.
	void SwapBuffer();

	/*molecule generating subfunctions*/

	//[DFS] Scan raw scores of each atom. [Do not use directly]
	void scan_score(int index, const Matrix& bondMatrix, const Arrayd& atomweight, Mark marks,
		Array& molRoots, Array& scores, Array& molofAtoms);
	
	//[DFS] Scan molecule chain. [Do not use directly]
	void scan_molecule(int index, const Matrix& bondMatrix, Mark marks, MatMolecule& molecule);
	
	//[DFS] Scan reaction between two frames. [Do not use directly]
	void create_reaction(int index, Reaction& reac, tsize_t bufferIndex_1, tsize_t bufferIndex_2,
		vector<bool>* marks, unsigned char flag);

	//generate score to raw score
	static int totmpscore(int weight, int nhcon, int terminalh);

	//make tmp score to smiles score
	static int tmpscore2score(int tmpscore);

	//test if an atom is terminal hydrogen
	static bool isterminalhydrogen(int tmpscore);

};

#endif // !REACTIONSTAT_H
