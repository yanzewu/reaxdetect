#pragma once
#include"includes.h"
#ifndef SMILES_H
#define SMILES_H

/* Written in 9/5/2016. Version 1.1. */

#define SMILESVERSION 101

#define maxhnum		0x08
#define maxnhcon	0x08
#define maxatmwht	0x80

//A uniquely representation of a molecule, canonical smiles
class smiles {
public:
	smiles() {}
	bool operator==(const smiles&)const;

	//kernel convert function
	void localconvert();
	
	//Write the string form in symbol
	void compile();

	/* interface for generating functions */

	//push a atom back with smiles score
	void push_back(int score) {
		atoms.push_back(snode(score));
	}

	//connect a parent with child
	void push_child(int parent, int index) {
		atoms[parent].children.push_back(index);
	}

	//return the number of atom with certain atom weight
	int has_element(int weight)const;

	//number of atoms
	size_t size() { return atoms.size(); }

	//output smiles. Just for compability.
	string toString()const { return symbol; }

	//output smiles. Valid only after compile();
	string toSmiles()const { return symbol; }

	//read smiles string. This function may contain error.
	void readSmiles(const string& sms);

	//written binary form. Do not use directly.
	unsigned int toBin(int* buffer)const;

	//read binary form. Do not use directly.
	void readBin(int* buffer);

protected:
	struct snode {
		int score;
		list<int> children;
		snode() { score = 0; }
		snode(int _score) { score = _score; }
		void swap(snode& other) {
			std::swap(score, other.score);
			children.swap(other.children);
		}
	};							//children tree node.
	vector<snode> atoms;		//atom list
	vector<int> looproots;		//roots of loop. See documentation for details.
	string symbol;				//string of chain.

	/* functions used in converting */

	//rank atoms by its score, and write into rank
	bool rankby(vector<int>& score, vector<int>& rank);

	//update score by scores of neighbours
	bool renew_score(const int index, const vector<int>& rank, const list<int>& _children_tmp, vector<int>& score);

	//label each atom
	void label(int index, int& crtlabel, vector<list<int> >& bonds, const vector<int>& rank, vector<int>& score);

	//recursive compare function. Do not use directly.
	bool compare(const smiles& sms, int myindex, int hisindex, bool* mark)const;

	/* functions used in compiling */

	//scans the loop of smiles, and write into looproots
	void scanloop(int index, int parent, bool* mark);

	//recursive compile function. Do not use directly.
	void compilepart(int index, const vector<string>& loopmark);

	//turn score to weight
	static int score2weight(int score);

	//turn score to hydrogen numbers
	static int score2hydro(int score);
};


#endif // !SMILES_H
