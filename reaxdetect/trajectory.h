// reading trajetory files.

#pragma once

#ifndef TRAJECTORY_H
#define TRAJECTORY_H

#include <map>
#include <vector>
#include <fstream>

#include "simulation.h"

using namespace std;


// reader of various trajectory file
class TrajReader {
public:

	struct Config {
		map<int, vector<double> > bondorder_cutoff;
		vector<double> bondorder_cutoff_default;

		int read_atompos;
		int count_bondorder;
		Config() : read_atompos(0) {
		}
	};

	struct atom{
		int id;
		double x, y, z, q;
	};

	struct bond{
		int id_1, id_2;
		int order;
		bond() : order(0) {}
		bond(int a, int b) : id_1(a), id_2(b), order(1) {
		}
	};

	//Data Structure recording information of atoms and bonds
	struct Frame {
		vector<atom> atoms;
		vector<bond> bonds;
	};	

	vector<int> atomTypes;		//Atom Types List
	vector<double> atomWeights;	//Atom Weights List
	
	map<int, vector<double> > bondorders; // Bond order statistics

	TrajReader() : config() {
	}
	TrajReader(const Config& c) : config(c){
	}

	// open file and read head
	int Open(const string& path);

	// read head
	int ReadTrjHead(Simulation* simulation);

	// read frame. no use directly
	int ReadTrjFrame(Frame&);

private:
	Config config;
	ifstream trjfile;
	int frameCount;
};
#endif // !TRAJECTORY_H
