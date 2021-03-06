// reading trajetory files.

#pragma once

#ifndef TRAJECTORY_H
#define TRAJECTORY_H

#include "simulation.h"

#include <map>
#include <vector>
#include <fstream>

// Stores atom position and charge.
struct atom {
    int id;
    double x, y, z, q;
};

// Stores bonded pair info.
struct bond {
    int id_1, id_2;
    int order;  /* notice: order 0 is a valid option! 
                If you need output order other than 0-3, please modify ``bond_symbol`` in smiles.cpp*/
    bond() : order(0) {}
    bond(int a, int b) : id_1(a), id_2(b), order(1) {
    }
};

// Stores all information in one frame.
struct TrajFrame {
    std::vector<atom> atoms;
    std::vector<bond> bonds;

    void clear() {
        atoms.clear();
        bonds.clear();
    }
};

// Interface of trajectory reader
// You may write your own trajectory reader, as long as the functions are implemented.
class TrajReader {
public:

    // Open the trajectory file. Return non-zero to halt.
    virtual int Open(const std::string&) = 0;

    /* Executed before read first frame. Return non-zero to halt.
    The following entries in Simulation class must be filled:
        - atomNumber: number of atom in the system;
        - atomWeights: list of atom weights of specified index.
    You may also modify ``timeStep`` and ``volume`` to overwritten command-line input.

     For trajectory that has no head, you can read first frame to get information and store it in
     a buffer, then output other information when ``ReadTrjFrame`` is first called.
    */
    virtual int ReadTrjHead(Simulation*) = 0;

    /* Executed in the read loop. Return non-zero to halt.
    If the system condition (atomWeights, timeStep, volume) is dynamically changed, you may remove the ``const`` before
    Simulation&. Some other places may also need to be changed. But the code should work as long as it can be compiled.
    */
    virtual int ReadTrjFrame(TrajFrame&, const Simulation&) = 0;
};


// reader of ReaxFF trajectory file generated by LAMMPS 'pair_style reax/c'
class ReaxTrajReader : public TrajReader {
public:

	struct Config {
		std::map<int, std::vector<double> > bondorder_cutoff;
		std::vector<double> bondorder_cutoff_default;

		int read_atompos;
		int count_bondorder;
		Config() : read_atompos(0) {
		}
	};
	
	std::map<int, std::vector<double> > bondorders; // Bond order statistics

	ReaxTrajReader() : config() {
	}
	ReaxTrajReader(const Config& c) : config(c){
	}

	int Open(const std::string& path);

	int ReadTrjHead(Simulation* simulation);

	int ReadTrjFrame(TrajFrame&, const Simulation&);

private:
	Config config;
	std::ifstream trjfile;
	int frameCount;
};

#define DefaultTrajReader ReaxTrajReader // <================= Change your own trajectory reader here
#define DefaultTrajConfig ReaxTrajReader::Config // <========== Change your own trajectory reader config here

#endif // !TRAJECTORY_H
