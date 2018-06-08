/* ReaxDetect 
 * A tool analyzing LAMMPS trajectory file with Reaction Force Field.
 * Written by Yanze Wu (https://github.com/flipboards)
*/

#pragma once
#ifndef REAXDETECT_H
#define REAXDETECT_H

#include "analyzer.h"
#include "reaxreader.h"
#include "simulation.h"
#include "util/serialize.h"

#include <string>
#include <map>

using namespace std;

#define MYVERSION 30101

// main entry for the program
class ReaxDetect {
public:
	ReaxDetect();

	// initialize arguments
	int init(int argc, const char** argv);

	// execute
	int exec();

private:

	int read_opt(int argc, const char** argv);

	void set_default_opt();

	void translate_opt();

	static void display_version();

	static int read_boc(const string& filename, map<int, Arrayd>&);

	Simulation simulation;
    Serializer config;
	ReaxTrajReader::Config config_traj;
	ReaxReader::Config config_reax;

	string input_path;
	string input_name;
	string config_path;

	const string default_config_path = "reacdetect.ini";

};

#endif // !REAXDETECT_H
