
#ifndef REAXDETECT_H
#define REAXDETECT_H

#include "analyzer.h"
#include "reaxreader.h"
#include "simulation.h"
#include "util/serialize.h"

#include <string>
#include <map>

using namespace std;

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

	void display_version()const;

	static int read_boc(const string& filename, map<int, Arrayd>&);

	Simulation simulation;
    Serializer config;
	DefaultTrajConfig config_traj;
	ReaxReader::Config config_reax;

	string input_path;
	string input_name;
	string config_path;

    const char* version_str = "2018v2 alpha";
};

#endif // !REAXDETECT_H
