// New class for ui
// Add in version 2.2.2; 

#pragma once
#ifndef REAXDETECT_H
#define REAXDETECT_H


#include "includes.h"
#include "confighandler.h"
#include "simulation.h"
#include "reaxreader.h"
#include "analyzer.h"

#define MYVERSION 30100

// main entry for the program
class ReaxDetect {
public:
	ReaxDetect();

	// initialize arguments
	int init(int argc, char** argv);

	// execute
	int exec();

private:

	int read_opt(int argc, char** argv);

	void set_default_opt();

	void translate_opt();

	void display_help();

	void display_version();

	static int error(int errorcode);

	Simulation simulation;
	ConfigReader cfg_reader;
	ReaxAnalyzer::Config config_analyzer;
	TrajReader::Config config_traj;
	ReaxReader::Config config_reax;

	string input_path;
	string input_name;
	string config_path;

	const string default_config_path = "reacdetect.ini";

	// options inside this program

	int read_line;

	enum {
		TRAJECTORY_FILE, RAWDATA_FILE
	} file_type;

	typedef struct {
		char write_rawdata;
		char write_fulldata;
		char write_kineticfile;
	} WriteOption;
	WriteOption write_option;
};

#endif // !REAXDETECT_H
