
#ifdef __unix__
#include <unistd.h>
#include <getopt.h>
#else
#include "getopt.h"
#endif

#include <stdlib.h>
#include "reaxdetect.h"
#include "reaxreader.h"
#include "analyzer.h"
#include "datawriter.h"
#include "filenamehandler.h"
#include "elementname.h"

#define ERROR_NO_INPUT	0x01
#define ERROR_BAD_INPUT 0x02
#define ERROR_BAD_ARGS	0x03
#define ERROR_BAD_OUTPUT	0x04
#define ERROR_NO_EXEC	0x08

ReaxDetect::ReaxDetect() {
}

int ReaxDetect::init(int argc, char** argv) {

	if (argc < 2) {
		display_version();
		exit(0);
		return 1;
	}

	set_default_opt();
	cfg_reader.read_data(path::split(argv[0])[0] + default_config_path);

	write_option.write_fulldata = 1;
	write_option.write_kineticfile = 1;
	write_option.write_rawdata = 1;

	simulation.timeStep = 0.0;
	simulation.volume = 1.0;

	int result = read_opt(argc, argv);
	if (result) {
		error(result);
	}

	if (!path::exists(input_path)) {
		error(ERROR_BAD_INPUT);
	}

	vector<string> input_path_split = path::split(input_path);
	input_name = input_path_split[0] + input_path_split[1];
	translate_opt();
	return 0;
}

void ReaxDetect::translate_opt()
{
	config_traj = TrajReader::Config();
	config_traj.read_atompos = stob(cfg_reader.get("ReadAtomPos", "false"));
	config_traj.count_bondorder = stob(cfg_reader.get("CountBondOrder", "false"));
	for (const auto& c : split(cfg_reader.get("BondOrderCutoffDefault", "0.3,0.5"), -1, ',')) {
		config_traj.bondorder_cutoff_default.push_back(stod(c));
	}
	string boc_path;
	if ((boc_path = cfg_reader.get("BondOrderCutoff", "default")) != "default") {
		read_boc(boc_path, config_traj.bondorder_cutoff);
	}

	config_reax.buffer_size = stoul(cfg_reader.get("FrameBufferSize", "2"));
	config_reax.recognize_interval = stoi(cfg_reader.get("RecognizeInterval", "1"));
	config_reax.recognize_limit = stoi(cfg_reader.get("RecognizeLimit", "-1"));
	config_reax.recognize_begin = stoi(cfg_reader.get("RecognizeBegin", "0"));

	config_analyzer = ReaxAnalyzer::Config();
	string sample_method_str = cfg_reader.get("SampleMethod", "fixint");
	if (sample_method_str == "fixint") {
		config_analyzer.sample_method = SAMPLE_FIXINT;
	}
	config_analyzer.sample_int = stoi(cfg_reader.get("SampleInterval", "1000"));
	config_analyzer.sample_range = stoi(cfg_reader.get("SampleRange", "1000"));
}

int ReaxDetect::exec() {

	printf("Reading Trajectory File...\n");
	ReaxDataWriter writer;
	ReaxReader reader(config_reax);

	TrajReader traj(config_traj);
	if (traj.Open(input_path) || traj.ReadTrjHead(&simulation)){
		error(ERROR_BAD_INPUT);
	}
	printf("Reading Frames...\n");
	reader.HandleData(traj, simulation);

	printf("Writing raw results...\n");
	if (write_option.write_fulldata && writer.Dump(input_name + "_full_dump.csv", input_name + "_full_reac.csv", reader)) {
		error(ERROR_BAD_OUTPUT);
	}
	writer.WriteConfig(input_name + "_config.txt", simulation, reader);
	if (cfg_reader.get("CountBondOrder", "false") == "true") {
		writer.WriteBondOrder(input_name + "_bondorder.csv", traj);
	}
	if (write_option.write_kineticfile) {
		printf("Analysing Kinetic Results...\n");
		ReaxAnalyzer analyzer(config_analyzer);
		analyzer.HandleData(reader, simulation);
		writer.WriteReport(input_name + "_full_report.csv", analyzer);
		writer.WriteRawReactionFreq(input_name + "_rawfreq.csv", analyzer);
		if (writer.WriteSample(input_name + "_sample.csv", analyzer)) {
			error(ERROR_BAD_OUTPUT);
		}
	}
	printf("Finished Successfully.\n");
	return 0;
}

int ReaxDetect::error(int errorcode) {
	switch (errorcode)
	{
	case ERROR_NO_INPUT:
		printf("Error: Please choose an input file!\n");
		break;
	case ERROR_BAD_INPUT:
		printf("Error: Bad input file.\n");
		break;
	case ERROR_BAD_ARGS:
		printf("Error: Please check your arguments.\n");
		break;
	case ERROR_BAD_OUTPUT:
		printf("Error: Cannot write output file.\n");
		break;
	case ERROR_NO_EXEC:
		exit(0);
		break;
	default:
		break;
	}
	exit(errorcode);
	return 1;
}

int ReaxDetect::read_boc(const string & filename, map<int, Arrayd>& bondorder_cutoff)
{
	ConfigReader boc_reader;
	boc_reader.read_data(filename);
	for (const auto& item : boc_reader._items) {
		vector<string> atom_names = split(item.first, 1, '-');
		int type1 = NameWeight.at(atom_names[0]), type2 = NameWeight.at(atom_names[1]);
		int bond_identifier = type1 >= type2 ? type1 * MAX_ATOM_TYPE + type2 : type2 * MAX_ATOM_TYPE + type1;
		for (const auto& c : split(item.second, -1, ',')) {
			bondorder_cutoff[bond_identifier].push_back(stof(c));
		}
	}
	return 0;
}

// bottom api

int ReaxDetect::read_opt(int argc, char** argv) {
	const option longOpts[] = {
		{ "version", no_argument, NULL, 0 },
		{ "help", no_argument, NULL, 'h' },
		{ "config", no_argument, NULL, 0},
		{ "dump", required_argument, NULL, 0},
		{ "", 0, 0, 0}
	};
	int longIndex;
	int opt = getopt_long(argc, argv, "c:v:s:t:b:f:l:h", longOpts, &longIndex);
	while (opt != -1) {
		switch (opt)
		{
			// config
		case 'c':
			config_path = string(optarg);
			if (path::exists(config_path)) {
				cfg_reader.read_data(config_path);
			}
			else {
				return ERROR_BAD_ARGS;
			}
			break;
			// volume
		case 'b':
			cfg_reader["FrameBufferSize"] = string(optarg);
			break;
		case 'f':
			cfg_reader["RecognizeInterval"] = string(optarg);
			break;
		case 'l':
			cfg_reader["RecognizeLimit"] = string(optarg);
		case 'v':
			simulation.volume = atof(optarg);
			break;
		case 's':
			cfg_reader["SampleInterval"] = string(optarg);
			cfg_reader["SampleRange"] = string(optarg);
			break;
		case 't':
			simulation.timeStep = atof(optarg);
			break;
		case '?':
		case 'h':
			display_help();
			return ERROR_NO_EXEC;
		case 0:
			if (string(longOpts[longIndex].name) == "version") {
				display_version();
				return ERROR_NO_EXEC;
			}
			else if (string(longOpts[longIndex].name) == "dump") {
				if (string(optarg) == "nodump") {
					write_option.write_fulldata = 0;
					write_option.write_rawdata = 0;
				}
				else if (string(optarg) == "full") {
					write_option.write_fulldata = 2;
					write_option.write_rawdata = 2;
				}
			}
			else if (string(longOpts[longIndex].name) == "config") {
				cfg_reader.fresh_data(path::split(argv[0])[0] + default_config_path);
				return ERROR_NO_EXEC;
			}
			break;
		default:
			break;
		}
		opt = getopt_long(argc, argv, "c:v:s:t:b:f:l:h", longOpts, &longIndex);
	}
	input_path = string(argv[argc - 1]);
	return 0;
}

void ReaxDetect::set_default_opt()
{
	cfg_reader["BondOrderCutoffDefault"] = "0.3,0.5";
	cfg_reader["BondOrderCutoff"] = "default";
	cfg_reader["SampleMethod"] = "fix_int";
	cfg_reader["SampleInterval"] = "1000";
	cfg_reader["SampleRange"] = "1000";
	cfg_reader["ReadAtomPos"] = "false";
	cfg_reader["FrameBufferSize"] = "2";
	cfg_reader["RecognizeInterval"] = "1";
	cfg_reader["RecognizeBegin"] = "0";
	cfg_reader["RecognizeLimit"] = "-1";
	cfg_reader["CountBondOrder"] = "false";
}

void ReaxDetect::display_version() {
	printf("Reaction Detector: Version %d.%d.%d\n", MYVERSION / 10000, (MYVERSION / 100) % 100, MYVERSION % 100);
}
void ReaxDetect::display_help() {
	printf("Reaction Detector:\n");
	printf("reaxdetect [-v volume] [-s samplerange] [-b buffersize] [-c config] [-h help] input_file\n");
	printf("--version, --help, --dump [dumpoption=nodump/full] --config\n");
}
