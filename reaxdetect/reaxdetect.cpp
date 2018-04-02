
#ifdef __unix__
#include <unistd.h>
#include <getopt.h>
#else
#include "getopt.h"
#endif

#include <stdlib.h>
#include "analyzer.h"
#include "datawriter.h"
#include "errors.h"
#include "reaxdetect.h"
#include "reaxreader.h"
#include "util/elements.h"
#include "util/path.h"
#include "util/strutil.h"


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

	read_opt(argc, argv);

	if (!path::exists(input_path)) {
        throw ReaxDetectError("Input file does not exist");
	}

	vector<string> input_path_split = path::split(input_path);
	input_name = input_path_split[0] + input_path_split[1];
	translate_opt();
	return 0;
}

void ReaxDetect::translate_opt()
{
	config_traj = TrajReader::Config();
	config_traj.read_atompos = stob(cfg_reader.at("ReadAtomPos"));
	config_traj.count_bondorder = stoi(cfg_reader.at("CountBondOrder"));
	for (const auto& c : split(cfg_reader.at("BondOrderCutoffDefault"), -1, ',')) {
		config_traj.bondorder_cutoff_default.push_back(stod(c));
	}
	string boc_path;
	if ((boc_path = cfg_reader.at("BondOrderCutoff")) != "default") {
		read_boc(boc_path, config_traj.bondorder_cutoff);
	}

    config_reax = ReaxReader::Config();
	config_reax.buffer_size = stoul(cfg_reader.at("FrameBufferSize"));
	config_reax.buffer_interval = stoi(cfg_reader.at("FrameBufferInterval"));
	config_reax.recognize_interval = stoi(cfg_reader.at("RecognizeInterval"));
	config_reax.recognize_limit = stoi(cfg_reader.at("RecognizeLimit"));
	config_reax.recognize_begin = stoi(cfg_reader.at("RecognizeBegin"));

	config_analyzer = ReaxAnalyzer::Config();
	string sample_method_str = cfg_reader.at("SampleMethod");
	if (sample_method_str == "fixint") {
		config_analyzer.sample_method = SAMPLE_FIXINT;
	}
	config_analyzer.sample_int = stod(cfg_reader.at("SampleInterval"));
	config_analyzer.sample_range = stod(cfg_reader.at("SampleRange"));
}

int ReaxDetect::exec() {

	printf("Reading Trajectory File...\n");
	ReaxDataWriter writer;
	ReaxReader reader(config_reax);

	TrajReader traj(config_traj);
    traj.Open(input_path);
    traj.ReadTrjHead(&simulation);

	printf("Reading Frames...\n");
	reader.HandleData(traj, simulation);

	printf("Writing raw results...\n");
	if (write_option.write_fulldata) {
        writer.Dump(input_name + "_full_dump.csv", input_name + "_full_reac.csv", reader);
	}
	writer.WriteConfig(input_name + "_config.txt", simulation, reader);
	if (cfg_reader.get("CountBondOrder", "0") != "0") {
		writer.WriteBondOrder(input_name + "_bondorder.csv", traj);
	}
	if (write_option.write_kineticfile) {
		printf("Analysing Kinetic Results...\n");
		simulation.timeStep *= config_reax.recognize_interval;
		ReaxAnalyzer analyzer(config_analyzer);
		analyzer.HandleData(reader, simulation);
		writer.WriteReport(input_name + "_full_report.csv", analyzer);
		writer.WriteRawReactionFreq(input_name + "_rawfreq.csv", analyzer);
        writer.WriteSample(input_name + "_sample.csv", analyzer);
	}
	printf("Finished Successfully.\n");
	return 0;
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
		{ "config", no_argument, NULL, 0},
		{ "dump", required_argument, NULL, 0},
		{ "help", no_argument, NULL, 'h' },
		{ "version", no_argument, NULL, 0 },
		{ "", 0, 0, 0}
	};
	int longIndex;
	int opt = getopt_long(argc, argv, "c:v:s:t:b:f:l:h", longOpts, &longIndex);
	while (opt != -1) {
		switch (opt)
		{
		case 'b':
			cfg_reader["FrameBufferSize"] = string(optarg);
			break;
		case 'c':
			config_path = string(optarg);
			if (path::exists(config_path)) {
				cfg_reader.read_data(config_path);
			}
			else {
                throw ReaxDetectError("Bad config path");
                return 1;
			}
			break;
		case 'f':
			cfg_reader["RecognizeInterval"] = string(optarg);
			break;
		case 'l':
			cfg_reader["RecognizeLimit"] = string(optarg);
		case 's':
			cfg_reader["SampleInterval"] = string(optarg);
			cfg_reader["SampleRange"] = string(optarg);
			break;
		case 't':
			simulation.timeStep = atof(optarg);
			break;
		case 'v':
			simulation.volume = atof(optarg);
			break;
		case '?':
		case 'h':
			display_help();
            exit(0);
		case 0:
			if (string(longOpts[longIndex].name) == "version") {
				display_version();
                exit(0);
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
            // write config file only
			else if (string(longOpts[longIndex].name) == "config") {
				cfg_reader.fresh_data(path::split(argv[0])[0] + default_config_path);
                exit(0);
			}
			break;
		default:
			break;
		}
		opt = getopt_long(argc, argv, "b:c:f:l:s:t:v:h", longOpts, &longIndex);
	}
	input_path = string(argv[argc - 1]);
	return 0;
}

void ReaxDetect::set_default_opt()
{
	cfg_reader["BondOrderCutoff"] = "default";
	cfg_reader["BondOrderCutoffDefault"] = "0.5,0.5";
	cfg_reader["CountBondOrder"] = "0";
	cfg_reader["FrameBufferSize"] = "2";
	cfg_reader["FrameBufferInterval"] = "1";
	cfg_reader["ReadAtomPos"] = "false";
	cfg_reader["RecognizeBegin"] = "0";
	cfg_reader["RecognizeInterval"] = "1";
	cfg_reader["RecognizeLimit"] = "-1";
	cfg_reader["SampleMethod"] = "fixint";
	cfg_reader["SampleInterval"] = "200.0";
	cfg_reader["SampleRange"] = "200.0";
}

void ReaxDetect::display_version() {
	printf("Reaction Detector: Version %d.%d.%d\n", MYVERSION / 10000, (MYVERSION / 100) % 100, MYVERSION % 100);
}
void ReaxDetect::display_help() {
	printf("Reaction Detector:\n"
	    "reaxdetect [-b buffersize] [-c config] [-f step] [-l limit]\n"
        "[-s samplerange] [-t timestep] [-v volume] [-h]\n"
	    "[--config] [--dump=nodump/full] [--help] [--version] trajectory\n");
}
