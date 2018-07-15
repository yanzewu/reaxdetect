
#ifdef __unix__
#include <unistd.h>
#include <getopt.h>
#else
#include "getopt.h"
#endif

#include "analyzer.h"
#include "datawriter.h"
#include "errors.h"
#include "reaxdetect.h"
#include "reaxreader.h"

#include "util/elements.h"
#include "util/path.h"
#include "util/strutil.h"
#include "util/serialize.h"

#include <stdlib.h>


ReaxDetect::ReaxDetect() {
}

int ReaxDetect::init(int argc, const char** argv) {

	set_default_opt();
	read_opt(argc, argv);
    translate_opt();

	if (!path::exists(input_path)) {
        throw ReaxDetectError("Input file does not exist");
	}

	vector<string> input_path_split = path::split(input_path);
	input_name = input_path_split[0] + input_path_split[1];
	
	return 0;
}

void ReaxDetect::translate_opt()
{
    simulation.timeStep = config.at<double>("timestep");
    simulation.volume = config.at<double>("volume");

    input_path = config.at<string>("file");

    /* NOTICE: If you want to use your own configer. You may change the code 
        about ``config_traj`` below in order to receive UI input.
    */
    config_traj = DefaultTrajConfig();
	config_traj.read_atompos = config.at<int>("ReadAtomPos");
	config_traj.count_bondorder = config.at<int>("CountBondOrder");
	for (const auto& c : split(config.at<string>("BondOrderCutoffDefault"), -1, ',')) {
		config_traj.bondorder_cutoff_default.push_back(stod(c));
	}
	string boc_path;
	if ((boc_path = config.at<string>("BondOrderCutoff")) != "default") {
		read_boc(boc_path, config_traj.bondorder_cutoff);
	}

    config_reax = ReaxReader::Config();
	config_reax.buffer_size = config.at<size_t>("FrameBufferSize");
	config_reax.recognize_interval = config.at<int>("RecognizeInterval");
	config_reax.recognize_limit = config.at<int>("RecognizeLimit");
	config_reax.recognize_begin = config.at<int>("RecognizeBegin");
}

int ReaxDetect::exec() {

	printf("Reading Trajectory File...\n");
	ReaxDataWriter writer;
	ReaxReader reader(config_reax);

    DefaultTrajReader traj(config_traj);
    traj.Open(input_path);
    traj.ReadTrjHead(&simulation);      // read necessary metadata

	printf("Reading Frames...\n");
	reader.HandleData(traj, simulation);    // main read loop & handle

	printf("Writing raw results...\n");
    writer.Dump(input_name + "_full_dump.csv", input_name + "_full_reac.csv", reader);  // Write species/reaction rates

	writer.WriteConfig(input_name + "_config.txt", simulation, reader); // Write conditional information (t. V)

	if (config.at<int>("CountBondOrder") != 0) {     // Used when writing PDF of bond order
		writer.WriteBondOrder(input_name + "_bondorder.csv", traj);
	}
	printf("Analysing Kinetic Results...\n");
	simulation.timeStep *= config_reax.recognize_interval;  // Scale time for correctly calculating lifetime

	ReaxAnalyzer analyzer;
	analyzer.HandleData(reader, simulation);
	writer.WriteReport(input_name + "_full_report.csv", analyzer);
	
    printf("Finished Successfully.\n");
	return 0;
}


int ReaxDetect::read_boc(const string & filename, map<int, Arrayd>& bondorder_cutoff)
{
	
    std::ifstream infile(filename);
    char _buffer[128];

    while (infile.getline(_buffer, 128)) {
        vector<string> sp = split(_buffer);
		vector<string> atom_names = split(sp[0], 1, '-');
		int type1 = NameWeight.at(atom_names[0]), type2 = NameWeight.at(atom_names[1]);
		int bond_identifier = type1 >= type2 ? type1 * MAX_ATOM_TYPE + type2 : type2 * MAX_ATOM_TYPE + type1;
		for (const auto& c : split(sp[1], -1, ',')) {
			bondorder_cutoff[bond_identifier].push_back(stof(c));
		}
    }

    infile.close();

	return 0;
}

// bottom api

int ReaxDetect::read_opt(int argc, const char** argv) {

    ArgParser parser;

    // Set command line argument here. The usage is similar to ``argparse`` in Python.
    parser.add_argument("-c", "--config", 1);
    parser.add_argument("-b", "--RecognizeBegin", 1);
    parser.add_argument("-s", "--RecognizeInterval", 1);
    parser.add_argument("-m", "--RecognizeLimit", 1);
    parser.add_argument("-t", "--timestep", 1);
    parser.add_argument("-v", "--volume", 1);
    parser.add_argument("-h", "--help", 0);
    parser.add_argument("--version", 0);
    parser.add_argument("file", 1);

    ArgList opts;
    try {
        opts = parser.parse(argc, argv);
    }
    catch (const ArgParseError& e) {
        std::cerr << e.what() << std::endl;
        parser.display_help(std::cerr);
        exit(1);
    }

    if (opts.count("help") > 0) {
        parser.display_help(std::cerr);
        exit(0);
    }
    if (opts.count("version") > 0) {
        display_version();
        exit(0);
    }
    if (opts.count("config") > 0) {
        config.read_file(opts.at("config")[0]);
    }

    if (opts.count("file") == 0) {
        std::cerr << "Missing argument: FILE" << std::endl;
        parser.display_help(std::cerr);
        exit(1);
    }

    config.read_arg(opts);

	return 0;
}

void ReaxDetect::set_default_opt()
{
    /* Set the default argument here. The input entry for argument is automatically created.
    */
    config.set_default({
        {"BondOrderCutoff", "default"},
        {"BondOrderCutoffDefault", "0.5,0.5"},
        {"CountBondOrder", 0},
        {"FrameBufferSize", 2},
        {"ReadAtomPos", 0},

        {"RecognizeBegin", 0},
        {"RecognizeInterval", 1},
        {"RecognizeLimit", (unsigned)-1},

        {"timestep", 1.0},
        {"volume", 1.0},
        {"file", ""}
    });
}

void ReaxDetect::display_version()const {
	printf("Reaction Detector: Version %s\n", version_str);
}

