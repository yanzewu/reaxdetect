
#include "datawriter.h"
#include "errors.h"
#include "util/elements.h"
#include "util/path.h"
#include "util/strutil.h"

#define FILENAME_LENGTH 37
#define EMPTY_BUFFER_COUNT	256


int ReaxDataWriter::Dump(const string& path, const string& path_reac, const ReaxReader& reader) {

	FILE* outfile_c = fopen(path.c_str(), "w");
	FILE* outfile_r = fopen(path_reac.c_str(), "w");

    if (!outfile_c)throw IOError(path);
    if (!outfile_r)throw IOError(path_reac);

	fprintf(outfile_c, "t,%s", join(reader.species).c_str());
	fprintf(outfile_r, "t");
	for (const auto& s : reader.reactions) {
		fprintf(outfile_r, ",%s,%s", s.to_string(reader.species).c_str(), (-s).to_string(reader.species).c_str());
	}

	for (const auto& fs : reader.fss) {
		fprintf(outfile_c, "\n%.2f", fs.t);
		for (const auto& c : fs.mol_freq) {
			fprintf(outfile_c, ",%d", c);
		}

		fprintf(outfile_r, "\n%.2f", fs.t);
		for (const auto& r : fs.reaction_freq) {
			fprintf(outfile_r, ",%d,%d", r.first, r.second);
		}

	}

	fclose(outfile_c);
	fclose(outfile_r);

	return 0;
}

int ReaxDataWriter::WriteReport(const string & path, const ReaxAnalyzer & analyzer)
{
	//ofstream outfile(path, ios_base::out);
	FILE* outfile = fopen(path.c_str(), "w");
    if (!outfile)throw IOError(path);

	fprintf(outfile, "Molecules:\nindex,name,life(ps)\n");
	for (size_t i = 0; i < analyzer.species.size(); i++) {
		fprintf(outfile, "%zd,%s,%.3f\n", i, analyzer.species[i].c_str(), analyzer.species_life[i]);
	}
	fprintf(outfile, "Reactions\nindex,reaction,freqplus,freqminus\n");
	for (size_t i = 0; i < analyzer.reactions.size(); i++) {
		fprintf(outfile, "%zd,%s,%d,%d\n", i, analyzer.reactions[i].c_str(), analyzer.rp[i], analyzer.rm[i]);
	}
	fclose(outfile);
	return 0;
}

int ReaxDataWriter::WriteBondOrder(const string & path, const TrajReader & reader)
{
	ofstream outfile(path, ios_base::out);
    if (!outfile.is_open()) throw IOError(path);

	for (const auto& entry : reader.bondorders) {
		outfile << WeightName[entry.first / MAX_ATOM_TYPE] << '-';
		outfile << WeightName[entry.first % MAX_ATOM_TYPE];
		for (const auto& b : entry.second) {
			outfile << ',' << b;
		}
		outfile << endl;
	}
	outfile.close();

	return 0;
}

int ReaxDataWriter::WriteConfig(const string & path, const Simulation & simulation, const ReaxReader& reader)
{
	ofstream outfile(path, ios_base::out);
    if (!outfile)throw IOError(path);

	outfile << "interval=" << reader.fss[1].t - reader.fss[0].t << "\n";
	outfile << "volume=" << simulation.volume << "\n";
	return 0;
}
