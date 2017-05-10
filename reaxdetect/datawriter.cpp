
#include "datawriter.h"
#include "filenamehandler.h"
#include "numarray.h"
#include "stringconvert.h"
#include "elementname.h"

#define FILENAME_LENGTH 37
#define EMPTY_BUFFER_COUNT	256


int ReaxDataWriter::Dump(const string& path, const string& path_reac, const ReaxReader& reader) {

	FILE* outfile_c = fopen(path.c_str(), "w");
	FILE* outfile_r = fopen(path_reac.c_str(), "w");

	if (!outfile_c || !outfile_r)return 1;

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

int ReaxDataWriter::WriteSample(const string& sample_path, const ReaxAnalyzer& analyzer)
{
	ofstream outfile(sample_path, ios_base::out);
	outfile << "SPECIES\nt,";
	outfile << join(analyzer.species) << "\n";
	for (const auto& sample : analyzer.samples) {
		outfile << sample.t << "," << join(sample.c) << "\n";
	}
	outfile << "END\nREACTIONS\nt,";
	outfile << join(analyzer.reactions) << "\n";
	for (const auto& sample : analyzer.samples) {
		outfile << sample.t << "," << join(sample.r) << "\n";
	}
	outfile << "END\nINIT\n";
	for (const auto& i : analyzer.init) {
		outfile << i.first << "," << i.second << "\n";
	}
	outfile << "END";
	outfile.close();
	return 0;
}

int ReaxDataWriter::WriteReport(const string & path, const ReaxAnalyzer & analyzer)
{
	//ofstream outfile(path, ios_base::out);
	FILE* outfile = fopen(path.c_str(), "w");
	if (!outfile)return 1;

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

int ReaxDataWriter::WriteRawReactionFreq(const string & path, const ReaxAnalyzer & analyzer)
{

	FILE* outfile = fopen(path.c_str(), "w");
	if (!outfile)return 1;

	fprintf(outfile, "Reaction,Freq_pro,Base_pro,Freq_con,Base_con\n");
	for (size_t i = 0; i < analyzer.reactions.size(); i++) {
		fprintf(outfile, "%s,%d,%u,%d,%u\n", 
			analyzer.reactions[i].c_str(), analyzer.rp[i], analyzer.sum_product_p[i], analyzer.rm[i], analyzer.sum_product_m[i]);
	}
	fclose(outfile);
	return 0;
}

int ReaxDataWriter::WriteBondOrder(const string & path, const ReaxReader & reader, const Simulation& simulation)
{
	ofstream outfile(path, ios_base::out);
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
