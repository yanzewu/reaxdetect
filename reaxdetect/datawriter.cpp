
#include <iomanip>

#include "datawriter.h"
#include "filenamehandler.h"
#include "numarray.h"
#include "stringconvert.h"

#define FILENAME_LENGTH 37
#define RWDXHEADER	"RAWDX2FF"
#define EMPTY_BUFFER_COUNT	256

const int BlockMark = 0xFFFF;
const int ElementMark = 0xFF00;

struct RwdxConfig {
	int dense_freq;
};

int ReaxDataWriter::ReadData(const string & path, ReaxReader & reader, Simulation * simulation)
{
	ifstream infile(path.c_str(), ios_base::in | ios_base::binary);
	if (!infile.is_open())return 1;
	int buffer[256]; tsize_t size;
	infile.read((char*)buffer, 8);
	if (string((char*)buffer, (char*)buffer + 8) != RWDXHEADER)return 1;//check file head

																		//read filehead
	infile.read(simulation->name, FILENAME_LENGTH);
	infile.read((char*)&simulation->timeStep, sizeof(simulation->timeStep));
	infile.read((char*)&simulation->volume, sizeof(simulation->volume));
	vread(infile, simulation->atomWeights);

	RwdxConfig rwdx;

	infile.read((char*)buffer, EMPTY_BUFFER_COUNT);
	memcpy(&rwdx, buffer, sizeof(rwdx));
	infile.read((char*)buffer, sizeof(int)); if (buffer[0] != BlockMark)return 1;

	//read molecule map
	reader.molecules.clear();
	infile.read((char*)&size, sizeof(size)); reader.molecules.resize(size);
	for (auto& s : reader.molecules) {
		infile.read((char*)&size, sizeof(size));
		infile.read((char*)buffer, size);
		s.read_bin((char*)buffer);
		infile.read((char*)buffer, sizeof(int)); if (buffer[0] != ElementMark)return 1;//endmark check
	}
	infile.read((char*)buffer, sizeof(int)); if (buffer[0] != BlockMark)return 1;//Blockmark check

																					//read Reaction map
	reader.reactions.clear();
	infile.read((char*)&size, sizeof(size)); reader.reactions.resize(size);
	for (auto& reaction : reader.reactions) {
		vread(infile, reaction.reagants);
		vread(infile, reaction.products);
		infile.read((char*)buffer, sizeof(int)); if (buffer[0] != ElementMark)return 1;
	}
	infile.read((char*)buffer, sizeof(int)); if (buffer[0] != BlockMark)return 1;

	//read framestatistic
	reader.fss.clear();
	infile.read((char*)&size, sizeof(size)); reader.fss.resize(size);
	for (auto& fs : reader.fss) {
		//read molecule frequency
		vread(infile, fs.mol_freq);
		infile.read((char*)buffer, sizeof(ElementMark)); if (buffer[0] != ElementMark)return 1;//endmark check
																								//read  reader.reactions
		if (rwdx.dense_freq) {
			vread(infile, fs.reaction_freq);
		}
		else {
			fs.reaction_freq.resize(reader.reactions.size());
			infile.read((char*)&size, sizeof(size));
			vector<diint> datas(size, diint(0, 0));
			vector<tsize_t> poses(size, 0);
			infile.read((char*)&buffer, size * sizeof(datas[0]));
			datas.assign((diint*)buffer, (diint*)buffer + size);
			infile.read((char*)&buffer, size * sizeof(poses[0]));
			poses.assign((tsize_t*)buffer, (tsize_t*)buffer + size);
			for (size_t i = 0; i < size; i++) {
				fs.reaction_freq[poses[i]] = datas[i];
			}
		}
		infile.read((char*)buffer, sizeof(ElementMark)); if (buffer[0] != ElementMark)return 1;
	}
	infile.close();
	for (auto& mol : reader.molecules) {
		reader.species.push_back(mol.to_smiles());
	}
	return 0;
}

int ReaxDataWriter::WriteData(const string& path, const ReaxReader& reader, const Simulation& simulation)
{
	string crtpath = path::get_path(path);

	ofstream outfile(crtpath.c_str(), ios_base::out | ios_base::binary);
	outfile.write(RWDXHEADER, 8);

	RwdxConfig rwdx;
	rwdx.dense_freq = 1;

	//write filehead [changeable section]
	outfile.write(simulation.name, FILENAME_LENGTH);	//37 bytes
	outfile.write((char*)&simulation.timeStep, sizeof(simulation.timeStep));	//8 bytes
	outfile.write((char*)&simulation.volume, sizeof(simulation.volume));	//8 bytes
	vwrite(outfile, simulation.atomWeights);	//vector

	char emptybuffer[EMPTY_BUFFER_COUNT];
	memset(&emptybuffer, 0, EMPTY_BUFFER_COUNT);
	memcpy(&emptybuffer, &rwdx, sizeof(rwdx));
	outfile.write(emptybuffer, EMPTY_BUFFER_COUNT);
	outfile.write((char*)&BlockMark, sizeof(BlockMark));

	//write molecule map
	tsize_t size = reader.molecules.size();
	outfile.write((char*)&size, sizeof(size));
	for (auto sms = reader.molecules.begin(); sms < reader.molecules.end(); sms++) {
		char buffer[255];
		size = sms->to_bin(buffer);
		outfile.write((char*)&size, sizeof(size));
		outfile.write(buffer, size);
		outfile.write((char*)&ElementMark, sizeof(ElementMark));
	}
	outfile.write((char*)&BlockMark, sizeof(int));

	//write reader.reactions
	size = reader.reactions.size();
	outfile.write((char*)&size, sizeof(size));

	for (auto reac = reader.reactions.begin(); reac < reader.reactions.end(); reac++) {
		vwrite(outfile, reac->reagants);
		vwrite(outfile, reac->products);
		outfile.write((char*)&ElementMark, sizeof(ElementMark));
	}
	outfile.write((char*)&BlockMark, sizeof(BlockMark));

	//write framestats
	size = reader.fss.size();
	outfile.write((char*)&size, sizeof(size));

	for (const auto& fs : reader.fss) {
		//write molecule frequency
		vwrite(outfile, fs.mol_freq);
		outfile.write((char*)&ElementMark, sizeof(ElementMark));
		vwrite(outfile, fs.reaction_freq);
		outfile.write((char*)&ElementMark, sizeof(ElementMark));
	}
	outfile.close();
	return 0;
}

int ReaxDataWriter::Dump(const string& path, const string& path_reac, const ReaxReader& reader, double interval) {
	ofstream outfile(path, ios_base::out);
	ofstream outfile_reac(path_reac, ios_base::out);

	outfile << "time," << join(reader.species);
	outfile_reac << "time,";
	for (const auto& s : reader.reactions) {
		outfile_reac << s.to_string(reader.species) << ",";
	}

	double t = 0;
	for (const auto& fs : reader.fss) {
		outfile << "\n" << t << "," << join(fs.mol_freq);
		outfile_reac << "\n" << t << "," << join(net(fs.reaction_freq));
		t += interval;
	}
	outfile.close();
	outfile_reac.close();
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
	outfile << "END";
	outfile.close();
	return 0;
}

int ReaxDataWriter::WriteReport(const string & path, const ReaxAnalyzer & analyzer)
{
	ofstream outfile(path, ios_base::out);
	outfile << "Molecules:\n";
	outfile << "index,\tname,\tlife(ps);" << endl;
	outfile << setprecision(4);
	for (size_t i = 0; i < analyzer.species.size(); i++) {
		outfile << i << "," << analyzer.species[i] << "," << analyzer.species_life[i] << endl;
	}
	outfile << "Reactions\n";
	outfile << "index,\treaction,\tfreqplus,\tfreqminus,\tkp_real,km_real,\tkp_min,kp_max,\tkm_min,km_max\n";
	for (size_t i = 0; i < analyzer.reactions.size(); i++) {
		outfile << i << "," << analyzer.reactions[i] << "," << analyzer.rp[i] << "," << analyzer.rm[i] << "," 
			<< analyzer.kp[i] << "," << analyzer.km[i] << ",";
		outfile << analyzer.kp_range[0][i] << "," << analyzer.kp_range[1][i] << ",";
		outfile << analyzer.km_range[0][i] << "," << analyzer.km_range[1][i] << endl;
	}

	outfile.close();

	return 0;
}

int ReaxDataWriter::WriteKineticFile(const string & path, const ReaxAnalyzer & analyzer)
{
	ofstream outfile(path);
	outfile << "# This file is generated by ReaxDetect.\n";
	outfile << "ELEMENTS\n" << join(analyzer.elements, "\n") << "\nEND\n";

	//species
	outfile << "SPECIES\n" << join(analyzer.species, "\n") << "\nEND\n";

	//reactions
	outfile << "REACTIONS" << endl;
	for (size_t i = 0; i < analyzer.reactions.size(); i++) {
		outfile << analyzer.reactions[i] << "\t" << analyzer.kp[i] << "\t" << analyzer.km[i] << endl;
	}
	outfile << "END" << endl;
	outfile.close();

	return 0;
}
