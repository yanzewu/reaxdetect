#pragma once

#ifndef DATAWRITER_H
#define DATAWRITER_H

#include "simulation.h"
#include "reaxreader.h"
#include "analyzer.h"

class ReaxDataWriter {
public:
	// read data from rwdx file
	static int ReadData(const string&, ReaxReader&, Simulation*);

	// write data into rwdx file
	static int WriteData(const string&, const ReaxReader&, const Simulation&);

	// dump full concentration
	static int Dump(const string& species_path, const string& reac_path, const ReaxReader&, double interval);

	// sample
	static int WriteSample(const string&, const ReaxAnalyzer&);

	static int WriteReport(const string& path, const ReaxAnalyzer&);

	static int WriteKineticFile(const string& path, const ReaxAnalyzer&);

	static int WriteRawReactionFreq(const string& path, const ReaxAnalyzer&);
};

#endif // !DATAWRITER_H
