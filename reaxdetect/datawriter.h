#pragma once

#ifndef DATAWRITER_H
#define DATAWRITER_H

#include "analyzer.h"
#include "simulation.h"
#include "reaxreader.h"

class ReaxDataWriter {
public:

	// dump full concentration
	static int Dump(const string& species_path, const string& reac_path, const ReaxReader&);

	// sample
	static int WriteSample(const string&, const ReaxAnalyzer&);

	static int WriteReport(const string& path, const ReaxAnalyzer&);

	static int WriteRawReactionFreq(const string& path, const ReaxAnalyzer&);

	static int WriteBondOrder(const string& path, const TrajReader&);

	static int WriteConfig(const string & path, const Simulation&, const ReaxReader&);
};

#endif // !DATAWRITER_H
