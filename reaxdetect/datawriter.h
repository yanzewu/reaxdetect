#pragma once

#ifndef DATAWRITER_H
#define DATAWRITER_H

#include "simulation.h"
#include "reaxreader.h"
#include "analyzer.h"

class ReaxDataWriter {
public:

	// dump full concentration
	static int Dump(const string& species_path, const string& reac_path, const ReaxReader&);

	// sample
	static int WriteSample(const string&, const ReaxAnalyzer&);

	static int WriteReport(const string& path, const ReaxAnalyzer&);

	static int WriteKineticFile(const string& path, const ReaxAnalyzer&);

	static int WriteRawReactionFreq(const string& path, const ReaxAnalyzer&);

	static int WriteBondOrder(const string& path, const ReaxReader&, const Simulation&);
};

#endif // !DATAWRITER_H
