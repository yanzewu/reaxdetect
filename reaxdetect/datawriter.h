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

	static int WriteReport(const string& path, const ReaxAnalyzer&);

	static int WriteBondOrder(const string& path, const ReaxTrajReader&);

	static int WriteConfig(const string & path, const Simulation&, const ReaxReader&);
};

#endif // !DATAWRITER_H
