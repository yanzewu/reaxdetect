#pragma once

#ifndef ANALYZER_H
#define ANALYZER_H

#include "simulation.h"
#include "reaxreader.h"

    // Basic analyzation and conclusion of simulation trajectory
class ReaxAnalyzer {
public:

	vector<string> elements;
	vector<string> species;
	vector<string> reactions;

	Arrayd species_life;
	Array freq_f, freq_b;   // Forward and backward occurences

	ReaxAnalyzer() {
	}

    // Handles data of ReaxReader, which must be invoked ReaxReader.HandleData first
	void HandleData(const ReaxReader&, const Simulation&);

private:

	void CalcMolLife(const ReaxReader&, double timeStep);

};

#endif // !ANALYZER_H
