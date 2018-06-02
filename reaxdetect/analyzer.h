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
	Array rp, rm;

	ReaxAnalyzer() {
	}

	void HandleData(const ReaxReader& rs, const Simulation&);

private:

	void CalcMolLife(const ReaxReader& rs, double timeStep);

};

#endif // !ANALYZER_H
