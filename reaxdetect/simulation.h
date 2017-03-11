#pragma once

#ifndef SIMULATION_H
#define SIMULATION_H


#include "includes.h"

// record arguments for a simulation
struct Simulation {
	char name[37];
	double temp;
	int atomNumber;
	double timeStep;
	double volume;
	vector<double> atomWeights;
};


#endif // SIMULATION_H
