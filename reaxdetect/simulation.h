#pragma once

#ifndef SIMULATION_H
#define SIMULATION_H

#include <vector>

// record arguments for a simulation
struct Simulation {
	char name[37];
	double temp;
	int atomNumber;
	double timeStep;
	double volume;
	std::vector<double> atomWeights;
};


#endif // SIMULATION_H
