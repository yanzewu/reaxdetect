
#ifndef SIMULATION_H
#define SIMULATION_H

#include <vector>

    // Metadata of simulation (either read from trajectory or input by user)
struct Simulation {
	int atomNumber;     // Number of atoms
	double timeStep;    // Frame timestep
	double volume;      // Volume
	std::vector<double> atomWeights;    // atom Weight List (index=>weight)
};


#endif // SIMULATION_H
