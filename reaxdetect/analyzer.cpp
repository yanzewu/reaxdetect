#include <math.h>
#include <algorithm>

#include "analyzer.h"
#include "util/elements.h"
#include "util/vecutil.h"

void ReaxAnalyzer::HandleData(const ReaxReader& reader, const Simulation& simulation)
{
	for (const auto& e : simulation.atomWeights) {
		elements.push_back(WeightNameL.at((int)(e + 0.1)));
	}
	printf("Elements are:");
	for (const auto& e : elements) {
		printf(" %s", e.c_str());
	}

    std::map<std::string, double> init;

	for (int i = 0; i < reader.species.size(); i++) {
		if (reader.fss[0].mol_freq[i] == 0)break;
		init[reader.species[i]] = reader.fss[0].mol_freq[i] / simulation.volume;
	}

	printf("\nConcentration of reactants are:\n");
	for (const auto& reactant : init) {
		printf("%s\t%.3e\n", reactant.first.c_str(), reactant.second);
	}

	printf("Encoding reactions...\n");
	species = reader.species;
	reactions.resize(reader.reactions.size());
	for (size_t i = 0; i < reader.reactions.size(); i++) {
		reactions[i] = reader.reactions[i].to_string(species);
	}

	printf("Calculating molecule lifetime...\n");
	CalcMolLife(reader, simulation.timeStep);
	
	printf("Calculating reaction product...\n");
	double interval = reader.fss[1].t - reader.fss[0].t;
	printf("Interval: %.3e\n", interval);
}

void ReaxAnalyzer::CalcMolLife(const ReaxReader& rs, double timeStep)
{
	Array c_integral(species.size());
	Array inflow(species.size(), 0);

	using namespace veccal;

	//the reagants need to be treated specially
	veccpy(rs.fss[0].mol_freq, inflow);

	rp.assign(reactions.size(), 0);
	rm.assign(reactions.size(), 0);

	//effective area
	for (const auto& fs : rs.fss) {
		c_integral += fs.mol_freq;
		for (size_t i = 0; i < reactions.size(); i++) {
			rp[i] += fs.reaction_freq[i].first;
			rm[i] += fs.reaction_freq[i].second;
		}
	}

	for (unsigned int j = 0; j < rs.reactions.size(); j++) {
		for (const auto& s : rs.reactions[j].products) {
			inflow[s] += rp[j];
		}
		for (const auto& s : rs.reactions[j].reagants) {
			inflow[s] += rm[j];
		}
	}

	species_life.resize(species.size());
	for (size_t i = 0; i < species.size(); i++) {
		if (inflow[i] > 0)species_life[i] = timeStep * c_integral[i] / inflow[i];
		else printf("Warning: inflow calc error: %zd\n", i);
	}
}
