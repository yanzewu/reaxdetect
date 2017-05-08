#include <math.h>
#include "analyzer.h"
#include "numarray.h"
#include "elementname.h"

void ReaxAnalyzer::HandleData(const ReaxReader& reader, const Simulation& simulation)
{
	for (const auto& e : simulation.atomWeights) {
		elements.push_back(WeightNameL.at((int)(e + 0.1)));
	}

	for (int i = 0; i < reader.species.size(); i++) {
		if (reader.fss[0].mol_freq[i] == 0)break;
		init[reader.species[i]] = reader.fss[0].mol_freq[i];
	}

	species = reader.species;
	reactions.resize(reader.reactions.size());
	for (size_t i = 0; i < reader.reactions.size(); i++) {
		reactions[i] = reader.reactions[i].to_string(species);
	}
	
	double interval = simulation.timeStep * (reader.fss[1].t - reader.fss[0].t);

	CalcMolLife(reader, simulation.timeStep);
	CountReaction(reader);

	if (config.sample_method == SAMPLE_FIXINT) {
		FixSample(reader, config.sample_int, config.sample_range, interval, simulation.volume);
	}
}

void ReaxAnalyzer::set_config(const Config & c)
{
	config = c;
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

void ReaxAnalyzer::CountReaction(const ReaxReader& rs)
{
	//version 2: confidence interval.
	using namespace veccal;

	//get frequency distribution
	auto size = rs.species.size();
	auto rsize = rs.reactions.size();

	sum_product_p.resize(rsize, 0);
	sum_product_m.resize(rsize, 0);
	for (auto fs = rs.fss.begin() + 1; fs < rs.fss.end(); fs++) {

		for (size_t i = 0; i < rsize; i++) {
			double product = 1.0;
			for (const auto& s : rs.reactions[i].reagants)
				product *= (double)(fs - 1)->mol_freq[s];
			sum_product_p[i] += product;

			product = 1.0;
			for (const auto& s : rs.reactions[i].products)
				product *= (double)(fs - 1)->mol_freq[s];
			sum_product_m[i] += product;

		}
	}
}

void ReaxAnalyzer::FixSample(const ReaxReader & reader, tsize_t sample_int, tsize_t range, double interval, double volume)
{
	using namespace veccal;

	for (size_t i = 0; i < reader.fss.size() - range; i += sample_int) {
		Array sum_c(reader.molecules.size(), 0), sum_r(reader.reactions.size(), 0);
		for (size_t j = 0; j < range; j++) {
			sum_c += reader.fss[i + j].mol_freq;
			sum_r += net(reader.fss[i + j].reaction_freq);
		}
		Sample sample;
		mul_into(sum_c, 1 / (double)(range), sample.c);
		mul_into(sum_r, 1 / (double)(range * interval), sample.r);
		sample.t = interval*(i + range / 2);
		samples.push_back(sample);
	}

}
