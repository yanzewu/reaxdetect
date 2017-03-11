#include "analyzer.h"
#include "numarray.h"
#include "elementname.h"

void ReaxAnalyzer::HandleData(const ReaxReader& reader, const Simulation& simulation)
{
	for (const auto& e : simulation.atomWeights) {
		elements.push_back(WeightName[(int)(e + 0.1)]);
	}

	species = reader.species;
	reactions.resize(reader.reactions.size());
	for (size_t i = 0; i < reader.reactions.size(); i++) {
		reactions[i] = reader.reactions[i].to_string(species);
	}
	
	CalcMolLife(reader, simulation.timeStep);
	CalculateRateConstant(reader, simulation.timeStep, simulation.volume, get_confidence(config.confidence), config.collide_prob);

	if (config.sample_method == SAMPLE_FIXINT) {
		FixSample(reader, config.sample_int, config.sample_range, simulation.timeStep, simulation.volume);
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

double ReaxAnalyzer::get_confidence(double confidence) {
	if (confidence == 0.90) return 1.645;
	else if (confidence == 0.95)return 1.96;
	else if (confidence == 0.99)return 2.576;
	else {
		printf("Warning: Confidence not valid. Using 95%% confidence.\n");
		return 1.96;
	}
}

void ReaxAnalyzer::CalculateRateConstant(const ReaxReader& rs, double timeStep, double volume, double z, double collideRatio)
{
	//version 2: confidence interval.
	using namespace veccal;

	//get frequency distribution
	auto size = rs.species.size();
	auto rsize = rs.reactions.size();

	Arrayd c(size, 0);
	Arrayd sum_product_p(rsize, 0), sum_product_m(rsize, 0);
	for (auto fs = rs.fss.begin() + 1; fs < rs.fss.end(); fs++) {

		mul_into((fs - 1)->mol_freq + fs->mol_freq, 0.5, c);

		//f_i = freq_i/\Proc(Ave(N_i))
		for (size_t i = 0; i < rsize; i++) {
			double product = 1.0;
			for (const auto& s : rs.reactions[i].reagants)
				product *= c[s];
			sum_product_p[i] += product;

			product = 1.0;
			for (const auto& s : rs.reactions[i].products)
				product *= c[s];
			sum_product_m[i] += product;

		}
	}

	kp.resize(rsize);
	km.resize(rsize);
	for (size_t i = 0; i < 2; i++) {
		kp_range[i].resize(rsize);
		km_range[i].resize(rsize);
	}

	for (size_t i = 0; i < rsize; i++) {
		if (sum_product_p[i] == 0.5)sum_product_p[i] = 1;
		if (sum_product_m[i] == 0.5)sum_product_m[i] = 1;

		double n = sum_product_p[i];
		double kp_raw = rp[i] / sum_product_p[i];
		double km_raw = rm[i] / sum_product_m[i];
		double x1 = kp_raw + z*z / (2.0*n);
		double x2 = z * sqrt(kp_raw * (1 - kp_raw) / n + z*z / (4.0 * n*n));

		kp_range[0][i] = (x1 - x2) / (1.0 + z*z / n);
		kp_range[1][i] = (x1 + x2) / (1.0 + z*z / n);

		n = sum_product_m[i];
		x1 = km_raw + z*z / (2.0*n);
		x2 = z * sqrt(km_raw * (1 - km_raw) / n + z*z / (4.0 * n*n));
		km_range[0][i] = (x1 - x2) / (1.0 + z*z / n);
		km_range[1][i] = (x1 + x2) / (1.0 + z*z / n);

		kp[i] = kp_raw * pow(volume, rs.reactions[i].reagants.size() - 1) / timeStep;
		km[i] = km_raw * pow(volume, rs.reactions[i].products.size() - 1) / timeStep;
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
		mul_into(sum_c, 1 / (volume * range), sample.c);
		mul_into(sum_r, 1 / (volume * range * interval), sample.r);
		sample.t = interval*(i + range / 2);
		samples.push_back(sample);
	}

}
