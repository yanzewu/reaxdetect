#pragma once

#ifndef ANALYZER_H
#define ANALYZER_H

#include "simulation.h"
#include "reaxreader.h"

#define SAMPLE_FIXINT	75
#define SAMPLE_SMART	77

// a middle class used in output of reactionstat data.
class ReaxAnalyzer {
public:

	struct Sample {
		double t;
		Arrayd c;
		Arrayd r;
	};
	vector<Sample> samples;

	vector<string> elements;
	vector<string> species;
	vector<string> reactions;

	Arrayd species_life;
	Arrayd kp, km;
	Arrayd kp_range[2], km_range[2];
	Array rp, rm;

	struct Config {
		double confidence;
		double collide_prob;
		char sample_method;
		tsize_t sample_int;
		tsize_t sample_range;
		Config() : confidence(0.90), collide_prob(1.0), sample_method(SAMPLE_FIXINT) {
		}
	};

	ReaxAnalyzer() : config() {
	}
	ReaxAnalyzer(const Config& c) : config(c) {
	}

	void HandleData(const ReaxReader& rs, const Simulation&);
	void set_config(const Config& c);

private:

	void CalcMolLife(const ReaxReader& rs, double timeStep);
	void CalculateRateConstant(const ReaxReader& rs, double timeStep, double volume, double confidence, double collideRatio);
	void FixSample(const ReaxReader&, tsize_t sample_int, tsize_t range, double interval, double volume);

	static double get_confidence(double confidence);

	Config config;
};

#endif // !ANALYZER_H
