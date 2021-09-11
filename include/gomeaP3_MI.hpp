#pragma once

#include <vector>
#include <unordered_map>
using namespace std;

#include "Config.hpp"
#include "Population_P3_MI.hpp"
#include "problems.hpp"
#include "shared.hpp"
#include "gomea.hpp"

class gomeaP3_MI: public GOMEA
{
public:
	int maximumNumberOfGOMEAs;
	int basePopulationSize, numberOfGOMEAs;

	vector<Population_P3_MI*> GOMEAs;

	gomeaP3_MI(Config *config_);
	~gomeaP3_MI();
	
	void initializeNewGOMEA();
	bool checkTermination();
	void GOMEAGenerationalSteps(int GOMEAIndex);
	void run();
};