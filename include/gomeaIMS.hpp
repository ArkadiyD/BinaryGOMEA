#pragma once

#include <vector>
#include <unordered_map>
using namespace std;

#include "Config.hpp"
#include "Population.hpp"
#include "problems.hpp"
#include "shared.hpp"
#include "gomea.hpp"

class gomeaIMS: public GOMEA
{
public:
	int maximumNumberOfGOMEAs;
	int generationsWithoutImprovement;
	int IMSsubgenerationFactor, basePopulationSize, numberOfGOMEAs, numberOfGenerationsIMS, minimumGOMEAIndex;

	vector<Population*> GOMEAs;

	gomeaIMS(Config *config_);
	~gomeaIMS();
	
	void initializeNewGOMEA();
	bool checkTermination();
	void generationalStepAllGOMEAs();
	bool checkTerminationGOMEA(int GOMEAIndex);
	void GOMEAGenerationalStepAllGOMEAsRecursiveFold(int indexSmallest, int indexBiggest);
	void run();
};