#pragma once

#include <cmath>
#include <iostream> 
#include <vector>
using namespace std;

#include "Individual.hpp"
#include "Config.hpp"
#include "shared.hpp"
#include "problems.hpp"
#include "FOS.hpp"
#include "PopulationGeneral.hpp"

class Population: public PopulationGeneral
{
public:

	Population(Config *config_, Problem *problemInstance_, sharedInformation *sharedInformationPointer_, size_t GOMEAIndex_, size_t populationSize_):
			PopulationGeneral(config_, problemInstance_, sharedInformationPointer_, GOMEAIndex_, populationSize_){};
	
	~Population();

	void makeOffspring();
	void generateOffspring();
	bool GOM(size_t offspringIndex, Individual *backup);
	bool conditionalGOM(size_t offspringIndex, Individual *backup, vector<vector<int> > &neighbors);	
};