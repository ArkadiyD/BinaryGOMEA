#pragma once

#include <cmath>
#include <iostream> 
#include <vector>
#include <csignal>
using namespace std;

#include "Individual.hpp"
#include "Config.hpp"
#include "shared.hpp"
#include "problems.hpp"
#include "FOS.hpp"
#include "PopulationGeneral.hpp"

class Population_P3_MI: public PopulationGeneral
{
public:
	
	Population_P3_MI(Config *config_, Problem *problemInstance_, sharedInformation *sharedInformationPointer_, size_t GOMEAIndex_, size_t populationSize_):
		PopulationGeneral(config_, problemInstance_, sharedInformationPointer_, GOMEAIndex_, populationSize_){};	
	~Population_P3_MI();

	void makeOffspring();
	void generateOffspring();
	bool GOM(size_t offspringIndex, Individual *backup);
	bool conditionalGOM(size_t offspringIndex, Individual *backup, vector<vector<int> > &neighbors);	
};