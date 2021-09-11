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

class PopulationGeneral
{
public:
	Config *config;
	Problem *problemInstance;
	sharedInformation *sharedInformationPointer;
	size_t GOMEAIndex;
	size_t populationSize;

	vector<Individual*> population;
	vector<Individual*> offspringPopulation;
	vector<int> noImprovementStretches;

	FOS *populationFOS;
	bool terminated;
	double averageFitness=-1e+308;
	size_t numberOfGenerations;
	int improved_restricted = 0;
	FOS *FOSInstance = NULL;
    vector<vector<double> > matrix;
   	size_t currentPyramidLevel=0;

	PopulationGeneral(Config *config_, Problem *problemInstance_, sharedInformation *sharedInformationPointer_, size_t GOMEAIndex_, size_t populationSize_);
	virtual ~PopulationGeneral(){};

	void tournamentSelection(int k, vector<Individual*> &population, vector<Individual*> &offspringPopulation);
	void hillClimberSingle(Individual *solution);	
	void hillClimberMultiple(Individual *solution);

	void calculateAverageFitness();	
	void copyOffspringToPopulation();
	void evaluateSolution(Individual *solution);
	void evaluateSolution(Individual *solution, Individual *solutionBefore, vector<int> &touchedGenes, double fitnessBefore);
	bool GOM(size_t offspringIndex, Individual *backup);
	void findNeighbors(vector<vector< int> > &neighbors);
	bool conditionalGOM(size_t offspringIndex, Individual *backup, vector<vector<int> > &neighbors);	
	bool FI(size_t offspringIndex, Individual *backup);
	bool conditionalFI(size_t offspringIndex, Individual *backup, vector<vector<int> > &neighbors);
	void updateElitistAndCheckVTR(Individual *solution);
	void checkTimeLimit();
};

