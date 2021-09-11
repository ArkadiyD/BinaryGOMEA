#pragma once

#include <vector>
#include <unordered_map>
#include <string>
#include <fstream>
#include <deque>
#include <algorithm>
#include <unordered_set>
#include <set>
#include <cassert>
using namespace std;

#include "Individual.hpp"
#include "utils.hpp"

class Config;
#include "Config.hpp"

class Problem
{
public:
	int numberOfVariables;
	int usePartialEvaluations;
	
	Problem(){};
	virtual ~Problem(){};
	virtual void initializeProblem(int numberOfVariables)=0;
	virtual double calculateFitness(Individual *solution)=0;
	virtual double calculateFitnessPartialEvaluations(Individual *solution, Individual *solutionBefore, vector<int> &touchedGenes, double fitnessBefore);
};

class oneMax:public Problem
{
public:
	oneMax(){cout<<"creating oneMax\n";}
	void initializeProblem(int numberOfVariables_)
	{
		numberOfVariables = numberOfVariables_;
	};
	double calculateFitness(Individual *solution);
	double calculateFitnessPartialEvaluations(Individual *solution, Individual *solutionBefore, vector<int> &touchedGenes, double fitnessBefore);
};

class concatenatedDeceptiveTrap:public Problem
{
	int k, s;
	bool bimodal;
	vector<vector<int> > trapsForVariable; //determines traps which each variables belongs to
public:
	concatenatedDeceptiveTrap(int k_, int s_, bool bimodal_): k(k_), s(s_), bimodal(bimodal_)
	{
		if (not bimodal_)
			cout<<"creating concatenated Deceptive Trap with trap size=" << k << " and shift=" << s << endl;
		else
		{
			if (k != 10 && k != 6)
			{
				cout << "Bimodal trap with k=" << k << " not implemented!" << endl;
				exit(0);
			}
			cout<<"creating bimodal concatenated Deceptive Trap with trap size=" << k << " and shift=" << s << endl;
		}
	}
	void initializeProblem(int numberOfVariables_);
	double calculateFitness(Individual *solution);
	double calculateFitnessPartialEvaluations(Individual *solution, Individual *solutionBefore, vector<int> &touchedGenes, double fitnessBefore);
};

struct NKSubfunction
{
	vector<int> variablesPositions;
	vector<double> valuesTable;
};

class ADF:public Problem
{
	string problemInstanceFilename;
	vector<NKSubfunction> subfunctions;
	vector<vector<int> > subfunctionsForVariable;

public:
	ADF(string problemInstanceFilename_): problemInstanceFilename(problemInstanceFilename_)
	{
		cout<<"creating ADF\n";
	}

	void initializeProblem(int numberOfVariables_);

	double calculateFitness(Individual *solution);
	double calculateFitnessPartialEvaluations(Individual *solution, Individual *solutionBefore, vector<int> &touchedGenes, double fitnessBefore);
};

class hierarchialDeceptiveTrap:public Problem
{
	int k;
	vector<int> transform;
public:
	hierarchialDeceptiveTrap(int k_): k(k_)
	{
		cout<<"creating hierarchialDeceptiveTrap\n";
	}

	void initializeProblem(int numberOfVariables_)
	{
		numberOfVariables = numberOfVariables_;
		if (!isPowerOfK(numberOfVariables, k))
		{
			cerr << "Number of bits should be a power of k! " << numberOfVariables << " is not a power of " << k << endl;
			exit(0);
		}
		transform.resize(numberOfVariables);		
	};

	double generalTrap(int unitation, double leftPeak, double rightPeak);
	double calculateFitness(Individual *solution);
	//double calculateFitnessPartialEvaluations(Individual *solution, Individual *solutionBefore, vector<double> &touchedGenes, double fitnessBefore);
};

class hierarchialIfAndOnlyIf:public Problem
{
public:
	hierarchialIfAndOnlyIf()
	{
		cout<<"creating hierarchialIfAndOnlyIf\n";
	}
	void initializeProblem(int numberOfVariables_)
	{
		numberOfVariables = numberOfVariables_;
		if (!isPowerOfK(numberOfVariables, 2))
		{
			cerr << "Number of bits should be a power of 2! " << numberOfVariables<< " is not a power of 2" << endl;
			exit(0);
		}
	};

	double calculateFitness(Individual *solution);
};

class maxCut:public Problem
{
	string problemInstanceFilename;
	vector<pair<pair<int, int>, double > > edges;
	vector<vector<int> > edgesForVariable;

public:
	maxCut(string problemInstanceFilename_): problemInstanceFilename(problemInstanceFilename_)
	{
		cout<<"creating maxCut\n";
	}
	void initializeProblem(int numberOfVariables_);
	double calculateFitness(Individual *solution);
	double calculateFitnessPartialEvaluations(Individual *solution, Individual *solutionBefore, vector<int> &touchedGenes, double fitnessBefore);
};

class leadingOnes:public Problem
{
public:
	leadingOnes()
	{
		cout<<"creating leadingOnes\n";
	}
	void initializeProblem(int numberOfVariables_)
	{
		numberOfVariables = numberOfVariables_;
	};

	double calculateFitness(Individual *solution);
	//double calculateFitnessPartialEvaluations(Individual *solution, Individual *solutionBefore, vector<double> &touchedGenes, double fitnessBefore);
};

class Clustering:public Problem
{
	string problemInstanceFilename;
	vector<vector<double> > points;
	vector<vector<double> > distances;
	
	int Dim;
public:
	Clustering(string problemInstanceFilename_): problemInstanceFilename(problemInstanceFilename_)
	{
		cout<<"creating Clustering\n";
	}
	void initializeProblem(int numberOfVariables_);
	double calculateFitness(Individual *solution);
};

class MAXSAT:public Problem
{
	string problemInstanceFilename;
	vector<vector<int> > subfunctions;
	vector<vector<int> > signs;
	
	vector<vector<int> > subfunctionsForVariable;

public:
	MAXSAT(string problemInstanceFilename_): problemInstanceFilename(problemInstanceFilename_)
	{
		cout<<"creating MAXSAT\n";
	}

	void initializeProblem(int numberOfVariables_);

	double calculateFitness(Individual *solution);
};

class SpinGlass:public Problem
{
	string problemInstanceFilename;
	vector<vector<int> > subfunctions;
	
	vector<vector<int> > subfunctionsForVariable;

public:
	SpinGlass(string problemInstanceFilename_): problemInstanceFilename(problemInstanceFilename_)
	{
		cout<<"creating SpinGlass\n";
	}

	void initializeProblem(int numberOfVariables_);	

	double calculateFitness(Individual *solution);
};

double deceptiveTrap(int unitation, int k);
double bimodalDeceptiveTrap(int unitation, int k);

void createProblemInstance(int problemIndex, int numberOfVariables, Config *config, Problem **problemInstance, string &instancePath, int k = 1, int s = 1);
bool problemNameByIndex(Config *config, string &problemName);

