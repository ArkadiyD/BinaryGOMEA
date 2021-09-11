#pragma once

#include <vector>
#include <string>
#include <algorithm>
#include <fstream>
#include <cmath>
#include <deque>
#include <random>
#include <bits/stdc++.h> 
#include <cassert>
using namespace std;

#include "Individual.hpp"
#include "utils.hpp"


class FOS
{	
public:
	vector<vector<int> > FOSStructure;	
	size_t numberOfVariables;
	size_t alphabetSize;
	
	vector<int> improvementCounters;
	vector<int> usageCounters;

	vector<int> single_counters;
	vector<vector<int> > pair_counters;
	
	FOS(size_t numberOfVariables_, size_t alphabetSize_): numberOfVariables(numberOfVariables_), alphabetSize(alphabetSize_)
	{}

	virtual ~FOS(){};

	size_t FOSSize()
	{
		return FOSStructure.size();
	}

	size_t FOSElementSize(int i)
	{
		return FOSStructure[i].size();
	}
	
	virtual void learnFOS(vector<Individual*> &population, vector<vector<int> > *VIG = NULL, mt19937 *rng = NULL) = 0;
	void writeToFileFOS(string folder, int populationIndex, int generation);
	void writeFOSStatistics(string folder, int populationIndex, int generation);
	void setCountersToZero();
	void shuffleFOS(vector<int> &indices, mt19937 *rng);
	void sortFOSAscendingOrder(vector<int> &indices);
	void sortFOSDescendingOrder(vector<int> &indices);
	void orderFOS(int orderingType, vector<int> &indices, mt19937 *rng);

	virtual void get_MI_Matrix(vector<vector<double> > &matrix){};
	virtual void addSolution(Individual *solution){};
	virtual void addSolutionTournamentSelection(Individual solution){};
	
	virtual void writeMIMatrixToFile(string folder, int populationIndex, int generation){};
};


class LTFOS: public FOS
{
private:
	vector<vector<double> > MI_Matrix;
	vector<vector<double> > S_Matrix;
	bool filtered;
	int similarityMeasure;
	int determineNearestNeighbour(int index, vector< vector< int > > &mpm);
	void computeMIMatrix(vector<Individual*> &population);
	void computeNMIMatrix(vector<Individual*> &population);
	void estimateParametersForSingleBinaryMarginal(vector<Individual*> &population, vector<size_t> &indices, size_t  &factorSize, vector<double> &result);

public:	
	LTFOS(size_t numberOfVariables_, size_t alphabetSize_, int similarityMeasure, bool filtered=false);
	~LTFOS(){};

	void learnFOS(vector<Individual*> &population, vector<vector<int> > *VIG = NULL, mt19937 *rng = NULL);
	void writeMIMatrixToFile(string folder, int populationIndex, int generation);
	void get_MI_Matrix(vector<vector<double> > &matrix)
	{
		matrix.resize(MI_Matrix.size());
		for (int i = 0; i < numberOfVariables; ++i)
		{
			matrix[i].resize(MI_Matrix[i].size());
			for (int j = 0; j < numberOfVariables; ++j)
			{
				matrix[i][j] = MI_Matrix[i][j];
				//cout << matrix[i][j] << " ";
			}
		}
	}

};

class LTFOSEfficient: public FOS
{
private:
	vector<vector<int> > single_counters;
	vector<vector<vector<int> > > pair_counters;

	vector<vector<double> > MI_Matrix;
	vector<vector<double> > S_Matrix;
	bool filtered;
	int similarityMeasure;
	int determineNearestNeighbour(int index, vector< vector< int > > &mpm);
	void computeMIMatrix(vector<Individual*> &population);
	void computeNMIMatrix(vector<Individual*> &population);
	void estimateParametersForSingleBinaryMarginal(vector<Individual*> &population, vector<size_t> &indices, size_t  &factorSize, vector<double> &result);
	vector<Individual> minHeap;
	vector<Individual> maxHeap;
	
public:	
	LTFOSEfficient(size_t numberOfVariables_, size_t alphabetSize_, int similarityMeasure, bool filtered=false);
	~LTFOSEfficient(){};

	void learnFOS(vector<Individual*> &population, vector<vector<int> > *VIG = NULL, mt19937 *rng = NULL);
	void buildGraph(double thresholdValue, mt19937 *rng);
	void writeMIMatrixToFile(string folder, int populationIndex, int generation);
	void addSolution(Individual *solution);
	void addSolutionTournamentSelection(Individual solution);
	
	void get_MI_Matrix(vector<vector<double> > &matrix)
	{
		matrix.resize(MI_Matrix.size());
		for (int i = 0; i < numberOfVariables; ++i)
		{
			matrix[i].resize(MI_Matrix[i].size());
			for (int j = 0; j < numberOfVariables; ++j)
			{
				matrix[i][j] = MI_Matrix[i][j];
				//cout << matrix[i][j] << " ";
			}
		}
	}

};


struct ascendingFitness{
  bool operator()(const Individual &a, const Individual &b) const{
    return a.fitness < b.fitness;
  }
};

struct descendingFitness{
  bool operator()(const Individual &a, const Individual &b) const{
    return a.fitness > b.fitness;
  }
};

bool FOSNameByIndex(size_t FOSIndex, string &FOSName);
void createFOSInstance(size_t FOSIndex, FOS **FOSInstance, size_t numberOfVariables, size_t alphabetSize, int similarityMeasure);
