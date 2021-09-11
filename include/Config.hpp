#pragma once

#include <string>
#include <cstdlib>
#include <iostream>
#include <random>
#include <unistd.h>
#include <getopt.h>

using namespace std;

#include "problems.hpp"
#include "FOS.hpp"

class Config
{
	void splitString(const string &str, vector<string> &splitted, char delim);
	bool isNumber(const string &s);

public:
	int usePartialEvaluations              = 0,                  
	    saveEvaluations                    = 0,
	    useForcedImprovements              = 0,
   		printHelp                          = 0;

   	double ratioVIG;
   	int conditionalGOM = 0;
   	double MI_threshold = 0.8;
   	double vtr = 1e+308;
	size_t problemIndex = 0, k = 1, s = 1,   
		FOSIndex = 0,              
	    orderFOS = 0,
	    numberOfVariables = 1;
	int similarityMeasure = 0;

	string folder = "test";
	string problemName,
		   FOSName;
	string problemInstancePath = "";

	long long timelimitMilliseconds = 1,
			  randomSeed = 42;
	
	size_t alphabetSize = 2;
	int populationScheme = 0;
	int maxArchiveSize = 1000000;
	int maximumNumberOfGOMEAs  = 100,
		IMSsubgenerationFactor = 4,
	    basePopulationSize     = 2,
	    populationSize = 1, 
	    maxGenerations = 200,
	    maxGenerationsWithoutImprovement = 10;
	int hillClimber = 0,
	    donorSearch = 0,
	    tournamentSelection = 0;

	mt19937 rng;

	bool parseCommandLine(int argc, char **argv);
	void checkOptions();
	void printUsage();
	void printOverview();
};