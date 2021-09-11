#pragma once

#include <vector>
#include <unordered_map>
using namespace std;

#include "utils.hpp"
#include "time.hpp"
#include "FOS.hpp"

struct sharedInformation
{
	double numberOfEvaluations;
	long long startTimeMilliseconds;
	double elitistSolutionHittingTimeMilliseconds,
	       elitistSolutionHittingTimeEvaluations;

	Individual elitist;
	solutionsArchive *evaluatedSolutions;
	bool firstEvaluationEver;

	Pyramid *pyramid;
	vector<FOS* > FOSes;
	
	sharedInformation(int maxArchiveSize)
	{
		numberOfEvaluations = 0;
		startTimeMilliseconds = getCurrentTimeStampInMilliSeconds();
		firstEvaluationEver = true;
		evaluatedSolutions = new solutionsArchive(maxArchiveSize);
		pyramid = new Pyramid();
	}

	~sharedInformation()
	{
		delete evaluatedSolutions;
		delete pyramid;
	}
};