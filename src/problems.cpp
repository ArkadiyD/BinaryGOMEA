#include "problems.hpp"

void createProblemInstance(int problemIndex, int numberOfVariables, Config *config, Problem **problemInstance, string &instancePath, int k, int s)
{
	switch (problemIndex)
	{
		case 0: *problemInstance = new oneMax(); break;
		case 1: *problemInstance = new concatenatedDeceptiveTrap(config->k, config->s, false); break;
		case 2: *problemInstance = new ADF(instancePath); break;
		case 3: *problemInstance = new maxCut(instancePath); break;		
		case 4: *problemInstance = new hierarchialIfAndOnlyIf(); break;
		case 5: *problemInstance = new leadingOnes(); break;
		case 6: *problemInstance = new hierarchialDeceptiveTrap(3); break;
		case 7: *problemInstance = new concatenatedDeceptiveTrap(config->k, config->s, true); break;
		case 8: *problemInstance = new MAXSAT(instancePath); break;
		case 9: *problemInstance = new SpinGlass(instancePath); break;

		default:
		{
			cerr << "No problem with index #" << problemIndex << " installed!\n";
			exit(0);
		};
	}
	(*problemInstance)->initializeProblem(numberOfVariables);
}

bool problemNameByIndex(Config *config, string &problemName)
{
	switch (config->problemIndex)
	{
		case 0: problemName = "OneMax"; break;
		case 1: problemName = "Concatenated Deceptive Trap with k=" + to_string(config->k) + " s=" + to_string(config->s); break;
		case 2: problemName = "ADF with k=" + to_string(config->k) + " s=" + to_string(config->s); break;
		case 3: problemName = "MaxCut"; break;
		case 4: problemName = "Hierarhical If-And-Only-If"; break;
		case 5: problemName = "Leading Ones"; break;		
		case 6: problemName = "Hierarhical Deceptive Trap-3"; break;
		case 7: problemName = "Concatenated Bimodal Deceptive Trap with k=" + to_string(config->k) + " s=" + to_string(config->s); break;
		case 8: problemName = "MAXSAT"; break;
		case 9: problemName = "SpinGlass"; break;

		default: return false; break;
	}
	return true;
}

double Problem::calculateFitnessPartialEvaluations(Individual *solution, Individual *solutionBefore, vector<int> &touchedGenes, double fitnessBefore)
{
	return 0.0;
}

double oneMax::calculateFitness(Individual *solution)
{
	double res = 0.0;
	for (size_t i = 0; i < solution->genotype.size(); ++i)
		res += solution->genotype[i];

	solution->fitness = res;
	return res;
}

double oneMax::calculateFitnessPartialEvaluations(Individual *solution, Individual *solutionBefore, vector<int> &touchedGenes, double fitnessBefore)
{
	double res = fitnessBefore;
	for (size_t i = 0; i < touchedGenes.size(); ++i)
	{
		int touchedGene = touchedGenes[i];
		res -= solutionBefore->genotype[touchedGene];
		res += solution->genotype[touchedGene];
	}
	solution->fitness = res;
	return res;
}

double deceptiveTrap(int unitation, int k)
{
	double result = 0.0;
	if (unitation != k)
		result = k - 1 - unitation;
	else
		result = unitation;
	
	return (double)result;
}

double bimodalDeceptiveTrap(int unitation, int k)
{
	double result = 0.0;
	if (k == 10)
	{
		if (unitation == 0 || unitation == 10)
			result = 10;
		else if (unitation == 5)
			result = 9;
		else if (unitation == 1 || unitation == 9)
			result = 5;
		else if (unitation == 2 || unitation == 8)
			result = 0;
		else if (unitation == 3 || unitation == 7)
			result = 3;
		else if (unitation == 4 || unitation == 6)
			result = 6;
	}
	else if (k == 6)
	{
		if (unitation == 0 || unitation == 6)
			result = 6;
		else if (unitation == 1 || unitation == 5)
			result = 0;
		else if (unitation == 2 || unitation == 4)
			result = 2;
		else if (unitation == 3)
			result = 5;
	}	
	else
	{
		cout << "Bimodal Trap with size " << k << " not implemented!";
		exit(0);
	}

	return (double)result;
}

void concatenatedDeceptiveTrap::initializeProblem(int numberOfVariables_)
{
	numberOfVariables = numberOfVariables_;
	trapsForVariable.resize(numberOfVariables);

	for (int i = 0; i < numberOfVariables; i += s)
	{
		for (int j = i; j < i + k; j++)			
			trapsForVariable[j % numberOfVariables].push_back(i);
	}
};
	
double concatenatedDeceptiveTrap::calculateFitness(Individual *solution)
{
	double res = 0.0;
	for (int i = 0; i < numberOfVariables; i += s)
	{
		//cout << i << s << endl;
		int unitation = 0;
		for (int j = i; j < i + k; j++)
			unitation += solution->genotype[j%numberOfVariables];

		if (!bimodal)
			res += deceptiveTrap(unitation, k);
		else if (bimodal)
			res += bimodalDeceptiveTrap(unitation, k);
	}

	solution->fitness = res;
	return res;
}

double concatenatedDeceptiveTrap::calculateFitnessPartialEvaluations(Individual *solution, Individual *solutionBefore, vector<int> &touchedGenes, double fitnessBefore)
{
	double res = fitnessBefore;

	unordered_set<int> touchedTraps; //starting positions of all touched traps
	for (int i = 0; i < touchedGenes.size(); i++)
	{
		int touchedGene = touchedGenes[i];
		for (int j = 0; j < trapsForVariable[touchedGene].size(); ++j)
		{
			touchedTraps.insert(trapsForVariable[touchedGene][j]);
		}
	}

	for (unordered_set<int>::iterator it = touchedTraps.begin(); it != touchedTraps.end(); ++it)
	{
		int unitation = 0, unitationBefore = 0;

		for (int j = *it; j < *it + k; j++)
		{
			unitationBefore += solutionBefore->genotype[j];
			unitation += solution->genotype[j];
		}
		res -= deceptiveTrap(unitationBefore, k);
		res += deceptiveTrap(unitation, k);
	}

	solution->fitness = res;
	return res;
}


void ADF::initializeProblem(int numberOfVariables_)
{
	numberOfVariables = numberOfVariables_;
	subfunctionsForVariable.resize(numberOfVariables);

	ifstream inFile(problemInstanceFilename, ifstream::in);
	if (inFile.fail())
	{
		cout << "Problem Instance File " << problemInstanceFilename << " does not exist!\n";
		exit(0);
	}
	int N, numFunctions;
	inFile >> N >> numFunctions;
	//cout << N << " " << numFunctions << endl;

	for (int i = 0; i < numFunctions; ++i)
	{
		NKSubfunction subfunction;

		int k;
		inFile >> k;

		for (int j = 0; j < k; ++j)
		{
			int var;
			inFile >> var;
			subfunction.variablesPositions.push_back(var);
		}
		int numCombinations = 1 << k;
		subfunction.valuesTable.resize(numCombinations);
		for (int j = 0; j < numCombinations; ++j)
		{
			string combination;
			inFile >> combination;
			int index = 0, pow = 1;
			for (size_t p = combination.size()-1; p > 0; p--)
			{
				if (combination[p] == '0' || combination[p] == '1')
				{
					index += ((int)(combination[p] == '1') * pow);
					pow *= 2;
				}
			}
			inFile >> subfunction.valuesTable[index];			

		}
		subfunctions.push_back(subfunction);

		//needed for partial evaluations
		for (size_t p = 0; p < subfunction.variablesPositions.size(); ++p)
		{
			subfunctionsForVariable[subfunction.variablesPositions[p]].push_back(subfunctions.size()-1);
		}

		//VIG.push_back(subfunction.variablesPositions);

	}
	inFile.close();
};


double ADF::calculateFitness(Individual *solution)
{
	long double res = 0.0;
	
	for (size_t i = 0; i < subfunctions.size(); ++i)
	{
		int index = 0, pow = 1;
		for (int j = subfunctions[i].variablesPositions.size()-1; j >= 0; --j)
		{
			int pos = subfunctions[i].variablesPositions[j];
			index += solution->genotype[pos] * pow;
			pow *= 2;
		}
		res += subfunctions[i].valuesTable[index];// * multiplier;
	}
	
	solution->fitness = (double)res;
	return (double)res;
}

double ADF::calculateFitnessPartialEvaluations(Individual *solution, Individual *solutionBefore, vector<int> &touchedGenes, double fitnessBefore)
{
	double res = fitnessBefore;

	unordered_set<int> touchedSubfunctions;
	for (int i = 0; i < touchedGenes.size(); i++)
	{
		int touchedGene = touchedGenes[i];
		for (int j = 0; j < subfunctionsForVariable[touchedGene].size(); ++j)
		{
			touchedSubfunctions.insert(subfunctionsForVariable[touchedGene][j]);
		}
	}

	for (unordered_set<int>::iterator it = touchedSubfunctions.begin(); it != touchedSubfunctions.end(); ++it)
	{
		int index = 0, indexBefore = 0, pow = 1;
		for (int j = subfunctions[*it].variablesPositions.size()-1; j >= 0; --j)
		{
			int pos = subfunctions[*it].variablesPositions[j];
			index += solution->genotype[pos] * pow;
			indexBefore += solutionBefore->genotype[pos] * pow;			
			pow *= 2;
		}
		
		res -= subfunctions[*it].valuesTable[indexBefore];		
		res += subfunctions[*it].valuesTable[index];

	}
	//cout << res << endl;

	solution->fitness = res;
	return res;
}

double hierarchialDeceptiveTrap::generalTrap(int unitation, double leftPeak, double rightPeak)
{
	if (unitation == -1)
		return 0; 

	if (unitation == k)
		return rightPeak;
  	
	return leftPeak * (1.0 - unitation / (k - 1.0));
};

double hierarchialDeceptiveTrap::calculateFitness(Individual *solution)
{
	double res = 0.0;
	int blockSize = k;
	double leftPeak, rightPeak;
	for (int i = 0; i < numberOfVariables; ++i)
		transform[i] = (int)solution->genotype[i];

	while (blockSize <= numberOfVariables)
	{
		for (int i = 0; i < numberOfVariables; i += blockSize)
		{
			int unitation = 0;
			for (int j = i; j < i + blockSize; j += blockSize / k)
			{
				if (transform[j] != 0 && transform[j] != 1)
				{
					unitation = -1;
					break;
				}
				unitation += transform[j];
			}
			
			if (unitation == 0)
				transform[i] = 0;
			else if (unitation == k)
				transform[i] = 1;
			else
				transform[i] = -1;

			if (blockSize < numberOfVariables)
			{
				leftPeak = 1.0;
				rightPeak = 1.0;
			}
			else
			{
				leftPeak = 0.9;
				rightPeak = 1.0;
			}
			res += generalTrap(unitation, leftPeak, rightPeak) * (blockSize / 3);
		}
		blockSize *= k;
	}
	solution->fitness = res;
	return res;
}


double hierarchialIfAndOnlyIf::calculateFitness(Individual *solution)
{
	double res = 0.0;
	int blockSize = 1;
	while (blockSize <= numberOfVariables)
	{
		for (int i = 0; i < numberOfVariables; i+=blockSize)
		{
			int unitation = 0;
			for (int j = i; j < i + blockSize; ++j)
				unitation += solution->genotype[j];

			if (unitation == blockSize || unitation == 0)
				res += blockSize;
			//cout << blockSize << "  " << res << endl;
		}
		blockSize *= 2;
	}
	solution->fitness = res;
	return res;
}

double leadingOnes::calculateFitness(Individual *solution)
{
	double res = 0.0;
	for (size_t i = 0; i < numberOfVariables; ++i)
	{
		if (solution->genotype[i] == 0)
			break;
		res += 1;
	}
	solution->fitness = res;
	return res;
}

void maxCut::initializeProblem(int numberOfVariables_)
{
	numberOfVariables = numberOfVariables_;
	edgesForVariable.resize(numberOfVariables);

	ifstream inFile(problemInstanceFilename, ifstream::in);
	if (inFile.fail())
	{
		cout << "Problem Instance File " << problemInstanceFilename << " does not exist!\n";
		exit(0);
	}
	int N, numEdges;
	inFile >> N >> numEdges;

	//cout << problemInstanceFilename << " " << N << "  " << numEdges << endl;

	for (int i = 0; i < numEdges; ++i)
	{
		int v1, v2;
		double w;
		inFile >> v1 >> v2 >> w;
		edges.push_back(make_pair(make_pair(v1-1, v2-1), w));
		//cout << v1 << " " << v2 << " " << w << endl;
		edgesForVariable[v1-1].push_back(i);
		edgesForVariable[v2-1].push_back(i);		
	}

	inFile.close();
}

double maxCut::calculateFitness(Individual *solution)
{
	double res = 0.0;
	for (size_t i = 0; i < edges.size(); ++i)
	{
		int v1 = edges[i].first.first;
		int v2 = edges[i].first.second;
		double w = edges[i].second;

		if (solution->genotype[v1] != solution->genotype[v2])
			res += w;
	}

	solution->fitness = res;
	return res;
}

double maxCut::calculateFitnessPartialEvaluations(Individual *solution, Individual *solutionBefore, vector<int> &touchedGenes, double fitnessBefore)
{
	double res = fitnessBefore;

	unordered_set<int> touchedEdges;
	for (int i = 0; i < touchedGenes.size(); i++)
	{
		int touchedGene = touchedGenes[i];
		for (int j = 0; j < edgesForVariable[touchedGene].size(); ++j)
		{
			touchedEdges.insert(edgesForVariable[touchedGene][j]);
		}
	}

	for (unordered_set<int>::iterator it = touchedEdges.begin(); it != touchedEdges.end(); ++it)
	{
		int v1 = edges[*it].first.first;
		int v2 = edges[*it].first.second;
		double w = edges[*it].second;

		if (solutionBefore->genotype[v1] != solutionBefore->genotype[v2])
			res -= w;

		if (solution->genotype[v1] != solution->genotype[v2])
			res += w;
	}

	solution->fitness = res;
	return res;
}

void Clustering::initializeProblem(int numberOfVariables_)
{
	numberOfVariables = numberOfVariables_;
	
	ifstream inFile(problemInstanceFilename, ifstream::in);
	if (inFile.fail())
	{
		cout << "Problem Instance File " << problemInstanceFilename << " does not exist!\n";
		exit(0);
	}
	int N;
	inFile >> N >> Dim;
	points.resize(N);
	for (int i = 0; i < N; ++i)
	{
		points[i].resize(Dim);
		for (int j = 0; j < Dim; ++j)
			inFile >> points[i][j];
		//cout << points[i][0];
	}

	inFile.close();
	distances.resize(numberOfVariables);
	for (size_t i = 0; i < numberOfVariables; ++i)
	{
		distances[i].resize(numberOfVariables);

		for (size_t j = i+1; j < numberOfVariables; ++j)
		{
			double dist = 0.0;

			for (int k = 0; k < Dim; ++k)
				dist += (points[i][k] - points[j][k])*(points[i][k] - points[j][k]);
		
			dist = pow(dist, 1.0/float(Dim));

			distances[i][j] = dist;
		}
	}
}

double Clustering::calculateFitness(Individual *solution)
{
	double res = 0.0;
	double intra_cluster_dist = 0.0, inter_cluster_dist = 0.0;
	int counter1=0, counter2=0;

	for (size_t i = 0; i < numberOfVariables; ++i)
	{
		for (size_t j = i+1; j < numberOfVariables; ++j)
		{
			if (solution->genotype[i] != solution->genotype[j])
			{
				inter_cluster_dist += distances[i][j];
				counter1 += 1;
				//if (distances[i][j] < min_dist)
				//	min_dist = distances[i][j];
			}
			else
			{
				intra_cluster_dist += distances[i][j];
				counter2 += 1;
				//if (distances[i][j] > max_dist)
				//	max_dist = distances[i][j];
			}
		}
	}

	inter_cluster_dist /= counter1;
	intra_cluster_dist /= counter2;
	
	if (counter1 == 0 or counter2 == 0)
		res = -1;
	else
		res = inter_cluster_dist / intra_cluster_dist;

	solution->fitness = res;
	return res;
}

void MAXSAT::initializeProblem(int numberOfVariables_)
{
	numberOfVariables = numberOfVariables_;
	subfunctionsForVariable.resize(numberOfVariables);

	ifstream inFile(problemInstanceFilename, ifstream::in);
	if (inFile.fail())
	{
		cout << "Problem Instance File " << problemInstanceFilename << " does not exist!\n";
		exit(0);
	}
	
	string line;
	for (int i = 0; i < 8; ++i)
		getline(inFile, line);

	while(true)
	{
		vector<int> subfunction;
		vector<int> subfunctionSigns;
		int var;
		inFile >> var;
		if (inFile.fail())
			break;
		if (var < 0)
		{
			subfunctionSigns.push_back(-1);
			var++;
		}
		else if (var > 0)
		{
			var--;
			subfunctionSigns.push_back(1);
		}
		var = abs(var);
		subfunction.push_back(var);

		while (true)
		{
			inFile >> var;
			if (var == 0)
				break;

			if (var < 0)
			{
				subfunctionSigns.push_back(-1);
				var++;
			}
			else if (var > 0)
			{
				var--;
				subfunctionSigns.push_back(1);
			}

			var = abs(var);
			subfunction.push_back(var);
		}
		for (int i = 0; i < subfunction.size(); ++i)
		{
			//cout << subfunction[i] << " | " << subfunctionSigns[i] << " ";
			subfunctionsForVariable[subfunction[i]].push_back(subfunctions.size());
		}
		//cout << endl;
		subfunctions.push_back(subfunction);
		signs.push_back(subfunctionSigns);

	}
	inFile.close();
};


double MAXSAT::calculateFitness(Individual *solution)
{
	long double res = 0.0;
	
	for (size_t i = 0; i < subfunctions.size(); ++i)
	{
		bool b = false;
		for (int j = 0; j < subfunctions[i].size(); ++j)
		{
			int var = subfunctions[i][j];
			int sign = signs[i][j];
			if (sign > 0 && solution->genotype[var] == 1)
				b = true;
			else if (sign < 0 && solution->genotype[var] == 0)
				b = true;			
		}
		if (b == false)
			res -= 1;		
	}
	
	solution->fitness = (double)res;
	return (double)res;
}



void SpinGlass::initializeProblem(int numberOfVariables_)
{
	numberOfVariables = numberOfVariables_;
	subfunctionsForVariable.resize(numberOfVariables);

	ifstream inFile(problemInstanceFilename, ifstream::in);
	if (inFile.fail())
	{
		cout << "Problem Instance File " << problemInstanceFilename << " does not exist!\n";
		exit(0);
	}
	
	double tmp1, tmp2;
	inFile >> tmp1 >> tmp2;
	
	while(true)
	{
		vector<int> subfunction;
		int var1, var2, sign;
		inFile >> var1 >> var2 >> sign;
		if (inFile.fail())
			break;
		
		var1--;
		var2--;

		subfunction.push_back(var1);
		subfunction.push_back(var2);
		subfunction.push_back(sign);		
		subfunctionsForVariable[var1].push_back(subfunctions.size());
		subfunctionsForVariable[var2].push_back(subfunctions.size());		
		subfunctions.push_back(subfunction);
		//cout << subfunction.size() << " " << var1 << " " << var2 << " " << sign << endl;
		
	}
	inFile.close();
};


double SpinGlass::calculateFitness(Individual *solution)
{
	long double res = 0.0;
	
	for (size_t i = 0; i < subfunctions.size(); ++i)
	{
		int cur_value = 1;
		for (int j = 0; j < 2; ++j)
		{
			if (solution->genotype[subfunctions[i][j]] == 1)
				cur_value *= 1;
			else
				cur_value *= -1;
		}
		cur_value *= subfunctions[i][2];
		res += cur_value;
	}
	
	res /= (double)numberOfVariables;
	solution->fitness = res;
	return (double)res;
}


