#include "PopulationGeneral.hpp"

PopulationGeneral::PopulationGeneral(Config *config_, Problem *problemInstance_, sharedInformation *sharedInformationPointer_, size_t GOMEAIndex_, size_t populationSize_): 
    config(config_), 
    problemInstance(problemInstance_),
    sharedInformationPointer(sharedInformationPointer_),
    GOMEAIndex(GOMEAIndex_), 
    populationSize(populationSize_)
{
    terminated = false;
    numberOfGenerations = 0;
    averageFitness = -1e+308;
    
    population.resize(populationSize);
    offspringPopulation.resize(populationSize);
    noImprovementStretches.resize(populationSize);
    
    vector<int> allGenes(config->numberOfVariables);
    iota(allGenes.begin(), allGenes.end(), 0);

    for (size_t i = 0; i < populationSize; ++i)
    {
      noImprovementStretches[i] = 0;

      population[i] = new Individual(config->numberOfVariables, config->alphabetSize);
      population[i]->randomInit(&config->rng);
      evaluateSolution(population[i], NULL, allGenes, 0.0);

      if (config->hillClimber == 1)
        hillClimberSingle(population[i]);
      else if (config->hillClimber == 2)
        hillClimberMultiple(population[i]);

      offspringPopulation[i] = new Individual(config->numberOfVariables, config->alphabetSize);
      *offspringPopulation[i] = *population[i];      
    }

    if (config->FOSIndex == 1 || config->FOSIndex == 2)
      createFOSInstance(config->FOSIndex, &FOSInstance, config->numberOfVariables, config->alphabetSize, config->similarityMeasure);
    else
    {
      if (sharedInformationPointer->FOSes.empty() || sharedInformationPointer->FOSes.size()-1 < currentPyramidLevel)
      {
        createFOSInstance(config->FOSIndex, &FOSInstance, config->numberOfVariables, config->alphabetSize, config->similarityMeasure);
        sharedInformationPointer->FOSes.push_back(FOSInstance);
      }
      FOSInstance = sharedInformationPointer->FOSes[currentPyramidLevel];
    }

    calculateAverageFitness();
}

void PopulationGeneral::calculateAverageFitness()
{
	averageFitness = 0.0;
	for (size_t i = 0; i < populationSize; ++i)
		averageFitness += population[i]->fitness;
	averageFitness /= populationSize;
}

void PopulationGeneral::copyOffspringToPopulation()
{
  for(size_t i = 0; i < populationSize; i++)
  {
  	*population[i] = *offspringPopulation[i];
  }
}

void PopulationGeneral::tournamentSelection(int k, vector<Individual*> &population, vector<Individual*> &offspringPopulation)
{
  int populationSize = population.size();

  vector<int> indices(populationSize * k);
  for (int i = 0; i < k; ++i)
  {
    for (int j = 0; j < populationSize; ++j)
      indices[populationSize*i + j] = j;

    shuffle(indices.begin() + populationSize*i, indices.begin() + populationSize*(i+1), config->rng);
  }
  for (int i = 0; i < populationSize; i++)
  {
    int winnerInd = 0;
    double winnerFitness = -1e+308;

    for (int j = 0; j < k; j++)
    {
      int challengerInd = indices[k*i+j];
      double challengerFitness = population[challengerInd]->fitness;
      //cout << i << " " << j << " " << challengerInd << endl;
      if (challengerFitness > winnerFitness)
      {
        winnerInd = challengerInd;
        winnerFitness = challengerFitness;
      }
    }

    *offspringPopulation[i] = *population[winnerInd];
  }
}

void PopulationGeneral::hillClimberSingle(Individual *solution)
{
	vector<int> positions(config->numberOfVariables);
  iota(positions.begin(), positions.end(), 0);   

  shuffle(positions.begin(), positions.end(), config->rng);

  for (int j = 0; j < positions.size(); ++j)
  {
    int curPos = positions[j];
    char curValue = solution->genotype[curPos];

  	for (char k = 0; k < config->alphabetSize; ++k)
  	{
    	if (k == curValue)
      		continue;

    	Individual backup = *solution;  
    	vector<int> touchedGenes(1, curPos);

    	solution->genotype[curPos] = k;

    	evaluateSolution(solution, &backup, touchedGenes, backup.fitness);

    	if (solution->fitness <= backup.fitness)
      		*solution = backup;
    }
  }
}

void PopulationGeneral::hillClimberMultiple(Individual *solution)
{
	vector<int> positions(config->numberOfVariables);
	iota(positions.begin(), positions.end(), 0);

	while (true)
	{
	  bool solutionImproved = false;

	  shuffle(positions.begin(), positions.end(), config->rng);

	  for (int j = 0; j < positions.size(); ++j)
	  {
	    int curPos = positions[j];
	    char curValue = solution->genotype[curPos];

	    for (char k = 0; k < config->alphabetSize; ++k)
	    {
	      if (k == curValue)
	        continue;

	      Individual backup = *solution;  
	      vector<int> touchedGenes(1, curPos);

	      solution->genotype[curPos] = k;

	      evaluateSolution(solution, &backup, touchedGenes, backup.fitness);

	      if (solution->fitness > backup.fitness)
	        solutionImproved = true;
	      else
	        *solution = backup;
	    }
	  }

	  if (!solutionImproved)
	    break;
	}
}

void PopulationGeneral::findNeighbors(vector<vector<int> > &neighbors)
{
  neighbors.resize(FOSInstance->FOSSize());   
  vector<vector<int> > neighborsMatrix(FOSInstance->FOSSize());
  vector<vector<double> > matrix;
  FOSInstance->get_MI_Matrix(matrix);

  for (size_t i = 0; i < FOSInstance->FOSSize(); i++)
  {
    neighbors[i].clear();
    neighborsMatrix[i].resize(config->numberOfVariables);
    fill(neighborsMatrix[i].begin(), neighborsMatrix[i].end(), 0);
    for (int p = 0; p < FOSInstance->FOSElementSize(i); ++p)
      neighborsMatrix[i][FOSInstance->FOSStructure[i][p]] = -1;

    if (FOSInstance->FOSElementSize(i) == 0 || FOSInstance->FOSElementSize(i) == config->numberOfVariables)
      continue;
    
    vector<pair<double, int> > score_for_variable(config->numberOfVariables);

    for (int p = 0; p < config->numberOfVariables; ++p)
    {
      score_for_variable[p].first = 0.0;
      score_for_variable[p].second = p;
      
      if (neighborsMatrix[i][p] < 0)
      {
        score_for_variable[p].first = -1;
        continue;
      }

      for (int q = 0; q < FOSInstance->FOSElementSize(i); ++q)
      {
        int var = FOSInstance->FOSStructure[i][q];
        //cout << p << " " << var << " " << matrix[p][var] << endl;
        score_for_variable[p].first += matrix[p][var];
      }
      score_for_variable[p].first /= FOSInstance->FOSElementSize(i);
      //cout << score_for_variable[p].first << " ";
    }
    sort(score_for_variable.begin(), score_for_variable.end(), [&](const pair<double, int> l, const pair<double,int> r){return l.first > r.first;});
    double top=0.0;
    for (int j = 0; j < score_for_variable.size(); ++j)
    {
      if (j == 0)
        top = score_for_variable[j].first;
      if (config->MI_threshold <= 1)
      {
        if (top <= 0.0 || score_for_variable[j].first/top < config->MI_threshold)
          break;
      }
      else
      {
        if (top <= 0.0 || j >= config->MI_threshold)
          break;      
      }

      int var = score_for_variable[j].second;

      if (neighborsMatrix[i][var] == 0)
      {
        neighborsMatrix[i][var] = 1;
        neighbors[i].push_back(var);
      }
    }   
  }
}

bool PopulationGeneral::FI(size_t offspringIndex, Individual *backup)
{
  vector<int> FOSIndices;
  FOSInstance->orderFOS(config->orderFOS, FOSIndices, &config->rng); 

  bool solutionHasChanged = 0;

  for (size_t i = 0; i < FOSInstance->FOSSize(); i++)
  {
    int ind = FOSIndices[i];

    if (FOSInstance->FOSElementSize(ind) == 0 || FOSInstance->FOSElementSize(ind) == config->numberOfVariables)
      continue;

    vector<int> touchedGenes;      
    bool donorEqualToOffspring = true;
    for(size_t j = 0; j < FOSInstance->FOSElementSize(ind); j++)
    {
      int variableFromFOS = FOSInstance->FOSStructure[ind][j];
      offspringPopulation[offspringIndex]->genotype[variableFromFOS] = sharedInformationPointer->elitist.genotype[variableFromFOS];
      touchedGenes.push_back(variableFromFOS);
      if (backup->genotype[variableFromFOS] != offspringPopulation[offspringIndex]->genotype[variableFromFOS])
        donorEqualToOffspring = false;
    }

    if (!donorEqualToOffspring)
    {
      evaluateSolution(offspringPopulation[offspringIndex], backup, touchedGenes, backup->fitness);

      if (offspringPopulation[offspringIndex]->fitness > backup->fitness)
      {
        *backup = *offspringPopulation[offspringIndex];
        solutionHasChanged = true;
      }
      else
      {
        *offspringPopulation[offspringIndex] = *backup;
      }
    }
    if (solutionHasChanged)
      break;
  }

  if (!solutionHasChanged)
  {
    *offspringPopulation[offspringIndex] = sharedInformationPointer->elitist;
  }

  return solutionHasChanged;
}

bool PopulationGeneral::conditionalFI(size_t offspringIndex, Individual *backup, vector<vector<int> > &neighbors)
{
  vector<bool> sampled(config->numberOfVariables, false);
  bool solutionHasChanged = 0;

  vector <bool> dependentMonitor(config->numberOfVariables, 0);

  vector<int> FOSIndices;
  FOSInstance->orderFOS(config->orderFOS, FOSIndices, &config->rng); 

  int already_sampled = 0;

  for (size_t i = 0; i < FOSInstance->FOSSize(); i++)
  {
    int ind = FOSIndices[i];

    if (FOSInstance->FOSElementSize(ind) == 0 || FOSInstance->FOSElementSize(ind) == config->numberOfVariables)
      continue;

    fill(dependentMonitor.begin(), dependentMonitor.end(), 0);
    vector<int> dependent;
    
    for(size_t j = 0; j < neighbors[ind].size(); j++)
    {
      int neighbor = neighbors[ind][j];
      
      if (sampled[neighbor])
      {
        if (dependentMonitor[neighbor] == 0)
          dependent.push_back(neighbor); 
        dependentMonitor[neighbor]=1; 
      }
    }

    for(size_t j = 0; j < FOSInstance->FOSElementSize(ind); j++)
    {
      int variableFromFOS = FOSInstance->FOSStructure[ind][j];    
      if (sampled[variableFromFOS] == 0)
        already_sampled++; 
      sampled[variableFromFOS]=1;  
    }
    
    bool canBeChosen = true;
    for (int k = 0; k < dependent.size(); ++k)
    {
      if (offspringPopulation[offspringIndex]->genotype[dependent[k]] != sharedInformationPointer->elitist.genotype[dependent[k]])
      {
        canBeChosen = false;
        break;
      }
    }
    if (!canBeChosen)
      continue;

    vector<int> touchedGenes;     
    bool donorEqualToOffspring = true; 
    for(size_t j = 0; j < FOSInstance->FOSElementSize(ind); j++)
    {
      int variableFromFOS = FOSInstance->FOSStructure[ind][j];
      offspringPopulation[offspringIndex]->genotype[variableFromFOS] = sharedInformationPointer->elitist.genotype[variableFromFOS];
      touchedGenes.push_back(variableFromFOS);
      if (backup->genotype[variableFromFOS] != offspringPopulation[offspringIndex]->genotype[variableFromFOS])
        donorEqualToOffspring = false;
    }

    if (!donorEqualToOffspring)
    {
      evaluateSolution(offspringPopulation[offspringIndex], backup, touchedGenes, backup->fitness);

      if (offspringPopulation[offspringIndex]->fitness > backup->fitness)
      {
        *backup = *offspringPopulation[offspringIndex];
        solutionHasChanged = true;
      }
      else
      {
        *offspringPopulation[offspringIndex] = *backup;
      }
    }
    if (solutionHasChanged)
      break;

  }

  if (!solutionHasChanged)
  {
    *offspringPopulation[offspringIndex] = sharedInformationPointer->elitist;
  }

  return solutionHasChanged;
}


void PopulationGeneral::evaluateSolution(Individual *solution, Individual *solutionBefore, vector<int> &touchedGenes, double fitnessBefore)
{  
  checkTimeLimit();

  /* Do the actual evaluation */
  archiveRecord searchResult;
  
  if (config->saveEvaluations)
    sharedInformationPointer->evaluatedSolutions->checkAlreadyEvaluated(solution->genotype, &searchResult);
  
  if (searchResult.isFound)
    solution->fitness = searchResult.value;
  else
  { 
    //cout << "before eval" << solution -> fitness << endl;
    if (config->usePartialEvaluations && solutionBefore != NULL)
    {
      problemInstance->calculateFitnessPartialEvaluations(solution, solutionBefore, touchedGenes, fitnessBefore);
      sharedInformationPointer->numberOfEvaluations += (double)touchedGenes.size() / config->numberOfVariables;
    }
    else
    {
      problemInstance->calculateFitness(solution);
      sharedInformationPointer->numberOfEvaluations += 1;
    }

    if (config->saveEvaluations)
      sharedInformationPointer->evaluatedSolutions->insertSolution(solution->genotype, solution->fitness);
  }

  updateElitistAndCheckVTR(solution);
}

void PopulationGeneral::checkTimeLimit()
{
  if (getMilliSecondsRunningSinceTimeStamp(sharedInformationPointer->startTimeMilliseconds) > config->timelimitMilliseconds)
  {
    cout << "TIME LIMIT REACHED!" << endl;
    throw customException("time");
  }
}

void PopulationGeneral::updateElitistAndCheckVTR(Individual *solution)
{
  /* Update elitist solution */
  if (sharedInformationPointer->firstEvaluationEver || (solution->fitness > sharedInformationPointer->elitist.fitness))
  {
    sharedInformationPointer->elitistSolutionHittingTimeMilliseconds = getMilliSecondsRunningSinceTimeStamp(sharedInformationPointer->startTimeMilliseconds);
    sharedInformationPointer->elitistSolutionHittingTimeEvaluations = sharedInformationPointer->numberOfEvaluations;

    sharedInformationPointer->elitist = *solution;
    //cout << ":" << sharedInformationPointer->elitist.fitness << endl;

    /* Check the VTR */
    if (solution->fitness >= config->vtr)
    {
      writeElitistSolutionToFile(config->folder, sharedInformationPointer->elitistSolutionHittingTimeEvaluations, sharedInformationPointer->elitistSolutionHittingTimeMilliseconds, solution);
      cout << "VTR HIT!\n";
      throw customException("vtr");
    }
  
    writeElitistSolutionToFile(config->folder, sharedInformationPointer->elitistSolutionHittingTimeEvaluations, sharedInformationPointer->elitistSolutionHittingTimeMilliseconds, solution);
  }

  sharedInformationPointer->firstEvaluationEver = false;
}
