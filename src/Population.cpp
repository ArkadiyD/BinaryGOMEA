#include "Population.hpp"


Population::~Population()
{
  for (size_t i = 0; i < populationSize; ++i)
  {
    delete population[i];
    delete offspringPopulation[i];
  }

  delete FOSInstance;
}

void Population::makeOffspring()
{
  vector<Individual*> populationForModelLearning(populationSize);
  for (int i = 0; i < populationSize; ++i)
  {
    populationForModelLearning[i] = new Individual(config->numberOfVariables, config->alphabetSize);
    *populationForModelLearning[i] = *population[i];
  }

  if (config->tournamentSelection)
    tournamentSelection(2, population, populationForModelLearning); //performs tournament selection and saves the winners to population array

  FOSInstance->learnFOS(populationForModelLearning, NULL, &config->rng);
  
  for (int i = 0; i < populationSize; ++i)
    delete populationForModelLearning[i];

  FOSInstance->setCountersToZero();

  generateOffspring();
}


void Population::generateOffspring()
{
  vector<vector<int> > neighbors;
  if (config->conditionalGOM > 0)
    findNeighbors(neighbors);

  for(size_t i = 0; i < populationSize; i++)
  {
      Individual backup = *population[i];  
      
      bool solutionHasChanged;
      if (config->conditionalGOM == 0 || numberOfGenerations == 0)
        solutionHasChanged = GOM(i, &backup);
      else
        solutionHasChanged = conditionalGOM(i, &backup, neighbors);

      /* Phase 2 (Forced Improvement): optimal mixing with elitist solution */
      if (config->useForcedImprovements)
      {
        if ((!solutionHasChanged) || (noImprovementStretches[i] > (1+(log(populationSize)/log(10)))))
        {
          if (config->conditionalGOM == 0 || numberOfGenerations == 0)
            FI(i, &backup);
          else
            conditionalFI(i, &backup, neighbors);
                      
           #ifdef DEBUG_GOM
             cout << "after FI " << i << " " << *offspringPopulation[i] << endl;
           #endif    
        }    
      }

    if(!(offspringPopulation[i]->fitness > population[i]->fitness))
      noImprovementStretches[i]++;
    else
      noImprovementStretches[i] = 0;

  }
}

bool Population::GOM(size_t offspringIndex, Individual *backup)
{
  size_t donorIndex;
  bool solutionHasChanged = false;
  bool thisIsTheElitistSolution = *offspringPopulation[offspringIndex] == sharedInformationPointer->elitist;//(sharedInformationPointer->elitistSolutionGOMEAIndex == GOMEAIndex) && (sharedInformationPointer->elitistSolutionOffspringIndex == offspringIndex);
  
  *offspringPopulation[offspringIndex] = *population[offspringIndex];

  vector<int> FOSIndices;
  FOSInstance->orderFOS(config->orderFOS, FOSIndices, &config->rng); 

  vector<int> donorIndices(populationSize);
  iota(donorIndices.begin(), donorIndices.end(), 0);

  for (size_t i = 0; i < FOSInstance->FOSSize(); i++)
  {
    int ind = FOSIndices[i];

    if (FOSInstance->FOSElementSize(ind) == 0 || FOSInstance->FOSElementSize(ind) == config->numberOfVariables)
      continue;

    bool donorEqualToOffspring = true;
    int indicesTried = 0;

    while (donorEqualToOffspring && indicesTried < donorIndices.size())
    {
      int j = config->rng() % (donorIndices.size() - indicesTried);
      //cout << "rand ind:" << j << " " << indicesTried+j << " " << donorIndices.size() << endl;
      swap(donorIndices[indicesTried], donorIndices[indicesTried + j]);
      donorIndex = donorIndices[indicesTried];
      indicesTried++;

      if (donorIndex == offspringIndex)
        continue;

      vector<int> touchedGenes;
      for(size_t j = 0; j < FOSInstance->FOSElementSize(ind); j++)
      {
        int variableFromFOS = FOSInstance->FOSStructure[ind][j]; 
        offspringPopulation[offspringIndex]->genotype[variableFromFOS] = population[donorIndex]->genotype[variableFromFOS];

        touchedGenes.push_back(variableFromFOS);

        if (backup->genotype[variableFromFOS] != offspringPopulation[offspringIndex]->genotype[variableFromFOS])
          donorEqualToOffspring = false;      
      }

      if (!donorEqualToOffspring)
      {
        evaluateSolution(offspringPopulation[offspringIndex], backup, touchedGenes, backup->fitness);

        // accept the change if this solution is not the elitist and the fitness is at least equally good (allows random walk in neutral fitness landscape)
        // however, if this is the elitist solution, only accept strict improvements, to avoid convergence problems
        if ((!thisIsTheElitistSolution && (offspringPopulation[offspringIndex]->fitness >= backup->fitness)) || 
            (thisIsTheElitistSolution && (offspringPopulation[offspringIndex]->fitness > backup->fitness)))   
        {       
          *backup = *offspringPopulation[offspringIndex];
          solutionHasChanged = true;
          FOSInstance->improvementCounters[ind]++;
        }
        else
        {
          *offspringPopulation[offspringIndex] = *backup;
        }

        FOSInstance->usageCounters[ind]++;
      }

      if (!config->donorSearch) //if not exhaustive donor search then stop searching anyway
        break;
    }
  }
  //if (thisIsTheElitistSolution)
  //  cout << "after:" << *offspringPopulation[offspringIndex] << endl;
  return solutionHasChanged;
}


bool Population::conditionalGOM(size_t offspringIndex, Individual *backup, vector<vector<int> > &neighbors)
{
  size_t donorIndex;
  bool solutionHasChanged = false;
  bool thisIsTheElitistSolution = *offspringPopulation[offspringIndex] == sharedInformationPointer->elitist;//(sharedInformationPointer->elitistSolutionGOMEAIndex == GOMEAIndex) && (sharedInformationPointer->elitistSolutionOffspringIndex == offspringIndex);
  
  *offspringPopulation[offspringIndex] = *population[offspringIndex];

  #ifdef DEBUG_GOM
  cout << "before gom " << offspringIndex << " " << *offspringPopulation[offspringIndex] << " " << thisIsTheElitistSolution << endl;  
  #endif 

  int already_sampled = 0;
  vector<bool> sampled(config->numberOfVariables, false);
  vector<int> donorIndices(populationSize);
  iota(donorIndices.begin(), donorIndices.end(), 0);  
  vector<int> FOSIndices;
  FOSInstance->orderFOS(config->orderFOS, FOSIndices, &config->rng); 

  vector <bool> dependentMonitor(config->numberOfVariables, 0);

  for (size_t i = 0; i < FOSInstance->FOSSize(); i++)
  {
    int ind = FOSIndices[i];
    fill(dependentMonitor.begin(), dependentMonitor.end(), 0);

    if (FOSInstance->FOSElementSize(ind) == 0 || FOSInstance->FOSElementSize(ind) == config->numberOfVariables)
      continue;

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
    // cout << "dependent:";
    // for (int j = 0; j < dependent.size(); ++j)
    // cout << dependent[j] << " ";
    // cout << endl;
    int indicesTried = 0;
    bool donorFound = false;

    for(size_t j = 0; j < FOSInstance->FOSElementSize(ind); j++)
    {
      int variableFromFOS = FOSInstance->FOSStructure[ind][j]; 
      if (sampled[variableFromFOS] == 0)
        already_sampled++;    
      sampled[variableFromFOS]=1;  
    }

    while (indicesTried < donorIndices.size())
    {
      int j = config->rng() % (donorIndices.size() - indicesTried);
      //cout << "rand ind:" << j << " " << indicesTried+j << " " << donorIndices.size() << endl;
      swap(donorIndices[indicesTried], donorIndices[indicesTried + j]);
      donorIndex = donorIndices[indicesTried];
      indicesTried++;

      if (donorIndex == offspringIndex)
        continue;
      
      bool dependentAreEqual = true;
      for (int k = 0; k < dependent.size(); ++k)
      {
        if (offspringPopulation[offspringIndex]->genotype[dependent[k]] != population[donorIndex]->genotype[dependent[k]])
        {
          dependentAreEqual = false;
          break;
        }
      }
      
      if (dependentAreEqual)
      {
        if (!config->donorSearch) // if not exhaustive search, stop
        {
          donorFound = true;
          break;
        }
        //if exhaustive donor search, check if donor is equal to offspring
        bool donorEqualToOffspring = true;
        for(size_t k = 0; k < FOSInstance->FOSElementSize(ind); k++)
        {
          int variableFromFOS = FOSInstance->FOSStructure[ind][k];      

          if (offspringPopulation[offspringIndex]->genotype[variableFromFOS] != population[donorIndex]->genotype[variableFromFOS])
          {
            donorEqualToOffspring = false;    
            break;
          }
        }
        if (!donorEqualToOffspring)
        {
          donorFound = true;
          break;
        }
      }

    }
    
    if (!donorFound)
      continue;

    bool donorEqualToOffspring = true;
    vector<int> touchedGenes;
    for(size_t j = 0; j < FOSInstance->FOSElementSize(ind); j++)
    {
      int variableFromFOS = FOSInstance->FOSStructure[ind][j];      
      offspringPopulation[offspringIndex]->genotype[variableFromFOS] = population[donorIndex]->genotype[variableFromFOS];

      touchedGenes.push_back(variableFromFOS);
      if (backup->genotype[variableFromFOS] != offspringPopulation[offspringIndex]->genotype[variableFromFOS])
        donorEqualToOffspring = false;    
    }

    if (!donorEqualToOffspring)
    {
      evaluateSolution(offspringPopulation[offspringIndex], backup, touchedGenes, backup->fitness);

      // accept the change if this solution is not the elitist and the fitness is at least equally good (allows random walk in neutral fitness landscape)
      // however, if this is the elitist solution, only accept strict improvements, to avoid convergence problems
      if ((!thisIsTheElitistSolution && (offspringPopulation[offspringIndex]->fitness >= backup->fitness)) || 
          (thisIsTheElitistSolution && (offspringPopulation[offspringIndex]->fitness > backup->fitness)))   
      {       
        *backup = *offspringPopulation[offspringIndex];
        
        solutionHasChanged = true;

        FOSInstance->improvementCounters[ind]++;
      }
      else
      {
        *offspringPopulation[offspringIndex] = *backup;
      }

      FOSInstance->usageCounters[ind]++;
    }

  }

  return solutionHasChanged;
}


