#include <iostream>
#include <fstream>
using namespace std;

#include "gomeaIMS.hpp"
#include "utils.hpp"

gomeaIMS::gomeaIMS(Config *config_): GOMEA(config_)
{
  prepareFolder(config->folder);
  initElitistFile(config->folder, config->populationScheme, config->populationSize);

  maximumNumberOfGOMEAs         = config->maximumNumberOfGOMEAs;
  IMSsubgenerationFactor        = config->IMSsubgenerationFactor;
  basePopulationSize            = config->basePopulationSize;
  numberOfGOMEAs                = 0;
  numberOfGenerationsIMS        = 0;
  minimumGOMEAIndex             = 0;
  numberOfGenerationsIMS        = 0;
  generationsWithoutImprovement = 0;

  createProblemInstance(config->problemIndex, config->numberOfVariables, config, &problemInstance, config->problemInstancePath);
  #ifdef DEBUG
    cout << "Problem Instance created! Problem number is " << config->problemIndex << endl;
  #endif

  sharedInformationInstance = new sharedInformation(config->maxArchiveSize);
  #ifdef DEBUG
    cout << "Shared Information instance created!\n";
  #endif
}

gomeaIMS::~gomeaIMS()
{
  for (int i = 0; i < numberOfGOMEAs; ++i)
    delete GOMEAs[i];

  delete problemInstance;
  delete sharedInformationInstance;
}


void gomeaIMS::run()
{
  while(!checkTermination())
  {
    if (numberOfGOMEAs < maximumNumberOfGOMEAs)
      initializeNewGOMEA();

    generationalStepAllGOMEAs();

    numberOfGenerationsIMS++;
  }
}

bool gomeaIMS::checkTermination()
{
  int i;

  if (numberOfGOMEAs == maximumNumberOfGOMEAs)
  {
    for (i = 0; i < maximumNumberOfGOMEAs; i++)
    {
      if (!GOMEAs[i]->terminated)
        return false;
    }

    return true;
  }
  
  return false;
}

void gomeaIMS::initializeNewGOMEA()
{
  #ifdef DEBUG
    cout << "Current number Of GOMEAs is " << numberOfGOMEAs << " | Creating New GOMEA!\n";
  #endif


  Population *newPopulation = NULL;

  if (numberOfGOMEAs == 0)
    newPopulation = new Population(config, problemInstance, sharedInformationInstance, numberOfGOMEAs, basePopulationSize);
  else
    newPopulation = new Population(config, problemInstance, sharedInformationInstance, numberOfGOMEAs, 2 * GOMEAs[numberOfGOMEAs-1]->populationSize);
  
  GOMEAs.push_back(newPopulation);
  numberOfGOMEAs++;
}

void gomeaIMS::generationalStepAllGOMEAs()
{
  int GOMEAIndexSmallest, GOMEAIndexBiggest;

  GOMEAIndexBiggest  = numberOfGOMEAs - 1;
  GOMEAIndexSmallest = 0;
  while(GOMEAIndexSmallest <= GOMEAIndexBiggest)
  {
    if (!GOMEAs[GOMEAIndexSmallest]->terminated)
      break;

    GOMEAIndexSmallest++;
  }

  GOMEAGenerationalStepAllGOMEAsRecursiveFold(GOMEAIndexSmallest, GOMEAIndexBiggest);
}

void gomeaIMS::GOMEAGenerationalStepAllGOMEAsRecursiveFold(int GOMEAIndexSmallest, int GOMEAIndexBiggest)
{
  int i, GOMEAIndex;

  for(i = 0; i < IMSsubgenerationFactor-1; i++)
  {
    for(GOMEAIndex = GOMEAIndexSmallest; GOMEAIndex <= GOMEAIndexBiggest; GOMEAIndex++)
    {
      if(!GOMEAs[GOMEAIndex]->terminated)
        GOMEAs[GOMEAIndex]->terminated = checkTerminationGOMEA(GOMEAIndex);

      if((!GOMEAs[GOMEAIndex]->terminated) && (GOMEAIndex >= minimumGOMEAIndex))
      {
        GOMEAs[GOMEAIndex]->calculateAverageFitness();
        double fitness_before = sharedInformationInstance->elitist.fitness;
        //cout << GOMEAIndex << " | avgFitness " << GOMEAs[GOMEAIndex]->averageFitness << endl;

        GOMEAs[GOMEAIndex]->makeOffspring();

        GOMEAs[GOMEAIndex]->copyOffspringToPopulation();

        GOMEAs[GOMEAIndex]->calculateAverageFitness();
        double fitness_after = sharedInformationInstance->elitist.fitness;
        //cout <<  GOMEAIndex << " " << GOMEAs[GOMEAIndex]->numberOfGenerations << " " << sharedInformationInstance->numberOfEvaluations << " " << GOMEAs[GOMEAIndex]->population.size() << " " << GOMEAs[GOMEAIndex]->averageFitness << endl;
        if (config->populationScheme==0 && sharedInformationInstance->numberOfEvaluations >= 100000000)
        {
         cout << "10^8 evals reached! Terminating...\n";
         GOMEAs[GOMEAIndex]->terminated = true;
        }
        GOMEAs[GOMEAIndex]->numberOfGenerations++;

        if (GOMEAs[GOMEAIndex]->numberOfGenerations >= config->maxGenerations)
          GOMEAs[GOMEAIndex]->terminated = true;
      }
    }

    for(GOMEAIndex = GOMEAIndexSmallest; GOMEAIndex < GOMEAIndexBiggest; GOMEAIndex++)
      GOMEAGenerationalStepAllGOMEAsRecursiveFold(GOMEAIndexSmallest, GOMEAIndex);
  }
}

bool gomeaIMS::checkTerminationGOMEA(int GOMEAIndex)
{
  for (int i = GOMEAIndex+1; i < numberOfGOMEAs; i++)
  {    
    if (GOMEAs[i]->averageFitness > GOMEAs[GOMEAIndex]->averageFitness)
    {
      minimumGOMEAIndex = GOMEAIndex+1;
      //cout << "" << i << ":" << GOMEAs[i]->averageFitness << " | " << GOMEAIndex << ":" << GOMEAs[GOMEAIndex]->averageFitness << endl;
      return true;
    }
  }

  for (size_t i = 1; i < GOMEAs[GOMEAIndex]->populationSize; i++)
  {
    for (size_t j = 0; j < config->numberOfVariables; j++)
    {
      if (GOMEAs[GOMEAIndex]->population[i]->genotype[j] != GOMEAs[GOMEAIndex]->population[0]->genotype[j])
      {
        return false;
      }
    }
  }

//cout << "CONVERGED!\n";
          

  return true;
}
