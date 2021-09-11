#include <iostream>
#include <fstream>
using namespace std;

#include "gomeaP3_MI.hpp"
#include "utils.hpp"

gomeaP3_MI::gomeaP3_MI(Config *config_): GOMEA(config_)
{
  prepareFolder(config->folder);
  initElitistFile(config->folder, config->populationScheme, config->populationSize);

  maximumNumberOfGOMEAs       = 1000000000;
  numberOfGOMEAs              = 0;
  
  createProblemInstance(config->problemIndex, config->numberOfVariables, config, &problemInstance, config->problemInstancePath);
  #ifdef DEBUG
    cout << "Problem Instance created! Problem number is " << config->problemIndex << endl;
  #endif

  sharedInformationInstance = new sharedInformation(config->maxArchiveSize);
  #ifdef DEBUG
    cout << "Shared Information instance created!\n";
  #endif
}

gomeaP3_MI::~gomeaP3_MI()
{
  for (int i = 0; i < numberOfGOMEAs; ++i)
    delete GOMEAs[i];

  delete problemInstance;
  delete sharedInformationInstance;
}

void gomeaP3_MI::run()
{
  while(!checkTermination())
  {
    if (numberOfGOMEAs < maximumNumberOfGOMEAs)
      initializeNewGOMEA();

    GOMEAGenerationalSteps(numberOfGOMEAs-1);
  }
}

bool gomeaP3_MI::checkTermination()
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

void gomeaP3_MI::initializeNewGOMEA()
{
  #ifdef DEBUG
    cout << "Current number Of GOMEAs is " << numberOfGOMEAs << " | Creating New GOMEA!\n";
  #endif


  Population_P3_MI *newPopulation = NULL;

  int populationSize;
  if (config->populationScheme == 2)
    populationSize = 1; //P3
  else if (config->populationScheme == 3)
    populationSize = (numberOfGOMEAs+1)*(numberOfGOMEAs+1); //P3-MI quadratic
  else if (config->populationScheme == 4)
    populationSize = numberOfGOMEAs+1; //P3-MI linear
  
  newPopulation = new Population_P3_MI(config, problemInstance, sharedInformationInstance, numberOfGOMEAs, populationSize);
  
  GOMEAs.push_back(newPopulation);
  numberOfGOMEAs++;
}

void gomeaP3_MI::GOMEAGenerationalSteps(int GOMEAIndex)
{
  while (true)
  {
    if(!GOMEAs[GOMEAIndex]->terminated)
    {
      GOMEAs[GOMEAIndex]->makeOffspring();

      GOMEAs[GOMEAIndex]->copyOffspringToPopulation();

      GOMEAs[GOMEAIndex]->currentPyramidLevel++;
    }
    else
      break;
  }
}


