#include "Config.hpp"

/**
 * Parses the command line.
 * For options, see printUsage.
 */
void Config::splitString(const string &str, vector<string> &splitted, char delim)
{
    size_t current, previous = 0;
    current = str.find(delim);
    while (current != string::npos)
    {
        splitted.push_back(str.substr(previous, current - previous));
        previous = current + 1;
        current = str.find(delim, previous);
    }
    splitted.push_back(str.substr(previous, current - previous));
}

bool Config::isNumber(const string &str)
{
	return !str.empty() && all_of(str.begin(), str.end(), ::isdigit);
}

bool Config::parseCommandLine(int argc, char **argv)
{
  const struct option longopts[] =
  {
	{"help",        no_argument,         0, 'h'},    
    {"FI",          no_argument,         0, 'f'},
    {"partial",     no_argument,         0, 'p'},
    {"saveEvals",   no_argument,         0, 's'},  
    {"donorSearch", no_argument,         0, 'd'},    
    {"tournamentSelection", no_argument, 0, 't'},       
    {"GOM",         required_argument,   0, 'C'},  
    {"hillClimber", required_argument,   0, 'H'},    
		{"problem",     required_argument,   0, 'P'},  
		{"L",           required_argument,   0, 'L'},  
		{"FOS",         required_argument,   0, 'F'},
		{"alphabet",    required_argument,   0, 'A'},
		{"scheme",      required_argument,   0, 'X'},
		{"instance",    required_argument,   0, 'I'},
		{"vtr",         required_argument,   0, 'V'},
		{"time",        required_argument,   0, 'T'},
		{"folder",      required_argument,   0, 'O'},
		{"threshold",   required_argument,   0, 'M'},        
    {"orderFOS",    required_argument,   0, 'B'}, 
    {"similarityMeasure", required_argument,   0, 'Z'}, 
    {"populationSize", required_argument,   0, 'Q'}, 
    {"seed",        required_argument,   0, 'S'},
		           
    {0,             0,                   0,  0 }
  };


  int c, index;
  while ((c = getopt_long(argc, argv, "h::c::f::p::s::d::t::C::P::F::L::O::T::S::V::A::M::I::X::B::Z::Q::H::", longopts, &index)) != -1)
  {
  	switch (c)
	{
		case 'h':
			printUsage();
			exit(0);
		case 'f':
			useForcedImprovements = 1;
			break;
		case 'p':
			usePartialEvaluations = 1;
			break;
		case 's':
			saveEvaluations = 1;
			break;
		case 'd':
			donorSearch = 1;
			break;
		case 't':
			tournamentSelection = 1;
			break;

		case 'C':
			conditionalGOM = atoi(optarg);
			break;
		case 'P':
			{
				const string optarg_str = string(optarg);
				if (isNumber(optarg_str)) 	
					problemIndex = atoi(optarg);
				else
				{
					vector<string> tmp;
					splitString(optarg_str, tmp, '_');
					cout << tmp[0] << " " << tmp[1] << " " << tmp[2] << endl;	
				
					problemIndex = atoi(tmp[0].c_str());
					k = atoi(tmp[1].c_str());
					s = atoi(tmp[2].c_str());
					cout << problemIndex << endl;
				}
				
			}
			break;
		case 'F':
			FOSIndex = atoi(optarg);
			break;
		case 'H':
			hillClimber = atoi(optarg);
			break;
		case 'L':
			numberOfVariables = atoi(optarg);
			break;
		case 'M':
			MI_threshold = atof(optarg);
			break;
		case 'O':
			folder= string(optarg);
			break;
		case 'T':
			timelimitMilliseconds = atoi(optarg) * 1000;
			break;
		case 'V':
			vtr = atof(optarg);
			break;
		case 'S':
			randomSeed = atoll(optarg);
			break;
		case 'A':
			alphabetSize = atoi(optarg);
			break;
		case 'I':
			problemInstancePath = string(optarg);
			break;
		case 'X':
			populationScheme = atoi(optarg);
			break;
		case 'Q':
			populationSize = atoi(optarg);
			break;
		case 'B':
			orderFOS = atoi(optarg);
			break;
		case 'Z':
			similarityMeasure = atoi(optarg);
			break;

		default:
			abort();
	}
  }

  if (populationScheme == 0) //fixed pop size
  {
  	maximumNumberOfGOMEAs = 1;
  	basePopulationSize = populationSize;
  }

  if (printHelp)
  {
  	printUsage();
  	exit(0);
  }
  return 1;
}

void Config::printUsage()
{
  cout << "Usage: GOMEA [-h] ...\n";
  cout << "   -h: Prints out this usage information.\n";

  cout << endl;
  cout << "General settings: \n";
  cout << "    --L: Number of variables. Default: 1\n";
  cout << "    --alphabet: Alphabet size. Default: 2 (binary optimization)\n";  
  cout << "    --problem: Index of optimization problem to be solved (maximization). Default: 0\n";
  cout << "    --instance: Problem instance. \n";  
  cout << "    --vtr: Value To Reach value. \n";  
  cout << "    --time: Time Limit (in seconds). \n";  
  cout << "    --folder: Folder where to save results. \n";  
  cout << "    --seed: Random seed. \n";  
  cout << "    -partial: Enables partial evaluations. Default: disabled.\n";
  cout << "    -saveEvals: Enables saving all evaluations in archive. Default: disabled\n";

  cout << endl;
  cout << "GOMEA Configuration settings: \n";  
  cout << "    --GOM: GOM type. Default: 0 (LT). 1 - conditionalGOM (CGOM) based on MI\n";
  cout << "    --threshold: Threshold value for CGOM. Default: 0.8\n";
  cout << "    --hillClimber: Hill Climber usage. Default: 0 (no HC). 1 - single HC. 2 - Exhaustive HC\n"; 
  cout << "    --FOS: FOS type. Default: 0 (LT). 1 - Filtered LT. 2 - Efficient implementation of LT for P3 (without Tournament Selection). 4 - Efficient implementation of filtered LT for P3 (without Tournament Selection).\n";
  cout << "    --orderFOS: FOS order. Default: 0 (randomly shuffled). 1 - sorted by the FOS elements size (ascending order).\n";
  cout << "    --similarityMeasure: FOS building similarity measure. Default: 0 (MI). 1 - NMI.\n";
  cout << "    --scheme: Population Management Scheme. Default: 0 (single population). 1 - IMS. 2 - P3. 3 - P3-MI (Quadratic). 4 - P3-MI (Linear).\n";
  cout << "    --populationSize: Population Size for single population run. \n";
  cout << "    -FI: Enables Forced Improvements. Default: disabled.\n";
  cout << "    -donorSearch: Enables Exhaustive Donor Search Default: disabled.\n";
  cout << "    -tournamentSelection: Enables Tournament Selection of size 2 prior to Linkage Model Learning\n";

  cout << endl;

}


void Config::printOverview()
{
  cout << "### Settings ######################################\n";
  cout << "#\n";
  cout << "# Use Forced Improvements : " << (useForcedImprovements ? "enabled" : "disabled")  << endl;
  cout << "# Use partial evaluations : " << (usePartialEvaluations ? "enabled" : "disabled")  << endl;
  cout << "# Save all evaluations : " << (saveEvaluations ? "enabled" : "disabled") << endl;
  string conditionalGOMDesc = "Unconditional";
  if (conditionalGOM == 1)
  	conditionalGOMDesc = "Conditional by MI, MI threshold="+to_string(MI_threshold);
  
  cout << "# conditionalGOM : " << (conditionalGOMDesc) << endl;
  string populationSchemeName = "Fixed Population Size";
  if (populationScheme == 1)
	populationSchemeName = "IMS";
  else if (populationScheme == 2)
  	populationSchemeName = "P3";  
  else if (populationScheme == 3)
  	populationSchemeName = "P3-Multiple Insertion Quadratic";  
  else if (populationScheme == 4)
  	populationSchemeName = "P3-Multiple Insertion Linear";

  cout << "# population scheme : " << populationSchemeName << endl;
  if (populationScheme == 0)
  	cout << "# population size : " << populationSize << endl;
  if (hillClimber == 0)
	  cout << "# use hill climber : " << "disabled" << endl;
  else if (hillClimber == 1)
	  cout << "# use hill climber : " << "single" << endl;
  else if (hillClimber == 2)
	  cout << "# use hill climber : " << "multiple" << endl;

  cout << "# use exhaustive donor search : " << (donorSearch ? "enabled" : "disabled") << endl;
  cout << "# use tournament selection : " << (tournamentSelection ? "enabled" : "disabled") << endl;
  cout << "# similarity measure : " << (similarityMeasure ? "normalized MI" : "MI") << endl;
  cout << "# FOS ordering : " << (orderFOS ? "ascending" : "random") << endl;
  
  cout << "#\n";
  cout << "###################################################\n";
  cout << "#\n";
  cout << "# Problem                      = " << problemName << endl;
  cout << "# Problem Instance Filename    = " << problemInstancePath << endl;
  cout << "# FOS                          = " << FOSName << endl;
  cout << "# Number of variables          = " << numberOfVariables << endl;
  cout << "# Alphabet size                = " << alphabetSize << endl;
  cout << "# Time Limit (seconds)         = " << timelimitMilliseconds/1000.0 << endl;
  cout << "# VTR                          = " << ((vtr < 1e+308) ? to_string(vtr) : "not set") << endl;
  cout << "# Random seed                  = " << randomSeed << endl;
  cout << "# Folder                       = " << folder << endl;
  cout << "#\n";
  cout << "### Settings ######################################\n";
}

void Config::checkOptions()
{
	if (!problemNameByIndex(this, problemName))
	{
		cerr << "No problem with index " << problemIndex << " installed!";
		exit(0);
	}

	if (!FOSNameByIndex(FOSIndex, FOSName))
	{
		cerr << "No FOS with index " << FOSIndex << " installed!";
		exit(0);
	}
}
