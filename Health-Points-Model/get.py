# Standard library imports
from collections import Counter
import importlib
from math import ceil, exp, log
from os.path import isfile
from pathlib import Path
from random import choice, choices, random
from statistics import mean

# Third-party imports
import numpy as np
from scipy.stats import genlogistic

import display

### FUNCTIONS ###
def parse_allele(allele_string, allele_name):
    ''' Parses the following number value out of file path string '''
    
    string_after_allele_name = allele_string[allele_string.index(allele_name)+len(allele_name):]
    index_after_allele_value = string_after_allele_name.find(next(filter(str.isalpha, string_after_allele_name)))
    if index_after_allele_value == -1:
        return string_after_allele_name
    else:
        return string_after_allele_name[:index_after_allele_value]


def create_config_file(file_path):
    '''Creates a file path with input info'''
    HP = "900" if not 'HP' in file_path else parse_allele(file_path,'HP')
    hp = "900" if not 'hp' in file_path else parse_allele(file_path,'hp')
    R = "10" if not 'R' in file_path else parse_allele(file_path,'R')
    r = "10" if not 'r' in file_path else parse_allele(file_path,'r')
    pop_cap_exponent = "1.2" if not 'pcE' in file_path else parse_allele(file_path,'pcE')
    pop_cap_mult = "0.1" if not 'pcM' in file_path else parse_allele(file_path,'pcM')
    TicksTillFertile = "1" if not 'TTF' in file_path else parse_allele(file_path,'TTF')
    InsultDistribution = "T7" if not 'IdT' in file_path else 'T'+parse_allele(file_path,'IdT')
    BirthProbability = "0.4" if not 'BP' in file_path else parse_allele(file_path, 'BP')
    LitterSize = "1" if not 'LS' in file_path else parse_allele(file_path,'LS')
    insult_recovery_ticks = "0" if not 'IRT' in file_path else parse_allele(file_path,'IRT')

    file_path = file_path.replace('.','pt')[:-4] + '.py'
    components = file_path.split('/')
        
    for i in range(1,len(components)):
        Path('./'+'/'.join(components[:i])).mkdir(parents=True, exist_ok=True)

    try:
        with open(file_path, 'a') as w:
            header = "import gen as c\n\n### CONFIG VARIABLES ###\n\nCONFIG = dict(c.CONFIG)\n"
            w.write(header)
            
            a = "\nCONFIG['HP'] = " + HP + "						  # the number of heath points that a member of the longer-lived species has at birth\n\
CONFIG['hp'] = " + hp + " 						  # the number of heath points that a member of the longer-lived species has at birth\n\
CONFIG['R'] = " + R + " 							  # The number of ticks between each age based decrease in HP for the longer-lived species\n\
CONFIG['r'] = " + r + " 							  # The number of ticks between each age based decrease in HP for the shorter-lived species\n\
CONFIG['LitterSize'] = " + LitterSize + " 					  # The number of NMRs born at a time\n\
CONFIG['UseLinearGrowth'] = " + str('Lin' in file_path) + "	 		  # True = linear population growth; False = exponential population growth\n\
CONFIG['InsultDistribution'] = '" + InsultDistribution + "' 		  # Choices are 'random', 'logistic', 'abs_logistic'\n\
CONFIG['BirthProbability'] = " + BirthProbability + " 		         # the probability of giving birth each time tick\n\
CONFIG['TicksTillFertile'] = " + TicksTillFertile + "		 \n\
CONFIG['insult_recovery_ticks'] = " + insult_recovery_ticks + " 		  # The number of ticks until your insult is healed\n\
CONFIG['pop_cap_exponent'] = " + pop_cap_exponent + "\n\
CONFIG['pop_cap_mult'] = " + pop_cap_mult + "\n\
CONFIG['LogorithmicPopulationCap'] = True 		# Decreases lifespan as population increases, should be used for exponential growth\n\n\
######## Generated Responses ########\n"
            w.write(a)
    except (IOError, OSError) as e:
        print(f'Error: Failed to create config file {file_path}: {e}')
        raise
    except Exception as e:
        print(f'Unexpected error creating config file {file_path}: {e}')
        raise
    
def configuration(run_args):
    ''' Import Configuration '''

    # total arguments
    num_args = len(run_args)
    # TODO: Standardize file path and name
    file_name = 'IQ-' + run_args[1] if 'UseImmortalQueenMode' in run_args and run_args.index('UseImmortalQueenMode') + 1 < len(run_args) and 'True' == run_args[run_args.index('UseImmortalQueenMode')+1] else run_args[1]
    if num_args > 1:
        if not isfile('config_files/' + file_name.replace('.', 'pt') + '.py'):
            create_config_file('config_files/' + file_name + '.py')
            
        file_path = 'config_files.' + file_name.replace('.', 'pt').replace('/', '.') 
        print(file_path)
        try:
            c = importlib.import_module(file_path)
        except (ImportError, ModuleNotFoundError) as e:
            print(f'Error: Failed to import config module {file_path}: {e}')
            raise
        except Exception as e:
            print(f'Unexpected error importing config module {file_path}: {e}')
            raise
    # else:
    #     file_path = 'config_files/wildLin.py' 
    #     import config_files.wildLin as c

    try:
        CONFIG = c.CONFIG
    except AttributeError as e:
        print(f'Error: Config module {file_path} does not have CONFIG attribute: {e}')
        raise
    CONFIG['file_path'] = 'config_files/' + file_name.replace('.', 'pt') + '.py'
    CONFIG['file_name'] =  file_name


# ######## Generated Responses ########

    CONFIG['ShowIndividualRunGraphs'] = False
    CONFIG['UseOverallRunStats'] = False		  # For testing, data analyisis and graphing
    CONFIG['GetIndividualRunData'] = True	  # For testing, data analyisis and graphing
    CONFIG['ShowIndividualRunStats'] = False	  # For testing, data analyisis and graphing
    CONFIG['StaticPopulationCap'] = False	  # NMRs only reproduce if population size is below the initial population
    CONFIG['NumRuns'] = 50
    CONFIG['StoreNMRData'] = False
    CONFIG['LongRun'] = False
    CONFIG['MaxPopulation'] = 10000			  # the maximum allowable population
    CONFIG['MaxTicks'] = 10000                # the maximum number of time ticks the simulation runs for
    CONFIG['GetAgeDistribution'] = False		  # Warning: This will take awhile to run if you don't lower max ticks to at least 2000
    CONFIG['GetLifespanDistribution'] = False  # For testing, data analyisis and graphing
    CONFIG['GetMortality'] = False             # Requires CONFIG['GetIndividualRunData'] = True to work
    CONFIG['GetEquilibriumPopulation'] = False             # Requires CONFIG['GetIndividualRunData'] = True to work
    CONFIG['UseEquilibriumForInitPop'] = True
            

    print("Total arguments passed:", num_args)
    for i in range(3, num_args, 2):
        print(run_args[i-1],'=',run_args[i], end = " ")

        # Makes sure things are properly typed
        if run_args[i] == 'True':
            CONFIG[run_args[i-1]] = True
        elif run_args[i] == 'False':
            CONFIG[run_args[i-1]] = False
        elif run_args[i] == 'None':
            CONFIG[run_args[i-1]] = None
        elif run_args[i].isnumeric():
            if '.' in run_args[i]:
                CONFIG[run_args[i-1]] = float(run_args[i])
            else:
                CONFIG[run_args[i-1]] = int(run_args[i])
        else:
            CONFIG[run_args[i-1]] = str(run_args[i])
        

    if CONFIG['LongRun']:
        CONFIG['GetAgeDistribution'] = True		  # Warning: This will take awhile to run if you don't lower max ticks to at least 2000
        CONFIG['GetLifespanDistribution'] = True  # For testing, data analyisis and graphing
        CONFIG['GetMortality'] = True             # Requires CONFIG['GetIndividualRunData'] = True to work
        CONFIG['GetEquilibriumPopulation'] = True             # Requires CONFIG['GetIndividualRunData'] = True to work
        CONFIG['UseOverallRunStats'] = True		  # For testing, data analyisis and graphing
        CONFIG['GetIndividualRunData'] = True
    elif CONFIG['GetAgeDistribution']:
        CONFIG['GetLifespanDistribution'] = False  # For testing, data analyisis and graphing
        CONFIG['GetIndividualRunData'] = False
        CONFIG['UseOverallRunStats'] = False		  # For testing, data analyisis and graphing

    ### Initialize Global Variables ###

    if CONFIG['GetEquilibriumPopulation']:         # Requires CONFIG['GetIndividualRunData'] = True to work
        CONFIG['GetIndividualRunData'] = True
    
    # Initialize colony reproduction mode defaults if not set
    if 'UseColonyReproductionMode' not in CONFIG:
        CONFIG['UseColonyReproductionMode'] = False
    if 'N0Colonies' not in CONFIG:
        CONFIG['N0Colonies'] = 5
    if 'MatingFlightInterval' not in CONFIG:
        CONFIG['MatingFlightInterval'] = 10
    if 'ThresholdColonySize' not in CONFIG:
        CONFIG['ThresholdColonySize'] = 30
    if 'MatingFlightSelection' not in CONFIG:
        CONFIG['MatingFlightSelection'] = 'random'
    if 'QueenHPMultiplier' not in CONFIG:
        CONFIG['QueenHPMultiplier'] = 2.0
    if 'BaseColonyFormationProbability' not in CONFIG:
        CONFIG['BaseColonyFormationProbability'] = 0.05
    if 'ColonyLogisticMult' not in CONFIG:
        CONFIG['ColonyLogisticMult'] = 0.1
    if 'ColonyLogisticExp' not in CONFIG:
        CONFIG['ColonyLogisticExp'] = 1.0
        
    # Sets display_factors 
    CONFIG['display_factors'] = ['Overall', 'w', 'W'] if CONFIG['UseOverallRunStats'] else ['w', 'W']

    print('\nCONFIG["display_factors"] = ', CONFIG['display_factors'])

    # Checks that we aren't trying to generate a new Age Distribution and have access to pre-generated distributions for all posible occurrances
    for pop in CONFIG['display_factors']:
        if not CONFIG['GetAgeDistribution'] and not pop in CONFIG['age_distribution']['CDF']:
            if 'DoNotUseAgeDistEver' not in CONFIG or not CONFIG['DoNotUseAgeDistEver']:
                print('Error: No Age Distribution Exists in Configuration, Generating Age Distribution (Warning: Will increase runtime)')
            if not CONFIG['LongRun']:
                CONFIG['NotUsingAgeDistOrLongRun'] = True


        
    if CONFIG['GetAgeDistribution']:
        print('Warning: Results may be skewed because age distributions are not in use,')
        print('         to fix the problem, re-run your code with the following flags')
        print("                 NumRuns 3 GetAgeDistribution False MaxTicks 3000")
        print("         Then copy the last few lines that are output, they should be in")
        print("         the form CONFIG['age_distribution']['CDF'][<pop_name>] = ___")
        print("         otherwise there could have been an error.")

    
    if CONFIG['UseEquilibriumForInitPop'] :
        CONFIG['InitialPopulation'] = '-' if not 'Overall' in CONFIG['equilibrium_population_size'] else CONFIG['equilibrium_population_size']['Overall']
        
    if not 'InitialPopulation' in CONFIG or CONFIG['InitialPopulation'] == '-':
        print('Warning: Results may be skewed because the population size is not set,')
        print('         to fix this, re-run your code with the following additional flags')
        print("                 GetEquilibriumPopulation True")
        print("         and then re-run the code.")

        CONFIG['InitialPopulation'] = 10  # the number of NMRs alive at start

    if CONFIG['GetAgeDistribution'] or ('NotUsingAgeDistOrLongRun' in CONFIG and CONFIG['NotUsingAgeDistOrLongRun'] == True) or CONFIG['DoNotUseAgeDistEver']:
        CONFIG['MaxAge'] =  ceil(CONFIG['HP']) if CONFIG['HP'] > CONFIG['hp'] else ceil(CONFIG['hp'])
    else:
        CONFIG['MaxAge'] = 2000
        for pop in CONFIG['display_factors']:
            CONFIG['MaxAge'] = len(CONFIG['age_distribution']['CDF'][pop]) + 100  if CONFIG['MaxAge'] < len(CONFIG['age_distribution']['CDF'][pop]) else ceil(CONFIG['MaxAge'])
    
    # print("CONFIG['MaxAge']",CONFIG['MaxAge'])
    # Adds comment explaining log cap equation
    if CONFIG['LogorithmicPopulationCap']:
        CONFIG['pop_cap_exponent'] = 1 if not 'pop_cap_exponent' in CONFIG else CONFIG['pop_cap_exponent']
        CONFIG['pop_cap_mult'] = .1 if not 'pop_cap_mult' in CONFIG else CONFIG['pop_cap_mult']

        if CONFIG['GetAgeDistribution']:
            data_setup = "\n\n# insult_severity = insult_severity * np.exp("+str(CONFIG['pop_cap_mult'])+" * population_size**"+ str(CONFIG['pop_cap_exponent'])+'\n'\
                + "# CONFIG['TicksTillFertile'] = " + str(CONFIG['TicksTillFertile'])

            display.write_to_file(msg = data_setup,
                        file_path = CONFIG['file_path'],
                        msg_type = 'logarithmic cap function',
                        header = '')

    print("CONFIG['InitialPopulation']",CONFIG['InitialPopulation'])
            

    return CONFIG


def random_parent(nmr_pool, tick, CONFIG, parent_sex):
  ''' Chooses a random fertile NMR from input nmr_pool '''
  
  if len(nmr_pool) > 0:
    default_parent = None
    possible_parents = choices(nmr_pool, k=len(nmr_pool))
    
    # Iterates through the randomized list of living NMRs to find the first eledgible parent
    for candidate in possible_parents:
        # checks if candidate is old enough and correct gender
        if candidate.sex == parent_sex:
            if CONFIG['TicksTillFertile'] <= tick - candidate.birth:
                return candidate
            elif default_parent is None or default_parent.birth < candidate.birth:
                # If no eligible parents are found, return the oldest candidate
                default_parent = candidate

    else:
      # This should only be used in rare instances
    #   print('Warning: No current eledgible parents found')
      return default_parent
  else:
    # print('Warning: No possible parents found')
    return 


def random_parent_alleles(nmr_pool, tick, CONFIG, parent_sex):
  ''' Chooses a random fertile NMR from input nmr_pool '''
  
  if len(nmr_pool) > 0:
    possible_parents = choices(nmr_pool, k=len(nmr_pool))
    
    # Iterates through the randomized list of living NMRs to find the first eledgible parent
    for candidate in possible_parents:
        # checks if candidate is old enough and correct gender
        if CONFIG['TicksTillFertile'] <= tick - candidate.birth and \
            candidate.sex == parent_sex:
            return candidate.w_pops['alleles']

    else:
    #   print('Warning: No eledgible parents of age found')
      return None
  else:
    # print('Warning: No possible parents found')
    return None


# Initializes Overarching Analytic Variables

def pop_death(run, tick, CONFIG, colonies=False):
    ''' Checks if a population has died out and updates if it has '''
    if CONFIG['UseImmortalQueenMode']:
      return False

                
    allele_count =  Counter(elem.w_pops['allele_name'] for elem in CONFIG['alive_NMRS']['female']) + \
                    Counter(elem.w_pops['allele_name'] for elem in CONFIG['alive_NMRS']['male']) 


    num_alive_NMRs = len(CONFIG['alive_NMRS']['female']) + len(CONFIG['alive_NMRS']['male'])
    # print("allele_count['ww']", allele_count['ww'], "allele_count['WW']", allele_count['WW'])
    if allele_count['ww'] == num_alive_NMRs or allele_count['WW'] == num_alive_NMRs:
        if CONFIG['overall_stats']['W']['pop_death_ticks'][run] == CONFIG['overall_stats']['w']['pop_death_ticks'][run]:
            if allele_count['WW'] == num_alive_NMRs:
                CONFIG['overall_stats']['w']['pop_death_count'] += 1 
                CONFIG['overall_stats']['w']['pop_death_ticks'][run] = tick
            else:
                CONFIG['overall_stats']['W']['pop_death_count'] += 1 
                CONFIG['overall_stats']['W']['pop_death_ticks'][run] = tick

            CONFIG['overall_stats']['pop_death_ticks'][run] = tick
            return True

    return False

def distribution_data( data_name, title, CONFIG):
  ''' Processes overall run data for use in any of the following distributions: Absolute, change_per_x_unit, hazard, PDF, CDF '''
  
  distributions = dict(
    Absolute = dict(),
    change_per_x_unit = dict(),
    hazard = dict(),
    PDF = dict(),
    CDF = dict()
  )

  # Scale Data Set
  for pop_name in CONFIG['display_factors']:
    CONFIG['overall_stats'][pop_name][data_name] = np.trim_zeros(CONFIG['overall_stats'][pop_name][data_name], trim='b')
    
  for pop_name in CONFIG['display_factors']:

    # Calculate PDF Data
    total_nmrs = sum(CONFIG['overall_stats'][pop_name][data_name])
    if total_nmrs > 0:

        # Save Absolute Data
        distributions['Absolute'][pop_name] = CONFIG['overall_stats'][pop_name][data_name]

        if data_name == 'equilibrium_lifespan_distribution' or data_name == 'lifespan_distribution':
            # Save change_per_x_unit Data
            try:
                distributions['change_per_x_unit'][pop_name] = [CONFIG['overall_stats'][pop_name][data_name][i] / sum(CONFIG['overall_stats'][pop_name][data_name][i:]) if sum(CONFIG['overall_stats'][pop_name][data_name][i:]) > 0 else 0 for i in range(len(CONFIG['overall_stats'][pop_name][data_name]))]
            except (ZeroDivisionError, IndexError) as e:
                print(f'Warning: Error calculating change_per_x_unit for {pop_name}: {e}')
                distributions['change_per_x_unit'][pop_name] = []

            # Save Hazard Function # To calculate it, we need to take -d log(S(t)) /dt, where S(t) is the number of individuals surviving at age t
            try:
                distributions['hazard'][pop_name] = [(log(sum(CONFIG['overall_stats'][pop_name][data_name][i:])) - log(sum(CONFIG['overall_stats'][pop_name][data_name][i+1:]))) / 1 if sum(CONFIG['overall_stats'][pop_name][data_name][i:]) > 0 and sum(CONFIG['overall_stats'][pop_name][data_name][i+1:]) > 0 else 0 for i in range(len(CONFIG['overall_stats'][pop_name][data_name]) - 1)]
            except (ValueError, IndexError) as e:
                print(f'Warning: Error calculating hazard function for {pop_name}: {e}')
                distributions['hazard'][pop_name] = []

        # Save PDF Data
        try:
            distributions['PDF'][pop_name] = [i / total_nmrs for i in CONFIG['overall_stats'][pop_name][data_name]]
        except ZeroDivisionError:
            print(f'Warning: total_nmrs is zero for {pop_name}, cannot calculate PDF')
            distributions['PDF'][pop_name] = []

        # Save CDF Data
        try:
            distributions['CDF'][pop_name] =  list(np.cumsum(distributions['PDF'][pop_name], dtype=float))
            distributions['CDF'][pop_name] = [float(i) for i in distributions['CDF'][pop_name]] # Stupid bit of code to fix numpy update bug
        except (ValueError, IndexError) as e:
            print(f'Warning: Error calculating CDF for {pop_name}: {e}')
            distributions['CDF'][pop_name] = []

    else:
        print('Error: there were no recorded instances of', title, 'data for the population', pop_name)

  return distributions

def insult_severity(CONFIG, colony_population_size, nmr=False, tick=False, after_pop_death = False, equilibrium_population_sizes = False):
    ''' Assigns Severity of insults using insult distribution set by config file'''

    if CONFIG['InsultDistribution'] == 'random':
        insult_severity = choice(range(CONFIG['HP'])) # Chooses a random number between 1 and the initial HP of the dominant allele
    # print(insult_severity)
    elif CONFIG['InsultDistribution'] == 'abs_logistic':
        insult_severity = abs(genlogistic.rvs(0.5799))
    elif CONFIG['InsultDistribution'] == 'logistic':
        insult_severity = genlogistic.rvs(.2, loc = 10)
        insult_severity = 0 if insult_severity < 1 else insult_severity
    if CONFIG['InsultDistribution'] == 'A1': # A1. Insult distribution 1:  y = x^(0.5)
        # Simulation 1:  Insult distribution 1:  y = x^(0.5)	(i.e., square root x, concave, decreasing slope).  Sample x from [1:1000000] so that max y = 1000.
        insult_severity = choice(range(1,1000000))**(0.5)
    elif CONFIG['InsultDistribution'] == 'A2': # A2. Insult distribution 2:  y = x 
        # Simulation 2:  Insult distribution 2:  y = x   (i.e., straight line).  Sample x from [1:1000] so that max y = 1000.
        insult_severity = choice(range(1,1000))
    elif CONFIG['InsultDistribution'] == 'A3': # A3. Insult distribution 3:  y = x^(1.5) 
        # Simulation 3:  Insult distribution 3:  y = x^(1.5)  (i.e., exponential, convex, increasing slope.  Sample x from [1: 100] so that max y = 1000.
        insult_severity = choice(range(1,100))**(1.5)
    elif CONFIG['InsultDistribution'] == 'B1': # B1. y = 1000 (1/x)	Sample from [1:100] 
        insult_severity = 1000 / choice(range(1,100))
    elif CONFIG['InsultDistribution'] == 'B2': # B2. y = 1000 (1/x)   Sample from [1:1000]
        insult_severity = 1000 / (random()*1000)#choice(range(1,1000))
    elif CONFIG['InsultDistribution'] == 'B3': # B3. y = 1000 (1/x)  Sample from {1:10000]
        insult_severity = 1000 / choice(range(1,10000))
    elif CONFIG['InsultDistribution'] == 'B4': # B4. y = 1000 (1/x)  Sample from [1:100000]
        insult_severity = 1000 / choice(range(1,100000))
    elif CONFIG['InsultDistribution'] == 'C1': #  C1. y = 1000 – e^(0.1 x)	Sample from [0:7], but sample every 0.01, so that there are 700 possible numbers, max y  = 999, min y is about 0.  This is concave, negative slope.  There might be a problem around x = 7 because the expression might go negative, I didn’t exactly multiply it out.
        insult_severity = 1000 - exp((choice(range(0,700)))/100)
    elif CONFIG['InsultDistribution'] == 'C2': #  C2. y = 1000 – x     Sample from [0:1000].  Linear and negative.
        insult_severity = 1000 - choice(range(1,1000))
    elif CONFIG['InsultDistribution'] == 'C3': #  C3. y = 1000 (1/x) 	Sample from [1:1000].  Convex and negative.
        insult_severity = 1000 / choice(range(1,1000))
    elif CONFIG['InsultDistribution'] == 'T1': # 2000 ( (1/x^2) (e^(-1/x)) ), sampling from [0.1 : 100] by units of 0.1
        rand = choice(range(1,1000)) / 10
        insult_severity = (2000 * exp(-1/rand)) / (rand**2)
    elif CONFIG['InsultDistribution'] == 'T2': # 2000 ( (1/x^2) (e^(-1/x)) ), sampling from [0.1 : 200] by units of 0.1
        rand = choice(range(1,2000)) / 10
        insult_severity = (2000 * exp(-1/rand)) / (rand**2)
    elif CONFIG['InsultDistribution'] == 'T3': 
        rand = random()*100 #choice(range(1,10000000)) / 100000
        insult_severity = (4*1000 * exp(-1/rand)) / (rand**2)
    elif CONFIG['InsultDistribution'] == 'T4': 
        rand = random()*1000 #1000 (1/x^2)  , sampling from [1:1000]
        insult_severity = 1000 / (rand**2)
    elif CONFIG['InsultDistribution'] == 'T5': 
        rand = random()*100 #1000 (1/x^2)  , sampling from [1:100]
        insult_severity = 1000 / (rand**2)
    elif CONFIG['InsultDistribution'] == 'T6': 
        rand = random()*50 #1000 (1/x^2)  , sampling from [1:50]
        insult_severity = 1000 / (rand**2)
    elif CONFIG['InsultDistribution'] == 'T7': 
        # print('T7')
        rand = random()*50000 
        insult_severity = 1000 / (rand**2)
    elif CONFIG['InsultDistribution'] == 'T8':
        init_hp = CONFIG[nmr.hp_type_name] + .1 # This is the biggest liberty I took, obviously dividing by log 1 is an issue but init_hp = either 1000 or 1200, I figured that .1 didn't make much difference there
        rand = random() # just a random number from 0 to 1
        alpha = .01 #The population size was too low with .1
        age = abs(nmr.birth - tick) + 1
        beta = -(log(nmr.current_health_points/init_hp)) / age # where I have assumed that health_points(t) = init_hp * exp(-beta t), which defines the parameter beta;
        dt = 1 # you do an insult every dt in time, defining the parameter dt: An insult is applied every tick
        insult_severity = (init_hp * alpha * (exp(beta * dt) - 1)) / (beta * log(1/rand))  
        # So, alpha is arbitrary, but beta and dt are determined by the other parameters of your simulation, and the correct values need to be used otherwise it won't work.
    elif CONFIG['InsultDistribution'] == 'T8B':
        init_hp = CONFIG[nmr.hp_type_name] + .1 # This is the biggest liberty I took, obviously dividing by log 1 is an issue but init_hp = either 1000 or 1200, I figured that .1 didn't make much difference there
        rand = random() # just a random number from 0 to 1
        alpha = .01 #The population size was too low with .1
        age = abs(nmr.birth - tick) + 1
        beta =  -log(2) / nmr.r_value_current # where I have assumed that health_points(t) = init_hp * exp(-beta t), which defines the parameter beta;
        insult_severity = (init_hp * alpha * (exp(beta ) - 1)) / (beta * log(1/rand))  
    elif CONFIG['InsultDistribution'] == 'T9': 
        rand = random()*50 
        insult_severity = 1000 / (rand**3)
    elif CONFIG['InsultDistribution'] == 'T99': 
        insult_severity = 5

    if CONFIG['LogorithmicPopulationCap']:
        if CONFIG['UseAvgForLogCap']:
          pop_size = mean(CONFIG['last_3_pop_sizes'])
        else:
          pop_size = colony_population_size if not after_pop_death and not colony_population_size == None else CONFIG['InitialPopulation']
          insult_severity = insult_severity * np.exp(CONFIG['pop_cap_mult'] * (pop_size**CONFIG['pop_cap_exponent']))

    if not 'insults' in CONFIG:
        CONFIG['insults'] = []
    CONFIG['insults'].append(insult_severity)  
      
    return insult_severity

def get_run_mode(CONFIG):
  ''' Returns the mode of the run as a string '''
  if not CONFIG['UseColonyReproductionMode']:
    mode = 'Individual'
  elif CONFIG['MatesWithinColony']:
    mode = 'Mole Rat'
  else:
    mode = 'Ant'
  return mode

def get_longevity_allele(CONFIG):
    ''' Returns the longevity allele as a string '''
    longevity = 'W' if CONFIG['R'] > CONFIG['r'] else 'w'
    return longevity

def get_longevity_or_vitality_allele(CONFIG, allele_name):
    ''' Returns the longevity or vitality allele as a string '''
    longevity = 'W' if CONFIG['R'] > CONFIG['r'] else 'w'
    if allele_name == longevity:
        return 'Longevity'
    else:
        return 'Vitality'