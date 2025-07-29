from genericpath import isfile
from math import ceil, log, exp
from pathlib import Path
from statistics import mean
import numpy as np
from collections import Counter
import importlib

from random import choices, random, choice
from scipy.stats import genlogistic

from display import write_to_file
from re import findall

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
    UseTradeoffAllele = True if 'w' in file_path else False
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

    w = open(file_path,'a')
    header = "import gen as c\n\n### CONFIG VARIABLES ###\n\nCONFIG = dict(c.CONFIG)\n"
    w.write(header)
    # components = file_path.split('/')[1:]
    
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
CONFIG['LogorithmicPopulationCap'] = True 		# Decreases lifespan as population increases, should be used for exponential growth\n\
CONFIG['UseTradeoffAllele'] = " + str(UseTradeoffAllele) + "\n\n\
######## Generated Responses ########\n"

    w.write(a)
    w.close()
    
def configuration(run_args):
    ''' Import Configuration '''

    # total arguments
    num_args = len(run_args)
    # CONFIG = dict()
    if num_args > 1:
        if not isfile('config_files/' + run_args[1].replace('.', 'pt') + '.py'):
            create_config_file('config_files/' + run_args[1] + '.py')
            
        file_path = 'config_files.' + run_args[1].replace('.', 'pt').replace('/', '.') 
        print(file_path)
        c = importlib.import_module(file_path)
    # else:
    #     file_path = 'config_files/wildLin.py' 
    #     import config_files.wildLin as c

    CONFIG = c.CONFIG
    CONFIG['file_path'] = 'config_files/' + run_args[1].replace('.', 'pt') + '.py'
    CONFIG['file_name'] =  run_args[1] 


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
    CONFIG['GetPopAFraction'] = False
    CONFIG['GetAgeDistribution'] = False		  # Warning: This will take awhile to run if you don't lower max ticks to at least 2000
    CONFIG['GetLifespanDistribution'] = True  # For testing, data analyisis and graphing
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
        elif '.' in run_args[i]:
            CONFIG[run_args[i-1]] = float(run_args[i])
        else:
            CONFIG[run_args[i-1]] = int(run_args[i])
        

    if CONFIG['LongRun']:
        CONFIG['GetAgeDistribution'] = True		  # Warning: This will take awhile to run if you don't lower max ticks to at least 2000
        CONFIG['GetLifespanDistribution'] = True  # For testing, data analyisis and graphing
        CONFIG['GetMortality'] = True             # Requires CONFIG['GetIndividualRunData'] = True to work
        CONFIG['GetEquilibriumPopulation'] = True             # Requires CONFIG['GetIndividualRunData'] = True to work
        CONFIG['GetPopAFraction'] = False           # Would take too long to run on long run
        CONFIG['UseOverallRunStats'] = True		  # For testing, data analyisis and graphing
        CONFIG['GetIndividualRunData'] = True
    elif CONFIG['GetAgeDistribution']:
        CONFIG['GetLifespanDistribution'] = False  # For testing, data analyisis and graphing
        CONFIG['GetIndividualRunData'] = False
        CONFIG['UseOverallRunStats'] = False		  # For testing, data analyisis and graphing
    elif CONFIG['GetPopAFraction']:
        CONFIG['GetLifespanDistribution'] = False  # For testing, data analyisis and graphing
    CONFIG['GetLifespanDistribution'] = True
    ### Initialize Global Variables ###

    if CONFIG['GetEquilibriumPopulation']:         # Requires CONFIG['GetIndividualRunData'] = True to work
        CONFIG['GetIndividualRunData'] = True
        
    # Sets display_factors 
    CONFIG['display_factors'] = ['Overall'] if CONFIG['UseOverallRunStats'] else []

    if CONFIG['HP'] != CONFIG['hp'] or CONFIG['R'] == CONFIG['r'] or CONFIG['UseTradeoffAllele']:
        CONFIG['show_hp'] = True
        CONFIG['display_factors'] += ['hp', 'HP']
    else:
        CONFIG['show_hp'] = False

    if CONFIG['R'] != CONFIG['r'] or CONFIG['HP'] == CONFIG['hp'] or CONFIG['UseTradeoffAllele']:
        CONFIG['show_r'] = True
        CONFIG['display_factors'] += ['r', 'R']
    else:
        CONFIG['show_r'] = False

    print('\nCONFIG["display_factors"] = ', CONFIG['display_factors'])

    # Checks that we aren't trying to generate a new Age Distribution and have access to pre-generated distributions for all posible occurrances
    for pop in CONFIG['display_factors']:
        if not CONFIG['GetAgeDistribution'] and not pop in CONFIG['age_distribution']['CDF']:
            print('Error: No Age Distribution Exists in Configuration, Generating Age Distribution (Warning: Will increase runtime)')
            CONFIG['GetAgeDistribution'] = True

        
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

    if CONFIG['GetAgeDistribution']:
        CONFIG['MaxAge'] =  ceil(CONFIG['HP']) if CONFIG['HP'] > CONFIG['hp'] else ceil(CONFIG['hp'])
    else:
        CONFIG['MaxAge'] = 500
        for pop in CONFIG['display_factors']:
            CONFIG['MaxAge'] = len(CONFIG['age_distribution']['CDF'][pop]) + 100  if CONFIG['MaxAge'] < len(CONFIG['age_distribution']['CDF'][pop]) else ceil(CONFIG['MaxAge'])

    # Adds comment explaining log cap equation
    if CONFIG['LogorithmicPopulationCap']:
        CONFIG['pop_cap_exponent'] = 1 if not 'pop_cap_exponent' in CONFIG else CONFIG['pop_cap_exponent']
        CONFIG['pop_cap_mult'] = .1 if not 'pop_cap_mult' in CONFIG else CONFIG['pop_cap_mult']

        if CONFIG['GetAgeDistribution']:
            data_setup = "\n\n# insult_severity = insult_severity * np.exp("+str(CONFIG['pop_cap_mult'])+" * population_size**"+ str(CONFIG['pop_cap_exponent'])+'\n'\
                + "# CONFIG['TicksTillFertile'] = " + str(CONFIG['TicksTillFertile'])

            write_to_file(msg = data_setup,
                        file_path = CONFIG['file_path'],
                        msg_type = 'logarithmic cap function')

    print("CONFIG['InitialPopulation']",CONFIG['InitialPopulation'])
            

    return CONFIG


def random_parent(nmr_pool, tick, CONFIG):
  ''' Chooses a random fertile NMR from input nmr_pool '''
  
  if len(nmr_pool) > 0:
    possible_parents = choices(nmr_pool)
    default_parent = possible_parents[0]

    while len(possible_parents) > 0 and CONFIG['TicksTillFertile'] >= tick \
        - possible_parents[0].birth:

  # Iterates through the randomized list of living NMRs to find the first eledgible parent

        possible_parents.pop()

    if len(possible_parents) > 0:
        return possible_parents[0]
    else:
      # This should only be used in rare instances
      return default_parent
  else:
    return 




# Initializes Overarching Analytic Variables

def pop_death(run, tick, CONFIG):
    ''' Checks if a population has died out and updates if it has '''
                
    allele_count = Counter(elem.hp_pops['allele_name'] for elem in CONFIG['alive_NMRS']['female']) + \
                    Counter(elem.r_pops['allele_name'] for elem in CONFIG['alive_NMRS']['female']) + \
                    Counter(elem.w_pops['allele_name'] for elem in CONFIG['alive_NMRS']['female']) + \
                    Counter(elem.hp_pops['allele_name'] for elem in CONFIG['alive_NMRS']['male']) + \
                    Counter(elem.r_pops['allele_name'] for elem in CONFIG['alive_NMRS']['male']) + \
                    Counter(elem.w_pops['allele_name'] for elem in CONFIG['alive_NMRS']['male']) 

    # print('allele_count',allele_count)

    num_alive_NMRs = len(CONFIG['alive_NMRS']['female']) + len(CONFIG['alive_NMRS']['male'])

    if CONFIG['UseTradeoffAllele']:
        if allele_count['ww'] == num_alive_NMRs or allele_count['WW'] == num_alive_NMRs:
            if CONFIG['overall_stats']['HP']['pop_death_ticks'][run] == CONFIG['overall_stats']['hp']['pop_death_ticks'][run]:
                if allele_count['WW'] == num_alive_NMRs:
                    CONFIG['overall_stats']['hp']['pop_death_count'] += 1 
                    CONFIG['overall_stats']['hp']['pop_death_ticks'][run] = tick
                    CONFIG['overall_stats']['r']['pop_death_count'] += 1 
                    CONFIG['overall_stats']['r']['pop_death_ticks'][run] = tick
                else:
                    CONFIG['overall_stats']['HP']['pop_death_count'] += 1 
                    CONFIG['overall_stats']['HP']['pop_death_ticks'][run] = tick
                    CONFIG['overall_stats']['R']['pop_death_count'] += 1 
                    CONFIG['overall_stats']['R']['pop_death_ticks'][run] = tick

                CONFIG['overall_stats']['pop_death_ticks'][run] = tick
                return True
    
    else:
        if CONFIG['show_hp'] and (allele_count['hphp'] == num_alive_NMRs or allele_count['HPHP'] == num_alive_NMRs):
            # checks that both populations still exist
            if CONFIG['overall_stats']['hp']['pop_death_ticks'][run] == CONFIG['overall_stats']['HP']['pop_death_ticks'][run]:

                if allele_count['HPHP'] == num_alive_NMRs:
                    CONFIG['overall_stats']['hp']['pop_death_count'] += 1 
                    CONFIG['overall_stats']['hp']['pop_death_ticks'][run] = tick
                else:
                    CONFIG['overall_stats']['HP']['pop_death_count'] += 1 
                    CONFIG['overall_stats']['HP']['pop_death_ticks'][run] = tick

                
                if not CONFIG['show_r']:
                    CONFIG['overall_stats']['pop_death_ticks'][run] = tick
                    return True    
                
        if CONFIG['show_r'] and (allele_count['rr'] == num_alive_NMRs or allele_count['RR'] == num_alive_NMRs):

            # checks that both populations still exist
            if CONFIG['overall_stats']['r']['pop_death_ticks'][run] == CONFIG['overall_stats']['R']['pop_death_ticks'][run]:

                if allele_count['RR'] == num_alive_NMRs:
                    CONFIG['overall_stats']['r']['pop_death_count'] += 1
                    CONFIG['overall_stats']['r']['pop_death_ticks'][run] = tick
                else:
                    CONFIG['overall_stats']['R']['pop_death_count'] += 1
                    CONFIG['overall_stats']['R']['pop_death_ticks'][run] = tick
                
                if not CONFIG['show_hp']:
                    CONFIG['overall_stats']['pop_death_ticks'][run] = tick
                    return True

        if CONFIG['show_hp'] and CONFIG['show_r']:
            if CONFIG['overall_stats']['hp']['pop_death_ticks'][run] != CONFIG['overall_stats']['HP']['pop_death_ticks'][run] and CONFIG['overall_stats']['r']['pop_death_ticks'][run] != CONFIG['overall_stats']['R']['pop_death_ticks'][run]:
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
            distributions['change_per_x_unit'][pop_name] = [CONFIG['overall_stats'][pop_name][data_name][i] / sum(CONFIG['overall_stats'][pop_name][data_name][i:]) for i in range(len(CONFIG['overall_stats'][pop_name][data_name]))]

            # Save Hazard Function # To calculate it, we need to take -d log(S(t)) /dt, where S(t) is the number of individuals surviving at age t
            distributions['hazard'][pop_name] = [(log(sum(CONFIG['overall_stats'][pop_name][data_name][i:])) - log(sum(CONFIG['overall_stats'][pop_name][data_name][i+1:]))) / 1 for i in range(len(CONFIG['overall_stats'][pop_name][data_name]) - 1)]

        # Save PDF Data
        distributions['PDF'][pop_name] = [i / total_nmrs for i in CONFIG['overall_stats'][pop_name][data_name]]

        # Save CDF Data
        distributions['CDF'][pop_name] =  list(np.cumsum(distributions['PDF'][pop_name], dtype=float))
        distributions['CDF'][pop_name] = [float(i) for i in distributions['CDF'][pop_name]] # Stupid bit of code to fix numpy update bug

    else:
        print('Error: there were no recorded instances of', title, 'data for the population', pop_name)

  return distributions

def insult_severity(CONFIG, nmr=False, tick=False, after_pop_death = False, equilibrium_population_sizes = False):
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
        rand = random()*50000 
        insult_severity = 1000 / (rand**2)
    elif CONFIG['InsultDistribution'] == 'T8':
        init_hp = CONFIG[nmr.hp_pops['type_name']] + .1 # This is the biggest liberty I took, obviously dividing by log 1 is an issue but init_hp = either 1000 or 1200, I figured that .1 didn't make much difference there
        rand = random() # just a random number from 0 to 1
        alpha = .01 #The population size was too low with .1
        age = abs(nmr.birth - tick) + 1
        beta = -(log(nmr.hp_pops['value']/init_hp)) / age # where I have assumed that health_points(t) = init_hp * exp(-beta t), which defines the parameter beta;
        dt = 1 # you do an insult every dt in time, defining the parameter dt: An insult is applied every tick
        insult_severity = (init_hp * alpha * (exp(beta * dt) - 1)) / (beta * log(1/rand))  
        # So, alpha is arbitrary, but beta and dt are determined by the other parameters of your simulation, and the correct values need to be used otherwise it won't work.
    elif CONFIG['InsultDistribution'] == 'T8B':
        init_hp = CONFIG[nmr.hp_pops['type_name']] + .1 # This is the biggest liberty I took, obviously dividing by log 1 is an issue but init_hp = either 1000 or 1200, I figured that .1 didn't make much difference there
        rand = random() # just a random number from 0 to 1
        alpha = .01 #The population size was too low with .1
        age = abs(nmr.birth - tick) + 1
        beta =  -log(2) / nmr.r_pops['value'] # where I have assumed that health_points(t) = init_hp * exp(-beta t), which defines the parameter beta;
        insult_severity = (init_hp * alpha * (exp(beta ) - 1)) / (beta * log(1/rand))  
    elif CONFIG['InsultDistribution'] == 'T9': 
        rand = random()*50 
        insult_severity = 1000 / (rand**3)
    elif CONFIG['InsultDistribution'] == 'T99': 
        insult_severity = 5

    if CONFIG['LogorithmicPopulationCap']:
        if CONFIG['UseAvgForLogCap']:
          population_size = mean(CONFIG['last_3_pop_sizes'])
        else:
          population_size = len(CONFIG['alive_NMRS']['female']) + len(CONFIG['alive_NMRS']['male']) if not after_pop_death else CONFIG['InitialPopulation']
          insult_severity = insult_severity * np.exp(CONFIG['pop_cap_mult'] * (population_size**CONFIG['pop_cap_exponent']))

    if not 'insults' in CONFIG:
        CONFIG['insults'] = []
    CONFIG['insults'].append(insult_severity)  
      
    return insult_severity
