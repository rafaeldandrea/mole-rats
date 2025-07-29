from time import asctime, time
from random import choices, random, choice
from Naked_Mole_Rat import Naked_Mole_Rat


### FUNCTIONS ###
def reset_overall_stats(CONFIG):
  '''Initializes overall_stats dict so that run data is clean'''
  
  CONFIG['overall_stats'] = dict(
      start_time = time(),
      pop_death_ticks = [0] * CONFIG['NumRuns'], # Number of ticks it took for a population to die each run
  )

  print('Starting Run at ', asctime())

  if CONFIG['UseOverallRunStats']:
     CONFIG['overall_stats']['Overall'] = dict(
        pop_death_count = 0, # Number of runs which the population died out
        avg_ls = [], # Average lifespan per run
        pop_death_ticks = [0] * CONFIG['NumRuns'], # Number of ticks it took for a population to die each run
        avg_equilibrium_population_size = [], # Average population size after equilibrium for each run that reaches equalibrium
      )
     if CONFIG['GetLifespanDistribution']:
        CONFIG['overall_stats']['Overall']['all_ls'] = [], # all lifespans

     
  if CONFIG['show_hp']:
      CONFIG['overall_stats']['HP'] =  dict(
        pop_death_count = 0, # Number of runs which the population died out
        avg_ls = [], # Average lifespan per run for the population
        pop_death_ticks = [0] * CONFIG['NumRuns'], # Number of ticks it took for a population to die each run
        avg_equilibrium_population_size = [], # Average population size after equilibrium for each run that reaches equalibrium
      )
      CONFIG['overall_stats']['hp'] =  dict(
        pop_death_count = 0, # Number of runs which the population died out
        avg_ls = [], # Average lifespan per run for the population
        pop_death_ticks = [0] * CONFIG['NumRuns'], # Number of ticks it took for a population to die each run
        avg_equilibrium_population_size = [], # Average population size after equilibrium for each run that reaches equalibrium
      )
      if CONFIG['GetLifespanDistribution']:
          CONFIG['overall_stats']['hp']['all_ls'] = [], # all lifespans
          CONFIG['overall_stats']['HP']['all_ls'] = [], # all lifespans

  if CONFIG['show_r']:
      CONFIG['overall_stats']['R'] =  dict(
        pop_death_count = 0, # Number of runs which the population died out
        avg_ls = [], # Average lifespan per run for the population      
        pop_death_ticks = [0] * CONFIG['NumRuns'], # Number of ticks it took for a population to die each run
        avg_equilibrium_population_size = [], # Average population size after equilibrium for each run that reaches equalibrium
      )
      CONFIG['overall_stats']['r'] =  dict(
        pop_death_count = 0, # Number of runs which the population died out
        avg_ls = [], # Average lifespan per run for the population
        pop_death_ticks = [0] * CONFIG['NumRuns'], # Number of ticks it took for a population to die each run
        avg_equilibrium_population_size = [], # Average population size after equilibrium for each run that reaches equalibrium
      )
      if CONFIG['GetLifespanDistribution']:
          CONFIG['overall_stats']['r']['all_ls'] = [], # all lifespans
          CONFIG['overall_stats']['R']['all_ls'] = [], # all lifespans
  return CONFIG['overall_stats']

def age_and_mortality_dist(CONFIG, pop_name, population, death_or_tick, dist_name):
    '''Updates values for age or mortality distribution'''
    for nmr in population:
      CONFIG['overall_stats'][pop_name][dist_name][abs(death_or_tick - nmr.birth) ] += 1

def NMR_birth(CONFIG, cur_tick, mom_w_alleles, dad_w_alleles, mom_hp_alleles, dad_hp_alleles, mom_r_alleles, dad_r_alleles, population = False, use_litter_size = True):
    '''Initializes new NMR'''
    new_NMRs = []
    litter_size = CONFIG['LitterSize'] if use_litter_size else 1
    for l in range(litter_size):
    # Choses diploid allele to be inherited from each parent
      parent_alleles = dict()
      
      if CONFIG['UseTradeoffAllele']:
        parent_alleles['mW'] = mom_w_alleles[choice([0,1])]
        parent_alleles['dW'] = dad_w_alleles[choice([0,1])]
      else:
        if CONFIG['show_hp']:
          parent_alleles['mHP'] = mom_hp_alleles[choice([0,1])]
          parent_alleles['dHP'] = dad_hp_alleles[choice([0,1])]
          
        if CONFIG['show_r']:
          parent_alleles['mR'] = mom_r_alleles[choice([0,1])]
          parent_alleles['dR'] = dad_r_alleles[choice([0,1])]

    # In this case only one NMR is born per litter
      new_NMR = Naked_Mole_Rat(CONFIG['num_NMRs'], cur_tick, parent_alleles, CONFIG)
      
      if population != False and CONFIG['GetIndividualRunData']:
          population.append(new_NMR)

    # add new NMRs to list of alive NMRs

      if new_NMR.sex == 'F':
        CONFIG['alive_NMRS']['female'].append(new_NMR)
      else:
        CONFIG['alive_NMRS']['male'].append(new_NMR)

      CONFIG['num_NMRs'] += 1
      new_NMRs.append(new_NMR)
    return new_NMRs
  
def set_random_NMR_age(active_pop_name, new_nmr, CONFIG, run):
    ''' Uses Binary Search To find a random NMR Age on the provided age distribution and ages nmr to that age '''
    
    # Uses Binary Search To find a random NMR Age on the provided age distribution
    random_num = random()
    lo, hi = 0, len(CONFIG['age_distribution']['CDF'][active_pop_name]) - 1
    best_ind = lo
    while lo <= hi:
        mid = int(lo + (hi - lo) // 2)
        if CONFIG['age_distribution']['CDF'][active_pop_name][mid] < random_num:
            lo = mid + 1
        elif CONFIG['age_distribution']['CDF'][active_pop_name][mid] > random_num:
            hi = mid - 1
        else:
            best_ind = mid
            break
          
        # check if mid is closer to val than best_ind
        if abs(CONFIG['age_distribution']['CDF'][active_pop_name][mid] - random_num) < abs(CONFIG['age_distribution']['CDF'][active_pop_name][best_ind] - random_num):
            best_ind = mid

    random_age = best_ind

    new_nmr.age_NMR(CONFIG, run, curr_tick=0, age=random_age)

def init_population(CONFIG, run, population = False):
    '''Initializes population to config specifications'''
    
    if CONFIG['UseAvgForLogCap']:
      CONFIG['last_3_pop_sizes'] = [CONFIG['InitialPopulation']] * 3
    
    while ((len(CONFIG['alive_NMRS']['female']) + len(CONFIG['alive_NMRS']['male'])) < CONFIG['InitialPopulation'] or len(CONFIG['alive_NMRS']['female']) < 1 or len(CONFIG['alive_NMRS']['male']) < 1):

        # Randomizes Parent Alleles
        mom_w_alleles = choices(['w','W'], weights=[CONFIG['InitialRecessiveAlleleFraction'], 1-CONFIG['InitialRecessiveAlleleFraction']], k = 2) if CONFIG['UseTradeoffAllele'] else ['w','w']
        dad_w_alleles = choices(['w','W'], weights=[CONFIG['InitialRecessiveAlleleFraction'], 1-CONFIG['InitialRecessiveAlleleFraction']], k = 2) if CONFIG['UseTradeoffAllele'] else ['w','w']

        mom_hp_alleles = choices(['hp','HP'], weights=[CONFIG['InitialRecessiveAlleleFraction'], 1-CONFIG['InitialRecessiveAlleleFraction']], k = 2) if CONFIG['show_hp'] else ['hp','hp']
        dad_hp_alleles = choices(['hp','HP'], weights=[CONFIG['InitialRecessiveAlleleFraction'], 1-CONFIG['InitialRecessiveAlleleFraction']], k = 2) if CONFIG['show_hp'] else ['hp','hp']
                  
        mom_r_alleles = choices(['r','R'], weights=[CONFIG['InitialRecessiveAlleleFraction'], 1-CONFIG['InitialRecessiveAlleleFraction']], k = 2) if CONFIG['show_r'] else ['r', 'r']
        dad_r_alleles = choices(['r','R'], weights=[CONFIG['InitialRecessiveAlleleFraction'], 1-CONFIG['InitialRecessiveAlleleFraction']], k = 2) if CONFIG['show_r'] else ['r', 'r']

        new_nmrs = NMR_birth(CONFIG, 0,  mom_w_alleles, dad_w_alleles, mom_hp_alleles, dad_hp_alleles, mom_r_alleles, dad_r_alleles, population, use_litter_size = False)
        # print(mom_w_or_hp_alleles)#Remove

        # Sets the age using the corresponding age distribution - please note, this is not implemented for use on both population types simultaneously
        for new_nmr in new_nmrs:
          if not CONFIG['GetAgeDistribution'] :

              # Gets the type that will be used to determine which age distribution is being used
              if not CONFIG['show_hp']:
                  pop_name = new_nmr.r_pops['type_name']
              elif not CONFIG['show_r']:
                  pop_name = new_nmr.hp_pops['type_name']
              else:
                  # In the case that both the allele types are in play, one is chosen at random
                  pop_name = choice([new_nmr.hp_pops['type_name'], new_nmr.r_pops['type_name']])
              
              set_random_NMR_age(pop_name, new_nmr, CONFIG, run)

          else:
            # Ages all NMRs to reproductive age without insults
            new_nmr.birth = -CONFIG['TicksTillFertile']

        # print((len(CONFIG['alive_NMRS']['female']) + len(CONFIG['alive_NMRS']['male'])) < CONFIG['InitialPopulation'] and (len(CONFIG['alive_NMRS']['female']) < 1 or len(CONFIG['alive_NMRS']['male']) < 1))

        
    # print("CONFIG['InitialPopulation']",CONFIG['InitialPopulation'])
    # print((len(CONFIG['alive_NMRS']['female']) + len(CONFIG['alive_NMRS']['male'])) < CONFIG['InitialPopulation'] and (len(CONFIG['alive_NMRS']['female']) < 1 or len(CONFIG['alive_NMRS']['male']) < 1))
    # print((len(CONFIG['alive_NMRS']['female']) + len(CONFIG['alive_NMRS']['male'])) < CONFIG['InitialPopulation'])
    # print( (len(CONFIG['alive_NMRS']['female']) < 1 or len(CONFIG['alive_NMRS']['male']) < 1))

    if CONFIG['UseLinearGrowth']:
        # print("len(CONFIG['alive_NMRS']['female'])",len(CONFIG['alive_NMRS']['female']))
        # print("len(CONFIG['alive_NMRS']['male'])",len(CONFIG['alive_NMRS']['male']))
        CONFIG['queen'] = choice(CONFIG['alive_NMRS']['female'])

    # print('end of init')#Remove

    return CONFIG['alive_NMRS']


def apply_insults(nmr_list, insult_tick, CONFIG, run):
  '''Obtains and applies insults for each NMR in list'''
  for nmr in list(nmr_list):
      nmr.apply_insult(CONFIG, insult_tick, run)