
from get import pop_death, random_parent
from update import init_population,apply_insults,NMR_birth,reset_overall_stats
from time import time
from random import random

def run_simulation(initial_pop_a_fraction, CONFIG):
    '''Runs simulation in a limited capacity'''
    
    if CONFIG['show_hp'] and CONFIG['show_r']:
      print('Warning, results may be off if both R and HP are being looked at')
      
    CONFIG['InitialRecessiveAlleleFraction'] = initial_pop_a_fraction
    print('When InitialRecessiveAlleleFraction = ', initial_pop_a_fraction)
    run_start_time = time()
    num_successful_runs = 0

    for run in range(CONFIG['NumRuns']):

        CONFIG['alive_NMRS'] = {'female': [], 'male': []}
        CONFIG['alive_NMRS'] = init_population(CONFIG, run)

        for tick in range(CONFIG['MaxTicks']):

            apply_insults(CONFIG['alive_NMRS']['male'], tick, CONFIG, run)
            apply_insults(CONFIG['alive_NMRS']['female'], tick, CONFIG, run)
            
            # if (random_parent(CONFIG['alive_NMRS']['male'], tick, CONFIG) == None or random_parent(CONFIG['alive_NMRS']['female'], tick, CONFIG) == None):
            #   if CONFIG['ShowIndividualRunStats']:
            #     print('Cant birth NMRs')
            #   break

            if pop_death(run, tick, CONFIG):
              num_successful_runs += 1
              break
            
            if (random_parent(CONFIG['alive_NMRS']['male'], tick, CONFIG) == None or random_parent(CONFIG['alive_NMRS']['female'], tick, CONFIG) == None):
                if CONFIG['ShowIndividualRunStats']:
                  print('Cant birth NMRs on tick ', tick)
            else:

              if CONFIG['UseLinearGrowth']:
                  if CONFIG['queen'].death:
                      CONFIG['queen'] = random_parent(CONFIG['alive_NMRS']['female'], tick, CONFIG)
                  elif random() < CONFIG['BirthProbability']:
                    if not CONFIG['StaticPopulationCap'] or (len(CONFIG['alive_NMRS']['female']) + len(CONFIG['alive_NMRS']['male'])) < CONFIG['InitialPopulation']:
                      dad = random_parent(CONFIG['alive_NMRS']['male'],  tick, CONFIG)
                      new_NMRs = NMR_birth(CONFIG, tick, CONFIG['queen'].w_pops['alleles'], dad.w_pops['alleles'], CONFIG['queen'].hp_pops['alleles'], dad.hp_pops['alleles'], CONFIG['queen'].r_pops['alleles'], dad.r_pops['alleles'])
                      # new_NMR.apply_insult(CONFIG, tick, run)
                      apply_insults(new_NMRs, tick, CONFIG, run)

              else:
                  for mom in CONFIG['alive_NMRS']['female']:
                      if CONFIG['TicksTillFertile'] <= tick - mom.birth and random() < CONFIG['BirthProbability']:
                        if not CONFIG['StaticPopulationCap'] or (len(CONFIG['alive_NMRS']['female']) + len(CONFIG['alive_NMRS']['male'])) < CONFIG['InitialPopulation']:
                          dad = random_parent(CONFIG['alive_NMRS']['male'], tick, CONFIG)
                          new_NMRs = NMR_birth(CONFIG, tick, mom.w_pops['alleles'], dad.w_pops['alleles'], mom.hp_pops['alleles'], dad.hp_pops['alleles'], mom.r_pops['alleles'], dad.r_pops['alleles'])
                          # new_NMR.apply_insult(CONFIG, tick, run)
                          apply_insults(new_NMRs, tick, CONFIG, run)


    print("That took", time() - run_start_time, 'seconds to run')
    run_start_time = time()
    print('num_successful_runs: ',num_successful_runs)

    if num_successful_runs  < 3:
      print('no successful runs')
      return -999, -1000
    elif CONFIG['show_hp'] and CONFIG['show_r']:
      print("CONFIG['overall_stats']['R']['pop_death_count'] / num_successful_runs: ", CONFIG['overall_stats']['R']['pop_death_count'] / num_successful_runs)
      print("CONFIG['overall_stats']['r']['pop_death_count'] / num_successful_runs: ", CONFIG['overall_stats']['r']['pop_death_count'] / num_successful_runs)
      print("CONFIG['overall_stats']['HP']['pop_death_count'] / num_successful_runs: ", CONFIG['overall_stats']['HP']['pop_death_count'] / num_successful_runs)
      print("CONFIG['overall_stats']['hp']['pop_death_count'] / num_successful_runs: ", CONFIG['overall_stats']['hp']['pop_death_count'] / num_successful_runs)
      return CONFIG['overall_stats']['HP']['pop_death_count'] / num_successful_runs, CONFIG['overall_stats']['hp']['pop_death_count'] / num_successful_runs,  CONFIG['overall_stats']['R']['pop_death_count'] / num_successful_runs, CONFIG['overall_stats']['r']['pop_death_count'] / num_successful_runs
    elif CONFIG['show_hp']:
      print("CONFIG['overall_stats']['HP']['pop_death_count'] / num_successful_runs: ", CONFIG['overall_stats']['HP']['pop_death_count'] / num_successful_runs)
      print("CONFIG['overall_stats']['hp']['pop_death_count'] / num_successful_runs: ", CONFIG['overall_stats']['hp']['pop_death_count'] / num_successful_runs)
      return CONFIG['overall_stats']['HP']['pop_death_count'] / num_successful_runs, CONFIG['overall_stats']['hp']['pop_death_count'] / num_successful_runs
    elif CONFIG['show_r']:
      print("CONFIG['overall_stats']['R']['pop_death_count'] / num_successful_runs: ", CONFIG['overall_stats']['R']['pop_death_count'] / num_successful_runs)
      print("CONFIG['overall_stats']['r']['pop_death_count'] / num_successful_runs: ", CONFIG['overall_stats']['r']['pop_death_count'] / num_successful_runs)
      return CONFIG['overall_stats']['R']['pop_death_count'] / num_successful_runs, CONFIG['overall_stats']['r']['pop_death_count'] / num_successful_runs

def find_equal_chance_initial_fraction(CONFIG):
    ''' Finds InitialRecessiveAlleleFraction which creates a 50-50 chance of either side winning '''

    best_pop_A_fraction = round(CONFIG['InitialRecessiveAlleleFraction'],2) # .5
    min_fraction = best_pop_A_fraction - .2 # 0
    max_fraction = best_pop_A_fraction + .2 # 1
    tolerance = 0.00005
    smallest_difference = 1
    retry_count = 0

    while True:
        mid_fraction = round((min_fraction + max_fraction) / 2, 15)
        if len(str(mid_fraction)) > 10 or mid_fraction < 0 or mid_fraction > 1:
              print('Error: This is getting too complicated')
              return best_pop_A_fraction

        deaths_HP, deaths_hp = run_simulation(mid_fraction, CONFIG)


        if abs(deaths_HP - deaths_hp) < smallest_difference:
          smallest_difference = abs(deaths_HP - deaths_hp)
          best_pop_A_fraction = mid_fraction
          

        reset_overall_stats(CONFIG)

        if deaths_HP == 'error' or abs(deaths_HP - deaths_hp) < tolerance or retry_count > 100:
            return best_pop_A_fraction
        elif deaths_HP < deaths_hp:
            min_fraction = mid_fraction
        elif deaths_hp < deaths_HP:
            max_fraction = mid_fraction
        else:
            print('lets try this again')
            retry_count += 1
