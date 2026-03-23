# Standard library imports
from collections import Counter
from math import ceil
from random import random
from statistics import mean, median
from time import strftime, time
from sys import argv

# Local imports
import display
import get
import update
from colony import Colony

#### Main Code ####

def track_colony_allele_frequencies(colonies):
    """
    Track allele frequencies for colony members, queens, and stored sperm 
    across all colonies in the landscape
    
    Returns:
        dict with 'W_allele', 'w_allele' 'total_alleles',
         'W_queen_allele', and 'w_queen_allele' frequencies
    """
    total_alleles_dict = {'W_allele': 0, 'w_allele': 0, 'total_alleles': 0,
                          'W_queen_allele': 0, 'w_queen_allele': 0}
    
    for colony in colonies:
        # Count Queen alleles
        if colony.queen is not None:
          queen_alleles = colony.queen.w_pops['alleles']
          for allele in queen_alleles:
              total_alleles_dict['total_alleles'] += 1
              if allele == 'W':  
                  total_alleles_dict['W_allele'] += 1
                  total_alleles_dict['W_queen_allele'] += 1
              else:  
                  total_alleles_dict['w_allele'] += 1
                  total_alleles_dict['W_queen_allele'] += 1

        # Count member alleles
        for member in colony.members:
            if member.death is not None:
                continue
            for allele in member.w_pops['alleles']:
                total_alleles_dict['total_alleles'] += 1
                if allele == 'W':  # HP-R
                    total_alleles_dict['W_allele'] += 1
                else:  # hp-r
                    total_alleles_dict['w_allele'] += 1
    
    return total_alleles_dict

def track_colony_times_fixed(colonies, run, tick, CONFIG):
    """
    Track fixation rates: number of colonies where the alleles have fixed
    
    Returns:
        dict with 'w_colonies_fixed' and 'W_colonies_fixed' and 'Not_fixed'
    """
    colony_fixation_counter = {'w_colonies_fixed': 0, 'W_colonies_fixed': 0, 'Not_fixed': 0}

    all_colonies_fixed = True
    first_col_w = None
    for colony in colonies:
          if colony.get_is_homozygous():
            # Check if colony has fixed alleles

            if first_col_w is None:
              first_col_w = colony.get_fixed_allele_name()
            else:
               if first_col_w != colony.get_fixed_allele_name():
                  all_colonies_fixed = False

            if colony.get_fixed_allele_name() == 'w':
              # Check if all alleles are w
              colony_fixation_counter['w_colonies_fixed'] += 1
            elif colony.get_fixed_allele_name() == 'W':
              # Check if all alleles are W
              colony_fixation_counter['W_colonies_fixed'] += 1
            else:
              print(f'Error: Allele type mismatch at tick {tick} for run {run}')
              if colony.queen is not None:
                print(f'Colony {colony.colony_id} has alleles {colony.queen.w_pops["alleles"]} and {colony.members[0].w_pops["alleles"]}')
              else:
                print(f'Colony {colony.colony_id} has no queen and {len(colony.members)} members')
              print(f'Colony {colony.colony_id} has fixed allele {colony.get_fixed_allele_name()}')
          else:
            all_colonies_fixed = False
            colony_fixation_counter['Not_fixed'] += 1
    
    if CONFIG['WriteDataToCSV']:
      display.figure4AtoB_csvs(CONFIG, tick, colony_fixation_counter['W_colonies_fixed'], colony_fixation_counter['w_colonies_fixed'], colony_fixation_counter['Not_fixed'], run)


    return all_colonies_fixed, colony_fixation_counter

CONFIG = get.configuration(list(argv))
CONFIG['overall_stats'] = update.reset_overall_stats(CONFIG)
CONFIG['num_NMRs'] = 0

# Initialize colony tracking if in colony mode
if CONFIG.get('UseColonyReproductionMode', False):
    CONFIG['colony_stats'] = {
        'colony_count': [],
        'num_runs_all_colonies_died': 0,
        'colony_deaths_per_run': [],
        'colony_age': [],
        'allele_frequencies': [],
        'fixation_rates': [],
        'colony_lifespan': {
          'Overall': [],
          'W': [],
          'w': []
        },
        'queen_lifespan': {
          'Overall': [],
          'W': [],
          'w': []
        },
        'num_nmrs_colony_contributed_to_mating_pool': {
          'Overall': [],
          'W': [],
          'w': []
        },
        'num_times_colony_contributed_to_mating_pool': {
          'Overall': [],
          'W': [],
          'w': []
        },
        'pct_of_colony_members_contributed_to_mating_pool_each_time': {
          'Overall': [],
          'W': [],
          'w': []
        },
    }

# Track ticks to reach equilibrium thresholds (50%, 75%, 85%) for Overall and HP and hp population
to_track_avg_thresholds = []
to_track_avg_thresholds += ['W'] if 'W' in CONFIG['equilibrium_population_size'] else []
to_track_avg_thresholds += ['w'] if 'w' in CONFIG['equilibrium_population_size'] else []
to_track_avg_thresholds += ['Overall'] if CONFIG['UseOverallRunStats'] and 'Overall' in CONFIG['equilibrium_population_size'] else []
equilibrium_average_health_points = [] 

for pop_name in CONFIG['display_factors']:
    # Initializes distribution arrays

    if CONFIG['GetLifespanDistribution']:
        CONFIG['overall_stats'][pop_name]['lifespan_distribution'] = [0] * CONFIG['MaxAge'] # Number of NMRs who died at each age - used to calculate age distribution

    if CONFIG['GetAgeDistribution'] :
        CONFIG['overall_stats'][pop_name]['age_distribution'] = [0] * CONFIG['MaxAge'] # Number of NMRs of each age at one time - used to calculate age distribution

    if CONFIG['GetMortality']:
        CONFIG['overall_stats'][pop_name]['equilibrium_lifespan_distribution'] = [0] * CONFIG['MaxAge'] # Number of NMRs who died at each age - used to calculate mortality


if CONFIG['StoreNMRData']:
  csv_name = 'csv/'+ argv[1]+ strftime('%b%d-%H')+'.csv'
  try:
    CONFIG['writer'] = open(csv_name, 'w')
    line =  'NMR ID, Age, Current Insults(oldest), Current Insults, Current Insults 1(most recent), Health Points, Birth Tick, Current Tick, run\n' 
    CONFIG['writer'].write(line)
  except (IOError, OSError) as e:
    print(f'Error: Failed to open CSV file {csv_name} for writing: {e}')
    CONFIG['StoreNMRData'] = False
  except IndexError as e:
    print(f'Error: Missing command line argument for CSV filename: {e}')
    CONFIG['StoreNMRData'] = False
  except Exception as e:
    print(f'Unexpected error opening CSV file {csv_name}: {e}')
    CONFIG['StoreNMRData'] = False

colony_has_disaster = []
colony_deaths_after_disaster = []
colony_deaths_not_after_disaster = []

for run in range(CONFIG['NumRuns']):

    run_start_time = time()
    
    CONFIG['alive_NMRS'] = dict(
        female = [],
        male = [],
        count = [],
    )
    
    equilibrium_population_size = dict()
    for pop_name in CONFIG['display_factors']:
      equilibrium_population_size[pop_name] = []
      CONFIG['max_experimental_age'][pop_name] = 0

    CONFIG['alive_NMRS']['ww'] = []
    CONFIG['alive_NMRS']['wW'] = []
    CONFIG['alive_NMRS']['WW'] = []

    if CONFIG['GetIndividualRunData']:
        population = []
    else:
        population = False

    # Initialize colonies or single population based on mode
    if CONFIG.get('UseColonyReproductionMode', False):
        colonies = update.init_colonies(CONFIG, run, population)
        next_colony_id = CONFIG.get('N0Colonies', 5)
        # Initialize colony death counter for this run
        colony_deaths_this_run = 0
        # Initialize alive_NMRS structure for tracking (used for some stats)
        CONFIG['alive_NMRS'] = dict(
            female=[],
            male=[],
            count=[],
            ww=[],
            wW=[],
            WW=[]
        )
        # Add all colony members to alive_NMRS for compatibility
        for colony in colonies:
            CONFIG['alive_NMRS']['female'].append(colony.queen)
            for member in colony.members:
                if member.sex == 'F':
                    CONFIG['alive_NMRS']['female'].append(member)
                else:
                    CONFIG['alive_NMRS']['male'].append(member)

    else:
        colonies = None
        next_colony_id = 0
        colony_deaths_this_run = 0
        CONFIG['alive_NMRS'] = update.init_population(CONFIG, run, population)
    
    # Initialize tracking for equilibrium thresholds (reset per run)
    equilibrium_thresholds_recorded = {
      'Overall': {'50pct': False, '75pct': False, '85pct': False, '100pct': False, '120pct': False}, 
      'W': {'50pct': False, '75pct': False, '85pct': False, '100pct': False, '120pct': False}, 
      'w': {'50pct': False, '75pct': False, '85pct': False, '100pct': False, '120pct': False}}
    
    equilibrium_health_points = []

    t = 0
    for tick in range(CONFIG['MaxTicks']):
        t = tick
        
        if CONFIG['WriteDataToCSV']:
          display.write_to_csv_each_tick(CONFIG, tick, run, colonies) 

        if CONFIG.get('UseColonyReproductionMode', False):
            # Colony reproduction mode
            
            # Conduct mating flight periodically
            if tick > 0 and tick % CONFIG.get('MatingFlightInterval', 10) == 0:
                fertilized_pioneer_candidates = update.conduct_mating_flight(colonies, tick, CONFIG, run, population)

                # Try to form new colonies from fertilized females
                num_colonies = len(colonies)
                for pioneer_group in fertilized_pioneer_candidates:
                    # Calculate formation probability
                    prob = update.calculate_colony_formation_probability(num_colonies, CONFIG)
                    if random() < prob:                        
                        # Form new colony
                        queen = pioneer_group[0]
                        colony_members = [] if len(pioneer_group) == 1 else pioneer_group[1:]
                        new_colony = update.form_new_colony(queen, tick, CONFIG, next_colony_id, colony_members, population)
                        colonies.append(new_colony)
                        next_colony_id += 1
                        num_colonies += 1
                    else:
                      # Pioneer candidates die if they fail to establish a colony
                      for pioneer in pioneer_group:
                        pioneer.kill_nmr(CONFIG, tick, run)
                 
            # Apply insults to all colony members and Queens
            for colony in colonies[:]:  # Use slice to avoid modification during iteration
                if not colony.is_alive(CONFIG, tick, run):
                  # Check if Queen or colony died (colony dies)
                  # Remove colony and kill all inhabitants
                  colonies.remove(colony)
                  colony_deaths_this_run += 1
                  continue
                else:

                  if random() < float(CONFIG['ColonyDisasterProbability']):
                    colony_has_disaster.append(True)
                    insult_multiplier = float(CONFIG['ColonyDisasterInsultMultiplier'])
                    # print(f'Colony disaster occurred at tick {tick} for run {run} and colony {colony.colony_id}')
                  else:
                    colony_has_disaster.append(False)
                    insult_multiplier = 1

                  deaths_this_time = 0
                  # Apply insults to members
                  for member in list(colony.members):
                      member.apply_insult(CONFIG, colony.get_size(), tick, run, insult_multiplier=insult_multiplier, is_queen=False)
                      # print(f'Member {member.id} applied insult {member.current_insults[-1]} with multiplier {insult_multiplier}')
                      # Remove dead members
                      if member.death is not None:
                          colony.remove_member(member)
                          if insult_multiplier > 1:
                            deaths_this_time += 1

                  if insult_multiplier > 1:
                    colony_deaths_after_disaster.append(deaths_this_time)
                  else:
                    colony_deaths_not_after_disaster.append(deaths_this_time)

                  # Queen produces offspring using stored sperm or male in colony
                  if random() < CONFIG['BirthProbability'] and colony.queen is not None:
                      # Chooses dad from stored sperm or heterozygous parent
                      dad_alleles = colony.queen.stored_sperm_alleles if not CONFIG['MatesWithinColony'] else get.random_parent_alleles(colony.members, tick, CONFIG, 'M')
                      if dad_alleles is not None:
                        new_NMRs = update.NMR_birth(CONFIG, tick, colony.queen.w_pops['alleles'], dad_alleles, population)
                        update.apply_insults(new_NMRs, tick, CONFIG, run, colony.get_size(), insult_multiplier=insult_multiplier)
                        # Add new NMRs to colony
                        for new_nmr in new_NMRs:
                            colony.add_member(new_nmr)
                  
                  if colony.queen is not None:
                    # Apply insults to Queen (with HP multiplier)
                    colony.queen.apply_insult(CONFIG, colony.get_size(), tick, run, insult_multiplier=insult_multiplier, is_queen=True)

                    if colony.queen.death is not None and CONFIG['TicksToQueenReplacement'] is not None:
                        # Track queen lifespan data
                        if CONFIG['UseOverallRunStats']:
                          CONFIG['colony_stats']['queen_lifespan']['Overall'].append(tick - colony.queen.birth)
                        CONFIG['colony_stats']['queen_lifespan'][colony.queen.w_pops['type_name']].append(tick - colony.queen.birth)
                        # Remove dead queen from colony
                        colony.queen = None
                        colony.last_queen_death_tick = tick

                  if colony.queen is None and CONFIG['TicksToQueenReplacement'] is not None:
                    if tick - colony.last_queen_death_tick >= int(CONFIG['TicksToQueenReplacement']):   
                      # Replace queen with a new one
                      colony.queen = get.random_parent(colony.get_females(), tick, CONFIG, 'F')
                      colony.remove_member(colony.queen)
                      colony.last_queen_death_tick = None

            if len(colonies) == 0:
                print(f'All colonies died out at tick {tick} for run {run}')
                CONFIG['colony_stats']['num_runs_all_colonies_died'] += 1
                if len(CONFIG['alive_NMRS']['female']) > 0 or len(CONFIG['alive_NMRS']['male']) > 0:
                    print(f'There are {len(CONFIG['alive_NMRS']['female'])} female NMRs and {len(CONFIG['alive_NMRS']['male'])} male NMRs left in alive nmrs')
                break
            
            
            # Track colony statistics and update allele counts for compatibility
            if CONFIG.get('UseColonyReproductionMode', False):
                CONFIG['colony_stats']['colony_count'].append(len(colonies))
                CONFIG['colony_stats']['colony_age'].append(mean([tick - colony.established_on for colony in colonies]))
                CONFIG['colony_stats']['allele_frequencies'].append(track_colony_allele_frequencies(colonies))
                all_colonies_fixed, fixation_rate = track_colony_times_fixed(colonies, run, tick, CONFIG)
                CONFIG['colony_stats']['fixation_rates'].append(fixation_rate)
                

                
                # Update allele counts for compatibility with existing stats tracking
                allele_count = Counter()

                # tmp_alive_nmrs = CONFIG['alive_NMRS']['female'][:] + CONFIG['alive_NMRS']['male'][:]
                for colony in colonies:
                  
                  for member in colony.members:
                      # tmp_alive_nmrs.remove(member)###
                      if member.death is None:
                          member_allele_name = member.w_pops['allele_name']
                          allele_count[member_allele_name] += 1
                      else:
                          print(f'Warning: Member {member.id} is dead {member.sex}')
                  
                  if colony.queen is not None and colony.queen.death is None:
                        queen_allele_name = colony.queen.w_pops['allele_name']
                        allele_count[queen_allele_name] += 1
                        # print(f'Queen {colony.queen} is {colony.queen.w_pops["allele_name"]}')
                        # tmp_alive_nmrs.remove(colony.queen)####
                if all_colonies_fixed:
                   if CONFIG['ShowIndividualRunStats']:
                     print (f'All colonies have the same fixed alleles at tick {tick} for run {run}')
                  #  break
                  #  print (f'All colonies have the same fixed alleles at tick {tick} for run {run}')

  
                  # #  TODO: Get rid of this, it is a temporary fix to avoid the population from getting stuck at a fixed allele
                  #  for nmr in tmp_alive_nmrs:
                  #   nmr.kill_nmr(CONFIG, tick, run)       

            current_pop = len(CONFIG['alive_NMRS']['female']) + len(CONFIG['alive_NMRS']['male'])
            
        else:
            # Original single population mode
            update.apply_insults(CONFIG['alive_NMRS']['male'], tick, CONFIG, run, len(CONFIG['alive_NMRS']['female']) + len(CONFIG['alive_NMRS']['male']))
            update.apply_insults(CONFIG['alive_NMRS']['female'], tick, CONFIG, run, len(CONFIG['alive_NMRS']['female']) + len(CONFIG['alive_NMRS']['male']))

          # Generate Births
            if get.random_parent_alleles(CONFIG['alive_NMRS']['male'], tick, CONFIG, 'M') == None or get.random_parent_alleles(CONFIG['alive_NMRS']['female'], tick, CONFIG, 'F') == None:
                if CONFIG['ShowIndividualRunStats']:
                  print('Can not birth NMRs on tick ', tick, '. There are ', len(CONFIG['alive_NMRS']['male']) + len(CONFIG['alive_NMRS']['female']),'NMRs total')
            else:
              if CONFIG['UseLinearGrowth']:
                  if not CONFIG['UseImmortalQueenMode'] and (not CONFIG['queen'] or not CONFIG['queen'].death == None):
                      CONFIG['queen'] = get.random_parent(CONFIG['alive_NMRS']['female'], tick, CONFIG, 'F')
                  elif random() < CONFIG['BirthProbability']:
                    if not CONFIG['StaticPopulationCap'] or (len(CONFIG['alive_NMRS']['female']) + len(CONFIG['alive_NMRS']['male'])) < CONFIG['InitialPopulation']:
                    # Choses random father
                      dad_alleles = get.random_parent_alleles(CONFIG['alive_NMRS']['male'], tick, CONFIG, 'M')
                      if dad_alleles is not None:
                        new_NMRs = update.NMR_birth(CONFIG, tick, CONFIG['queen'].w_pops['alleles'], dad_alleles, population)
                        update.apply_insults(new_NMRs, tick, CONFIG, run, len(CONFIG['alive_NMRS']['female']) + len(CONFIG['alive_NMRS']['male']))
              else:
                  for mom in CONFIG['alive_NMRS']['female']:
                    if not CONFIG['StaticPopulationCap'] or (len(CONFIG['alive_NMRS']['female']) + len(CONFIG['alive_NMRS']['male'])) < CONFIG['InitialPopulation']:
                      if CONFIG['TicksTillFertile'] <= tick - mom.birth \
                              and random() < CONFIG['BirthProbability']:

                        # Choses random father
                          dad_alleles = get.random_parent_alleles(CONFIG['alive_NMRS']['male'], tick, CONFIG, 'M')
                          if dad_alleles is not None:
                            new_NMRs = update.NMR_birth(CONFIG, tick, mom.w_pops['alleles'], dad_alleles, population)
                            update.apply_insults(new_NMRs, tick, CONFIG, run, len(CONFIG['alive_NMRS']['female']) + len(CONFIG['alive_NMRS']['male']))


        allele_count = Counter(elem.w_pops['alleles'][0] + elem.w_pops['alleles'][1] for elem in CONFIG['alive_NMRS']['female']) + \
                       Counter(elem.w_pops['alleles'][0] + elem.w_pops['alleles'][1] for elem in CONFIG['alive_NMRS']['male'])

        # print(f'Allele count: {allele_count}')
    # Set population Counts
        if CONFIG['UseAvgForLogCap']:
          CONFIG['last_3_pop_sizes'] = CONFIG['last_3_pop_sizes'][1:] + [len(CONFIG['alive_NMRS']['female']) + len(CONFIG['alive_NMRS']['male']) ]

          print(CONFIG['last_3_pop_sizes'])

        CONFIG['alive_NMRS']['ww'].insert(tick, allele_count['ww'])
        CONFIG['alive_NMRS']['WW'].insert(tick, allele_count['WW'])
        CONFIG['alive_NMRS']['wW'].insert(tick, allele_count['Ww'] + allele_count['wW'])

        CONFIG['alive_NMRS']['count'].insert(tick, len(CONFIG['alive_NMRS']['female']) + len(CONFIG['alive_NMRS']['male']))

        if len(CONFIG['alive_NMRS']['female']) + len(CONFIG['alive_NMRS']['male']) > CONFIG['MaxPopulation']:
            print ('The maximum population has been exceeded at ', tick, ' ticks on run', run)
            break
        if ((len(CONFIG['alive_NMRS']['female']) < 1 or len(CONFIG['alive_NMRS']['male']) < 1)):
            lost_pop = 'female' if len(CONFIG['alive_NMRS']['female']) < 1 else 'male'
            print ('There are no',lost_pop, 'NMRs left at', tick, ' ticks on run', run)
            break
        if CONFIG['overall_stats']['pop_death_ticks'][run] == 0:
          get.pop_death(run, tick, CONFIG, colonies)

        elif not CONFIG['GetAgeDistribution'] and not CONFIG['GetMortality'] and not CONFIG['GetEquilibriumPopulation']:
          break
        else:
          winning_pops = ['Overall'] if CONFIG['UseOverallRunStats'] else []
          if len(CONFIG['alive_NMRS']['female']) > 0:
            winning_pops += [CONFIG['alive_NMRS']['female'][0].w_pops['type_name']]

          if CONFIG['GetAgeDistribution']:
            for population_name in winning_pops:
                update.age_and_mortality_dist(CONFIG, population_name, CONFIG['alive_NMRS']['female'], tick, 'age_distribution')
                update.age_and_mortality_dist(CONFIG, population_name, CONFIG['alive_NMRS']['male'], tick, 'age_distribution')

          if CONFIG['GetEquilibriumPopulation']:
            for population_name in winning_pops:
              equilibrium_population_size[population_name].append(len(CONFIG['alive_NMRS']['male']) + len(CONFIG['alive_NMRS']['female']))
        

        current_pop = len(CONFIG['alive_NMRS']['female']) + len(CONFIG['alive_NMRS']['male'])

        # Track ticks to reach equilibrium thresholds (50%, 75%, 85%, 120%) for Overall and HP and hp population
    
        for population_name in to_track_avg_thresholds:

          if equilibrium_thresholds_recorded[population_name]['100pct']:
            if equilibrium_thresholds_recorded[population_name]['120pct']:
              # At equilibrium population size, what is the average number of Health Points in the population?
              total_health_points = sum(elem.current_health_points for elem in CONFIG['alive_NMRS']['female']) + sum(elem.current_health_points for elem in CONFIG['alive_NMRS']['male'])
              # total_health_points = sum(elem.current_health_points for elem in CONFIG['alive_NMRS']['female'] if elem.hp_type_name == pop) + sum(elem.current_health_points for elem in CONFIG['alive_NMRS']['male'] if elem.hp_type_name == pop)
              equilibrium_health_points.append(total_health_points)
            elif len(CONFIG['overall_stats'][population_name]['ticks_to_100pct_equilibrium']) > 0 and \
              tick == (CONFIG['overall_stats'][population_name]['ticks_to_100pct_equilibrium'][-1] + ceil(0.2 * CONFIG['equilibrium_population_size'][population_name])):

              stat_name = '100plus20'
              CONFIG['overall_stats'][population_name]['ticks_to_100plus20pct_equilibrium'].append(tick)
              equilibrium_thresholds_recorded[population_name]['120pct'] = True

              if population_name == 'Overall':
                fhp_female_count = len(CONFIG['alive_NMRS']['female']) * 2
                fhp_male_count = len(CONFIG['alive_NMRS']['male']) if CONFIG['UseImmortalQueenMode'] else len(CONFIG['alive_NMRS']['male']) * 2               # In immortal queen mode, the males are considered haploid and the females are diploid so we need to adjust the counts accordingly

              else:
                fhp_female_count = sum(elem.w_pops['alleles'].count(population_name) for elem in CONFIG['alive_NMRS']['female'])
                fhp_male_count = sum(1 for elem in CONFIG['alive_NMRS']['male'] if elem.w_pops['alleles'][0] == population_name) if CONFIG['UseImmortalQueenMode'] else sum(elem.w_pops['alleles'].count(population_name) for elem in CONFIG['alive_NMRS']['male'])               # In immortal queen mode, the males are considered haploid and the females are diploid so we need to adjust the counts accordingly

              
              CONFIG['overall_stats'][population_name]['NMRs_at_equilibrium_at']['female']['percent_'+stat_name].append(fhp_female_count)
              CONFIG['overall_stats'][population_name]['NMRs_at_equilibrium_at']['male']['percent_'+stat_name].append(fhp_male_count)
              CONFIG['overall_stats'][population_name]['NMRs_at_equilibrium_at']['Overall']['percent_'+stat_name].append(fhp_male_count + fhp_female_count)

          else:   
              # Check for pre-calculated equilibrium thresholds
              if CONFIG['equilibrium_population_size'].get(population_name):
                  # Check each threshold (only record once per run)
                  for stat_name in ['50', '75', '85', '100']:
                    if not equilibrium_thresholds_recorded[population_name][stat_name+'pct'] and current_pop >= (int(stat_name) / 100) * CONFIG['equilibrium_population_size'].get(population_name):
                      CONFIG['overall_stats'][population_name]['ticks_to_'+stat_name+'pct_equilibrium'].append(tick)
                      equilibrium_thresholds_recorded[population_name][stat_name+'pct'] = True
                      if population_name == 'Overall':
                        fhp_female_count = len(CONFIG['alive_NMRS']['female'])
                        fhp_male_count = len(CONFIG['alive_NMRS']['male'])
                      else:
                        fhp_female_count = sum(1 for elem in CONFIG['alive_NMRS']['female'] if elem.hp_type_name == population_name)
                        fhp_male_count = sum(1 for elem in CONFIG['alive_NMRS']['male'] if elem.hp_type_name == population_name)
                      
                      CONFIG['overall_stats'][population_name]['NMRs_at_equilibrium_at']['female']['percent_'+stat_name].append(fhp_female_count)
                      CONFIG['overall_stats'][population_name]['NMRs_at_equilibrium_at']['male']['percent_'+stat_name].append(fhp_male_count)
                      CONFIG['overall_stats'][population_name]['NMRs_at_equilibrium_at']['Overall']['percent_'+stat_name].append(fhp_male_count + fhp_female_count)
              else:
                print(f'Warning: No equilibrium population size found for {population_name}, re-run with GetEquilibriumPopulation True to calculate equilibrium thresholds then this feature will be available')
 
    if CONFIG['overall_stats']['pop_death_ticks'][run] >= CONFIG['MaxTicks'] - 1:
      winner = 'MaxTicksExceeded'
    elif CONFIG['overall_stats']['W']['pop_death_ticks'][run] > 0:
      if get.get_longevity_allele(CONFIG) == 'W':
        winner = 'Vitality'
      else:
        winner = 'Longevity'
    elif CONFIG['overall_stats']['w']['pop_death_ticks'][run] > 0:
      if get.get_longevity_allele(CONFIG) == 'W':
        winner = 'Longevity'
      else:
        winner = 'Vitality'
    else:
      winner = 'Other'

    if t >= CONFIG['MaxTicks'] - 1 and winner == 'Other':
      print(f'WARNING: Max ticks exceeded without fixation on run {run}, consider increasing max ticks')



    if CONFIG['WriteDataToCSV']:
      display.suppfig_3F_and_4D_csv(CONFIG, winner, run)

    if CONFIG['ShowIndividualRunStats']:
        print ('\nRUN NUMBER: ', run)
        print ('Runs took ', time() - run_start_time,
               ' seconds to run\n')
        print('Num ticks: ', tick)

        if len(CONFIG['alive_NMRS']['female']) > 0:
          print(f'If it fixed, the winning population was {CONFIG["alive_NMRS"]["female"][0].w_pops["alleles"][0]}')

    if CONFIG['GetIndividualRunData']:
      display.individual_run_stats(CONFIG, equilibrium_population_size, population, run, t)
      t = 0

    if CONFIG['ShowIndividualRunGraphs']:
      
      display.individual_run_graph(CONFIG, run, tick, 'population')
      if CONFIG['UseColonyReproductionMode']:
        display.individual_run_graph(CONFIG, run, tick, 'colony')
            
    if len(equilibrium_health_points) > 0:
      equilibrium_average_health_points.append(mean(equilibrium_health_points))
    
    # Store colony deaths for this run
    if CONFIG.get('UseColonyReproductionMode', False):
        CONFIG['colony_stats']['colony_deaths_per_run'].append(colony_deaths_this_run)
      
        
          
# Showing overall statistics data
print ('\n------- Overall Statistics -------\n')

print (CONFIG['NumRuns'], ' runs took ', time() - CONFIG['overall_stats']['start_time'], ' seconds to run\n')


dom_info = 'The gomperz curve of population W had a half life of ', CONFIG['R'], ' and began with ', CONFIG['HP'], ' health points'
rec_info = 'The gomperz curve of population w had a half life of ', CONFIG['r'], ' and began with ', CONFIG['hp'], ' health points'
display.simulation_stats('W','w', dom_info, rec_info, CONFIG)

if len(CONFIG['overall_stats']['pop_death_ticks']) > 0:
    print ('\nOn average a run lasted ', mean(CONFIG['overall_stats']['pop_death_ticks']), 'ticks\n')

if CONFIG['GetIndividualRunData']:
    display.run_stats_using_individual_data(CONFIG)
    
    # Print ticks until equilibrium threshold statistics        
    print('\n--- Equilibrium Population Threshold Statistics ---')

    longevity = get.get_longevity_allele(CONFIG)
    vitality = 'w' if longevity == 'W' else 'W'
    print(f'Longevity is {longevity} and Vitality is {vitality}')
    
    for stat_name in ['50', '75', '85', '100', '100plus20']:
      
      for population_name in to_track_avg_thresholds:
          
        if len(CONFIG['overall_stats'][population_name]['ticks_to_'+stat_name+'pct_equilibrium']) > 0:
          avg_ticks_to_pct_equilibrium = mean(CONFIG['overall_stats'][population_name]['ticks_to_'+stat_name+'pct_equilibrium'])
          print(f'Average ticks to reach {stat_name}% of equilibrium using the {population_name} population equilibrium size: {avg_ticks_to_pct_equilibrium:.2f} (from {len(CONFIG["overall_stats"][population_name]["ticks_to_"+stat_name+"pct_equilibrium"])} runs)')

        print()

      for gender in ['female', 'male', 'Overall']:
        if len(CONFIG['overall_stats'][longevity]['NMRs_at_equilibrium_at'][gender]['percent_'+stat_name]) > 0 and len(CONFIG['overall_stats'][vitality]['NMRs_at_equilibrium_at'][gender]['percent_'+stat_name]) > 0 and mean(CONFIG['overall_stats'][vitality]['NMRs_at_equilibrium_at'][gender]['percent_'+stat_name]) != 0:
          print(f'The {gender} ratio of L to V alleles at {stat_name}% of equilibrium is {mean(CONFIG['overall_stats'][longevity]['NMRs_at_equilibrium_at'][gender]['percent_'+stat_name]) / mean(CONFIG['overall_stats'][vitality]['NMRs_at_equilibrium_at'][gender]['percent_'+stat_name])} or {mean(CONFIG['overall_stats'][longevity]['NMRs_at_equilibrium_at'][gender]['percent_'+stat_name])} : { mean(CONFIG['overall_stats'][vitality]['NMRs_at_equilibrium_at'][gender]['percent_'+stat_name])}')
          print()
    
    print('--------------------------------')

    if CONFIG['UseOverallRunStats'] and len(CONFIG['overall_stats']['W']['avg_ls']) > 0 and len(CONFIG['overall_stats']['w']['avg_ls']) > 0 and len(CONFIG['overall_stats']['Overall']['avg_ls']) > 0 and mean(CONFIG['overall_stats']['Overall']['avg_ls']) != 0:
        print ('The average lifespan of population W was ', 100
                * (mean(CONFIG['overall_stats']['W']['avg_ls'])
                - mean(CONFIG['overall_stats']['w']['avg_ls']))
                / mean(CONFIG['overall_stats']['Overall']['avg_ls']), ' percent longer')

    
if CONFIG['GetLifespanDistribution']:
  if CONFIG['MakeOverallGraphs']:
    lifespan_distribution = get.distribution_data('lifespan_distribution', 'Distribution of Lifespans', CONFIG)
    display.distribution_data( lifespan_distribution, 
                                'lifespan_distribution',
                                'Distribution of Lifespans',
                                CONFIG,
                                # compare_types = ['PDF','change_per_x_unit', 'hazard'],
                                # graph_types=['PDF'],
                                # graph_types = ['PDF','change_per_x_unit', 'hazard'],
                                # write_types = ['PDF','change_per_x_unit', 'hazard'],
                                # write_types=['PDF'],
                                # print_types = ['hazard'],
                                compare_types=['PDF'],
                                # graph_types=['PDF'],
                                graph_types=['hazard','PDF'],
                                write_types=['hazard','PDF'],
                                setAxis = False,
                              )
  
if CONFIG['GetMortality']:
  if CONFIG['MakeOverallGraphs']:
    equilibrium_lifespan_distribution = get.distribution_data('equilibrium_lifespan_distribution', 'Equilibrium Distribution of Lifespans', CONFIG)
    display.distribution_data( equilibrium_lifespan_distribution, 
                                'equilibrium_lifespan_distribution',
                                'Change in Mortality',
                                CONFIG,
                                # compare_types = ['PDF','change_per_x_unit', 'hazard'],
                                compare_types=['PDF'],
                                # graph_types = ['PDF','change_per_x_unit', 'hazard', 'CDF', 'Absolute'],
                                graph_types=['PDF'],
                                write_types=['PDF'],
                                setAxis = False
                              )

if CONFIG['GetAgeDistribution']:
  if CONFIG['MakeOverallGraphs']:
    age_distribution = get.distribution_data('age_distribution', 'Age Distribution', CONFIG)
    display.distribution_data( age_distribution,
                              'age_distribution', 
                              'Age Distribution',
                              CONFIG,
                              # graph_types = ['PDF', 'hazard'],
                              # graph_types = ['PDF', 'Absolute'],
                              # print_types = [ 'CDF'],
                              # write_types = ['CDF'],
                              compare_types=['PDF'],
                              graph_types=['PDF'],
                              write_types=['PDF','CDF'],
                              setAxis = False,
                              # compare_types = ['PDF']
                              )

      
    display.write_log_comparison_to_csv(CONFIG, age_distribution)
# print(CONFIG)

print('\n\nmax_age rec = '+    str(CONFIG['max_experimental_age']['w'] ))
print('max_age dom = '+    str(CONFIG['max_experimental_age']['W'] ))

if CONFIG['UseOverallRunStats']:
  print('max_age overall = '+str(CONFIG['max_experimental_age']['Overall'] ))

run_type = 'Linear' if CONFIG['UseLinearGrowth'] else 'Exponential'
header = "Population,	Allele,	Trial,	HealthPoints,	HalfLife,	Fixation,	Lifespan-Avg, Lifespan-Max,	Probability"
line = run_type + ',HPR,1,' +  str(CONFIG['HP']) + ',' + str(CONFIG['R'])+ ',' + \
  str(CONFIG['overall_stats']['w']['pop_death_count'])+ ',' + \
    str(mean(CONFIG['overall_stats']['W']['avg_ls']))+ ','+\
        str(CONFIG['max_experimental_age']['W']) + ',' + \
          str(CONFIG['overall_stats']['w']['pop_death_count'] / CONFIG['NumRuns']) + '\n'
display.write_to_file(line, 'data/fixation-data.csv', header)
line = run_type + ',hpr,1,' +  str(CONFIG['hp']) + ',' + str(CONFIG['r'])+ ',' + \
  str(CONFIG['overall_stats']['W']['pop_death_count'])+ ',' + \
    str(mean(CONFIG['overall_stats']['w']['avg_ls']))+ ','+ \
      str(CONFIG['max_experimental_age']['w']) + ',' + \
        str(CONFIG['overall_stats']['W']['pop_death_count'] / CONFIG['NumRuns']) + '\n'
display.write_to_file(line, 'data/fixation-data.csv', header)

# Print colony reproduction mode statistics
if CONFIG.get('UseColonyReproductionMode', False):
    print('\n------- Colony Reproduction Mode Statistics -------\n')
    if 'colony_stats' in CONFIG and len(CONFIG['colony_stats']['colony_count']) > 0:
        final_colony_count = CONFIG['colony_stats']['colony_count'][-1]
        avg_colony_count = mean(CONFIG['colony_stats']['colony_count'])
        print(f'Final colony count: {final_colony_count}')
        print(f'Average colony count: {avg_colony_count:.2f}')
        print(f'Average colony age: {mean(CONFIG['colony_stats']['colony_age']):.2f}')

        for factor in CONFIG['display_factors']:
          print()
          if len(CONFIG['colony_stats']['queen_lifespan'][factor]) > 0:
            print(f'The average {factor} queen lived {mean(CONFIG['colony_stats']['queen_lifespan'][factor])} ticks')
            print(f'The average {factor} colony lived {mean(CONFIG['colony_stats']['colony_lifespan'][factor])} ticks')
            print(f'The average {factor} colony contributed {mean(CONFIG['colony_stats']['num_nmrs_colony_contributed_to_mating_pool'][factor])} nmrs to the mating pool')
            print(f'The average {factor} colony contributed to the mating pool {mean(CONFIG['colony_stats']['num_times_colony_contributed_to_mating_pool'][factor])} times')
            print(f'The average {factor} colony contributed {mean(CONFIG['colony_stats']['pct_of_colony_members_contributed_to_mating_pool_each_time'][factor])} percent of their colony to the mating pool each time')

        if len(CONFIG['colony_stats']['allele_frequencies']) > 0:
        #  This has the following info: 'W_allele', 'w_allele' 'total_alleles', 'W_queen_allele', and 'w_queen_allele' frequencies
            w_alleles = [a['w_allele'] for a in CONFIG['colony_stats']['allele_frequencies']]
            W_alleles = [a['W_allele'] for a in CONFIG['colony_stats']['allele_frequencies']]
            
            print(f'\nMean allele frequencies:')
            print(f'  W alleles: {mean(W_alleles):.4f}')
            print(f'  w allele: {mean(w_alleles):.4f}')
        

        if len(CONFIG['colony_stats']['fixation_rates']) > 0:
            w_colonies_fixed = [a['w_colonies_fixed'] for a in CONFIG['colony_stats']['fixation_rates']]
            W_colonies_fixed = [a['W_colonies_fixed'] for a in CONFIG['colony_stats']['fixation_rates']]
            not_fixed = [a['Not_fixed'] for a in CONFIG['colony_stats']['fixation_rates']]
            total_colonies = [w_colonies_fixed[t] + W_colonies_fixed[t] + not_fixed[t] for t in range(len(CONFIG['colony_stats']['fixation_rates']))]
            print(f'\nMean fixation rates:')
            print(f'  Homozygous w: {mean(w_colonies_fixed):.4f}')
            print(f'  Homozygous W: {mean(W_colonies_fixed):.4f}')
            print(f'  Heterozygous: {mean(not_fixed):.4f}')
            print(f'\nPercent fixation rates:')
            print(f'  Homozygous w: {mean([w_colonies_fixed[t]/total_colonies[t] for t in range(len(total_colonies))]):.4f}')
            print(f'  Homozygous W: {mean([W_colonies_fixed[t]/total_colonies[t] for t in range(len(total_colonies))]):.4f}')
            print(f'  Heterozygous: {mean([not_fixed[t]/total_colonies[t] for t in range(len(total_colonies))]):.4f}')

    print(f'\nPercentage of runs where all colonies died: {CONFIG["colony_stats"]["num_runs_all_colonies_died"] / CONFIG["NumRuns"]}')
    if len(CONFIG['colony_stats']['colony_deaths_per_run']) > 0:
        avg_colony_deaths = mean(CONFIG['colony_stats']['colony_deaths_per_run'])
        print(f'Average number of colony deaths per run: {avg_colony_deaths:.2f}')
        if len(CONFIG['overall_stats']['pop_death_ticks']) == len(CONFIG['colony_stats']['colony_deaths_per_run']) and avg_colony_deaths != 0:
          avg_run_ticks =  mean(CONFIG['overall_stats']['pop_death_ticks'])
          print(f'Average number of ticks per colony death: {avg_run_ticks / avg_colony_deaths}')

        # Get indices (run numbers) where W fixed (all WW)
        W_fixed_runs = [i for i, tick in enumerate(CONFIG['overall_stats']['w']['pop_death_ticks']) if tick != 0]
        if len(W_fixed_runs) > 0:
          num_W_win_colony_deaths = [CONFIG['colony_stats']['colony_deaths_per_run'][i] for i in W_fixed_runs]
          print(f'\nAverage number of colony deaths per run where W fixed: {mean(num_W_win_colony_deaths)}')
          W_avg_ticks_to_fixation = sum(CONFIG['overall_stats']['w']['pop_death_ticks']) / len(W_fixed_runs)
          print(f'Average number of ticks per colony death where W fixed: {W_avg_ticks_to_fixation / mean(num_W_win_colony_deaths)}')

        # Get indices (run numbers) where w fixed (all ww)
        w_fixed_runs = [i for i, tick in enumerate(CONFIG['overall_stats']['W']['pop_death_ticks']) if tick != 0]
        if len(w_fixed_runs) > 0:
          num_w_win_colony_deaths = [CONFIG['colony_stats']['colony_deaths_per_run'][i] for i in w_fixed_runs]
          print(f'\nAverage number of colony deaths per run where w fixed: {mean(num_w_win_colony_deaths)}')
          w_avg_ticks_to_fixation = sum(CONFIG['overall_stats']['W']['pop_death_ticks']) / len(w_fixed_runs)
          print(f'Average number of ticks per colony death where w fixed: {w_avg_ticks_to_fixation / mean(num_w_win_colony_deaths)}')

if len(colony_has_disaster) > 0:
  disasters = [i for i in colony_has_disaster if i]
  print(f'Percentage of ticks where a colony had a disaster: {len(disasters) / len(colony_has_disaster)}')

  if len(colony_deaths_after_disaster) > 0:
    print(f'Average number of colony deaths per run after a disaster: {mean(colony_deaths_after_disaster)}')
  if len(colony_deaths_not_after_disaster) > 0:
    print(f'Average number of colony deaths per run not after a disaster: {mean(colony_deaths_not_after_disaster)}')
  # print('avg insult',mean(CONFIG['insults']))
  # print('avg old insult',mean(CONFIG['old_insults']))
else:
  print('No disasters occurred')