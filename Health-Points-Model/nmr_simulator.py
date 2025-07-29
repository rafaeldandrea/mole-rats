
### Imports ###
from os.path import isfile
from random import random
from statistics import mean
from time import time, strftime
from collections import Counter
import sys
import display
import get 
import allele_fraction_finder
import update 


#### Main Code ####

CONFIG = get.configuration(list(sys.argv))
CONFIG['overall_stats'] = update.reset_overall_stats(CONFIG)
CONFIG['num_NMRs'] = 0


for pop_name in CONFIG['display_factors']:
    # Initializes distribution arrays

    if CONFIG['GetLifespanDistribution']:
        CONFIG['overall_stats'][pop_name]['lifespan_distribution'] = [0] * CONFIG['MaxAge'] # Number of NMRs who died at each age - used to calculate age distribution

    if CONFIG['GetAgeDistribution'] :
        CONFIG['overall_stats'][pop_name]['age_distribution'] = [0] * CONFIG['MaxAge'] # Number of NMRs of each age at one time - used to calculate age distribution

    if CONFIG['GetMortality']:
        CONFIG['overall_stats'][pop_name]['equilibrium_lifespan_distribution'] = [0] * CONFIG['MaxAge'] # Number of NMRs who died at each age - used to calculate mortality

if CONFIG['GetPopAFraction']:
    # Find the value of InitialRecessiveAlleleFraction where the chance of population HP dying out is the same as the chance of population hp dying out
    equal_chance_fraction = allele_fraction_finder.find_equal_chance_initial_fraction(CONFIG)
    print("InitialRecessiveAlleleFraction for equal chances:", equal_chance_fraction)
    new_init_pop = "\nCONFIG['InitialRecessiveAlleleFraction'] = " + str(equal_chance_fraction) + "  # This was found in " + str(CONFIG['NumRuns']) + 'runs'
    display.write_to_file(msg = new_init_pop,
              file_path = CONFIG['file_path'],
              msg_type = 'new initial population')

    CONFIG['InitialRecessiveAlleleFraction'] = equal_chance_fraction
    print("I'll Prove it below")
    

if CONFIG['StoreNMRData']:
  csv_name = 'csv/'+ sys.argv[1]+ strftime('%b%d-%H')+'.csv'
  CONFIG['writer'] = open(csv_name,'w')
  line =  'NMR ID,' + 'Age,' +'Current Insults(oldest),Current Insults ,Current Insults 1(most recent),' + 'Health Points, Birth Tick, Current Tick, run\n' 
  CONFIG['writer'].write(line)

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

    if CONFIG['show_hp']:
      CONFIG['alive_NMRS']['HPHP'] = []
      CONFIG['alive_NMRS']['hpHP'] = []
      CONFIG['alive_NMRS']['hphp'] = []

    if CONFIG['show_r']:
      CONFIG['alive_NMRS']['rr'] = []
      CONFIG['alive_NMRS']['RR'] = []
      CONFIG['alive_NMRS']['rR'] = []

    if CONFIG['UseTradeoffAllele']:
      CONFIG['alive_NMRS']['ww'] = []
      CONFIG['alive_NMRS']['wW'] = []
      CONFIG['alive_NMRS']['WW'] = []

    if CONFIG['GetIndividualRunData']:
        population = []
    else:
        population = False

    CONFIG['alive_NMRS'] = update.init_population(CONFIG, run, population)
    t = 0
    for tick in range(CONFIG['MaxTicks']):
        t = tick
        update.apply_insults(CONFIG['alive_NMRS']['male'], tick, CONFIG, run)
        update.apply_insults(CONFIG['alive_NMRS']['female'], tick, CONFIG, run)

      # Generate Births
        if (get.random_parent(CONFIG['alive_NMRS']['male'], tick, CONFIG) == None or get.random_parent(CONFIG['alive_NMRS']['female'], tick, CONFIG) == None):
            if CONFIG['ShowIndividualRunStats']:
              print('Cant birth NMRs on tick ', tick, '. There are ', len(CONFIG['alive_NMRS']['male']) + len(CONFIG['alive_NMRS']['female']),'NMRs total')
        else:
          if CONFIG['UseLinearGrowth']:
              if not CONFIG['queen'] or not CONFIG['queen'].death == None:
                  # print(len(CONFIG['alive_NMRS']['female']))
                  CONFIG['queen'] = get.random_parent(CONFIG['alive_NMRS']['female'], tick, CONFIG)
              elif random() < CONFIG['BirthProbability']:
                if not CONFIG['StaticPopulationCap'] or (len(CONFIG['alive_NMRS']['female']) + len(CONFIG['alive_NMRS']['male'])) < CONFIG['InitialPopulation']:
                # Choses random father
                  dad = get.random_parent(CONFIG['alive_NMRS']['male'], tick, CONFIG)
                  new_NMRs = update.NMR_birth(CONFIG, tick, CONFIG['queen'].w_pops['alleles'], dad.w_pops['alleles'], CONFIG['queen'].hp_pops['alleles'], dad.hp_pops['alleles'], CONFIG['queen'].r_pops['alleles'], dad.r_pops['alleles'], population)
                  # new_NMR.apply_insult(CONFIG, tick, run)
                  update.apply_insults(new_NMRs, tick, CONFIG, run)

          else:
              for mom in CONFIG['alive_NMRS']['female']:
                if not CONFIG['StaticPopulationCap'] or (len(CONFIG['alive_NMRS']['female']) + len(CONFIG['alive_NMRS']['male'])) < CONFIG['InitialPopulation']:
                  if CONFIG['TicksTillFertile'] <= tick - mom.birth \
                          and random() < CONFIG['BirthProbability']:

                    # Choses random father
                      dad = get.random_parent(CONFIG['alive_NMRS']['male'], tick, CONFIG)
                      new_NMRs = update.NMR_birth(CONFIG, tick, mom.w_pops['alleles'], dad.w_pops['alleles'], mom.hp_pops['alleles'], dad.hp_pops['alleles'], mom.r_pops['alleles'], dad.r_pops['alleles'], population)
                      # new_NMR.apply_insult(CONFIG, tick, run)
                      update.apply_insults(new_NMRs, tick, CONFIG, run)


        allele_count = Counter(elem.hp_pops['alleles'][0] + elem.hp_pops['alleles'][1] for elem in CONFIG['alive_NMRS']['female']) + \
                       Counter(elem.r_pops['alleles'][0] + elem.r_pops['alleles'][1] for elem in CONFIG['alive_NMRS']['female']) + \
                       Counter(elem.hp_pops['alleles'][0] + elem.hp_pops['alleles'][1] for elem in CONFIG['alive_NMRS']['male']) + \
                       Counter(elem.r_pops['alleles'][0] + elem.r_pops['alleles'][1] for elem in CONFIG['alive_NMRS']['male']) + \
                       Counter(elem.w_pops['alleles'][0] + elem.w_pops['alleles'][1] for elem in CONFIG['alive_NMRS']['female']) + \
                       Counter(elem.w_pops['alleles'][0] + elem.w_pops['alleles'][1] for elem in CONFIG['alive_NMRS']['male'])

    # Set population Counts
        if CONFIG['UseAvgForLogCap']:
          CONFIG['last_3_pop_sizes'] = CONFIG['last_3_pop_sizes'][1:] + [len(CONFIG['alive_NMRS']['female']) + len(CONFIG['alive_NMRS']['male']) ]

          print(CONFIG['last_3_pop_sizes'])
        if CONFIG['show_hp']:
          CONFIG['alive_NMRS']['HPHP'].insert(tick, allele_count['HPHP'])
          CONFIG['alive_NMRS']['hphp'].insert(tick, allele_count['hphp'])
          CONFIG['alive_NMRS']['hpHP'].insert(tick, allele_count['HPhp'] + allele_count['hpHP'])

        if CONFIG['show_r']:
          CONFIG['alive_NMRS']['RR'].insert(tick, allele_count['RR'])
          CONFIG['alive_NMRS']['rr'].insert(tick, allele_count['rr'])
          CONFIG['alive_NMRS']['rR'].insert(tick, allele_count['Rr'] + allele_count['rR'])

        if CONFIG['UseTradeoffAllele']:
          CONFIG['alive_NMRS']['ww'].insert(tick, allele_count['ww'])
          CONFIG['alive_NMRS']['WW'].insert(tick, allele_count['WW'])
          CONFIG['alive_NMRS']['wW'].insert(tick, allele_count['Ww'] + allele_count['wW'])

        CONFIG['alive_NMRS']['count'].insert(tick, len(CONFIG['alive_NMRS']['female']) + len(CONFIG['alive_NMRS']['male']))

        if len(CONFIG['alive_NMRS']['female']) + len(CONFIG['alive_NMRS']['male']) > CONFIG['MaxPopulation']:
            print ('The maximum population has been exceeded at ', tick, ' ticks')
            break
        if len(CONFIG['alive_NMRS']['female']) < 1 or len(CONFIG['alive_NMRS']['male']) < 1:
            lost_pop = 'female' if len(CONFIG['alive_NMRS']['female']) < 1 else 'male'
            print ('There are no',lost_pop, 'NMRs left at', tick, ' ticks')
            break
        if CONFIG['overall_stats']['pop_death_ticks'][run] == 0:
          get.pop_death(run, tick, CONFIG)
        elif not CONFIG['GetAgeDistribution'] and not CONFIG['GetMortality'] and not CONFIG['GetEquilibriumPopulation']:
          break
        else:
          winning_pops = ['Overall'] if CONFIG['UseOverallRunStats'] else []
          if len(CONFIG['alive_NMRS']['female']) > 0:
            if CONFIG['show_r']:
              winning_pops += [CONFIG['alive_NMRS']['female'][0].r_pops['type_name']]
            if CONFIG['show_hp']:
              winning_pops += [CONFIG['alive_NMRS']['female'][0].hp_pops['type_name']]

          if CONFIG['GetAgeDistribution']:
            for pop in winning_pops:
                update.age_and_mortality_dist(CONFIG, pop, CONFIG['alive_NMRS']['female'], tick, 'age_distribution')
                update.age_and_mortality_dist(CONFIG, pop, CONFIG['alive_NMRS']['male'], tick, 'age_distribution')

          if CONFIG['GetEquilibriumPopulation']:
            for pop in winning_pops:
              equilibrium_population_size[pop].append(len(CONFIG['alive_NMRS']['male']) + len(CONFIG['alive_NMRS']['female']))



    if CONFIG['ShowIndividualRunStats']:
        print ('\nRUN NUMBER: ', run)
        print ('Runs took ', time() - run_start_time,
               ' seconds to run\n')
        print('Num ticks: ', tick)

        if CONFIG['UseTradeoffAllele']:
          display.print_allele_counts(CONFIG, 'w')
          # if CONFIG['alive_NMRS']['ww'][-1] == 0:
          #   print('No Williams alleles')
          # elif CONFIG['alive_NMRS']['WW'][-1] == 0:
          #   print('No Williams alleles')  
          # else:
          #   print('no winners')
          # totalwAlleles = [(2*CONFIG['alive_NMRS']['ww'][i]) + CONFIG['alive_NMRS']['wW'][i] for i in range(len(CONFIG['alive_NMRS']['wW']))]
          # totalWAlleles = [(2*CONFIG['alive_NMRS']['WW'][i]) + CONFIG['alive_NMRS']['wW'][i] for i in range(len(CONFIG['alive_NMRS']['wW']))]
          # print('Recessive williams alleles(w): ', totalwAlleles)
          # print('Dominant williams alleles(W): ', totalWAlleles)

        else:
          if CONFIG['show_hp']:
            totalhpAlleles = [2*CONFIG['alive_NMRS']['hphp'][i] + CONFIG['alive_NMRS']['hpHP'][i] for i in range(len(CONFIG['alive_NMRS']['hpHP']))],
            totalHPAlleles = [2*CONFIG['alive_NMRS']['HPHP'][i] + CONFIG['alive_NMRS']['hpHP'][i] for i in range(len(CONFIG['alive_NMRS']['hpHP']))],
            print('Recessive HP alleles(hp): ', totalhpAlleles)
            print('Dominant HP alleles(HP): ', totalHPAlleles)

          if CONFIG['show_r']:

            totalRAlleles = [2*CONFIG['alive_NMRS']['RR'][i] + CONFIG['alive_NMRS']['rR'][i] for i in range(len(CONFIG['alive_NMRS']['rR']))],
            totalrAlleles = [2*CONFIG['alive_NMRS']['rr'][i] + CONFIG['alive_NMRS']['rR'][i] for i in range(len(CONFIG['alive_NMRS']['rR']))],
            print('Recessive half-life alleles(r): ', totalrAlleles)
            print('Dominant half-life alleles(R): ', totalRAlleles)


    if CONFIG['GetIndividualRunData']:
        # print('At the end of the run there were ', len(CONFIG['alive_NMRS']['male'])+ len(CONFIG['alive_NMRS']['female']), 'nmrs alive')
        display.individual_run_stats(CONFIG, equilibrium_population_size, population, run, t)
        t = 0

    if CONFIG['ShowIndividualRunGraphs']:
        display.individual_run_graph(CONFIG, run)
        
          
# Showing overall statistics data
print ('\n------- Overall Statistics -------\n')

print (CONFIG['NumRuns'], ' runs took ', time() - CONFIG['overall_stats']['start_time'], ' seconds to run\n')


if CONFIG['show_hp']:
    dom_info = 'Population HP began with ', CONFIG['HP'], ' health points\n'
    rec_info = 'Population hp began with ', CONFIG['hp'], ' health points'
    display.simulation_stats('HP','hp', dom_info, rec_info, CONFIG)

if CONFIG['show_r']:

    dom_info = 'The gomperz curve of population R had a half life of ', CONFIG['R']
    rec_info = 'The gomperz curve of population r had a half life of ', CONFIG['r']
    display.simulation_stats('R','r', dom_info, rec_info, CONFIG)

if len(CONFIG['overall_stats']['pop_death_ticks']) > 0:
    print ('\nOn average a run lasted ', mean(CONFIG['overall_stats']['pop_death_ticks']), 'ticks\n')

if CONFIG['GetIndividualRunData']:
    display.run_stats_using_individual_data(CONFIG)
    if CONFIG['show_hp']:
        if CONFIG['UseOverallRunStats']:
            print ('The average lifespan of population HP was ', 100
                    * (mean(CONFIG['overall_stats']['HP']['avg_ls'])
                    - mean(CONFIG['overall_stats']['hp']['avg_ls']))
                    / mean(CONFIG['overall_stats']['Overall']['avg_ls']), ' percent longer')

    
    if CONFIG['show_r']:
        if CONFIG['UseOverallRunStats']:
            print ('The average lifespan of population R was ', 100
                    * (mean(CONFIG['overall_stats']['R']['avg_ls'])
                    - mean(CONFIG['overall_stats']['r']['avg_ls']))
                    / mean(CONFIG['overall_stats']['Overall']['avg_ls']), ' percent longer')
    
if CONFIG['GetLifespanDistribution']:
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
                              graph_types=['PDF'],
                              write_types=['PDF'],
                              setAxis = False,
                            )
  
if CONFIG['GetMortality']:
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

print('\n\nmax_age overall = '+str(CONFIG['max_experimental_age']['Overall'] ))
print('max_age rec = '+    str(CONFIG['max_experimental_age']['hp'] ))
print('max_age dom = '+    str(CONFIG['max_experimental_age']['HP'] ))

if CONFIG['show_hp'] and CONFIG['show_r']:

  run_type = 'Linear' if CONFIG['UseLinearGrowth'] else 'Exponential'
  header = "Population,	Allele,	Trial,	HealthPoints,	HalfLife,	Fixation,	Lifespan-Avg, Lifespan-Max,	Probability"
  line = run_type + ',HPR,1,' +  str(CONFIG['HP']) + ',' + str(CONFIG['R'])+ ',' + \
    str(CONFIG['overall_stats']['hp']['pop_death_count'])+ ',' + \
      str(mean(CONFIG['overall_stats']['HP']['avg_ls']))+ ','+\
          str(CONFIG['max_experimental_age']['R']) + ',' + \
            str(CONFIG['overall_stats']['hp']['pop_death_count'] / CONFIG['NumRuns']) + '\n'
  display.write_to_file(line, 'data/fixation-data.csv', header)
  line = run_type + ',hpr,1,' +  str(CONFIG['hp']) + ',' + str(CONFIG['r'])+ ',' + \
    str(CONFIG['overall_stats']['HP']['pop_death_count'])+ ',' + \
      str(mean(CONFIG['overall_stats']['hp']['avg_ls']))+ ','+ \
        str(CONFIG['max_experimental_age']['r']) + ',' + \
          str(CONFIG['overall_stats']['HP']['pop_death_count'] / CONFIG['NumRuns']) + '\n'
  display.write_to_file(line, 'data/fixation-data.csv', header)

# print('avg insult',mean(CONFIG['insults']))
# print('avg old insult',mean(CONFIG['old_insults']))
