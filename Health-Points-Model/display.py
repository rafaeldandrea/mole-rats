# Standard library imports
from datetime import datetime
from importlib import import_module
from math import exp, log
from os.path import isfile
from statistics import mean

# Third-party imports
from matplotlib.pyplot import close, figtext, savefig, show
from numpy import trim_zeros
from pandas import DataFrame

# Local imports
from get import get_run_mode, get_longevity_allele, get_longevity_or_vitality_allele

def run_stats_using_individual_data(CONFIG):
    ''' Shows overall run stats only collected using individual run data '''
    
    if CONFIG['GetEquilibriumPopulation']:
      data_setup = "\n\nCONFIG['InitialPopulation'] = " 
      if CONFIG['UseOverallRunStats'] and len(CONFIG['overall_stats']['Overall']['avg_equilibrium_population_size']) > 0:
          data_setup += str(round(mean(CONFIG['overall_stats']['Overall']['avg_equilibrium_population_size']))) 
          CONFIG['InitialPopulation'] = round(mean(CONFIG['overall_stats']['Overall']['avg_equilibrium_population_size']))
      else:
          data_setup += "'-'" 
          CONFIG['InitialPopulation'] = '-'

      write_to_file(msg = data_setup,
                    file_path = CONFIG['file_path'],
                    msg_type = 'Failed to population size to file')
                  
    for pop_name in CONFIG['display_factors']:
        if len(CONFIG['overall_stats'][pop_name]['avg_ls']) > 0:
          print ('The average lifespan of a NMR in population ', pop_name, ' was ',
                    mean(CONFIG['overall_stats'][pop_name]['avg_ls']))

          data_setup = "\nCONFIG['avg_lifespan']['"+pop_name+"'] = " + str(mean(CONFIG['overall_stats'][pop_name]['avg_ls'])) 
          write_to_file(msg = data_setup,
                        file_path = CONFIG['file_path'],
                        msg_type = 'Failed to average lifespan to file')
        else:
          print('No average lifespan data exists')

        if CONFIG['GetEquilibriumPopulation'] and len(CONFIG['overall_stats'][pop_name]['avg_equilibrium_population_size']) > 0:

            avg_pop = mean(CONFIG['overall_stats'][pop_name]['avg_equilibrium_population_size'])
            print('The average equilibrium population size in population ',pop_name,' was ', avg_pop)

            data_setup = "\nCONFIG['" + 'equilibrium_population_size' + "']['" + pop_name + "'] = " + str(avg_pop)
            write_to_file(msg = data_setup,
                          file_path = CONFIG['file_path'],
                          msg_type = 'age distribution')

def simulation_stats(dom_name, rec_name, dom_info, rec_info, CONFIG):
    ''' Prints overall simulation stats '''

    print (
        rec_name, ' was a recessive allele, that began as ', 
        float(CONFIG['InitialRecessiveAlleleFraction']) * 100, ' percent of the population\n',
        rec_info
    )
    print (
        dom_name, ' was a dominant allele, that began as ', 
        abs(1 - float(CONFIG['InitialRecessiveAlleleFraction'])) * 100, ' percent of the population\n',
        dom_info
    )


    print (
    'In ',
    CONFIG['NumRuns'],
    ' runs, population ', rec_name, ' fixed ',
    CONFIG['overall_stats'][dom_name]['pop_death_count'],
    " times. That's ",
    100 * CONFIG['overall_stats'][dom_name]['pop_death_count'] / CONFIG['NumRuns'],
    ' percent of the time.',
    )

    print (
          'In ',
          CONFIG['NumRuns'],
          ' runs, population ', dom_name, ' fixed ',
          CONFIG['overall_stats'][rec_name]['pop_death_count'],
          " times. That's ",
          100 * CONFIG['overall_stats'][rec_name]['pop_death_count'] / CONFIG['NumRuns'],
          ' percent of the time.',
          )
    
    # for pop in [rec_name, dom_name]:
    #   print (
    #       'In ',
    #       CONFIG['NumRuns'],
    #       ' runs, population ', pop, ' died out ',
    #       CONFIG['overall_stats'][pop]['pop_death_count'],
    #       " times. That's ",
    #       100 * CONFIG['overall_stats'][pop]['pop_death_count'] / CONFIG['NumRuns'],
    #       ' percent of the time.',
    #       )
      
    print (
      'In ',
      CONFIG['NumRuns'],
      ' runs, neither population died out ',
      abs(CONFIG['NumRuns'] - CONFIG['overall_stats'][rec_name]['pop_death_count'] - CONFIG['overall_stats'][dom_name]['pop_death_count']),
      " times. That's ",
      100 * abs(CONFIG['NumRuns'] - CONFIG['overall_stats'][rec_name]['pop_death_count'] - CONFIG['overall_stats'][dom_name]['pop_death_count']) / CONFIG['NumRuns'],
      ' percent of the time.',
      )

    try:
      CONFIG['pcnt_run_failures'] = abs(CONFIG['NumRuns'] - CONFIG['overall_stats'][rec_name]['pop_death_count'] - CONFIG['overall_stats'][dom_name]['pop_death_count']) / CONFIG['NumRuns']
    except ZeroDivisionError:
      print('Warning: NumRuns is zero, cannot calculate run failure percentage')
      CONFIG['pcnt_run_failures'] = 0
    except KeyError as e:
      print(f'Warning: Missing key in CONFIG for run failure calculation: {e}')
      CONFIG['pcnt_run_failures'] = 0
    msg = "\nCONFIG['pcnt_run_failures'] = " + str(CONFIG['pcnt_run_failures'])
    write_to_file(msg = msg,
              file_path = CONFIG['file_path'],
              msg_type = 'distribution')

  
def individual_run_stats(CONFIG, equilibrium_population_size, population, run, end_tick):
  ''' Prints individual run stats and stores additional insight into overall trends for use in run_stats_using_individual_data()'''
  for pop_name in CONFIG['display_factors']:

    ages = []
    insults = []
    deaths = [0] * (CONFIG['MaxTicks'] + CONFIG['MaxAge'])

    if CONFIG['GetEquilibriumPopulation'] and len(equilibrium_population_size[pop_name]) > 0:
        CONFIG['overall_stats'][pop_name]['avg_equilibrium_population_size'].append( mean(equilibrium_population_size[pop_name]))

    # ages_at_tick_2000 = [0] *  CONFIG['MaxAge']
    
    if not 'old_insults' in CONFIG:
      CONFIG['old_insults'] = []    
    
    for nmr in population:
        if pop_name == 'Overall' or nmr.w_pops['type_name'] == pop_name:
          if nmr.death == None:
            nmr.age_NMR(CONFIG, run, end_tick, colony_population_size=len(CONFIG['alive_NMRS']['female']) + len(CONFIG['alive_NMRS']['male']))
          if nmr.death > 0:
            total_insult_on_death = sum(nmr.current_insults)
            ages.append(abs(abs(nmr.death - nmr.birth-1)))
            insults.append(total_insult_on_death)
            deaths[nmr.death] += 1
            CONFIG['max_experimental_age'][pop_name] = CONFIG['max_experimental_age'][pop_name] if nmr.death-nmr.birth < CONFIG['max_experimental_age'][pop_name] else nmr.death-nmr.birth

        # if nmr.death and nmr.death >= 2000 and nmr.birth <= 2000:
        #   nmr_age  =  2000 - nmr.birth 
          # ages_at_tick_2000[nmr_age] += 1 


    if len(ages) > 0:

      CONFIG['overall_stats'][pop_name]['avg_ls'].append(mean(ages))
      if CONFIG['GetLifespanDistribution']:
        if not isinstance(CONFIG['overall_stats'][pop_name]['all_ls'], list):
          CONFIG['overall_stats'][pop_name]['all_ls'] = list(CONFIG['overall_stats'][pop_name]['all_ls']) + ages if len(CONFIG['overall_stats'][pop_name]['all_ls']) > 1 else ages
        else:
          CONFIG['overall_stats'][pop_name]['all_ls'].extend(ages)

      # # Generating ages file
      # if pop_name == 'W' or pop_name == 'w':
      #   ages_at_tick_2000 = trim_zeros(ages_at_tick_2000, 'b')
      #   header = 'age, count, allele, run, tick number(2000), type\n'

      #   for x in range(len(ages_at_tick_2000)):
      #     run_type = 'Linear' if CONFIG['UseLinearGrowth'] else 'Exponential'
      #     line = str(x)+ ',' + str(ages_at_tick_2000[x])+ ',' + pop_name+ ',' + str(run) + ',2000,' + run_type + '\n'
      #     line = header + line if run == 0 and x==0 else line
      #     write_to_file(line, 'data/ages-data.csv', "writing ages to file", header)


def write_to_file(msg, file_path, msg_type = 'message', header = ''):
  """Write message to file with proper error handling."""
  try: 
    if not isfile(file_path):
      with open(file_path, "x") as myfile:
          myfile.write(header + '\n')
          myfile.write(msg)
    else:
      with open(file_path, "a") as myfile:
        myfile.write(msg)
  except (IOError, OSError) as e:
      print(f'Unexpected error appending {msg_type} to file {file_path}: {e}')

      

  
def distribution_data(current_distribution, 
                              data_name,
                              title,
                              CONFIG,
                              graph_types = [],
                              print_types = [],
                              write_types = [],
                              compare_types = [],
                              setAxis = True
                              ):
  
  '''Writes, prints or graphs distibutions'''

  # print("# Now displaying the", title, "data which was just generated using the configuration ",CONFIG,"\n")

  # Writes Data to Config File
  for pop_name in CONFIG['display_factors']:
    for data_type in list(set(write_types) | set(print_types) | set(compare_types)):
      if pop_name in current_distribution[data_type]:
          
        data_setup =  str(list(current_distribution[data_type][pop_name])) + " # Data Generated in " + str(CONFIG['NumRuns']) + " Runs"

        if data_type in print_types:
          print("CONFIG['" + data_name + "']['"+data_type+"']['" + pop_name + "'] = " + data_setup)
        if data_type in write_types:
          writen_data_setup = '\n' + "CONFIG['" + data_name + "']['"+data_type+"']['" + pop_name + "'] = " + data_setup 
          write_to_file(msg = writen_data_setup,
                        file_path = CONFIG['file_path'],
                        msg_type = 'distribution')

        if data_type in compare_types and (CONFIG['file_path'].endswith('Lin.py') or CONFIG['file_path'].endswith('Exp.py')) :
          if 'Exp.py' in CONFIG['file_path']:# or 'Lin.py' in CONFIG['file_path']:
            comparison_data_setup = "\nCONFIG['ComparisonDistibution']['" + data_name + "']['" + data_type + "']['" + pop_name + "'] = " + data_setup
            comparison_path = CONFIG['file_path'][:-6] + 'Exp.py' if CONFIG['UseLinearGrowth'] else CONFIG['file_path'][:-6] + 'Lin.py' 
            write_to_file(msg = comparison_data_setup,
                          file_path = comparison_path,
                          msg_type = 'comparison distribution')
  
      else:
        print('# Error: Could not retrieve ',data_type,' data for', pop_name, '\nTry rerunning') 


  style_dict = dict(
     W = ['mx-', 'm--'],
     w  = ['bs-', 'b--'], 
     Overall = ['k-','k--'],
     theoretical_gompertz_R = ['mx-', 'm--'],
     theoretical_gompertz_r = ['gd-','g--'], 
  )
  graph_styles = dict()

  current_pop_prefix = 'Linear population, pure ' if CONFIG['UseLinearGrowth'] else 'Exponential population, pure '
  comparison_pop_prefix = 'Linear population, pure ' if not CONFIG['UseLinearGrowth'] else 'Exponential population, pure '
  distributions = dict()


  # Loads and formats distributions data for graphs

  for distribution_type in graph_types:
      
      distributions[distribution_type] = dict()
      use_comparison_dist = data_name in CONFIG['ComparisonDistibution'] and distribution_type in CONFIG['ComparisonDistibution'][data_name]

      # Scale Data Set
      if setAxis:
        max_array_len = 60#120
      else:
        max_array_len = 0

      for pop_name in CONFIG['display_factors']:
          dist = [] if not pop_name in current_distribution[distribution_type] else current_distribution[distribution_type][pop_name]
          graph_styles[current_pop_prefix + pop_name] = style_dict[pop_name][0] #+ 'o-' if CONFIG['UseLinearGrowth'] else  style_dict[pop_name] + 'x-' 

          if not setAxis:
            max_array_len = len(dist) if len(dist) > max_array_len else max_array_len
          distributions[distribution_type][current_pop_prefix + pop_name] = list (dist )

          if use_comparison_dist and pop_name in CONFIG['ComparisonDistibution'][data_name][distribution_type]:
            graph_styles[comparison_pop_prefix + pop_name] = style_dict[pop_name][1] #+ 'x-' if CONFIG['UseLinearGrowth'] else  style_dict[pop_name] + 'o-' 
            if not setAxis:
              array_len = len(CONFIG['ComparisonDistibution'][data_name][distribution_type][pop_name])
              max_array_len = array_len if array_len > max_array_len else max_array_len
            distributions[distribution_type][comparison_pop_prefix + pop_name] = list(CONFIG['ComparisonDistibution'][data_name][distribution_type][pop_name] )
          # else:
          #    print('Could not find ComparisonDistibution for ',distribution_type,'in', pop_name)
      for dist_key in distributions[distribution_type]:
         if len(distributions[distribution_type][dist_key]) > max_array_len:
           distributions[distribution_type][dist_key] = distributions[distribution_type][dist_key][:max_array_len]
         else:
           distributions[distribution_type][dist_key] += [0] * (max_array_len - len(distributions[distribution_type][dist_key]))

  # Create Graph Labels

  yLabels = dict(
     CDF = 'Cumulative Distribution',
     Absolute = 'Number of NMRs',
     change_per_x_unit = 'Change per x unit',
     hazard = 'Hazard (1/ticks)',
     PDF = 'Proportion of Population'
  )
  
  full_title = title + ' for ' 

  full_title += ' Using Distribution ' + CONFIG['InsultDistribution']
  full_title += ' Using log cap multiplier ' + str(CONFIG['pop_cap_mult']) + ' and exponent ' + str(CONFIG['pop_cap_exponent'])
  prefix_full_title = True if 'Absolute' in graph_types else False

  # Create Graphs
  fig_name = "graphs/"+ data_name + 'HP' + str(CONFIG['HP']) + 'hp' + str(CONFIG['hp']) + 'R' + str(CONFIG['R']) + 'r' + str(CONFIG['r'])
  fig_name += 'Lin' if CONFIG['UseLinearGrowth'] else 'Exp'
  fig_name += CONFIG['InsultDistribution']
  fig_name += 'NoCap' if not CONFIG['LogorithmicPopulationCap'] else ''

  for distribution_type in graph_types:

    use_logy=True if distribution_type == 'change_per_x_unit' or distribution_type == 'hazard' else False
    
    for pop_name in ['r', 'R']:
        alpha = .01 #The population size was too low with .1
        beta =  -log(2) / CONFIG[pop_name]    
        dt = 1 # you do an insult every dt in time, defining the parameter dt: An insult is applied every tick
        M = 48 # unknown constant I can play with

        if distribution_type == 'hazard' and (CONFIG['InsultDistribution'] == 'T8' or CONFIG['InsultDistribution'] == 'T8B'):
          # we should try plotting the theoretical Gompertz on the hazard plot. For T8, it should be hazard(t) = alpha * exp(-beta * t).
          distributions[distribution_type]['theoretical_gompertz_' + pop_name] = [alpha * exp(-beta * t) for t in range(max_array_len)]
        elif distribution_type == 'hazard' and CONFIG['InsultDistribution'] == 'T6':
          # we should try plotting the theoretical Gompertz on the hazard plot. For T6, it should be hazard(t) = exp(beta * t / 2) / (M * dt).
          distributions[distribution_type]['theoretical_gompertz_' + pop_name] = [exp(-beta * t / 2) / (M * dt) for t in range(max_array_len)]
          # d = [distributions[distribution_type]['theoretical_gompertz_' + pop_name][t]/ distributions[distribution_type][current_pop_prefix+pop_name][t]  for t in range(1,30)]
          # # print('diff = ', d)#][t]- distributions[distribution_type][pop_name][t]  for t in range(max_array_len)])
          # print(mean(d)) ## Use this to get M
        elif distribution_type == 'change_per_x_unit' and CONFIG['InsultDistribution'] == 'T6':
          # we should try plotting the theoretical Gompertz on the mortality plot. For T6, it should be D(t) = exp(beta * t / 2) / M 
          distributions[distribution_type]['theoretical_gompertz_' + pop_name] = [exp(-beta * t / 2) / (M * dt) - distributions[distribution_type][current_pop_prefix+pop_name][1] for t in range(max_array_len)]

    # if distribution_type == 'change_per_x_unit':
    #   DataFrame(distributions[distribution_type]) \
    #           .plot(title = full_title, xlabel='Age (Ticks)', ylabel=yLabels[distribution_type], style=graph_styles, logy=False)
    #   savefig(fig_name + distribution_type +'NotLog'+ '.png')
    
    if prefix_full_title:
      full_title = 'Absolute ' + full_title if distribution_type == 'Absolute' else 'Relative ' + full_title

    print(title,'graphing ',distribution_type)
    if distribution_type == 'hazard': 
      pop_key = current_pop_prefix + 'W'
      ax1 = DataFrame(distributions[distribution_type][pop_key]) \
              .plot(xlabel='Age (Ticks)', ylabel=yLabels[distribution_type], style='k', logy=use_logy, legend=False)
      if 'PDF' in graph_types and pop_key in distributions['PDF']:
        ax2 = ax1.twinx()
        DataFrame(distributions['PDF'][pop_key]).plot(ax=ax2, style='b', legend=False)
        ax2.set_ylabel(yLabels['PDF'], color='b')
        ax2.tick_params(axis='y', labelcolor='b')
        ax1.get_figure().subplots_adjust(right=0.85)
    else:
      DataFrame(distributions[distribution_type]) \
              .plot(title = full_title, xlabel='Age (Ticks)', ylabel=yLabels[distribution_type], style=graph_styles, logy=use_logy)
      figtext(.1,0,s='** This graph was created using data from ' + str(CONFIG['NumRuns'])+ ' runs')
    # figtext(.1,0.1,s='** equalibrium population size is around ' + str(CONFIG['InitialPopulation']))

    # show()
    datetime_string = datetime.now().strftime("%Y%m%d_%H%M%S")
    savefig(fig_name + distribution_type + '_' + datetime_string + '.png')

def individual_run_graph(CONFIG, run, end_tick, graph_type='population'):
    ''' Graphs stats from single run '''
    
    fig_name = "graphs/LinearRun"  if CONFIG['UseLinearGrowth'] else 'graphs/ExponentialRun'
    fig_name += CONFIG['InsultDistribution']
  
    full_title = 'Linear Growth Plot' if CONFIG['UseLinearGrowth'] else 'Exponential Growth Plot'
    style_dict = dict(
      W = 'r',
      w = 'b', 
      wW = 'm',
    )
    graph_styles = dict()

    graph_content = {}
    r_name_end = False
    hp_name_end = False
    
    # Scale Data Set
    # max_array_len = 0
    # array_len = len(current_distribution[distribution_type][pop_name])
    # max_array_len = array_len if array_len > max_array_len else max_array_len
    if CONFIG['UseColonyReproductionMode'] and graph_type == 'colony':
        graph_content = CONFIG['colony_stats']['fixation_rates'][-end_tick:]
        w_name_end = 'ColonyFixationRates.png'

    elif CONFIG['overall_stats']['w']['pop_death_ticks'][run] != CONFIG['overall_stats']['W']['pop_death_ticks'][run]:
      # Only Graphs population if someone won
      # fig_name += 'W' + str(CONFIG['W']) + 'w' + str(CONFIG['w']) 
      graph_content['Homozygous W'] = CONFIG['alive_NMRS']['WW']
      graph_content['Homozygous w'] = CONFIG['alive_NMRS']['ww']
      graph_content['Heterozygous W-w'] = CONFIG['alive_NMRS']['wW']

      graph_styles['Homozygous W'] = style_dict['W']
      graph_styles['Homozygous w'] = style_dict['w'] 
      graph_styles['Heterozygous W-w'] = style_dict['wW']
      
      w_name_end = 'WWin.png' if CONFIG['overall_stats']['w']['pop_death_ticks'][run] > 0 else 'wWin.png'

    # print(graph_content)

    if len(graph_content) > 0:
      DataFrame(graph_content).plot(title = full_title, xlabel = 'Time (ticks)', ylabel = 'Population', style=graph_styles)
      figtext(.25,0,s='** This graph was created using ' + str(CONFIG['NumRuns'])+ ' runs')
      show()
      
      if w_name_end:
        savefig(fig_name + w_name_end)

      close()
      
def write_log_comparison_to_csv(CONFIG, age_distribution):

    ''' For use in comparing log cap multiplier and exponent, 
        NOTE: import may mess with CONFIG values '''

    for pop_name in CONFIG['display_factors']:
      # Adds Column labels to file if they dont exist yet
      # max_age = "'-'" if not pop_name in age_distribution['CDF'] else str(len(age_distribution['CDF'][pop_name]))
      max_age_str = "CONFIG['max_experimental_age']['"+pop_name+"'] = "+ str(CONFIG['max_experimental_age'][pop_name])
      write_to_file(msg='\n'+max_age_str,
                            file_path=CONFIG['file_path'],
                            msg_type='maximimum age')
      # print("CONFIG['max_experimental_age']['R']1",CONFIG['max_experimental_age']['R'])

    if 'Exp.py' in CONFIG['file_path'] and 'Overall' in CONFIG['display_factors']:
      # Only write comparisons if Lin has already run
      comparison_path = CONFIG['file_path'][:-6] + 'Lin' 
      comparison_path = 'config_files.' + comparison_path[13:]
      
      try:
        avg_ls = float(CONFIG['avg_lifespan']['Overall']) if 'Overall' in CONFIG['avg_lifespan'] else '-'
        comparison_CONFIG = import_module(comparison_path)
      except (ImportError, ModuleNotFoundError) as e:
        print(f'Warning: Could not import comparison config from {comparison_path}: {e}')
        return
      except KeyError as e:
        print(f'Warning: Missing key in CONFIG: {e}')
        return
      except Exception as e:
        print(f'Unexpected error loading comparison config: {e}')
        return

      csv_name = 'csv/'+ str(CONFIG['InsultDistribution'])+'.csv'
      header = 'pop_cap_exponent, pop_cap_mult, HP, hp, R, r, ticks until fertile, population size lin, population size exp, maximum age lin, maximum age exp, percent run failures lin, percent run failures exp, avg lifespan lin, avg lifespan exp'
      
      try:
        if not isfile(csv_name):
          with open(csv_name, 'w') as w:
            w.write(header)
        
        line = '\n' +str(CONFIG['pop_cap_exponent']) + ',' + str(CONFIG['pop_cap_mult']) + ',' + str(CONFIG['TicksTillFertile']) + ','\
          + str(CONFIG['HP']) + ',' + str(CONFIG['hp']) + ',' + str(CONFIG['R']) + ',' + str(CONFIG['r']) + ',' \
          + str(comparison_CONFIG.CONFIG['InitialPopulation'])  + ',' + str(CONFIG['InitialPopulation']) + ',' \
          + str(comparison_CONFIG.CONFIG['max_experimental_age']['Overall']) + ',' + str(CONFIG['max_experimental_age'][pop_name]) + ',' \
          + str(comparison_CONFIG.CONFIG['pcnt_run_failures'])  + ',' + str(CONFIG['pcnt_run_failures']) + ',' \
          + str(comparison_CONFIG.CONFIG['avg_lifespan']['Overall'])  + ',' + str(avg_ls)

        with open(csv_name, 'a') as w:
          w.write(line)
      except (IOError, OSError) as e:
        print(f'Failed to write CSV comparison data to {csv_name}: {e}')
      except KeyError as e:
        print(f'Missing required CONFIG key for CSV comparison: {e}')
      except Exception as e:
        print(f'Unexpected error writing CSV comparison: {e}')
              
              
def print_allele_counts(CONFIG, allele_name):
    recessive_allele_name = allele_name.lower() + allele_name.lower()
    dominant_allele_name = allele_name.upper() + allele_name.upper()
    heterozygous_allele_name = allele_name.lower() + allele_name.upper()
    
    if CONFIG['alive_NMRS'][recessive_allele_name][-1] == 0:
      print('The recessive allele(',allele_name.lower(),') is extinct')
    elif CONFIG['alive_NMRS'][dominant_allele_name][-1] == 0:
      print('The dominant allele(',allele_name.upper(),') is extinct')  
    else:
      print('Time limit exceeded')
      
    totalwAlleles = [(2*CONFIG['alive_NMRS'][recessive_allele_name][i]) + CONFIG['alive_NMRS'][heterozygous_allele_name][i] for i in range(len(CONFIG['alive_NMRS'][heterozygous_allele_name]))]
    totalWAlleles = [(2*CONFIG['alive_NMRS'][dominant_allele_name][i]) + CONFIG['alive_NMRS'][heterozygous_allele_name][i] for i in range(len(CONFIG['alive_NMRS'][heterozygous_allele_name]))]
    print('Recessive alleles(', allele_name.lower(), '): ', totalwAlleles)
    print('Dominant alleles (', allele_name.upper(), '): ', totalWAlleles)


def get_scenario_data_setup(CONFIG, run):
  ''' Generates the data setup for the scenario '''


  """
  Scenario 1:
  Linear Individual Mode with R9r10HP1000hp900pcE1.1pcM0.07BP0.3TTF3IRT2
  Scenario 2:
  Linear Individual Mode with R9r9.2HP1000hp900pcE1.1pcM0.07BP0.3TTF3IRT2
  Scenario 3:
  Linear Individual Mode with R9r9.3HP1000hp600pcE1.1pcM0.07BP0.3TTF3IRT2
  Scenario 4:
  Exponential Individual Mode with R9r10HP1000hp900pcE1.1pcM0.07BP0.3TTF3IRT2
  Scenario 5:
  Exponential Individual Mode with R9r9.2HP1000hp900pcE1.1pcM0.07BP0.3TTF3IRT2
  Scenario 6:
  Exponential Individual Mode with R9r9.3HP1000hp600pcE1.1pcM0.07BP0.3TTF3IRT2
  Scenario 7:
  Linear Ant Mode with R9r10HP1000hp900pcE1.1pcM0.07BP0.3TTF3IRT2
  Scenario 8:
  Linear Ant Mode with R9r9.2HP1000hp900pcE1.1pcM0.07BP0.3TTF3IRT2
  Scenario 9:
  Linear Ant Mode with R9r9.3HP1000hp600pcE1.1pcM0.07BP0.3TTF3IRT2
  Scenario 10:
  Linear Mole Rat Mode with R9r10HP1000hp900pcE1.1pcM0.07BP0.3TTF3IRT2
  Scenario 11:
  Linear Mole Rat Mode with R9r9.2HP1000hp900pcE1.1pcM0.07BP0.3TTF3IRT2
  Scenario 12:
  Linear Mole Rat Mode with R9r9.3HP1000hp600pcE1.1pcM0.07BP0.3TTF3IRT2
  Scenario -999:
  Everything else
  """

  scenario = -999
  mode = get_run_mode(CONFIG)

  if mode == 'Individual':
    if CONFIG['file_name'] == 'R9r10HP1000hp900pcE1.1pcM0.07BP0.3TTF3IRT2Lin':
      scenario = 1
    elif CONFIG['file_name'] == 'R9r9.2HP1000hp900pcE1.1pcM0.07BP0.3TTF3IRT2Lin':
      scenario = 2
    elif CONFIG['file_name'] == 'R9r9.3HP1000hp600pcE1.1pcM0.07BP0.3TTF3IRT2Lin':
      scenario = 3
    elif CONFIG['file_name'] == 'R9r10HP1000hp900pcE1.1pcM0.07BP0.3TTF3IRT2Exp':
      scenario = 4
    elif CONFIG['file_name'] == 'R9r9.2HP1000hp900pcE1.1pcM0.07BP0.3TTF3IRT2Exp':
      scenario = 5
    elif CONFIG['file_name'] == 'R9r9.3HP1000hp600pcE1.1pcM0.07BP0.3TTF3IRT2Exp':
      scenario = 6
  if mode == 'Ant':
    if CONFIG['file_name'] == 'R9r10HP1000hp900pcE1.1pcM0.07BP0.3TTF3IRT2Lin':
      scenario = 7
    elif CONFIG['file_name'] == 'R9r9.2HP1000hp900pcE1.1pcM0.07BP0.3TTF3IRT2Lin':
      scenario = 8
    elif CONFIG['file_name'] == 'R9r9.3HP1000hp600pcE1.1pcM0.07BP0.3TTF3IRT2Lin':
      scenario = 9
  if mode == 'Mole Rat':
    if CONFIG['file_name'] == 'R9r10HP1000hp900pcE1.1pcM0.07BP0.3TTF3IRT2Lin':
      scenario = 10
    elif CONFIG['file_name'] == 'R9r9.2HP1000hp900pcE1.1pcM0.07BP0.3TTF3IRT2Lin':
      scenario = 11
    elif CONFIG['file_name'] == 'R9r9.3HP1000hp600pcE1.1pcM0.07BP0.3TTF3IRT2Lin':
      scenario = 12


  
    # Called at the end of each run
  full_scenario_header = "Scenario, num_runs, date_time, mode"
  ''' "Scenario,num_runs, init_rec_allele_fraction, init_pop, \
LitterSize, BirthProbability, reproduction_type, insult_recovery_ticks, \
pop_cap_exponent, pop_cap_mult, ticks until fertile, HP, hp, R, r, max_ticks, \
use_immortal_queen_mode, use_colony_reproduction_mode, n0_colonies, \
threshold_colony_size, mating_flight_selection, queen_hp_multiplier, \
base_colony_formation_probability, colony_logistic_mult, colony_logistic_exp, \
pioneer_group_size, colony_disaster_probability, colony_disaster_insult_multiplier, \
ticks_to_queen_replacement, mates_within_colony, date_time, mode"

'''
  time_stamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
  reproduction_type = 'Linear' if CONFIG['UseLinearGrowth'] else 'Exponential'
  full_scenario_data_setup = str(scenario) + ',' + str(run) + ',' + str(time_stamp) + ',' + mode
  
  '''
    full_scenario_data_setup = str(scenario) + ',' + str(CONFIG['NumRuns']) + ',' + str(CONFIG['InitialRecessiveAlleleFraction']) + ',' + str(CONFIG['InitialPopulation']) + ',' \
  + str(CONFIG['LitterSize']) + ',' + str(CONFIG['BirthProbability']) + ',' + reproduction_type + ',' + str(CONFIG['insult_recovery_ticks']) + ',' \
  + str(CONFIG['pop_cap_exponent']) + ',' + str(CONFIG['pop_cap_mult']) + ',' + str(CONFIG['TicksTillFertile']) + ',' + str(CONFIG['HP']) + ',' + str(CONFIG['hp']) + ',' \
  + str(CONFIG['R']) + ',' + str(CONFIG['r']) + ',' + str(CONFIG['MaxTicks']) + ',' + str(CONFIG['UseImmortalQueenMode']) + ',' + str(CONFIG['UseColonyReproductionMode']) + ',' \
  + str(CONFIG['N0Colonies']) + ',' + str(CONFIG['ThresholdColonySize']) + ',' + str(CONFIG['MatingFlightSelection']) + ',' + str(CONFIG['QueenHPMultiplier']) + ',' \
  + str(CONFIG['BaseColonyFormationProbability']) + ',' + str(CONFIG['ColonyLogisticMult']) + ',' + str(CONFIG['ColonyLogisticExp']) + ',' + str(CONFIG['PioneerGroupSize']) + ',' \
  + str(CONFIG['ColonyDisasterProbability']) + ',' + str(CONFIG['ColonyDisasterInsultMultiplier']) + ',' + str(CONFIG['TicksToQueenReplacement']) + ',' \
  + str(CONFIG['MatesWithinColony'])+',' + str(time_stamp) + ',' + mode'''
  return full_scenario_header, full_scenario_data_setup, scenario, run, mode
  
def suppfig_3F_and_4D_csv(CONFIG, winner, run):
  ''' Generates CSV files for the supplementary figure '''
  """
Supplementary Figure (parameter sweep of colony reproduction)

This will be a bar plot showing the frequency with which each allele wins in each mode under each parameter scenario for a relevant parameter. We discussed that this would be the threshold colony size for reproduction. We can also sweep something like the number of starting health points and half-life of the longevity and vitality alleles. 

This is the only data set (data_suppfig.csv) where I need data for more than one parameter scenario.

The headings would be the following:
Scenario [integer]: which parameter scenario the data point refers to. This will encompass parameters like colony threshold size and the Gompertz parameters
Run: as above
Mode: as above
Winner: as above
  """
  mode = get_run_mode(CONFIG)


  supplementary_figure_header = "Scenario, Run, Mode, Winner"

  full_scenario_header, full_scenario_data_setup, scenario, run_num, mode = get_scenario_data_setup(CONFIG, run)
  supplementary_figure_data_setup = str(scenario) + ',' + str(run) + ',' + str(mode) + ',' + winner + '\n'

  write_to_file(msg = supplementary_figure_data_setup,
              file_path = "csv/data_suppfig.csv",
              msg_type = 'Failed to write supplementary figure data to file on run {Run}',
              header = supplementary_figure_header)

  if mode == 'Individual':
    figure3F_csv(msg = supplementary_figure_data_setup,
              header = supplementary_figure_header)
  else:
    figure4D_csv(msg = supplementary_figure_data_setup,
              header = supplementary_figure_header)

def write_to_csv_each_tick(CONFIG, tick, run, colonies):
  ''' Writes to CSV files each tick '''

  full_scenario_header, full_scenario_data_setup, scenario, run_num, mode = get_scenario_data_setup(CONFIG, run)
  full_message = ""
  header = ""
  file_path = ""

  if CONFIG['UseColonyReproductionMode']:
    header = "Scenario, Run, Time, Mode, Allele, Age, Sex, Dominant Allele"
    file_path = "csv/data_fig4C.csv"
    longevity_allele_name = get_longevity_allele(CONFIG)

    # data_fig4C.csv 
    for colony in colonies:
      fixed_allele_name = colony.get_fixed_allele_name()
      alleles = ""
      if fixed_allele_name == longevity_allele_name:
        alleles = "Longevity"
      elif fixed_allele_name == None:
        alleles = "Mix"
      else:
        alleles = "Vitality"
      
      if CONFIG['UseLinearGrowth'] and colony.queen is not None:
        full_message += str(scenario) + ',' + str(run_num) + ',' + str(tick) + ',' + str(mode) + ',' + alleles + ',' + str(tick - colony.queen.birth) + ',' + 'queen'+ ',' + get_longevity_or_vitality_allele(CONFIG, colony.queen.w_pops['type_name']) + '\n'

      for member in colony.members:
        sex = 'male' if member.sex == 'M' else 'female'
        full_message += str(scenario) + ',' + str(run_num) + ',' + str(tick) + ',' + str(mode) + ',' + alleles + ',' + str(tick - member.birth) + ',' + sex + ',' + get_longevity_or_vitality_allele(CONFIG, member.w_pops['type_name']) + '\n'


  else:
    data_3AtoE_header = "Scenario, Run, Time, Population, Age, Sex, Allele1, Allele2, DominantAllele"
    file_path = "csv/data_fig3AtoE.csv"
    header = data_3AtoE_header
    population = 'Linear' if CONFIG['UseLinearGrowth'] else 'Exponential'
    queen = CONFIG['queen'] if CONFIG['UseLinearGrowth'] else None
    for nmr in CONFIG['alive_NMRS']['female'] + CONFIG['alive_NMRS']['male']:
      sex = nmr.sex if nmr is not queen else 'queen'
      full_message += str(scenario) + ',' + str(run) + ',' + str(tick) + ',' + population + ',' + str(tick - nmr.birth) + ',' + sex + ',' + get_longevity_or_vitality_allele(CONFIG, nmr.w_pops['alleles'][0]) + ',' + get_longevity_or_vitality_allele(CONFIG, nmr.w_pops['alleles'][1]) + ',' + get_longevity_or_vitality_allele(CONFIG, nmr.w_pops['type_name']) + '\n'  
  
  write_to_file(msg = full_message,
              file_path = file_path,
              msg_type = f'Failed to write {file_path} data to file',
              header = header)

def figure4D_csv(msg, header):
  """
  In the data frame for panel D (data_fig4D.csv), each row pertains to a simulation run of a given scenario and a given mode. The headings would be
  Scenario[integer]: which Gompertz scenario was used. 
  Run: as above
  Mode: as above
  Winner [character]: which allele fixed in all populations by the end of the simulation ("Longevity", "Vitality", or "Neither")

  """
  write_to_file(msg = msg,
              file_path = "csv/data_fig4D.csv",
              msg_type = 'Failed to write figure 4D data to file',
              header = header)


def figure3F_csv(msg, header):
  """
  The second data frame (data_fig3F.csv) would give me a summary of the outcome for all runs for Fig 3F. That data frame would have the following headings:
  Scenario [integer]: which parameter scenario the data refers to 
  Run [integer]: which simulation replicate the data refers to
  Population [character]: "Linear" or "Exponential"
  Ticks [integer]: number of ticks the simulation was run for
  Winner [character]: Identity of the allele that fixed in the population. Either "Longevity", "Vitality", or "died out" or "timed out" (in case neither allele has fixed by the end of the simulation)
  """
  
  write_to_file(msg = msg,
              file_path = "csv/data_fig3F.csv",
              msg_type = 'Failed to write figure 3F data to file',
              header = header)

def figure4AtoB_csvs(CONFIG, tick, num_W_colonies, num_w_colonies, num_mix_colonies, run):
    ''' Generates CSV files for the figure '''

    """
    In the data frame for panels A and B (data_fig4AB.csv), each row pertains to a time point for a given run of a given scenario and a given mode. The headings would be 
    Scenario: as above
    Run: as above
    Time: as above
    Mode [character]: "Ant", "Mole Rat", or "Individual"
    Longevity [character]: how many colonies carry only the longevity allele
    Vitality [character]: how many colonies carry only the vitality allele
    Mix [character]: how many colonies carry both alleles
    """

    data_4AtoB_header = "Scenario, Run, Time, Mode, Longevity, Vitality, Mix"
    longevity_allele_name = get_longevity_allele(CONFIG)
    vitality = num_w_colonies if longevity_allele_name == 'W' else num_W_colonies
    longevity = num_W_colonies if longevity_allele_name == 'W' else num_w_colonies
 

    full_scenario_header, full_scenario_data_setup, scenario, run_num, mode = get_scenario_data_setup(CONFIG, run)
    data_4AtoB_msg = str(scenario) + ',' + str(run_num) + ',' + str(tick) + ',' + str(mode) + ',' + str(longevity) + ',' + str(vitality) + ',' + str(num_mix_colonies) + '\n'

    write_to_file(msg = data_4AtoB_msg,
              file_path = "csv/data_fig4AB.csv",
              msg_type = 'Failed to write figure 4AtoB data to file',
              header =  data_4AtoB_header)

