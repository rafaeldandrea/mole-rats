
from genericpath import isfile
import importlib
import math
from statistics import mean, median
from matplotlib.pyplot import show, savefig, close, figtext
from numpy import trim_zeros
from pandas import DataFrame

def run_stats_using_individual_data(CONFIG):
    ''' Shows overall run stats only collected using individual run data '''
    
    if CONFIG['GetEquilibriumPopulation']:
      data_setup = "\n\nCONFIG['InitialPopulation'] = " 
      if CONFIG['UseOverallRunStats'] and len(CONFIG['overall_stats']['Overall']['avg_equilibrium_population_size']) > 0:
          data_setup += str(round(mean(CONFIG['overall_stats']['Overall']['avg_equilibrium_population_size']))) 
      # elif not CONFIG['UseOverallRunStats']:
          CONFIG['InitialPopulation'] = str(round(mean(CONFIG['overall_stats']['Overall']['avg_equilibrium_population_size']))) 
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
        CONFIG['InitialRecessiveAlleleFraction'] * 100, ' percent of the population\n',
        rec_info
    )
    print (
        dom_name, ' was a dominant allele, that began as ', 
        abs(1 - CONFIG['InitialRecessiveAlleleFraction']) * 100, ' percent of the population\n',
        dom_info
    )

    '''In  100  runs, population  r  died out  63  times. That's  63.0  percent of the time.

In  100  runs, population  R  died out  36  times. That's  36.0  percent of the time.

 

I would like you to reverse the statement, so it would say:

 

In  100  runs, population  r  fixed 36  times. That's  36.0  percent of the time.

In  100  runs, population  R  fixed  63  times. That's  63.0  percent of the time.

'''

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

    CONFIG['pcnt_run_failures'] = abs(CONFIG['NumRuns'] - CONFIG['overall_stats'][rec_name]['pop_death_count'] - CONFIG['overall_stats'][dom_name]['pop_death_count']) / CONFIG['NumRuns']
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
        # print('adding',  mean(equilibrium_population_size[pop_name]), ' to ', pop_name)

    ages_at_tick_2000 = [0] * 500
    death_insults_for_old_nmrs = []
    if not 'old_insults' in CONFIG:
      CONFIG['old_insults'] = []    
    # if not 'old_insults2' in CONFIG:
    #   CONFIG['old_insults2'] = []
    # instultsat1000 = []
    # print(len(population))
    # print("CONFIG['MaxAge']",CONFIG['MaxAge'])
    # print("CONFIG['MaxTicks'] ",CONFIG['MaxTicks'] )
    for nmr in population:
        if pop_name == 'Overall' or nmr.hp_pops['type_name'] == pop_name or nmr.r_pops['type_name'] == pop_name:
          if nmr.death == None:
            nmr.age_NMR(CONFIG, run, end_tick)
          if nmr.death > 0:
            total_insult_on_death = sum(nmr.current_insults)
            ages.append(abs(abs(nmr.death - nmr.birth-1)))
            insults.append(total_insult_on_death)
            deaths[nmr.death] += 1
            CONFIG['max_experimental_age'][pop_name] = CONFIG['max_experimental_age'][pop_name] if nmr.death-nmr.birth < CONFIG['max_experimental_age'][pop_name] else nmr.death-nmr.birth
        #     if pop_name == 'hp' and nmr.death - nmr.birth > 98:
        #       death_insults_for_old_nmrs.append(total_insult_on_death)#
        #       CONFIG['old_insults'].append(nmr.current_insults[-1])
        #       # CONFIG['old_insults'].append(nmr.current_insults[-1])
        #     elif pop_name == 'HP' and nmr.death - nmr.birth > 89:
        #       death_insults_for_old_nmrs.append(total_insult_on_death)
        #       CONFIG['old_insults'].append(nmr.current_insults[-1])
        if nmr.death and nmr.death >= 2000 and nmr.birth <= 2000:
          nmr_age  =  2000 - nmr.birth 
          ages_at_tick_2000[nmr_age] += 1 
          # if len(nmr.current_insults) > 0:
          #   instultsat1000.append(sum(nmr.current_insults))


    if len(ages) > 0:
      CONFIG['overall_stats'][pop_name]['avg_ls'].append(mean(ages))
      if CONFIG['GetLifespanDistribution']:
        if CONFIG['overall_stats'][pop_name]['all_ls'] is not list:
          CONFIG['overall_stats'][pop_name]['all_ls'] = list(CONFIG['overall_stats'][pop_name]['all_ls']) + ages if len(CONFIG['overall_stats'][pop_name]['all_ls']) > 1 else ages
        else:
          CONFIG['overall_stats'][pop_name]['all_ls'].extend(ages)

      # Generating ages file
      if pop_name == 'R' or pop_name == 'r':
        ages_at_tick_2000 = trim_zeros(ages_at_tick_2000, 'b')
        for x in range(len(ages_at_tick_2000)):
          header = 'age, count, allele, run, tick number(2000), type\n'
          # print('writing ages to file ages')
          # write_to_file('ages', 'ages.txt') # population, exp or lin, 
          # write_to_file(ages, 'ages.txt')
          run_type = 'Linear' if CONFIG['UseLinearGrowth'] else 'Exponential'
          line = str(x)+ ',' + str(ages_at_tick_2000[x])+ ',' + pop_name+ ',' + str(run) + ',2000,' + run_type + '\n'
          line = header + line if run == 0 and x==0 else line
          write_to_file(line, 'data/ages-data.csv', header)
          # write_to_file(line, 'data/ages-' + CONFIG['file_name'] + '.csv', header)


        # print('avg killing insult', mean(insults))
        # print('avg killing insult 1000', mean(instultsat1000))
        # print('avg old killing insult', mean(death_insults_for_old_nmrs))
        # print('pop size at 1000', sum(ages_at_tick_2000))
        # print('avg insult',mean(CONFIG['insults']))
        # print()
          
        

    # if CONFIG['ShowIndividualRunStats']:
      
    #   deaths = trim_zeros(deaths, 'b')
    #   # print('The Death Array for Population', pop_name, 'is', deaths)

    #   if len(ages) > 0:
    #     print ('Average Population ', pop_name,' Lifespan: ',
    #             mean(ages))
    #     print ('Median Population ', pop_name,' Lifespan: ',
    #             median(ages))
    #   if len(insults) > 0:
    #     print ('Average Population ', pop_name,' Killing Insult: ',
    #             mean(insults))
    #     print ('Median Population ', pop_name,' Killing Insult: ',
    #             median(insults), '\n')

def write_to_file(msg, file_path, msg_type='message', header=''):
  try: 
    with open(file_path, "a") as myfile:
      myfile.write(msg)
  except:
      try: 
        with open(file_path, "w") as myfile:
          myfile.write(header)
          myfile.write(msg)
      except:
        print('Failed to write ' + msg_type + ' to file ' + file_path)

      

  
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
          if 'Lin.py' in CONFIG['file_path'] or 'Exp.py' in CONFIG['file_path']:
            comparison_data_setup = "\nCONFIG['ComparisonDistibution']['" + data_name + "']['" + data_type + "']['" + pop_name + "'] = " + data_setup
            comparison_path = CONFIG['file_path'][:-6] + 'Exp.py' if CONFIG['UseLinearGrowth'] else CONFIG['file_path'][:-6] + 'Lin.py' 
            write_to_file(msg = comparison_data_setup,
                          file_path = comparison_path,
                          msg_type = 'comparison distribution')
  
      else:
        print('# Error: Could not retrieve ',data_type,' data for', pop_name, '\nTry rerunning') 


  style_dict = dict(
     R = ['ro-', 'r--'],
     HP = ['mx-', 'm--'],
     r = ['bs-', 'b--'], 
     hp = ['gd-','g--'],     
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
     hazard = 'hazard',
     PDF = 'Proportion of Population'
  )
  
  full_title = title + ' for ' 

  # if CONFIG['show_r'] and CONFIG['show_hp']:
  #   full_title += 'Gompertz Allele (R or r) and Health Points Allele (HP or hp)' 
  # else:
  #   full_title += 'Health Points allele (HP or hp)' if CONFIG['show_hp'] else 'Gompertz Allele (R or r)' 

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
        beta =  -math.log(2) / CONFIG[pop_name]    
        dt = 1 # you do an insult every dt in time, defining the parameter dt: An insult is applied every tick
        M = 48 # unknown constant I can play with
        if distribution_type == 'hazard' and (CONFIG['InsultDistribution'] == 'T8' or CONFIG['InsultDistribution'] == 'T8B'):
          # we should try plotting the theoretical Gompertz on the hazard plot. For T8, it should be hazard(t) = alpha * math.exp(-beta * t).
          distributions[distribution_type]['theoretical_gompertz_' + pop_name] = [alpha * math.exp(-beta * t) for t in range(max_array_len)]
        elif distribution_type == 'hazard' and CONFIG['InsultDistribution'] == 'T6':
          # we should try plotting the theoretical Gompertz on the hazard plot. For T6, it should be hazard(t) = exp(beta * t / 2) / (M * dt).
          distributions[distribution_type]['theoretical_gompertz_' + pop_name] = [math.exp(-beta * t / 2) / (M * dt) for t in range(max_array_len)]
          # d = [distributions[distribution_type]['theoretical_gompertz_' + pop_name][t]/ distributions[distribution_type][current_pop_prefix+pop_name][t]  for t in range(1,30)]
          # # print('diff = ', d)#][t]- distributions[distribution_type][pop_name][t]  for t in range(max_array_len)])
          # print(mean(d)) ## Use this to get M
        elif distribution_type == 'change_per_x_unit' and CONFIG['InsultDistribution'] == 'T6':
          # we should try plotting the theoretical Gompertz on the mortality plot. For T6, it should be D(t) = exp(beta * t / 2) / M 
          distributions[distribution_type]['theoretical_gompertz_' + pop_name] = [math.exp(-beta * t / 2) / (M * dt) - distributions[distribution_type][current_pop_prefix+pop_name][1] for t in range(max_array_len)]

    # if distribution_type == 'change_per_x_unit':
    #   DataFrame(distributions[distribution_type]) \
    #           .plot(title = full_title, xlabel='Age (Ticks)', ylabel=yLabels[distribution_type], style=graph_styles, logy=False)
    #   savefig(fig_name + distribution_type +'NotLog'+ '.png')
    
    if prefix_full_title:
      full_title = 'Absolute ' + full_title if distribution_type == 'Absolute' else 'Relative ' + full_title

    print(title,'graphing ',distribution_type)
    DataFrame(distributions[distribution_type]) \
            .plot(title = full_title, xlabel='Age (Ticks)', ylabel=yLabels[distribution_type], style=graph_styles, logy=use_logy)
    figtext(.1,0,s='** This graph was created using data from ' + str(CONFIG['NumRuns'])+ ' runs')
    # figtext(.1,0.1,s='** equalibrium population size is around ' + str(CONFIG['InitialPopulation']))

    # show()
    savefig(fig_name + distribution_type + '.png')

def individual_run_graph(CONFIG, run):
    ''' Graphs stats from single run '''
    
    fig_name = "graphs/LinearRun"  if CONFIG['UseLinearGrowth'] else 'graphs/ExponentialRun'
    fig_name += CONFIG['InsultDistribution']
  
    full_title = 'Linear Growth Plot' if CONFIG['UseLinearGrowth'] else 'Exponential Growth Plot'
    style_dict = dict(
      R = 'r',
      HP = 'm',
      r = 'b', 
      hp = 'g',     
      hpHP = 'k',
      rR = 'y'
    )
    graph_styles = dict()

    graph_content = {}
    r_name_end = False
    hp_name_end = False
    
    # Scale Data Set
    # max_array_len = 0
    # array_len = len(current_distribution[distribution_type][pop_name])
    # max_array_len = array_len if array_len > max_array_len else max_array_len
        
    if CONFIG['show_hp'] and CONFIG['overall_stats']['hp']['pop_death_ticks'][run] != CONFIG['overall_stats']['HP']['pop_death_ticks'][run]:
      # Only Graphs population if someone won
      fig_name += 'HP' + str(CONFIG['HP']) + 'hp' + str(CONFIG['hp']) 
      graph_content['Homozygous HP'] = CONFIG['alive_NMRS']['HPHP']
      graph_content['Homozygous hp'] = CONFIG['alive_NMRS']['hphp']
      graph_content['Heterozygous HP-hp'] = CONFIG['alive_NMRS']['hpHP']

      graph_styles['Homozygous HP'] = style_dict['HP']
      graph_styles['Homozygous hp'] = style_dict['hp'] 
      graph_styles['Heterozygous HP-hp'] = style_dict['hpHP']
      
      hp_name_end = 'HPWin.png' if CONFIG['overall_stats']['hp']['pop_death_ticks'][run] > 0 else 'hpWin.png'

    if CONFIG['show_r'] and CONFIG['overall_stats']['r']['pop_death_ticks'][run] != CONFIG['overall_stats']['R']['pop_death_ticks'][run]:

      fig_name += 'R' + str(CONFIG['R']) + 'r' + str(CONFIG['r'])
      graph_content['Homozygous R'] = CONFIG['alive_NMRS']['RR']
      graph_content['Homozygous r'] = CONFIG['alive_NMRS']['rr']
      graph_content['Heterozygous R-r'] = CONFIG['alive_NMRS']['rR']

      graph_styles['Homozygous R'] = style_dict['R']
      graph_styles['Homozygous r'] = style_dict['r'] 
      graph_styles['Heterozygous R-r'] = style_dict['rR']
      
      r_name_end = 'rWin.png' if CONFIG['overall_stats']['r']['pop_death_ticks'][run] == 0 else 'RWin.png'

    # print(graph_content)

    if len(graph_content) > 0:
      DataFrame(graph_content).plot(title = full_title, xlabel = 'Time (ticks)', ylabel = 'Population', style=graph_styles)
      # figtext(.25,0,s='** This graph was created using ' + str(CONFIG['NumRuns'])+ ' runs')
      show()
      
      if hp_name_end:
        savefig(fig_name + hp_name_end)
      if r_name_end:
        savefig(fig_name + r_name_end)

      close()
      
def write_log_comparison_to_csv(CONFIG, age_distribution):

    ''' For use in comparing log cap multiplier and exponent, 
        NOTE: import may mess with CONFIG values '''

    for pop_name in CONFIG['display_factors']:
      # Adds Column labels to file if they dont exist yet
      # max_age = "'-'" if not pop_name in age_distribution['CDF'] else str(len(age_distribution['CDF'][pop_name]))
      max_age_str = "CONFIG['max_experimental_age']['"+pop_name+"'] = "+ CONFIG['max_experimental_age'][pop_name]
      write_to_file(msg='\n'+max_age_str,
                            file_path=CONFIG['file_path'],
                            msg_type='maximimum age')
      # print("CONFIG['max_experimental_age']['R']1",CONFIG['max_experimental_age']['R'])

    if 'Exp.py' in CONFIG['file_path'] and 'Overall' in CONFIG['display_factors']:
      # Only write comparisons if Lin has already run
      comparison_path = CONFIG['file_path'][:-6] + 'Lin' 
      comparison_path = 'config_files.' + comparison_path[13:]
      # print(len('config_files'))
      # comparison_path[12] = '.'
      # print('comparison_path',comparison_path)
      # c_CONFIG =  dict(CONFIG)
      
      avg_ls = float(CONFIG['avg_lifespan']['Overall']) if 'Overall' in CONFIG['avg_lifespan'] else '-'
      comparison_CONFIG = importlib.import_module(comparison_path)   
                

      csv_name = 'csv/'+ str(CONFIG['InsultDistribution'])+'.csv'
      if not isfile(csv_name):
        w = open(csv_name,'w')
        header = 'pop_cap_exponent, pop_cap_mult, ticks until fertile, HP, hp, R, r, population size lin, population size exp, maximum age lin, maximum age exp, percent run failures lin, percent run failures exp, avg lifespan lin, avg lifespan exp'
        # header = 'pop_cap_exponent, pop_cap_mult, population size lin, population size exp, maximum age lin R, maximum age exp R, maximum age lin r, maximum age exp r'
        w.write(header)
        w.close()
      # print("CONFIG['max_experimental_age']['R']",CONFIG['max_experimental_age']['R'])
      # print("comparison_CONFIG.CONFIG['max_experimental_age']['R']",comparison_CONFIG.CONFIG['max_experimental_age']['R'])
        header = 'pop_cap_exponent, pop_cap_mult, HP, hp, R, r, ticks until fertile, population size lin, population size exp, maximum age lin, maximum age exp, percent run failures lin, percent run failures exp, avg lifespan lin, avg lifespan exp'

      w = open(csv_name,'a')
      line = '\n' +str(CONFIG['pop_cap_exponent']) + ',' + str(CONFIG['pop_cap_mult']) + ',' + str(CONFIG['TicksTillFertile']) + ','\
        + str(CONFIG['HP']) + ',' + str(CONFIG['hp']) + ',' + str(CONFIG['R']) + ',' + str(CONFIG['r']) + ',' \
        + str(comparison_CONFIG.CONFIG['InitialPopulation'])  + ',' + str(CONFIG['InitialPopulation']) + ',' \
        + str(comparison_CONFIG.CONFIG['max_experimental_age']['Overall']) + ',' + str(CONFIG['max_experimental_age'][pop_name]) + ',' \
        + str(comparison_CONFIG.CONFIG['pcnt_run_failures'])  + ',' + str(CONFIG['pcnt_run_failures']) + ',' \
        + str(comparison_CONFIG.CONFIG['avg_lifespan']['Overall'])  + ',' + str(avg_ls) 
        #, percent run failures lin, percent run failures exp, 
        # avg lifespan lin, avg lifespan exp

      # print('max_age',max_age)
      w.write(line)
      w.close()
              
      # print(max_age_str)

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
