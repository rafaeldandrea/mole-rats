# Standard library imports
from random import choice, choices, random, shuffle
from math import exp
from time import asctime, time

# Local imports
from Naked_Mole_Rat import Naked_Mole_Rat
from colony import Colony


### FUNCTIONS ###
def _create_nmrs_at_equilibrium_template():
    '''Creates the NMRs_at_equilibrium_at dictionary structure template'''
    percent_keys = ['50', '75', '85', '100', '100plus20']
    genders = ['female', 'male', 'Overall']
    return {
        gender: {f'percent_{pct}': [] for pct in percent_keys}
        for gender in genders
    }

def _create_population_stats_template(num_runs, include_equilibrium_tracking=False):
    '''Creates a template for population statistics dictionary
    
    Args:
        num_runs: Number of simulation runs
        include_equilibrium_tracking: Whether to include equilibrium tracking fields
    
    Returns:
        Dictionary template for population statistics
    '''
    template = {
        'pop_death_count': 0,
        'avg_ls': [],
        'pop_death_ticks': [0] * num_runs,
        'avg_equilibrium_population_size': [],
    }
    
    if include_equilibrium_tracking:
        template.update({
            'ticks_to_50pct_equilibrium': [],
            'ticks_to_75pct_equilibrium': [],
            'ticks_to_85pct_equilibrium': [],
            'ticks_to_100pct_equilibrium': [],
            'ticks_to_100plus20pct_equilibrium': [],
            'NMRs_at_equilibrium_at': _create_nmrs_at_equilibrium_template(),
        })
    
    return template

def reset_overall_stats(CONFIG):
  '''Initializes overall_stats dict so that run data is clean'''
  
  CONFIG['overall_stats'] = dict(
      start_time = time(),
      pop_death_ticks = [0] * CONFIG['NumRuns'], # Number of ticks it took for a population to die each run
  )

  print('Starting Run at ', asctime())

  CONFIG['overall_stats']['W'] = _create_population_stats_template(
        CONFIG['NumRuns'],
        include_equilibrium_tracking=True
  )
  CONFIG['overall_stats']['w'] = _create_population_stats_template(
        CONFIG['NumRuns'],
        include_equilibrium_tracking=True
  )

  if CONFIG['UseOverallRunStats']:
    CONFIG['overall_stats']['Overall'] = _create_population_stats_template(
         CONFIG['NumRuns'], 
         include_equilibrium_tracking=True
    )
    if CONFIG['GetLifespanDistribution']:
        CONFIG['overall_stats']['Overall']['all_ls'] = [] # all lifespans


  if CONFIG['GetLifespanDistribution']:
    CONFIG['overall_stats']['w']['all_ls'] = [] # all lifespans
    CONFIG['overall_stats']['W']['all_ls'] = [] # all lifespans

  return CONFIG['overall_stats']

def age_and_mortality_dist(CONFIG, pop_name, population, death_or_tick, dist_name):
    '''Updates values for age or mortality distribution in the overall_stats dictionary
    
    Args:
        CONFIG: Configuration dictionary
        pop_name: Name of the population (w, W, Overall)
        population: List of NMRs
        death_or_tick: Nmr death tick or current tick
        dist_name: Name of the distribution (age, mortality)

    '''
    
    for nmr in population:
      CONFIG['overall_stats'][pop_name][dist_name][abs(death_or_tick - nmr.birth) ] += 1

def NMR_birth(CONFIG, cur_tick, mom_w_alleles, dad_w_alleles, population = False, use_litter_size = True):
    '''Initializes new NMR'''
    new_NMRs = []
    litter_size = CONFIG['LitterSize'] if use_litter_size else 1

    if 'UseImmortalQueenMode' in CONFIG and CONFIG['UseImmortalQueenMode']:
      # In immortal queen mode, the mother and father have heterozygus alleles for female offspring 
        mom_w_alleles = ['W', 'w']
        dad_w_alleles = ['W', 'w']

    
    for l in range(litter_size):
    # Choses diploid allele to be inherited from each parent
      sex = 'F' if CONFIG['num_NMRs'] % 2 == 0 else 'M'
      parent_alleles = dict()
      

      parent_alleles['mW'] = mom_w_alleles[choice([0,1])]
      parent_alleles['dW'] = dad_w_alleles[choice([0,1])] if not (sex == 'M' and CONFIG['UseImmortalQueenMode']) else parent_alleles['mW']       # In immortal queen mode, the male NMRs will be haploid and have a 50:50 chance of being V or L


      new_NMR = Naked_Mole_Rat(CONFIG['num_NMRs'], cur_tick, parent_alleles, CONFIG)
      


    # add new NMRs to list of alive NMRs

      if new_NMR.sex == 'F':
        CONFIG['alive_NMRS']['female'].append(new_NMR)
      else:
        CONFIG['alive_NMRS']['male'].append(new_NMR)
      
      if population != False and CONFIG.get('GetIndividualRunData', False):
          population.append(new_NMR)

      CONFIG['num_NMRs'] += 1
      new_NMRs.append(new_NMR)

    return new_NMRs
  
def set_random_NMR_age(active_pop_name, new_nmr, CONFIG, run):
    ''' 
    Uses Binary Search To find a random NMR Age on the provided age distribution and ages nmr to that age 
    
    Args:
        active_pop_name: Name of the population
        new_nmr: NMR object
        CONFIG: Configuration dictionary
        run: Run number
        
    '''
    
    # Uses Binary Search To find a random NMR Age on the provided age distribution
    if len(CONFIG['age_distribution']['CDF'][active_pop_name]) == 0:
      # handle empty case
      print(f'Warning: No age distribution data found for {active_pop_name}')
      
    else:
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

      new_nmr.age_NMR(CONFIG, run, curr_tick=0, age=random_age, colony_population_size=len(CONFIG['alive_NMRS']['female']) + len(CONFIG['alive_NMRS']['male']))

def init_population(CONFIG, run, population = False):
    '''
    Initializes population to config specifications
    
    Args:
        CONFIG: Configuration dictionary
        run: Run number
        population: Optional population list
        
    Returns:
        List of NMRs
    '''
    
    if CONFIG['UseAvgForLogCap']:
      CONFIG['last_3_pop_sizes'] = [CONFIG['InitialPopulation']] * 3
    
    while ((len(CONFIG['alive_NMRS']['female']) + len(CONFIG['alive_NMRS']['male'])) < CONFIG['InitialPopulation'] or len(CONFIG['alive_NMRS']['female']) < 1 or len(CONFIG['alive_NMRS']['male']) < 1):

        # Randomizes Parent Alleles
        mom_w_alleles = choices(['w','W'], weights=[float(CONFIG['InitialRecessiveAlleleFraction']), 1-float(CONFIG['InitialRecessiveAlleleFraction'])], k = 2)
        dad_w_alleles = choices(['w','W'], weights=[float(CONFIG['InitialRecessiveAlleleFraction']), 1-float(CONFIG['InitialRecessiveAlleleFraction'])], k = 2) 

        new_nmrs = NMR_birth(CONFIG, 0,  mom_w_alleles, dad_w_alleles, population, use_litter_size = False)

        # Sets the age using the corresponding age distribution - please note, this is not implemented for use on both population types simultaneously
        for new_nmr in new_nmrs:
          if not CONFIG['GetAgeDistribution'] and ('NotUsingAgeDistOrLongRun' not in CONFIG or CONFIG['NotUsingAgeDistOrLongRun'] == False) and not CONFIG['DoNotUseAgeDistEver']:

              # Gets the type that will be used to determine which age distribution is being used
              pop_name = new_nmr.w_pops['type_name']

              set_random_NMR_age(pop_name, new_nmr, CONFIG, run)

          else:
            # Ages all NMRs to reproductive age without insults
            new_nmr.birth = -CONFIG['TicksTillFertile']


    if CONFIG['UseLinearGrowth']:
        CONFIG['queen'] = choice(CONFIG['alive_NMRS']['female'])


    return CONFIG['alive_NMRS']


def apply_insults(nmr_list, insult_tick, CONFIG, run, colony_population_size, insult_multiplier=1):
  '''Obtains and applies insults for each NMR in list'''
  for nmr in list(nmr_list):
      nmr.apply_insult(CONFIG, colony_population_size, insult_tick, run, insult_multiplier=insult_multiplier)

def init_colonies(CONFIG, run, population=False):
    """
    Initialize N0 colonies with heterozygous Queens
    
    Each original colony starts with a heterozygous Queen:
    - Genotype: (HP-R)/(hp-r) = (W)/(w)
    - Stored sperm: 50% (W), 50% (w)
    
    Args:
        CONFIG: Configuration dictionary
        run: Run number
        population: Optional population list for tracking
        
    Returns:
        List of Colony objects
    """
    colonies = []
    n0_colonies = CONFIG.get('N0Colonies', 5)

    # This means one allele is W (dominant, HP-R) and one is w (recessive, hp-r)
    heterozygous_alleles = dict(
        mW='W',  # Mother allele: W (W -> Use HP and R values from CONFIG)
        dW='w'   # Father allele: w (w -> Use hp and r values from CONFIG)
    )    
    for colony_id in range(n0_colonies):
        # Create heterozygous Queen: (HP-R)/(hp-r) = (W)/(w)

        # Create Queen NMR
        # The queen should be female and the id should be even
        if CONFIG['num_NMRs'] % 2 != 0:
            CONFIG['num_NMRs'] += 1

        queen = Naked_Mole_Rat(CONFIG['num_NMRs'], 0, heterozygous_alleles, CONFIG)

        # CONFIG['num_NMRs'?] += 1
        if not CONFIG['MatesWithinColony']:
            queen.stored_sperm_alleles = ['W', 'w']  # 50% (W), 50% (w)
        # CONFIG['num_NMRs'] += 1 # Increment the number of NMRs #TODO: is this being done twice?
        
        # Add Queen to population tracking if needed
        if population != False and CONFIG.get('GetIndividualRunData', False):
            population.append(queen)
        
        # Create initial colony members (start with a few workers)
        # For original colonies, we'll start with some initial members
        initial_members = []
        initial_colony_size = CONFIG.get('InitialPopulation', 10) // n0_colonies
        initial_colony_size = max(initial_colony_size, 2)  # At least 2 members
        
        for i in range(initial_colony_size - 1):  # -1 because Queen is already counted

            # For original colonies, members are heterozygous            
            CONFIG['num_NMRs'] += 1
            member = Naked_Mole_Rat(CONFIG['num_NMRs'], 0, heterozygous_alleles, CONFIG)

            # Set age if not using age distribution
            if not CONFIG.get('GetAgeDistribution', False) and ('NotUsingAgeDistOrLongRun' not in CONFIG or CONFIG['NotUsingAgeDistOrLongRun'] == False) and not CONFIG['DoNotUseAgeDistEver']:
                pop_name = member.w_pops['type_name']
                set_random_NMR_age(pop_name, member, CONFIG, run)
            else:
                member.birth = -CONFIG['TicksTillFertile']
            
            initial_members.append(member)
            
            if population != False and CONFIG.get('GetIndividualRunData', False):
                population.append(member)
        
        # Create colony
        colony = Colony(queen, initial_members, colony_id, 0)
        colonies.append(colony)
    
    return colonies


def calculate_colony_formation_probability(num_colonies, CONFIG):
    """
    Calculate probability of new colony formation using logistic equation
    
    Args:
        num_colonies: Current number of colonies
        CONFIG: Configuration dictionary
        
    Returns:
        Probability (0.0 to 1.0)
    """
    base_prob = CONFIG.get('BaseColonyFormationProbability', 0.05)
    logistic_mult = CONFIG.get('ColonyLogisticMult', 0.1)
    logistic_exp = CONFIG.get('ColonyLogisticExp', 1.0)
    
    # Logistic equation: probability = base_prob / exp(logistic_mult * num_colonies^logistic_exp)
    probability = float(base_prob) / exp(float(logistic_mult) * (float(num_colonies) ** float(logistic_exp)))
    
    # print('num cols at calculation of formation prob', num_colonies, 'prob', probability)
    return min(probability, 1.0)  # Cap at 1.0

def conduct_mating_flight(colonies, tick, CONFIG, run, population=False):
    """
    Conduct mating flight: colonies contribute individuals who mate randomly
    
    Args:
        colonies: List of Colony objects
        tick: Current tick
        CONFIG: Configuration dictionary
        population: Optional population list
        
    Returns:
        List of fertilized females (potential new Queens)
    """
    threshold = CONFIG.get('ThresholdColonySize', 30)
    selection_method = CONFIG.get('MatingFlightSelection', 'random')
    
    # Collect contributors from all colonies
    all_contributors = []
    contributor_colonies = {}  # Track which colony each contributor belongs to

    for colony in colonies:
        colony_size = colony.get_size()
        num_to_contribute = max(0, colony_size - threshold)
        
        if num_to_contribute > 0:
            contributors = colony.get_contributors(num_to_contribute, selection_method)
            all_contributors.extend(contributors)
            
            # Track contributors for removal
            for contributor in contributors:
                contributor_colonies[contributor] = colony

            # Update Stats
            colony.num_nmrs_contributed_to_mating_pool += num_to_contribute
            colony.num_times_contributed_to_mating_pool += 1

        colony.pct_of_colony_members_contributed_to_mating_pool_each_time.append(num_to_contribute / colony_size)
    
    # Separate into males and females
    females = [c for c in all_contributors if c.sex == 'F']
    males = [c for c in all_contributors if c.sex == 'M']
    
    # Random mating: pair random females with random males
    shuffle(females)
    shuffle(males)
     # Pair up (limited by whichever sex runs out first)
    num_pairs = min(len(females), len(males))
    # print(f'Mating, there are {len(contributor_colonies)} contributing colonies and {len(all_contributors)} nmrs mating')
    
    use_stored_sperm = True if CONFIG['PioneerGroupSize'] == 1  or not CONFIG['MatesWithinColony'] else False # In this case, the queen will start the colony alone using stored sperm
    size_of_pioneer_group = CONFIG['PioneerGroupSize']

    if CONFIG['PioneerGroupSize'] % 2 != 0:
        if not use_stored_sperm:
            print(f'Warning: PioneerGroupSize is not even, it will be rounded up to the nearest even number')
        size_of_pioneer_group = CONFIG['PioneerGroupSize'] + 1

    num_pioneer_groups = num_pairs // ( size_of_pioneer_group // 2 )
    num_per_gender = size_of_pioneer_group // 2
    pioneer_groups = []   
    all_pioneer_group_members = []
    
    for i in range(num_pioneer_groups):
        colony_females = females[i* num_per_gender:(i+1)* num_per_gender]
        colony_males = males[i* num_per_gender:(i+1)* num_per_gender]

        # This represents the sperm she received from the mating flight
        if use_stored_sperm:
          colony_females[0].stored_sperm_alleles = colony_males[0].w_pops['alleles']
          pioneer_groups.append(colony_females)
          all_pioneer_group_members.extend(colony_females)
        else:
          pioneer_groups.append(colony_females + colony_males)
          all_pioneer_group_members.extend(colony_females + colony_males)

    # Remove contributors from their colonies
    for contributor, colony in contributor_colonies.items():
        colony.remove_member(contributor)
        #   If they have not joined a pioneer group, they die
        if contributor not in all_pioneer_group_members:
            contributor.kill_nmr(CONFIG, tick, run)

        elif use_stored_sperm and CONFIG['PioneerGroupSize'] == 1:
            # If we use stored sperm, the male's die after they have mated
            if contributor.sex == 'M':
                contributor.kill_nmr(CONFIG, tick, run)

    return pioneer_groups

def form_new_colony(fertilized_female, tick, CONFIG, colony_id, colony_members, population=False):
    """
    Form a new colony from a fertilized female
    
    Args:
        fertilized_female: Naked_Mole_Rat object (new Queen)
        tick: Current tick
        CONFIG: Configuration dictionary
        colony_id: Unique colony ID
        colony_members: List of NMRs to add to the colony
        population: Optional population list
        
    Returns:
        Colony object or None if formation failed
    """
    # Create new colony with Queen
    colony = Colony(fertilized_female, colony_members, colony_id, tick)
    
    return colony