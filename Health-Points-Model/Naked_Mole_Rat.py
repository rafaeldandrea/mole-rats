# Local imports
from get import insult_severity

class Naked_Mole_Rat:

    def get_dominant_trait(self, dominant_trait, recessive_trait, mom_allele, dad_allele):
      ''' Returns dominant allele if the dominant trait exists '''
      return dominant_trait if mom_allele == dominant_trait or dad_allele == dominant_trait else recessive_trait

    def age_NMR(self, CONFIG, run, curr_tick, colony_population_size, age=False):
      ''' Ages NMR up to either input age or death'''
      if age:
        self.birth = self.birth - age
        for tick in range(age):
            if self.death == None:
                self.apply_insult(CONFIG, colony_population_size, -tick, run, after_pop_death=True)
      else:
        while self.death == None:
          self.apply_insult(CONFIG, colony_population_size, curr_tick, run, after_pop_death=True)
          curr_tick += 1
          
    def kill_nmr(self, CONFIG, insult_tick, run):
        '''Sets death to insult_tick and updates corresponding stats'''

        sex = 'female' if self.sex == 'F' else 'male'

        if self.id in [nmr.id for nmr in CONFIG['alive_NMRS'][sex]]:
          CONFIG['alive_NMRS'][sex].remove(self)
        else:
          if CONFIG['ShowIndividualRunStats']:
            print(f'Warning: NMR {self.id} not found in CONFIG["alive_NMRS"][{sex}] on tick {insult_tick}')
                
        if self.death is not None:
          if CONFIG['ShowIndividualRunStats']:
            print(f'Warning: NMR {self.id} is already dead on tick {insult_tick}')
          return
        self.death = insult_tick

        if (CONFIG['GetLifespanDistribution'] or CONFIG['GetMortality']) and abs(self.death - self.birth) > CONFIG['MaxAge']:
          print('Warning: RAISE MAX AGE')
          
        if CONFIG['GetLifespanDistribution'] and self.death > 0 and self.birth > 0:

          if CONFIG['UseOverallRunStats']:
              CONFIG['overall_stats']['Overall']['lifespan_distribution'][abs(self.death - self.birth)] += 1
              
          CONFIG['overall_stats'][self.w_pops['type_name']]['lifespan_distribution'][abs(self.death - self.birth)] += 1
              
        if CONFIG['GetMortality'] and CONFIG['overall_stats']['pop_death_ticks'][run] > 0 and self.birth > 100 and self.birth + 100 < CONFIG['MaxTicks']:

          if CONFIG['UseOverallRunStats']:
              CONFIG['overall_stats']['Overall']['equilibrium_lifespan_distribution'][abs(self.death - self.birth)] += 1
          
          CONFIG['overall_stats'][self.w_pops['type_name']]['equilibrium_lifespan_distribution'][abs(self.death - self.birth)] += 1



    def apply_insult(self, CONFIG, colony_population_size, insult_tick, run, insult_multiplier=1, after_pop_death=False, is_queen=False):
      ''' Gets and applies insult to NMR 
      
      Args:
          CONFIG: Configuration dictionary
          insult_tick: Current tick
          run: Run number
          after_pop_death: Whether population has died
          is_queen: Whether this NMR is a Queen (applies HP multiplier)
      '''
      total_insults = sum(self.current_insults)
    
      if CONFIG['StoreNMRData'] and self.id % 100 == 0:
        
        id_str = str(self.id) +',' + str(abs(self.birth -insult_tick)) +',' +str(self.current_insults) +',' +str(self.current_health_points)+',' +str(self.birth)+','+str(insult_tick)+ ','+str(run)+ ','+str(self.r_value_current)+'\n'
        CONFIG['writer'].write(id_str)

      # Updates Health Points
      if self.current_health_points - total_insults <= 0:
        # Removes Dead NMRs from the Population who died of an insult or naturally
          self.kill_nmr(CONFIG, insult_tick, run)
          
      # Simulates increased frailty due to aging reducing health points
      age = abs(int(self.birth) - int(insult_tick))
      try:
        if self.r_value_current == 0:
          print(f'Warning: r value is zero for NMR {self.id}, using default value of 1')
          self.r_value_current = 1
        self.current_health_points = CONFIG[self.hp_type_name] * (.5**(age / self.r_value_current))
      except ZeroDivisionError:
        print(f'Error: Division by zero when calculating health points for NMR {self.id}')
        self.current_health_points = CONFIG[self.hp_type_name]
      except KeyError as e:
        print(f'Error: Missing CONFIG key {e} when calculating health points for NMR {self.id}')
        self.current_health_points = 0

      if is_queen:
        self.current_health_points *= CONFIG['QueenHPMultiplier']

      # Removes recovered insults from current insults list and adds newly generated insult
      self.current_insults = self.current_insults[1:] + [insult_severity(CONFIG, colony_population_size, nmr=self, tick=insult_tick, after_pop_death=after_pop_death) * insult_multiplier]
      

    def store_allele_info(self, dom_name, rec_name, inherited_alleles):
      '''Makes allele info into accessable format'''

      mom_allele = dom_name if inherited_alleles['mW'] == 'W' else rec_name
      dad_allele = dom_name if inherited_alleles['dW'] == 'W' else rec_name
    
      active_allele = self.get_dominant_trait(dom_name, rec_name, mom_allele, dad_allele)
      allele_info = dict(
        alleles = [mom_allele, dad_allele],
        allele_name = mom_allele + dad_allele if mom_allele == dad_allele else rec_name + dom_name,
        type_name = active_allele
        )

      return allele_info

    def __str__(self):
        ''' This is used if you try to print an NMR'''
        id_str = 'NMR ' + str(self.id) + ' ' + str(self.sex) + ' ' + str(self.colony_ids) + \
            '\nBorn: ' + str(self.birth) + \
            '\nDied: ' + str(self.death ) + \
            '\nCurrent Insults: ' + str(self.current_insults) + \
            '\nHealth Points: ' + str(self.current_health_points)
        return id_str

    def __init__(
        self,
        id,
        birth,
        inherited_alleles,
        CONFIG
        ):
        '''Initializes NMR'''
        self.id = id
        self.sex = 'F' if id % 2 == 0 else 'M'
        self.birth = birth
        self.death = None  # When death is determined it will be set to a positive number
        self.current_insults = [0] * (CONFIG['insult_recovery_ticks'] + 1)  # This will be a list of healing insults where the index coresponds to the number of remaining ticks

        self.w_pops = self.store_allele_info('W', 'w', inherited_alleles)

        self.hp_type_name = 'HP' if self.w_pops['type_name'] == 'W' else 'hp'
        self.current_health_points = CONFIG[self.hp_type_name]

        r_type_name = 'R' if self.w_pops['type_name'] == 'W' else 'r'
        self.r_value_current = CONFIG[r_type_name]
        
        # This is set when creating a Queen, not in __init__
        self.stored_sperm_alleles = None  # Will be ['W', 'w'] for original colonies
        self.colony_ids = [] # Will be a list of colony ids that this NMR was a member of, the last one will be the one it is currently in



