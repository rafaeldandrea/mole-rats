
### Classes ###
from get import insult_severity


class Naked_Mole_Rat:

    def get_dominant_trait(self, dominant_trait, recessive_trait, mom_allele, dad_allele):
      ''' Returns dominant allele if the dominant trait exists '''
      return dominant_trait if mom_allele == dominant_trait or dad_allele == dominant_trait else recessive_trait

    def age_NMR(self, CONFIG, run, curr_tick, age=False):
      ''' Ages NMR up to either input age or death'''
      if age:
        self.birth = self.birth - age
        for tick in range(age):
            if self.death == None:
                # print(tick, age, self)
                self.apply_insult(CONFIG, -tick, run, after_pop_death=True)
      else:
        while self.death == None:
          # print('age_nmr')
          self.apply_insult(CONFIG, curr_tick, run, after_pop_death=True)
          curr_tick += 1
          
    def kill_nmr(self, CONFIG, insult_tick, run):
        '''Sets death to insult_tick and updates corresponding stats'''
        self.death = insult_tick
        sex = 'female' if self.sex == 'F' else 'male'
        # print (sex, len(CONFIG['alive_NMRS'][sex]))
        # print ('    ', [type(i) == Naked_Mole_Rat for i in CONFIG['alive_NMRS'][sex]])
        # print ('    ', [i.id for i in CONFIG['alive_NMRS'][sex]])
        # print ('    ', [i.birth for i in CONFIG['alive_NMRS'][sex]])
        # print ('    ', [i.death for i in CONFIG['alive_NMRS'][sex]])

        CONFIG['alive_NMRS'][sex].remove(self)
        # print(self)
        # if (abs(insult_tick - self.birth) < 20):
        #   print(abs(self.death - self.birth) , '=',self.death, '-' ,self.birth, '   ', insult_tick)
        if (CONFIG['GetLifespanDistribution'] or CONFIG['GetMortality']) and abs(self.death - self.birth) > CONFIG['MaxAge']:
          print('Warning: RAISE MAX AGE')
          
        if CONFIG['GetLifespanDistribution'] and self.death > 0 and self.birth > 0:

          if CONFIG['UseOverallRunStats']:
              # print(abs(self.death - self.birth) , '=',self.death, '-' ,self.birth, '   ', insult_tick)
              CONFIG['overall_stats']['Overall']['lifespan_distribution'][abs(self.death - self.birth)] += 1
              
          if CONFIG['show_hp']:
              CONFIG['overall_stats'][self.hp_pops['type_name']]['lifespan_distribution'][abs(self.death - self.birth)] += 1

          if CONFIG['show_r']:
              CONFIG['overall_stats'][self.r_pops['type_name']]['lifespan_distribution'][abs(self.death - self.birth)] += 1
              
        if CONFIG['GetMortality'] and CONFIG['overall_stats']['pop_death_ticks'][run] > 0 and self.birth > 100 and self.birth + 100 < CONFIG['MaxTicks']:

          if CONFIG['UseOverallRunStats']:
              CONFIG['overall_stats']['Overall']['equilibrium_lifespan_distribution'][abs(self.death - self.birth)] += 1

          if CONFIG['show_hp']:
              CONFIG['overall_stats'][self.hp_pops['type_name']]['equilibrium_lifespan_distribution'][abs(self.death - self.birth)] += 1

          if CONFIG['show_r']:
              CONFIG['overall_stats'][self.r_pops['type_name']]['equilibrium_lifespan_distribution'][abs(self.death - self.birth)] += 1



    def apply_insult(self, CONFIG, insult_tick, run, after_pop_death=False):
      ''' Gets and applies insult to NMR '''
      total_insults = sum(self.current_insults)
    
      if CONFIG['StoreNMRData'] and self.id % 100 == 0:
        
        id_str = str(self.id) +',' + str(abs(self.birth -insult_tick)) +',' +str(self.current_insults) +',' +str(self.hp_pops['value'])+',' +str(self.birth)+','+str(insult_tick)+ ','+str(run)+ ','+str(self.r_pops['value'])+'\n'
        CONFIG['writer'].write(id_str)
        # print(id_str)
      # print('insult_tick',type(insult_tick))
      # print('Killing ', type(self.birth))       
      # Updates Health Points

      if self.hp_pops['value'] - total_insults <= 0:
        # Removes Dead NMRs from the Population who died of an insult or naturally
          self.kill_nmr(CONFIG, insult_tick, run)
          
      # Simulates increased frailty due to aging reducing health points
      age = abs(int(self.birth) - int(insult_tick))
      self.hp_pops['value'] = CONFIG[self.hp_pops['type_name']] * (.5**(age / self.r_pops['value']))
      #   self.hp_pops['value'] = self.hp_pops['value'] - self.r_pops['value']

      # Removes recovered insults from current insults list and adds newly generated insult
      self.current_insults = self.current_insults[1:] + [insult_severity(CONFIG, nmr=self, tick=insult_tick, after_pop_death=after_pop_death)]
      
  # print(count)

    def store_allele_info(self, dom_name, rec_name, inherited_alleles, use_parents, use_williams_allele):
      '''Makes allele info into accessable format'''
      if use_parents:
        if use_williams_allele:
          mom_allele = dom_name if inherited_alleles['mW'] == 'W' else rec_name
          dad_allele = dom_name if inherited_alleles['dW'] == 'W' else rec_name
        else:
          mom_allele = inherited_alleles['m'+dom_name]         
          dad_allele = inherited_alleles['d'+dom_name] 
      
      if use_parents:
        active_allele = self.get_dominant_trait(dom_name, rec_name, mom_allele, dad_allele)
        allele_info = dict(
          alleles = [mom_allele, dad_allele],
          allele_name = mom_allele + dad_allele if mom_allele == dad_allele else rec_name + dom_name,
          type_name = active_allele
          )
      else: 
        allele_info = dict(
          alleles = [rec_name, rec_name],
          allele_name = rec_name + rec_name,
          type_name = rec_name
          )
        
      return allele_info

    def __str__(self):
        ''' This is used if you try to print an NMR'''
        id_str = 'NMR ' + str(self.id) + \
            '\nBorn: ' + str(self.birth) + \
            '\nDied: ' + str(self.death ) + \
            '\nCurrent Insults: ' + str(self.current_insults) + \
            '\nHealth Points: ' + str(self.hp_pops['value'])
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

        self.w_pops = self.store_allele_info('W', 'w', inherited_alleles, CONFIG['UseTradeoffAllele'], CONFIG['UseTradeoffAllele'])

        self.hp_pops = self.store_allele_info('HP', 'hp', inherited_alleles, CONFIG['show_hp'], CONFIG['UseTradeoffAllele'])
        self.hp_pops['value'] = CONFIG[self.hp_pops['type_name']]

        self.r_pops = self.store_allele_info('R', 'r', inherited_alleles, CONFIG['show_r'], CONFIG['UseTradeoffAllele'])
        self.r_pops['value'] = CONFIG[self.r_pops['type_name']]




