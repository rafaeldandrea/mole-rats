# Standard library imports
from random import choice, sample, shuffle
from statistics import mean

# Local imports
from Naked_Mole_Rat import Naked_Mole_Rat

class Colony:
    """Represents a colony with a Queen and members"""
    
    def __init__(self, queen, members=None, colony_id=0, established_on=0):
        """
        Initialize a colony
        
        Args:
            queen: Naked_Mole_Rat object representing the Queen
            members: List of Naked_Mole_Rat objects (excluding Queen)
            colony_id: Unique identifier for the colony
        """
        self.queen = queen
        self.members = []
        self.colony_id = colony_id
        self.established_on = established_on

        self.num_nmrs_contributed_to_mating_pool = 0
        self.num_times_contributed_to_mating_pool = 0
        self.pct_of_colony_members_contributed_to_mating_pool_each_time = []
        self.is_homozygous = False

        for member in members:
            self.add_member(member)

        self.last_queen_death_tick = None
        
    def get_size(self):
        """Get total colony size including Queen"""
        return len(self.members) + 1
    
    def is_alive(self, CONFIG, tick, run):
        """Check if colony is alive and kill it if it isn't"""
        ouput = True

        if self.queen is None:
            if CONFIG['TicksToQueenReplacement'] is None:
                ouput = False
            elif len(self.get_females()) == 0:
                ouput = False
            
            if CONFIG['MatesWithinColony']:
                if len(self.get_males()) == 0:
                    ouput = False
        elif self.queen.death is not None:
            ouput = False
        elif self.queen.stored_sperm_alleles == None and len(self.get_males()) == 0:
            ouput = False
        
        if not ouput:
            self.kill_colony(CONFIG, tick, run)

        return ouput

    def get_contributors(self, num_to_contribute, selection_method='random'):
        """
        Select individuals to contribute to mating flight
        
        Args:
            num_to_contribute: Number of individuals to select
            selection_method: 'random', 'youngest', or 'oldest'
            
        Returns:
            List of selected NMRs
        """
        if num_to_contribute <= 0 or len(self.members) == 0:
            return []
        
        num_to_contribute = min(num_to_contribute, len(self.members))
        
        if selection_method == 'random':
            return sample(self.members, num_to_contribute)
        elif selection_method == 'youngest':
            # Sort by birth tick (higher birth = younger)
            sorted_members = sorted(self.members, key=lambda x: x.birth, reverse=True)
            return sorted_members[:num_to_contribute]
        elif selection_method == 'oldest':
            # Sort by birth tick (lower birth = older)
            sorted_members = sorted(self.members, key=lambda x: x.birth)
            return sorted_members[:num_to_contribute]
        else:
            # Default to random
            return sample(self.members, num_to_contribute)
    
    def remove_contributors(self, contributors):
        """
        Remove contributors from colony members
        
        Args:
            contributors: List of NMRs to remove
        """
        for contributor in contributors:
            if contributor in self.members:
                self.members.remove(contributor)
    
    def get_females(self):
        """Get all female members (excluding Queen)"""
        return [m for m in self.members if m.sex == 'F']
    
    def get_males(self):
        """Get all male members (excluding Queen)"""
        return [m for m in self.members if m.sex == 'M']
    
    def add_member(self, nmr):
        """Add a member to the colony"""
        if nmr != self.queen:
            self.members.append(nmr)
            nmr.colony_ids.append(self.colony_id)
    
    def remove_member(self, nmr):
        """Remove a member from the colony"""
        if nmr in self.members:
            self.members.remove(nmr)
    
    def get_is_homozygous(self):
        '''Check if the colony is homozygous, having the same allele in all members
        
        Args:
            None
        
        Returns:
            True if the colony is homozygous, False otherwise
        '''
        if self.is_homozygous:
            return True
        else:
            if len(self.members) == 0 and self.queen is None:
                return False
            comparison_nmr = self.queen if self.queen is not None else self.members[0]
            if comparison_nmr.w_pops['alleles'][0] == comparison_nmr.w_pops['alleles'][1]:
                if comparison_nmr.stored_sperm_alleles is not None and self.queen is not None:
                    if not comparison_nmr.w_pops['alleles'] == comparison_nmr.stored_sperm_alleles:
                        return False
                
                for member in self.members:
                    if not comparison_nmr.w_pops['alleles'] == member.w_pops['alleles']:
                        return False
                        
                self.is_homozygous = True
                return True
            else:
                return False
    
    def get_fixed_allele_name(self):
        if self.get_is_homozygous():
            if self.queen is not None:
                return self.queen.w_pops['alleles'][0]
            else:
                if len(self.members) > 0:
                    return self.members[0].w_pops['alleles'][0]
                else:
                    return None
        else:
            return None
    
    def kill_colony(self, CONFIG, tick, run):
        for factor in CONFIG['display_factors']:
            if factor == 'Overall' or (self.queen is not None and factor == self.queen.w_pops['type_name']):
                if self.queen is not None and CONFIG['TicksToQueenReplacement'] is None:
                    CONFIG['colony_stats']['queen_lifespan'][factor].append(tick - self.queen.birth)
                CONFIG['colony_stats']['colony_lifespan'][factor].append(tick - self.established_on)
                
                CONFIG['colony_stats']['num_nmrs_colony_contributed_to_mating_pool'][factor].append(self.num_nmrs_contributed_to_mating_pool)
                CONFIG['colony_stats']['num_times_colony_contributed_to_mating_pool'][factor].append(self.num_times_contributed_to_mating_pool)
                if len(self.pct_of_colony_members_contributed_to_mating_pool_each_time ) > 1:
                    CONFIG['colony_stats']['pct_of_colony_members_contributed_to_mating_pool_each_time'][factor].append(mean(self.pct_of_colony_members_contributed_to_mating_pool_each_time ))

        # Kill the queen (ensure she's removed from alive_NMRS)
        if self.queen and self.queen.death is None:
            self.queen.kill_nmr(CONFIG, tick, run)
        
        # Kill all members
        for member in list(self.members):
            member.kill_nmr(CONFIG, tick, run)
        
    def __str__(self):
        return f"Colony {self.colony_id}: Queen {self.queen.id}, {len(self.members)} members, Size: {self.get_size()}"

