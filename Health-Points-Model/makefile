#### Makefile ####
# NMR Simulator - Quick Reference
#   make runDefault   Run default: Individual, Ant, Mole-Rat modes (R9r9.2, R9r9.3, R9r10)
#   make quick        Fast test run (2 runs, Lin only)
#   make help         Show all targets and options
#   Override options: make run DATA_COMBINATIONS=R10r9 REPRODUCTION_TYPE=Exp NUM_RUNS=20

.PHONY: default run runDefault setRunDetails help quick GetInitialDistributionData GetAgeDistributionData cleanCopy

help:
	@echo "NMR Simulator - Available targets:"
	@echo ""
	@echo "  make default             Default: run Individual, Ant, and Mole-Rat modes"
	@echo "  make runDefault    Run Individual, Ant, and Mole-Rat modes (R9r9.2, R9r9.3, R9r10 configs)"
	@echo "  make quick        Fast test (2 runs, Lin only) - use to verify setup"
	@echo "  make GetInitialDistributionData  Equilibrium + age + lifespan (LONG_RUN=True)"
	@echo "  make GetAgeDistributionData      Age distribution only"
	@echo ""
	@echo "  make run          Direct run (requires NUM_RUNS=, REPRODUCTION_TYPE=)"
	@echo "  make cleanCopy    Copy project to ../NMR_SIM_Copy"
	@echo ""
	@echo "Common overrides:"
	@echo "  DATA_COMBINATIONS=R9r10HP1000hp900pcE1.1pcM0.07BP0.3TTF3IRT2"
	@echo "  REPRODUCTION_TYPE=Lin | Exp   (omit for both)"
	@echo "  NUM_RUNS=50  RUN_MULTIPLIER=10  LONG_RUN=True"

RUN_MULTIPLIER ?= 1
LONG_RUN ?= True
GET_AGE_DIST ?= False
DO_NOT_USE_AGE_DIST_EVER ?= True
USE_OVERALL_RUN_STATS = False
WRITE_DATA_TO_CSV = False
USE_AVG_POP_SIZE_FOR_LOG_CAP = False

SHOW_INDIV_RUN_STATS ?= True
SHOW_INDIV_RUN_GRAPH ?= False
MAKE_DISTRIBUTION_GRAPHS ?= True
INIT_REC_ALLELE_FRACTION = 0.5 #.75 means that 75% of alleles will be recessive in the beginning, .285 means 28.8% and so on.
SET_INIT_POP_TO_EQUILIBRIUM = False
INIT_POP = 16 # Only if SET_INIT_POP_TO_EQUILIBRIUM = False will this be used
DATA_COMBINATIONS ?= R9r9.2HP1000hp900pcE1.1pcM0.07BP0.3TTF3IRT2 #R9r9.3HP1000hp600pcE1.1pcM0.07BP0.3TTF3IRT2 # R9r10HP1000hp900pcE1.1pcM0.07BP0.3TTF3IRT2 
MAX_TICKS = 10000
IMMORTAL_QUEEN_MODE = False

# Colony Reproduction Mode settings
USE_COLONY_REPRODUCTION_MODE = False
# The following settings are only used if USE_COLONY_REPRODUCTION_MODE is True
########################################################
N0_COLONIES = 16
MATING_FLIGHT_INTERVAL = 6
THRESHOLD_COLONY_SIZE = 23
MATING_FLIGHT_SELECTION = 'youngest' # 'random' or 'youngest' or 'oldest'
QUEEN_HP_MULTIPLIER = 1
BASE_COLONY_FORMATION_PROBABILITY = 1
COLONY_LOGISTIC_MULTIPLIER = 0.06
COLONY_LOGISTIC_EXP = 1.1

# Starting Inbreeding and HP Simulation Settings
# This should have the following settings:
# The following settings are only used if USE_COLONY_REPRODUCTION_MODE is True
PIONEER_GROUP_SIZE = 1 # The number of NMRs in the pioneer group, this should be either 1 (if the queen has stored sperm and you want her to start the colony alone) or an even number to make it half female, half male
COLONY_DISASTER_PROBABILITY = 0 # The probability of a colony disaster occurring for a colony each tick
COLONY_DISASTER_INSULT_MULTIPLIER = 1 # The multiplier for the insult severity on the tick of a colony disaster
TICKS_TO_QUEEN_REPLACEMENT = None # If None, the queen will not be replaced and on the queen's death, all members of the the colony will die
MATES_WITHIN_COLONY = False # Whether to allow mates within the colony or not (insead of using stored sperm), this is counted as True if USE_COLONY_REPRODUCTION_MODE is False


########################################################


# Ex. R10r9
# This will eiter create or use a config file with the default configuration

    # HP = "900" 
    # hp = "900" 
    # R = "10" 
    # r = "10" 
    # LitterSize = "1" (LS)
    # pop_cap_exponent = "1.2" (pcE)
    # pop_cap_mult = "0.1" (pcM)  
    # TicksTillFertile = "1" (TTF)
    # InsultDistribution = "T7" (IdT)
    # BirthProbability = "0.4" (BP)
    # insult_recovery_ticks = 0 (IRT)

# To use a different configuration, use the corresponding letter code 
# Ex. TTF2HP290pcE2IdT4LS5IRT2 will create and run the following configuration

    # HP = "290" 
    # hp = "900" 
    # R = "10" 
    # r = "10" 
    # LitterSize = "5" (LS)
    # pop_cap_exponent = "2" (pcE)
    # pop_cap_mult = "0.1" (pcM)  
    # TicksTillFertile = "2" (TTF)
    # InsultDistribution = "T4" (Id)
    # BirthProbability = "0.4" (BP)
    # insult_recovery_ticks = 2 (IRT)

# Ex. R10r9HP900hp1000

    # HP = "900" 
    # hp = "1000" 
    # R = "10" 
    # r = "9" 
    # LitterSize = "1" (LS)
    # pop_cap_exponent = "1.2" (pcE)
    # pop_cap_mult = "0.1" (pcM)  
    # TicksTillFertile = "1" (TTF)
    # InsultDistribution = "T7" (IdT)
    # BirthProbability = "0.4" (BP)
    # insult_recovery_ticks = 0 (IRT)


# Default target
default: runDefault

# Run simulation with default makefile settings
runSim:
	$(MAKE) setRunDetails LIN_RUN_NUM?=1 EXP_RUN_NUM?=1;

# Fast test run (2 runs, Lin only) — use to verify setup
quick:
	$(MAKE) setRunDetails LIN_RUN_NUM=2 EXP_RUN_NUM=0 REPRODUCTION_TYPE=Lin USE_COLONY_REPRODUCTION_MODE=False

run:
	$(foreach var,$(DATA_COMBINATIONS),\
		python3 nmr_simulator.py $(var)$(REPRODUCTION_TYPE) \
			UseColonyReproductionMode $(USE_COLONY_REPRODUCTION_MODE)\
			NumRuns $(NUM_RUNS) \
			LongRun $(LONG_RUN)\
			WriteDataToCSV $(WRITE_DATA_TO_CSV)\
			UseOverallRunStats $(USE_OVERALL_RUN_STATS)\
			MakeOverallGraphs $(MAKE_DISTRIBUTION_GRAPHS)\
			GetAgeDistribution $(GET_AGE_DIST)\
			DoNotUseAgeDistEver $(DO_NOT_USE_AGE_DIST_EVER)\
			ShowIndividualRunStats $(SHOW_INDIV_RUN_STATS)\
			ShowIndividualRunGraphs $(SHOW_INDIV_RUN_GRAPH)\
			UseAvgForLogCap $(USE_AVG_POP_SIZE_FOR_LOG_CAP)\
			UseEquilibriumForInitPop $(SET_INIT_POP_TO_EQUILIBRIUM)\
			InitialPopulation $(INIT_POP)\
			MaxTicks $(MAX_TICKS)\
			UseImmortalQueenMode $(IMMORTAL_QUEEN_MODE)\
			InitialRecessiveAlleleFraction $(INIT_REC_ALLELE_FRACTION)\
			N0Colonies $(N0_COLONIES) \
			MatingFlightInterval $(MATING_FLIGHT_INTERVAL)\
			ThresholdColonySize $(THRESHOLD_COLONY_SIZE)\
			MatingFlightSelection $(MATING_FLIGHT_SELECTION)\
			QueenHPMultiplier $(QUEEN_HP_MULTIPLIER)\
			BaseColonyFormationProbability $(BASE_COLONY_FORMATION_PROBABILITY)\
			ColonyLogisticMult $(COLONY_LOGISTIC_MULTIPLIER)\
			ColonyLogisticExp $(COLONY_LOGISTIC_EXP)\
			PioneerGroupSize $(PIONEER_GROUP_SIZE)\
			ColonyDisasterProbability $(COLONY_DISASTER_PROBABILITY)\
			ColonyDisasterInsultMultiplier $(COLONY_DISASTER_INSULT_MULTIPLIER)\
			MatesWithinColony $(MATES_WITHIN_COLONY)\
			TicksToQueenReplacement $(TICKS_TO_QUEEN_REPLACEMENT);\
		)

setRunDetails:

	@if [ "$(REPRODUCTION_TYPE)" = "" ]; then \
		$(MAKE) run REPRODUCTION_TYPE=Lin NUM_RUNS?=`echo $$(($(LIN_RUN_NUM) * $(RUN_MULTIPLIER)))`; \
		$(MAKE) run REPRODUCTION_TYPE=Exp NUM_RUNS?=`echo $$(($(EXP_RUN_NUM) * $(RUN_MULTIPLIER)))`; \
	elif [ "$(REPRODUCTION_TYPE)" = "Lin" ]; then \
		$(MAKE) run NUM_RUNS?=`echo $$(($(LIN_RUN_NUM) * $(RUN_MULTIPLIER)))`; \
	elif [ "$(REPRODUCTION_TYPE)" = "Exp" ]; then \
		$(MAKE) run NUM_RUNS?=`echo $$(($(EXP_RUN_NUM) * $(RUN_MULTIPLIER)))`; \
	fi ;\


runMoleRatMode:
	$(MAKE) setRunDetails DATA_COMBINATIONS="R9r10HP1000hp900pcE1.1pcM0.07BP0.3TTF3IRT2 \
			R9r9.3HP1000hp600pcE1.1pcM0.07BP0.3TTF3IRT2 R9r9.2HP1000hp900pcE1.1pcM0.07BP0.3TTF3IRT2" \
			USE_COLONY_REPRODUCTION_MODE=True N0_COLONIES=16 THRESHOLD_COLONY_SIZE=23 \
			MATING_FLIGHT_SELECTION='youngest' QUEEN_HP_MULTIPLIER=1 \
			COLONY_LOGISTIC_MULTIPLIER=0.06 COLONY_LOGISTIC_EXP=1.1 \
			PIONEER_GROUP_SIZE=8 COLONY_DISASTER_PROBABILITY=0.03 \
			COLONY_DISASTER_INSULT_MULTIPLIER=200000 TICKS_TO_QUEEN_REPLACEMENT=20 \
			MATES_WITHIN_COLONY=True LIN_RUN_NUM=10 REPRODUCTION_TYPE=Lin \
			DO_NOT_USE_AGE_DIST_EVER=True

runIndividualMode:
	$(MAKE) setRunDetails DATA_COMBINATIONS="R9r9.2HP1000hp900pcE1.1pcM0.07BP0.3TTF3IRT2 \
			R9r9.3HP1000hp600pcE1.1pcM0.07BP0.3TTF3IRT2 R9r10HP1000hp900pcE1.1pcM0.07BP0.3TTF3IRT2" \
			USE_COLONY_REPRODUCTION_MODE=False IMMORTAL_QUEEN_MODE=False \
			LIN_RUN_NUM=10 EXP_RUN_NUM=10 DO_NOT_USE_AGE_DIST_EVER=True
			
runAntMode:
	$(MAKE) setRunDetails DATA_COMBINATIONS="R9r10HP1000hp900pcE1.1pcM0.07BP0.3TTF3IRT2 \
	        R9r9.3HP1000hp600pcE1.1pcM0.07BP0.3TTF3IRT2 R9r9.2HP1000hp900pcE1.1pcM0.07BP0.3TTF3IRT2" \
			USE_COLONY_REPRODUCTION_MODE=True N0_COLONIES=16 THRESHOLD_COLONY_SIZE=23 \
			MATING_FLIGHT_SELECTION='youngest' QUEEN_HP_MULTIPLIER=1 \
			COLONY_LOGISTIC_MULTIPLIER=0.06 COLONY_LOGISTIC_EXP=1.1 \
			PIONEER_GROUP_SIZE=1 COLONY_DISASTER_PROBABILITY=0 \
			COLONY_DISASTER_INSULT_MULTIPLIER=1 TICKS_TO_QUEEN_REPLACEMENT=None \
			MATES_WITHIN_COLONY=False LIN_RUN_NUM=10 REPRODUCTION_TYPE=Lin \
			DO_NOT_USE_AGE_DIST_EVER=True

# Run Individual, Ant, and Mole-Rat modes with default configs (Lin+Exp, 10 runs each)
runDefault: runIndividualMode runAntMode runMoleRatMode



# This will get equilibrium population size, age ditribution and lifespan distributions 
GetInitialDistributionData:
	$(MAKE) setRunDetails LIN_RUN_NUM=10 EXP_RUN_NUM=10 LONG_RUN=True;

# This will get the age ditribution
GetAgeDistributionData:
	$(MAKE) setRunDetails LIN_RUN_NUM=10 EXP_RUN_NUM=10 GET_AGE_DIST=True;


cleanCopy:
	mkdir ../NMR_SIM_Copy
	mkdir ../NMR_SIM_Copy/config_files
	mkdir ../NMR_SIM_Copy/data
	mkdir ../NMR_SIM_Copy/csv
	mkdir ../NMR_SIM_Copy/graphs
	cp makefile ../NMR_SIM_Copy/
	cp *.py ../NMR_SIM_Copy/
	cp requirements.txt ../NMR_SIM_Copy/
	cp README.md ../NMR_SIM_Copy/

