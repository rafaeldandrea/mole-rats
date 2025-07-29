

#### Makefile ####
RUN_MULTIPLIER ?= 30
LONG_RUN ?= False
POP_A_FRAC ?= False
SHOW_INDIV_RUN_STATS ?= False
SHOW_INDIV_RUN_GRAPH ?= False
GET_AGE_DIST ?= False
SET_INIT_POP_TO_EQUILIBRIUM = True
INIT_REC_ALLELE_FRACTION =0.5 #.75 means that 75% of alleles will be recessive in the beginning, .285 means 28.8% and so on.
INIT_POP = 104 # Only if SET_INIT_POP_TO_EQUILIBRIUM = False will this be used
DATA_COMBINATIONS ?=  wR9r10HP1000hp900pcE1.1pcM0.07BP0.3TTF3IRT2 wR10r9HP900hp1000pcE1.1pcM0.07BP0.3TTF3IRT2 # #wR9r10HP1000hp900pcE1.1pcM0.07BP0.3TTF3 #LS4wR10r9HP900hp1000 #HP290pcE2IdT3BP.9 #HP1200hp1000 HP1500hp1000 # R.72r.7 R.73r.7 R.71r.7 R5.5r5 R6r5 R7r5 HP1100hp1000 #WA-T7-HP1000hp1200R6r5 #T7-R6r5 #B2-R6r5 #B2-R6r5 T6-HP1200hp1000 #T6-HP1200hp1000R6r5 T6-R6r5  #T7-R6r5 #T5-R6r5 T6-R6r5 # B2-R6r5 #T4-R6r5 T3-R6r5  
MAX_TICKS = 5001

# REPRODUCTION_TYPE=Lin
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
    # UseTradeoffAllele = False (w)
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
    # UseTradeoffAllele = False
    # insult_recovery_ticks = 2 (IRT)

# Ex. wR10r9HP900hp1000

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
    # UseTradeoffAllele = True
    # insult_recovery_ticks = 0 (IRT)

USE_AVG_POP_SIZE_FOR_LOG_CAP = False
# REPRODUCTION_TYPE ?=BasicLin
LIN_DATA_COMBINATIONS = test6 #williamsAlleleT7 #T9 #T6 #T8B T8 #T8B T9 #T6 T3 A2 logistic abs_logistic T8 # T5 T3 #B2 T4  T2 #A1 A2 A3 #B1 B2 B3 B4 C1 C2 C3 T1 


# Default target
default: runSim 
	# GetInitialDistributionData 
	# make GetPopAFraction 

run:
	$(foreach var,$(DATA_COMBINATIONS),\
		python3 nmr_simulator.py $(var)$(REPRODUCTION_TYPE) \
			NumRuns $(NUM_RUNS) \
			LongRun $(LONG_RUN)\
			GetPopAFraction $(POP_A_FRAC)\
			GetAgeDistribution $(GET_AGE_DIST)\
			ShowIndividualRunStats $(SHOW_INDIV_RUN_STATS)\
			ShowIndividualRunGraphs $(SHOW_INDIV_RUN_GRAPH)\
			UseAvgForLogCap $(USE_AVG_POP_SIZE_FOR_LOG_CAP)\
			UseEquilibriumForInitPop $(SET_INIT_POP_TO_EQUILIBRIUM)\
			InitialPopulation $(INIT_POP)\
			MaxTicks $(MAX_TICKS)\
			InitialRecessiveAlleleFraction $(INIT_REC_ALLELE_FRACTION);\
		)

setRunDetails:

	@if [ "$(REPRODUCTION_TYPE)" = "" ]; then \
		$(MAKE) run REPRODUCTION_TYPE=Lin NUM_RUNS=`echo $$(($(LIN_RUN_NUM) * $(RUN_MULTIPLIER)))`; \
		$(MAKE) run REPRODUCTION_TYPE=Exp NUM_RUNS=`echo $$(($(EXP_RUN_NUM) * $(RUN_MULTIPLIER)))`; \
	elif [ "$(REPRODUCTION_TYPE)" = "BasicLin" ]; then \
		$(MAKE) run REPRODUCTION_TYPE='' NUM_RUNS=`echo $$(($(LIN_RUN_NUM) * $(RUN_MULTIPLIER)))` DATA_COMBINATIONS="$(LIN_DATA_COMBINATIONS)"; \
	elif [ "$(REPRODUCTION_TYPE)" = "Lin" ]; then \
		$(MAKE) run NUM_RUNS=`echo $$(($(LIN_RUN_NUM) * $(RUN_MULTIPLIER)))`; \
	elif [ "$(REPRODUCTION_TYPE)" = "Exp" ]; then \
		$(MAKE) run NUM_RUNS=`echo $$(($(EXP_RUN_NUM) * $(RUN_MULTIPLIER)))`; \
	fi ;\

# This will get lifespan distribution
runSim:
	$(MAKE) setRunDetails LIN_RUN_NUM=15 EXP_RUN_NUM=15

# This will get equilibrium population size, age ditribution and lifespan distributions 
GetInitialDistributionData:
	$(MAKE) setRunDetails LIN_RUN_NUM=10 EXP_RUN_NUM=5 LONG_RUN=True;

# This will get the age ditribution
GetAgeDistributionData:
	$(MAKE) setRunDetails LIN_RUN_NUM=10 EXP_RUN_NUM=5 GET_AGE_DIST=True;

# This gets the initial recessive allele fraction 
GetPopAFraction:
	$(MAKE) setRunDetails LIN_RUN_NUM=100 EXP_RUN_NUM=10 POP_A_FRAC=True;

GetInitialPopAFraction:
	$(MAKE) setRunDetails LIN_RUN_NUM=100 EXP_RUN_NUM=10 POP_A_FRAC="True InitialRecessiveAlleleFraction .5 GetLifespanDistribution False";

changeNumRuns:
	find ./config_files -type f -exec sed -i 's/NumRuns = 50,/NumRuns = 5000,/g' {} \;

removeAllGeneratedResults%:
	find ./config_files/$*.py -type f -exec sed -i "/######## Generated Responses ########/q" {} \;
	# rm data/*;
	# rm csv/*;
	# rm graphs/*;
	# mkdir data;
	# mkdir csv;
	# mkdir graphs;

POSSIBLE_pop_cap_exponent ?= 1.9  3 #1.2 1.1 #1.3 1.6 1 #.9 1.1 1.2 
POSSIBLE_pop_cap_mult ?= 1 2 #.01 .05 .1 .5 1 .001 


runForLogCap1:
	$(foreach var,$(POSSIBLE_pop_cap_exponent),\
		make runForLogCap POP_CAP_EXPONENT=$(var);\
		)

runForLogCap:
	$(foreach var,$(POSSIBLE_pop_cap_mult),\
		python3 nmr_simulator.py T7-R6r5Lin \
			NumRuns 100 \
			LongRun True\
			pop_cap_exponent $(POP_CAP_EXPONENT)\
			pop_cap_mult $(var)\
			ShowIndividualRunStats $(SHOW_INDIV_RUN_STATS);\
		python3 nmr_simulator.py T7-R6r5Exp \
			NumRuns 50 \
			LongRun True\
			pop_cap_exponent $(POP_CAP_EXPONENT)\
			pop_cap_mult $(var)\
			ShowIndividualRunStats $(SHOW_INDIV_RUN_STATS);\
		)


# default: GetInitialDistributionData
# default: GetInitialDistributionData GetInitialPopAFraction runSim
# default:  GetAgeDistributionData runSim
# default:  runSim 
# default:  GetInitialPopAFraction
# virtualenv venv                    
# source venv/bin/activate