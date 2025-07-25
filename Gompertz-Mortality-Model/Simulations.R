# Required libraries: tidyverse, furrr, pracma
library(tidyverse)
library(furrr)

# Local parameters
test = TRUE   # Change to TRUE for the full set of simulations
save_address = "[LOCAL ADDRESS]"  # Change to appropriate address for saving simulation data

# Load project functions
source("https://github.com/rafaeldandrea/mole-rats/raw/refs/heads/main/Gompertz-Mortality-Model/Functions.R") 

# Directory to save data
save_directory = 
  paste0(
    save_address, 
    Sys.Date(),
    '/'
  )

# Create directory to save data
dir.create(save_directory, showWarnings = FALSE)

# Parallel jobs
plan(multicore, workers = availableCores() - 2)

# Fixed parameters
parameters = 
  tibble(
    mutant_type = c('longevity', 'vitality'),
    carrying_capacity = 1e2,
    interval = 3 / carrying_capacity,
    intrinsic_fecundity = 2 / (1 - interval),
    N0 = carrying_capacity,
    burnin = 100,
    maxtime = 1100,
    save = FALSE, # if TRUE, saves transient states
    m0_wt = .1,
    m0_mutant = c(.1, .05),
    mp0_wt = .02,
    mp0_mutant = c(.015, .02)
  ) |>
  expand_grid(
    model = c('linear', 'exponential'),
    p0 = seq(0, 20, by = 2) / 200,
    dominance = 1,
    seed = 1:5000
  ) |>
  filter((model == 'exponential') | (seed <= 500)) |>
  rowid_to_column(var = 'ID') |>
  select(
    ID,
    model,
    N0,
    p0,
    mutant_type,
    dominance,
    carrying_capacity,
    interval,
    intrinsic_fecundity,
    m0_wt,
    m0_mutant,
    mp0_wt,
    mp0_mutant,
    burnin,
    maxtime,
    seed,
    save
  )

# If running on personal computer, test code on a single run
if(test == TRUE){
  parameters =
    parameters |>
    filter(
      model == 'linear',
      mutant_type == 'longevity',
      seed <= 2,
      p0 == .1
    )
}


# Execute code
out = 
   parameters |>
   Simulate(
      save.data = TRUE,
      dir.path = save_directory
   )
