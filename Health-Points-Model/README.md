# NMR Simulator

A  simulation of population dynamics, modeling allele competition between longevity (HP-R) and vitality (hp-r) genotypes. The simulator supports both individual (single population) and colony-based reproduction modes, with configurable Gompertz aging, insult distributions, and population growth models.

## Table of Contents

- [Requirements](#requirements)
- [Installation](#installation)
- [Project Structure](#project-structure)
- [Quick Start](#quick-start)
- [Running the Simulation](#running-the-simulation)
- [Makefile Reference](#makefile-reference)
- [Configuration](#configuration)
- [Data Combination Codes](#data-combination-codes)
- [Output](#output)
- [Modes of Operation](#modes-of-operation)

---

## Requirements

- **Python 3.x** (3.9+ recommended)
- **NumPy**
- **SciPy**
- **Matplotlib**
- **pandas**

---

## Installation

1. Clone or download this repository.

2. Create and activate a virtual environment (recommended):

   ```bash
   python3 -m venv venv
   source venv/bin/activate   # On macOS/Linux
   # or: venv\Scripts\activate   # On Windows
   ```

3. Install dependencies:

   ```bash
   pip install numpy scipy matplotlib pandas
   ```

4. Create required directories (they will be auto-created if missing, but you may need them for output):

   ```bash
   mkdir -p config_files csv graphs data
   ```

---

## Project Structure

```
NMR_Sim/
├── nmr_simulator.py   # Main entry point
├── makefile           # Build and run targets
├── get.py             # Configuration parsing, parent selection, insult severity
├── update.py          # Population/colony updates, births, mating flights
├── display.py         # Statistics, graphs, CSV output
├── gen.py             # Base CONFIG template
├── Naked_Mole_Rat.py  # NMR agent with aging and insult mechanics
├── colony.py          # Colony class (queen + members)
├── config_files/      # Auto-generated per-run configs
├── data/              # Fixation and aggregate data
├── csv/               # Per-tick and per-run CSV exports
└── graphs/            # Generated plots
```

---

## Quick Start

**Fast test run (verify setup):**

```bash
make quick
```

**Default full simulation (Individual + Ant + Mole-Rat modes as seen in the research paper, 10 runs each):**

```bash
make
```

**Run with custom configuration:**

```bash
make runSim DATA_COMBINATIONS=R9r10HP1000hp900pcE1.1pcM0.07BP0.3TTF3IRT2 REPRODUCTION_TYPE=Exp NUM_RUNS=20
```

---

## Running the Simulation

### Using the Makefile (Recommended)

| Command | Description |
|--------|-------------|
| `make` | Default: runs Individual, Ant, and Mole-Rat modes (R9r9.2, R9r9.3, R9r10; 10 runs each) |
| `make runDefault` | Same as `make`: Individual + Ant + Mole-Rat modes |
| `make quick` | Fast test (2 Individual runs, Lin only) — use to verify setup |
| `make run` | Runs as according to makefile setting or additions (such as `NUM_RUNS=4`) |
| `make runIndividualMode` | Individual mode: R9r9.2, R9r9.3, R9r10 (Lin + Exp, 10 runs each as shown in the paper) |
| `make runAntMode` | Ant colony mode: stored sperm, no disasters (Lin, 10 runs as shown in the paper) |
| `make runMoleRatMode` | Mole-rat colony mode: in-colony mating, disasters, queen replacement (Lin, 10 runs as shown in the paper) |
| `make GetInitialDistributionData` | Equilibrium + age + lifespan (sets `LONG_RUN=True`) |
| `make GetAgeDistributionData` | Age distribution only |
| `make cleanCopy` | Copy project to `../NMR_SIM_Copy` without any generated data |
| `make help` | List all targets and options |

### Running Directly with Python

The main script is `nmr_simulator.py`. It expects at minimum a **data combination string** and key-value overrides:

```bash
python3 nmr_simulator.py <DATA_COMBINATION><REPRODUCTION_TYPE> \
  NumRuns 50 \
  UseColonyReproductionMode False \
  WriteDataToCSV True \
  MaxTicks 100000
```

**Required arguments:**

- First positional argument: `DATA_COMBINATION` + `REPRODUCTION_TYPE` (e.g., `R9r9.2HP1000hp900pcE1.1pcM0.07BP0.3TTF3IRT2Lin`).

**Common overrides:**

- `NumRuns 50` — number of simulation runs
- `UseColonyReproductionMode True|False`
- `MaxTicks 1000000`
- `WriteDataToCSV True|False`
- `LongRun True` — enables age distribution, lifespan distribution, mortality, equilibrium
- `GetAgeDistribution True`
- `ShowIndividualRunStats True`
- `ShowIndividualRunGraphs True`

### Config File Generation

If no config file exists for the given combination, one is created automatically under `config_files/` from the data combination string. The simulator will then load and use it.

---

## Makefile Reference

### Main Variables

| Variable | Default | Description |
|----------|---------|-------------|
| `DATA_COMBINATIONS` | `R9r9.2HP1000hp900pcE1.1pcM0.07BP0.3TTF3IRT2` | Config code(s) for run(s) |
| `REPRODUCTION_TYPE` | (both) | `Lin`, `Exp`, or empty for both |
| `NUM_RUNS` | derived from `LIN_RUN_NUM`/`EXP_RUN_NUM` × `RUN_MULTIPLIER` | Total runs per config |
| `LIN_RUN_NUM` | (set per target) | Number of Linear runs |
| `EXP_RUN_NUM` | (set per target) | Number of Exponential runs |
| `RUN_MULTIPLIER` | 1 | Multiplier for Lin/Exp run counts if using target `setRunDetails` |
| `LONG_RUN` | False | Enable the collection and storing of equilibrium, lifespan, mortality, etc. (Warning: will slow runtime)|
| `MAX_TICKS` | 1000000 | Max simulation ticks per run |
| `WRITE_DATA_TO_CSV` | True | Write per-tick/per-run CSVs |
| `USE_OVERALL_RUN_STATS` | False | Use overall run statistics |
| `MAKE_DISTRIBUTION_GRAPHS` | True | Makes graphs for age distribution, lifespan distribution and mortality and saves them to graphs directiory (Warning: this will only occur durring long runs or runs that generate age distribution) |
| `USE_AVG_POP_SIZE_FOR_LOG_CAP` | False | Use average population size for logistic cap (Warning: to use this feature you must have previously generated an average population size by running the simulation with the current configuration using LONG_RUN=True) |
| `SHOW_INDIV_RUN_STATS` | False | Show statistics for each individual run |
| `SHOW_INDIV_RUN_GRAPH` | False | Show graphs for each individual run |
| `INIT_REC_ALLELE_FRACTION` | 0.5 | Initial fraction of recessive alleles on simulation comencement (e.g. 0.75 = 75%) |
| `SET_INIT_POP_TO_EQUILIBRIUM` | False | Start from equilibrium population size (Warning: to use this feature you must have previously generated an average population size by running the simulation with the current configuration using LONG_RUN=True) |
| `INIT_POP` | 16 | Initial population (used if `SET_INIT_POP_TO_EQUILIBRIUM=False` or no equilibrium population is available) |
| `IMMORTAL_QUEEN_MODE` | False | Queen does not age or die from insults (only applicable in colony simulation)|
| `GET_AGE_DIST` | False | Get age distribution data with run (Warning: will slow runtime) |
| `DO_NOT_USE_AGE_DIST_EVER` | True | Disable use of age distribution from config |

### Colony Mode Variables (only in proper use when `USE_COLONY_REPRODUCTION_MODE=True` and `REPRODUCTION_TYPE=Lin`)

| Variable | Default | Description |
|----------|---------|-------------|
| `N0_COLONIES` | 16 | Initial number of colonies |
| `THRESHOLD_COLONY_SIZE` | 23 | Colony size threshold for participation in mating flights |
| `MATING_FLIGHT_INTERVAL` | 6 | Ticks between mating flights |
| `MATING_FLIGHT_SELECTION` | `youngest` | `random`, `youngest`, or `oldest` |
| `QUEEN_HP_MULTIPLIER` | 1 | Multiplier for queen health points (to reflect difference in lifestyle) |
| `BASE_COLONY_FORMATION_PROBABILITY` | 1 | Base probability of colony formation (Should be a number between 0 and 1)|
| `COLONY_LOGISTIC_MULTIPLIER` | 0.06 | Logistic growth multiplier for colony formation probability (to simulate scarcity of suitable colony establishment locations) |
| `COLONY_LOGISTIC_EXP` | 1.1 | Logistic growth exponent for colony formation probability (to simulate scarcity of suitable colony establishment locations) |
| `PIONEER_GROUP_SIZE` | 1 | Members of pioneer group necesary to establish colony (1 for stored-sperm or an even number for colonies that mate internally) |
| `COLONY_DISASTER_PROBABILITY` | 0.1 | Probability of disaster per tick per colony |
| `COLONY_DISASTER_INSULT_MULTIPLIER` | 200000 | Insult severity multiplier on disaster tick |
| `TICKS_TO_QUEEN_REPLACEMENT` | 20 | Replace queen after N ticks; `None` = no replacement |
| `MATES_WITHIN_COLONY` | False | Allow in-colony mating (vs stored sperm from mating pool) |

---

## Configuration

### Data Combination Codes

The first argument encodes simulation parameters. Format: `R{r1}r{r2}HP{hp1}hp{hp2}pcE{exp}pcM{mult}BP{bp}TTF{ttf}IRT{irt}`.

| Code | Parameter | Default | Description |
|------|-----------|---------|-------------|
| `R` | R | 10 | Half-life ticks for longevity genotype (W) |
| `r` | r | 10 | Half-life ticks for vitality genotype (w) |
| `HP` | HP | 900 | Initial health points for longevity (W) |
| `hp` | hp | 900 | Initial health points for vitality (w) |
| `pcE` | pop_cap_exponent | 1.2 | Population cap exponent |
| `pcM` | pop_cap_mult | 0.1 | Population cap multiplier |
| `BP` | BirthProbability | 0.4 | Birth probability per tick |
| `TTF` | TicksTillFertile | 1 | Ticks until fertile |
| `IRT` | insult_recovery_ticks | 0 | Insult recovery time |

**Examples:**

- `R9r10HP1000hp900pcE1.1pcM0.07BP0.3TTF3IRT2` — R=9, r=10, HP=1000, hp=900, etc.
- `R9r9.2HP1000hp900pcE1.1pcM0.07BP0.3TTF3IRT2` — Use `.` for decimals (e.g., 9.2 → `9pt2` in filename).

### Reproduction Types

- **Lin** — Linear population growth (single queen, all births from queen).
- **Exp** — Exponential growth (each female will reproduce after `TicksTillFertile`).

---

## Output

- **Terminal:** Fixation counts, average lifespans, equilibrium stats, colony stats (if colony mode).
- **`data/`:** `fixation-data.csv` and other aggregate outputs.
- **`csv/`:** Per-tick and per-run CSVs when `WriteDataToCSV` is True.
- **`graphs/`:** Plots ( lifespan distributions, age distributions, etc.).
- **`config_files/`:** Generated configs with equilibrium, lifespan, and age data appended for reuse.

---

## Modes of Operation

### Individual Mode (`UseColonyReproductionMode False`)

- Single colony population.
- Linear: one queen; all births from queen (subject to fertility and population cap).
- Exponential: each female reproduces independently (subject to fertility and population cap).

### Colony Mode (`UseColonyReproductionMode True`)

- Multiple colonies, each with a queen and workers.
- Periodic mating flights; selected individuals form new colonies with some probability.
- Colony disasters, queen replacement, and in-colony vs. stored-sperm mating are configurable.
- **Ant mode:** `MATES_WITHIN_COLONY False` (stored sperm).
- **Mole-rat mode:** `MATES_WITHIN_COLONY True` (in-colony mates).

---

## Example Commands

```bash
# Quick sanity check
make quick

# Default full run
make

# Custom data combo and run count
make runSim DATA_COMBINATIONS=R9r10HP1000hp900pcE1.1pcM0.07BP0.3TTF3IRT2 NUM_RUNS=50

# Lin only, 20 runs
make run REPRODUCTION_TYPE=Lin NUM_RUNS=20

# Long run (equilibrium, age, lifespan)
make GetInitialDistributionData

# Age distribution only
make GetAgeDistributionData

# Copy project for backup
make cleanCopy
```

