library(tidyverse)
library(pracma)


# Better sampling function
resample = \(x, ...) x[sample.int(length(x), ...)]

# Mortality follows the Gompertz curve
Mortality = \(
  age, 
  N, 
  density_dependence, 
  mu0, # mortality at birth (age = 0) when N << K
  mup0 # derivative of mortality at birth when N << K
){
  
   beta = 1 / log(1 / (1 - mu0)) * mup0 / (1 - mu0)
   alpha = mup0 / (1 - mu0) * 1 / (exp(beta) - 1)
 
  return(
    1 - 
      exp(-density_dependence * N) * 
      exp(-alpha / beta * exp(beta * age) * (exp(beta) - 1))
  )
}

# Number of offspring per year is K in linear population, f0 * N in exponential population
Fecundity = \(model, N, K, f0){
  
   if(model == 'linear') fecundity = K
   if(model == 'exponential') fecundity = f0
   return(fecundity)  
} 

# Determine mortality phenotype based on alleles
Phenotype = \(L1, L2, dominance){
   phenotype = 
      (dominance == 0) * ifelse(any(L1 == 1, L2 == 1), 1, 2) +
      (dominance == 1) * ifelse(any(L1 == 2, L2 == 2), 2, 1) +
      (dominance == .5) * ifelse(all(L1 == 1, L2 == 1), 1, ifelse(all(L1 == 2, L2 == 2), 2, 3))
   
   return(phenotype)
}

# Call Exponential or Linear population model
Model = \(model, ...){
   if(model == 'exponential') return(Exponential(model, ...))
   if(model == 'linear') return(Linear(model, ...))
}

# Exponential Model
Exponential = \(
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
){
   
   set.seed(seed)
   
   c1 = c(m0_wt, m0_mutant, (m0_wt + m0_mutant) / 2)
   c2 = c(mp0_wt, mp0_mutant, (mp0_wt + mp0_mutant) / 2)
   
   density_dependence = 1 / carrying_capacity
   
   time = 0
   step = 0
   burnin_steps = floor(burnin / interval)
   
   # initial conditions
   population = 
      tibble(
         Time = time,
         age = 0,
         sex = rep(c('F', 'M'), N0 / 2),
         L1 = 1,
         L2 = 1
      ) |>
      rowwise() |>
      mutate(phenotype = Phenotype(L1, L2, dominance)) |>
      ungroup() |>
      mutate(
         mu0 = c1[phenotype],
         mup0 = c2[phenotype],
         mortality = 
           Mortality(
             age = 0,
             N = n(),
             mu0 = mu0, 
             mup0 = mup0,
             density_dependence = density_dependence
            ),
         fecundity = 
            ifelse(
               sex == 'F',
               Fecundity(
                  model = 'exponential', 
                  N = n(), 
                  K = carrying_capacity, 
                  f0 = intrinsic_fecundity
               ),
               0
            )
      )
   
   record = population
   
   while(time < maxtime){
      
      # update time
      time = time + interval
      step = step + 1
      
      # deaths
      death = rbinom(n = nrow(population), size = 1, p = interval * population$mortality)
      population = population[death == 0, ]
      
      # stop if population collapses
      if(nrow(population) == 0) break
      
      # stop if mutant is lost or fixes
      if(all(step > burnin_steps, p0 > 0, population$L1 == 1, population$L2 == 1)) break
      if(all(step > burnin_steps, p0 < 1, population$L1 == 2, population$L2 == 2)) break
      
      # update fecundity due to change in population size,
      # update mortality due to aging
      population =
         population |>
         mutate(
            Time = time,
            age = interval + age
         ) |>
         rowwise() |>
         mutate(phenotype = Phenotype(L1, L2, dominance)) |>
         ungroup() |>
         mutate(
            mu0 = c1[phenotype],
            mup0 = c2[phenotype],
            mortality = 
              Mortality(
                age = age,
                N = n(),
                mu0 = mu0, 
                mup0 = mup0,
                density_dependence = density_dependence
              ),
            fecundity = 
               ifelse(
                  sex == 'F',
                  Fecundity(
                     model = 'exponential', 
                     N = n(), 
                     K = carrying_capacity, 
                     f0 = intrinsic_fecundity
                  ),
                  0
               )
         )
      
      if(all(population$sex == 'M')) next
      if(all(population$sex != 'M')) next
      
      # births
      males = 
         population |>
         filter(sex == 'M')
      
      females = 
         population |>
         filter(sex == 'F') |>
         mutate(
            mate = resample(seq(nrow(males)), replace = TRUE, size = n()),
            mate_L1 = males$L1[mate],
            mate_L2 = males$L2[mate],
            number_of_offspring = rpois(n = n(), lambda = interval * fecundity)
         ) 
      
      offspring = 
         females |>
         rowwise() |>
         reframe(
            Time = time, 
            age = 0,
            sex = resample(c('F', 'M'), replace = TRUE, size = number_of_offspring),
            L1 = resample(c(L1, L2), replace = TRUE, size = number_of_offspring),
            L2 = resample(c(mate_L1, mate_L2), replace = TRUE, size = number_of_offspring),
            phenotype = Phenotype(L1, L2, dominance)
         ) |>
         ungroup() |>
         mutate(
           mu0 = c1[phenotype],
           mup0 = c2[phenotype],
           mortality = 
             Mortality(
               age = 0,
               N = n() + nrow(population),
               mu0 = mu0, 
               mup0 = mup0,
               density_dependence = density_dependence
             ),
           fecundity = 
             ifelse(
               sex == 'F',
               Fecundity(
                 model = 'exponential', 
                 N = n() + nrow(population), 
                 K = carrying_capacity, 
                 f0 = intrinsic_fecundity
               ),
               0
             )
         )
      
      # add offspring to the population
      population =
         population |>
         bind_rows(offspring)
      
      # insert mutants at proportion p0 when time == burnin
      if(step == burnin_steps){
        index = resample(1:nrow(population), round(p0 * nrow(population)))
        
        mutants = 
          population[index, ] |>
          mutate(
            age = 0,
            L1 = 2,
            L2 = 2
          ) |>
          rowwise() |>
          mutate(phenotype = Phenotype(L1, L2, dominance)) |>
          ungroup() |>
          mutate(
            mu0 = c1[phenotype],
            mup0 = c2[phenotype],
            mortality = 
              Mortality(
                age = 0,
                N = n(),
                mu0 = mu0, 
                mup0 = mup0,
                density_dependence = density_dependence
              )
          )
        
        population[index, ] = mutants
      }
      
      # record outcome
      if(save == TRUE & time %% 10 == 0){
         record =
            record |>
            bind_rows(population)
      }
   }
   
   parms = 
      tibble(
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
   
   if(save){
      return(
         parms |> 
            bind_cols(
               record |> 
                  bind_rows(population) |>
                  rename(time = Time)
            )
      )
   }
   
   if(!save){
      return(
         parms |> 
            bind_cols(
               population |> 
                  rename(time = Time)
            )
      )
   }
      
   
}


# Linear Model
Linear = \(
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
){
   set.seed(seed)
   
   c1 = c(m0_wt, m0_mutant, (m0_wt + m0_mutant) / 2)
   c2 = c(mp0_wt, mp0_mutant, (mp0_wt + mp0_mutant) / 2)
   
   density_dependence = 1 / carrying_capacity
   
   time = 0
   step = 0
   burnin_steps = floor(burnin / interval)
   
   population = 
      tibble(
         Time = time,
         age = 0,
         sex = rep(c('F', 'M'), N0 / 2),
         L1 = 1,
         L2 = 1
      ) |>
      rowwise() |>
      mutate(phenotype = Phenotype(L1, L2, dominance)) |>
      ungroup() |>
      mutate(
         mu0 = c1[phenotype],
         mup0 = c2[phenotype],
         mortality = 
           Mortality(
             age = age,
             N = n(),
             mu0 = mu0, 
             mup0 = mup0,
             density_dependence = density_dependence
           ),
         fecundity = 0
      )
   
   queen_ID = with(population, resample(which(sex == 'F'), size = 1))
   queen_phenotype = population$phenotype[queen_ID]
   population$sex[queen_ID] = 'Q'
   population$fecundity[queen_ID] = 
      Fecundity(
         model = 'linear', 
         N = nrow(population), 
         K = carrying_capacity, 
         f0 = intrinsic_fecundity
      )
   
   record = population
   
   while(time < maxtime){
      
      # update time
      time = time + interval
      step = step + 1
      
      # deaths
      death = rbinom(n = nrow(population), size = 1, p = interval * population$mortality)
      population = population[death == 0, ]
      
      # stop if population collapses
      if(nrow(population) == 0) break
      
      # stop if mutant is lost or fixes
      if(all(step > burnin_steps, p0 > 0, population$L1 == 1, population$L2 == 1)) break
      if(all(step > burnin_steps, p0 < 1, population$L1 == 2, population$L2 == 2)) break
      
      # update mortality due to aging
      population =
         population |>
         mutate(
            Time = time,
            age = interval + age
         ) |>
         rowwise() |>
         mutate(phenotype = Phenotype(L1, L2, dominance)) |>
         ungroup() |>
         mutate(
            mu0 = c1[phenotype],
            mup0 = c2[phenotype],
            mortality = 
              Mortality(
                age = age,
                N = nrow(population),
                mu0 = mu0, 
                mup0 = mup0,
                density_dependence = density_dependence
              )
         )
      
      # don't bother with reproduction if only males are left or if all males are gone
      if(all(population$sex == 'M')) next
      if(all(population$sex != 'M')) next
      
      # has the queen died?
      if(all(population$sex != 'Q')){
         queen_ID = with(population, resample(which(sex == 'F'), size = 1))
         population$sex[queen_ID] = 'Q'
      }
      
      # update queen fecundity changes due to change in population size
      population = 
         population |>
         mutate(
            fecundity = 
               ifelse(
                  sex == 'Q',
                  Fecundity(
                    model = 'linear', 
                    N = nrow(population), 
                    K = carrying_capacity, 
                    f0 = intrinsic_fecundity
                  ),
                  0
               )
         )
      
      # births
      males = 
         population |>
         filter(sex == 'M')
      
      queen = 
         population |>
         filter(sex == 'Q') |>
         mutate(
            mate = resample(seq(nrow(males)), size = 1),
            mate_L1 = males$L1[mate],
            mate_L2 = males$L2[mate],
            number_of_offspring = rpois(n = 1, lambda = interval * fecundity)
         ) 
      
      offspring = 
         queen |>
         rowwise() |>
         reframe(
            Time = time, 
            age = 0,
            sex = resample(c('F', 'M'), replace = TRUE, size = number_of_offspring),
            L1 = resample(c(L1, L2), replace = TRUE, size = number_of_offspring),
            L2 = resample(c(mate_L1, mate_L2), replace = TRUE, size = number_of_offspring)
         ) |>
         mutate(phenotype = Phenotype(L1, L2, dominance)) |>
         ungroup() |>
         mutate(
            mu0 = c1[phenotype],
            mup0 = c2[phenotype],
            mortality = 
              Mortality(
                age = age,
                N = n() + nrow(population),
                mu0 = mu0, 
                mup0 = mup0,
                density_dependence = density_dependence
              ),
            fecundity = 0
         )
         
      # add offspring to the population
      population =
         population |>
         bind_rows(offspring)
      
      # insert mutants at proportion p0 when time == burnin
      if(step == burnin_steps){
        index = resample(1:nrow(population), round(p0 * nrow(population)))
        
        mutants = 
          population[index, ] |> 
          mutate(
            age = 0,
            L1 = 2,
            L2 = 2
          ) |>
          rowwise() |>
          mutate(phenotype = Phenotype(L1, L2, dominance)) |>
          ungroup() |>
          mutate(
            mu0 = c1[phenotype],
            mup0 = c2[phenotype],
            mortality = 
              Mortality(
                age = 0,
                N = n(),
                mu0 = mu0, 
                mup0 = mup0,
                density_dependence = density_dependence
              )
          )
        
        population[index, ] = mutants
      }
      
      # record outcome
      if(save == TRUE & time %% 10 == 0){
         record =
            record |>
            bind_rows(population)
      }
   }
   
   parms = 
      tibble(
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
   
   if(save){
      return(
         parms |> 
            bind_cols(
               record |> 
                  bind_rows(population) |>
                  rename(time = Time)
            )
      )
   }
   
   if(!save){
      return(
         parms |> 
            bind_cols(
               population |> 
                  rename(time = Time)
            )
      )
   }

}


# Execute simulation runs specified in parameters
Simulate = 
   \(
      parameters, 
      save.data, 
      dir.path
   ){
      result = 
         parameters |>
         future_pmap_dfr(
            .f = Model,
            .options = furrr_options(seed = NULL)
         )
      
      if(save.data){
         filename = paste0(dir.path, 'results.RDS')
         saveRDS(result, filename)
      } 
      
      if(!save.data) return(result)
   }

