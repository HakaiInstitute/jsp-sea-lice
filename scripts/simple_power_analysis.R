library(tidyverse)
library(here)
library(glmmTMB)
library(parallel)
library(foreach)
library(doParallel)

num_cores <- detectCores() - 1
registerDoParallel(cores = num_cores)


sealice_abund_long <- read_csv("https://raw.githubusercontent.com/HakaiInstitute/jsp-data/develop/supplemental_materials/tidy_data/sealice_all_stages_ts_long.csv")

main_lice <- sealice_abund_long |> 
  mutate(spp_stage = paste(louse_species, life_stage, sep = "_")) |> 
  pivot_wider(id_cols = c("year", "survey_date", "seine_id", "ufn", "species"), 
              names_from = spp_stage, values_from = count, values_fill = 0) |> 
  dplyr::select(year:species, cal_att = "Caligus clemensi_attached", 
                cal_mot = "Caligus clemensi_motile", 
                lep_att = "Lepeophtheirus salmonis_attached", 
                lep_mot = "Lepeophtheirus salmonis_motile") |> 
  filter(year != 2022) |> 
  # Removed 2022 because it was the first year farms were removed and I don't want to simulate data based on that.
  dplyr::select(year, species, cal_mot, cal_att, lep_mot, lep_att, seine_id)


# Create simulated dataset that have cal_mot counts that have been reduced 
sim_lice <- main_lice |> 
  select(-c(lep_mot, lep_att, cal_att)) |> 
  mutate(farm = 1,
         year = year + 7,
         seine_id = as_factor(paste0(year, seine_id)),
         #NOTE: modify the below probs to adjust the effect size !
         rand_int = sample(c(0,1), size = nrow(main_lice), replace = TRUE, prob = c(0.5, 0.5)),
         #cal_mot = if_else(cal_mot > 0, cal_mot -1, cal_mot)
         #cal_mot = if_else(cal_mot > 3, cal_mot -1, cal_mot) #-2.7222
         #cal_mot = if_else(cal_mot > 3, cal_mot -2, cal_mot) #-5.4%
         #cal_mot = if_else(cal_mot > 2, cal_mot -1, cal_mot) #-7.3%
         #cal_mot = if_else(cal_mot > 2, cal_mot -2, cal_mot) #-15%
         #cal_mot = if_else(cal_mot > 1, cal_mot -1, cal_mot) #-21.4%
         #cal_mot = if_else(cal_mot > 0, cal_mot - sample(c(1,0), 1, prob = c(0.3, 0.7)), cal_mot)
         # cal_mot = case_when(cal_mot == 2 ~ 1, #-44.4
         #                      cal_mot == 3 ~ 2, #-58.5
         #                      TRUE ~ cal_mot)
         cal_mot = case_when(cal_mot == 1 ~ rand_int,
                           cal_mot >= 2 ~ cal_mot - rand_int, #-44.4
                            TRUE ~ cal_mot)
         )

effect <- (mean(sim_lice$cal_mot) - mean(main_lice$cal_mot)) / mean(main_lice$cal_mot) * 100
effect

main_lice <- main_lice |> 
  mutate(farm = 0)

n_iterations <- 10
sample_size <- 0.5


p_values <- foreach(i = 1:n_iterations, .combine = rbind) %dopar% {
  
  sim_lice_adj <- sim_lice |> 
    dplyr::group_by(year, seine_id, species) |> 
    slice_sample(prop = sample_size, replace = TRUE) |> 
    ungroup()
 
  sim_main_lice <- sim_lice_adj |> 
    bind_rows(main_lice) |> 
    mutate(farm = as_factor(farm),
           year = as_factor(year))

  with_farm <- glmmTMB(cal_mot ~ species + farm + (1 | year/seine_id),
                           ziformula = ~1,
                           data = sim_main_lice, family = nbinom2, 
                           na.action = na.omit)
  
  no_farm <- glmmTMB(cal_mot ~ species + (1 | year/seine_id),
                         ziformula = ~1,
                         data = sim_main_lice, family = nbinom2, 
                         na.action = na.omit)
  
  anova_comp <- anova(with_farm, no_farm)
  
  p_value <- anova_comp[2,"Pr(>Chisq)"]
  mean_change <- (mean(sim_lice_adj$cal_mot) - mean(main_lice$cal_mot)) / mean(main_lice$cal_mot) * 100
  
  tibble(p_value, mean_change)
}

p_values

mean(p_values$p_value, na.rm = TRUE)
mean(p_values$mean_change, na.rm = TRUE)
sum(p_values$p_value < 0.05) / (n_iterations - sum(is.na(p_values))) 




# Sample sizes that work:
# -66% cal_mot = if_else(cal_mot > 0, cal_mot -1, cal_mot)
# -21% cal_mot = if_else(cal_mot > 1, cal_mot -1, cal_mot)
#-15% cal_mot = if_else(cal_mot > 2, cal_mot -2, cal_mot) 
