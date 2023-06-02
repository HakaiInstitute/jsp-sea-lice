# Sealice power analysis
library(tidyverse)
library(here)
library(rsample)
library(parallel)
library(foreach)
library(doParallel)
library(glmmTMB)
library(beepr)
library(lmtest)

# Set up a parallel backend
num_cores <- detectCores() - 1
registerDoParallel(cores = num_cores)

sealice_abund_long <- read_csv("https://raw.githubusercontent.com/HakaiInstitute/jsp-data/develop/supplemental_materials/tidy_data/sealice_all_stages_ts_long.csv")

# wrangle data so every row is one fish
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


# Create resample function grouped by hierarchy of year, seine_id, and species
# so that sealice counts are properly nested wrt exp. design
re_sample <- function(df, sample_size) {
  df |> dplyr::group_by(year, seine_id, species) |>
    slice_sample(prop = sample_size, replace = TRUE) |> 
    ungroup()
}

# set effect size
sequence <- c(1, 0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1, 0)

# Set the proportion of the original sample size (of each speices in each seine) to simulate
sample_sizes <- c(0.6, 0.8, 1, 1.2, 1.4)
# create function to simulate data that uses effect sizes and sample sizes 
# and the column of sea lice species and stage you'd like to simulate 
simulate_data_v3 <- function(sample_sizes, column) {
  main_lice <- dplyr::select(main_lice, year, species, seine_id, column)
  simulated_datasets <- list()
  
  for (i in sequence) {
   #p_changes <- list()
    for (sample_size in sample_sizes) {
      
      #re sample sim_lice dataset (2022-2028) at prescribed sample size
      simulated_data <- re_sample(df = main_lice, sample_size)
      
      simulated_data$rand_int <- sample(c(0,1), size = nrow(simulated_data), 
                                                   replace = TRUE, 
                                                   prob = c(i, 1-i))
      
      simulated_data[[column]] <- case_when(simulated_data[[column]] == 1 ~ simulated_data$rand_int,
                                            simulated_data[[column]] > 1 ~ simulated_data[[column]] - simulated_data$rand_int,
                                            TRUE ~ simulated_data[[column]])
      
      #p_change <- (mean(simulated_data[[column]]) - mean(main_lice[[column]])) / mean(main_lice[[column]]) * 100
      #p_changes[sample_size] <- p_change
      
      #simulated_dataset$p_change <- p_change
      simulated_data$year <- simulated_data$year + 7
      
      
      #resample original data (with 2015-2021 years)
      past_data <- re_sample(df = main_lice, 1) # Should I always keep sample size = 1 for the past data?
      
      #add farm effect binary variable
      simulated_data$farm <- 1
      past_data$farm <- 0
      
      simulated_data <- bind_rows(simulated_data, past_data)
      simulated_data$farm <- as_factor(simulated_data$farm)
      simulated_data$seine_id <- factor(paste0(simulated_data$year, 
                                               simulated_data$seine_id))
      simulated_data$year <- factor(simulated_data$year)
      
      # create name of simulated dataset based on effect and sample size 
      simulated_datasets[[paste0("mean_", i, "_n_",
                                 sample_size)]] <- simulated_data

    }
   #print(mean(unlist(p_changes)))
  }
  return(simulated_datasets)
}

# Create function to fit each combined dataset (future simulated data with farm effect, 
# and past data without farm effect) with one model that includes the farm effect
# and one that doesn't, compare the fit, and store the p_value from MLE (anova) comparison

compare_models <- function(simulated_datasets, column) {
  # Initialize list to store model comparison results
  model_comparison_results <- list()
  
  #Loop through each dataset in the list
  for (i in seq_along(simulated_datasets)) {
    dataset <- simulated_datasets[[i]]
    
    # Had to remove year as fixed effect due to colinearity with farm effect
    with_farm <- glmmTMB(as.formula(paste0(column, "~ species + farm + (1 | year/seine_id)")),
                         ziformula = ~1,
                         data = dataset, family = nbinom2, na.action = na.omit)
    
    no_farm   <- glmmTMB(as.formula(paste0(column, "~ species + (1 | year/seine_id)")),
                         ziformula = ~1,
                         data = dataset, family = nbinom2, na.action = na.omit)
    
    # Compare models using MLE via the anova function
    #anova_comp <- anova(with_farm, no_farm, na.action = na.omit)
    anova_comp <- anova(no_farm, with_farm, na.action = na.omit)
    
    # Store the result in the list
    model_comparison_results[[i]] <- list(
      dataset_name = names(simulated_datasets)[i],
      #extract p_value
      p_value = anova_comp[2,8]
    )
  }
  
  # Return the list of model comparison results
  return(model_comparison_results)
}

# Now combine all the functions into a function that iterates through every
# simulated dataset n_iteration times
calculate_p_proportions <- function(sample_sizes, column, n_iterations = 10) {
  #tic("All iterations")
  
  # calculate and store number of model comparisons that need to be computed
  n_models <- length(sequence) * length(sample_sizes)
  
  # Initialize a matrix to store p-values for each model comparison
  p_values <- matrix(0, nrow = n_iterations, ncol = n_models)
  
  # Create parallel loop to simulate datasets, compare models and store p values,
  # and store the percent of iterations that resulted in a p_value < 0.05
  p_values <- foreach(i = 1:n_iterations, .combine = "rbind") %dopar% {
    
    # Simulate datasets
    simulated_datasets <- simulate_data_v3(sample_sizes, column = column)
    #p_change <- simulated_datasets[["mean_1_n_2"]][["p_change"]][[1]]
    # Compare models and store the results
    model_comparison_results <- compare_models(simulated_datasets, column)
  
    # Extract p-values from model_comparison_results
    p_value_vector <- sapply(model_comparison_results, function(x) x$p_value)
    #dataset_names <- sapply(model_comparison_results, function(x) x$dataset_name)
    # Return the p_value_vector
    #names(p_value_vector) <- dataset_names
    #tibble(p_value_vector, p_change)
    p_value_vector

  }
  
 
  
  # Calculate the proportion of times p is < 0.05 for each model comparison
  p_less_than_0_05 <- colSums(p_values < 0.05, na.rm = TRUE) / n_iterations
  #p_less_than_0_05 <- summarise(p_value$p_value_vector < 0.05) / n()
  
  # Generate a data frame with sample_size and effect columns to store p_values
  # expand grid will vary sample size first while keeping effect constant which
  # is the same order the nested for loops in the simulate data function are nested
  sample_size_effect <- expand.grid(sample_size = sample_sizes, 
                                    effect = c(-44.4, -42.7, -41, -38.2, -33.5, -31.6, -29.8, -29.4, -27.25, -23.7, -21.4),
                                    stringsAsFactors = FALSE)
  
  
  
  # Add the p_proportions column to the sample_size_effect data frame
  sample_size_effect$p_proportions <- p_less_than_0_05
  
  #toc()
  print(sample_size_effect)
  return(sample_size_effect)
  
}

# Define number of iterations to run
n_iterations <- 100

# call the function to create a data frame of p_value proportions for each combination
# of effect size(mean_multiplier) and sample_size
p_proportions <- calculate_p_proportions(sample_sizes = sample_sizes, 
                                         column = "cal_mot",
                                         n_iterations = n_iterations
                                         )

beepr::beep(2)

print(p_proportions)
# #convert multiplier to percent reduction 
p_proportions$effect <- abs(p_proportions$effect)
#p_proportions$effect <- (1 - p_proportions$effect) * 100
#p_proportions$effect <- (100 - p_proportions$effect) / 100
# 
# #convert sample size fraction to actual number of samples in dataset
# p_proportions$sample_prop <- p_proportions$sample_size 
# p_proportions$sample_size <- p_proportions$sample_size * 4084


# Create a heatmap of p_proportions as a function of effect and sample_size
heatmap_plot <- ggplot(p_proportions, aes(x = factor(effect), y = factor(sample_size), fill = p_proportions)) +
  geom_tile() +
  scale_fill_gradient2(low = "grey", high = "black", limits = c(0, 1), name = "Power") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  labs(x = "% decrease in sealice abundance", y = "Sample Size")

# Print the heatmap
print(heatmap_plot)
beepr::beep(5)





