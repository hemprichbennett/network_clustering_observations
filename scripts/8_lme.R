#### Header ####
## Project: network clustering (real)
## Script purpose: Running mixed effects models
## Date: 2019-04-17
## Author: Dave Hemprich-Bennett (hemprich.bennett@gmail.com)
## Notes


##### Setup #####

library(lme4)
library(lmerTest)
library(tidyverse)
library(broom)

# Function for capitalising strings
firstup <- function(x) {
  substr(x, 1, 1) <- toupper(substr(x, 1, 1))
  x
}


wide_df <- read.csv("results/clumped_values.csv",
  header = T,
  stringsAsFactors = F
)

wide_df$iteration <- NULL
long_df <- reshape2::melt(wide_df,
  id.vars = c(
    "threshold", "dataset",
    "net_level", "taxonomic_level"
  )
)

long_df <- rename(long_df, metric = variable)
long_df$metric <- gsub("H2", "H2'", long_df$metric)

long_df$metric <- as.character(long_df$metric)



##### Remove datapoints and sets which will crash the model #####

long_df <- long_df %>%
  group_by(metric) %>%
  mutate(
    n_datasets = n_distinct(dataset),
    n_values = n_distinct(value)
  ) %>%
  filter(
    !is.na(value) & # Remove na values
      n_datasets > 1 & # Remove metrics which only have one valid dataset
      n_values > 1 # Remove metrics which only have one valid value
  )




# Remove qualitative metrics from quantitative comparisons ----------------

metric_types <- read.csv('data/metric_types.csv',
                         stringsAsFactors = F, row.names = NULL, header = T)

quant_metrics <- metric_types %>%
  # isolate the desired column
  pull(quantitative) %>%
  # turn the names into the ones expected in desired_combinations
  str_replace('_', '.') %>%
  str_replace('higher', 'HL') %>%
  str_replace('lower', 'LL')



long_df <- long_df %>%
  filter(!(dataset %in% c('canary', 'galapagos') & metric %in% quant_metrics))

# check which values have been retained
retained <- long_df %>% 
  group_by(dataset, metric) %>% 
  summarise(n = n()) %>% 
  pivot_wider(names_from = metric, values_from = n)

##### do the modelling #####

list_lmer <- function(met) {
  desired_subset <- filter(long_df, metric == met)


    mod <- lmer(value ~ threshold +
      (1 | dataset) +
      net_level +
      taxonomic_level,
    data = desired_subset
    )



  return(mod)
}


# kill compartment diversity and isa, as they're both always 1
long_df <- long_df %>%
  filter(!metric %in% c(
    "compartment.diversity",
    "interaction.strength.asymmetry"
  ))


# Make fixed models -------------------------------------------------------


fixed_mods <- lapply(
  unique(long_df$metric),
  function(x) list_lmer(met = x)
)

names(fixed_mods) <- unique(long_df$metric)

fixed_success <- tibble(
  metric = names(fixed_mods),
  success = NA
)

for (i in 1:length(fixed_mods)) {
  if (TRUE %in% grepl("failed", summary(fixed_mods[[i]]))) {
    fixed_success$success[i] <- FALSE
  } else {
    fixed_success$success[i] <- TRUE
  }
}




# compare model type successes --------------------------------------------
group_by(fixed_success, success) %>% summarise(n = n())
#group_by(random_success, success) %>% summarise(n = n())

filter(fixed_success, success == FALSE)




# Try and make something useful out of the outputs... ---------------------

summary_list <- list()

for(i in 1:length(fixed_mods)){
  broom_tibble <- broom.mixed::tidy(fixed_mods[[i]]) %>%
    # Remove the intercept, its not so informative for plotting
    filter(term != '(Intercept)') %>%
    mutate(metric = names(fixed_mods)[i],
           term = gsub('threshold', 'Clustering threshold', term),
           term = gsub('net_levelupper', 'Network level (upper)', term),
           term = gsub('taxonomic_level', 'Taxonomic level: ', term),
           term = gsub('sd__\\(Intercept\\)', 'Random effect: dataset', term),
           term = gsub('sd__Observation', 'Residual', term),
           term = gsub('\\(Intercept\\)', 'Intercept', term),
           effect = gsub('ran_pars', 'Random', effect),
           effect = gsub('fixed', 'Fixed', effect)) 
    
    
  for_ms <- broom_tibble %>% 
    rename_all(firstup) %>%
    mutate_if(is.double, round, 2) %>%
    mutate_all(~ gsub('^0$', '<0.001', .x))
  
  write.csv(for_ms,
    paste0("results/model_outputs/", names(fixed_mods)[i],".csv"),
    row.names = F)
  
  summary_list[[i]] <- broom_tibble
  }

summary_df <- do.call(rbind, summary_list)
write_csv(summary_df, 'results/model_outputs/all_model_outputs.csv')
# ggplot(summary_df,
#        aes(y = metric, x = term, fill = estimate)) +
#   geom_tile()


# Function for making INDIVIDUAL barplots
barplot <- function(met, save_file = F){

    outplot <- ggplot(filter(summary_df, metric == met)
         , aes(x = term, y = estimate)) +
      ggtitle(met)+
    geom_bar(stat = 'identity') +
    coord_flip() + theme_bw() +
      labs(x = 'Estimate', y = 'Term')
    if(save_file == T){
      ggsave(paste0('figures/lm_plots/',met, '.pdf'), outplot)
    }else{
      return(outplot)
    }
  
  
}

firstup <- function(x) {
  substr(x, 1, 1) <- toupper(substr(x, 1, 1))
  x
} # A function to capitalise the metric names when making plots


barplot('connectance')
lapply(unique(summary_df$metric), function(x) barplot(met = x, save_file = T))


# Now isolate a few metrics of specific interest, plot them ---------------

desired_metrics <- c(
  "Functional complementarity, higher level", "Modularity", "Nestedness", 
  "NODF", "Robustness, higher level", "H2'"
)

desired_subset <- summary_df %>%
  # Improve the metric names
  mutate(metric = gsub('\\.HL', ', higher level', metric),
    metric = gsub('\\.', ' ', metric),
    metric = firstup(metric)) %>%
  # Filter so we only have the desired metrics
  filter(metric %in% desired_metrics) %>%
  # Rename fun comp, for nicer plotting
  mutate(metric = gsub("Functional complementarity, higher level", 
                       "Functional complementarity,\nhigher level", metric)) %>%
  # Reorder the factors to z-a, as ggplot plots them in an annoying order otherwise
  mutate(term = as.factor(term),
           term = factor(term, levels=rev(levels(term))))

facet_plot <- ggplot(desired_subset, aes(x = term, y = estimate)) +
  geom_bar(stat = 'identity') +
  coord_flip() + theme_bw() +
  labs(y = 'Estimate', x = 'Term') +
  facet_wrap(. ~ metric, scales = 'free_x', ncol = 2)
facet_plot

ggsave('figures/grid_lm.jpg', facet_plot, height = 8, width = 8)
