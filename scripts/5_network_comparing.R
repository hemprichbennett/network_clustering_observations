#### Header ####
## Project: network_resolution
## Script purpose: Making figures for the clustering data
## Date: 2019-04-05
## Author: Dave Hemprich-Bennett (hemprich.bennett@gmail.com)
## Notes This array script is used to generate all the possible outcomes of
## clumping, for plotting with a later script
##################################################

##### set up session #####

library(tidyverse)
library(tidyverse)

original_values <- read_csv("results/original_values.csv")

clumped_df <- read_csv("results/clumped_values.csv")

i <- 1

#### Find which comparison is being done in this iteration #####
possible_pairs <- expand.grid(
  net1_ID = unique(original_values$dataset),
  net2_ID = unique(original_values$dataset),
  net1_clust = unique(clumped_df$threshold),
  net2_clust = unique(clumped_df$threshold),
  net_level = unique(clumped_df$net_level),
  taxonomic_level = unique(clumped_df$taxonomic_level),
  metric = colnames(original_values)[6:length(colnames(original_values))]
)

# Remove the self-comparisons
possible_pairs <- possible_pairs[-which(possible_pairs$net1_ID == possible_pairs$net2_ID), ]


# remove quantitative metrics involving canary or galapagos networks, as they're binary

metric_types <- read.csv('data/metric_types.csv',
                         stringsAsFactors = F, row.names = NULL, header = T)

quant_metrics <- metric_types %>%
  # isolate the desired column
  pull(quantitative) %>%
  # turn the names into the ones expected in desired_combinations
  str_replace('_', '.') %>%
  str_replace('higher', 'HL') %>%
  str_replace('lower', 'LL')


possible_pairs <- possible_pairs %>%
# filter out canary and galapagos networks from comparisons
    filter(! (net1_ID %in% c('canary', 'galapagos') & 
           metric %in% quant_metrics)) %>%
    filter(! (net2_ID %in% c('canary', 'galapagos') & 
           metric %in% quant_metrics)) %>%
  # filter out gen1 and hrat from order-level bird comparisons at 1 similarity, 
  # as theres only one order (passeriformes)
  filter(! (net1_ID %in% c('gen1', 'hrat') & 
              net1_clust == 1 & 
              net_level == 'upper' & 
              taxonomic_level == 'order')) %>%
  filter(! (net2_ID %in% c('gen1', 'hrat') & 
              net2_clust == 1 & 
              net_level == 'upper' &
              taxonomic_level == 'order'))


conclusion_check <- function(i) {
  ##### Set the basic parameters ####
  desired_combination <- possible_pairs[i, ]

  net1_ID <- as.character(desired_combination$net1_ID)
  net2_ID <- as.character(desired_combination$net2_ID)
  net1_clust <- desired_combination$net1_clust
  net2_clust <- desired_combination$net2_clust
  net_level <- as.character(desired_combination$net_level)
  taxonomic_level <- as.character(desired_combination$taxonomic_level)
  metric <- as.character(desired_combination$metric)

  ###### Find the original 'real' values #####

  net1_trueval <- original_values[which(original_values$dataset == net1_ID), metric]
  net2_trueval <- original_values[which(original_values$dataset == net2_ID), metric]
  if (is.na(net1_trueval) | is.na(net2_trueval)) {
    biggest <- "neither"
  } else if (net1_trueval > net2_trueval) {
    biggest <- "net1"
  } else if (net2_trueval > net1_trueval) {
    biggest <- "net2"
  } else if (net1_trueval == net2_trueval) {
    biggest <- "neither"
  }

  ##### Find the desired values #####
  net1_subset <- clumped_df[which(clumped_df$dataset == net1_ID &
    clumped_df$threshold == net1_clust &
    clumped_df$net_level == net_level &
    clumped_df$taxonomic_level == taxonomic_level), metric]

  net2_subset <- clumped_df[which(clumped_df$dataset == net2_ID &
    clumped_df$threshold == net2_clust &
    clumped_df$net_level == net_level &
    clumped_df$taxonomic_level == taxonomic_level), metric]

  # Generate all possible pairs of simulated networks
  all_pairs <- expand.grid(net1_subset, net2_subset)
  colnames(all_pairs) <- c("net1", "net2")

  # get the indicies of all erroneous conclusions
  if (biggest == "net1") {
    erroneous_conclusions <- which(all_pairs$net2 > all_pairs$net1)
  } else if (biggest == "net2") {
    erroneous_conclusions <- which(all_pairs$net1 > all_pairs$net2)
  } else if (biggest == "neither") {
    erroneous_conclusions <- which(all_pairs$net1 != all_pairs$net2)
  }

  # So was a bad comparison possible?
  if (length(erroneous_conclusions) > 0) {
    false_conclusion_possible <- T
  } else {
    false_conclusion_possible <- F
  }

  # Make an outstring
  outstr <- cbind(desired_combination, false_conclusion_possible)

  # Return it
  return(outstr)
}
#for(i in 1:nrow(possible_pairs)){ print(conclusion_check(i))}



# Now do the comparisons
start_time <- Sys.time()
out_list <- lapply(seq(1, nrow(possible_pairs)), function(x) conclusion_check(x))
end_time <- Sys.time()

cat("duration was", end_time - start_time)
# Combine the (MANY) outputs into a df
out_df <- data.table::rbindlist(out_list)

# Save
write.csv(out_df, "results/possible_conclusions.csv")
