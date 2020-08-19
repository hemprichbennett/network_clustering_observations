#### Header ####
## Project: network clustering (real)
## Script purpose: generating the values for plotting, by analysing how many 
## of each pair of clustering threshold can give an erroneous conclusion
## Date: 2019-04-16
## Author: Dave Hemprich-Bennett (hemprich.bennett@gmail.com)
## Notes
##################################################

##### Setup ####
library(tidyverse)
# Read in data
possible_conclusions <- read.csv("results/possible_conclusions.csv",
  stringsAsFactors = F, row.names = 1
)

# Make the df to use as input parameters for the function
desired_combinations <- unique(possible_conclusions[, 3:7])



###### Make function #####

# This function takes all of the comparisons done in the previous script and
# finds what proportion of them were able to give dodgy conclusions

percentage_calc <- function(i) {
  # Isolate the desired values
  this_combination <- desired_combinations[i, ]
  net1_clust <- this_combination$net1_clust
  net2_clust <- this_combination$net2_clust
  net_level <- this_combination$net_level
  taxonomic_level <- this_combination$taxonomic_level
  metric <- this_combination$metric

  # Find the matching results
  this_comparison <- possible_conclusions[
    which(possible_conclusions$net1_clust == net1_clust &
      possible_conclusions$net2_clust == net2_clust &
      possible_conclusions$net_level == net_level &
      possible_conclusions$taxonomic_level == taxonomic_level &
      possible_conclusions$metric == metric),
  ]

  # Calculate the percentage
  bad_comparisons <- length(which(this_comparison$false_conclusion_possible == T))

  bad_percentage <- bad_comparisons / nrow(this_comparison) * 100

  # output the findings
  this_combination$percentage <- bad_percentage
  return(this_combination)
}


# Remove quantitative comparisons which shouldn't be made -----------------

# because the canary and galapagos networks were binary, they should not be 
# used in any quantitative comparisons. We here remove them from the 
# possible_conclusions data frame before using the comparison_calc function


metric_types <- read.csv('data/metric_types.csv',
                         stringsAsFactors = F, row.names = NULL, header = T)

quant_metrics <- metric_types %>%
  # isolate the desired column
  pull(quantitative) %>%
  # turn the names into the ones expected in desired_combinations
  str_replace('_', '.') %>%
  str_replace('higher', 'HL') %>%
  str_replace('lower', 'LL')


# possible_conclusions <- possible_conclusions[-
#   which(possible_conclusions$net1_ID %in% c('canary', 'galapagos') &
#         possible_conclusions$metric %in% quant_metrics),]
# 
# possible_conclusions <- possible_conclusions[-
#   which(possible_conclusions$net2_ID %in% c('canary', 'galapagos') &
#           possible_conclusions$metric %in% quant_metrics),]

####### Analysis #####
start_time <- Sys.time()
out_list <- lapply(seq(1, nrow(desired_combinations)), function(x) percentage_calc(x))
out_df <- data.table::rbindlist(out_list)
end_time <- Sys.time()
cat("Analysis took", end_time - start_time)


write.csv(out_df, "results/percentage_errors.csv")
