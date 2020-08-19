#### Header ####
## Project: network resolution
## Script purpose: calculating the original values for each of the observation
## networks
## Date: 2019-04-05
## Author: Dave Hemprich-Bennett (hemprich.bennett@gmail.com)
## Notes
##################################################

library(bipartite)

###### read in and format data #####
inpath <- "data/input_data/"
filenames <- list.files(path = inpath, pattern = "matrix")

all_files <- lapply(paste0(inpath, filenames), function(f)
  read.csv(f, header = T, stringsAsFactors = F, row.names = 1))

names(all_files) <- gsub("matrix_|\\.csv", "", filenames)

###### Calculations #####

# Make a function so we can vectorise the calculations
clumped_stats <- function(input_matrix) {
  outstr <- c(
    networklevel(input_matrix),
    slot(bipartite::computeModules(web = input_matrix), "likelihood")
  )
  # Modularity needs naming, it doesn't have one by default
  names(outstr)[length(outstr)] <- "modularity" 


  out_df <- data.frame(
    dataset = NA, iteration = 0,
    threshold = 0, net_level = NA,
    taxonomic_level = NA, t(outstr)
  )
  return(out_df)
}



# Calculate
all_stats <- lapply(all_files, clumped_stats)

stats_df <- do.call(rbind, all_stats)

stats_df$dataset <- row.names(stats_df)


write.csv(stats_df, "results/original_values.csv", row.names = F)
