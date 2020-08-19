#### Header ####
## Project: network_resolution
## Script purpose: combining together the ridiculous number of clumping outputs
## Date: 2019-04-05
## Author: Dave Hemprich-Bennett (hemprich.bennett@gmail.com)
## Notes
##################################################

###### Read in data #####
inpath <- "data/clumped_output_data/"
filenames <- list.files(path = inpath)

filenames <- paste0(inpath, filenames)

all_files <- lapply(filenames, function(f)
  read.csv(f, header = T, stringsAsFactors = F, row.names = 1))

##### Combine and save #####
giant_df <- do.call(rbind, all_files)

write.csv(giant_df, "results/clumped_values.csv", row.names = F)
