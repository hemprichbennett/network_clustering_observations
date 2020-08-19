#### Header ####
## Project: network_clustering_real
## Script purpose: obtaining metadata for all networks
## Date: 2019-03-20
## Author: Dave Hemprich-Bennett (hemprich.bennett@gmail.com)
## Notes
##################################################

###### setup #####

library(taxize)
library(plyr)



dirpath <- "data/input_data/"
filenames <- list.files(path = dirpath, pattern = "*.csv", include.dirs = T)

if (TRUE %in% grepl("metadata", filenames)) {
  filenames <- filenames[-grep("metadata", filenames)]
}



filepaths <- paste0(dirpath, filenames)


###### get and save taxonomy #####
# Its an ugly nested loop, but as it needs user input its the easiest way

# for both taxonomic levels for each file, we need to obtain and save the
# metadata
for (chosen_file in 1:length(filenames)) {
  net <- read.csv(filepaths[chosen_file], header = T, stringsAsFactors = F, row.names = 1)
  netname <- filenames[chosen_file]
  netname <- gsub("\\.csv", "", netname)
  netname <- gsub("matrix_", "", netname)
  network_levels <- c("upper", "lower")

  # selected_level <- 1
  for (selected_level in 1:2) {
    target_level <- network_levels[selected_level]
    cat('working on ',netname, ', ', target_level, '\n' )

    if (target_level == "upper") {
      target_sp <- colnames(net)
    } else {
      target_sp <- rownames(net)
    }
    target_sp <- gsub("\\.", " ", target_sp)


    taxonomic_info <- classification(target_sp, db = "itis")

    desired_levels <- c("family", "order", "genus", "species")
    out_file <- matrix(nrow = 0, ncol = length(desired_levels))
    colnames(out_file) <- desired_levels

    # Keep only the taxonomic information which we're after, as there are different
    # numbers of fields for some species
    for (species in 1:length(taxonomic_info)) {
        sp_name <- names(taxonomic_info)[species]
      if (is.null(dim(taxonomic_info[[species]]))) {
        out_vec <- c(NA, NA, NA, sp_name)
        cat("No taxonomy found for ", sp_name, " in ", netname, " ", "target_level")
      }
      else {
        desired_rows <- which(taxonomic_info[[species]]$rank %in% desired_levels)
        
        # some species don't have all 4 values, so need to be worked on manually
        # later
        
        if(length(desired_rows) <4){
          out_vec <- c(NA, NA, NA, sp_name)
          cat("No taxonomy found for ", sp_name, " in ", netname, " ", "target_level")
          next()
        }

        desired_vals <- taxonomic_info[[species]][desired_rows, c(1, 2)]

        # out_file <- rbind(out_file, rep(NA, ncol(out_file)))

        out_vec <- rep(NA, ncol(out_file))

        # As the desired_vals file isn't in the right order, the following orders the
        # values taxonomically before putting them into out_file
        for (column in 1:ncol(out_file)) {
          column_name <- colnames(out_file)[column]
          out_vec[column] <- desired_vals$name[which(desired_vals$rank == column_name)]
        }
      }

      out_file <- rbind(out_file, out_vec)
    }
    out_name <- paste0(
      "data/input_data/metadata_", netname, "_",
      target_level, ".csv"
    )
    print(out_name)
    write.csv(out_file, out_name, row.names = FALSE)
  }
  
}
