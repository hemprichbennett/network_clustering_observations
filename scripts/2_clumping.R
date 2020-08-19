#### Header ####
## Project: Network_OTUs
## Script purpose: Clumping taxonomically seven different networks,
## in a HPC array
## Date: 2019-04-04
## Author: Dave Hemprich-Bennett (hemprich.bennett@gmail.com)
## Notes
##################################################

start_time <- Sys.time()


##### Setup #####

dir <- getwd()
basedir <- strsplit(dir, split = "/")[[1]][2]
print(basedir)
if (grepl("data", basedir)) {
  loc <- 'apocrita'
  library(here, lib.loc = "/data/home/btw863/r_packages/")
  require(methods)
  library(crayon, lib.loc = "/data/home/btw863/r_packages/")
  library(vctrs, lib.loc = "/data/home/btw863/r_packages/")
  library(permute, lib.loc = "/data/home/btw863/r_packages/")
  library(lattice, lib.loc = "/data/home/btw863/r_packages/")
  library(vegan, lib.loc = "/data/home/btw863/r_packages/")
  library(statnet.common, lib.loc = "/data/home/btw863/r_packages/")
  library(network, lib.loc = "/data/home/btw863/r_packages/")
  library(sna, lib.loc = "/data/home/btw863/r_packages/")
  library(bipartite, lib.loc = "/data/home/btw863/r_packages/")
  library(stringr, lib.loc = "/data/home/btw863/r_packages/")

  in_arg <- commandArgs(trailingOnly = TRUE)
  
  in_arg <- as.numeric(in_arg[1])
  cat('in_arg is', in_arg, '\n')
} else if (grepl("home", basedir)){
  loc <- 'ARCUS'
  setwd('/home/zool2291/projects/network_clustering_real')
  library(here, lib.loc = "/home/zool2291/r_packages/3.5_and_others/")
  require(methods)
  library(crayon, lib.loc = "/home/zool2291/r_packages/3.5_and_others/")
  library(vctrs, lib.loc = "/home/zool2291/r_packages/3.5_and_others/")
  library(permute, lib.loc = "/home/zool2291/r_packages/3.5_and_others/")
  #library(lattice, lib.loc = "/home/zool2291/r_packages/3.5_and_others/")
  library(vegan, lib.loc = "/home/zool2291/r_packages/3.5_and_others/")
  library(statnet.common, lib.loc = "/home/zool2291/r_packages/3.5_and_others/")
  library(network, lib.loc = "/home/zool2291/r_packages/3.5_and_others/")
  library(sna, lib.loc = "/home/zool2291/r_packages/3.5_and_others/")
  library(bipartite, lib.loc = "/home/zool2291/r_packages/3.5_and_others/")
  library(stringr, lib.loc = "/home/zool2291/r_packages/3.5_and_others/")
  
  in_arg <- commandArgs(trailingOnly = TRUE)
  in_arg <- as.numeric(in_arg[1])
  cat('in_arg is', in_arg, '\n')
}else {
  library("here")
  library(vegan)
  library(bipartite)
  library(stringr)
  library(reshape2)
  loc <- 'local'

  in_arg <- commandArgs(trailingOnly = TRUE)
  in_arg <- as.numeric(in_arg[1])
  cat('in_arg is', in_arg, '\n')
}



getwd()
inpath <- "data/input_data/"
# The datasets to use
matrixes_names <- list.files(path = inpath, pattern = "matrix")
print(matrixes_names)
metadata_names <- list.files(path = inpath, pattern = "metadata")
network_levels <- c("upper", "lower") # The network levels we can clump
iterations <- seq(1, 100, 1) # The number of iterations we want
thresholds <- seq(0.1, 1, 0.1) # The thresholds we can clump by
taxonomic_levels <- c("order", "family", "genus")

# For the random elements, lets set a seed based on the iteration being run
set.seed(Sys.time())


param_list <- list(
  matrixes_names, network_levels,
  iterations, thresholds, taxonomic_levels
)
names(param_list) <- c(
  "matrixes_names", "network_levels", "iterations",
  "thresholds", "taxonomic_levels"
)
combinations <- do.call(expand.grid, param_list)

##### Define clumping function #####


clumper <- function(input_matrix, clump_to_level, portion_to_clump, network_level, metadata) {
  # cat('clump_to_level is ',clump_to_level,'\n')
  
  # Find if network is weighted or binary
  weighting <- ifelse(sum(range(input_matrix)) == 1, 'binary', 'weighted')
  
  if (network_level == "upper") {
    new_matrix <- matrix(nrow = nrow(input_matrix), ncol = 0)
    rownames(new_matrix) <- rownames(input_matrix)
    nodes_to_clump <- sample(1:ncol(input_matrix), round(ncol(input_matrix) * portion_to_clump))
    a <- 1
    for (i in 1:ncol(input_matrix)) {
      if (i %in% nodes_to_clump) {
        c <- gsub("\\.", " ", colnames(input_matrix)[i])
        r <- which(metadata$species == c)
        newname <- metadata[r, clump_to_level]
        #cat("i is", i, "c is", c, "r is ", r, "newname is", newname, "\n")
        if (ncol(new_matrix) > 0) {
          print(c)
          print(newname)
          if (newname %in% colnames(new_matrix)) {
            to_merge <- which(colnames(new_matrix) == newname)
            new_matrix[, to_merge] <- new_matrix[, to_merge] + input_matrix[, i]
          } else {
            new_matrix <- cbind(new_matrix, as.numeric(input_matrix[, i]))
            colnames(new_matrix)[a] <- newname
            a <- a + 1
          }
        } else {
          new_matrix <- cbind(new_matrix, as.numeric(input_matrix[, i]))
          colnames(new_matrix)[a] <- newname
          a <- a + 1
        }
      } else {
        new_matrix <- cbind(new_matrix, as.numeric(input_matrix[, i]))
        colnames(new_matrix)[a] <- colnames(input_matrix)[i]
        a <- a + 1
      }
      # print(new_matrix)
    }
  }
  if (network_level == "lower") {
    new_matrix <- matrix(ncol = ncol(input_matrix), nrow = 0)
    colnames(new_matrix) <- colnames(input_matrix)
    nodes_to_clump <- sample(1:nrow(input_matrix), round(nrow(input_matrix) * portion_to_clump))
    a <- 1
    for (i in 1:nrow(input_matrix)) {
      if (i %in% nodes_to_clump) {
        c <- gsub("\\.", " ", rownames(input_matrix)[i])
        r <- which(metadata$species == c)
        newname <- metadata[r, clump_to_level]
        # cat('i is', i, 'c is', c, 'r is',r,  'newname is', newname,'\n')
        if (nrow(new_matrix) > 0) {
          
          if (newname %in% rownames(new_matrix)) {
            to_merge <- which(rownames(new_matrix) == newname)
            #   print(to_merge)
            #    print(new_matrix[to_merge, ] + as.numeric(input_matrix[i, ]))
            new_matrix[to_merge, ] <- new_matrix[to_merge, ] + as.numeric(input_matrix[i, ])
            #    print(new_matrix)
          } else {
            new_matrix <- rbind(new_matrix, as.numeric(input_matrix[i, ]))
            #   cat('else 1, newname is', newname, '\n')
            rownames(new_matrix)[a] <- newname
            a <- a + 1
            #  print(new_matrix)
          }
        } else {
          new_matrix <- rbind(new_matrix, as.numeric(input_matrix[i, ]))
          # cat('else 2, newname is', newname, '\n')
          rownames(new_matrix)[a] <- newname
          a <- a + 1
        }
      } else {
        new_matrix <- rbind(new_matrix, as.numeric(input_matrix[i, ]))
        # cat('else 3, nextname is', rownames(input_matrix)[i], '\n')
        rownames(new_matrix)[a] <- rownames(input_matrix)[i]
        a <- a + 1
      }
      # print(new_matrix)
    }
  }
  
  if(weighting == 'binary'){
    new_matrix <- ifelse(new_matrix == 0, 0, 1) # Make it binary, as the input data was binary.
  }
  
  
  new_matrix <- as.matrix(new_matrix) # Make sure its a matrix
  return(new_matrix)
}


clumping_analysis <- function(input_argument){
  
  # Isolate the variables we're after in this iteration
  desired_combo <- combinations[input_argument, ]
  print(desired_combo)
  dataset <- gsub("matrix_|\\.csv", "", desired_combo$matrixes_names)
  net_level <- as.character(desired_combo$network_levels)
  iteration <- desired_combo$iterations
  threshold <- desired_combo$thresholds
  taxonomic_level <- as.character(desired_combo$taxonomic_levels)
  
  # These iterations crash as theres only one order
  if (dataset == "gen1" & net_level == "upper" & taxonomic_level == "order" & threshold == "1") {
    print("skipping gen1 order")
    return()
  }
  if (dataset == "hrat" & net_level == "upper" & taxonomic_level == "order" & threshold == "1") {
    print("skipping hrat order")
    return()
  }
  # Read in the matrix and appropriate metadata
  original_matrix <- read.csv(paste0(inpath, desired_combo$matrixes_names),
                              header = T, stringsAsFactors = F, row.names = 1
  )
  
  colnames(original_matrix) <- gsub("\\.", " ", colnames(original_matrix))
  rownames(original_matrix) <- gsub("\\.", " ", rownames(original_matrix))
  
  metadata_name <- metadata_names[grep(dataset, metadata_names)]
  metadata_name <- metadata_name[grep(net_level, metadata_name)]
  metadata_file <- read.csv(paste0(inpath, metadata_name),
                            header = T, stringsAsFactors = F, row.names = NULL
  )
  
  
  # Clump the matrix by the given parameters
  clumped_mat <- clumper(
    input_matrix = original_matrix,
    portion_to_clump = threshold,
    network_level = net_level,
    clump_to_level = taxonomic_level,
    metadata = metadata_file
  )
  out_mat_name <- paste0('./data/clumped_matrixes/', iteration, '_',
                        dataset, '_', threshold, '_', net_level, '_',
                        taxonomic_level,'.csv')
  print(out_mat_name)
  
  write.csv(clumped_mat, file = out_mat_name)


  # Calculate the metrics for the clumped matrix
  # Sadly we can't do modularity for tiny matrices

  # this checks to see how many rows and columns are left after excluding all 
	# that are just full of ones: modularity excludes them by default but breaks if
  # it tries to analyse a tiny matrix

	n_row <- nrow(clumped_mat)
	n_col <- ncol(clumped_mat)

	rowSums(clumped_mat)

	# the number of rows which are not filled with 1 or 0 (these are removed and
	# can lead to a tiny matrix, which modularity understandably hates)
	n_good_rows <- length((rowSums(clumped_mat) < n_col) & 
	                        (rowSums(clumped_mat) >0))
	
	# same but for columns
	n_good_cols <- length((colSums(clumped_mat) < n_row) & 
	                        (colSums(clumped_mat) >0))
	  

  if (n_good_cols > 3 & n_good_rows > 3) {
    clumped_stats <- c(
      networklevel(clumped_mat),
      slot(bipartite::computeModules(web = clumped_mat), "likelihood")
    )
  } else {
    clumped_stats <- c(
      networklevel(clumped_mat),
      NA
    )
  }
  
  
  
  
  # Format the output
  names(clumped_stats)[length(clumped_stats)] <- "modularity" # Modularity needs naming, it doesn't have one by default
  
  
  clumped_df <- data.frame(
    dataset = dataset, iteration = iteration,
    threshold = threshold, net_level = net_level,
    taxonomic_level = taxonomic_level, t(clumped_stats)
  )
  
  
  # Save the output file
  outfilename <- paste(dataset, iteration, threshold, net_level, taxonomic_level, sep = "_")
  
  write.csv(clumped_df, paste0("./data/clumped_output_data/", outfilename, ".csv"))
  
}

# Use the function
if(loc == 'ARCUS'){
  cat('running on ARCUS, in_arg is', in_arg, '\n')
  start_pos <- (in_arg * 54) -53
  end_pos <- start_pos + 53
  cat('start_pos is', start_pos, 'end_pos is', end_pos, '\n')
  for(i in start_pos:end_pos){clumping_analysis(input_argument = i)}
}else if(loc == 'apocrita'){
  clumping_analysis(input_argument = in_arg)
}else if(loc == 'local'){
  lapply(seq(1:nrow(combinations)), function(i) clumping_analysis(input_argument = i))
}

#clumping_analysis(14800)
#end_time <- Sys.time()

#cat("Iteration took ", end_time - start_time, " of whatever unit of system time this is\n")
