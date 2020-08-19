#### Header ####
## Project: network clustering (real)
## Script purpose: ploting the results from the previous script
## Date: 2019-04-16
## Author: Dave Hemprich-Bennett (hemprich.bennett@gmail.com)
## Notes
##################################################

##### Setup ####
library(ggplot2)
library(scales)
library(magrittr)
library(dplyr)

firstup <- function(x) {
  substr(x, 1, 1) <- toupper(substr(x, 1, 1))
  x
} # A function to capitalise the metric names when making plots


all_df <- read.csv("results/percentage_errors.csv",
  stringsAsFactors = F,
  header = T, row.names = 1
)

# Capitalise variables for plotting
all_df$net_level <- firstup(all_df$net_level)
all_df$taxonomic_level <- firstup(all_df$taxonomic_level)

# change network level naming: upper and lower isn't that informative!
all_df$net_level <- ifelse(all_df$net_level == 'Upper', 'Seed consumer', 'Plant')

# reorder factors
all_df$taxonomic_level <- factor(all_df$taxonomic_level,
  levels = c("Order", "Family", "Genus")
)

all_df$metric <- factor(all_df$metric, 
                        levels = unique(all_df$metric)[
                          order(unique(all_df$metric))
                          ])

##### Improve the metric names ####
all_df$metric <- gsub("\\.", " ", all_df$metric)

##### Rename multiple metrics for better plotting #####


all_df$metric <- gsub(" LL", ", lower level", all_df$metric)
all_df$metric <- gsub(" HL", ", higher level", all_df$metric)

all_df$metric <- firstup(all_df$metric)

# all_df$metric <- gsub(
#   "Alatalo interaction evenness",
#   "Alatalo\ninteraction\nevenness",
#   all_df$metric
# )
# 
# all_df$metric <- gsub(
#   "Cluster coefficient",
#   "Cluster\ncoefficient",
#   all_df$metric
# )
# 
# all_df$metric <- gsub(
#   "Compartment diversity",
#   "Compartment\ndiversity",
#   all_df$metric
# )
# 
# all_df$metric <- gsub(
#   "Extinction slope",
#   "Extinction\nslope",
#   all_df$metric
# )
# 
all_df$metric <- gsub(
  "H2",
  "H2'",
  all_df$metric
)

# all_df$metric <- gsub(
#   "Interaction evenness",
#   "Interaction\nevenness",
#   all_df$metric
# )
# 
# all_df$metric <- gsub(
#   "Interaction strength asymmetry",
#   "Interaction\nstrength\nasymmetry",
#   all_df$metric
# )
# 
# all_df$metric <- gsub(
#   "Linkage density",
#   "Linkage\ndensity",
#   all_df$metric
# )
# 
# all_df$metric <- gsub(
#   "Links per species",
#   "Links\nper species",
#   all_df$metric
# )
# 
# all_df$metric <- gsub(
#   "Mean number of shared partners",
#   "Mean number\nof shared\npartners",
#   all_df$metric
# )
# 
# all_df$metric <- gsub(
#   "Number of compartments",
#   "Number of\ncompartments",
#   all_df$metric
# )
# 
# all_df$metric <- gsub(
#   "Number of species",
#   "Number of\nspecies",
#   all_df$metric
# )
# 
# all_df$metric <- gsub(
#   "Partner diversity",
#   "Partner\ndiversity",
#   all_df$metric
# )
# 
# all_df$metric <- gsub(
#   "Shannon diversity",
#   "Shannon\ndiversity",
#   all_df$metric
# )
# 
# all_df$metric <- gsub(
#   "Specialisation asymmetry",
#   "Specialisation\nasymmetry",
#   all_df$metric
# )
# 
# all_df$metric <- gsub(
#   "Web asymmetry",
#   "Web\nasymmetry",
#   all_df$metric
# )
# 
# all_df$metric <- gsub(
#   "Weighted cluster coefficient",
#   "Weighted\ncluster\ncoefficient",
#   all_df$metric
# )
# 
# all_df$metric <- gsub(
#   "Weighted ",
#   "Weighted\n",
#   all_df$metric
# )




desired_metrics <- c(
  "Functional complementarity, higher level", "Modularity", "Nestedness", 
  "NODF", "Robustness, higher level", "H2'"
)

# Reorder the metrics so they're plotted alphabetically
desired_metrics <- desired_metrics[order(desired_metrics)]

# Reorder the network level, so that seed consumers are plotted above plants
all_df$net_level <- factor(all_df$net_level, levels = c('Seed consumer', 'Plant'))




# Function for making the basic plot --------------------------------------


tileplot_function <- function(df, metric, legend_desired = F, save_file = F) {
  tileplot <- ggplot(df, aes(x = net1_clust, y = net2_clust)) +
    geom_tile(aes(fill = percentage)) +
    scale_fill_viridis_c(
      option = "cividis", direction = -1,
      name = "Erroneous findings possible",
      limits = c(0, 100),
      breaks = c(0, 100),
      labels = c("0%", "100%") # Make the legend show values as a percentage
    ) +
    facet_grid(net_level ~ taxonomic_level) +
    theme(
      panel.background = element_blank(), # Get rid of the annoying background formatting
      strip.background = element_rect(fill = "white"), # Get rid of the grey facet labels
      legend.position = "bottom",
      strip.text = element_text(size = 11),
      axis.title.x = element_blank(),
      axis.title.y = element_blank()
    ) +
    ggtitle(bquote(~ underline(.(metric)))) # bquote and .( are needed to use the value of metric, rather than a string 'metric'


  if (legend_desired == F) {
    tileplot <- tileplot + theme(
      legend.position = "none"
    )
  }
  if(save_file == T){
    filename <- paste0("figures/individual_metrics/", metric, ".jpg")
    ggsave(filename, tileplot, 
         units = "in", height = 8,
         width = 16, dpi = 300
    )
  }

  return(tileplot)
}


###### Make a plot of a subset of desired metrics, across clumping types #####

# Make the desired plots in a simple for loop
plot_list <- list()
for (i in 1:length(desired_metrics)) {
  lev <- desired_metrics[i]
  df <- filter(all_df, metric == lev)

  plot_list[[lev]] <- tileplot_function(df, metric = lev, legend_desired = F)
}



# Function for extracting the legend
g_legend <- function(a.gplot) {
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}

# Extract the legend
plot_legend <- g_legend(tileplot_function(df, metric = lev, legend_desired = T))


# plot layout
lay <- rbind(c(1,2,3),
             c(4,5, 6),
             c(7,7,7),
             c(8, 8, 8))

# Save the plot
jpeg("figures/grid.jpg",
  units = "in", height = 8,
  width = 16, res = 300
)
grid_plot <- gridExtra::grid.arrange(plot_list[[1]], plot_list[[2]], plot_list[[3]],
  plot_list[[4]], plot_list[[5]], plot_list[[6]], 
  grid::textGrob('Proportion of Network 1 nodes collapsed'), 
  plot_legend,
  layout_matrix = lay,
  heights = unit(c(3.3,3.3,0.5, 1), c('in', 'in', 'in')),
  
  left = "Proportion of Network 2 nodes collapsed"
)
dev.off()



##### Plot each metric, for SI #####


# all_df$metric <- gsub(
#   "Functional complementarity",
#   "Functional\ncomplementarity",
#   all_df$metric
# ) # Changed here as newline would have been annoying earlier
# 
# all_df$metric <- gsub('higher', '\nhigher', all_df$metric)
# all_df$metric <- gsub('lower', '\nlower', all_df$metric)

tax_levels <- unique(all_df$taxonomic_level)
net_levels <- unique(all_df$net_level)



# Use the plot making function on all metrics
lapply(unique(all_df$metric), function(x)
  tileplot_function(filter(all_df, metric == x), metric = x, legend_desired = T,
                    save_file = T))
