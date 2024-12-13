## Load necessary libraries
library(ggplot2)
library(paletteer)
library(argparse)

## Parsers
parser <- ArgumentParser(description="Example R script using argparse")
parser$add_argument("-csv", "--csv_file", type="character", help="Csv file with coordinated info", required=TRUE)
parser$add_argument("-inf", "--infectors_file", type="character", help="Csv file with the infector of each 'infector' event", required=TRUE)
args <- parser$parse_args()

# Read the data from the CSV file
coords_data <- read.csv(args$csv_file, sep = ",", header = TRUE)
infectors_data <- read.csv(args$infectors_file, sep = ",", header = TRUE)

# Extracting the number from the CSV filename
base_name <- basename(args$csv_file)       # Get the filename only (without the path)
name_part <- tools::file_path_sans_ext(base_name)  # Remove the file extension
number_part <- unlist(strsplit(name_part, "_"))[length(unlist(strsplit(name_part, "_")))]  # Get the number part

## Function to plot scatter plot from CSV files
plot_coordinates <- function(coords_data, infectors_data) {

  # Assigning labels for every mutation value
  mutation_labels <- c('Healthy individuals','Normal spreaders','Super spreaders')
  mutation_colors <- c('#FD7901FF','#0E7175FF','#C35BCAFF')

  # Create a data frame with mutation and corresponding labels and colors
  coords_data$Mutation_label <- factor(coords_data$Mutation, 
                                       levels = c(0.0, 1.0, 2.0), 
                                       labels = mutation_labels)

  # Initialize an Event.Type column for coords_data with 'Individuals' for all individuals
  coords_data$Event.Type <- 'Individuals'

  # Get the row corresponding to the current event in infectors_data
  event_row <- infectors_data[infectors_data$Event == as.integer(number_part), ]

  # Check if the event is 'infection' and get the infector
  if (event_row$Event.Type[1] == 'infection') {
    infector <- event_row$Individual[1]
    # Assign 'Infector' label to the respective individual
    coords_data$Event.Type[infector + 1] <- 'Infector'  # R is 1-indexed
  }

  # Count occurrences for each mutation type
  mutation_counts <- table(coords_data$Mutation_label)

  # Modify mutation labels to include counts
  mutation_labels_with_count <- paste0(mutation_labels, " (", mutation_counts[mutation_labels], ")")

  # Convert Event.Type into a factor to use in the legend
  coords_data$Event.Type <- factor(coords_data$Event.Type, 
                                   levels = c('Individuals', 'Infector'))

  # Create the base plot with ggplot2
  ScatterPlot <- ggplot(coords_data, aes(x = x, y = y)) +
    geom_point(aes(color = Mutation_label, shape = Event.Type), size = 1.5) +
    scale_color_manual(name = 'Colors', 
                       values = mutation_colors,
                       labels = mutation_labels_with_count) +
    scale_shape_manual(name = 'Shapes', values = c(16, 8)) +  # 16 for 'Individuals', 8 for 'Infector'
    labs(title = 'Scatter Plot of xy coordinates', x = 'x', y = 'y') +
    theme_minimal() +
    guides(color = guide_legend(order = 1), 
           shape = guide_legend(order = 2))  # Control order of legends

  return(ScatterPlot)
}


## Generic function for saving plots as png files
save_plot <- function(plot, filename, width = 8, height = 6, dpi = 300) {
  # Apply white background to the plot
  plot <- plot +
    theme(panel.background = element_rect(fill = "white", color = NA),
          plot.background = element_rect(fill = "white", color = NA))
  # Save the plot as PNG
  ggsave(filename, plot = plot, width = width, height = height, dpi = dpi)
}

# Example usage
plot_coords <- plot_coordinates(coords_data, infectors_data)
## Save plot as png file
scatter_plot_name <- paste0('plotCoords_', number_part, '.png')
save_plot(plot_coords, scatter_plot_name)
