## Load necessary libraries
library(ggplot2)
library(gridExtra)
library(paletteer)
library(argparse)
library(rlang)

## Parsers
parser <- ArgumentParser(description="Example R script using argparse")
parser$add_argument("-csv", "--csv_file", type="character", help="Csv file with statistics", required=TRUE)
parser$add_argument("-s", "--sample_size", type="integer", help="Sample size", required=TRUE)
parser$add_argument("-l", "--genome_length", type="integer", help="Genome length size", required=TRUE)
parser$add_argument("-m", "--mutation_rate", type="double", help="Mutation rate", required=TRUE)
parser$add_argument("-pop", "--population_type", type="character", help="Type of population (mix, super strains, normal strains)", required=TRUE)
args <- parser$parse_args()

csv_file <- args$csv_file
sample_size <- args$sample_size
genome_length <- args$genome_length
mutation_rate <- args$mutation_rate
suffix <- sub(".*(_\\d{2}_\\d{2}_\\d{4}_\\d{2}_\\d{2}).*", "\\1", csv_file)

data <- read.csv(csv_file, sep = ",", header = TRUE) # Read the data from the CSV file

################### Plotting ###################

line_plot <- function(data, x_var, y_var, title, subtitle, x_label, y_label, 
                      line_color = "black", line_size = 1, line_type = "solid",
                      smooth_line = FALSE, smooth_color = "red", smooth_size = 1, smooth_type = "solid") 
                      {
  p <- ggplot(data, aes({{ x_var }}, {{ y_var }})) +
    geom_line(color = line_color, size = line_size, linetype = line_type) +  
    labs(title = title,
         subtitle = subtitle,
         x = x_label,
         y = y_label) +
    theme_minimal()
  
  if (smooth_line) {
    p <- p + geom_smooth(method = "loess", color = smooth_color, size = smooth_size, 
                         linetype = smooth_type, se = FALSE)}

  return(p)
}

plot.Ne <- line_plot(data = data, 
                     x_var = Time,
                     y_var = Theta_w/(2*mutation_rate*genome_length),
                     title = "Approximate effective population size",
                     subtitle = paste("Sample =", sample_size, ", Genome length =", genome_length, ", Mutation rate =", mutation_rate),
                     x_label = "Simulation Time",
                     y_label = "Ne",
                     line_color = "#04225CFF",
                     line_size = 0.5,
                     line_type = "solid",  
                     smooth_line = TRUE,
                     smooth_color = "#DC3B34FF",
                     smooth_size = 1.0,
                     smooth_type = "solid")  

## Generic function for saving plots as png files
save_plot <- function(plot, filename, width = 8, height = 6, dpi = 300) {
  # Apply white background to the plot
  plot <- plot +
    theme(panel.background = element_rect(fill = "white", color = NA),
          plot.background = element_rect(fill = "white", color = NA))
  # Save the plot as PNG
  ggsave(filename, plot = plot, width = width, height = height, dpi = dpi)
}

save_plot(plot.Ne, paste0("plotNe", suffix, ".png"))