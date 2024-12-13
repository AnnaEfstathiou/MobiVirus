## Load necessary libraries
library(ggplot2)
library(paletteer)
library(argparse)

## Parsers
parser <- ArgumentParser(description="Example R script using argparse")
parser$add_argument("-csv", "--csv_file", type="character", help="Csv file with strain info", required=TRUE)
args <- parser$parse_args()

csv_file <- args$csv_file

# Read the data from the CSV file
data <- read.csv(csv_file, sep = ",", header = TRUE)

plot.NormalStrains <- ggplot(data, aes(x = Time, y = Normal.spreaders)) +
  geom_line(color = "#DDAA33FF") +
  labs(x = 'Simulation Time', y = 'Number of strains', title = 'Number of normal strains over Time') +
  theme_minimal()

## Generic function for saving plots as png files
save_plot <- function(plot, filename, width = 8, height = 6, dpi = 300) {
  # Apply white background to the plot
  plot <- plot +
    theme(panel.background = element_rect(fill = "white", color = NA),
          plot.background = element_rect(fill = "white", color = NA))
  # Save the plot as PNG
  ggsave(filename, plot = plot, width = width, height = height, dpi = dpi)
}
## Save plot as png files
save_plot(plot.NormalStrains, paste0("plotStrains", suffix, ".png"))