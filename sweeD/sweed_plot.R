## Load necessary libraries
library(ggplot2)
library(paletteer)
library(argparse)

## Parsers
parser <- ArgumentParser(description="Example R script using argparse")
parser$add_argument("-f", "--file", type="character", help="Csv file with strain info", required=TRUE)
parser$add_argument("-s", "--suffix", type="character", help="Suffix for file naming", required=TRUE)
args <- parser$parse_args()

data_file <- args$file
suffix <- args$suffix

# Read the data from the CSV file
data <- read.table(data_file, sep = "\t", header = TRUE, skip = 2) # skip the first 2 lines because the table starts in line 3

# Use ggplot directly on the dataframe without creating a separate variable
plot.SelectiveSweep <- ggplot(data, aes(x = Position, y = Likelihood)) +
  geom_line(color = "#492406") +
  labs(x = 'Positions', y = 'Likelihood', title = 'Selective Sweep') +
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
save_plot(plot.SelectiveSweep, paste0("plotSelSw_", suffix, ".png"))