## Load necessary libraries
library(ggplot2)
library(argparse)

## Parsers
parser <- ArgumentParser(description="Example R script using argparse")
parser$add_argument("-csv", "--csv_file", type="character", help="Csv file with statistics", required=TRUE)
parser$add_argument("-s", "--sample_size", type="integer", help="Sample size", required=TRUE)
parser$add_argument("-l", "--genome_length", type="integer", help="Genome length size", required=TRUE)
parser$add_argument("-pop", "--population_type", type="character", help="Type of population (mix, super strains, normal strains)", required=TRUE)
args <- parser$parse_args()

event_number <- function(csv_file) {
  # Use regular expression to find the number at the end of the file name before '.csv'
  match <- regmatches(csv_file, regexec("_(\\d+)\\.csv$", csv_file))
  
  if (length(match[[1]]) > 1) {
    number <- match[[1]][2] # Extract the captured number
  } else {
    print("No valid number found in the file name.")
    number <- NULL
  }
  
  return(number)
}

csv_file <- args$csv_file
sample_size <- args$sample_size
genome_length <- args$genome_length
suffix <- event_number(csv_file)

# Read the data from the CSV file
data <- read.csv(csv_file, sep = ",", header = TRUE)

plot.SelectiveSweep <- ggplot(data, aes(x = Window)) +
  geom_line(aes(y = Tajima.s.D,
                color = "Tajima's D")) +
  geom_line(aes(y = Pi,
                color = "Pi")) +
  geom_line(aes(y = Theta.W,
                color = "Theta Waterson")) +
  labs(title = "Selective Sweep",
       subtitle = paste("Sample =", sample_size, ", Genome length =", genome_length),
       x = "Window",
       y = "Value") +
  scale_color_manual(values = c("Theta Waterson" = "#084351", "Pi" = "#df6f46", "Tajima's D" = "#72874EFF")) +
  theme_minimal()

### Save plots ##
save_plot <- function(plot, prefix, suffix, extension = ".png", width = 8, height = 6, dpi = 300) {
  # Construct the full filename
  filename <- paste0(prefix, suffix, extension)
  
  # Apply white background to the plot
  plot <- plot +
    theme(panel.background = element_rect(fill = "white", color = NA),
          plot.background = element_rect(fill = "white", color = NA))
  
  # Save the plot as PNG
  ggsave(filename, plot = plot, width = width, height = height, dpi = dpi)
}
save_plot(plot.SelectiveSweep, "plotSelSw_", suffix)