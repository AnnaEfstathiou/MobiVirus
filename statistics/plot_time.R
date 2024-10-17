library(ggplot2)
library(argparse)


## Parsers
parser <- ArgumentParser(description="Example R script using argparse")
parser$add_argument("-csv", "--csv_file", type="character", help="Csv file with statistics", required=TRUE)
parser$add_argument("-png", "--save_png", action="store", help="Save the PNG file")
args <- parser$parse_args()

csv_file <- args$csv_file

time_related_data <- read.csv(csv_file, colClasses = c("numeric", "numeric"))

# Use ggplot directly on the dataframe without creating a separate variable
plot.TimeEvents <- ggplot(time_related_data, aes(x = Time, y = Event)) +
  geom_line(color = "#492406") +
  labs(x = 'Simulation Time', y = 'Events', title = 'Events over Simulation Time') +
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
save_plot(plot.TimeEvents, "plotTimeEvents.png")