## Load necessary libraries
library(ggplot2)
library(paletteer)
library(argparse)

## Parsers
parser <- ArgumentParser(description="Example R script using argparse")
parser$add_argument("-csv", "--csv_file", type="character", help="Csv file with strain info", required=TRUE)
args <- parser$parse_args()

csv_file <- args$csv_file
suffix <- sub(".*(_\\d{2}_\\d{2}_\\d{4}_\\d{2}_\\d{2}).*", "\\1", csv_file)

# Read the data from the CSV file
data <- read.csv(csv_file, sep = ",", header = TRUE)

## Plot the number of viral strains
plot.ViralStrains <- ggplot(data, aes(x = Time)) +
  geom_line(aes(y = Total.infected,
                color = "All strains")) +
  geom_line(aes(y = Super.spreaders,
                color = "Super strains")) +
  geom_line(aes(y = Normal.spreaders,
                color = "Normal strains")) +
  labs(title = "Number of viral strains over Time",
       x = "Simulation Time",
       y = "Number of strains") +
  scale_colour_paletteer_d("khroma::highcontrast") +
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
save_plot(plot.ViralStrains, paste0("plotStrains", suffix, ".png"))