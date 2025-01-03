## Load necessary libraries
library(ggplot2)
library(paletteer)
library(argparse)
library(dplyr)

## Parsers
parser <- ArgumentParser(description="Example R script for multiple files with infection rates in one plot")
parser$add_argument("-csvs", "--csv_files", nargs="+", type="character", help="CSV files with strain info followed by their infection rates", required=TRUE)
args <- parser$parse_args()

# Split arguments into file names and infection rates
csv_args <- args$csv_files
num_files <- length(csv_args) / 2
files <- csv_args[seq(1, length(csv_args), by = 2)]
infection_rates <- as.numeric(csv_args[seq(2, length(csv_args), by = 2)])

# Initialize an empty data frame to store combined data
combined_data <- data.frame()
last_points <- data.frame()  # To store last x-values for each infection rate

## Loop through each file and infection rate to read and label data
for (i in 1:num_files) {
  file <- files[i]
  infection_rate <- infection_rates[i]
  
  # Read data from each CSV file
  data <- read.csv(file, sep = ",", header = TRUE)
  
  # Add a column for infection rate as a factor (to be used as a label in the plot)
  data$InfectionRate <- as.factor(infection_rate)
  
  # Find the last time point for this file
  last_time <- max(data$Time)
  
  # Append the last time point to last_points data frame
  last_points <- bind_rows(last_points, data.frame(InfectionRate = as.factor(infection_rate), LastTime = last_time))
  
  # Select relevant columns and bind them to the combined data frame
  combined_data <- bind_rows(combined_data, data)
}

## Plot combined data with each line colored by infection rate
plot_combined <- ggplot(combined_data, aes(x = Time, y = Normal.spreaders, color = InfectionRate)) +
  geom_line() +
  labs(title = "Normal Strains Over Time by Infection Rate",
       x = "Simulation Time",
       y = "Number of Normal Strains",
       color = "Infection Rate") +
  scale_colour_paletteer_d("ggthemes::Classic_20") +
  theme_minimal()

# Add vertical lines for the last time point of each infection rate with matching colors
plot_combined <- plot_combined + 
  geom_vline(data = last_points, aes(xintercept = LastTime, color = InfectionRate), linetype = "dashed", show.legend = FALSE) 

## Save the combined plot as a PNG file
save_plot <- function(plot, filename, width = 8, height = 6, dpi = 300) {
  # Apply white background to the plot
  plot <- plot +
    theme(panel.background = element_rect(fill = "white", color = NA),
          plot.background = element_rect(fill = "white", color = NA))
  # Save the plot as PNG
  ggsave(filename, plot = plot, width = width, height = height, dpi = dpi)
}

# Save the combined plot
save_plot(plot_combined, "timeVSrate.png")