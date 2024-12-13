# Load necessary libraries
library(ggplot2)
library(viridis)
library(argparse)

## Parsers
parser <- ArgumentParser(description="Example R script using argparse")
parser$add_argument("-csv", "--csv_file", type="character", help="Csv file with statistics", required=TRUE)
parser$add_argument("-l", "--genome_length", type="integer", help="Genome length size", required=TRUE)
parser$add_argument("-png", "--save_png", action="store", help="Save the PNG file")
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
genome_length <- args$genome_length
suffix <- event_number(csv_file)

ld_data <- read.csv(csv_file, sep = ",", header = TRUE)

ld_data$i_scaled <- as.numeric(as.factor(ld_data$i * genome_length))
ld_data$j_scaled <- as.numeric(as.factor(ld_data$j * genome_length))

# LD Plot for rsq
plot.RSQ <- ggplot(ld_data, aes(x = i_scaled, y = j_scaled, fill = rsq)) +
  geom_tile() +
  labs(title = "LD Plot (rsq)", x = "pos i", y = "pos j") +
  scale_fill_viridis(discrete = FALSE) +
  theme_minimal()

# LD Plot for D
plot.D <- ggplot(ld_data, aes(x = i_scaled, y = j_scaled, fill = D)) +
  geom_tile() +
  labs(title = "LD Plot (D)", x = "pos i", y = "pos j") +
  scale_fill_viridis(discrete = FALSE) +
  theme_minimal()

# LD Plot for Dprime
plot.Dprime <- ggplot(ld_data, aes(x = i_scaled, y = j_scaled, fill = Dprime)) +
  geom_tile() +
  labs(title = "LD Plot (D')", x = "pos i", y = "pos j") +
  scale_fill_viridis(discrete = FALSE) +
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

## Saving plots as png files
if (!is.null(args$save_png)) {
  save_plot(plot.RSQ, paste0("plotLD_rsq_", suffix, ".png"))
  save_plot(plot.D, paste0("plotTLD_d_", suffix, ".png"))
  save_plot(plot.Dprime, paste0("plotLD_dprime_", suffix, ".png"))
} else {
  # Save all plots in a single PDF file
  pdf("PlotsLD.pdf", width = 8, height = 6)
  print(plot.RSQ)
  print(plot.D)
  print(plot.Dprime)
  dev.off()
}