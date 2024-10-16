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
parser$add_argument("-pop", "--population_type", type="character", help="Type of population (mix, super strains, normal strains)", required=TRUE)
parser$add_argument("-png", "--save_png", action="store", help="Save the PNG file")
parser$add_argument("-pdf", "--save_pdf", action="store", help="Save the PNG file")
args <- parser$parse_args()

csv_file <- args$csv_file
sample_size <- args$sample_size
genome_length <- args$genome_length

# Read the data from the CSV file
data <- read.csv(csv_file, sep = ",", header = TRUE)

################################################
################### Plotting ###################
################################################

## Generic function for creating single line plots
line_plot <- function(data, x_var, y_var, title, subtitle, x_label, y_label, line_color) {
  ggplot(data, aes({{ x_var }}, {{ y_var }})) +
    geom_line(color = line_color) +
    labs(title = title,
         subtitle = subtitle,
         x = x_label,
         y = y_label) +
    theme_minimal()
}

## Tajima's D plot
plot.TajimasD <- line_plot(data = data, 
                                  x_var = Time,
                                  y_var = TajimasD,
                                  title = "Tajima's D over Time",
                                  subtitle = paste("Sample =", sample_size, ", Genome length =", genome_length),
                                  x_label = "Simulation Time",
                                  y_label = "Tajima's D",
                                  line_color = "#C35BCAFF")

## Haplotype Diversity plot
plot.HaplotypeDiv <- line_plot(data = data,
                                      x_var = Time,
                                      y_var = Hp_Diversity,
                                      title = "Haplotype Diversity over Time",
                                      subtitle = paste("Sample =", sample_size, ", Genome length =", genome_length), 
                                      x_label = "Simulation Time",
                                      y_label = "Haplotype Diversity",
                                      line_color = "#004488FF")

## Unique haplotypes plot
plot.UniqueSeqs <- line_plot(data = data,
                                    x_var = Time,
                                    y_var = Unique_seqs,
                                    title = "Number of unique haplotypes over Time", 
                                    subtitle = paste("Sample =", sample_size, ", Genome length =", genome_length), 
                                    x_label = "Simulation Time",
                                    y_label = "# unique haplotypes",
                                    line_color = "#BB5566FF")

## Pi & Theta Waterson Estimators Plot
plot.ThetawPi <- ggplot(data, aes(x = Time)) +
  geom_line(aes(y = Pi, color = "Pi")) +
  geom_line(aes(y = Theta_w, color = "Theta Waterson")) +
  labs(title = "Pi & Theta Waterson Estimators over Time",
       subtitle = paste("Sample =", sample_size, ", Genome length =", genome_length),
       x = "Simulation Time",
       y = "Values") +
  scale_color_manual(values = c("Pi" = "#0E7175FF", "Theta Waterson" = "#FD7901FF")) +
  theme_minimal()

## Fst Plot
plot.Fst <- ggplot(data, aes(x = Time)) +
  geom_line(aes(y = Fst_hsm,
                color = "HSM Fst")) +
  geom_line(aes(y = Fst_slatkin,
                color = "Slatkin Fst")) +
  geom_line(aes(y = Fst_hbk,
                color = "HBK Fst")) +
  labs(title = "Fst over Time",
       subtitle = paste("Sample =", sample_size, ", Genome length =", genome_length),
       x = "Simulation Time",
       y = "Value") +
  scale_colour_paletteer_d("ltc::trio3") +
  theme_minimal()


## Save plots ##

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
  save_plot(plot.TajimasD, "plotTajimasD.png")
  save_plot(plot.ThetawPi, "plotThetawPi.png")
  save_plot(plot.HaplotypeDiv, "plotHaplotypeDiv.png")
  save_plot(plot.UniqueSeqs, "plotUniqueSeqs.png")
  save_plot(plot.Fst, "plotFst.png")
}

## Saving plots in pdf file
if (!is.null(args$save_pdf)) {
  pdf("PlotStatistics.pdf", width = 8, height = 10)
  grid.arrange(plot.TajimasD, plot.ThetawPi,
              ncol = 1, nrow = 2,
              padding = unit(1, "lines"))
  grid.arrange(plot.HaplotypeDiv, plot.UniqueSeqs,
              ncol = 1, nrow = 2,
              padding = unit(1, "lines"))
  grid.arrange(plot.Fst,
              ncol = 1, nrow = 1,
              padding = unit(1, "lines"))
  dev.off()
}