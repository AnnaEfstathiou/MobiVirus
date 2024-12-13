library(ggplot2)
library(argparse)

## Parsers
parser <- ArgumentParser(description = "Example R script using argparse")
parser$add_argument("-csv", "--csv_file", type="character", help="Csv file with statistics", required=TRUE)
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
suffix <- event_number(csv_file)
sfs_data <- read.csv(csv_file, sep = ",", header = TRUE)

plot.SFS <- ggplot(sfs_data, aes(x = factor(SNP.Count), y = Frequency)) +
  geom_bar(stat = "identity", fill = "#204035FF", width = 0.2) +
  xlab("SNP Count") +
  ylab("Frequency") +
  ggtitle("Site Frequency Spectrum") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 8)) +
  scale_x_discrete(breaks = seq(0, max(sfs_data$SNP.Count), by = 5))  

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
save_plot(plot.SFS, "plotSFS_", suffix)