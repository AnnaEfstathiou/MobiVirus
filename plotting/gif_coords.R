library(ggplot2)
library(dplyr)
library(gganimate)
library(gifski)



# Function to load coordinates files
load_coords_files <- function(directory) {

  # List all files that match the pattern coords_*.csv
  files <- list.files(path = directory, pattern = "^coords_[0-9]+\\.csv$", full.names = TRUE)

  # Initialize an empty list to store data frames
  data_list <- list()

  # Loop over each file
  for (file in files) {

    # Extract the frame number from the file name using regex
    event_number <- as.numeric(sub("coords_([0-9]+)\\.csv", "\\1", basename(file)))

    # Read the CSV file
    df <- read.csv(file)

    # Add the frame number as a new column
    df$event <- event_number

    # Append the data frame to the list
    data_list[[length(data_list) + 1]] <- df
  }

  # Combine all data frames into one
  coords_data <- bind_rows(data_list)

  # Ensure the data is sorted by the frame column
  coords_data <- arrange(coords_data, event)

  return(coords_data)
}



# Load the coordinate data from files
coords_data <- load_coords_files("/home/anna/mobivirus/simulation/simulation_10_10_2024_12_07/samples")



# Define the labels and colors for mutations
mutation_labels <- c('Healthy individuals', 'Normal spreaders', 'Super spreaders')
mutation_colors <- c('#1a001c', '#FD7901FF', '#0E7175FF')

# Add the mutation labels without counting
coords_data$Mutation_label <- factor(coords_data$Mutation,
                                     levels = c(0, 1, 2),
                                     labels = mutation_labels)

# Plot with gganimate
Plot.movement <- ggplot(coords_data, aes(x = x, y = y, color = Mutation_label)) +
  geom_point() +
  theme_bw() +
  labs(title = 'Event: {closest_state}', x = 'x', y = 'y') +
  scale_color_manual(name = 'Colors', values = mutation_colors) +
  transition_states(event, transition_length = 20, state_length = 10) +
  ease_aes('linear')

# Animate the plot and make it slower
anim <- animate(Plot.movement, nframes = 100, fps = 20, width = 800, height = 600, renderer = gifski_renderer())

# Save the animated GIF to a file
anim_save("plot_simulation.gif", animation = anim)