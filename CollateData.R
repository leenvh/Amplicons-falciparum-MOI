require(dplyr)
library(tidyr)
library(ggplot2)
library(reshape2)
library(patchwork)
library(gridExtra)
library(stringr)
library(cowplot)
library(future.apply)
library(furrr)


setwd("~/Google Drive/My Drive/PhD_LSHTM/Projects/AmpliconSequencing/Collate")
TRAP <- read.csv("~/Google Drive/My Drive/PhD_LSHTM/Projects/AmpliconSequencing/Haplotypr_Analysis/Matching_Haplotypes_TRAP_Molten.csv",header=TRUE,sep=',')
CSP <- read.csv("~/Google Drive/My Drive/PhD_LSHTM/Projects/AmpliconSequencing/Haplotypr_Analysis/Matching_Haplotypes_CSP_Molten.csv",header=TRUE,sep=',')
DrugRes <- read.csv("~/Google Drive/My Drive/PhD_LSHTM/Projects/AmpliconSequencing/MinION/merged_data.csv",header=TRUE,sep=',')


#Data cleaning and formatting
TRAP<-TRAP[,-c(1,13)]
CSP<-CSP[,-c(1,13)]
DrugRes<-DrugRes[,-1]
colnames(DrugRes)[2]<-"Individual"

human_TRAP <- TRAP %>% 
  filter(species == "human") %>%
  mutate(Percentage = as.numeric(Percentage),
         Day = factor(Day, levels = c(0, 2, 7, 14, 21, 28, 35)))
human_CSP <- CSP %>% 
  filter(species == "human") %>%
  mutate(Percentage = as.numeric(Percentage),
         Day = factor(Day, levels = c(0, 2, 7, 14, 21, 28, 35)))
mosquito_TRAP <- TRAP %>% filter(species == "mosquito")%>%
  mutate(Percentage = as.numeric(Percentage),
         Day = factor(Day, levels = c(0, 2, 7, 14, 21, 28, 35)))
mosquito_CSP <- CSP %>% filter(species == "mosquito")%>%
  mutate(Percentage = as.numeric(Percentage),
         Day = factor(Day, levels = c(0, 2, 7, 14, 21, 28, 35)))


# Function to create stacked plot
create_stacked_plot <- function(df, day, is_mosquito = FALSE) {
  # Extract the string after the last underscore for SampleID
  df$ShortSampleID <- str_extract(df$SampleID, "(?<=_)[^_]+$")
  
  # Determine the width of the bars based on whether it's a mosquito plot
  # and how many mosquito samples there are for that day
  bar_width <- if (is_mosquito) {
    if (length(unique(df$SampleID)) > 1) {
      (1 / length(unique(df$SampleID)) * 2)  # Adjusted width for multiple mosquito samples
    } else {
      1  # Full width if only one mosquito sample
    }
  } else {
    1  # Full width for human samples
  }
  
  # Create the plot using the shortened SampleID
  plot <- ggplot(df, aes(x = ShortSampleID, y = Percentage, fill = Haplotype)) +
    geom_bar(stat = "identity", position = "stack", width = bar_width) +
    theme_classic() +
    theme(legend.position = "none",axis.title.x = element_blank())
  
  # Conditionally remove the y-axis elements if not Day 0 human
  if (day != 0) {
    plot <- plot + 
      theme(axis.text.y = element_blank(),
            axis.ticks.y = element_blank(),
            axis.title.y = element_blank(),
            axis.line.y = element_blank())
  }
  
  # Conditionally add the y-axis title if Day 0 
  if (day == 0) {
    plot <- plot + 
      labs(y = "Haplotype %")
  }
  
  # Remove the x-axis line and annotations for mosquito plots
  if (is_mosquito) {
    plot <- plot + 
      theme(axis.text.x = element_blank(),
            axis.ticks.x = element_blank(),
            axis.title.x = element_blank(),
            axis.line.x = element_blank())
  }
  
  return(plot)
}




# Define a function to create a combined plot for each day
create_combined_day_plot <- function(human_df, mosquito_df, day) {
  # Create the human plot for the day
  human_plot <- create_stacked_plot(human_df %>% filter(Day == day), day, is_mosquito = FALSE)
  
  # Create the mosquito plot for the day
  mosquito_plot <- create_stacked_plot(mosquito_df %>% filter(Day == day), day, is_mosquito = TRUE)
  
  # Combine the mosquito and human plots for the day
  combined_day_plot <- mosquito_plot / human_plot
  
  return(combined_day_plot)
}


# Function to expand dataset to include all days
expand_dataset <- function(df, all_days) {
  df %>% 
    complete(Day = all_days, fill = list(Percentage = 0))
}


#Function to Process Individual
process_individual <- function(human_df, mosquito_df, ordered_days, individual) {
  combined_plots_per_day <- list()
  
  for (day in ordered_days) {
    human_data_for_day <- human_df %>% filter(Individual == individual, Day == day)
    mosquito_data_for_day <- mosquito_df %>% filter(Individual == individual, Day == day)
    
    # Call the create_combined_day_plot function
    combined_day_plot <- create_combined_day_plot(human_data_for_day, mosquito_data_for_day, day)
    combined_plots_per_day[[as.character(day)]] <- combined_day_plot
  }
  
  # Remove NULL values from the list
  combined_plots_per_day <- combined_plots_per_day[!sapply(combined_plots_per_day, is.null)]
  
  # Combine all daily plots into a single plot for the individual
  combined_plot <- wrap_plots(combined_plots_per_day, nrow = 1)
  print(paste("Processing individual:", individual))
  return(combined_plot)
}

#Function to Process all Individual
process_all_individuals <- function(human_data, mosquito_data, ordered_days) {
  all_individuals <- unique(c(human_data$Individual, mosquito_data$Individual))
  all_days <- unique(c(human_data$Day, mosquito_data$Day))
  
  # Expand both datasets to include all days
  expanded_human <- expand_dataset(human_data, all_days)
  expanded_mosquito <- expand_dataset(mosquito_data, all_days)
  
  combined_plots <- list()
  for (individual in all_individuals) {
    combined_plots[[individual]] <- process_individual(expanded_human, expanded_mosquito, ordered_days, individual)
  }
  
  return(combined_plots)
}


#Executing the Functions for TRAP and CSP Data
# Assuming you have loaded your TRAP and CSP data into human_TRAP, mosquito_TRAP, human_CSP, and mosquito_CSP
ordered_days <- c(0, 2, 7, 14, 21, 28, 35)

# Process TRAP data
combined_TRAP_plots <- process_all_individuals(human_TRAP, mosquito_TRAP, ordered_days)
print(combined_TRAP_plots[["PQ-04-018"]])
# Process CSP data
combined_CSP_plots <- process_all_individuals(human_CSP, mosquito_CSP, ordered_days)
print(combined_CSP_plots[["PQ-04-003"]])


# Set up parallel processing
plan(multisession, workers = 4)  # Choose the number of workers based on your CPU

# Function to process individual plots
process_individual_plots <- function(individual) {
  trap_plot <- combined_TRAP_plots[[individual]]
  csp_plot <- combined_CSP_plots[[individual]]
  
  # Create an annotation plot for Mosquito and Human labels
  annotation_plot <- ggplot() +
    annotate("text", x = 0, y = 100, label = "Mosquito", vjust = -27, size = 5) +
    annotate("text", x = 0, y = 100, label = "Human", vjust = 5, size = 5) +
    theme_void() +
    theme(plot.margin = margin(0, 0, 0, 0))
  
  # Combine the TRAP and CSP plots with equal widths
  combined_plot <- annotation_plot  + csp_plot + trap_plot +
    plot_layout(ncol = 3, widths = c(1,8,8))
  
  combined_plot_annotated <- ggdraw() +
    draw_plot(combined_plot) +
    draw_label(individual, size = 20, fontface = 'bold',
               x = 0.5, y = 1, hjust = 0.5, vjust = 1)+
    draw_label("TRAP", size = 15, x = 0.1, y = 1, hjust = 0.1, vjust = 1.5)+
    draw_label("CSP", size = 15, x = 0.58, y = 1, hjust = 0.58, vjust = 1.5)
  
  return(combined_plot_annotated)
}

# List of individuals
individuals <- names(combined_TRAP_plots)

# Process all individuals in parallel
combined_plots_list <- future_map(individuals, process_individual_plots)
names(combined_plots_list) <- individuals

print(combined_plots_list[["PQ-04-058"]])