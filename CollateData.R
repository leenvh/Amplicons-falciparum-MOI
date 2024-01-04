require(dplyr)
library(tidyr)
library(ggplot2)
library(reshape2)
library(patchwork)
library(gridExtra)
library(stringr)


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
mosquito_TRAP <- TRAP %>% filter(species == "mosquito")%>%
  mutate(Percentage = as.numeric(Percentage),
         Day = factor(Day, levels = c(0, 2, 7, 14, 21, 28, 35)))

write.csv(human_TRAP, "human_TRAP.csv")
write.csv(mosquito_TRAP, "mosquito_TRAP.csv")

# Function to expand dataset to include all days
expand_dataset <- function(df, all_days) {
  df %>% 
    complete(Day = all_days, fill = list(Percentage = 0))
}


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


# Identify all unique days
all_days <- unique(c(human_TRAP$Day, mosquito_TRAP$Day))

# Expand both datasets to include all days
expanded_human_TRAP <- expand_dataset(human_TRAP, all_days)
expanded_mosquito_TRAP <- expand_dataset(mosquito_TRAP, all_days)

# Define a function to create a combined plot for each day
create_combined_day_plot <- function(human_df, mosquito_df, day) {
  # Create the human plot for the day
  human_plot <- create_stacked_plot(human_df %>% filter(Day == day))
  
  # Create the mosquito plot for the day
  mosquito_plot <- create_stacked_plot(mosquito_df %>% filter(Day == day))
  
  # Combine the mosquito and human plots for the day
  combined_day_plot <- mosquito_plot / human_plot
  
  return(combined_day_plot)
}

# Define the order of days
ordered_days <- c(0, 2, 7, 14, 21, 28,35)

# Loop through each individual and create combined plots
combined_TRAP_plots <- list()
for (individual in unique(expanded_human_TRAP$Individual)) {
  combined_plots_per_day <- list()
  
  for (day in ordered_days) {
    human_data_for_day <- expanded_human_TRAP %>% filter(Individual == individual, Day == day)
    mosquito_data_for_day <- expanded_mosquito_TRAP %>% filter(Individual == individual, Day == day)
    
    # Create the human plot for the day or a blank plot with the day label if no data, or a blank plot with day label and y axis if no data at day0
    human_plot <- if (nrow(human_data_for_day) > 0) {
      create_stacked_plot(human_data_for_day, day)
    } else {
      blank_plot <- ggplot(data.frame(SampleID = paste("Day", day)), aes(x = SampleID)) +
        geom_blank() +
        scale_x_discrete(limits = paste("Day", day)) +
        labs(x = "") +
        theme_minimal() +
        theme(
          axis.text.x = element_text(size = 8),
          axis.text.y = element_text(size = 8),
          axis.ticks.y = element_line(),
          axis.title.y = element_text(),
          panel.grid = element_blank(),
          plot.margin = margin(5.5, 5.5, 5.5, 5.5),
          axis.line.x = element_line(color="black"),
          axis.ticks.x = element_line(color="black")
        )
      
      # Add y-axis label for Day 0 only
      if (day == 0) {
        blank_plot <- blank_plot + labs(y = "Haplotype %") + 
          scale_y_continuous(breaks = c(0, 25, 50, 75, 100), limits = c(0, 100)) +
          theme(
            axis.line.y = element_line(color="black"),
            axis.ticks.y = element_line(color="black"),
            axis.text.y = element_text(size = 8),
          )
      } else {
        blank_plot <- blank_plot + theme(axis.title.y = element_blank())
      }
      
      blank_plot
    }
    
    # Create the mosquito plot for the day or a blank plot if no data
    mosquito_plot <- if (nrow(mosquito_data_for_day) > 0) {
      create_stacked_plot(mosquito_data_for_day, day, is_mosquito = TRUE)
    } else {
      blank_plot_mosq <- ggplot() +
        theme_classic() + 
        theme(
          axis.text.x = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.title.y = element_blank(),
          axis.title.x = element_blank(),
          axis.line.y = element_blank(),
          axis.line.x = element_blank(),
          panel.grid = element_blank()
        ) 
    # Add y-axis label for Day 0 only
      if (day == 0) {
        blank_plot_mosq <- blank_plot_mosq + labs(y = "Haplotype %") + 
          scale_y_continuous(breaks = c(0, 25, 50, 75, 100), limits = c(0, 100)) +
          theme(
            axis.line.y = element_line(color="black"),
            axis.ticks.y = element_line(color="black"),
            axis.text.y = element_text(size = 8),
          )
      } else {
        blank_plot_mosq <- blank_plot_mosq + theme(axis.title.y = element_blank())
      }
      
      blank_plot_mosq
    }
    
    # Combine the mosquito and human plots for the day if mosquito data is present
    combined_day_plot <- if (!is.null(mosquito_plot)) {
      mosquito_plot / human_plot
    } else {
      human_plot
    }
    
    combined_plots_per_day[[as.character(day)]] <- combined_day_plot
  }
  
  # Remove NULL values from the list
  combined_plots_per_day <- combined_plots_per_day[!sapply(combined_plots_per_day, is.null)]
  
  # Create an annotation plot for Mosquito and Human labels
  annotation_plot <- ggplot() +
    annotate("text", x = 0, y = 100, label = "Mosquito", vjust = -27, size = 5) +
    annotate("text", x = 0, y = 100, label = "Human", vjust = 5, size = 5) +
    theme_void() +
    theme(plot.margin = margin(0, 0, 0, 0))
  
  # Combine all daily plots into a single plot for the individual
  if (length(combined_plots_per_day) > 0) {
    combined_TRAP_plot <- wrap_plots(list(annotation_plot, wrap_plots(combined_plots_per_day, nrow = 1)), ncol = 2, widths = c(0.1, 1)) +
      plot_annotation(title = individual,theme = theme(
                        plot.title = element_text(hjust = 0.5, size = 14)))
    
    combined_TRAP_plots[[individual]] <- combined_TRAP_plot
  } else {
    # If there are no plots, create a placeholder for the individual
    combined_TRAP_plots[[individual]] <- ggplot() + theme_void() +
      labs(title = paste("No data for individual", individual))
  }
}


print(combined_TRAP_plots[["PQ-04-018"]])
print(combined_TRAP_plots[["PQ-04-035"]])


# Determine the layout of the grid
number_of_individuals <- length(combined_TRAP_plots)
number_of_columns <- 5 # or however many columns you want
number_of_rows <- ceiling(number_of_individuals / number_of_columns)

# Create the large plot
large_plot <- wrap_plots(combined_TRAP_plots, ncol = number_of_columns, nrow = number_of_rows)
ggsave("large_combined_plot_TRAP.pdf", large_plot, width = 49, height = 49, units = "in")

