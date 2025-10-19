# ui.R
source("global.R", local = TRUE)

# Define the theme
my_theme <- bs_theme(
  version = 5,
  bg = "#FFFFFF",
  fg = "#000000",
  primary = "#2C3E50",
  secondary = "#18BC9C"
)

ui <- page_navbar(
  theme = my_theme,
  title = "Microbiome Data Explorer",
  
  # Tab 1: Individual Samples
  nav_panel(
    "Individual Samples",
    individual_samples_ui("individual_samples")
  ),
  
  # Tab 2: Experiment Analysis
  nav_panel(
    "Experiment Analysis", 
    experiment_analysis_ui("experiment_analysis")
  ),
  
  # Tab 3: File Downloads
  nav_panel(
    "Data Download",
    data_download_ui("data_download")
  )
)