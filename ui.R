# ui.R
source("01_modules/sample_explorer.R")
source("01_modules/experiment_aggregator.R")
source("01_modules/file_downloader.R")

ui <- page_navbar(
  title = "Microbiome Data Explorer",
  theme = bslib::bs_theme(
    version = 5,
    bootswatch = "flatly",  # Choose your preferred theme
    primary = "#2c3e50",
    secondary = "#34495e"
  ),
  
  nav_panel(
    "Sample Explorer",
    sample_explorer_ui("sample_tab")
  ),
  
  nav_panel(
    "Experiment Aggregator", 
    experiment_aggregator_ui("experiment_tab")
  ),
  
  nav_panel(
    "Download Data",
    file_downloader_ui("download_tab")
  ),
  
  # Optional: Add footer or additional info
  footer = div(
    class = "text-center p-3",
    "Microbiome Data Explorer | Built with R Shiny",
    br(),
    "© 2025 Keith Goddard"
  )
)