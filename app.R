# 1. Load packages (usually from global.R)
source("global.R")
source("01_modules/taxon_explorer.R")
source("01_modules/sample_explorer.R")
source("01_modules/experiment_explorer.R")
source("01_modules/sample_compare.R")


ui <- page_navbar(
  title = "Microbiome Explorer",
  theme = bs_theme(bootswatch = "flatly"),
  navbar_options = list(class = "bg-primary", theme = "dark"),
  
  nav_panel("Taxon Explorer", taxon_explorer_ui("taxon")),
  nav_panel("Sample Explorer", sample_explorer_ui("sample")),
  nav_panel("Experiment Explorer", experiment_explorer_ui("experiment")),
  nav_panel("Before After Explorer", sample_compare_ui("compare")),
  
  # Inject custom CSS
  tags$head(
    tags$style(
      HTML("
        /* Add space between app title and tab links */
        .navbar-brand + .navbar-nav {
          margin-left: 3rem !important;
        }
      ")
    )
  )
  
  
)

server <- function(input, output, session) {
  # Pass the phyloseq object to each module
  taxon_explorer_server("taxon", ps = ps)
  sample_explorer_server("sample", ps = ps)
  experiment_explorer_server("experiment", ps = ps)
  sample_compare_server("compare", ps = ps)
}


shinyApp(ui = ui, server = server)
