# server.R
# server.R
source("01_modules/taxon_explorer.R")
source("01_modules/sample_explorer.R")
source("01_modules/experiment_explorer.R")

server <- function(input, output, session) {
  # Pass the phyloseq object to each module
  taxon_explorer_server("taxon", ps = ps)
  sample_explorer_server("sample", ps = ps)
  experiment_explorer_server("experiment", ps = ps)
}