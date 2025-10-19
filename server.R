# server.R
source("global.R", local = TRUE)

server <- function(input, output, session) {
  # Call module servers
  individual_samples_server("individual_samples", ps)
  experiment_analysis_server("experiment_analysis", ps)
  data_download_server("data_download", ps)
}