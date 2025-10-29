# server.R
server <- function(input, output, session) {
  # Load your phyloseq object here
  # For example:
  # physeq_obj <- reactive({
  #   readRDS("data/your_phyloseq_object.rds")
  # })
  
  # For this example, I'll assume you have your physeq_obj available
  # You may want to make it reactive if loading from file
  
  # Call module servers
  callModule(sample_explorer_server, "sample_tab", physeq_obj)
  callModule(experiment_aggregator_server, "experiment_tab", physeq_obj)
  callModule(file_downloader_server, "download_tab", physeq_obj)
}