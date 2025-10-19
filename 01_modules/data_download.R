# modules/data_download.R

# UI for Data Download Tab
data_download_ui <- function(id) {
  ns <- NS(id)
  
  tagList(
    card(
      card_header("Download Processed Data"),
      selectInput(ns("download_tax_level"), "Taxonomic Level for Download:",
                  choices = taxonomic_levels, selected = "Genus"),
      selectInput(ns("download_exp"), "Experiment(s) for Download:",
                  choices = c("All", experiments), multiple = TRUE,
                  selected = "All"),
      selectInput(ns("file_format"), "File Format:",
                  choices = c("CSV", "TSV", "Excel")),
      downloadButton(ns("download_data"), "Download Data")
    )
  )
}

# Server for Data Download Tab
data_download_server <- function(id, phyloseq_obj) {
  moduleServer(id, function(input, output, session) {
    
    # Download handler
    output$download_data <- downloadHandler(
      filename = function() {
        paste("microbiome_data_", input$download_tax_level, "_", 
              format(Sys.time(), "%Y%m%d_%H%M%S"), 
              switch(input$file_format,
                     "CSV" = ".csv",
                     "TSV" = ".tsv",
                     "Excel" = ".xlsx"), 
              sep = "")
      },
      content = function(file) {
        # Filter data based on download criteria
        if ("All" %in% input$download_exp) {
          download_ps <- phyloseq_obj
        } else {
          download_ps <- filter_phyloseq_by_experiment(phyloseq_obj, input$download_exp)
        }
        
        # Prepare data at specified taxonomic level
        final_data <- get_abundance_table(download_ps, input$download_tax_level)
        
        # Write file based on selected format
        switch(input$file_format,
               "CSV" = write.csv(final_data, file, row.names = FALSE),
               "TSV" = write.table(final_data, file, sep = "\t", row.names = FALSE),
               "Excel" = writexl::write_xlsx(final_data, file))
      }
    )
  })
}