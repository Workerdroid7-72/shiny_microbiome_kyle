# modules/sample_explorer.R
sample_explorer_ui <- function(id) {
  ns <- NS(id)
  
  tagList(
    fluidRow(
      column(3,
             card(
               title = "Filters",
               selectInput(
                 ns("tax_level"),
                 "Taxonomic Level:",
                 choices = TAX_LEVELS,
                 selected = "Genus"
               ),
               selectInput(
                 ns("sample_id"),
                 "Select Sample:",
                 choices = NULL  # Will be populated by server
               ),
               sliderInput(
                 ns("top_n"),
                 "Show Top N taxa:",
                 min = 5, max = 50, value = 10
               )
             )
      ),
      column(9,
             card(
               title = "Sample Abundance Data",
               DT::dataTableOutput(ns("abundance_table")),
               downloadButton(ns("download_table"), "Download Table")
             )
      )
    ),
    card(
      title = "Abundance Visualization",
      plotlyOutput(ns("abundance_plot"))
    )
  )
}

sample_explorer_server <- function(input, output, session, physeq_obj) {
  ns <- session$ns
  
  # Reactive: Get available samples
  available_samples <- reactive({
    # Assuming physeq_obj is available
    sample_names(physeq_obj())
  })
  
  # Update sample selector
  observe({
    updateSelectInput(
      session,
      "sample_id",
      choices = available_samples(),
      selected = if(length(available_samples()) > 0) available_samples()[1] else NULL
    )
  })
  
  # Reactive: Get abundance data for selected sample
  sample_abundance_data <- reactive({
    req(input$sample_id, input$tax_level, physeq_obj())
    
    # Extract abundance data for selected sample
    otu_mat <- otu_table(physeq_obj())
    tax_mat <- tax_table(physeq_obj())
    
    # Get data for selected sample
    sample_data <- as.numeric(otu_mat[, input$sample_id])
    
    # Create data frame with taxonomy
    df <- data.frame(
      Taxon = rownames(tax_mat),
      Abundance = sample_data,
      stringsAsFactors = FALSE
    )
    
    # Add taxonomic level
    df$Level <- tax_mat[, input$tax_level]
    
    # Aggregate by taxonomic level
    result <- df %>%
      group_by(Level) %>%
      summarise(Abundance = sum(Abundance), .groups = 'drop') %>%
      filter(Level != "", !is.na(Level)) %>%
      arrange(desc(Abundance))
    
    # Limit to top N
    result <- head(result, input$top_n)
    
    return(result)
  })
  
  # Output: Abundance table
  output$abundance_table <- DT::renderDataTable({
    req(sample_abundance_data())
    
    sample_abundance_data() %>%
      DT::datatable(
        options = list(
          pageLength = 15,
          lengthMenu = list(c(10, 15, 25, -1), c(10, 15, 25, "All"))
        )
      )
  })
  
  # Output: Abundance plot
  output$abundance_plot <- renderPlotly({
    req(sample_abundance_data())
    
    p <- ggplot(sample_abundance_data(), aes(x = reorder(Level, Abundance), y = Abundance)) +
      geom_col(fill = "#3498db", alpha = 0.7) +
      coord_flip() +
      labs(
        title = paste("Abundance at", input$tax_level, "Level"),
        x = input$tax_level,
        y = "Abundance"
      ) +
      theme_minimal()
    
    ggplotly(p, tooltip = c("x", "y"))
  })
  
  # Download handler
  output$download_table <- downloadHandler(
    filename = function() {
      paste("sample_abundance_", input$sample_id, "_", Sys.Date(), ".csv", sep = "")
    },
    content = function(file) {
      write.csv(sample_abundance_data(), file, row.names = FALSE)
    }
  )
}