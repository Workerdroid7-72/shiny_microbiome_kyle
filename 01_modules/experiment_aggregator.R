# modules/experiment_aggregator.R
experiment_aggregator_ui <- function(id) {
  ns <- NS(id)
  
  tagList(
    fluidRow(
      column(3,
             card(
               title = "Filters",
               selectInput(
                 ns("tax_level_exp"),
                 "Taxonomic Level:",
                 choices = TAX_LEVELS,
                 selected = "Genus"
               ),
               selectInput(
                 ns("experiment_var"),
                 "Group by Experiment Variable:",
                 choices = NULL  # Will be populated by server
               ),
               checkboxInput(
                 ns("normalize"),
                 "Normalize by Total Counts",
                 value = TRUE
               ),
               sliderInput(
                 ns("top_n_exp"),
                 "Show Top N taxa:",
                 min = 5, max = 50, value = 10
               )
             )
      ),
      column(9,
             card(
               title = "Aggregated Abundance Data",
               DT::dataTableOutput(ns("agg_abundance_table"))
             )
      )
    ),
    card(
      title = "Aggregated Visualization",
      plotlyOutput(ns("agg_plot"))
    )
  )
}

experiment_aggregator_server <- function(input, output, session, physeq_obj) {
  ns <- session$ns
  
  # Reactive: Get available metadata variables
  available_vars <- reactive({
    meta <- sample_data(physeq_obj())
    vars <- colnames(meta)
    # Filter out any variables that might not be suitable for grouping
    vars[!vars %in% c("Sample", "sample_id")]  # Adjust based on your metadata
  })
  
  # Update experiment variable selector
  observe({
    updateSelectInput(
      session,
      "experiment_var",
      choices = available_vars(),
      selected = if(length(available_vars()) > 0) available_vars()[1] else NULL
    )
  })
  
  # Reactive: Get aggregated abundance data
  agg_abundance_data <- reactive({
    req(input$tax_level_exp, input$experiment_var, physeq_obj())
    
    # Extract data
    otu_mat <- otu_table(physeq_obj())
    tax_mat <- tax_table(physeq_obj())
    meta <- sample_data(physeq_obj())
    
    # Convert to data frames
    otu_df <- as.data.frame(otu_mat)
    otu_df$Taxon <- rownames(otu_mat)
    
    # Add taxonomic level
    tax_df <- as.data.frame(tax_mat)
    tax_df$Taxon <- rownames(tax_mat)
    
    # Merge
    merged_df <- merge(otu_df, tax_df[, c("Taxon", input$tax_level_exp)], by = "Taxon")
    
    # Add sample metadata
    meta_long <- meta
    meta_long$Sample <- rownames(meta)
    
    # Reshape to long format
    library(tidyr)
    long_df <- merged_df %>%
      pivot_longer(
        cols = -c(Taxon, input$tax_level_exp),
        names_to = "Sample",
        values_to = "Abundance"
      )
    
    # Add metadata
    final_df <- merge(long_df, meta_long, by = "Sample")
    
    # Group and aggregate
    agg_df <- final_df %>%
      group_by(.data[[input$experiment_var]], .data[[input$tax_level_exp]]) %>%
      summarise(Abundance = sum(Abundance), .groups = 'drop') %>%
      filter(.data[[input$tax_level_exp]] != "", !is.na(.data[[input$tax_level_exp]]))
    
    # Normalize if requested
    if(input$normalize) {
      agg_df <- agg_df %>%
        group_by(.data[[input$experiment_var]]) %>%
        mutate(Abundance = Abundance / sum(Abundance) * 100) %>%
        ungroup()
    }
    
    # Get top N taxa overall
    top_taxa <- agg_df %>%
      group_by(.data[[input$tax_level_exp]]) %>%
      summarise(total_abundance = sum(Abundance), .groups = 'drop') %>%
      arrange(desc(total_abundance)) %>%
      slice_head(n = input$top_n_exp) %>%
      pull(.data[[input$tax_level_exp]])
    
    # Filter data to top taxa
    result <- agg_df %>%
      filter(.data[[input$tax_level_exp]] %in% top_taxa) %>%
      arrange(desc(Abundance))
    
    return(result)
  })
  
  # Output: Aggregated abundance table
  output$agg_abundance_table <- DT::renderDataTable({
    req(agg_abundance_data())
    
    agg_abundance_data() %>%
      DT::datatable(
        options = list(
          pageLength = 15,
          lengthMenu = list(c(10, 15, 25, -1), c(10, 15, 25, "All"))
        )
      )
  })
  
  # Output: Aggregated plot
  output$agg_plot <- renderPlotly({
    req(agg_abundance_data())
    
    p <- ggplot(agg_abundance_data(), 
                aes(x = .data[[input$experiment_var]], 
                    y = Abundance, 
                    fill = .data[[input$tax_level_exp]])) +
      geom_col(position = "dodge") +
      labs(
        title = paste("Aggregated Abundance by", input$experiment_var),
        x = input$experiment_var,
        y = if(input$normalize) "Relative Abundance (%)" else "Abundance",
        fill = input$tax_level_exp
      ) +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
    
    ggplotly(p, tooltip = c("x", "y", "fill"))
  })
}