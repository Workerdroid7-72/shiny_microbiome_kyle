# modules/experiment_analysis.R

# UI for Experiment Analysis Tab
experiment_analysis_ui <- function(id) {
  ns <- NS(id)
  
  tagList(
    card(
      card_header("Aggregate Analysis by Experiment"),
      layout_sidebar(
        sidebar = sidebar(
          width = 300,
          checkboxGroupInput(ns("exp_select"), "Select Experiments:",
                             choices = experiments,
                             selected = experiments[1]),
          selectInput(ns("tax_level_aggregate"), "Taxonomic Level:",
                      choices = taxonomic_levels, selected = "Genus"),
          numericInput(ns("top_n_agg"), "Show Top N Taxa:",
                       value = 10, min = 5, max = 20),
          actionButton(ns("update_agg"), "Update Analysis")
        ),
        navset_card_underline(
          nav_panel("Table", DTOutput(ns("aggregate_table"))),
          nav_panel("Plot", plotlyOutput(ns("aggregate_plot")))
        )
      )
    )
  )
}

# Server for Experiment Analysis Tab
experiment_analysis_server <- function(id, phyloseq_obj) {
  moduleServer(id, function(input, output, session) {
    
    # Reactive expression for aggregate data
    aggregate_data <- reactive({
      req(input$exp_select)
      
      # Filter phyloseq by selected experiments
      ps_filtered <- filter_phyloseq_by_experiment(phyloseq_obj, input$exp_select)
      
      # Get abundance table
      abun_table <- get_abundance_table(ps_filtered, input$tax_level_aggregate)
      
      return(abun_table)
    }) %>% bindEvent(input$update_agg)
    
    # Aggregate Table
    output$aggregate_table <- renderDT({
      data <- aggregate_data()
      
      # Calculate mean abundances across samples
      summary_data <- data %>%
        select(-SampleID) %>%
        summarise(across(everything(), mean, na.rm = TRUE)) %>%
        pivot_longer(cols = everything(), names_to = "Taxa", values_to = "Mean_Abundance") %>%
        arrange(desc(Mean_Abundance))
      
      datatable(
        summary_data,
        options = list(
          pageLength = 10,
          scrollX = TRUE
        ),
        rownames = FALSE
      )
    })
    
    # Aggregate Plot
    output$aggregate_plot <- renderPlotly({
      data <- aggregate_data()
      
      # Prepare data for plotting
      plot_data <- data %>%
        select(-SampleID) %>%
        summarise(across(everything(), mean, na.rm = TRUE)) %>%
        pivot_longer(cols = everything(), names_to = "Taxa", values_to = "Mean_Abundance") %>%
        arrange(desc(Mean_Abundance)) %>%
        slice_head(n = input$top_n_agg)
      
      p <- ggplot(plot_data, aes(x = reorder(Taxa, Mean_Abundance), y = Mean_Abundance, fill = Taxa)) +
        geom_col() +
        coord_flip() +
        labs(x = input$tax_level_aggregate, y = "Mean Abundance",
             title = paste("Experiments:", paste(input$exp_select, collapse = ", "))) +
        theme_minimal() +
        theme(legend.position = "none")
      
      ggplotly(p)
    })
  })
}