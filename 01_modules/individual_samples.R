# modules/individual_samples.R

# UI for Individual Samples Tab
individual_samples_ui <- function(id) {
  ns <- NS(id)
  
  tagList(
    card(
      card_header("Sample Selection and Visualization"),
      layout_sidebar(
        sidebar = sidebar(
          width = 300,
          selectInput(ns("sample_select"), "Select Sample:", 
                      choices = NULL),
          selectInput(ns("tax_level_individual"), "Taxonomic Level:",
                      choices = taxonomic_levels, selected = "Genus"),
          numericInput(ns("top_n"), "Show Top N Taxa:", 
                       value = 10, min = 5, max = 20)
        ),
        navset_card_underline(
          nav_panel("Table", DTOutput(ns("individual_table"))),
          nav_panel("Plot", plotlyOutput(ns("individual_plot")))
        )
      )
    )
  )
}

# Server for Individual Samples Tab
individual_samples_server <- function(id, phyloseq_obj) {
  moduleServer(id, function(input, output, session) {
    
    # Update sample choices based on phyloseq object
    observe({
      sample_choices <- sample_names(phyloseq_obj)
      updateSelectInput(session, "sample_select", choices = sample_choices)
    })
    
    # Reactive expression for individual sample data
    individual_data <- reactive({
      req(input$sample_select, input$tax_level_individual)
      
      # Get abundance table at specified taxonomic level
      abun_table <- get_abundance_table(phyloseq_obj, input$tax_level_individual)
      
      # Filter for selected sample
      sample_data <- abun_table %>%
        filter(SampleID == input$sample_select) %>%
        select(-SampleID) %>%
        t() %>%
        as.data.frame() %>%
        rownames_to_column("Taxa") %>%
        setNames(c("Taxa", "Abundance")) %>%
        filter(Abundance > 0) %>%
        arrange(desc(Abundance))
      
      return(sample_data)
    })
    
    # Individual Sample Table
    output$individual_table <- renderDT({
      data <- individual_data()
      datatable(
        data,
        options = list(
          pageLength = 10,
          scrollX = TRUE,
          dom = 'tip'
        ),
        rownames = FALSE
      )
    })
    
    # Individual Sample Plot
    output$individual_plot <- renderPlotly({
      data <- individual_data() %>%
        arrange(desc(Abundance)) %>%
        slice_head(n = input$top_n)
      
      p <- ggplot(data, aes(x = reorder(Taxa, Abundance), y = Abundance, fill = Taxa)) +
        geom_col() +
        coord_flip() +
        labs(x = input$tax_level_individual, y = "Abundance", 
             title = paste("Sample:", input$sample_select)) +
        theme_minimal() +
        theme(legend.position = "none")
      
      ggplotly(p)
    })
  })
}