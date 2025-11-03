# modules/experiment_explorer.R

# experiment_explorer.R

experiment_explorer_ui <- function(id) {
  ns <- NS(id)
  
  layout_sidebar(
    sidebar = sidebar(
      width = 280,
      selectInput(
        ns("tax_level"),
        "Taxonomic Level",
        choices = c("Family", "Genus", "Species"),
        selected = "Genus"
      ),
      selectInput(
        ns("top_n"),
        "Number of Top Taxa",
        choices = c(
          "Top 5" = 5,
          "Top 10" = 10,
          "Top 20" = 20
        ),
        selected = 10
      ),
      selectInput(ns("experiment"), "Select Experiment", choices = NULL),
      uiOutput(ns("sample_ui")),
      actionButton(ns("apply_btn"), "Show Results", class = "btn-primary")
    ),
    # plotOutput(NS(id, "taxa_plot"), height = "400px"),
    # DTOutput(NS(id, "taxa_table"))
    plotOutput(ns("taxa_plot"), height = "400px"),
    # Namespaced
    DTOutput(ns("taxa_table"))                      # Namespaced
  )
}

experiment_explorer_server <- function(id, ps) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns
    
    # Ensure 'experiment' column exists
    if (!"experiment" %in% names(sample_data(ps))) {
      stop("Your phyloseq object's sample_data must contain a column named 'experiment'.")
    }
    
    
    
    # Populate experiment dropdown once (ensure character)
    experiments <- sort(unique(as.character(sample_data(ps)$experiment_name)))
    updateSelectInput(session, "experiment", choices = experiments)
    
    # Render UI for sample selector (static structure)
    output$sample_ui <- renderUI({
      selectInput(ns("sample"),
                  "Select Sample",
                  choices = character(0),
                  multiple = FALSE)
    })
    
    # Update sample choices when experiment_name changes
    observeEvent(input$experiment, {
      req(input$experiment)
      
      # Safely get sample data
      sample_df <- as.data.frame(sample_data(ps), stringsAsFactors = FALSE)
      
      # Ensure both are character and clean whitespace
      sample_df$experiment_name <- trimws(as.character(sample_df$experiment_name))
      selected_exp <- trimws(as.character(input$experiment))
      
      # Debug (remove after testing)
      # cat("Selected:", selected_exp, "| Available:", paste(unique(sample_df$experiment_name), collapse = ", "), "\n")
      
      # Find matching samples
      samples_in_exp <- rownames(sample_df)[sample_df$experiment_name == selected_exp]
      
      # Update sample selector
      if (length(samples_in_exp) == 0) {
        updateSelectInput(session,
                          "sample",
                          choices = character(0),
                          selected = character(0))
      } else {
        updateSelectInput(
          session,
          "sample",
          choices = sort(samples_in_exp),
          selected = samples_in_exp[1]
        )
      }
    })
    
    
    
    # Add this to see if button is firing
    observe({
      cat("DEBUG: input$apply_btn value:", input$apply_btn, "\n")
    })
    
    # Reactive values to store results and status
    results <- reactiveValues(plot_data = NULL,
                              table_data = NULL,
                              ready = FALSE)
    
    # Compute and display results when button is clicked
    observeEvent(input$apply_btn, {
      req(input$tax_level,
          input$top_n,
          input$experiment,
          input$sample)
      
      tax_level <- input$tax_level
      top_n <- as.numeric(input$top_n)
      sample_id <- input$sample
      
      cat("DEBUG: Processing sample:",
          sample_id,
          "for tax level:",
          tax_level,
          "\n")
      
      # Wrap everything in tryCatch to catch errors
      tryCatch({
        # Validate tax level exists in phyloseq
        if (!tax_level %in% rank_names(ps)) {
          showModal(
            modalDialog(
              title = "Error",
              "Selected taxonomic level not found in phyloseq taxonomy table.",
              easyClose = TRUE
            )
          )
          return()
        }
        
        # Get counts for this sample from original ps
        sample_otu <- otu_table(ps)[, sample_id, drop = FALSE]
        tax_df <- as.data.frame(tax_table(ps))
        
        # Merge OTU + taxonomy
        df <- data.frame(
          taxon = rownames(sample_otu),
          relative_abundance = as.numeric(sample_otu[, 1]),
          # Use correct column name
          stringsAsFactors = FALSE
        )
        df <- cbind(df, tax_df[match(df$taxon, rownames(tax_df)), , drop = FALSE])
        
        # Handle missing taxonomy gracefully
        df[is.na(df)] <- ""
        
        # Create display label based on taxonomic level
        if (tax_level == "Species") {
          # Combine Genus and Species, but handle cases where Genus is missing
          df$display_label <- trimws(paste(df$Genus, df$Species))
          # Clean up if both are empty or just whitespace
          df$display_label[df$display_label == ""] <- "Unknown"
          group_var <- "display_label"
        } else {
          group_var <- tax_level
          df$display_label <- df[[tax_level]]
        }
        
        # Aggregate by display_label
        df_agg <- df %>%
          filter(relative_abundance > 0) %>%  # Use correct column name
          group_by(!!sym("display_label")) %>%
          summarise(relative_abundance = sum(relative_abundance),
                    .groups = "drop") %>%  # Use correct column name
          filter(display_label != "" &
                   display_label != "Unknown") %>%
          arrange(desc(relative_abundance))  # Use correct column name
        
        # Top N
        df_agg <- head(df_agg, top_n)
        
        # Rename column for consistency
        colnames(df_agg)[colnames(df_agg) == "display_label"] <- "TaxonLabel"
        
        cat("DEBUG: After processing - nrow:", nrow(df_agg), "\n")
        if (nrow(df_agg) > 0) {
          cat("DEBUG: First few TaxonLabels:",
              paste(head(df_agg$TaxonLabel, 3), collapse = ", "),
              "\n")
        }
        
        # Store results in reactive values
        results$plot_data <- df_agg
        results$table_data <- df_agg
        results$ready <- TRUE
        results$sample_id <- sample_id
        results$tax_level <- tax_level
        
        cat("DEBUG: Stored results. plot_data nrow:",
            nrow(results$plot_data),
            "\n")
        
      }, error = function(e) {
        cat("ERROR in observeEvent: ", conditionMessage(e), "\n")
        showModal(modalDialog(
          title = "Error",
          paste("An error occurred:", conditionMessage(e)),
          easyClose = TRUE
        ))
      })
    })
    
    # Render plot using the namespaced output ID
    output$taxa_plot <- renderPlot({
      cat("DEBUG: renderPlot called. req(results$ready)...\n")
      req(results$ready)  # Only depends on results$ready
      cat("DEBUG: renderPlot - data available. nrow:",
          nrow(results$plot_data),
          "\n")
      
      df_clean <- results$plot_data
      sample_id <- results$sample_id
      tax_level <- results$tax_level
      
      p <- ggplot(df_clean, aes(
        x = reorder(TaxonLabel, relative_abundance),
        y = relative_abundance
      )) +  # Use correct column name
        geom_col(fill = "#2c7bb6") +
        coord_flip() +
        labs(
          x = "Taxon",
          y = "Relative Abundance",
          # Update y-axis label
          title = paste("Top", nrow(df_clean), tax_level, "in Sample:", sample_id)
        ) +
        theme_minimal(base_size = 10)
      
      print(p)
    })
    
    # Render table using the namespaced output ID
    output$taxa_table <- renderDT({
      cat("DEBUG: renderDT called. req(results$ready)...\n")
      req(results$ready)  # <-- NEW: Require ready flag
      cat("DEBUG: renderDT - data available. nrow:",
          nrow(results$table_data),
          "\n")
      
      df_clean <- results$table_data
      
      datatable(
        df_clean,
        options = list(pageLength = min(10, nrow(df_clean)), scrollX = TRUE),
        rownames = FALSE
      )
    })
    
  })
}