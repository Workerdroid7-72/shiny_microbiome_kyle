# sample_explorer.R

sample_explorer_ui <- function(id) {
  ns <- NS(id)
  
  layout_sidebar(
    sidebar = sidebar(
      width = 280,
      textInput(ns("sample_id"), "Enter 4-digit Sample ID:", value = ""),
      selectInput(
        ns("tax_level"), 
        "Taxonomic Level:", 
        choices = c("Family", "Genus", "Species")
      ),
      selectInput(
        ns("top_n"), 
        "Show top taxa:", 
        choices = c("Top 5" = 5, "Top 10" = 10, "Top 20" = 20),  # ✅ Updated choices
        selected = 10  # ✅ Default to "Top 10"
      ),
      actionButton(ns("go"), "Show Results", class = "btn-primary")
      
    ),
    plotOutput(ns("taxa_plot")),
    DTOutput(ns("taxa_table"))
  )
  
}


sample_explorer_server <- function(id, ps) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns
    
    sample_taxa <- eventReactive(input$go, {
      short_id <- input$sample_id
      
      # Validate: must be 4 digits
      if (nchar(short_id) != 4 || !grepl("^[0-9]{4}$", short_id)) {
        showNotification("Please enter exactly 4 digits (e.g., 1234).", type = "warning")
        return(NULL)
      }
      
      # 🔑 Reconstruct full sample ID
      full_id <- paste0("ESSE.222.", short_id)
      
      # Check if full_id exists in phyloseq object
      if (!full_id %in% sample_names(ps)) {
        showNotification(paste("Sample", full_id, "not found!"), type = "error")
        return(NULL)
      }
      
      # Get counts for this sample
      sample_otu <- otu_table(ps)[, full_id, drop = FALSE]
      tax_df <- as.data.frame(tax_table(ps))
      
      # Merge OTU + taxonomy
      df <- data.frame(
        taxon = rownames(sample_otu),
        count = as.numeric(sample_otu[, 1]),
        stringsAsFactors = FALSE
      )
      df <- cbind(df, tax_df[match(df$taxon, rownames(tax_df)), , drop = FALSE])
      
      # Handle missing taxonomy gracefully
      df[is.na(df)] <- ""
      
      # Create display label based on taxonomic level
      if (input$tax_level == "Species") {
        # Combine Genus and Species, but handle cases where Genus is missing
        df$display_label <- trimws(paste(df$Genus, df$Species))
        # Optional: clean up if both are empty or just whitespace
        df$display_label[df$display_label == ""] <- "Unknown"
        group_var <- "display_label"
      } else {
        group_var <- input$tax_level
        df$display_label <- df[[group_var]]
      }
      
      # Aggregate by display_label
      df_agg <- df %>%
        group_by(display_label) %>%
        summarise(count = sum(count), .groups = "drop") %>%
        filter(display_label != "" & !is.na(display_label)) %>%
        arrange(desc(count))
      
      # Top N
      df_agg <- head(df_agg, as.numeric(input$top_n))
      
      # Store the taxonomic level used for this result
      df_agg$used_tax_level <- input$tax_level
      
      # Rename column for consistency in plot/table labels
      colnames(df_agg)[colnames(df_agg) == "display_label"] <- input$tax_level
      
      df_agg
    })
    
    output$taxa_plot <- renderPlot({
      req(sample_taxa())
      
      # Use the taxonomic level that was used when the data was generated
      current_tax_level <- sample_taxa()$used_tax_level[1]  # All rows have the same level
      
      # If current_tax_level is NULL or different from input$tax_level, don't render
      # Only render if the tax level matches what was used to generate the data
      if(is.null(current_tax_level) || current_tax_level != input$tax_level) {
        # Don't render anything until the button is pressed again
        return(NULL)
      }
      
      ggplot(sample_taxa(), aes(x = reorder(!!sym(input$tax_level), count), y = count)) +
        geom_col(fill = "#45B29D") +
        coord_flip() +
        labs(
          x = ifelse(input$tax_level == "Species", "Genus Species", input$tax_level),
          y = "Relative Abundance",
          title = paste("Top Taxa in Sample")
        ) +
        theme_minimal()
    })
    
    output$taxa_table <- renderDT({
      req(sample_taxa())
      
      # Only render if the tax level matches what was used to generate the data
      current_tax_level <- sample_taxa()$used_tax_level[1]  # All rows have the same level
      
      if(is.null(current_tax_level) || current_tax_level != input$tax_level) {
        # Return an empty table until the button is pressed again
        return(datatable(data.frame(), options = list(pageLength = 10)))
      }
      
      datatable(sample_taxa(), options = list(pageLength = 10)) %>%
        formatRound(columns = "count", digits = 4)
    })
  })
}