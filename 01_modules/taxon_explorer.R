# modules/taxon_explorer.R


taxon_explorer_ui <- function(id) {
  ns <- NS(id)
  
  layout_sidebar(
    sidebar = sidebar(
      selectInput(ns("taxon_level"), "Taxon Level:", choices = c("Genus", "Species")),
      uiOutput(ns("taxon_selector")),
      checkboxInput(
        ns("include_experiments"),
        "Include samples from experiments",
        value = TRUE
      ),
      width = 280
    ),
    DTOutput(ns("sample_table")),
    fill = TRUE,
    fillable = TRUE
  )
}

taxon_explorer_server <- function(id, ps) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns
    
    # Precompute unique taxon choices with proper labels
    tax_df <- as.data.frame(tax_table(ps))
    
    # Ensure required columns exist
    stopifnot("Genus" %in% colnames(tax_df), "Species" %in% colnames(tax_df))
    
    # For Genus: just unique genus names (remove NAs)
    unique_genus <- sort(unique(tax_df$Genus[!is.na(tax_df$Genus)]))
    
    # For Species: create "Genus species" labels, but only where both are non-NA
    species_rows <- !is.na(tax_df$Genus) & !is.na(tax_df$Species)
    species_labels <- paste(tax_df$Genus[species_rows], tax_df$Species[species_rows])
    # Keep them in the same order as rownames for later matching
    names(species_labels) <- rownames(tax_df)[species_rows]  # map back to OTU IDs
    
    # Remove duplicates in label (though shouldn't happen if data is clean)
    # But we keep all OTUs that match a label — handled in matching logic
    unique_species_labels <- sort(unique(species_labels))
    
    # UI for taxon selector (with server-side selectize)
    output$taxon_selector <- renderUI({
      selectizeInput(
        ns("taxon_selector"),
        "Select Taxon:",
        choices = character(0),  # populated dynamically
        options = list(
          placeholder = "Start typing...",
          searchConjunction = "and",
          maxOptions = 10000
        ),
        multiple = FALSE,
        selected = NULL
      )
    })
    
    # Update choices based on taxon level using server-side rendering
    observe({
      level <- input$taxon_level
      if (is.null(level)) return()
      
      if (level == "Genus") {
        choices <- setNames(unique_genus, unique_genus)  # value = label
      } else if (level == "Species") {
        choices <- setNames(unique_species_labels, unique_species_labels)
      } else {
        choices <- character(0)
      }
      
      # Use updateSelectizeInput with server = TRUE for performance
      updateSelectizeInput(
        session,
        "taxon_selector",
        choices = choices,
        selected = NULL,
        server = TRUE  # 👈 THIS enables server-side filtering
      )
    })
    
    matching_samples <- reactive({
      req(input$taxon_selector)
      taxon_label <- input$taxon_selector
      level <- input$taxon_level
      include_experiments <- input$include_experiments
      
      if (is.null(taxon_label) || taxon_label == "") return(data.frame())
      
      rel_abund_mat <- as.matrix(otu_table(ps))
      
      # Find matching OTUs based on level and label
      if (level == "Genus") {
        matching_otus <- rownames(tax_df)[
          !is.na(tax_df$Genus) & 
            tax_df$Genus == taxon_label
        ]
      } else if (level == "Species") {
        # Match on combined "Genus Species" label
        # Reconstruct labels for all rows (with NAs handled)
        full_labels <- ifelse(
          !is.na(tax_df$Genus) & !is.na(tax_df$Species),
          paste(tax_df$Genus, tax_df$Species),
          NA_character_
        )
        matching_otus <- names(species_labels)[species_labels == taxon_label]
        # Alternative (slower but clearer):
        # matching_otus <- rownames(tax_df)[full_labels == taxon_label]
      } else {
        return(data.frame())
      }
      
      if (length(matching_otus) == 0) return(data.frame())
      
      # Compute abundance per sample
      rel_abund_per_sample <- colSums(rel_abund_mat[matching_otus, , drop = FALSE])
      sample_ids <- names(rel_abund_per_sample)[rel_abund_per_sample > 0]
      rel_abund_vals <- rel_abund_per_sample[rel_abund_per_sample > 0]
      
      if (length(sample_ids) == 0) return(data.frame())
      
      # Extract metadata
      meta_df <- as.data.frame(sample_data(ps))
      meta_df$sample_id <- rownames(meta_df)
      
      valid_samples <- sample_ids[sample_ids %in% meta_df$sample_id]
      if (length(valid_samples) == 0) return(data.frame())
      
      meta_subset <- meta_df[meta_df$sample_id %in% valid_samples, , drop = FALSE]
      meta_subset <- meta_subset[match(valid_samples, meta_subset$sample_id), , drop = FALSE]
      
      # Apply experiment filter
      exp_col <- meta_subset$experiment  # adjust if column name differs
      
      if (include_experiments) {
        keep_rows <- !is.na(exp_col)
      } else {
        keep_rows <- is.na(exp_col)
      }
      
      filtered_samples <- valid_samples[keep_rows]
      if (length(filtered_samples) == 0) return(data.frame())
      
      filtered_abund <- rel_abund_vals[filtered_samples]
      filtered_meta <- meta_subset[keep_rows, , drop = FALSE]
      
      # Validate metadata columns
      required_cols <- c("seq_type", "experiment_name", "notes")
      missing_cols <- setdiff(required_cols, colnames(filtered_meta))
      if (length(missing_cols) > 0) {
        warning("Missing metadata columns: ", paste(missing_cols, collapse = ", "))
        return(data.frame())
      }
      
      result <- data.frame(
        sample_id = filtered_samples,
        relative_abundance = filtered_abund,
        seq_type = filtered_meta$seq_type,
        experiment_name = filtered_meta$experiment_name,
        notes = filtered_meta$notes,
        row.names = NULL
      )
      
      result <- result[order(-result$relative_abundance), , drop = FALSE]
      result$relative_abundance <- round(result$relative_abundance, 6)
      
      if (nrow(result) == 0) return(data.frame())
      result
    })
    
    output$sample_table <- renderDT({
      df <- matching_samples()
      req(nrow(df) > 0)
      
      datatable(
        df,
        rownames = FALSE,
        width = "100%",
        options = list(
          pageLength = 25,
          responsive = TRUE,
          autoWidth = FALSE,
          scrollX = FALSE,
          columnDefs = list(
            list(className = "dt-center", targets = 1)
          )
        ),
        class = "compact"
      ) %>%
        DT::formatRound("relative_abundance", digits = 4)
    })
  })
}