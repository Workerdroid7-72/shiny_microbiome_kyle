# sample_explorer.R

sample_explorer_ui <- function(id) {
  ns <- NS(id)
  
  layout_sidebar(
    sidebar = sidebar(
      width = 280,
      textInput(ns("sample_id"), "Enter Sample ID (3 or 4 digits):", value = ""),
      selectInput(
        ns("tax_level"), 
        "Taxonomic Level:", 
        choices = c("Family", "Genus", "Species")
      ),
      selectInput(
        ns("top_n"), 
        "Show top taxa:", 
        choices = c("Top 5" = 5, "Top 10" = 10, "Top 20" = 20),
        selected = 10
      ),
      actionButton(ns("go"), "Show Results", class = "btn-primary"),
      hr(),  # optional separator
      htmlOutput(ns("sample_metadata")),
      
      
      # Hidden input to receive clicked taxon — properly hidden
      div(
        textInput(ns("clicked_taxon"), label = NULL, value = ""),
        style = "display: none;"
      )
    ),
    plotOutput(ns("taxa_plot")),
    DTOutput(ns("taxa_table"))
  )
}

# sample_explorer_server <- function(id, ps) {
#   moduleServer(id, function(input, output, session) {
#     ns <- session$ns
# 
#     current_tax_level_reactive <- reactiveVal(NULL)
# 
#     sample_taxa <- eventReactive(input$go, {
#       short_id <- input$sample_id
# 
#       # NEW validation: allow 3 or 4 digits
#       if (!grepl("^[0-9]{3,4}$", short_id)) {
#         showNotification("Please enter 3 or 4 digits (e.g., 123 or 1234).", type = "warning")
#         return(NULL)
#       }
# 
#       full_id <- paste0("ESSE.222.", short_id)
# 
#       if (!full_id %in% sample_names(ps)) {
#         showNotification(paste("Sample", full_id, "not found!"), type = "error")
#         return(NULL)
#       }
# 
#       sample_otu <- otu_table(ps)[, full_id, drop = FALSE]
#       tax_df <- as.data.frame(tax_table(ps))
# 
#       df <- data.frame(
#         taxon = rownames(sample_otu),
#         relative_abundance = as.numeric(sample_otu[, 1]),
#         stringsAsFactors = FALSE
#       )
#       df <- cbind(df, tax_df[match(df$taxon, rownames(tax_df)), , drop = FALSE])
#       df[is.na(df)] <- ""
# 
#       if (input$tax_level == "Species") {
#         df$display_label <- trimws(paste(df$Genus, df$Species))
#         df$display_label[df$display_label == ""] <- "Unknown"
#         group_var <- "display_label"
#       } else {
#         group_var <- input$tax_level
#         df$display_label <- df[[group_var]]
#       }
# 
#       df_agg <- df %>%
#         group_by(display_label) %>%
#         summarise(relative_abundance = sum(relative_abundance), .groups = "drop") %>%
#         filter(display_label != "" & !is.na(display_label)) %>%
#         arrange(desc(relative_abundance))
# 
#       df_agg <- head(df_agg, as.numeric(input$top_n))
# 
#       df_agg$used_tax_level <- input$tax_level
#       colnames(df_agg)[colnames(df_agg) == "display_label"] <- input$tax_level
# 
#       current_tax_level_reactive(input$tax_level)
# 
#       df_agg
#     })
# 
#     sample_metadata_text <- eventReactive(input$go, {
#       short_id <- input$sample_id
#       if (!grepl("^[0-9]{3,4}$", short_id)) return(NULL)
# 
#       full_id <- paste0("ESSE.222.", short_id)
#       if (!full_id %in% sample_names(ps)) return(NULL)
# 
#       # Extract sample metadata correctly:
#       meta_df <- as.data.frame(sample_data(ps))
#       # Ensure rownames match sample names (they usually do, but let's be safe)
#       rownames(meta_df) <- sample_names(ps)
# 
#       # Now safely extract the row
#       if (!full_id %in% rownames(meta_df)) return(NULL)
# 
#       row_meta <- meta_df[full_id, , drop = FALSE]  # keep as 1-row data.frame
# 
#       # Safely extract values (as character, handling NA/empty)
#       get_val <- function(col) {
#         if (col %in% names(row_meta)) {
#           val <- as.character(row_meta[[col]])
#           if (is.na(val) || val == "" || val == "NA") "(none)" else val
#         } else {
#           "(column missing)"
#         }
#       }
# 
#       exp_name <- get_val("experiment_name")
#       notes_val <- get_val("notes")
#       seq_type_val <- get_val("seq_type")
# 
#       # Format as HTML
#       glue::glue(
#         "<strong>Experiment:</strong> {exp_name}<br/>",
#         "<strong>Notes:</strong> {notes_val}<br/>",
#         "<strong>Seq Type:</strong> {seq_type_val}"
#       )
#     })
# 
#     output$taxa_plot <- renderPlot({
#       req(sample_taxa())
#       current_tax_level <- sample_taxa()$used_tax_level[1]
#       if (is.null(current_tax_level) || current_tax_level != input$tax_level) {
#         return(NULL)
#       }
# 
#       # Precompute the label (optional but cleaner)
#       # You may want to round or format the values for readability
#       plot_data <- sample_taxa() %>%
#         mutate(
#           label = round(relative_abundance, digits = 3),  # adjust digits as needed
#           taxon = reorder(!!sym(input$tax_level), relative_abundance)
#         )
# 
#       ggplot(plot_data, aes(x = taxon, y = relative_abundance)) +
#         geom_col(fill = "#45B29D") +
#         geom_text(
#           aes(label = ifelse(relative_abundance > 0.03, label, "")),  # only show if >2%
#           hjust = 1.1,          # places text just inside the right edge of each bar
#           color = "white",      # contrast well with your bar color
#           fontface = "bold",
#           size = 3.5
#         ) +
#         coord_flip() +
#         labs(
#           x = ifelse(input$tax_level == "Species", "Genus Species", input$tax_level),
#           y = "Relative Abundance",
#           title = "Top Taxa in Sample"
#         ) +
#         theme_minimal() +
#         theme(
#           axis.text.x = element_text(size = 9),
#           axis.text.y = element_text(size = 9)
#         )
#     })
# 
#     output$taxa_table <- renderDT({
#       req(sample_taxa())
#       current_tax_level <- sample_taxa()$used_tax_level[1]
#       if (is.null(current_tax_level) || current_tax_level != input$tax_level) {
#         return(datatable(data.frame(), options = list(pageLength = 10)))
#       }
# 
#       df_clean <- sample_taxa() %>%
#         select(-used_tax_level)
# 
#       # Add "Compare" button column as HTML
#       df_clean$compare <- sprintf(
#         '<button type="button" class="btn btn-xs btn-outline-secondary" data-taxon="%s">%s Compare</button>',
#         gsub('"', '\\"', df_clean[[input$tax_level]]),
#         as.character(shiny::icon("search"))
#       )
# 
#       datatable(
#         df_clean,
#         escape = FALSE,
#         options = list(
#           dom = 'Bfrtip',
#           pageLength = 10,
#           columnDefs = list(
#             list(orderable = FALSE, targets = ncol(df_clean))  # last column = button
#           )
#         ),
#         callback = JS(paste0("
#           table.on('click', 'button', function() {
#             var taxon = $(this).data('taxon');
#             Shiny.setInputValue('", ns("clicked_taxon"), "', taxon);
#           });
#         "))
#       ) %>%
#         formatRound(columns = "relative_abundance", digits = 6)
#     })
# 
#     # Observe clicks via hidden input
#     
#     observeEvent(input$clicked_taxon, {
#       req(input$clicked_taxon != "")
# 
#       local_taxon <- input$clicked_taxon
#       tax_level <- isolate(current_tax_level_reactive())
# 
#       if (is.null(tax_level)) {
#         showNotification("Please run the analysis first.", type = "warning")
#         return()
#       }
# 
#       # Map OTUs to taxonomy
#       tax_df_full <- as.data.frame(tax_table(ps), stringsAsFactors = FALSE)
#       rownames(tax_df_full) <- taxa_names(ps)
# 
#       # Find matching OTUs
#       if (tax_level == "Species") {
#         combined <- trimws(paste(tax_df_full$Genus, tax_df_full$Species))
#         matches <- combined == local_taxon
#         if (local_taxon == "Unknown") {
#           matches <- (tax_df_full$Genus == "" | is.na(tax_df_full$Genus)) &
#             (tax_df_full$Species == "" | is.na(tax_df_full$Species))
#         }
#       } else {
#         matches <- tax_df_full[[tax_level]] == local_taxon
#       }
# 
#       matching_otus <- rownames(tax_df_full)[matches]
#       matching_otus <- matching_otus[matching_otus %in% taxa_names(ps)]
# 
#       if (length(matching_otus) == 0) {
#         showNotification("No OTUs match this taxon.", type = "warning")
#         return()
#       }
# 
#       # Get relative abundances across all samples
#       otu_mat <- otu_table(ps)
#       rel_abund_per_sample <- colSums(otu_mat[matching_otus, , drop = FALSE])
#       result_df <- data.frame(
#         sample_id = names(rel_abund_per_sample),
#         relative_abundance = as.numeric(rel_abund_per_sample),
#         stringsAsFactors = FALSE
#       ) %>%
#         filter(relative_abundance > 0) %>%
#         arrange(desc(relative_abundance))
# 
#       if (nrow(result_df) == 0) {
#         showNotification("This taxon is absent from all samples.", type = "info")
#         return()
#       }
# 
#       showModal(modalDialog(
#         title = paste("Samples containing:", local_taxon),
#         size = "l",
#         DTOutput(ns("compare_table")),
#         footer = modalButton("Close")
#       ))
# 
#       output$compare_table <- renderDT({
#         datatable(result_df, options = list(dom = 'Bfrtip',pageLength = 10)) %>%
#           formatRound(columns = "relative_abundance", digits = 6)
#       })
#     })
#     
#     
#     # Use bindEvent to trigger the modal display only when compare_data_reactive() returns a valid result
#     bindEvent(compare_data_reactive(), {
#       req(compare_data_reactive())  # Only proceed if we have valid data
#       
#       showModal(modalDialog(
#         title = paste("Samples containing:", input$clicked_taxon),
#         size = "l",
#         DTOutput(ns("compare_table")),
#         footer = modalButton("Close")
#       ))
#     })
#     
#     # Render the comparison table using the reactive data
#     output$compare_table <- renderDT({
#       req(compare_data_reactive())
#       datatable(compare_data_reactive(), options = list(dom = 'Bfrtip', pageLength = 10)) %>%
#         formatRound(columns = "relative_abundance", digits = 6)
#     })
# 
#     output$sample_metadata <- renderUI({
#       req(sample_metadata_text())
#       HTML(sample_metadata_text())
#     })
#   })
# }


sample_explorer_server <- function(id, ps) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns
    
    current_tax_level_reactive <- reactiveVal(NULL)
    
    sample_taxa <- eventReactive(input$go, {
      short_id <- input$sample_id
      
      # NEW validation: allow 3 or 4 digits
      if (!grepl("^[0-9]{3,4}$", short_id)) {
        showNotification("Please enter 3 or 4 digits (e.g., 123 or 1234).", type = "warning")
        return(NULL)
      }
      
      full_id <- paste0("ESSE.222.", short_id)
      
      if (!full_id %in% sample_names(ps)) {
        showNotification(paste("Sample", full_id, "not found!"), type = "error")
        return(NULL)
      }
      
      sample_otu <- otu_table(ps)[, full_id, drop = FALSE]
      tax_df <- as.data.frame(tax_table(ps))
      
      df <- data.frame(
        taxon = rownames(sample_otu),
        relative_abundance = as.numeric(sample_otu[, 1]),
        stringsAsFactors = FALSE
      )
      df <- cbind(df, tax_df[match(df$taxon, rownames(tax_df)), , drop = FALSE])
      df[is.na(df)] <- ""
      
      if (input$tax_level == "Species") {
        df$display_label <- trimws(paste(df$Genus, df$Species))
        df$display_label[df$display_label == ""] <- "Unknown"
        group_var <- "display_label"
      } else {
        group_var <- input$tax_level
        df$display_label <- df[[group_var]]
      }
      
      df_agg <- df %>%
        group_by(display_label) %>%
        summarise(relative_abundance = sum(relative_abundance), .groups = "drop") %>%
        filter(display_label != "" & !is.na(display_label)) %>%
        arrange(desc(relative_abundance))
      
      df_agg <- head(df_agg, as.numeric(input$top_n))
      
      df_agg$used_tax_level <- input$tax_level
      colnames(df_agg)[colnames(df_agg) == "display_label"] <- input$tax_level
      
      current_tax_level_reactive(input$tax_level)
      
      df_agg
    })
    
    sample_metadata_text <- eventReactive(input$go, {
      short_id <- input$sample_id
      if (!grepl("^[0-9]{3,4}$", short_id)) return(NULL)
      
      full_id <- paste0("ESSE.222.", short_id)
      if (!full_id %in% sample_names(ps)) return(NULL)
      
      # Extract sample metadata correctly:
      meta_df <- as.data.frame(sample_data(ps))
      # Ensure rownames match sample names (they usually do, but let's be safe)
      rownames(meta_df) <- sample_names(ps)
      
      # Now safely extract the row
      if (!full_id %in% rownames(meta_df)) return(NULL)
      
      row_meta <- meta_df[full_id, , drop = FALSE]  # keep as 1-row data.frame
      
      # Safely extract values (as character, handling NA/empty)
      get_val <- function(col) {
        if (col %in% names(row_meta)) {
          val <- as.character(row_meta[[col]])
          if (is.na(val) || val == "" || val == "NA") "(none)" else val
        } else {
          "(column missing)"
        }
      }
      
      exp_name <- get_val("experiment_name")
      notes_val <- get_val("notes")
      seq_type_val <- get_val("seq_type")
      
      # Format as HTML
      glue::glue(
        "<strong>Experiment:</strong> {exp_name}<br/>",
        "<strong>Notes:</strong> {notes_val}<br/>",
        "<strong>Seq Type:</strong> {seq_type_val}"
      )
    })
    
    output$taxa_plot <- renderPlot({
      req(sample_taxa())
      current_tax_level <- sample_taxa()$used_tax_level[1]
      if (is.null(current_tax_level) || current_tax_level != input$tax_level) {
        return(NULL)
      }
      
      # Precompute the label (optional but cleaner)
      # You may want to round or format the values for readability
      plot_data <- sample_taxa() %>%
        mutate(
          label = round(relative_abundance, digits = 3),  # adjust digits as needed
          taxon = reorder(!!sym(input$tax_level), relative_abundance)
        )
      
      ggplot(plot_data, aes(x = taxon, y = relative_abundance)) +
        geom_col(fill = "#45B29D") +
        geom_text(
          aes(label = ifelse(relative_abundance > 0.03, label, "")),  # only show if >2%
          hjust = 1.1,          # places text just inside the right edge of each bar
          color = "white",      # contrast well with your bar color
          fontface = "bold",
          size = 3.5
        ) +
        coord_flip() +
        labs(
          x = ifelse(input$tax_level == "Species", "Genus Species", input$tax_level),
          y = "Relative Abundance",
          title = "Top Taxa in Sample"
        ) +
        theme_minimal() +
        theme(
          axis.text.x = element_text(size = 9),
          axis.text.y = element_text(size = 9)
        )
    })
    
    output$taxa_table <- renderDT({
      req(sample_taxa())
      current_tax_level <- sample_taxa()$used_tax_level[1]
      if (is.null(current_tax_level) || current_tax_level != input$tax_level) {
        return(datatable(data.frame(), options = list(pageLength = 10)))
      }
      
      df_clean <- sample_taxa() %>%
        select(-used_tax_level)
      
      # Add "Compare" button column as HTML
      df_clean$compare <- sprintf(
        '<button type="button" class="btn btn-xs btn-outline-secondary" data-taxon="%s">%s Compare</button>',
        gsub('"', '\\"', df_clean[[input$tax_level]]),
        as.character(shiny::icon("search"))
      )
      
      datatable(
        df_clean,
        escape = FALSE,
        options = list(
          dom = 'Bfrtip',
          pageLength = 10,
          columnDefs = list(
            list(orderable = FALSE, targets = ncol(df_clean))  # last column = button
          )
        ),
        callback = JS(paste0("
          table.on('click', 'button', function() {
            var taxon = $(this).data('taxon');
            Shiny.setInputValue('", ns("clicked_taxon"), "', taxon);
          });
        "))
      ) %>%
        formatRound(columns = "relative_abundance", digits = 6)
    })
    
    # Refactored section using bindEvent - final corrected version
    compare_data_reactive <- eventReactive(input$clicked_taxon, {
      req(input$clicked_taxon != "")
      
      local_taxon <- input$clicked_taxon
      tax_level <- isolate(current_tax_level_reactive())
      
      if (is.null(tax_level)) {
        showNotification("Please run the analysis first.", type = "warning")
        return(NULL)
      }
      
      # Map OTUs to taxonomy
      tax_df_full <- as.data.frame(tax_table(ps), stringsAsFactors = FALSE)
      rownames(tax_df_full) <- taxa_names(ps)
      
      # Find matching OTUs
      if (tax_level == "Species") {
        combined <- trimws(paste(tax_df_full$Genus, tax_df_full$Species))
        matches <- combined == local_taxon
        if (local_taxon == "Unknown") {
          matches <- (tax_df_full$Genus == "" | is.na(tax_df_full$Genus)) &
            (tax_df_full$Species == "" | is.na(tax_df_full$Species))
        }
      } else {
        matches <- tax_df_full[[tax_level]] == local_taxon
      }
      
      matching_otus <- rownames(tax_df_full)[matches]
      matching_otus <- matching_otus[matching_otus %in% taxa_names(ps)]
      
      if (length(matching_otus) == 0) {
        showNotification("No OTUs match this taxon.", type = "warning")
        return(NULL)
      }
      
      # Get relative abundances across all samples
      otu_mat <- otu_table(ps)
      rel_abund_per_sample <- colSums(otu_mat[matching_otus, , drop = FALSE])
      result_df <- data.frame(
        sample_id = names(rel_abund_per_sample),
        relative_abundance = as.numeric(rel_abund_per_sample),
        stringsAsFactors = FALSE
      ) %>%
        filter(relative_abundance > 0) %>%
        arrange(desc(relative_abundance))
      
      if (nrow(result_df) == 0) {
        showNotification("This taxon is absent from all samples.", type = "info")
        return(NULL)
      }
      
      return(result_df)
    })
    
    # Use observeEvent to handle the click event (which is often clearer than bindEvent for side effects)
    observeEvent(input$clicked_taxon, {
      req(input$clicked_taxon != "")
      
      # Get the data (this will trigger the eventReactive if needed)
      data_to_show <- compare_data_reactive()
      req(data_to_show) # Stop if no data
      
      # Create the table output *inside* the observeEvent
      output$compare_table <- renderDT({
        datatable(data_to_show, options = list(dom = 'Bfrtip', pageLength = 10)) %>%
          formatRound(columns = "relative_abundance", digits = 6)
      })
      
      # Show the modal
      showModal(modalDialog(
        title = paste("Samples containing:", input$clicked_taxon),
        size = "l",
        DTOutput(ns("compare_table")),
        footer = modalButton("Close")
      ))
    })
    
    output$sample_metadata <- renderUI({
      req(sample_metadata_text())
      HTML(sample_metadata_text())
    })
  })
}