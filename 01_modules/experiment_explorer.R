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
      actionButton(ns("apply_btn"), "Show Results", class = "btn-primary"),
      uiOutput(ns("download_ui")),  # dynamic button
      hr(),
      htmlOutput(ns("sample_metadata")),
      
      # Hidden input to receive clicked taxon
      div(
        textInput(ns("clicked_taxon"), label = NULL, value = ""),
        style = "display: none;"
      )
    ),
    
    # NEW: Tabbed interface for visualizations
    navset_card_tab(
      title = "Visualization",
      nav_panel("Bar", plotOutput(ns("plot_bar"), height = 400)),
      nav_panel("Pie", plotOutput(ns("plot_pie"), height = 400)),
      #nav_panel("Treemap", plotOutput(ns("plot_treemap"), height = 400))
      nav_panel("Stacked", plotOutput(ns("plot_stacked"), height = 150))
    ),
    DTOutput(ns("taxa_table"))                      
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
      
      sample_df <- as.data.frame(sample_data(ps), stringsAsFactors = FALSE)
      sample_df$experiment_name <- trimws(as.character(sample_df$experiment_name))
      selected_exp <- trimws(as.character(input$experiment))
      
      samples_in_exp <- rownames(sample_df)[sample_df$experiment_name == selected_exp]
      
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
    
    # NEW: Reactive values to store the last *clicked* sample (for reset logic)
    last_clicked_sample <- reactiveVal(NULL)
    
    # NEW: Reactive values to store results (instead of eventReactive)
    results <- reactiveValues(
      plot_data = NULL,
      table_data = NULL,
      sample_id = NULL,
      tax_level = NULL,
      metadata_html = NULL
    )
    
    # NEW: Reset outputs when key inputs change (prevents stale data)
    observeEvent({
      input$experiment
      input$sample
    }, {
      # Reset all results when inputs change
      results$plot_data <- NULL
      results$table_data <- NULL
      results$sample_id <- NULL
      results$tax_level <- NULL
      results$metadata_html <- NULL
      last_clicked_sample(NULL)
    }, ignoreNULL = FALSE)
    
    # NEW: Process data when button is clicked
    observeEvent(input$apply_btn, {
      req(input$tax_level, input$top_n, input$experiment, input$sample)
      
      tax_level <- input$tax_level
      top_n <- as.numeric(input$top_n)
      sample_id <- input$sample
      
      # Remember which sample was clicked
      last_clicked_sample(sample_id)
      
      # Wrap everything in tryCatch to catch errors
      tryCatch({
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
        
        sample_otu <- otu_table(ps)[, sample_id, drop = FALSE]
        tax_df <- as.data.frame(tax_table(ps))
        
        df <- data.frame(
          taxon = rownames(sample_otu),
          relative_abundance = as.numeric(sample_otu[, 1]),
          stringsAsFactors = FALSE
        )
        df <- cbind(df, tax_df[match(df$taxon, rownames(tax_df)), , drop = FALSE])
        df[is.na(df)] <- ""
        
        if (tax_level == "Species") {
          df$display_label <- trimws(paste(df$Genus, df$Species))
          df$display_label[df$display_label == ""] <- "Unknown"
          group_var <- "display_label"
        } else {
          group_var <- tax_level
          df$display_label <- df[[tax_level]]
        }
        
        df_agg <- df %>%
          dplyr::filter(relative_abundance > 0) %>%
          dplyr::group_by(!!sym("display_label")) %>%
          dplyr::summarise(relative_abundance = sum(relative_abundance), .groups = "drop") %>%
          dplyr::filter(display_label != "" & display_label != "Unknown") %>%
          dplyr::arrange(dplyr::desc(relative_abundance))
        
        df_agg <- head(df_agg, top_n)
        colnames(df_agg)[colnames(df_agg) == "display_label"] <- "TaxonLabel"
        
        # Store results
        results$plot_data <- df_agg
        results$table_data <- df_agg
        results$sample_id <- sample_id
        results$tax_level <- tax_level
        
        # Also process metadata
        full_id <- input$sample
        if (full_id %in% sample_names(ps)) {
          meta_df <- as.data.frame(sample_data(ps))
          rownames(meta_df) <- sample_names(ps)
          
          if (full_id %in% rownames(meta_df)) {
            row_meta <- meta_df[full_id, , drop = FALSE]
            
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
            
            results$metadata_html <- HTML(glue::glue(
              "<strong>Experiment:</strong> {exp_name}<br/>",
              "<strong>Notes:</strong> {notes_val}<br/>",
              "<strong>Seq Type:</strong> {seq_type_val}"
            ))
          }
        }
      }, error = function(e) {
        showModal(modalDialog(
          title = "Error",
          paste("An error occurred:", conditionMessage(e)),
          easyClose = TRUE
        ))
      })
    })
    
    # 3️⃣ Render metadata
    output$sample_metadata <- renderUI({
      if (is.null(last_clicked_sample())) {
        HTML("<em>Select a sample and click 'Show Results' to display metadata.</em>")
      } else {
        req(results$metadata_html)
        results$metadata_html
      }
    })
    
    # NEW: Bar plot renderer
    output$plot_bar <- renderPlot({
      if (is.null(results$plot_data)) {
        return(NULL)  # Clear plot
      }
      
      df_clean <- results$plot_data
      sample_id <- results$sample_id
      tax_level <- results$tax_level
      
      if (nrow(df_clean) == 0) {
        return(ggplot2::ggplot() + 
                 ggplot2::annotate("text", x = 0.5, y = 0.5, label = "No data to display", size = 5) +
                 ggplot2::theme_void())
      }
      
      df_clean <- df_clean %>%
        dplyr::mutate(
          label = round(relative_abundance, digits = 3),
          TaxonLabel_ordered = reorder(TaxonLabel, relative_abundance)
        )
      
      p <- ggplot2::ggplot(df_clean, ggplot2::aes(
        x = TaxonLabel_ordered,
        y = relative_abundance,
        label = label
      )) +
        ggplot2::geom_col(fill = "#2c7bb6") +
        ggplot2::geom_text(
          ggplot2::aes(label = ifelse(relative_abundance >= 0.03, label, "")),
          hjust = 1.05,
          color = "white",
          fontface = "bold",
          size = 3,
          na.rm = TRUE
        ) +
        ggplot2::coord_flip() +
        ggplot2::labs(
          x = "Taxon",
          y = "Relative Abundance",
          title = paste("Top", nrow(df_clean), tax_level, "in Sample:", sample_id)
        ) +
        ggplot2::theme_minimal(base_size = 10) +
        ggplot2::scale_x_discrete(expand = ggplot2::expansion(mult = c(0, 0.15)))
      
      print(p)
    })
    
    # NEW: Pie chart renderer
    output$plot_pie <- renderPlot({
      if (is.null(results$plot_data)) {
        return(NULL)  # Clear plot
      }
      
      df_clean <- results$plot_data
      sample_id <- results$sample_id
      tax_level <- results$tax_level
      
      if (nrow(df_clean) == 0) {
        return(ggplot2::ggplot() + 
                 ggplot2::annotate("text", x = 0.5, y = 0.5, label = "No data to display", size = 5) +
                 ggplot2::theme_void())
      }
      
      # Create a pie chart using coord_polar
      p <- ggplot2::ggplot(df_clean, ggplot2::aes(
        x = factor(1),  # Single bar
        y = relative_abundance,
        fill = TaxonLabel
      )) +
        ggplot2::geom_bar(stat = "identity", width = 1) +
        ggplot2::coord_polar(theta = "y") +
        ggplot2::labs(
          title = paste("Top", nrow(df_clean), tax_level, "in Sample:", sample_id),
          fill = tax_level
        ) +
        ggplot2::theme_void() +
        ggplot2::theme(
          axis.title = element_blank(),
          axis.text = element_blank(),
          axis.ticks = element_blank(),
          plot.title = element_text(hjust = 0.5)
        )
      
      print(p)
    })
    
    # NEW: Stacked bar renderer (horizontal, compact)
    # output$plot_stacked <- renderPlot({
    #   if (is.null(results$plot_data)) {
    #     return(NULL)  # Clear plot
    #   }
    #   
    #   df_clean <- results$plot_data
    #   sample_id <- results$sample_id
    #   tax_level <- results$tax_level
    #   
    #   if (nrow(df_clean) == 0) {
    #     return(ggplot2::ggplot() + 
    #              ggplot2::annotate("text", x = 0.5, y = 0.5, label = "No data to display", size = 5) +
    #              ggplot2::theme_void())
    #   }
    #   
    #   # For stacked bar, we'll use the same data but format differently
    #   p <- ggplot2::ggplot(df_clean, ggplot2::aes(
    #     x = factor(1),  # Single x position
    #     y = relative_abundance,
    #     fill = TaxonLabel
    #   )) +
    #     ggplot2::geom_bar(stat = "identity", width = 0.5) +
    #     ggplot2::coord_flip() +  # Horizontal orientation
    #     ggplot2::labs(
    #       title = paste("Composition -", tax_level, "Level"),
    #       y = "Relative Abundance",
    #       fill = tax_level
    #     ) +
    #     ggplot2::theme_minimal() +
    #     ggplot2::theme(
    #       axis.title.y = element_blank(),
    #       axis.text.y = element_blank(),
    #       axis.ticks.y = element_blank(),
    #       plot.title = element_text(hjust = 0.5)
    #     )
    #   
    #   print(p)
    # })
    
    output$plot_stacked <- renderPlot({
      if (is.null(results$plot_data)) {
        return(NULL)  # Clear plot
      }
      
      df_clean <- results$plot_data
      sample_id <- results$sample_id
      tax_level <- results$tax_level
      
      if (nrow(df_clean) == 0) {
        return(ggplot2::ggplot() + 
                 ggplot2::annotate("text", x = 0.5, y = 0.5, label = "No data to display", size = 5) +
                 ggplot2::theme_void())
      }
      
      # Generate colors: use a palette that can handle many categories
      n_taxa <- nrow(df_clean)
      if (n_taxa <= 12) {
        # Use brewer palette for smaller numbers
        colors <- RColorBrewer::brewer.pal(
          n = if (n_taxa <= 3) 3 else if (n_taxa <= 8) 8 else 12,
          name = if (n_taxa <= 8) "Set2" else "Paired"
        )[1:n_taxa]
      } else {
        # Use viridis for larger numbers (works well with many categories)
        colors <- viridisLite::viridis(n_taxa)
      }
      
      # For stacked bar, we'll use the same data but format differently
      p <- ggplot2::ggplot(df_clean, ggplot2::aes(
        x = factor(1),  # Single x position
        y = relative_abundance,
        fill = TaxonLabel
      )) +
        ggplot2::geom_bar(
          stat = "identity", 
          width = 0.5,
          color = "white",  # White lines between sections
          size = 0.5        # Line thickness
        ) +
        ggplot2::coord_flip() +  # Horizontal orientation
        ggplot2::scale_fill_manual(values = colors) +
        ggplot2::labs(
          title = paste("Composition -", tax_level, "Level"),
          y = "Relative Abundance",
          fill = tax_level
        ) +
        ggplot2::theme_minimal() +
        ggplot2::theme(
          axis.title.y = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          plot.title = element_text(hjust = 0.5)
        )
      
      print(p)
    })
    
    # 5️⃣ Render table
    output$taxa_table <- DT::renderDT({
      if (is.null(results$table_data)) {
        return(DT::datatable(data.frame()))  # Clear table
      }
      
      df_clean <- results$table_data
      
      DT::datatable(
        df_clean,
        options = list(
          dom = 'Bfrtip',
          lengthChange = FALSE
        ),
        rownames = FALSE
      ) %>%
        DT::formatRound(
          columns = "relative_abundance", 
          digits = 4,
          mark = ""
        )
    })
    
  })
}

# experiment_explorer_server <- function(id, ps) {
#   moduleServer(id, function(input, output, session) {
#     ns <- session$ns
#     
#     # Ensure 'experiment' column exists
#     if (!"experiment" %in% names(sample_data(ps))) {
#       stop("Your phyloseq object's sample_data must contain a column named 'experiment'.")
#     }
#     
#     # Populate experiment dropdown once (ensure character)
#     experiments <- sort(unique(as.character(sample_data(ps)$experiment_name)))
#     updateSelectInput(session, "experiment", choices = experiments)
#     
#     # Render UI for sample selector (static structure)
#     output$sample_ui <- renderUI({
#       selectInput(ns("sample"),
#                   "Select Sample",
#                   choices = character(0),
#                   multiple = FALSE)
#     })
#     
#     # Update sample choices when experiment_name changes
#     observeEvent(input$experiment, {
#       req(input$experiment)
#       
#       sample_df <- as.data.frame(sample_data(ps), stringsAsFactors = FALSE)
#       sample_df$experiment_name <- trimws(as.character(sample_df$experiment_name))
#       selected_exp <- trimws(as.character(input$experiment))
#       
#       samples_in_exp <- rownames(sample_df)[sample_df$experiment_name == selected_exp]
#       
#       if (length(samples_in_exp) == 0) {
#         updateSelectInput(session,
#                           "sample",
#                           choices = character(0),
#                           selected = character(0))
#       } else {
#         updateSelectInput(
#           session,
#           "sample",
#           choices = sort(samples_in_exp),
#           selected = samples_in_exp[1]
#         )
#       }
#     })
#     
#     # 1️⃣ NEW: Reactive to store sample metadata (only when button is pressed)
#     sample_metadata_html <- eventReactive(input$apply_btn, {
#       req(input$sample)
#       
#       full_id <- input$sample
#       if (!full_id %in% sample_names(ps)) {
#         return(HTML("<em>No metadata found for sample.</em>"))
#       }
#       
#       meta_df <- as.data.frame(sample_data(ps))
#       rownames(meta_df) <- sample_names(ps)
#       
#       if (!full_id %in% rownames(meta_df)) {
#         return(HTML("<em>Sample not found in metadata.</em>"))
#       }
#       
#       row_meta <- meta_df[full_id, , drop = FALSE]
#       
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
#       HTML(glue::glue(
#         "<strong>Experiment:</strong> {exp_name}<br/>",
#         "<strong>Notes:</strong> {notes_val}<br/>",
#         "<strong>Seq Type:</strong> {seq_type_val}"
#       ))
#     })
#     
#     
#     
#     # 2️⃣ NEW: Reactive to store taxa processing results (only when button is pressed)
#     processed_data <- eventReactive(input$apply_btn, {
#       req(input$tax_level, input$top_n, input$experiment, input$sample)
#       
#       tax_level <- input$tax_level
#       top_n <- as.numeric(input$top_n)
#       sample_id <- input$sample
#       
#       # Wrap everything in tryCatch to catch errors
#       tryCatch({
#         if (!tax_level %in% rank_names(ps)) {
#           showModal(
#             modalDialog(
#               title = "Error",
#               "Selected taxonomic level not found in phyloseq taxonomy table.",
#               easyClose = TRUE
#             )
#           )
#           return(NULL)
#         }
#         
#         sample_otu <- otu_table(ps)[, sample_id, drop = FALSE]
#         tax_df <- as.data.frame(tax_table(ps))
#         
#         df <- data.frame(
#           taxon = rownames(sample_otu),
#           relative_abundance = as.numeric(sample_otu[, 1]),
#           stringsAsFactors = FALSE
#         )
#         df <- cbind(df, tax_df[match(df$taxon, rownames(tax_df)), , drop = FALSE])
#         df[is.na(df)] <- ""
#         
#         if (tax_level == "Species") {
#           df$display_label <- trimws(paste(df$Genus, df$Species))
#           df$display_label[df$display_label == ""] <- "Unknown"
#           group_var <- "display_label"
#         } else {
#           group_var <- tax_level
#           df$display_label <- df[[tax_level]]
#         }
#         
#         df_agg <- df %>%
#           dplyr::filter(relative_abundance > 0) %>%
#           dplyr::group_by(!!sym("display_label")) %>%
#           dplyr::summarise(relative_abundance = sum(relative_abundance), .groups = "drop") %>%
#           dplyr::filter(display_label != "" & display_label != "Unknown") %>%
#           dplyr::arrange(dplyr::desc(relative_abundance))
#         
#         df_agg <- head(df_agg, top_n)
#         colnames(df_agg)[colnames(df_agg) == "display_label"] <- "TaxonLabel"
#         
#         # Return a list of both plot and table data
#         list(
#           plot_data = df_agg,
#           table_data = df_agg,
#           sample_id = sample_id,
#           tax_level = tax_level
#         )
#       }, error = function(e) {
#         showModal(modalDialog(
#           title = "Error",
#           paste("An error occurred:", conditionMessage(e)),
#           easyClose = TRUE
#         ))
#         return(NULL)
#       })
#     })
#     
#     
#     
#     # 3️⃣ Render metadata
#     output$sample_metadata <- renderUI({
#       sample_metadata_html()
#     })
#     
#     # 4️⃣ Render plot
#     output$taxa_plot <- renderPlot({
#       req(processed_data())  # Wait for data to be ready
#       
#       data <- processed_data()
#       df_clean <- data$plot_data
#       sample_id <- data$sample_id
#       tax_level <- data$tax_level
#       
#       if (nrow(df_clean) == 0) {
#         return(ggplot2::ggplot() + 
#                  ggplot2::annotate("text", x = 0.5, y = 0.5, label = "No data to display", size = 5) +
#                  ggplot2::theme_void())
#       }
#       
#       df_clean <- df_clean %>%
#         dplyr::mutate(
#           label = round(relative_abundance, digits = 3),
#           TaxonLabel_ordered = reorder(TaxonLabel, relative_abundance)
#         )
#       
#       p <- ggplot2::ggplot(df_clean, ggplot2::aes(
#         x = TaxonLabel_ordered,
#         y = relative_abundance,
#         label = label
#       )) +
#         ggplot2::geom_col(fill = "#2c7bb6") +
#         ggplot2::geom_text(
#           ggplot2::aes(label = ifelse(relative_abundance >= 0.03, label, "")),
#           hjust = 1.05,
#           color = "white",
#           fontface = "bold",
#           size = 3,
#           na.rm = TRUE
#         ) +
#         ggplot2::coord_flip() +
#         ggplot2::labs(
#           x = "Taxon",
#           y = "Relative Abundance",
#           title = paste("Top", nrow(df_clean), tax_level, "in Sample:", sample_id)
#         ) +
#         ggplot2::theme_minimal(base_size = 10) +
#         ggplot2::scale_x_discrete(expand = ggplot2::expansion(mult = c(0, 0.15)))
#       
#       print(p)
#     })
#     
#     # 5️⃣ Render table
#     output$taxa_table <- DT::renderDT({
#       req(processed_data())
#       
#       df_clean <- processed_data()$table_data
#       
#       DT::datatable(
#         df_clean,
#         options = list(
#           dom = 'Bfrtip',
#           lengthChange = FALSE
#         ),
#         rownames = FALSE
#       ) %>%
#         DT::formatRound(
#           columns = "relative_abundance", 
#           digits = 4,
#           mark = ""
#         )
#     })
#     
#   })
# }