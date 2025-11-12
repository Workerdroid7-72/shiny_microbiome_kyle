# sample_compare.R

sample_compare_ui <- function(id) {
  ns <- NS(id)
  
  layout_sidebar(
    sidebar = sidebar(
      width = 300,
      textInput(ns("sample_id1"), "Before Sample (3 or 4 digits):", value = ""),
      textInput(ns("sample_id2"), "After Sample (3 or 4 digits):", value = ""),
      
      selectInput(
        ns("tax_level"), 
        "Taxonomic Level:", 
        choices = c("Family", "Genus", "Species")
      ),
      selectInput(
        ns("top_n"), 
        "Show top taxa (combined):", 
        choices = c("Top 5" = 5, "Top 10" = 10, "Top 20" = 20),
        selected = 10
      ),
      actionButton(ns("go"), "Compare Samples", class = "btn-primary"),
      hr(),
      htmlOutput(ns("sample1_metadata")),
      hr(),
      htmlOutput(ns("sample2_metadata"))
    ),
    plotOutput(ns("paired_plot"), height = 400),
    DTOutput(ns("compare_table"))
  )
}

sample_compare_server <- function(id, ps) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns
    
    cat("📊 sample_explorer: otu_table sum =", sum(otu_table(ps)), 
        "| sample 0012 sum =", sum(otu_table(ps)[, "ESSE.222.0012"]), "\n")
    
    current_tax_level_reactive <- reactiveVal(NULL)
    
    # Helper: validate and expand short ID
    expand_sample_id <- function(short_id) {
      if (!grepl("^[0-9]{3,4}$", short_id)) {
        return(NULL)
      }
      paste0("ESSE.222.", short_id)
    }
    
    # Fetch taxa data for a single sample (returns df with tax_level col + rel_abund)
    get_sample_taxa_df <- function(full_id, tax_level) {
      # ✅ Copy from sample_explorer — proven to work
      sample_otu <- otu_table(ps)[, full_id, drop = FALSE]  # keep as 1-col matrix
      tax_df <- as.data.frame(tax_table(ps))
      
      df <- data.frame(
        taxon = rownames(sample_otu),
        relative_abundance = as.numeric(sample_otu[, 1]),  # safe: 1 column
        stringsAsFactors = FALSE
      )
      df <- cbind(df, tax_df[match(df$taxon, rownames(tax_df)), , drop = FALSE])
      df[is.na(df)] <- ""
      
      # Rest identical to sample_explorer...
      if (tax_level == "Species") {
        df$display_label <- trimws(paste(df$Genus, df$Species))
        df$display_label[df$display_label == ""] <- "Unknown"
        group_var <- "display_label"
      } else {
        group_var <- tax_level
        df$display_label <- df[[group_var]]
      }
      
      df_agg <- df %>%
        dplyr::group_by(display_label) %>%
        dplyr::summarise(relative_abundance = sum(relative_abundance), .groups = "drop") %>%
        dplyr::filter(display_label != "" & !is.na(display_label)) %>%
        dplyr::arrange(dplyr::desc(relative_abundance))
      
      df_agg
    }
    
    # Reactive: compute aligned taxon table for both samples
    paired_taxa <- eventReactive(input$go, {
      cat("🔍 [DEBUG] Button 'go' clicked at", Sys.time(), "\n")
      cat("  Input IDs: '", input$sample_id1, "' | '", input$sample_id2, "'\n", sep = "")
      
      # Validate and expand IDs
      expand_id <- function(short) {
        if (!grepl("^[0-9]{3,4}$", short)) return(NULL)
        paste0("ESSE.222.", short)
      }
      
      id1 <- expand_id(input$sample_id1)
      id2 <- expand_id(input$sample_id2)
      
      cat("  Expanded IDs: '", id1, "' | '", id2, "'\n", sep = "")
      
      if (is.null(id1) || is.null(id2)) {
        showNotification("Please enter 3 or 4 digits for both samples.", type = "error")
        return(NULL)
      }
      
      if (!id1 %in% sample_names(ps)) {
        showNotification(paste("Sample 1 (", id1, ") not found.", sep = ""), type = "error")
        return(NULL)
      }
      if (!id2 %in% sample_names(ps)) {
        showNotification(paste("Sample 2 (", id2, ") not found.", sep = ""), type = "error")
        return(NULL)
      }
      
      cat(" ✅ IDs valid and found. Proceeding...\n")
      
      # 🔁 Reuse sample_explorer's aggregation logic (proven to work)
      get_taxa_df <- function(full_id) {
        # Extract OTU data for sample
        sample_otu <- otu_table(ps)[, full_id, drop = FALSE]
        tax_df <- as.data.frame(tax_table(ps))
        
        df <- data.frame(
          taxon = rownames(sample_otu),
          relative_abundance = as.numeric(sample_otu[, 1]),
          stringsAsFactors = FALSE
        )
        # Merge taxonomy
        df <- cbind(df, tax_df[match(df$taxon, rownames(tax_df)), , drop = FALSE])
        df[is.na(df)] <- ""
        
        # Build display label
        if (input$tax_level == "Species") {
          df$display_label <- trimws(paste(df$Genus, df$Species))
          df$display_label[df$display_label == ""] <- "Unknown"
        } else {
          # Safely get tax level column (case-insensitive fallback)
          tax_col <- input$tax_level
          if (!tax_col %in% colnames(df)) {
            # Try lowercase
            lower_match <- grep(tolower(tax_col), tolower(colnames(df)), value = TRUE)
            if (length(lower_match) > 0) {
              tax_col <- lower_match[1]
            } else {
              warning("Taxonomic level '", input$tax_level, "' not found.")
              df$display_label <- "Unknown"
            }
          }
          if (tax_col %in% colnames(df)) {
            df$display_label <- df[[tax_col]]
          }
          df$display_label[is.na(df$display_label) | df$display_label == ""] <- "Unknown"
        }
        
        # Aggregate
        df_agg <- df %>%
          dplyr::group_by(display_label) %>%
          dplyr::summarise(relative_abundance = sum(relative_abundance), .groups = "drop") %>%
          dplyr::filter(display_label != "Unknown") %>%  # keep Unknown if you prefer
          dplyr::arrange(dplyr::desc(relative_abundance))
        
        df_agg
      }
      
      # Get data for both samples
      df1 <- get_taxa_df(id1)
      df2 <- get_taxa_df(id2)
      
      cat("   df1: ", nrow(df1), "taxa | df2: ", nrow(df2), "taxa\n")
      if (nrow(df1) == 0 && nrow(df2) == 0) {
        showNotification("No taxa found for either sample at this level.", type = "info")
        return(NULL)
      }
      
      # Merge on union of taxa
      all_taxa <- union(df1$display_label, df2$display_label)
      result <- tibble::tibble(display_label = all_taxa) %>%
        dplyr::left_join(
          df1 %>% dplyr::rename(rel_abund_1 = relative_abundance),
          by = "display_label"
        ) %>%
        dplyr::left_join(
          df2 %>% dplyr::rename(rel_abund_2 = relative_abundance),
          by = "display_label"
        ) %>%
        dplyr::mutate(
          rel_abund_1 = dplyr::coalesce(rel_abund_1, 0),
          rel_abund_2 = dplyr::coalesce(rel_abund_2, 0),
          difference = rel_abund_2 - rel_abund_1
        ) %>%
        dplyr::arrange(dplyr::desc(rel_abund_1)) %>%
        head(as.numeric(input$top_n))
      
      # Rename display_label to tax_level for consistency with sample_explorer
      result <- result %>%
        dplyr::rename(!!sym(input$tax_level) := display_label)
      
      cat(" ✅ Final result: ", nrow(result), " taxa. Columns: ", paste(colnames(result), collapse = ", "), "\n")
      
      # Attach sample IDs for plot title
      attr(result, "sample_ids") <- c(id1, id2)
      result
    })
    
    # Metadata helpers
    get_metadata_html <- function(full_id) {
      meta_df <- as.data.frame(sample_data(ps))
      rownames(meta_df) <- sample_names(ps)
      if (!full_id %in% rownames(meta_df)) return(NULL)
      
      row_meta <- meta_df[full_id, , drop = FALSE]
      get_val <- function(col) {
        if (col %in% names(row_meta)) {
          val <- as.character(row_meta[[col]])
          if (is.na(val) || val == "" || val == "NA") "(none)" else val
        } else {
          "(column missing)"
        }
      }
      glue::glue(
        "<strong>{full_id}</strong><br/>",
        "<em>Experiment:</em> {get_val('experiment_name')}<br/>",
        "<em>Notes:</em> {get_val('notes')}<br/>",
        "<em>Seq Type:</em> {get_val('seq_type')}"
      )
    }
    
    metadata1 <- eventReactive(input$go, {
      id1 <- expand_sample_id(input$sample_id1)
      if (is.null(id1)) return(NULL)
      get_metadata_html(id1)
    })
    
    metadata2 <- eventReactive(input$go, {
      id2 <- expand_sample_id(input$sample_id2)
      if (is.null(id2)) return(NULL)
      get_metadata_html(id2)
    })
    
    # Plot: paired bars (dodged)
    output$paired_plot <- renderPlot({
      data <- paired_taxa()
      cat("📈 paired_plot called. paired_taxa() is", class(data), "| nrows =", 
          ifelse(is.null(data), "NULL", nrow(data)), "\n")
      
      if (is.null(data)) {
        plot(1, type = "n", xlab = "", ylab = "", main = "No data")
        text(1, 1, "paired_taxa() returned NULL", col = "red")
        return()
      }
      
      if (nrow(data) == 0) {
        plot(1, type = "n", xlab = "", ylab = "", main = "Empty data")
        text(1, 1, "paired_taxa() has 0 rows", col = "orange")
        return()
      }
      
      req(paired_taxa())
      df <- paired_taxa()
      tax_level <- attr(df, "sample_ids", exact = TRUE)
      if (is.null(tax_level)) return(NULL)
      
      sample_ids <- attr(df, "sample_ids")
      
      # Melt for ggplot
      plot_df <- df %>%
        select(!!sym(input$tax_level), rel_abund_1, rel_abund_2) %>%
        tidyr::pivot_longer(
          cols = c(rel_abund_1, rel_abund_2),
          names_to = "sample",
          values_to = "relative_abundance"
        ) %>%
        mutate(
          sample = factor(sample,
                          levels = c("rel_abund_1", "rel_abund_2"),
                          labels = c(sample_ids[1], sample_ids[2]))
        ) %>%
        mutate(
          !!sym(input$tax_level) := reorder(!!sym(input$tax_level), -relative_abundance)
        )
      
      ggplot(plot_df, aes(
        x = !!sym(input$tax_level),
        y = relative_abundance,
        fill = sample
      )) +
        geom_col(position = "dodge", alpha = 0.85, width = 0.75) +
        geom_text(
          aes(label = ifelse(
            relative_abundance > 0.05,  # higher threshold for inside
            sprintf("%.2f", relative_abundance),
            ""
          )),
          position = position_dodge(width = 0.75),
          hjust = 1.0,           # right-align inside bar
          vjust = 0.5,
          size = 3,
          color = "white",
          fontface = "bold"
        ) +
        # Add a slight margin to prevent cutoff
        scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +  # extra space on right
        coord_flip() +
        scale_fill_manual(values = c("#45B29D", "#3f7be2")) +
        labs(
          x = ifelse(input$tax_level == "Species", "Genus Species", input$tax_level),
          y = "Relative Abundance",
          title = paste("Paired Taxa Abundance:",
                        sample_ids[1], "vs", sample_ids[2]),
          fill = "Sample"
        ) +
        theme_minimal(base_size = 10) +
        theme(
          legend.position = "top",
          axis.text = element_text(size = 9)
        )
    })
    
    # DT: with color-coded difference
    # output$compare_table <- renderDT({
    #   cat("🧾 compare_table: entering renderDT\n")
    #   
    #   data <- paired_taxa()
    #   cat("   paired_taxa() returned:", 
    #       if (is.null(data)) "NULL" else paste(class(data), "with", nrow(data), "rows"), "\n")
    #   
    #   if (is.null(data)) {
    #     return(datatable(data.frame(Error = "No data available."), 
    #                      options = list(dom = 't')))
    #   }
    #   
    #   if (nrow(data) == 0) {
    #     return(datatable(data.frame(Warning = "No taxa to display."), 
    #                      options = list(dom = 't')))
    #   }
    #   
    #   # Proceed only if valid
    #   tax_level <- input$tax_level
    #   colnames(data)[colnames(data) == tax_level] <- "Taxon"
    #   
    #   df_clean <- data %>%
    #     dplyr::rename(
    #       `Sample 1` = rel_abund_1,
    #       `Sample 2` = rel_abund_2,
    #       Difference = difference
    #     ) %>%
    #     dplyr::select(Taxon, `Sample 1`, `Sample 2`, Difference) %>%
    #     dplyr::mutate(dplyr::across(dplyr::where(is.numeric), ~ round(., 6)))
    #   
    #   cat("   Building DT with", nrow(df_clean), "rows, cols:", colnames(df_clean), "\n")
    #   
    #   # Ensure df_clean is a data.frame
    #   if (!is.data.frame(df_clean)) {
    #     cat("   ❌ df_clean is not a data.frame! Class:", class(df_clean), "\n")
    #     df_clean <- as.data.frame(df_clean)
    #   }
    #   
    #   dt_obj <- datatable(
    #     df_clean,
    #     options = list(pageLength = min(10, nrow(df_clean)), dom = 't'),
    #     class = "compact"
    #   )
    #   
    #   cat("   DT object created. Applying formatting...\n")
    #   
    #   dt_obj %>%
    #     formatStyle(
    #       "Difference",
    #       backgroundColor = styleColorBar(colour = c("#E74C3C", "#2ECC71")),
    #       color = "white",
    #       fontWeight = "bold"
    #     ) %>%
    #     formatRound(columns = c("Sample 1", "Sample 2", "Difference"), digits = 6)
    # })
    
    output$compare_table <- renderDT({
      data <- paired_taxa()
      req(data)
      
      tax_level <- input$tax_level
      colnames(data)[colnames(data) == tax_level] <- "Taxon"
      
      df_clean <- data %>%
        dplyr::rename(
          `Sample 1` = rel_abund_1,
          `Sample 2` = rel_abund_2,
          Difference = difference
        ) %>%
        dplyr::select(Taxon, `Sample 1`, `Sample 2`, Difference) %>%
        dplyr::mutate(dplyr::across(dplyr::where(is.numeric), ~ round(., 4)))
      
      datatable(df_clean, options = list(dom = 't', pageLength = min(10, nrow(df_clean)))) %>%
        formatStyle(
          "Difference",
          backgroundColor = styleInterval(c(0), c("#E74C3C", "#2ECC71")),
          color = "white",
          fontWeight = "bold"
        ) %>%
        formatRound(columns = c("Sample 1", "Sample 2", "Difference"), digits = 4)
    })
    
    # Metadata outputs
    output$sample1_metadata <- renderUI({ HTML(metadata1()) })
    output$sample2_metadata <- renderUI({ HTML(metadata2()) })
    
  })
}