# helpers/data_helpers.R

# Function to aggregate phyloseq object at specific taxonomic level
aggregate_phyloseq <- function(phyloseq_obj, level = "Genus") {
  # Check if the taxonomic level exists
  if (!level %in% rank_names(phyloseq_obj)) {
    stop("Taxonomic level not found in phyloseq object")
  }
  
  # Merge taxa at the specified level
  ps_agg <- tax_glom(phyloseq_obj, taxrank = level, NArm = FALSE)
  
  return(ps_agg)
}

# Function to convert phyloseq to data frame for plotting
ps_to_df <- function(phyloseq_obj, level = "Genus") {
  # Aggregate first
  ps_agg <- aggregate_phyloseq(phyloseq_obj, level)
  
  # Melt for plotting
  melted <- psmelt(ps_agg)
  
  return(melted)
}

# Function to get abundance table from phyloseq
get_abundance_table <- function(phyloseq_obj, level = "Genus") {
  ps_agg <- aggregate_phyloseq(phyloseq_obj, level)
  
  # Extract abundance table
  abun_table <- otu_table(ps_agg)
  tax_table <- tax_table(ps_agg)
  
  # Combine with taxonomy
  rownames(abun_table) <- tax_table[, level]
  
  # Transpose to have samples as rows
  abun_table <- as.data.frame(t(abun_table))
  abun_table$SampleID <- rownames(abun_table)
  
  return(abun_table)
}

# Function to filter phyloseq by experiment
filter_phyloseq_by_experiment <- function(phyloseq_obj, experiments) {
  ps_filtered <- subset_samples(phyloseq_obj, Experiment %in% experiments)
  return(ps_filtered)
}