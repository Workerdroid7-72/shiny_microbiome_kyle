# global.R
# Load required packages
library(shiny)
library(bslib)
library(phyloseq)
library(dplyr)
library(tidyr)
library(ggplot2)
library(DT)


# Define your phyloseq object (replace with your actual object)
ps <- readRDS("00_data/03_cleaned_data/phyloseq_object.rds")
stopifnot(nsamples(ps) > 0, ntaxa(ps) > 0) # checking to make sure the ps object is correct

meta_df <- as.data.frame(sample_data(ps))
tax_table_df <- as.data.frame(tax_table(ps))
unique_family <- sort(unique(tax_table_df$Family))
unique_genus <- sort(unique(tax_table_df$Genus))
unique_species <- sort(unique(tax_table_df$Species))

# Remove NA or blank entries
unique_genus <- unique_genus[!is.na(unique_genus) & unique_genus != ""]
unique_species <- unique_species[!is.na(unique_species) & unique_species != ""]
