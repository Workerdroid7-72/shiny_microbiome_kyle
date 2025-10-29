# global.R
# Load required packages
library(shiny)
library(bslib)
library(phyloseq)
library(dplyr)
library(tibble)
library(plotly)
library(ggplot2)
library(DT)
library(readr)
library(data.table)

# Define taxonomic levels
TAX_LEVELS <- c("Order" = "Order", "Family" = "Family", "Genus" = "Genus")

# Define your phyloseq object (replace with your actual object)
ps <- readRDS("01_cleaned_data/phyloseq_object.rds")
our_sample_table <- phyloseq::sample_data(ps)
our_taxonomic_table <- phyloseq::tax_table(ps)
our_otu_table <- phyloseq::otu_table(ps)
