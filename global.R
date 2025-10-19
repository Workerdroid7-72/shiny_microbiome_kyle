# global.R
library(shiny)
library(bslib)
library(tidyverse)
library(DT)
library(plotly)
library(phyloseq)
library(reshape2)


# Load the phyloseq object
#ps <- load_phyloseq_data()

# Global variables
taxonomic_levels <- c("Order", "Family", "Genus")
experiments <- unique(sample_data(ps)$Experiment)

# Source helper functions
source("01_helpers/data_helpers.R")

# Source modules
source("01_modules/individual_samples.R")
source("01_modules/experiment_analysis.R")
source("01_modules/data_download.R")