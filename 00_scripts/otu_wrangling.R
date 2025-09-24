library(tidyverse)
library(readr)
library(purrr)
library(dplyr)
library(tidyr)
library(fs)
library(here)

rm(list = ls())

# PROBLEM 1 - otu data is contained in tab delimited text files, and not csv files.
# - this means there are no discrete columns, and these will in fact have to be deduced/calculated 
# - for each file when it is imported

# Set the path to your CSV files
file_path <- "00_otu_data/"

# Get list of text files (adjust pattern as needed)
text_files <- dir_ls(file_path, regexp = "\\.(txt|tsv|tab)$")

# Function to read tab-delimited files and standardize format
read_otu_text_file <- function(file) {
  # Read as tab-delimited
  read_delim(file, delim = "\t", show_col_types = FALSE) %>%
    # Ensure first column is named OTU (handles #OTU or other names)
    rename(OTU = 1) %>%
    mutate(OTU = as.character(OTU)) %>%
    # Convert to long format for merging
    pivot_longer(cols = -OTU, 
                 names_to = "Sample", 
                 values_to = "Count") %>%
    filter(!is.na(OTU), !is.na(Sample))
}

# 00_ Read and merge all files ----
merged_otu <- text_files %>%
  map(read_otu_text_file) %>%
  bind_rows() %>%
  group_by(OTU, Sample) %>%
  summarize(Count = sum(Count, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = Sample, 
              values_from = Count, 
              values_fill = 0)




# 01_ construct taxonomic table by breaking up OTU column ----

# 02_ construct otu_table ----

# 03_ construct sample_table ----