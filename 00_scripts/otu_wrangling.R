library(tidyverse)
library(readr)
library(purrr)
library(dplyr)
library(tidyr)
library(fs)
library(here)

rm(list = ls())

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

#_______________________________________________________________________________
# 00_ Read and merge all files ----
merged_otu <- text_files %>%
  map(read_otu_text_file) %>%
  bind_rows() %>%
  group_by(OTU, Sample) %>%
  summarize(Count = sum(Count, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = Sample, 
              values_from = Count, 
              values_fill = 0)

merged_otu <- merged_otu %>%
  rename_with(~ sub("^ESSE\\.", "", .x), .cols = starts_with("ESSE.222."))



#add in a unique identifier column
merged_otu <- merged_otu %>% 
  mutate(otu_id = row_number()) %>% 
  relocate(otu_id, .before = 1)



#_______________________________________________________________________________
#finding values for a specific species [this is for kyle]
spp_df <- merged_otu %>% filter(otu_id=='516')

df_long <- spp_df %>%
  select(-c("otu_id")) %>% 
  pivot_longer(
    cols = starts_with("222."),       # selects all sample columns
    names_to = "SampleID",            # new column for sample names
    values_to = "Abundance"           # new column for abundance values
  )

for_kyle <- df_long %>% filter(Abundance > 0.001) %>% 
  arrange(desc(Abundance)) 


#_______________________________________________________________________________
# 01_ construct taxonomic table by breaking up OTU column ----

split_taxa <- merged_otu %>%
  slice(-1) %>% 
  select(otu_id, OTU) 

#will need to retain a link btwn a OTU identifier and the taxonomic table
  
split_taxa <- split_taxa %>%   
separate(
    OTU,
    into = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"),
    sep = ";",
    remove = TRUE
  )

#clean up the values
split_taxa <- split_taxa %>%
  mutate(Kingdom = str_extract(Kingdom, "(?<=__).*")) %>% 
  mutate(Phylum = str_extract(Phylum, "(?<=__).*")) %>% 
  mutate(Class = str_extract(Class, "(?<=__).*")) %>% 
  mutate(Order = str_extract(Order, "(?<=__).*")) %>% 
  mutate(Family = str_extract(Family, "(?<=__).*")) %>% 
  mutate(Genus = str_extract(Genus, "(?<=__).*")) %>% 
  mutate(Species = str_extract(Species, "(?<=__).*")) 

for_tax_df <- as.data.frame(split_taxa)
for_tax_df_unique <- for_tax_df %>% 
  distinct()

for_tax_df_unique <- for_tax_df_unique %>%
  tibble::column_to_rownames("otu_id")

write_csv(for_tax_df_unique, "00_cleaned_data/tax_table.csv")

tax_mat <- as.matrix(for_tax_df_unique)
tax_phy <- phyloseq::tax_table(tax_mat)  # for use in constructing our phyloseq obj

#_______________________________________________________________________________
# 02_ construct otu_table ----
merged_otu <- merged_otu %>% 
  select(-OTU)
otu_df <- merged_otu %>%
  tibble::column_to_rownames("otu_id")
otu_mat <- as.matrix(otu_df)
class(otu_mat) <- "numeric"
otu_rel <- phyloseq::otu_table(otu_mat, taxa_are_rows = TRUE) # for use in 
                                                # constructing our phyloseq obj

# 03_ construct sample_table ----