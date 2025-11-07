library(tidyverse)
library(readr)
library(purrr)
library(dplyr)
library(tidyr)
library(fs)
library(here)
library(stringr)

rm(list = ls())

# Set the path to your CSV files
file_path <- "00_data/01_otu_data/"

# Get list of text files (adjust pattern as needed)
text_files <- dir_ls(file_path, regexp = "_spp.*\\.(txt|tab)$")
text_files2 <- dir_ls(file_path, regexp = "_spp.*\\.(tsv)$")

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
merged_otu_16s <- text_files %>%
  map(read_otu_text_file) %>%
  bind_rows() %>%
  group_by(OTU, Sample) %>%
  summarize(Count = sum(Count, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = Sample, 
              values_from = Count, 
              values_fill = 0)

merged_otu_shotgun <- text_files2 %>%
  map(read_otu_text_file) %>%
  bind_rows() %>%
  group_by(OTU, Sample) %>%
  summarize(Count = sum(Count, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = Sample, 
              values_from = Count, 
              values_fill = 0)

#make sure sample swab numbers have a CONSISTENT naming convention
merged_otu_16s <- merged_otu_16s %>%
  rename_with(~ ifelse(startsWith(.x, "ESSE."), .x, paste0("ESSE.", .x)))
merged_otu_shotgun <- merged_otu_shotgun %>%
  rename_with(~ ifelse(startsWith(.x, "ESSE."), .x, paste0("ESSE.", .x)))

# rename the OTU column from ESSE.OTU to OTU
merged_otu_16s <- merged_otu_16s %>% 
  rename('OTU' = 'ESSE.OTU')
merged_otu_shotgun <- merged_otu_shotgun %>% 
  rename('OTU' = 'ESSE.OTU')


#******************************************************************************
#* lots of tidying and fixing required BEFORE we merge these 2 DFs into 1
#* 

#step 1 - deal with Domain vs Kingdom
merged_otu_shotgun <- merged_otu_shotgun %>% 
  mutate(OTU = str_replace(OTU, "d__Bacteria;", "k__Bacteria;"))
merged_otu_shotgun <- merged_otu_shotgun %>% 
  mutate(OTU = str_replace(OTU, "d__Eukaryota;", "k__Eukaryota;"))

#step 2 - deal with the species name including the genus
merged_otu_shotgun <- merged_otu_shotgun %>% 
  mutate(OTU = str_replace(OTU, "s__[^ ]+ ", "s__"))

#step 3 : try to deal with the missing CLASS values from the shotgun data
#approach will be to try and find the correct class from the 16s data,
# and us it to "fill in" the right value.

# Step 1: Parse df_16s into a lookup table
lookup <- merged_otu_16s %>%
  mutate(
    family = str_extract(OTU, "f__[^;]+"),
    genus  = str_extract(OTU, "g__[^;]+"),
    species = str_extract(OTU, "s__[^;]+"),
    class   = str_extract(OTU, "c__[^;]+")
  ) %>%
  select(family, genus, species, class)

# Step 2: Define the row-wise fixer function
fix_otu <- function(otu_string, lookup_df) {
  fam <- str_extract(otu_string, "f__[^;]+")
  gen <- str_extract(otu_string, "g__[^;]+")
  spe <- str_extract(otu_string, "s__[^;]+")
  
  match_row <- lookup_df %>%
    filter(family == fam, genus == gen, species == spe)
  
  if (nrow(match_row) == 1) {
    class_val <- match_row$class[1]
    str_replace(otu_string, "(p__[^;]+;)", paste0("\\1", class_val, ";"))
  } else {
    str_replace(otu_string, "(p__[^;]+;)", "\\1c__NA;")
  }
}


# Step 3: Apply to df_shotgun using purrr::map
df_shotgun <- merged_otu_shotgun %>%
  mutate(OTU = map_chr(OTU, ~fix_otu(.x, lookup)))


#merge these two dataframes now : the df_shotgun and the merged_otu_16s 

# Step 1: Pivot both dataframes to long format
long_16s <- merged_otu_16s %>%
  pivot_longer(cols = starts_with("ESSE.222"), names_to = "sample_id", values_to = "abundance")

long_shotgun <- df_shotgun %>%
  pivot_longer(cols = starts_with("ESSE.222"), names_to = "sample_id", values_to = "abundance")


# Step 2: Combine and aggregate
combined_long <- bind_rows(long_16s, long_shotgun) %>%
  group_by(OTU, sample_id) %>%
  summarise(abundance = sum(abundance, na.rm = TRUE), .groups = "drop")

# Step 4: Pivot back to wide format
merged_df <- combined_long %>%
  pivot_wider(names_from = sample_id, values_from = abundance, values_fill = 0)



#add in a unique identifier column
merged_otu <- merged_df %>% 
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

split_taxa <- split_taxa %>% 
  mutate(OTU = str_replace(OTU, "d__Viruses;", "k__Viruses;"))

#will need to retain a link btwn a OTU identifier and the taxonomic table
# deal with thte fact that the TSV files have a Domain instead of a Kingdom


  
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

#deal with NA values at the GENUS level
for_tax_df_unique <- for_tax_df_unique %>% 
  mutate(
    Genus = case_when(
      is.na(Genus) | toupper(trimws(Genus)) == "NA" ~ paste0(Family, "(F)_NA"),
      TRUE ~ as.character(Genus)
    )
  )

#now do the same for Family
for_tax_df_unique <- for_tax_df_unique %>% 
  mutate(
    Family = case_when(
      is.na(Family) | toupper(trimws(Family)) == "NA" ~ paste0(Order, "(O)_NA"),
      TRUE ~ as.character(Family)
    )
  )



for_tax_df_unique <- for_tax_df_unique %>%
  tibble::column_to_rownames("otu_id")

write_csv(for_tax_df_unique, "00_data/03_cleaned_data/tax_table.csv")

  # for use in constructing our phyloseq obj

#_______________________________________________________________________________
# 02_ construct otu_table ----
merged_otu <- merged_otu %>% 
  select(-OTU)
otu_df <- merged_otu %>%
  tibble::column_to_rownames("otu_id")
write_csv(otu_df, "00_data/03_cleaned_data/otu_table.csv")
 # for use in 
                                                # constructing our phyloseq obj

#_______________________________________________________________________________
# 03_ construct sample_table ----
# this is the metadata table: swab number, and then other cols. 
# for now, only 'real row' is the experimental group
#below code is NOT the correct, but an outlien of what is required

# for the metadata, we have to deal with 2 seperate files
# FILE 1: metadata for 16s samples
# FILE 2: metadata for shotgun samples
# PLAN: import them seperatley, tidy them up and make sure columns are correctly uniform
# then, 

swab_details <- read_csv("00_data/02_metadata/metadata.csv")
swab_details <- swab_details %>% 
  janitor::clean_names()
#remove uneeded columns
swab_details <- swab_details %>% 
  select("corrected_id", "description", "description1", "description2", "description3", "sample_id")
#change column names to more meaningful ones
swab_details <- swab_details %>% 
  rename('swab_type' = "description") %>% 
  rename('experiment' = 'description1') %>% 
  rename('experiment_name' = 'description2') %>% 
  rename('notes' = 'description3')

swab_details <- swab_details %>% 
  mutate(seq_type ="16s")


#create a new batch name, from PART of the sample_id
swab_details <- swab_details %>% 
  separate(sample_id, into = c("batch_id", "rest"), sep = "_", remove = FALSE)
swab_details <- swab_details %>% 
  select(-c(rest))



swab_details <- swab_details %>%
  group_by(corrected_id) %>%
  filter(!(corrected_id == "ESSE.222.0032" & row_number() == 1)) %>%
  ungroup()


#_______________________________________________________________________________
#now do the same kind of thing for the shotgun samples
swab_details_shotgun <- read_csv("00_data/02_metadata/Shotgun.csv")
swab_details_shotgun <- swab_details_shotgun %>% 
   janitor::clean_names()
#remove unneccesary columns 
swab_details_shotgun <- swab_details_shotgun %>% 
  select("sample_lable", "description", "description1", "description2", "description3", "internal_id")

swab_details_shotgun <- swab_details_shotgun %>% 
  rename('swab_type' = "description") %>% 
  rename('experiment' = 'description1') %>% 
  rename('experiment_name' = 'description2') %>% 
  rename('notes' = 'description3') %>% 
  rename('sample_id' = 'internal_id')

swab_details_shotgun <- swab_details_shotgun %>% 
  mutate(seq_type ="shotgun")

#create a new batch name, from PART of the sample_id
swab_details_shotgun <- swab_details_shotgun %>% 
  separate(sample_id, into = c("batch_id", "rest"), sep = "_", remove = FALSE)
swab_details_shotgun <- swab_details_shotgun %>% 
  select(-c(rest))


swab_details_shotgun <- swab_details_shotgun %>% 
  rename('corrected_id' = 'sample_lable')

#MERGE the two into 1 merged_swab_meta
merged_metadata <- rbind(swab_details, swab_details_shotgun)

merged_metadata <- merged_metadata %>% 
  rename('otu_id' = 'corrected_id')

  
#prep this DF for sample_table
swab_df <- merged_metadata %>%
  as.data.frame() %>%              # Convert tibble to plain data frame
  column_to_rownames("otu_id")

#we have a duplicate row for ESSE.222.0032
write_csv(swab_df, "00_data/03_cleaned_data/sample_table.csv")


# now create the RDS from all of the above

otu_mat <- as.matrix(otu_df)
class(otu_mat) <- "numeric"
otu_rel <- phyloseq::otu_table(otu_mat, taxa_are_rows = TRUE)


tax_mat <- as.matrix(for_tax_df_unique)
tax_phy <- phyloseq::tax_table(tax_mat)

sample_phy <- phyloseq::sample_data(swab_df)


microbiome_phyloseq_obj <- phyloseq::phyloseq(otu_rel, tax_phy, sample_phy)

##SAVE this now
saveRDS(microbiome_phyloseq_obj, file = "00_data/03_cleaned_data/phyloseq_object.rds")
