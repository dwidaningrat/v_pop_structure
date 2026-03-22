# CONFIG #####

setwd("/Users/waodedwidaningrat/viridans/capsular_genes_flanking/")
library(dplyr)
library(readr)
library(stringr)
library(purrr)
library(tidyr)
library(dplyr)
library(dplyr)
library(purrr)
library(tidyr)
library(readr)
library(tibble)

# Read files
fileA <- read.csv("all_cps_genes_blast_rbh_14jan2026_modified.csv", stringsAsFactors = FALSE)
fileB <- read.csv("cps_gene_pathway.csv", stringsAsFactors = FALSE)
# If file B has unique query_id, this will work directly.
# If file B has duplicates, take the first pathway for each ID:
fileB_unique <- fileB %>%
  group_by(query_id) %>%
  summarise(pathway = first(pathway), .groups = "drop")

fileA$query_id <- sub("_1$", "", fileA$query_id)

# Join to add pathway into file A
merged_data <- fileA %>%
  left_join(fileB_unique, by = "query_id")

# Inspect
head(merged_data)

# Save output
write.csv(merged_data, "cps_genes_blast_allgenomes_pathway_rbh.csv", row.names = FALSE)


#PARSING BLAST RESULTS TO SEE DISTRIBUTION OF CPS GENES AND PATHWAYS PRESENCE ABSENCE #####
# Input / output files
infile  <- "cps_genes_blast_allgenomes_pathway_rbh.csv"
outfile <- "contiguity_cps-genes_with_pathways_rbh.csv"

ignore_case <- TRUE


# Read CSV
df <- read_csv(infile, show_col_types = FALSE)

df_best <- df %>%
  # 1. Filter out weak hits
  dplyr::filter(percent_identity >= 50) %>%
  
  # 2. Rank remaining hits by alignment quality
  dplyr::arrange(
    dplyr::desc(bit_score),
    dplyr::desc(percent_identity),
    dplyr::desc(align_length)
  ) %>%
  
  # 3. Keep one hit per accession × gene
  dplyr::group_by(accession, query_id) %>%
  dplyr::slice_head(n = 1) %>%
  dplyr::ungroup() %>%
  
  # 4. Optional normalization
  dplyr::mutate(
    query_up = if (ignore_case) toupper(query_id) else query_id
  )



collapse_safe <- function(vec, sep = ";") {
  vec <- vec[!is.na(vec) & vec != ""]
  if (length(vec) == 0) return(NA_character_)
  paste(unique(vec), collapse = sep)
}

# CPSA–D gene set
cpsabcd_set <- c("CPSA", "CPSB", "CPSC", "CPSD")

# Function to get contig that has ALL CPSA–D genes (single CPSABCD contig)
get_contig_target <- function(ci, genes_up) {
  present <- ci %>% filter(query_up %in% genes_up)
  if (nrow(present) == 0) return(NA_character_)
  
  contig_genes <- present %>%
    group_by(contig) %>%
    summarise(n_genes = n_distinct(query_up), .groups = "drop") %>%
    filter(n_genes == length(genes_up))  # require all genes_up
  
  if (nrow(contig_genes) == 0) return(NA_character_)
  
  contig_genes$contig[1]
}

# Function to get all genes in a contig (by name of contig)
get_all_genes_in_contig <- function(ci, contig_name) {
  if (is.na(contig_name)) return(NA_character_)
  genes <- ci %>% filter(contig == contig_name) %>% pull(query_id)
  collapse_safe(genes)
}

# Function to check if all given genes are on the same contig
check_same_contig_fun <- function(ci, genes) {
  present <- ci %>% filter(query_id %in% genes)
  if (nrow(present) <= 1) return(NA_character_)
  uniq_contigs <- unique(present$contig)
  if (length(uniq_contigs) == 1) "yes" else "no"
}

# Function to get genes in a given range around CPSA on the CPSABCD contig
get_genes_around_cpsa <- function(ci, cps_contig, cpsabcd_contiguos,
                                  upstream = 2000, downstream = 30298) {
  cps_contig <- cps_contig[1]
  cpsabcd_contiguos <- cpsabcd_contiguos[1]
  
  # Only proceed if actually have a full CPSABCD cluster on a single contig
  if (is.na(cps_contig) || cpsabcd_contiguos != "yes") return(character(0))
  
  # Restrict to that CPSABCD contig
  tbl <- ci %>% dplyr::filter(contig == cps_contig)
  
  start_cpsa <- tbl %>%
    dplyr::filter(query_up == "CPSA") %>%
    pull(start)
  if (length(start_cpsa) == 0) return(character(0))
  
  min_coord <- max(0, start_cpsa - upstream)
  max_coord <- start_cpsa + downstream
  
  tbl %>%
    dplyr::filter(start >= min_coord & start <= max_coord) %>%
    dplyr::arrange(start) %>%
    pull(query_id) %>%
    unique()
}

# summarise per genome
summary_df <- df_best %>%
  group_by(accession) %>%
  summarise(
    gene_count = n(),
    gene_name = collapse_safe(query_id),
    ci = list(tibble(
      query_id = query_id,
      query_up = query_up,
      contig   = contig,
      start    = start,
      pathway  = pathway
    )),
    .groups = "drop"
  ) %>%
  mutate(
    # single contig that contains ALL CPSABCD genes
    contig_name_cpsabcd = map_chr(ci, ~ get_contig_target(.x, cpsabcd_set)),
    
    # genome-wide cpsABCD count
    gene_count_cpsabcd_genome = map_int(
      ci,
      ~ n_distinct(.x$query_up[.x$query_up %in% cpsabcd_set])
    ),
    
    genes_contiguos_cpsabcd = map2_chr(ci, contig_name_cpsabcd, get_all_genes_in_contig),
    
    # count distinct CPSABCD genes on that contig only
    gene_count_cpsabcd = map2_int(
      ci, contig_name_cpsabcd,
      ~ if (is.na(.y)) 0L else
        .x %>%
        filter(contig == .y, query_up %in% cpsabcd_set) %>%
        summarise(n = n_distinct(query_up)) %>%
        pull(n)
    ),
    
    # CPSABCD considered "contiguous" only if all 4 are present on that single contig
    cpsabcd_contiguos = if_else(gene_count_cpsabcd == length(cpsabcd_set),
                                "yes", "no")
  )

# Pathways
pathways <- df_best$pathway %>% unique() %>% na.omit()

# Build pathway columns with CPSA neighborhood info

pathway_summary <- summary_df %>%
  mutate(
    pathway_info = pmap(
      list(ci, cpsabcd_contiguos, contig_name_cpsabcd),
      function(tbl, cps_flag, cps_contig) {
        out <- list()
        
        # thresholds per pathway
        contig_thresholds <- list(
          "initialtransferase" = 1, "rhamnosyltranferase" = 1, "galactopyranose" = 1,
          "phosphotransferase" = 1, "acetyltransferase" = 1, "pyruvyltransferase" = 1,
          "transposase" = 1, "lipoprotein" = 1, "polymerase" = 1, "flippase" = 1, "glucugalacturonic" = 2,
          "pseudogenes" = 1, "hypothetical" = 1, "glycerolcdp" =1, "glycerolndp" =3,
          "arabinitol" = 2, "mannac" = 1, "glycosyltransferase" = 2,
          "fucnac" = 3, "mannitol" = 2, "ribofuranose" = 1,
          "rhamnose" = 4
        )
        
        for (pw in pathways) {
          pw_genes <- tbl %>% filter(pathway == pw) %>% pull(query_id)
          
          # basic counts / contiguity (ignoring CPSABCD)
          out[[paste0(pw, "_gene_count")]] <- length(pw_genes)
          out[[paste0(pw, "_gene_names")]] <- collapse_safe(pw_genes)
          out[[paste0(pw, "_contiguous")]] <- check_same_contig_fun(tbl, pw_genes)
          
          #same contig as CPSABCD: ALL genes on that single contig
          if (cps_flag == "yes" && !is.na(cps_contig) && length(pw_genes) > 0) {
            pw_contigs_unique <- tbl %>%
              filter(query_id %in% pw_genes) %>%
              pull(contig) %>%
              unique()
            
            out[[paste0(pw, "_samecontig_as_cpsABCD")]] <-
              if (length(pw_contigs_unique) == 1 && pw_contigs_unique == cps_contig) "yes" else "no"
          } else {
            out[[paste0(pw, "_samecontig_as_cpsABCD")]] <- NA_character_
          }
          
          #same contig as CPSABCD (with thresholds)
          if (cps_flag == "yes" && !is.na(cps_contig) && length(pw_genes) > 0) {
            pw_contigs <- tbl %>%
              filter(query_id %in% pw_genes) %>%
              pull(contig)
            
            # count genes on *the CPSABCD contig only*
            n_same_contig <- sum(pw_contigs == cps_contig)
            threshold <- contig_thresholds[[pw]]
            if (is.null(threshold)) threshold <- 1
            
            out[[paste0(pw, "_contig_cpsABCD")]] <-
              if (n_same_contig >= threshold) "yes" else "no"
          } else {
            out[[paste0(pw, "_contig_cpsABCD")]] <- NA_character_
          }
          
          #genes within CPSA ranges (restricted to CPSABCD contig)
          for (downstream in c(2000, 3000, 4000, 5000, 10000, 15000, 20000, 25000, 30298)) {
            genes_in_range <- get_genes_around_cpsa(
              ci = tbl,
              cps_contig = cps_contig,
              cpsabcd_contiguos = cps_flag,
              upstream = 2000,
              downstream = downstream
            )
            
            # only genes in this pathway AND in the CPSABCD-contig window
            genes_in_range_pw <- intersect(pw_genes, genes_in_range)
            
            # names
            colname_base <- paste0(pw, "_within_2k_", downstream)
            out[[colname_base]] <- collapse_safe(genes_in_range_pw)
            
            # count
            out[[paste0(colname_base, "_count")]] <- length(genes_in_range_pw)
          }
          
          #same contig AND within CPSA window (2k up, 30298 down; with thresholds)
          if (cps_flag == "yes" && !is.na(cps_contig) && length(pw_genes) > 0) {
            
            # genes of this pathway on CPSABCD contig within [CPSA-2000, CPSA+30298]
            genes_window <- get_genes_around_cpsa(
              ci                 = tbl,
              cps_contig         = cps_contig,
              cpsabcd_contiguos  = cps_flag,
              upstream           = 2000,
              downstream         = 30298
            )
            
            # count of pathway’s genes in that window
            n_window <- length(intersect(pw_genes, genes_window))
            
            # use same custom thresholds as for *_contig_cpsABCD
            threshold <- contig_thresholds[[pw]]
            if (is.null(threshold)) threshold <- 1
            
            out[[paste0(pw, "_range_cpsABCD")]] <-
              if (n_window >= threshold) "yes" else "no"
            
          } else {
            out[[paste0(pw, "_range_cpsABCD")]] <- NA_character_
          }
        }
        
        out
      }
    )
  ) %>%
  unnest_wider(pathway_info)

# Remove ci column, write output
final_df <- pathway_summary %>% select(-ci)
write_csv(final_df, outfile, na = "")
message("Done! Output written to: ", outfile)



#CREATE TABLE DISTRIBUTION ACROSS SPECIES #####

library(openxlsx)
# Read files
cps_genes <- read.csv("contiguity_cps-genes_with_pathways_rbh.csv", stringsAsFactors = FALSE)
species <- read.csv("accession_species_group.csv", stringsAsFactors = FALSE)


# Join to add pathway info
merged_data <- cps_genes %>%
  left_join(species, by = "accession")

# Inspect
head(merged_data)

# Save output
write.csv(merged_data, "contiguity_cps-genes_with_pathways_species_rbh.csv", row.names = FALSE)


# Input / output files
infile  <- "contiguity_cps-genes_with_pathways_species_rbh.csv"

library(dplyr)

# Read CSV
dataset <- read_csv(infile, show_col_types = FALSE)


# Convert all blank strings to NA
dataset[dataset == ""] <- NA

# Trim whitespace from character columns
dataset <- data.frame(lapply(dataset, function(x) if(is.character(x)) trimws(x) else x))

# define species order by vgs and build lookup
if (!all(c("species", "vgs") %in% names(dataset))) {
  stop("Both 'species' and 'vgs' columns must be present in the dataset.")
}

species_vgs_lookup <- dataset %>%
  select(Species = species, vgs) %>%
  distinct()

species_order <- species_vgs_lookup %>%
  arrange(vgs, Species) %>%
  pull(Species)

# ensure species is an ordered factor
dataset$species <- factor(dataset$species, levels = species_order)

# Vector of column names to tabulate against species_coregenome
cols_to_tab <- c(
  "cpsabcd_contiguos", "gene_count_cpsabcd", "gene_count_cpsabcd_genome",
  "arabinitol_gene_count", "arabinitol_contiguous", "arabinitol_contig_cpsABCD", "arabinitol_range_cpsABCD",
  "glycosyltransferase_gene_count", "glycosyltransferase_contiguous", "glycosyltransferase_contig_cpsABCD", "glycosyltransferase_range_cpsABCD",
  "chlorohydrolase_gene_count", "chlorohydrolase_contiguous", "chlorohydrolase_contig_cpsABCD", "chlorohydrolase_range_cpsABCD",
  "lipoprotein_gene_count", "lipoprotein_contiguous", "lipoprotein_contig_cpsABCD", "lipoprotein_range_cpsABCD",
  "acetyltransferase_gene_count", "acetyltransferase_contiguous", "acetyltransferase_contig_cpsABCD", "acetyltransferase_range_cpsABCD",
  "fucnac_gene_count", "fucnac_contiguous", "fucnac_contig_cpsABCD", "fucnac_range_cpsABCD",
  "glucugalacturonic_gene_count", "glucugalacturonic_contiguous", "glucugalacturonic_contig_cpsABCD", "glucugalacturonic_range_cpsABCD",
  "galactopyranose_gene_count", "galactopyranose_contiguous", "galactopyranose_contig_cpsABCD", "galactopyranose_range_cpsABCD",
  "glycerolcdp_gene_count", "glycerolcdp_contiguous", "glycerolcdp_contig_cpsABCD", "glycerolcdp_range_cpsABCD",
  "glycerolndp_gene_count", "glycerolndp_contiguous", "glycerolndp_contig_cpsABCD", "glycerolndp_range_cpsABCD",
  "mannitol_gene_count", "mannitol_contiguous", "mannitol_contig_cpsABCD", "mannitol_range_cpsABCD",
  "mannac_gene_count", "mannac_contiguous", "mannac_contig_cpsABCD", "mannac_range_cpsABCD",
  "rhamnose_gene_count", "rhamnose_contiguous", "rhamnose_contig_cpsABCD", "rhamnose_range_cpsABCD",
  "phosphotransferase_gene_count", "phosphotransferase_contiguous", "phosphotransferase_contig_cpsABCD", "phosphotransferase_range_cpsABCD",
  "initialtransferase_gene_count", "initialtransferase_contiguous", "initialtransferase_contig_cpsABCD", "initialtransferase_range_cpsABCD",
  "ribofuranose_gene_count", "ribofuranose_contiguous", "ribofuranose_contig_cpsABCD", "ribofuranose_range_cpsABCD",
  "rhamnosyltranferase_gene_count", "rhamnosyltranferase_contiguous", "rhamnosyltranferase_contig_cpsABCD", "rhamnosyltranferase_range_cpsABCD",
  "hypothetical_gene_count", "hypothetical_contiguous", "hypothetical_contig_cpsABCD", "hypothetical_range_cpsABCD",
  "transposase_gene_count", "transposase_contiguous", "transposase_contig_cpsABCD", "transposase_range_cpsABCD",
  "flippase_gene_count", "flippase_contiguous", "flippase_contig_cpsABCD", "flippase_range_cpsABCD",
  "pyruvyltransferase_gene_count", "pyruvyltransferase_contiguous", "pyruvyltransferase_contig_cpsABCD", "pyruvyltransferase_range_cpsABCD",
  "polymerase_gene_count", "polymerase_contiguous", "polymerase_contig_cpsABCD", "polymerase_range_cpsABCD"
)

# Columns that are present in the dataset
present_cols <- cols_to_tab[cols_to_tab %in% names(dataset)]
present_cols

# Columns that are missing
missing_cols <- cols_to_tab[!cols_to_tab %in% names(dataset)]
missing_cols

# Print summary
cat("Number of columns to tabulate:", length(cols_to_tab), "\n")
cat("Present in dataset:", length(present_cols), "\n")
cat("Missing from dataset:", length(missing_cols), "\n")



# Create all contingency tables #####
tables_list <- lapply(cols_to_tab, function(col) {
  table(dataset$species, dataset[[col]])
})


# Name the list elements for easy access
names(tables_list) <- cols_to_tab

# Example: print first contingency table
tables_list[[1]]


# Create a new workbook
wb <- createWorkbook()

make_unique_sheet_name <- function(base_name, wb, max_len = 31) {
  # existing sheet names (case-insensitive)
  existing <- tolower(names(wb))
  
  # initial truncated candidate
  base_trunc <- substr(base_name, 1, max_len)
  candidate <- base_trunc
  
  if (!(tolower(candidate) %in% existing)) {
    return(candidate)
  }
  
  # if exists, append _1, _2, ... while respecting max_len
  i <- 1
  repeat {
    suffix <- paste0("_", i)
    cut_len <- max_len - nchar(suffix)
    candidate <- paste0(substr(base_trunc, 1, cut_len), suffix)
    
    if (!(tolower(candidate) %in% existing)) {
      return(candidate)
    }
    
    i <- i + 1
  }
}


# Loop over the tables and add each as a separate sheet
for (col_name in names(tables_list)) {
  tbl <- tables_list[[col_name]]
  
  # Skip empty tables
  if (nrow(tbl) == 0 || ncol(tbl) == 0) {
    warning(paste("Skipping empty table for column:", col_name))
    next
  }
  
  sheet_name <- make_unique_sheet_name(col_name, wb)
  
  df <- as.data.frame.matrix(tbl)
  
  # order rows by species_order
  df <- df[species_order[species_order %in% rownames(df)], , drop = FALSE]
  
  # add Species and vgs
  species_vec <- rownames(df)
  vgs_vec <- species_vgs_lookup$vgs[match(species_vec, species_vgs_lookup$Species)]
  
  df <- cbind(
    Species = species_vec,
    vgs     = vgs_vec,
    df
  )
  rownames(df) <- NULL
  
  addWorksheet(wb, sheetName = sheet_name)
  writeData(wb, sheet = sheet_name, x = df)
}


# Save the workbook
saveWorkbook(wb, file = "distribution_cps-pathway_species_rbh.xlsx", overwrite = TRUE)


# Create a new workbook
wb <- createWorkbook()

# Loop over the tables and add each as a separate sheet
for (col_name in names(tables_list)) {
  tbl <- tables_list[[col_name]]
  
  if (nrow(tbl) == 0 || ncol(tbl) == 0) {
    warning(paste("Skipping empty table for column:", col_name))
    next
  }
  
  sheet_name <- make_unique_sheet_name(col_name, wb)
 
  df <- as.data.frame.matrix(tbl)
  
  # order rows by vgs-sorted species_order
  df <- df[species_order[species_order %in% rownames(df)], , drop = FALSE]
  
  df_percent <- round(df / rowSums(df) * 100, 3)
  
  species_vec <- rownames(df)
  vgs_vec <- species_vgs_lookup$vgs[match(species_vec, species_vgs_lookup$Species)]
  
  df_combined <- data.frame(
    Species = species_vec,
    vgs     = vgs_vec
  )
  
  for (j in 1:ncol(df)) {
    count_col   <- df[, j]
    percent_col <- df_percent[, j]
    
    df_combined[[paste0(colnames(df)[j], "_n")]]       <- count_col
    df_combined[[paste0(colnames(df)[j], "_percent")]] <- percent_col
  }
  
  rownames(df_combined) <- NULL
  
  addWorksheet(wb, sheetName = sheet_name)
  writeData(wb, sheet = sheet_name, x = df_combined)
}


# Save the workbook
saveWorkbook(wb, file = "distribution_cps-pathway_species_with_rowpercent_rbh.xlsx", overwrite = TRUE)


# Identify columns ending with "contig_cpsabcd" that exist in tables_list
cpsabcd_cols <- grep("contig_cpsABCD$", names(tables_list), value = TRUE)

# Initialize a list to store processed tables
df_list <- list()

for (col_name in cpsabcd_cols) {
  tbl <- tables_list[[col_name]]
  
  # Skip empty tables
  if (nrow(tbl) == 0 || ncol(tbl) == 0) next
  
  # Convert to data frame
  df <- as.data.frame.matrix(tbl)
  
  # Calculate row-wise percentages
  df_percent <- round(df / rowSums(df) * 100, 3)
  
  # Combine count and percentage as separate columns
  df_combined <- data.frame(Species = rownames(df))
    #, vgs     = species_vgs_lookup$vgs[match(rownames(df), species_vgs_lookup$Species)])
  
  for (j in 1:ncol(df)) {
    count_col <- df[, j]
    percent_col <- df_percent[, j]
    
    df_combined[[paste0(col_name, "_", colnames(df)[j], "_n")]] <- count_col
    df_combined[[paste0(col_name, "_", colnames(df)[j], "_percent")]] <- percent_col
  }
  
  df_list[[col_name]] <- df_combined
}

# Only proceed if df_list is not empty
if (length(df_list) == 0) {
  stop("No valid tables found for columns ending with 'contig_cpsabcd'. Check your dataset and tables_list.")
}

# Merge all data frames by Species
big_table <- df_list[[1]]
if (length(df_list) > 1) {
  for (i in 2:length(df_list)) {
    big_table <- full_join(big_table, df_list[[i]], by = "Species")
  }
}

# Replace NA with 0
big_table[is.na(big_table)] <- 0

# add vgs column and order by vgs, then Species
big_table <- big_table %>%
  mutate(
    vgs = species_vgs_lookup$vgs[match(Species, species_vgs_lookup$Species)]
  ) %>%
  relocate(vgs, .after = Species) %>%
  arrange(vgs, Species)


# order Species by vgs-based species_order
big_table$Species <- factor(big_table$Species, levels = species_order)
big_table <- big_table[order(big_table$Species), ]
big_table$Species <- as.character(big_table$Species)


# Save as Excel
wb <- createWorkbook()
addWorksheet(wb, "cpsABCD_combined")
writeData(wb, "cpsABCD_combined", big_table)
saveWorkbook(wb, "cpsABCD_combined_table_rbh.xlsx", overwrite = TRUE)

# View first rows
head(big_table)


# Identify columns ending with "range_cpsABCD" that exist in tables_list
cpsabcd2_cols <- grep("range_cpsABCD$", names(tables_list), value = TRUE)

# Initialize a list to store processed tables
df_list <- list()

for (col_name in cpsabcd2_cols) {
  tbl <- tables_list[[col_name]]
  
  # Skip empty tables
  if (nrow(tbl) == 0 || ncol(tbl) == 0) next
  
  # Convert to data frame
  df <- as.data.frame.matrix(tbl)
  
  # Calculate row-wise percentages
  df_percent <- round(df / rowSums(df) * 100, 3)
  
  # Combine count and percentage as separate columns
  df_combined <- data.frame(Species = rownames(df))
  #, vgs     = species_vgs_lookup$vgs[match(rownames(df), species_vgs_lookup$Species)])
  
  for (j in 1:ncol(df)) {
    count_col <- df[, j]
    percent_col <- df_percent[, j]
    
    df_combined[[paste0(col_name, "_", colnames(df)[j], "_n")]] <- count_col
    df_combined[[paste0(col_name, "_", colnames(df)[j], "_percent")]] <- percent_col
  }
  
  df_list[[col_name]] <- df_combined
}

# Only proceed if df_list is not empty
if (length(df_list) == 0) {
  stop("No valid tables found for columns ending with 'range_cpsabcd'. Check your dataset and tables_list.")
}

# Merge all data frames by Species
big_table2 <- df_list[[1]]
if (length(df_list) > 1) {
  for (i in 2:length(df_list)) {
    big_table2 <- full_join(big_table2, df_list[[i]], by = "Species")
  }
}

# Replace NA with 0
big_table2[is.na(big_table2)] <- 0

# add vgs column and order by vgs, then Species
big_table2 <- big_table2 %>%
  mutate(
    vgs = species_vgs_lookup$vgs[match(Species, species_vgs_lookup$Species)]
  ) %>%
  relocate(vgs, .after = Species) %>%
  arrange(vgs, Species)


# order Species by vgs-based species_order
big_table2$Species <- factor(big_table2$Species, levels = species_order)
big_table2 <- big_table2[order(big_table2$Species), ]
big_table2$Species <- as.character(big_table2$Species)


# Save as Excel
wb <- createWorkbook()
addWorksheet(wb, "cpsABCD_combined_inrange")
writeData(wb, "cpsABCD_combined_inrange", big_table2)
saveWorkbook(wb, "cpsABCD_combined_inrange_table_rbh.xlsx", overwrite = TRUE)

# View first rows
head(big_table2)



# Identify all tables in tables_list that are gene_count tables
gene_count_cols <- grep("_gene_count$", names(tables_list), value = TRUE, ignore.case = TRUE)

# Initialize a list to store data frames
df_list <- list()

for (col_name in gene_count_cols) {
  tbl <- tables_list[[col_name]]
  
  # Skip empty tables
  if (nrow(tbl) == 0 || ncol(tbl) == 0) next
  
  # Convert table to data frame (rows = species, columns = categories)
  df <- as.data.frame.matrix(tbl)
  
  # Calculate row-wise percentages
  df_percent <- round(df / rowSums(df) * 100, 3)
  
  # Create new data frame with species as first column
  df_combined <- data.frame(Species = rownames(df))
  #,vgs     = species_vgs_lookup$vgs[match(rownames(df), species_vgs_lookup$Species)])
  
  # Add count and percent columns for each category
  for (j in 1:ncol(df)) {
    df_combined[[paste0(col_name, "_", colnames(df)[j], "_n")]] <- df[, j]
    df_combined[[paste0(col_name, "_", colnames(df)[j], "_percent")]] <- df_percent[, j]
  }
  
  df_list[[col_name]] <- df_combined
}

# Only proceed if df_list is not empty
if (length(df_list) == 0) {
  stop("No valid gene_count tables found in tables_list.")
}

# Merge all gene_count tables by Species
big_gene_count_table <- df_list[[1]]
if (length(df_list) > 1) {
  for (i in 2:length(df_list)) {
    big_gene_count_table <- full_join(big_gene_count_table, df_list[[i]], by = "Species")
  }
}

# Replace NA with 0
big_gene_count_table[is.na(big_gene_count_table)] <- 0

big_gene_count_table <- big_gene_count_table %>%
  mutate(
    vgs = species_vgs_lookup$vgs[match(Species, species_vgs_lookup$Species)]
  ) %>%
  relocate(vgs, .after = Species) %>%
  arrange(vgs, Species)


# order Species by vgs-based species_order
big_gene_count_table$Species <- factor(big_gene_count_table$Species, levels = species_order)
big_gene_count_table <- big_gene_count_table[order(big_gene_count_table$Species), ]
big_gene_count_table$Species <- as.character(big_gene_count_table$Species)


# Save as Excel
wb <- createWorkbook()
addWorksheet(wb, "gene_count_with_percent")
writeData(wb, "gene_count_with_percent", big_gene_count_table)
saveWorkbook(wb, "gene_count_with_percent_table_rbh.xlsx", overwrite = TRUE)

# View first rows
head(big_gene_count_table)



#PLOTTING #####
library(dplyr)
library(tidyr)
library(ggplot2)
library(readr)   # read_csv
library(grid)    # for unit()
library(stringr) # for str_wrap if needed


# CONFIG
species_csv <- "species_order_cps_plot.csv"
pathway_csv <- "pathway_order_plot.csv"   # CSV that defines x-axis category order and decorated names

# Exclude pattern
exclude_pattern <- "flippase|polymerase|lipoprotein|hypothetical|chlorohydrolase"

apply_exclude <- TRUE

show_species_counts <- FALSE

show_group_counts <- FALSE

group_label_wrap_width <- 18

# Read species + group order from CSV
if (!file.exists(species_csv)) stop("CSV file not found: ", species_csv)
sp_df <- read_csv(species_csv, show_col_types = FALSE)

# assume first column = Species, second column = Group, take rows 2:32
if (nrow(sp_df) < 32) stop("species_order_cps_plot.csv appears to have fewer than 32 rows.")
sp_sub <- sp_df[2:32, ]
species_col <- names(sp_sub)[1]
group_col   <- names(sp_sub)[2]

species_order <- as.character(sp_sub[[species_col]])
group_lookup <- sp_sub %>%
  select(!!sym(species_col), !!sym(group_col)) %>%
  rename(Species = !!sym(species_col), Group = !!sym(group_col)) %>%
  mutate(Species = as.character(Species), Group = as.character(Group))

# Read pathway order CSV (defines x axis order + decorated names)
if (!file.exists(pathway_csv)) stop("pathway_order_plot.csv not found: ", pathway_csv)
path_df <- read_csv(pathway_csv, show_col_types = FALSE)
if (ncol(path_df) < 1) stop("pathway_order_plot.csv must have at least one column with pathway short names")

# assume first column is the short pathway id,
# and if present, second column is the decorated pathway_name to display on the axis
path_short <- as.character(path_df[[1]])
if (ncol(path_df) >= 2) {
  path_name <- as.character(path_df[[2]])
} else {
  # fallback: if no decorated name column, use the short name as the display name
  path_name <- path_short
}

# build a lookup table for short -> decorated name
path_map <- tibble(source_short = path_short, pathway_name = path_name)

# Prepare long counts from big_table2
n_cols <- grep("_n$", names(big_table2), value = TRUE)
if (length(n_cols) == 0) stop("No '*_n' columns found in big_table2")

counts_long <- big_table2 %>%
  select(Species, vgs, all_of(n_cols)) %>%
  pivot_longer(
    cols = -c(Species, vgs),
    names_to = c("source", "level"),
    names_pattern = "^(.*)_(.*)_n$",
    values_to = "n"
  )

counts_long <- counts_long %>%
  group_by(Species, vgs, source) %>%
  mutate(N = sum(n, na.rm = TRUE)) %>%
  ungroup()

yes_long <- counts_long %>%
  filter(tolower(level) == "yes") %>%
  mutate(prop = ifelse(N > 0, n / N, NA_real_),
         percent = round(prop * 100, 1)) %>%
  # shorten x-axis display label
  mutate(source_short = sub("_.*$", "", source))

# Apply exclude pattern
if (apply_exclude) {
  yes_long <- yes_long %>%
    filter(!grepl(exclude_pattern, source_short, ignore.case = TRUE))
} else {
  message("apply_exclude = FALSE: not filtering categories by exclude_pattern")
}

# Compute species-level totals and group totals
# species_N as max(N) observed for that species
species_totals <- yes_long %>%
  group_by(Species) %>%
  summarize(species_N = max(N, na.rm = TRUE), .groups = "drop") %>%
  left_join(group_lookup, by = "Species")  # add Group column from CSV

# group_N as sum of species_N for species present
group_totals <- species_totals %>%
  group_by(Group) %>%
  summarize(group_N = sum(species_N, na.rm = TRUE), .groups = "drop")

species_totals <- species_totals %>%
  left_join(group_totals, by = "Group")

# Build plotting labels
species_totals <- species_totals %>%
  mutate(
    Species_label = if (show_species_counts) {
      sprintf("%s (%d)", Species, ifelse(is.na(species_N), 0L, species_N))
    } else {
      as.character(Species)
    }
  )

# Build Group_label
species_totals <- species_totals %>%
  mutate(
    Group_clean = ifelse(is.na(Group), "Unknown", Group),
    Group_label_plain = if (show_group_counts) {
      paste0(Group_clean, " (", ifelse(is.na(group_N), 0L, group_N), ")")
    } else {
      Group_clean
    },
    # wrap the plain label here as a single string with \n inserted by str_wrap
    Group_label_wrapped = str_wrap(Group_label_plain, width = group_label_wrap_width)
  )

# join labels back to long data
yes_long <- yes_long %>%
  left_join(species_totals %>% select(Species, Species_label, Group, Group_label_wrapped), by = "Species") %>%
  rename(Group_label = Group_label_wrapped)

# Apply species order from CSV and create ordered Species_label factor
# Keep only species present both in CSV and data
present_species <- intersect(species_order, unique(yes_long$Species))
if (length(present_species) == 0) stop("None of the species from CSV appear in the data.")

# Build ordered vector of Species_label according to CSV order
label_map <- species_totals %>%
  filter(Species %in% present_species) %>%
  mutate(Species = factor(Species, levels = species_order)) %>%
  arrange(Species)

species_label_levels <- label_map$Species_label

# convert Species_label to factor with levels in desired order
yes_long <- yes_long %>%
  filter(Species %in% present_species) %>%
  mutate(
    Species_label = factor(Species_label, levels = species_label_levels),
    Group_label = factor(Group_label, levels = unique(label_map$Group_label_wrapped[order(match(label_map$Species, species_order))]))
  )

# using pathway_order_plot.csv mapping and use decorated pathway_name for labels
# distinct source -> short mapping
source_map <- yes_long %>%
  distinct(source, source_short)

# attach pathway position  and decorated name by matching source_short to path_map
source_map <- source_map %>%
  left_join(path_map, by = "source_short") %>%
  mutate(pathway_pos = match(source_short, path_short))

# For sources not in the CSV, give them positions after the listed pathways (preserve alphabetical within)
max_pos <- ifelse(all(is.na(source_map$pathway_pos)), 0, max(source_map$pathway_pos, na.rm = TRUE))
source_map <- source_map %>%
  mutate(pathway_pos = ifelse(is.na(pathway_pos), max_pos + row_number(), pathway_pos)) %>%
  arrange(pathway_pos, source)

# final factor levels for 'source' left -> right
source_levels <- source_map %>% pull(source)

# mapping from source -> decorated label; fallback to short if no decorated name
labels_map <- source_map %>%
  mutate(pathway_name = ifelse(is.na(pathway_name), source_short, pathway_name)) %>%
  select(source, pathway_name) %>%
  deframe()  # named character vector: names = source, values = pathway_name

# apply ordering
yes_long <- yes_long %>% mutate(source = factor(source, levels = source_levels))

# prepare tile text (blank for NA)
yes_long_ori <- yes_long %>%
  mutate(plot_label = ifelse(is.na(prop), "", sprintf("%.1f%%\n(%d/%d)", percent, n, N)),
         text_col = ifelse(is.na(prop), "black", ifelse(prop > 0.5, "white", "black")))

yes_long <- yes_long %>%
  mutate(plot_label = case_when(
    is.na(prop) ~ "",
    n == 0 & N > 0 ~ sprintf("0.0%%\n(%d/%d)", n, N),
    prop > 0 & (prop * 100) < 0.05 ~ sprintf("<0.1%%\n(%d/%d)", n, N),
    TRUE ~ sprintf("%.1f%%\n(%d/%d)", round(prop * 100, 1), n, N)
    ),
    text_col = ifelse(is.na(prop), "black", ifelse(prop > 0.5, "white", "black"))
    )

# helper to create parseable expressions for y labels so species is italic and "Streptococcus" → "S. "
make_y_expr <- function(lbl) {
  # lbl is like "Species" or "Species (123)"
  lbl <- as.character(lbl)
  m <- regmatches(lbl, regexec("^(.+?) \\((\\d+)\\)$", lbl))[[1]]
  if (length(m) == 0) {
    sp_full <- lbl
    # abbreviate Streptococcus -> S. (only when at the start)
    sp_abbrev <- sub("^Streptococcus\\s+", "S. ", sp_full)
    sp_esc <- gsub('"', '\\"', sp_abbrev, fixed = TRUE)
    return(paste0('italic("', sp_esc, '")'))
  } else {
    sp_full <- m[2]
    num <- m[3]
    sp_abbrev <- sub("^Streptococcus\\s+", "S. ", sp_full)
    sp_esc <- gsub('"', '\\"', sp_abbrev, fixed = TRUE)
    # italic species name, plain count in parentheses
    return(paste0('italic("', sp_esc, '")*plain(" (', num, ')")'))
  }
}




# Build the 3 plots (Layer1, Layer2, Layer3) #####

library(dplyr); library(tidyr); library(ggplot2); library(stringr); library(grid)

# Denominators
# 1) total genomes per species (from accession-level 'dataset')
if (!exists("dataset")) stop("`dataset` (accession-level table) not found")
# allow both 'species' and 'Species'
species_colname <- if ("species" %in% names(dataset)) "species" else if ("Species" %in% names(dataset)) "Species" else stop("No species column found in dataset")
species_totals_tbl <- dataset %>%
  count(!!sym(species_colname), name = "N_total") %>%
  rename(Species = !!sym(species_colname))

N_total_vec <- setNames(as.integer(species_totals_tbl$N_total), species_totals_tbl$Species)

# 2) genomes per species that have contiguous cpsABCD
if (!"cpsabcd_contiguos" %in% names(dataset)) stop("`cpsabcd_contiguos` column missing in dataset")
species_cpsabcd_tbl <- dataset %>%
  filter(!is.na(cpsabcd_contiguos) & tolower(cpsabcd_contiguos) == "yes") %>%
  count(!!sym(species_colname), name = "N_cpsabcd_yes") %>%
  rename(Species = !!sym(species_colname))

# ensure species with 0 are represented
species_all <- species_totals_tbl$Species
species_cpsabcd_full <- tibble(Species = species_all) %>%
  left_join(species_cpsabcd_tbl, by = "Species") %>%
  mutate(N_cpsabcd_yes = ifelse(is.na(N_cpsabcd_yes), 0L, as.integer(N_cpsabcd_yes)))

N_cpsabcd_vec <- setNames(as.integer(species_cpsabcd_full$N_cpsabcd_yes), species_cpsabcd_full$Species)

# Helper to create the standardized long table used for plotting
# input_table: data.frame (Species + many *_n columns)
# pattern: regex to match *_n columns to pivot
# filter_level: if provided, keep only rows where level == filter_level (e.g., "yes")
# special_for_presence: if TRUE, treat any category != "0" as presence and sum them (for big_gene_count_table)
build_yes_long <- function(input_table, pattern = "_n$", filter_level = NULL, special_for_presence = FALSE) {
  if (!("Species" %in% names(input_table))) {
    # try first column as Species
    input_table <- input_table %>% rename(Species = 1)
  }
  # find relevant *_n columns
  ncols <- grep(pattern, names(input_table), value = TRUE, ignore.case = TRUE)
  if (length(ncols) == 0) return(NULL)
  
  long <- input_table %>%
    select(Species, all_of(ncols)) %>%
    pivot_longer(cols = -Species, names_to = c("source", "level"), names_pattern = "^(.*)_(.*)_n$", values_to = "n") %>%
    mutate(source_short = sub("_.*$", "", source)) %>%
    # normalize level to lowercase
    mutate(level = tolower(level))
  
  if (special_for_presence) {
    pres_long <- long %>%
      mutate(pres_flag = ifelse(level == "0", 0L, 1L)) %>%
      group_by(Species, source, source_short) %>%
      summarise(n = sum(n[pres_flag == 1], na.rm = TRUE), .groups = "drop")
    pres_long <- pres_long %>% mutate(level = "presence")
    return(pres_long)
  } else {
    if (!is.null(filter_level)) {
      long <- long %>% filter(level == filter_level)
    }
    return(long)
  }
}

# Layer 1: any gene present in genome (use big_gene_count_table)
if (!exists("big_gene_count_table")) stop("big_gene_count_table not found; cannot make Layer1")
layer1_long <- build_yes_long(big_gene_count_table, pattern = "_n$", filter_level = NULL, special_for_presence = TRUE)
# attach denominators N_total and compute prop
layer1_long <- layer1_long %>%
  mutate(N = as.integer(N_total_vec[Species])) %>%
  mutate(prop = ifelse(!is.na(N) & N > 0, n / N, NA_real_), percent = round(prop * 100, 1))


# Layer 2: pathway on same contig as cpsABCD (use big_table, level == "yes")
if (!exists("big_table")) stop("big_table not found; cannot make Layer2")
layer2_long <- build_yes_long(big_table, pattern = "_n$", filter_level = "yes", special_for_presence = FALSE)
# denominators = number genomes with contiguous cpsABCD per species
layer2_long <- layer2_long %>%
  mutate(N = as.integer(N_cpsabcd_vec[Species])) %>%
  mutate(prop = ifelse(!is.na(N) & N > 0, n / N, NA_real_), percent = round(prop * 100, 1))

# Layer 3: pathway on same contig AND within range (use big_table2, level == "yes")
if (!exists("big_table2")) stop("big_table2 not found; cannot make Layer3")
layer3_long <- build_yes_long(big_table2, pattern = "_n$", filter_level = "yes", special_for_presence = FALSE)
layer3_long <- layer3_long %>%
  mutate(N = as.integer(N_cpsabcd_vec[Species])) %>%
  mutate(prop = ifelse(!is.na(N) & N > 0, n / N, NA_real_), percent = round(prop * 100, 1))

# suppress (blank) Layer3 cells where Layer2 had n == 0
suppress_layer3_when_layer2_zero <- TRUE

if (suppress_layer3_when_layer2_zero) {
  
  # ensure layer3 has plot_label and text_col
  layer3_long <- layer3_long %>%
    mutate(
      plot_label = case_when(
        is.na(prop) ~ "",
        n == 0 & N > 0 ~ sprintf("0.0%%\n(%d/%d)", n, N),
        prop > 0 & (prop * 100) < 0.05 ~ sprintf("<0.1%%\n(%d/%d)", n, N),
        TRUE ~ sprintf("%.1f%%\n(%d/%d)", round(prop * 100, 1), n, N)
      ),
      text_col = ifelse(is.na(prop), "black", ifelse(prop > 0.5, "white", "black")),
      Species = as.character(Species),
      source_short = as.character(source_short)
    )
  
  # prepare layer2 summary keyed by Species + source_short
  layer2_j <- layer2_long %>%
    mutate(
      Species = as.character(Species),
      source_short = as.character(source_short),
      n_layer2 = as.integer(n),
      N_layer2 = as.integer(N)
    ) %>%
    select(Species, source_short, n_layer2, N_layer2) %>%
    distinct()
  
  # diagnostics: unique keys
  n_keys_l2 <- nrow(distinct(layer2_j, Species, source_short))
  n_keys_l3 <- nrow(distinct(layer3_long, Species, source_short))
  message("Unique (Species, source_short) keys: layer2=", n_keys_l2, " ; layer3=", n_keys_l3)
  
  # left-join layer2 info into layer3
  layer3_joined <- layer3_long %>%
    left_join(layer2_j, by = c("Species", "source_short"))
  
  # report rows without matching layer2 entry
  n_na_join <- sum(is.na(layer3_joined$n_layer2))
  if (n_na_join > 0) {
    message("Note: ", n_na_join, " layer3 rows had no matching layer2 key (n_layer2 is NA). ",
            "This can happen if 'source_short' names differ between big_table and big_table2.")
    print(head(distinct(layer3_joined %>% filter(is.na(n_layer2)) %>% select(Species, source_short)), 8))
  }
  
  # treat NA n_layer2 as 0 for suppression decision
  layer3_supp <- layer3_joined %>%
    mutate(
      n_layer2 = ifelse(is.na(n_layer2), 0L, as.integer(n_layer2)),
      N_layer2 = ifelse(is.na(N_layer2), 0L, as.integer(N_layer2))
    )
  
  n_to_suppress <- sum(layer3_supp$n_layer2 == 0, na.rm = TRUE)
  total_rows <- nrow(layer3_supp)
  message("Layer3 suppression active: ", n_to_suppress, " / ", total_rows, " rows will be blanked because layer2 n == 0")
  
  # now blank (set to NA or empty) only where layer2 n == 0
  layer3_final <- layer3_supp %>%
    mutate(
      prop = ifelse(n_layer2 == 0, NA_real_, prop),
      percent = ifelse(n_layer2 == 0, NA_real_, percent),
      plot_label = ifelse(n_layer2 == 0, "", plot_label),
      text_col = ifelse(is.na(prop), "black", ifelse(prop > 0.5, "white", "black"))
    ) %>%
    select(-n_layer2, -N_layer2)
  
  # overwrite layer3_long with suppressed version
  layer3_long <- layer3_final
  message("Finished: layer3_long updated with suppression (if any).")
  
} else {
  message("Layer3 suppression disabled; no changes made to layer3_long.")
}

# Function to make plot from a long df (follows your template exactly) 
# Titles / subtitles / captions for each layer
plot_title_1   <- "LAYER 1: Proportion of genomes within each species with pathway gene detected anywhere in the genome"
plot_subtitle_1 <- ""
plot_caption_1 <- paste0(
  "Description:\n",
  "n = number of genomes within each species having at least 1 gene of the corresponding pathway, could be located anywhere in the genome.\n",
  "N = total number of genomes for each species in the dataset.\n"
)

plot_title_2   <- "LAYER 2: Proportion of genomes within each species with pathway gene detected in the same contig of cpsABCD"
plot_subtitle_2 <- ""
plot_caption_2 <- paste0(
  "Description:\n",
  "n = number of genomes within each species detected with the corresponding pathway with threshold applied for each pathway\n",
  "(e.g. at least 4 rhamnose genes need to be in the same contig as cpsABCD to be considered 'present' for rhamnose pathway).\n",
  "N = number of genomes for each species with contiguous cpsABCD.\n"
)

plot_title_3   <- "LAYER 3: Proportion of genomes within each species with pathway gene detected in the same contig of cpsABCD and within cps-like range"
plot_subtitle_3 <- "specified cps-like range (cps-like regions): within range of 2,000bp upstream and 30,298bp downstream of cpsA where cpsABCD is contiguous"
plot_caption_3 <- paste0(
  "Description:\n",
  "n = number of genomes within each species detected with the corresponding pathway with threshold applied for each pathway\n",
  "(e.g. at least 4 rhamnose genes need to be in the same contig as cpsABCD to be considered 'present' for rhamnose pathway),\n",
  "and those pathway gene(s) must be within the specified cps-like range.\n",
  "N = number of genomes for each species with contiguous cpsABCD.\n"
)

#  Replace function: make_plot_from_yes_long now accepts title, subtitle, caption
make_plot_from_yes_long <- function(yes_long, title, subtitle = "", caption = "") {
  if (is.null(yes_long) || nrow(yes_long) == 0) {
    message("No data for ", title); return(NULL)
  }
  
  # optional exclude
  if (exists("apply_exclude") && apply_exclude) {
    yes_long <- yes_long %>% filter(!grepl(exclude_pattern, source_short, ignore.case = TRUE))
  }
  
  # species / group totals for labels
  species_tot_local <- yes_long %>%
    group_by(Species) %>%
    summarise(species_N = max(N, na.rm = TRUE), .groups = "drop") %>%
    left_join(group_lookup, by = "Species")
  group_tot_local <- species_tot_local %>% group_by(Group) %>% summarise(group_N = sum(species_N, na.rm = TRUE), .groups = "drop")
  species_tot_local <- species_tot_local %>% left_join(group_tot_local, by = "Group")
  
  species_tot_local <- species_tot_local %>%
    mutate(
      Species_label = if (exists("show_species_counts") && show_species_counts) sprintf("%s (%d)", Species, ifelse(is.na(species_N), 0L, species_N)) else as.character(Species),
      Group_clean = ifelse(is.na(Group), "Unknown", Group),
      Group_label_plain = if (exists("show_group_counts") && show_group_counts) paste0(Group_clean, " (", ifelse(is.na(group_N), 0L, group_N), ")") else Group_clean,
      Group_label_wrapped = str_wrap(Group_label_plain, width = ifelse(exists("group_label_wrap_width"), group_label_wrap_width, 18))
    )
  
  yes_long <- yes_long %>% left_join(species_tot_local %>% select(Species, Species_label, Group, Group_label_wrapped), by = "Species") %>% rename(Group_label = Group_label_wrapped)
  
  # ordering species
  present_species <- intersect(species_order, unique(yes_long$Species))
  if (length(present_species) == 0) stop("None of species from CSV appear in data for plot ", title)
  label_map <- species_tot_local %>% filter(Species %in% present_species) %>% mutate(Species = factor(Species, levels = species_order)) %>% arrange(Species)
  species_label_levels <- label_map$Species_label
  yes_long <- yes_long %>% filter(Species %in% present_species) %>%
    mutate(Species_label = factor(Species_label, levels = species_label_levels),
           Group_label = factor(Group_label, levels = unique(label_map$Group_label_wrapped[order(match(label_map$Species, species_order))])))
  
  # order x axis using path_map
  source_map_local <- yes_long %>% distinct(source, source_short) %>% left_join(path_map, by = "source_short")
  source_map_local <- source_map_local %>% mutate(pathway_pos = match(source_short, path_short))
  max_pos_local <- ifelse(all(is.na(source_map_local$pathway_pos)), 0, max(source_map_local$pathway_pos, na.rm = TRUE))
  source_map_local <- source_map_local %>% mutate(pathway_pos = ifelse(is.na(pathway_pos), max_pos_local + row_number(), pathway_pos)) %>% arrange(pathway_pos, source)
  source_levels_local <- source_map_local$source
  labels_map_local <- source_map_local %>% mutate(pathway_name = ifelse(is.na(pathway_name), source_short, pathway_name)) %>% select(source, pathway_name) %>% deframe()
  yes_long <- yes_long %>% mutate(source = factor(source, levels = source_levels_local))
  
  # tile labels
  yes_long <- yes_long %>%
    mutate(plot_label = case_when(
      is.na(prop) ~ "",
      n == 0 & N > 0 ~ sprintf("0.0%%\n(%d/%d)", n, N),
      prop > 0 & (prop * 100) < 0.05 ~ sprintf("<0.1%%\n(%d/%d)", n, N),
      TRUE ~ sprintf("%.1f%%\n(%d/%d)", round(prop * 100, 1), n, N)
    ),
    text_col = ifelse(is.na(prop), "black", ifelse(prop > 0.5, "white", "black")))
  
  # y label parser uses make_y_expr from your script
  p <- ggplot(yes_long, aes(x = source, y = Species_label, fill = prop)) +
    geom_tile(color = "grey40", linewidth = 0.1, height = 0.9, width = 0.9) +
    geom_text(aes(label = plot_label, colour = text_col), size = 2.1, lineheight = 0.9) +
    scale_x_discrete(labels = function(x) labels_map_local[x]) +
    scale_y_discrete(labels = function(x) parse(text = vapply(x, make_y_expr, FUN.VALUE = character(1)))) +
    facet_grid(rows = vars(Group_label), scales = "free_y", space = "free", switch = "y",
               labeller = labeller(Group_label = label_wrap_gen(width = ifelse(exists("group_label_wrap_width"), group_label_wrap_width, 18)))) +
    scale_fill_gradientn(name = "Proportion\n(0–1)", colours = c("#fff5f0", "#fb6a4a", "#67000d"), limits = c(0,1), na.value = "grey95") +
    scale_colour_identity(guide = "none") +
    labs(title = title, subtitle = subtitle, caption = caption, x = "Sugar Biosynthetic Pathways, Transferases, Transposase", y = NULL) +
    theme_minimal(base_size = 12) +
    theme(
      strip.placement = "outside",
      strip.text.y.left = element_text(angle = 0, hjust = 1, size = 11, face = "bold"),
      strip.background = element_rect(fill = "grey95", colour = NA),
      panel.grid = element_blank(),
      panel.spacing = unit(0.15, "lines"),
      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
      axis.text.y = element_text(size = 9),
      legend.position = "right",
      plot.caption.position = "plot",
      plot.caption = element_text(size = 10, face = "plain", hjust = 0, lineheight = 0.95),
      plot.subtitle = element_text(size = 10, face = "plain", margin = margin(b = 6))
    )
  
  p
}

# Build & print the three plots with unique labels
p_layer1 <- make_plot_from_yes_long(layer1_long, title = plot_title_1, subtitle = plot_subtitle_1, caption = plot_caption_1)
p_layer2 <- make_plot_from_yes_long(layer2_long, title = plot_title_2, subtitle = plot_subtitle_2, caption = plot_caption_2)
p_layer3 <- make_plot_from_yes_long(layer3_long, title = plot_title_3, subtitle = plot_subtitle_3, caption = plot_caption_3)

if (!is.null(p_layer1)) print(p_layer1)
if (!is.null(p_layer2)) print(p_layer2)
if (!is.null(p_layer3)) print(p_layer3)



#PLOT 4, BASICALLY THIS IS LAYER 3 PLOT, BUT THE N IS NUMBER OF CONTIGUOUS CPSABCD GENOMES THE THE PATHWAY GENE, NOT ALL GENOMES WITH CONTIGUOUS CPSABCD
#THE n is number of genomes with the corresponding pathway present in the same contig AND WITHIN RANGE

plot_subtitle <- "specified cps-like range (cps-like regions): within range of 2,000bp upstream and 30,298bp downstream of cpsA where cpsABCD is contiguous"

# Plot title/subtitle/caption
plot_title   <- "LAYER 4: Proportion of genomes within each species with pathway gene detected in the cps-like regions"
plot_subtitle <- "specified cps-like range (cps-like regions): within range of 2,000bp upstream and 30,298bp downstream of cpsA where cpsABCD is contiguous"

# explicit caption
plot_caption <- paste0(
  "Three layers:\n",
  "1) pathway gene(s) detected in the genome ; 2) pathway gene(s) located in the same contig as cpsABCD; 3) pathway gene(s) within the specified range of cpsABCD\n",
  "This plot shows proportions only for genomes with contiguous cpsABCD (2 and 3).\n",
  "Within each cell: N = number of genomes having the corresponding pathway AND located on the same contig as cpsABCD; ",
  "n = number of those genomes with the pathway located within the specified cps-like range.\n",
  "Empty cell mean 0 genome in that species having the corresponding pathway gene (s) anywhere in the genome or if any they not in the same contig as cpsABCD\n",
  "or cpsABCD not contiguos in that genome (see 1st and 2nd layer plots for details).\n",
  "Coloured cells show the proportion of genomes (n/N) with pathway gene(s) located on the same cpsABCD contig and within the specified range."
)

# PLOT using facets by Group_label on the left
p_yes_prop_groups <- ggplot(yes_long, aes(x = source, y = Species_label, fill = prop)) +
  geom_tile(color = "grey40", linewidth = 0.1, height = 0.9, width = 0.9) +
  geom_text(aes(label = plot_label, colour = text_col), size = 2.1, lineheight = 0.9) +
  # use decorated pathway_name labels (fall back to short if no decorated name)
  scale_x_discrete(labels = function(x) labels_map[x]) +
  # make y labels italic for the species part, and abbreviate Streptococcus -> S.
  scale_y_discrete(labels = function(x) parse(text = vapply(x, make_y_expr, FUN.VALUE = character(1)))) +
  # use label_wrap_gen to wrap Group_label; Group_label is plain wrapped text already but label_wrap_gen ensures further wrapping
  facet_grid(rows = vars(Group_label), scales = "free_y", space = "free", switch = "y",
             labeller = labeller(Group_label = label_wrap_gen(width = group_label_wrap_width))) +
  scale_fill_gradientn(
    name = "Proportion\n(0–1)",
    colours = c("#fff5f0", "#fb6a4a", "#67000d"),
    limits = c(0, 1),
    na.value = "grey95"
  ) +
  scale_colour_identity(guide = "none") +
  labs(
    title = plot_title,
    subtitle = plot_subtitle,
    caption = plot_caption,
    x = "Sugar Biosynthetic Pathways, Transferases, Transposase",
    y = NULL
  ) +
  theme_minimal(base_size = 12) +
  theme(
    # place the strip (group label) on the left; make it bold
    strip.placement = "outside",
    strip.text.y.left = element_text(angle = 0, hjust = 1, size = 11, face = "bold"),
    strip.background = element_rect(fill = "grey95", colour = NA),
    panel.grid = element_blank(),
    panel.spacing = unit(0.15, "lines"),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    axis.text.y = element_text(size = 9),
    legend.position = "right",
    plot.caption.position = "plot",          # position caption relative to plot area
    plot.caption          = element_text(size = 10, face = "plain", hjust = 0, lineheight = 0.95),
    # subtitle styling
    plot.subtitle         = element_text(size = 10, face = "plain", margin = margin(b = 6))
  )

p_yes_prop_groups










# Plot distribution of pathway_count (number of pathways in-range per genome) #####

library(dplyr)
library(tidyr)
library(ggplot2)
library(readr)
library(stringr)
library(grid)

species_csv <- "species_order_cps_plot.csv"
accession_file <- "contiguity_cps-genes_with_pathways_species_rbh.csv" # produced earlier

# Controls
show_species_counts <- FALSE   # include (Species (N)) on y-axis labels
show_group_counts   <- FALSE
group_label_wrap_width <- 18

# Read accession-level data
if (!exists("dataset")) {
  if (!file.exists(accession_file)) stop("Accession-level file not found: ", accession_file)
  dataset <- read_csv(accession_file, show_col_types = FALSE)
}

# normalize column names
names(dataset) <- make.names(names(dataset), unique = TRUE)

# find the *_range_cpsABCD columns
range_cols <- grep("range_cpsABCD$", names(dataset), value = TRUE, ignore.case = TRUE)
if (length(range_cols) == 0) stop("No '*_range_cpsABCD' columns found in dataset")

# compute pathway_count per accession: count of columns == "yes"
dataset2 <- dataset %>%
  mutate_at(vars(all_of(range_cols)), ~ ifelse(is.na(.), NA_character_, tolower(as.character(.)))) %>%
  rowwise() %>%
  mutate(
    pathway_count = sum(c_across(all_of(range_cols)) == "yes", na.rm = TRUE)
  ) %>%
  ungroup()

# sanity: pathway_count should be integer >= 0
dataset2$pathway_count <- as.integer(dataset2$pathway_count)

# compute denominators (total genomes per species) 
species_colname <- if ("species" %in% names(dataset2)) "species" else if ("Species" %in% names(dataset2)) "Species" else stop("No species column found in dataset")

species_totals_tbl <- dataset2 %>%
  count(!!sym(species_colname), name = "N_total") %>%
  rename(Species = !!sym(species_colname))

# ensure Species column exists in dataset2 as character for joins
dataset2 <- dataset2 %>% mutate(Species = as.character(!!sym(species_colname)))

#  tabulate pathway_count by species
pathcount_by_species <- dataset2 %>%
  group_by(Species, pathway_count) %>%
  summarise(n = n(), .groups = "drop") %>%
  right_join(species_totals_tbl, by = "Species") %>%   # ensure species with zero counts for some pathway_count preserved
  mutate(n = ifelse(is.na(n), 0L, n),
         N = N_total,
         prop = ifelse(N > 0, n / N, NA_real_),
         percent = round(prop * 100, 1))

# make sure pathway_count factor levels from 0..max (so x-axis consistent)
max_pc <- max(pathcount_by_species$pathway_count, na.rm = TRUE)
pc_levels <- as.character(0:max_pc)
pathcount_by_species <- pathcount_by_species %>%
  mutate(pathway_count = factor(as.character(pathway_count), levels = pc_levels))

#  prepare species order & group lookup like in the other plots
if (!file.exists(species_csv)) stop("CSV file not found: ", species_csv)
sp_df <- read_csv(species_csv, show_col_types = FALSE)
if (nrow(sp_df) < 2) stop("species_order_cps_plot.csv seems too small")

# assume first two columns: Species, Group
species_col <- names(sp_df)[1]
group_col   <- names(sp_df)[2]

sp_sub <- sp_df[2:32, ]
species_order <- as.character(sp_sub[[species_col]])
group_lookup <- sp_sub %>%
  select(!!sym(species_col), !!sym(group_col)) %>%
  rename(Species = !!sym(species_col), Group = !!sym(group_col)) %>%
  mutate(Species = as.character(Species), Group = as.character(Group))

# add group information to pathcount table
pathcount_by_species <- pathcount_by_species %>%
  left_join(group_lookup, by = "Species")

# compute species & group labels for plotting
species_totals_plot <- pathcount_by_species %>%
  group_by(Species, Group) %>%
  summarise(species_N = max(N, na.rm = TRUE), .groups = "drop") %>%
  group_by(Group) %>%
  mutate(group_N = sum(species_N, na.rm = TRUE)) %>%
  ungroup() %>%
  mutate(
    Species_label = if (show_species_counts) sprintf("%s (%d)", Species, species_N) else as.character(Species),
    Group_label_plain = if (show_group_counts) paste0(ifelse(is.na(Group), "Unknown", Group), " (", ifelse(is.na(group_N), 0L, group_N), ")") else ifelse(is.na(Group), "Unknown", Group),
    Group_label_wrapped = str_wrap(Group_label_plain, width = group_label_wrap_width)
  )

pathcount_by_species <- pathcount_by_species %>%
  left_join(species_totals_plot %>% select(Species, Species_label, Group, Group_label_wrapped), by = "Species") %>%
  rename(Group_label = Group_label_wrapped)

# filter to species in CSV order
present_species <- intersect(species_order, unique(pathcount_by_species$Species))
if (length(present_species) == 0) stop("None of the species from CSV appear in the data.")

label_map <- species_totals_plot %>%
  filter(Species %in% present_species) %>%
  mutate(Species = factor(Species, levels = species_order)) %>%
  arrange(Species)

species_label_levels <- label_map$Species_label

pathcount_by_species <- pathcount_by_species %>%
  filter(Species %in% present_species) %>%
  mutate(
    Species_label = factor(Species_label, levels = species_label_levels),
    Group_label = factor(Group_label, levels = unique(label_map$Group_label_wrapped[order(match(label_map$Species, species_order))]))
  )

# helper for italic species and S. abbreviation
make_y_expr <- function(lbl) {
  lbl <- as.character(lbl)
  m <- regmatches(lbl, regexec("^(.+?) \\((\\d+)\\)$", lbl))[[1]]
  if (length(m) == 0) {
    sp_full <- lbl
    sp_abbrev <- sub("^Streptococcus\\s+", "S. ", sp_full)
    sp_esc <- gsub('"', '\\"', sp_abbrev, fixed = TRUE)
    return(paste0('italic("', sp_esc, '")'))
  } else {
    sp_full <- m[2]
    num <- m[3]
    sp_abbrev <- sub("^Streptococcus\\s+", "S. ", sp_full)
    sp_esc <- gsub('"', '\\"', sp_abbrev, fixed = TRUE)
    return(paste0('italic("', sp_esc, '")*plain(" (', num, ')")'))
  }
}

# Create tile label text and text color
pathcount_by_species <- pathcount_by_species %>%
  mutate(plot_label = ifelse(is.na(prop), "", sprintf("%.1f%%\n(%d/%d)", percent, n, N)),
         text_col = ifelse(is.na(prop), "black", ifelse(prop > 0.5, "white", "black")))

# Plot title/subtitle/caption
plot_title_pc <- "Distribution of number of pathways located in cpsABCD-range per genome"
plot_subtitle_pc <- "pathway_count = number of distinct pathways with threshold met and within cpsABCD range (per genome)"
plot_caption_pc <- paste0(
  "Each cell: percent (n/N)\n",
  "n = number of genomes in the species with pathway_count = X (i.e. X pathways met the 'in-range' criteria)\n",
  "N = total number of genomes for that species in the dataset\n",
  "Pathway_count is computed by counting pathways for which *_range_cpsABCD == 'yes' for that accession."
)

p_pathcount <- ggplot(pathcount_by_species, aes(x = pathway_count, y = Species_label, fill = prop)) +
  geom_tile(color = "grey40", linewidth = 0.15, height = 0.9, width = 0.9) +
  geom_text(aes(label = plot_label, colour = text_col), size = 2.1, lineheight = 0.9) +
  scale_x_discrete(expand = c(0,0)) +
  scale_y_discrete(labels = function(x) parse(text = vapply(x, make_y_expr, FUN.VALUE = character(1)))) +
  facet_grid(rows = vars(Group_label), scales = "free_y", space = "free", switch = "y",
             labeller = labeller(Group_label = label_wrap_gen(width = group_label_wrap_width))) +
  scale_fill_gradientn(name = "Proportion\n(0–1)", colours = c("#fff5f0", "#fb6a4a", "#67000d"), limits = c(0,1), na.value = "grey95") +
  scale_colour_identity(guide = "none") +
  labs(title = plot_title_pc, subtitle = plot_subtitle_pc, caption = plot_caption_pc, x = "pathway_count (number of pathways in-range per genome)", y = NULL) +
  theme_minimal(base_size = 12) +
  theme(
    strip.placement = "outside",
    strip.text.y.left = element_text(angle = 0, hjust = 1, size = 11, face = "bold"),
    strip.background = element_rect(fill = "grey95", colour = NA),
    panel.grid = element_blank(),
    panel.spacing = unit(0.15, "lines"),
    axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0.5),
    axis.text.y = element_text(size = 9),
    legend.position = "right",
    plot.caption.position = "plot",
    plot.caption = element_text(size = 9, hjust = 0, lineheight = 0.95),
    plot.subtitle = element_text(size = 10, margin = margin(b = 6))
  )

# print the plot
print(p_pathcount)


# Save species x pathway_count tables to Excel #####
library(dplyr)
library(tidyr)
library(openxlsx)

# check required objects
if (!exists("dataset2")) stop("dataset2 not found (accession-level table with Species).")
if (!exists("pathcount_by_species")) stop("pathcount_by_species not found (species × pathway_count summary).")

# ensure Species column name consistent
if ("species" %in% names(dataset2) && !("Species" %in% names(dataset2))) {
  dataset2 <- dataset2 %>% rename(Species = species)
}
if ("species" %in% names(pathcount_by_species) && !("Species" %in% names(pathcount_by_species))) {
  pathcount_by_species <- pathcount_by_species %>% rename(Species = species)
}

# Ensure Species exists now
if (!"Species" %in% names(dataset2)) stop("dataset2 must contain a 'Species' column.")
if (!"Species" %in% names(pathcount_by_species)) stop("pathcount_by_species must contain a 'Species' column.")

# make sure pathway_count is present
if (!"pathway_count" %in% names(pathcount_by_species)) stop("pathcount_by_species must contain 'pathway_count' column.")

# Compute species-level totals (N_total) from dataset2
species_Ntbl <- dataset2 %>%
  group_by(Species) %>%
  summarise(N_total = n(), .groups = "drop")

# Optionally capture vgs if present
vgs_lookup <- NULL
if ("vgs" %in% names(dataset2)) {
  vgs_lookup <- dataset2 %>% distinct(Species, vgs)
}

# Normalize pathcount_by_species: ensure expected columns exist
# expected minimal columns: Species, pathway_count, n, N, prop, percent
required_cols <- c("Species", "pathway_count", "n", "N", "prop", "percent")
for (cl in required_cols) {
  if (!(cl %in% names(pathcount_by_species))) {
    pathcount_by_species[[cl]] <- NA
  }
}

# keep any optional metadata columns if present
meta_cols <- intersect(c("Group", "vgs", "Species_label", "Group_label"), names(pathcount_by_species))

# coerce pathway_count numeric for sorting
pathcount_by_species <- pathcount_by_species %>%
  mutate(pathway_count = suppressWarnings(as.integer(as.character(pathway_count))))

# Prepare long_out
long_out <- pathcount_by_species %>%
  select(any_of(c("Species", meta_cols, "pathway_count", "n", "N", "prop", "percent"))) %>%
  arrange(Species, pathway_count)

# Build wide counts (Species x pathway_count)
counts_wide <- long_out %>%
  select(Species, pathway_count, n) %>%
  mutate(pathway_count = ifelse(is.na(pathway_count), "NA", as.character(pathway_count))) %>%
  pivot_wider(names_from = pathway_count, values_from = n, values_fill = 0) %>%
  # attach N_total and vgs
  left_join(species_Ntbl, by = "Species")

if (!is.null(vgs_lookup)) {
  counts_wide <- counts_wide %>% left_join(vgs_lookup, by = "Species")
}

# convert counts to percent of N_total
# identify pathway_count columns
exclude_cols <- c("Species", "N_total", intersect(names(counts_wide), c("vgs")))
pathway_cols <- setdiff(names(counts_wide), exclude_cols)

percent_wide <- counts_wide
for (col in pathway_cols) {
  pct_col <- paste0(col, "_pct")
  percent_wide[[pct_col]] <- ifelse(!is.na(percent_wide$N_total) & percent_wide$N_total > 0,
                                    round(100 * percent_wide[[col]] / percent_wide$N_total, 3),
                                    NA_real_)
}

# Select percent_wide_out with Species, vgs, N_total and percent columns
percent_cols <- grep("_pct$", names(percent_wide), value = TRUE)
percent_wide_out <- percent_wide %>%
  select(Species, any_of("vgs"), N_total, all_of(percent_cols))

# Write to Excel workbook
out_xlsx <- "species_pathway_count_distribution_rbh.xlsx"
wb <- createWorkbook()

addWorksheet(wb, "counts")
writeData(wb, "counts", counts_wide)
setColWidths(wb, "counts", cols = 1:ncol(counts_wide), widths = "auto")

addWorksheet(wb, "percent")
writeData(wb, "percent", percent_wide_out)
setColWidths(wb, "percent", cols = 1:ncol(percent_wide_out), widths = "auto")

addWorksheet(wb, "long")
writeData(wb, "long", long_out)
setColWidths(wb, "long", cols = 1:ncol(long_out), widths = "auto")

# README
addWorksheet(wb, "README")
readme_text <- c(
  "species_pathway_count_distribution_rbh.xlsx",
  "",
  "Sheets:",
  " - counts : wide counts table (Species x pathway_count) plus N_total and optional vgs",
  " - percent: percentage columns for each pathway_count (columns end with _pct)",
  " - long   : long format table with Species, pathway_count, n, N, prop, percent and optional metadata",
  "",
  "Notes:",
  " - pathway_count = integer count of pathways meeting the specified condition per accession.",
  " - N_total is computed from dataset2 (accession-level table)."
)
writeData(wb, "README", as.data.frame(readme_text), colNames = FALSE)
setColWidths(wb, "README", cols = 1, widths = 120)

saveWorkbook(wb, file = out_xlsx, overwrite = TRUE)
message("Saved species × pathway_count distribution workbook to: ", out_xlsx)

write.xlsx(dataset2, file = "contiguity_cps-genes_with_pathways_species_count_rbh.xlsx")


##### DISTRIBUTION OF CPSABCD PERCENTAGE IDENTITY ACROSS DIFFERENT SPECIES #####
library(dplyr)
library(readr)
library(tidyr)
library(openxlsx)

blast_file   <- "all_cps_genes_blast_rbh_14jan2026_modified.csv"
species_file <- "accession_species_group.csv"
outfile      <- "cpsABCD_percentage_identity_by_species_rbh.xlsx"

blast <- read_csv(blast_file, show_col_types = FALSE)
species <- read_csv(species_file, show_col_types = FALSE)

blast$query_id <- sub("_1$", "", blast$query_id)

genes_of_interest <- c("CPSA", "CPSB", "CPSC", "CPSD")


blast <- blast %>%
  filter(query_id %in% genes_of_interest)

blast_best <- blast %>%
  arrange(desc(bit_score)) %>%
  group_by(accession, query_id) %>%
  slice_head(n = 1) %>%
  ungroup()

st <- blast <- read_csv("/Users/waodedwidaningrat/viridans/fastq_ena_check/final_metadata_dphil_22dec2025.csv",
                        show_col_types = FALSE)

cols_to_keep <- c("accession_ori", "species_coregenome", "vgs_coregenome", "insilico_serotype", "serobav2_serotype", "pfaster_serotype", "pfaster_probability",
                  "pfaster_note", "serobav2", "serobav2_genetic_variant", "gpsc", "cluster", "serotype_phenotypic_fix")

st <- st[, cols_to_keep]

colnames(st)[colnames(st) == "accession_ori"] <- "accession"

blast_best_st <- st %>%
  left_join(blast_best, by= "accession")

blast_best_st <- blast_best_st %>%
  mutate(
    pid_category = ifelse(percent_identity >= 50, "≥50%", "<50%")
  )

# Save output
write.csv(blast_best_st, "cpsabcd_blast_allgenomes_species_rbh.csv", row.names = FALSE)

blast_best <- blast_best %>%
  left_join(species, by = "accession")

# safety check
if (!"species" %in% names(blast_best)) {
  stop("Column 'species' not found after join. Check accession_species_group.csv")
}


blast_best <- blast_best %>%
  mutate(
    pid_category = ifelse(percent_identity >= 50, "≥50%", "<50%")
  )

species_totals <- species %>%
  count(species, name = "N_species")

# COUNT ≥50% HITS PER GENE × SPECIES
hits_50 <- blast_best %>%
  filter(pid_category == "≥50%") %>%
  count(query_id, species, name = "n_ge50")

# BUILD FULL GENE × SPECIES GRID

full_grid <- expand_grid(
  gene    = genes_of_interest,
  species = species_totals$species
)

# MERGE COUNTS + COMPUTE <50% BY COMPLEMENT

summary_table <- full_grid %>%
  left_join(hits_50, by = c("gene" = "query_id", "species")) %>%
  left_join(species_totals, by = "species") %>%
  mutate(
    n_ge50 = ifelse(is.na(n_ge50), 0L, n_ge50),
    n_lt50 = N_species - n_ge50
  ) %>%
  pivot_longer(
    cols = c(n_ge50, n_lt50),
    names_to = "pid_category",
    values_to = "n"
  ) %>%
  mutate(
    pid_category = recode(pid_category,
                          n_ge50 = "≥50%",
                          n_lt50 = "<50%"),
    proportion = round(n / N_species, 3)
  ) %>%
  arrange(gene, species, pid_category)


# WRITE EXCEL

wb <- createWorkbook()

# Combined sheet
addWorksheet(wb, "ALL_GENES")
writeData(wb, "ALL_GENES", summary_table)
setColWidths(wb, "ALL_GENES", cols = 1:ncol(summary_table), widths = "auto")

# One sheet per gene
for (g in genes_of_interest) {
  
  gene_df <- summary_table %>%
    filter(gene == g) %>%
    select(-gene)
  
  addWorksheet(wb, substr(g, 1, 31))
  writeData(wb, substr(g, 1, 31), gene_df)
  setColWidths(wb, substr(g, 1, 31), cols = 1:ncol(gene_df), widths = "auto")
}

saveWorkbook(wb, outfile, overwrite = TRUE)

message("Done. Output written to: ", outfile)












##### PLOTTING HEATMAPS #####

library(dplyr)
library(tidyr)
library(readr)
library(ComplexHeatmap)
library(circlize)
library(grid)

infile  <- "cps_genes_blast_allgenomes_pathway_rbh.csv"

# Read CSV
df <- read_csv(infile, show_col_types = FALSE)

ignore_case <- TRUE

df_best <- df %>%
  # 1. Filter out weak hits
  dplyr::filter(percent_identity >= 50) %>%
  
  # 2. Rank remaining hits by alignment quality
  dplyr::arrange(
    dplyr::desc(bit_score),
    dplyr::desc(percent_identity),
    dplyr::desc(align_length)
  ) %>%
  
  # 3. Keep one hit per accession × gene
  dplyr::group_by(accession, query_id) %>%
  dplyr::slice_head(n = 1) %>%
  dplyr::ungroup() %>%
  
  # 4. Optional normalization
  dplyr::mutate(
    query_up = if (ignore_case) toupper(query_id) else query_id
  )

pid_long <- df_best %>%
  select(
    gene = query_up,      # or query_id
    accession,
    percent_identity
  )

pid_matrix_df <- pid_long %>%
  pivot_wider(
    names_from  = accession,
    values_from = percent_identity
  )

pid_matrix <- pid_matrix_df %>%
  column_to_rownames("gene") %>%
  as.matrix()

dim(pid_matrix)      # genes × accessions
pid_matrix[1:5, 1:5] # inspect


#anywhere in the genome
pid_matrix_anywhere <- pid_matrix

#same contig as cpsABCD
pid_samecontig <- df_best %>%
  left_join(
    summary_df %>% select(accession, contig_name_cpsabcd, cpsabcd_contiguos),
    by = "accession"
  ) %>%
  filter(
    cpsabcd_contiguos == "yes",
    contig == contig_name_cpsabcd
  ) %>%
  select(
    gene = query_up,
    accession,
    percent_identity
  ) %>%
  pivot_wider(
    names_from = accession,
    values_from = percent_identity
  ) %>%
  column_to_rownames("gene") %>%
  as.matrix()

#within CPS-like region
pid_inrange <- df_best %>%
  left_join(
    summary_df %>% select(accession, contig_name_cpsabcd, cpsabcd_contiguos, ci),
    by = "accession"
  ) %>%
  filter(cpsabcd_contiguos == "yes") %>%
  rowwise() %>%
  mutate(
    in_range = query_id %in%
      get_genes_around_cpsa(
        ci = ci,
        cps_contig = contig_name_cpsabcd,
        cpsabcd_contiguos = cpsabcd_contiguos,
        upstream = 2000,
        downstream = 30298
      )
  ) %>%
  ungroup() %>%
  filter(in_range) %>%
  select(
    gene = query_up,
    accession,
    percent_identity
  ) %>%
  pivot_wider(
    names_from = accession,
    values_from = percent_identity
  ) %>%
  column_to_rownames("gene") %>%
  as.matrix()

#mapping file
cps_annot <- read_csv("cps_gene_pathway.csv")
acc_annot <- read_csv("accession_species_group_heatmap.csv")
species_hex <- read_csv("species_hex_codes.csv")

gene_order_csv <- cps_annot %>%
  pull(query_id) %>%
  unique()

enforce_gene_order <- function(mat, gene_order_csv) {
  ordered_genes <- intersect(gene_order_csv, colnames(mat))
  mat[, ordered_genes, drop = FALSE]
}

pid_matrix_anywhere_plot <- t(pid_matrix_anywhere)
pid_samecontig_plot <- t(pid_samecontig)
pid_inrange_plot <- t(pid_inrange)

pid_matrix_anywhere_plot <- enforce_gene_order(pid_matrix_anywhere_plot, gene_order_csv)
pid_samecontig_plot      <- enforce_gene_order(pid_samecontig_plot, gene_order_csv)
pid_inrange_plot         <- enforce_gene_order(pid_inrange_plot, gene_order_csv)

#Column labels
gene_labels_map <- cps_annot %>%
  distinct(query_id, proper_gene) %>%
  deframe()

make_col_labels <- function(mat, gene_labels_map) {
  cn <- colnames(mat)
  labels <- cn
  hit <- cn %in% names(gene_labels_map)
  labels[hit] <- gene_labels_map[cn[hit]]
  labels
}

col_labels_anywhere   <- make_col_labels(pid_matrix_anywhere_plot, gene_labels_map)
col_labels_samecontig <- make_col_labels(pid_samecontig_plot, gene_labels_map)
col_labels_inrange    <- make_col_labels(pid_inrange_plot, gene_labels_map)

#accession ordering + wrapped group names:

prepare_rows <- function(mat, acc_annot, wrap_width = 10) {
  
  acc_order <- acc_annot %>%
    left_join(species_hex, by = "species") %>%
    distinct(accession, species, group, hex) %>%
    arrange(group, species) %>%
    filter(accession %in% rownames(mat))
  
  mat <- mat[acc_order$accession, , drop = FALSE]
  
  row_split <- factor(
    str_wrap(acc_order$group, width = wrap_width),
    levels = unique(str_wrap(acc_order$group, width = wrap_width))
  )
  
  row_ha <- rowAnnotation(
    Species = acc_order$species,
    col = list(Species = setNames(acc_order$hex, acc_order$species)),
    show_annotation_name = FALSE
  )
  
  list(mat = mat, row_split = row_split, row_ha = row_ha)
}

anywhere_rows   <- prepare_rows(pid_matrix_anywhere_plot, acc_annot)
samecontig_rows <- prepare_rows(pid_samecontig_plot, acc_annot)
inrange_rows    <- prepare_rows(pid_inrange_plot, acc_annot)

#colour schemes
col_fun <- colorRamp2(
  c(50, 70, 90, 100),
  c("#f7fbff", "#6baed6", "#2171b5", "#08306b")
)

col_fun <- colorRamp2(
  c(0, 10, 20, 30, 40, 50, 100),
  c("white", "lightyellow", "yellow", "orange", "orangered", "red", "darkred"))

col_fun <- colorRamp2(
  c(50, 60, 70, 80, 90, 100),
  c("lightyellow", "yellow", "orange", "orangered", "red", "darkred"))


column_grid_fun <- function(j, i, x, y, width, height, fill) {
  grid.lines(
    x = unit(c(x - width/2, x - width/2), "npc"),
    y = unit(c(y - height/2, y + height/2), "npc"),
    gp = gpar(col = "grey85", lwd = 0.1)
  )
}

column_grid_fun <- function(j, i, x, y, width, height, fill) {
  if (i == 1) {
    grid.lines(
      x = unit(c(x - width/2, x - width/2), "npc"),
      y = unit(c(0, 1), "npc"),
      gp = gpar(col = "grey85", lwd = 0.1)
    )
  }
}

column_grid_fun <- function(j, i, x, y, width, height, fill) {
  grid.lines(
    x = unit(c(x - width/2, x - width/2), "npc"),
    y = unit(c(0, 1), "npc"),
    gp = gpar(col = "grey85", lwd = 0.1)
  )
}

#anywhere in the genome
Heatmap(
  anywhere_rows$mat,
  name = "Percentage Identity",
  col = col_fun,
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  row_split = anywhere_rows$row_split,
  show_row_names = FALSE,
  show_column_names = TRUE,
  column_labels = col_labels_anywhere,
  column_names_gp = gpar(fontsize = 6),
  left_annotation = anywhere_rows$row_ha,
  row_title_gp = gpar(fontsize = 5),
  na_col = "white",
  heatmap_legend_param = list(
    title = "Percentage identity (%)",
    title_gp = gpar(fontsize = 6, fontface = "bold"),
    labels_gp = gpar(fontsize = 8)
  ),
  column_title = "Percentage identity of capsular genes anywhere in the genome",
  column_title_gp = gpar(fontsize = 16, fontface = "bold"),
  cell_fun = column_grid_fun
)


#same contig
Heatmap(
  samecontig_rows$mat,
  name = "Percentage Identity",
  col = col_fun,
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  row_split = samecontig_rows$row_split,
  show_row_names = FALSE,
  show_column_names = TRUE,
  column_labels = col_labels_samecontig,
  column_names_gp = gpar(fontsize = 6),
  left_annotation = samecontig_rows$row_ha,
  row_title_gp = gpar(fontsize = 5),
  na_col = "white",
  heatmap_legend_param = list(
    title = "Percentage identity (%)",
    title_gp = gpar(fontsize = 6, fontface = "bold"),
    labels_gp = gpar(fontsize = 8)
  ),
  column_title = "Percentage identity of capsular genes located in the same contig as cpsABCD",
  column_title_gp = gpar(fontsize = 16, fontface = "bold"),
  cell_fun = column_grid_fun
)



Heatmap(
  inrange_rows$mat,
  name = "Percentage Identity",
  col = col_fun,
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  row_split = inrange_rows$row_split,
  show_row_names = FALSE,
  show_column_names = TRUE,
  column_labels = col_labels_inrange,
  column_names_gp = gpar(fontsize = 6),
  left_annotation = inrange_rows$row_ha,
  row_title_gp = gpar(fontsize = 5),
  na_col = "white",
  heatmap_legend_param = list(
    title = "Percentage identity (%)",
    title_gp = gpar(fontsize = 6, fontface = "bold"),
    labels_gp = gpar(fontsize = 8)
  ),
  column_title = "Percentage identity of capsular genes located within CPS-like regions",
  column_title_gp = gpar(fontsize = 16, fontface = "bold"),
  cell_fun = column_grid_fun
)






#####PLOTTING HEATMAPS anywhere in the genome filter by bit_score only #####

library(dplyr)
library(tidyr)
library(readr)
library(ComplexHeatmap)
library(circlize)
library(grid)

infile  <- "cps_genes_blast_allgenomes_pathway_rbh.csv"

# Read CSV
df <- read_csv(infile, show_col_types = FALSE)

ignore_case <- TRUE

#Best hit per accession × gene
df_best <- df %>%
  arrange(desc(bit_score)) %>%
  group_by(accession, query_id) %>%
  slice_head(n = 1) %>%
  ungroup() %>%
  mutate(query_up = if(ignore_case) toupper(query_id) else query_id)


pid_long <- df_best %>%
  select(
    gene = query_up,      # or query_id
    accession,
    percent_identity
  )

pid_matrix_df <- pid_long %>%
  pivot_wider(
    names_from  = accession,
    values_from = percent_identity
  )

pid_matrix <- pid_matrix_df %>%
  column_to_rownames("gene") %>%
  as.matrix()

dim(pid_matrix)      # genes × accessions
pid_matrix[1:5, 1:5] # inspect


#anywhere in the genome
pid_matrix_anywhere <- pid_matrix

#same contig as cpsABCD
pid_samecontig <- df_best %>%
  left_join(
    summary_df %>% select(accession, contig_name_cpsabcd, cpsabcd_contiguos),
    by = "accession"
  ) %>%
  filter(
    cpsabcd_contiguos == "yes",
    contig == contig_name_cpsabcd
  ) %>%
  select(
    gene = query_up,
    accession,
    percent_identity
  ) %>%
  pivot_wider(
    names_from = accession,
    values_from = percent_identity
  ) %>%
  column_to_rownames("gene") %>%
  as.matrix()

#within CPS-like region
pid_inrange <- df_best %>%
  left_join(
    summary_df %>% select(accession, contig_name_cpsabcd, cpsabcd_contiguos, ci),
    by = "accession"
  ) %>%
  filter(cpsabcd_contiguos == "yes") %>%
  rowwise() %>%
  mutate(
    in_range = query_id %in%
      get_genes_around_cpsa(
        ci = ci,
        cps_contig = contig_name_cpsabcd,
        cpsabcd_contiguos = cpsabcd_contiguos,
        upstream = 2000,
        downstream = 30298
      )
  ) %>%
  ungroup() %>%
  filter(in_range) %>%
  select(
    gene = query_up,
    accession,
    percent_identity
  ) %>%
  pivot_wider(
    names_from = accession,
    values_from = percent_identity
  ) %>%
  column_to_rownames("gene") %>%
  as.matrix()

#mapping file
cps_annot <- read_csv("cps_gene_pathway.csv")
acc_annot <- read_csv("accession_species_group_heatmap.csv")
species_hex <- read_csv("species_hex_codes.csv")

gene_order_csv <- cps_annot %>%
  pull(query_id) %>%
  unique()

enforce_gene_order <- function(mat, gene_order_csv) {
  ordered_genes <- intersect(gene_order_csv, colnames(mat))
  mat[, ordered_genes, drop = FALSE]
}

pid_matrix_anywhere_plot <- t(pid_matrix_anywhere)
pid_samecontig_plot <- t(pid_samecontig)
pid_inrange_plot <- t(pid_inrange)

pid_matrix_anywhere_plot <- enforce_gene_order(pid_matrix_anywhere_plot, gene_order_csv)
pid_samecontig_plot      <- enforce_gene_order(pid_samecontig_plot, gene_order_csv)
pid_inrange_plot         <- enforce_gene_order(pid_inrange_plot, gene_order_csv)

#Column labels
gene_labels_map <- cps_annot %>%
  distinct(query_id, proper_gene) %>%
  deframe()

make_col_labels <- function(mat, gene_labels_map) {
  cn <- colnames(mat)
  labels <- cn
  hit <- cn %in% names(gene_labels_map)
  labels[hit] <- gene_labels_map[cn[hit]]
  labels
}

col_labels_anywhere   <- make_col_labels(pid_matrix_anywhere_plot, gene_labels_map)
col_labels_samecontig <- make_col_labels(pid_samecontig_plot, gene_labels_map)
col_labels_inrange    <- make_col_labels(pid_inrange_plot, gene_labels_map)

#accession ordering + wrapped group names:

prepare_rows <- function(mat, acc_annot, wrap_width = 10) {
  
  acc_order <- acc_annot %>%
    left_join(species_hex, by = "species") %>%
    distinct(accession, species, group, hex) %>%
    arrange(group, species) %>%
    filter(accession %in% rownames(mat))
  
  mat <- mat[acc_order$accession, , drop = FALSE]
  
  row_split <- factor(
    str_wrap(acc_order$group, width = wrap_width),
    levels = unique(str_wrap(acc_order$group, width = wrap_width))
  )
  
  row_ha <- rowAnnotation(
    Species = acc_order$species,
    col = list(Species = setNames(acc_order$hex, acc_order$species)),
    show_annotation_name = FALSE
  )
  
  list(mat = mat, row_split = row_split, row_ha = row_ha)
}

anywhere_rows   <- prepare_rows(pid_matrix_anywhere_plot, acc_annot)
samecontig_rows <- prepare_rows(pid_samecontig_plot, acc_annot)
inrange_rows    <- prepare_rows(pid_inrange_plot, acc_annot)

#colour schemes
col_fun <- colorRamp2(
  c(50, 70, 90, 100),
  c("#f7fbff", "#6baed6", "#2171b5", "#08306b")
)

col_fun <- colorRamp2(
  c(0, 10, 20, 30, 40, 50, 100),
  c("white", "lightyellow", "yellow", "orange", "orangered", "red", "darkred"))


column_grid_fun <- function(j, i, x, y, width, height, fill) {
  grid.lines(
    x = unit(c(x - width/2, x - width/2), "npc"),
    y = unit(c(y - height/2, y + height/2), "npc"),
    gp = gpar(col = "grey85", lwd = 0.1)
  )
}

column_grid_fun <- function(j, i, x, y, width, height, fill) {
  if (i == 1) {
    grid.lines(
      x = unit(c(x - width/2, x - width/2), "npc"),
      y = unit(c(0, 1), "npc"),
      gp = gpar(col = "grey85", lwd = 0.1)
    )
  }
}

column_grid_fun <- function(j, i, x, y, width, height, fill) {
  grid.lines(
    x = unit(c(x - width/2, x - width/2), "npc"),
    y = unit(c(0, 1), "npc"),
    gp = gpar(col = "grey85", lwd = 0.1)
  )
}

#anywhere in the genome
Heatmap(
  anywhere_rows$mat,
  name = "Percentage Identity",
  col = col_fun,
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  row_split = anywhere_rows$row_split,
  show_row_names = FALSE,
  show_column_names = TRUE,
  column_labels = col_labels_anywhere,
  column_names_gp = gpar(fontsize = 6),
  left_annotation = anywhere_rows$row_ha,
  row_title_gp = gpar(fontsize = 5),
  na_col = "white",
  heatmap_legend_param = list(
    title = "Percentage identity (%)",
    title_gp = gpar(fontsize = 6, fontface = "bold"),
    labels_gp = gpar(fontsize = 8)
  ),
  column_title = "Percentage identity of capsular genes anywhere in the genome",
  column_title_gp = gpar(fontsize = 16, fontface = "bold"),
  cell_fun = column_grid_fun
)
