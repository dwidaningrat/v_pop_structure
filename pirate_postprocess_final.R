setwd("/Users/waodedwidaningrat/viridans/pirate_31July2024/results_pirate_rmduplicate_4/summary_pirate_4/analysis_in_R/")

if (!requireNamespace("ggvenn", quietly = TRUE)) install.packages("ggvenn")
if (!requireNamespace("ggplot2", quietly = TRUE)) install.packages("ggplot2")
if (!requireNamespace("stringr", quietly = TRUE)) install.packages("stringr")
if (!requireNamespace("eulerr", quietly = TRUE)) install.packages("eulerr")
if (!requireNamespace("ComplexUpset", quietly = TRUE)) install.packages("ComplexUpset")

packageVersion("ggvenn")
packageDescription("ggvenn")[c("Package","Version","Title","Maintainer","Author")]


library(ComplexUpset)
library(ggvenn)
library(ggplot2)
library(stringr)
library(eulerr)
library(data.table)
library(dplyr)
library(openxlsx)

# Input
pirate_csv    <- "PIRATE_gene_families_forR_tree_accession.csv"

mapping_csv  <- "accession_species_group_tree_accession.csv"
species_hex_csv <- "species_hex_codes.csv"
group_hex_csv <- "group_hex.csv"


# Read Mapping
map_df <- fread(mapping_csv, data.table = FALSE)
names(map_df) <- tolower(names(map_df))
stopifnot(all(c("accession","species","group") %in% names(map_df)))

map_df$accession <- trimws(as.character(map_df$accession))
map_df$species   <- trimws(as.character(map_df$species))
map_df$group     <- trimws(as.character(map_df$group))
map_df <- map_df %>% filter(!is.na(accession), !is.na(species), !is.na(group)) %>% distinct()

# Read PIRATE results
hdr <- names(fread(pirate_csv, nrows = 0, check.names = FALSE))
hdr[1] <- "gene"

mapped_genomes <- intersect(hdr[-1], map_df$accession)
if (length(mapped_genomes) == 0) stop("No accession IDs in PIRATE results matched mapping data.")

pir <- fread(
  pirate_csv,
  select = c(hdr[1], mapped_genomes),
  data.table = FALSE,
  check.names = FALSE
)
names(pir)[1] <- "gene"
pir$gene <- as.character(pir$gene)

# Presence/absence matrix
presence_mat <- as.data.frame(
  lapply(pir[, -1, drop = FALSE], function(x) {
    x <- as.character(x)
    !is.na(x) & nzchar(x)
  }),
  check.names = FALSE
)
colnames(presence_mat) <- colnames(pir)[-1]
presence_mat$gene <- pir$gene

# Species -> genomes (accessions) 
species_list <- split(map_df$accession, map_df$species)
species_list <- lapply(species_list, function(acc) intersect(acc, mapped_genomes))
species_list <- species_list[lengths(species_list) > 0]

# Species -> gene-family set
gene_sets_species <- lapply(names(species_list), function(sp) {
  acc <- species_list[[sp]]
  present_any <- rowSums(presence_mat[, acc, drop = FALSE]) > 0
  presence_mat$gene[present_any]
})
names(gene_sets_species) <- names(species_list)


# Read species hex colours
col_df <- fread(species_hex_csv, data.table = FALSE)
nms <- tolower(names(col_df))
sp_i <- which(nms %in% c("species","sp","taxon") | grepl("species|taxon", nms))[1]
hx_i <- which(nms %in% c("hex","color","colour") | grepl("hex|colou?r|color", nms))[1]
if (is.na(sp_i) || is.na(hx_i)) stop("Could not detect species/hex columns. Headers are: ",
                                     paste(names(col_df), collapse = ", "))

col_df <- col_df[, c(sp_i, hx_i)]
names(col_df) <- c("species", "hex")
col_df$species <- trimws(as.character(col_df$species))
col_df$hex     <- trimws(as.character(col_df$hex))

get_species_cols <- function(species_vec) {
  cols <- col_df$hex[match(species_vec, col_df$species)]
  cols[is.na(cols)] <- "#BDBDBD"
  cols
}



# helper
sets_to_membership_df <- function(sets) {
  all_genes <- sort(unique(unlist(sets, use.names = FALSE)))
  mem <- as.data.frame(sapply(sets, function(v) all_genes %in% v))
  mem$gene <- all_genes
  mem
}

# group -> vector of species
group_to_species <- map_df %>%
  distinct(group, species) %>%
  group_by(group) %>%
  summarise(species = list(sort(unique(species))), .groups = "drop")

group_to_species %>%
  mutate(
    n_species = lengths(species),
    n_nonempty = vapply(species, function(spv) sum(lengths(gene_sets_species[spv]) > 0), integer(1))
  ) %>%
  select(group, n_species, n_nonempty) %>%
  print()

grp_col_df <- fread(group_hex_csv, data.table = FALSE)
nms <- tolower(names(grp_col_df))
gr_i <- which(nms %in% c("group","grp") | grepl("group", nms))[1]
hx_i <- which(nms %in% c("hex","color","colour") | grepl("hex|colou?r|color", nms))[1]
if (is.na(gr_i) || is.na(hx_i)) stop("Could not detect group/hex columns in group_hex.csv. Headers are: ",
                                     paste(names(grp_col_df), collapse = ", "))

grp_col_df <- grp_col_df[, c(gr_i, hx_i)]
names(grp_col_df) <- c("group", "hex")
grp_col_df$group <- trimws(as.character(grp_col_df$group))
grp_col_df$hex   <- trimws(as.character(grp_col_df$hex))

get_group_cols <- function(group_vec) {
  cols <- grp_col_df$hex[match(group_vec, grp_col_df$group)]
  cols[is.na(cols)] <- "#BDBDBD"
  cols
}

fmtN <- function(x) format(x, big.mark = ",", scientific = FALSE, trim = TRUE)

make_quant_labels <- function(fit, min_percent = 3, digits = 1) {
  vals <- fit$original.values
  pct  <- 100 * vals / sum(vals)
  
  count_str <- ifelse(abs(vals - round(vals)) < 1e-9,
                      as.character(as.integer(round(vals))),
                      format(vals, trim = TRUE))
  
  pct_str <- sprintf(paste0("%.", digits, "f"), pct)
  
  lbl <- ifelse(pct >= min_percent, paste0(count_str, " [", pct_str, "%]"), "")
  names(lbl) <- names(vals)
  lbl
}

# plot venn (2–8) or UpSet (>8)
plot_group_species <- function(sets, set_cols, title_text) {
  sets <- as.list(sets)
  
  # drop empty sets (safety)
  keep <- lengths(sets) > 0
  sets <- sets[keep]
  set_cols <- set_cols[keep]
  
  set_names <- names(sets)
  k <- length(set_names)
  cat(title_text, " | k =", k, "\n")
  if (k < 2) return(invisible(NULL))
  
  # membership df (genes x sets)
  mem <- sets_to_membership_df(sets)  # cols = set_names + gene
  totalN <- nrow(mem)
  
  # set sizes (per species)
  set_sizes <- lengths(sets)
  
  # wrap + append N to labels
  wrapped <- stringr::str_wrap(set_names, width = 14)
  wrappedN <- paste0(wrapped, "\nN=", fmtN(set_sizes))
  
  # rename membership columns to wrappedN
  names(mem)[match(set_names, names(mem))] <- wrappedN
  
  if (k <= 5) {
    p <- ggvenn::ggvenn(
      mem,
      columns = wrappedN,
      fill_color = set_cols,
      fill_alpha = 1,
      show_percentage = TRUE,
      digits = 1,
      set_name_size = 5,
      text_size = 4,
      padding = 0.12
    ) +
      ggplot2::labs(
        title = title_text,
        subtitle = paste0("Total N = ", fmtN(totalN))
      ) +
      ggplot2::theme(
        plot.title    = ggplot2::element_text(face = "bold", margin = ggplot2::margin(b = 2)),
        plot.subtitle = ggplot2::element_text(margin = ggplot2::margin(b = 4)),
        plot.margin   = ggplot2::margin(2, 2, 2, 2)
      )
    print(p)
    
  } else {
    fit <- eulerr::euler(mem[, wrappedN, drop = FALSE], shape = "ellipse")
    q_labs <- make_quant_labels(fit, min_percent = 3, digits = 1)
    
    legend_labels <- wrappedN  # already has "Species\nN=..."
    
    p <- plot(
      fit,
      fills = list(fill = setNames(set_cols, wrappedN), alpha = 1),
      quantities = q_labs,          # your "count [x%]" with <3% hidden
      labels = FALSE,               # no species names on ellipses
      legend = list(labels = legend_labels, ncol = 2),
      main = paste0(title_text, "\nTotal N = ", fmtN(totalN))
    )
    print(p)
  }
}




# 8 diagrams: within each GROUP, sets are SPECIES

# species -> hex
get_species_cols <- function(species_vec) {
  cols <- col_df$hex[match(species_vec, col_df$species)]
  cols[is.na(cols)] <- "#BDBDBD"
  cols
}

# count how many species per group
print(group_to_species %>% mutate(n_species = lengths(species)) %>% select(group, n_species))

for (g in group_to_species$group) {
  sp_vec <- group_to_species$species[[which(group_to_species$group == g)]]
  sp_vec <- intersect(sp_vec, names(gene_sets_species))  # safety
  
  sets_g <- gene_sets_species[sp_vec]
  cols_g <- get_species_cols(names(sets_g))
  
  plot_group_species(
    sets = sets_g,
    set_cols = cols_g,
    title_text = paste0("Gene-family overlap within group: ", g)
  )
}


# one diagram: sets are GROUPS
group_names_raw <- sort(unique(map_df$group))

# gene-family sets per group (union across species in group)
gene_sets_group <- lapply(group_names_raw, function(g) {
  sp_vec <- sort(unique(map_df$species[map_df$group == g]))
  sp_vec <- intersect(sp_vec, names(gene_sets_species))
  unique(unlist(gene_sets_species[sp_vec], use.names = FALSE))
})
names(gene_sets_group) <- group_names_raw

# colours from group_hex.csv
group_cols_raw <- get_group_cols(group_names_raw)

names(gene_sets_group) <- group_names_raw

group_cols_raw <- get_group_cols(group_names_raw)

# totals
group_sizes <- lengths(gene_sets_group)
all_genes_g <- sort(unique(unlist(gene_sets_group, use.names = FALSE)))
totalN_g <- length(all_genes_g)

# wrap AFTER colours/sizes are matched
grp_wrapped <- stringr::str_wrap(group_names_raw, width = 14)
grp_wrappedN <- paste0(grp_wrapped, "\nN=", fmtN(group_sizes))

names(gene_sets_group) <- grp_wrappedN
group_cols_wrapped <- setNames(group_cols_raw, grp_wrappedN)

# membership table
mem_g <- as.data.frame(sapply(gene_sets_group, function(v) all_genes_g %in% v))
names(mem_g) <- grp_wrappedN

fit_g <- eulerr::euler(mem_g, shape = "ellipse")
q_labs_g <- make_quant_labels(fit_g, min_percent = 3, digits = 1)

p <- plot(
  fit_g,
  fills = list(fill = group_cols_wrapped, alpha = 1),
  quantities = q_labs_g,
  labels = FALSE,
  legend = list(labels = grp_wrappedN, ncol = 2),
  main = paste0("Gene-family overlap by GROUP (union across species)\nTotal N = ", fmtN(totalN_g))
)
print(p)






# tables: group/species + exclusives

# all genomes columns present in presence_mat
genomes <- setdiff(colnames(presence_mat), "gene")
stopifnot(length(genomes) > 0)


genes_union_over_genomes <- function(acc_vec, presence_mat, genomes_all) {
  acc_vec <- intersect(acc_vec, genomes_all)
  if (length(acc_vec) == 0) return(character(0))
  present_any <- rowSums(as.matrix(presence_mat[, acc_vec, drop = FALSE])) > 0
  presence_mat$gene[present_any]
}

# helper: intersection gene set over a set of genomes (core genes)
genes_intersection_over_genomes <- function(acc_vec, presence_mat, genomes_all) {
  acc_vec <- intersect(acc_vec, genomes_all)
  if (length(acc_vec) == 0) return(character(0))
  present_all <- rowSums(as.matrix(presence_mat[, acc_vec, drop = FALSE])) == length(acc_vec)
  presence_mat$gene[present_all]
}


# build gene sets per species (union across genomes in species)
map_df2 <- map_df %>%
  mutate(
    accession = trimws(as.character(accession)),
    species   = trimws(as.character(species)),
    group     = trimws(as.character(group))
  ) %>%
  filter(!is.na(accession), !is.na(species), !is.na(group)) %>%
  distinct()

species_list <- split(map_df2$accession, map_df2$species)
gene_sets_species2 <- lapply(species_list, genes_union_over_genomes, presence_mat = presence_mat, genomes_all = genomes)
names(gene_sets_species2) <- names(species_list)

# build gene sets per group (union across genomes in group)
group_list <- split(map_df2$accession, map_df2$group)
gene_sets_group2 <- lapply(group_list, genes_union_over_genomes, presence_mat = presence_mat, genomes_all = genomes)
names(gene_sets_group2) <- names(group_list)

# build CORE gene sets per species (present in ALL genomes of species)
core_gene_sets_species <- lapply(
  species_list,
  genes_intersection_over_genomes,
  presence_mat = presence_mat,
  genomes_all = genomes
)
names(core_gene_sets_species) <- names(species_list)


# helper: >= X% presence gene set over a set of genomes
genes_fraction_over_genomes <- function(acc_vec, presence_mat, genomes_all, frac = 0.95) {
  acc_vec <- intersect(acc_vec, genomes_all)
  if (length(acc_vec) == 0) return(character(0))
  min_n <- ceiling(frac * length(acc_vec))
  present_n <- rowSums(as.matrix(presence_mat[, acc_vec, drop = FALSE]))
  presence_mat$gene[present_n >= min_n]
}


# build CORE gene sets per group (present in ALL genomes of group)
core_gene_sets_group <- lapply(
  group_list,
  genes_intersection_over_genomes,
  presence_mat = presence_mat,
  genomes_all = genomes
)
names(core_gene_sets_group) <- names(group_list)

# counts tables
genes_per_group <- tibble(
  group = names(gene_sets_group2),
  n_genes_union = lengths(gene_sets_group2)
) %>% arrange(desc(n_genes_union), group)

genes_per_species <- tibble(
  species = names(gene_sets_species2),
  n_genes_union = lengths(gene_sets_species2)
) %>% arrange(desc(n_genes_union), species)

# exclusives: present in entity, absent from all other entities
exclusive_count_table <- function(gene_sets_named, entity_col = "entity") {
  entities <- names(gene_sets_named)
  all_union <- unique(unlist(gene_sets_named, use.names = FALSE))
  
  excl_counts <- lapply(entities, function(e) {
    others <- setdiff(entities, e)
    others_union <- unique(unlist(gene_sets_named[others], use.names = FALSE))
    excl <- setdiff(gene_sets_named[[e]], others_union)
    tibble(
      !!entity_col := e,
      n_exclusive = length(excl),
      n_total_union = length(gene_sets_named[[e]]),
      n_shared = length(gene_sets_named[[e]]) - length(excl)
    )
  }) %>% bind_rows()
  
  excl_counts %>% arrange(desc(n_exclusive), .data[[entity_col]])
}

exclusive_by_group <- exclusive_count_table(gene_sets_group2, entity_col = "group")
exclusive_by_species <- exclusive_count_table(gene_sets_species2, entity_col = "species")

# >=95% gene sets per species
softcore_gene_sets_species <- lapply(
  species_list,
  genes_fraction_over_genomes,
  presence_mat = presence_mat,
  genomes_all = genomes,
  frac = 0.95
)
names(softcore_gene_sets_species) <- names(species_list)

# >=95% gene sets per group
softcore_gene_sets_group <- lapply(
  group_list,
  genes_fraction_over_genomes,
  presence_mat = presence_mat,
  genomes_all = genomes,
  frac = 0.95
)
names(softcore_gene_sets_group) <- names(group_list)

# core gene distributions

# genes shared by ALL genomes within each species
core_genes_per_species <- tibble(
  species = names(core_gene_sets_species),
  n_core_genes = lengths(core_gene_sets_species)
) %>% arrange(desc(n_core_genes), species)

# genes shared by ALL genomes within each group
core_genes_per_group <- tibble(
  group = names(core_gene_sets_group),
  n_core_genes = lengths(core_gene_sets_group)
) %>% arrange(desc(n_core_genes), group)

# genes shared by ALL genomes across the entire dataset
core_genes_all_genomes <- genes_intersection_over_genomes(
  acc_vec = genomes,
  presence_mat = presence_mat,
  genomes_all = genomes
)

core_genes_dataset <- tibble(
  scope = "all_genomes",
  n_core_genes = length(core_genes_all_genomes)
)

# >=95% (soft-core) gene distributions

# genes shared by >=95% of genomes within each species
softcore_genes_per_species <- tibble(
  species = names(softcore_gene_sets_species),
  n_softcore_genes = lengths(softcore_gene_sets_species)
) %>% arrange(desc(n_softcore_genes), species)

# genes shared by >=95% of genomes within each group
softcore_genes_per_group <- tibble(
  group = names(softcore_gene_sets_group),
  n_softcore_genes = lengths(softcore_gene_sets_group)
) %>% arrange(desc(n_softcore_genes), group)

# genes shared by >=95% of genomes across the entire dataset
softcore_genes_all_genomes <- genes_fraction_over_genomes(
  acc_vec = genomes,
  presence_mat = presence_mat,
  genomes_all = genomes,
  frac = 0.95
)

softcore_genes_dataset <- tibble(
  scope = "all_genomes",
  threshold = ">=95%",
  n_softcore_genes = length(softcore_genes_all_genomes)
)

# write workbook
out_xlsx <- "gene_presence_distributions.xlsx"
wb <- createWorkbook()

addWorksheet(wb, "1_genes_per_group")
writeDataTable(wb, "1_genes_per_group", genes_per_group)

addWorksheet(wb, "2_genes_per_species")
writeDataTable(wb, "2_genes_per_species", genes_per_species)

addWorksheet(wb, "3_exclusive_by_group")
writeDataTable(wb, "3_exclusive_by_group", exclusive_by_group)

addWorksheet(wb, "4_exclusive_by_species")
writeDataTable(wb, "4_exclusive_by_species", exclusive_by_species)

addWorksheet(wb, "5_core_genes_per_species")
writeDataTable(wb, "5_core_genes_per_species", core_genes_per_species)

addWorksheet(wb, "6_core_genes_per_group")
writeDataTable(wb, "6_core_genes_per_group", core_genes_per_group)

addWorksheet(wb, "7_core_genes_all_genomes")
writeDataTable(wb, "7_core_genes_all_genomes", core_genes_dataset)

addWorksheet(wb, "8_softcore_genes_per_species")
writeDataTable(wb, "8_softcore_genes_per_species", softcore_genes_per_species)

addWorksheet(wb, "9_softcore_genes_per_group")
writeDataTable(wb, "9_softcore_genes_per_group", softcore_genes_per_group)

addWorksheet(wb, "10_softcore_genes_all_genomes")
writeDataTable(wb, "10_softcore_genes_all_genomes", softcore_genes_dataset)

saveWorkbook(wb, out_xlsx, overwrite = TRUE)

message("Wrote: ", out_xlsx)


#SANITY CHECK
# soft-core should be >= core
stopifnot(all(softcore_genes_per_species$n_softcore_genes >= core_genes_per_species$n_core_genes))
stopifnot(all(softcore_genes_per_group$n_softcore_genes >= core_genes_per_group$n_core_genes))
stopifnot(softcore_genes_dataset$n_softcore_genes >= core_genes_dataset$n_core_genes)

# core should be <= union
stopifnot(all(core_genes_per_species$n_core_genes <= genes_per_species$n_genes_union))
stopifnot(all(core_genes_per_group$n_core_genes <= genes_per_group$n_genes_union))






# Heatmaps

# packages
if (!requireNamespace("ape", quietly = TRUE)) install.packages("ape")
if (!requireNamespace("data.table", quietly = TRUE)) install.packages("data.table")
if (!requireNamespace("dplyr", quietly = TRUE)) install.packages("dplyr")
if (!requireNamespace("ggplot2", quietly = TRUE)) install.packages("ggplot2")
if (!requireNamespace("stringr", quietly = TRUE)) install.packages("stringr")

library(ape); library(data.table); library(dplyr); library(ggplot2); library(stringr)

# inputs (edit if needed)
tree_file <- "core_genome_tree_microreact_midpoint_decreasing.nwk"
mapping_csv <- "accession_species_group_tree_accession.csv"  
order_csv   <- "species_group_order.csv"

# checks & read mapping
stopifnot(exists("presence_mat"))
stopifnot("gene" %in% colnames(presence_mat))

map_df_file <- data.table::fread(mapping_csv, data.table = FALSE)
names(map_df_file) <- tolower(names(map_df_file))

# detection of columns
stopifnot(any(c("accession") %in% names(map_df_file)))
stopifnot(any(c("species")     %in% names(map_df_file)))
stopifnot(any(c("group")       %in% names(map_df_file)))

# normalize
map_df2 <- map_df_file %>%
  transmute(
    accession = trimws(as.character(accession)),
    species   = trimws(as.character(species)),
    group     = trimws(as.character(group))
  ) %>%
  filter(!is.na(accession), !is.na(species), !is.na(group)) %>%
  distinct()

# genomes present in presence_mat (columns)
genomes_all <- setdiff(colnames(presence_mat), "gene")
if (length(genomes_all) == 0) stop("No genome/accession columns found in presence_mat.")

# restrict mapping to only genomes present
map_df2 <- map_df2 %>% filter(accession %in% genomes_all)

# tree order helper
if (!requireNamespace("ape", quietly = TRUE)) install.packages("ape")
if (!requireNamespace("phytools", quietly = TRUE)) install.packages("phytools")

library(ape)
library(phytools)

get_tree_order_midpoint <- function(tree_path, genomes_vec) {
  tr <- ape::read.tree(tree_path)
  
  # Midpoint rooting needs edge lengths; if missing, give all edges length 1
  if (is.null(tr$edge.length)) {
    tr$edge.length <- rep(1, nrow(tr$edge))
  }
  
  # midpoint root, then ladderize
  tr <- phytools::midpoint.root(tr)
  tr <- ape::ladderize(tr, right = TRUE)
  
  tips <- tr$tip.label
  
  in_tree   <- intersect(tips, genomes_vec)
  not_in_tr <- setdiff(genomes_vec, tips)
  
  c(in_tree, sort(not_in_tr))
}

genomes_all <- setdiff(colnames(presence_mat), "gene")

tree_order_all <- get_tree_order_midpoint(tree_file, genomes_all)

tree_rank <- setNames(seq_along(tree_order_all), tree_order_all)


# species/group order
read_species_group_order_auto <- function(path, map_df_reference) {
  raw <- tryCatch(data.table::fread(path, data.table = FALSE, stringsAsFactors = FALSE),
                  error = function(e) stop("Cannot read order CSV: ", e$message))
  nms <- tolower(names(raw))
  sp_i <- which(nms %in% c("species","sp","taxon"))[1]
  gr_i <- which(nms %in% c("group","grp","clade"))[1]
  if (is.na(sp_i) || is.na(gr_i)) {
    ref_sp <- unique(as.character(map_df_reference$species))
    ref_gr <- unique(as.character(map_df_reference$group))
    col_vals <- lapply(raw, function(x) unique(na.omit(as.character(x[x != ""]))))
    sp_scores <- sapply(col_vals, function(vals) mean(vals %in% ref_sp))
    gr_scores <- sapply(col_vals, function(vals) mean(vals %in% ref_gr))
    sp_i <- which.max(sp_scores); gr_i <- which.max(gr_scores)
    if (sp_i == gr_i && length(col_vals) > 1) gr_i <- order(gr_scores, decreasing = TRUE)[2]
  }
  order_df <- data.frame(group = as.character(raw[[gr_i]]), species = as.character(raw[[sp_i]]), stringsAsFactors = FALSE)
  order_df <- order_df[!(is.na(order_df$group) & is.na(order_df$species)), , drop = FALSE]
  order_df$row_order <- seq_len(nrow(order_df))
  list(order_df = order_df, group_levels = unique(order_df$group), species_levels = unique(order_df$species))
}
ord <- NULL
if (file.exists(order_csv)) {
  ord <- read_species_group_order_auto(order_csv, map_df2)
  order_df <- ord$order_df
  group_levels <- ord$group_levels
} else {
  order_df <- data.frame(group = unique(map_df2$group), species = NA_character_, row_order = seq_along(unique(map_df2$group)), stringsAsFactors = FALSE)
  group_levels <- unique(map_df2$group)
}

# ---- order genomes by CSV (group -> species), tree as tiebreak ----
order_genomes_by_csv <- function(map_df_ref, genomes_vec, order_df, tree_order = NULL) {
  m <- map_df_ref %>% filter(accession %in% genomes_vec)
  key <- order_df %>% group_by(group, species) %>% summarise(row_order = min(row_order, na.rm = TRUE), .groups = "drop")
  m2 <- left_join(m, key, by = c("group","species"))
  m2$row_order[is.na(m2$row_order)] <- 1e9
  if (!is.null(tree_order)) {
    tree_rank <- setNames(seq_along(tree_order), tree_order)
    m2$tree_r <- ifelse(m2$accession %in% names(tree_rank), tree_rank[m2$accession], 1e9)
  } else {
    m2$tree_r <- 1e9
  }
  m2 <- m2 %>% arrange(row_order, tree_r, accession)
  missing_map <- setdiff(genomes_vec, m2$accession)
  c(m2$accession, sort(missing_map))
}

# gene order by prevalence in a genome subset (returns only genes with >0 presence)
get_gene_order <- function(presence_mat, genomes_subset) {
  genomes_subset <- intersect(genomes_subset, setdiff(colnames(presence_mat), "gene"))
  if (length(genomes_subset) == 0) return(character(0))
  prev <- rowSums(as.matrix(presence_mat[, genomes_subset, drop = FALSE]))
  present_any <- prev > 0
  genes_present <- presence_mat$gene[present_any]
  prev_present <- prev[present_any]
  genes_present[order(prev_present, genes_present, decreasing = TRUE)]
}

# helper to create non-repeating species labels for y axis
# For a given genome_order (vector of accessions), returns a named vector of labels (one per accession)
# Labels are empty "" for most accessions; only one accession per species gets the label (the first in genomic order).
make_species_labels_for_axis <- function(genome_order, map_df_ref) {
  # genome_order: vector of accession IDs in the plotted order (top->bottom order should be reversed outside)
  # map_df_ref: mapping df with accession,species,group
  species_for_acc <- map_df_ref$species[match(genome_order, map_df_ref$accession)]
  group_for_acc   <- map_df_ref$group[match(genome_order, map_df_ref$accession)]
  labels <- rep("", length(genome_order))
  names(labels) <- genome_order
  # find first occurrence per species
  uniq_sp <- unique(species_for_acc[!is.na(species_for_acc)])
  for (sp in uniq_sp) {
    idx <- which(species_for_acc == sp)
    if (length(idx) == 0) next
    first_i <- idx[1]
    grp <- group_for_acc[first_i]
    labels[first_i] <- ifelse(is.na(grp) || grp == "", sp, paste0(sp, " (", grp, ")"))
  }
  labels
}

# full label per accession
make_accession_labels <- function(accession_levels, map_df_ref) {
  sp <- map_df_ref$species[match(accession_levels, map_df_ref$accession)]
  gr <- map_df_ref$group[match(accession_levels, map_df_ref$accession)]
  lbl <- paste0(accession_levels, " | ", sp, " | ", gr)
  lbl[is.na(sp)] <- accession_levels[is.na(sp)]  # fallback
  names(lbl) <- accession_levels
  lbl
}

# y positions where species changes (for horizontal lines)
species_break_positions <- function(genome_order, map_df_ref) {
  sp <- map_df_ref$species[match(genome_order, map_df_ref$accession)]
  # break after each run of species
  idx <- which(sp[-1] != sp[-length(sp)])
  idx + 0.5
}

# y positions between every accession row (species plots)
accession_break_positions <- function(genome_order) {
  if (length(genome_order) <= 1) return(numeric(0))
  seq_len(length(genome_order) - 1) + 0.5
}


# plot function
plot_presence_heatmap2 <- function(presence_mat, genomes_subset, genome_order, gene_order,
                                   title_text, out_file, show_gene_labels = FALSE, map_df_ref) {
  genomes_subset <- intersect(genomes_subset, setdiff(colnames(presence_mat), "gene"))
  if (length(genomes_subset) == 0) return(invisible(NULL))
  
  # FORCE tree order
  genome_order <- intersect(genome_order, genomes_subset)
  if (length(genome_order) == 0) return(invisible(NULL))
  
  gene_order <- intersect(gene_order, presence_mat$gene)
  if (length(gene_order) == 0) return(invisible(NULL))
  
  dt <- as.data.table(presence_mat[, c("gene", genome_order), drop = FALSE])
  for (g in genome_order) dt[, (g) := as.logical(get(g))]
  long_dt <- melt(dt, id.vars = "gene", variable.name = "accession", value.name = "present")
  
  long_dt[, gene := factor(gene, levels = gene_order)]
  long_dt[, accession := factor(accession, levels = rev(genome_order))]
  
  # y labels: one species label per species (first row of species)
  labels_for_genome_order <- make_species_labels_for_axis(genome_order, map_df_ref)
  plot_levels <- rev(genome_order)
  y_labels <- labels_for_genome_order[plot_levels]
  y_labels[is.na(y_labels)] <- ""
  names(y_labels) <- plot_levels
  
  # horizontal lines between species
  breaks <- species_break_positions(genome_order, map_df_ref)
  
  n <- length(genome_order)
  breaks_plot <- n - breaks + 1
  
  p <- ggplot(long_dt, aes(x = gene, y = accession, fill = present)) +
    geom_raster() +
    scale_fill_manual(values = c(`TRUE` = "#08306B", `FALSE` = "white"), na.value = "white") +
    scale_y_discrete(labels = y_labels) +
    labs(
      title = title_text,
      subtitle = paste0("Genomes: ", length(genome_order), " | Genes: ", length(gene_order)),
      x = "Genes (sorted by prevalence within plotted subset)",
      y = "",
      fill = "Present"
    ) +
    theme_minimal(base_size = 11) +
    theme(
      panel.grid = element_blank(),
      axis.text.y = element_text(size = 7, hjust = 1),
      axis.title.x = element_text(margin = margin(t = 6)),
      plot.subtitle = element_text(size = 9),
      plot.title = element_text(face = "bold")
    )
  
  if (!show_gene_labels) {
    p <- p + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
  } else {
    p <- p + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 3))
  }
  
  n_genes <- length(gene_order)
  n_genomes <- length(genome_order)
  w_in <- max(10, min(70, n_genes / 150))
  h_in <- max(6,  min(45, n_genomes / 25))
  
  ggsave(out_file, p, width = w_in, height = h_in, dpi = 300, limitsize = FALSE)
  invisible(p)
}


# output
dir.create("heatmaps", showWarnings = FALSE)
dir.create(file.path("heatmaps", "groups"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path("heatmaps", "species"), showWarnings = FALSE, recursive = TRUE)

# Plot 1: ALL genomes, y-order by tree, x by prevalence
# For the all-genomes plots, keep all genes
# global gene order (all genes, including those absent in some subsets)
gene_order_global <- presence_mat$gene[order(rowSums(as.matrix(presence_mat[, genomes_all, drop = FALSE])), decreasing = TRUE)]

genome_order_p1 <- tree_order_all
gene_order_p1 <- gene_order_global

plot_presence_heatmap2(
  presence_mat = presence_mat,
  genomes_subset = genomes_all,
  genome_order = genome_order_p1,
  gene_order = gene_order_p1,
  title_text = "Plot 1: Presence/absence (genomes ordered by tree)",
  out_file = file.path("heatmaps", "plot1_tree_order_all_genomes.pdf"),
  show_gene_labels = FALSE,
  map_df_ref = map_df2
)

# Plot 2: ALL genomes, y-order by species_group_order.csv (group->species) within species: tree order (when available)
# Uses accession_species_group.csv mapping for labels & ordering

genome_order_p2 <- tree_order_all
gene_order_p2 <- gene_order_global

plot_presence_heatmap2(
  presence_mat = presence_mat,
  genomes_subset = genomes_all,
  genome_order = genome_order_p2,
  gene_order = gene_order_p2,
  title_text = "Plot 2: Presence/absence (genomes ordered by CSV species->group)",
  out_file = file.path("heatmaps", "plot2_species_group_order_all_genomes.pdf"),
  show_gene_labels = FALSE,
  map_df_ref = map_df2
)


# Group plots: one plot per group
# each plot uses ONLY genes that are present in that group's genomes
# - y-axis: species (one label per species; label format "species (group)")

groups_in_data <- sort(unique(map_df2$group))
groups_plot_order <- unique(c(intersect(if (exists("group_levels")) group_levels else character(0), groups_in_data), setdiff(groups_in_data, if (exists("group_levels")) group_levels else character(0))))

for (g in groups_plot_order) {
  acc_g <- map_df2$accession[map_df2$group == g]
  acc_g <- intersect(acc_g, genomes_all)
  if (length(acc_g) < 1) next
  
  # order genomes by CSV within this group (species order), then tree within species
  map_g <- map_df2 %>% filter(group == g, accession %in% acc_g)
  #genome_order_g <- order_genomes_by_csv(map_g, acc_g, order_df, tree_order = tree_order_all)
  genome_order_g <- intersect(tree_order_all, acc_g)
  
  
  # gene order restricted to genes present in this group's genomes
  gene_order_g <- get_gene_order(presence_mat, acc_g)
  if (length(gene_order_g) == 0) next
  
  safe_g <- gsub("[^A-Za-z0-9._-]+", "_", g)
  plot_presence_heatmap2(
    presence_mat = presence_mat,
    genomes_subset = acc_g,
    genome_order = genome_order_g,
    gene_order = gene_order_g,
    title_text = paste0("Group: ", g, " — genes present in group (N=", length(gene_order_g), ")"),
    out_file = file.path("heatmaps", "groups", paste0("group_", safe_g, "_presence_absence.pdf")),
    show_gene_labels = FALSE,
    map_df_ref = map_df2
  )
}

# Species plots: one plot per species
# each plot shows accession IDs on the y-axis (one per genome)
# each plot uses ONLY genes that are present in that species' genomes

species_in_data <- sort(unique(map_df2$species))

for (sp in species_in_data) {
  acc_sp <- map_df2$accession[map_df2$species == sp]
  acc_sp <- intersect(acc_sp, genomes_all)
  if (length(acc_sp) < 1) next
  
  # genome order: keep tree order for this species where possible, else accession sort
  #genome_order_sp <- c(intersect(tree_order_all, acc_sp), setdiff(acc_sp, intersect(tree_order_all, acc_sp)))
  genome_order_sp <- intersect(tree_order_all, acc_sp)
  
  # genes present in this species (union across its genomes), ordered by prevalence within species
  gene_order_sp <- get_gene_order(presence_mat, acc_sp)
  if (length(gene_order_sp) == 0) next
  
  safe_sp <- gsub("[^A-Za-z0-9._-]+", "_", sp)
  
  # function that keeps accession IDs as y-axis labels
  plot_presence_heatmap_species <- function(presence_mat, genomes_subset, genome_order, gene_order,
                                            title_text, out_file, show_gene_labels = FALSE, map_df_ref) {
    genomes_subset <- intersect(genomes_subset, setdiff(colnames(presence_mat), "gene"))
    if (length(genomes_subset) == 0) return(invisible(NULL))
    
    # FORCE tree order
    genome_order <- intersect(genome_order, genomes_subset)
    if (length(genome_order) == 0) return(invisible(NULL))
    
    gene_order <- intersect(gene_order, presence_mat$gene)
    if (length(gene_order) == 0) return(invisible(NULL))
    
    dt <- as.data.table(presence_mat[, c("gene", genome_order), drop = FALSE])
    for (g in genome_order) dt[, (g) := as.logical(get(g))]
    long_dt <- melt(dt, id.vars = "gene", variable.name = "accession", value.name = "present")
    
    long_dt[, gene := factor(gene, levels = gene_order)]
    long_dt[, accession := factor(accession, levels = rev(genome_order))]
    
    # full labels per accession
    plot_levels <- rev(genome_order)
    acc_labels <- make_accession_labels(plot_levels, map_df_ref)
    
    # horizontal line between every row
    n <- length(genome_order)
    breaks_plot <- accession_break_positions(plot_levels)  # already in plot order coords
    
    p <- ggplot(long_dt, aes(x = gene, y = accession, fill = present)) +
      geom_raster() +
      scale_fill_manual(values = c(`TRUE` = "#08306B", `FALSE` = "white"), na.value = "white") +
      scale_y_discrete(labels = acc_labels) +
      labs(
        title = title_text,
        subtitle = paste0("Genomes: ", length(genome_order), " | Genes: ", length(gene_order)),
        x = "Genes (sorted by prevalence within species)",
        y = "",
        fill = "Present"
      ) +
      theme_minimal(base_size = 11) +
      theme(
        panel.grid = element_blank(),
        axis.text.y = element_text(size = 6),
        axis.title.x = element_text(margin = margin(t = 6)),
        plot.subtitle = element_text(size = 9),
        plot.title = element_text(face = "bold")
      )
    
    if (!show_gene_labels) {
      p <- p + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
    } else {
      p <- p + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 3))
    }
    
    n_genes <- length(gene_order)
    n_genomes <- length(genome_order)
    w_in <- max(10, min(70, n_genes / 150))
    h_in <- max(6,  min(45, n_genomes / 25))
    
    ggsave(out_file, p, width = w_in, height = h_in, dpi = 300, limitsize = FALSE)
    invisible(p)
  }
  
  plot_presence_heatmap_species(
    presence_mat = presence_mat,
    genomes_subset = acc_sp,
    genome_order = genome_order_sp,
    gene_order = gene_order_sp,
    title_text = paste0("Species: ", sp, " — accessions (N=", length(acc_sp), "); genes present in species (N=", length(gene_order_sp), ")"),
    out_file = file.path("heatmaps", "species", paste0("species_", safe_sp, "_presence_absence_by_accession.pdf")),
    show_gene_labels = FALSE,
    map_df_ref = map_df2
  )
}


message("Done. Plots printed to the plotting device and PDFs written into 'heatmaps/' (groups & species directories).")













# Heatmaps (GROUP/SPECIES ordered; NO TREE)

# packages
if (!requireNamespace("data.table", quietly = TRUE)) install.packages("data.table")
if (!requireNamespace("dplyr", quietly = TRUE)) install.packages("dplyr")
if (!requireNamespace("ggplot2", quietly = TRUE)) install.packages("ggplot2")
if (!requireNamespace("stringr", quietly = TRUE)) install.packages("stringr")

library(data.table)
library(dplyr)
library(ggplot2)
library(stringr)

# inputs

mapping_csv <- "accession_species_group_tree_accession.csv"

# checks & read mapping
stopifnot(exists("presence_mat"))
stopifnot("gene" %in% colnames(presence_mat))

map_df_file <- data.table::fread(mapping_csv, data.table = FALSE)
names(map_df_file) <- tolower(names(map_df_file))
stopifnot(all(c("accession", "species", "group") %in% names(map_df_file)))

# normalize
map_df2 <- map_df_file %>%
  transmute(
    accession = trimws(as.character(accession)),
    species   = trimws(as.character(species)),
    group     = trimws(as.character(group))
  ) %>%
  filter(!is.na(accession), !is.na(species), !is.na(group), nzchar(accession), nzchar(species), nzchar(group)) %>%
  distinct(accession, .keep_all = TRUE)

# genomes present in presence_mat (columns)
genomes_all <- setdiff(colnames(presence_mat), "gene")
if (length(genomes_all) == 0) stop("No genome/accession columns found in presence_mat.")

# restrict mapping to only genomes present
map_df2 <- map_df2 %>% filter(accession %in% genomes_all)

# define group/species ranks by FIRST APPEARANCE in mapping file
group_levels <- unique(map_df2$group)

map_df2 <- map_df2 %>%
  mutate(group_rank = match(group, group_levels)) %>%
  group_by(group) %>%
  mutate(species_rank = match(species, unique(species))) %>%
  ungroup()

# global genome order: GROUP -> SPECIES -> accession (alphabetical within species)
genome_order_all <- map_df2 %>%
  arrange(group_rank, species_rank, accession) %>%
  pull(accession)

# append any genomes not in mapping at the end
genome_order_all <- c(genome_order_all, setdiff(genomes_all, genome_order_all))

# gene order helpers
get_gene_order <- function(presence_mat, genomes_subset) {
  genomes_subset <- intersect(genomes_subset, setdiff(colnames(presence_mat), "gene"))
  if (length(genomes_subset) == 0) return(character(0))
  prev <- rowSums(as.matrix(presence_mat[, genomes_subset, drop = FALSE]))
  present_any <- prev > 0
  genes_present <- presence_mat$gene[present_any]
  prev_present <- prev[present_any]
  genes_present[order(prev_present, genes_present, decreasing = TRUE)]
}

# global gene order (all genes; prevalence across ALL genomes)
gene_order_global <- presence_mat$gene[
  order(rowSums(as.matrix(presence_mat[, genomes_all, drop = FALSE])), decreasing = TRUE)
]

# axis label helpers

# hierarchical labels:
# show GROUP once at start of its block
# show SPECIES once per species
make_group_species_labels_for_axis <- function(genome_order, map_df_ref) {
  grp <- map_df_ref$group[match(genome_order, map_df_ref$accession)]
  sp  <- map_df_ref$species[match(genome_order, map_df_ref$accession)]
  
  labels <- rep("", length(genome_order))
  names(labels) <- genome_order
  
  # groups in the plotted order
  grp_order <- unique(grp[!is.na(grp)])
  
  for (g in grp_order) {
    idx_g <- which(grp == g)
    if (length(idx_g) == 0) next
    
    # species in the plotted order within group
    sp_order <- unique(sp[idx_g][!is.na(sp[idx_g])])
    
    for (j in seq_along(sp_order)) {
      s <- sp_order[j]
      idx_s <- idx_g[sp[idx_g] == s]
      if (length(idx_s) == 0) next
      
      first_row <- idx_s[1]
      if (j == 1) {
        # group label + first species in group on the same row
        labels[first_row] <- paste0(g, "\n  ", s)
      } else {
        labels[first_row] <- paste0("  ", s)
      }
    }
  }
  
  labels
}

# per-accession labels for species plots
make_accession_labels <- function(accession_levels, map_df_ref) {
  sp <- map_df_ref$species[match(accession_levels, map_df_ref$accession)]
  gr <- map_df_ref$group[match(accession_levels, map_df_ref$accession)]
  lbl <- paste0(accession_levels, " | ", sp, " | ", gr)
  lbl[is.na(sp)] <- accession_levels[is.na(sp)]  # fallback
  names(lbl) <- accession_levels
  lbl
}

# breaks
species_break_positions <- function(genome_order, map_df_ref) {
  sp <- map_df_ref$species[match(genome_order, map_df_ref$accession)]
  idx <- which(sp[-1] != sp[-length(sp)])
  idx + 0.5
}

accession_break_positions <- function(genome_order) {
  if (length(genome_order) <= 1) return(numeric(0))
  seq_len(length(genome_order) - 1) + 0.5
}

# plotting functions 

plot_presence_heatmap2 <- function(presence_mat, genomes_subset, genome_order, gene_order,
                                   title_text, out_file, show_gene_labels = FALSE, map_df_ref) {
  genomes_subset <- intersect(genomes_subset, setdiff(colnames(presence_mat), "gene"))
  if (length(genomes_subset) == 0) return(invisible(NULL))
  
  # enforce ordering (GROUP/SPECIES), by filtering genome_order
  genome_order <- genome_order[genome_order %in% genomes_subset]
  if (length(genome_order) == 0) return(invisible(NULL))
  
  gene_order <- intersect(gene_order, presence_mat$gene)
  if (length(gene_order) == 0) return(invisible(NULL))
  
  dt <- as.data.table(presence_mat[, c("gene", genome_order), drop = FALSE])
  for (g in genome_order) dt[, (g) := as.logical(get(g))]
  long_dt <- melt(dt, id.vars = "gene", variable.name = "accession", value.name = "present")
  
  long_dt[, gene := factor(gene, levels = gene_order)]
  long_dt[, accession := factor(accession, levels = rev(genome_order))]
  
  # hierarchical y labels: GROUP + species
  labels_for_order <- make_group_species_labels_for_axis(genome_order, map_df_ref)
  plot_levels <- rev(genome_order)
  y_labels <- labels_for_order[plot_levels]
  y_labels[is.na(y_labels)] <- ""
  names(y_labels) <- plot_levels
  
  # species breaks
  breaks <- species_break_positions(genome_order, map_df_ref)
  n <- length(genome_order)
  breaks_plot <- n - breaks + 1
  
  p <- ggplot(long_dt, aes(x = gene, y = accession, fill = present)) +
    geom_raster() +
    scale_fill_manual(values = c(`TRUE` = "#08306B", `FALSE` = "white"), na.value = "white") +
    scale_y_discrete(labels = y_labels) +
    labs(
      title = title_text,
      subtitle = paste0("Genomes: ", length(genome_order), " | Genes: ", length(gene_order)),
      x = "Genes (sorted by prevalence within plotted subset)",
      y = "",
      fill = "Present"
    ) +
    theme_minimal(base_size = 11) +
    theme(
      panel.grid = element_blank(),
      axis.text.y = element_text(size = 7, hjust = 1),
      axis.title.x = element_text(margin = margin(t = 6)),
      plot.subtitle = element_text(size = 9),
      plot.title = element_text(face = "bold")
    )
  
  if (!show_gene_labels) {
    p <- p + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
  } else {
    p <- p + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 3))
  }
  
  n_genes <- length(gene_order)
  n_genomes <- length(genome_order)
  w_in <- max(10, min(70, n_genes / 150))
  h_in <- max(6,  min(45, n_genomes / 25))
  
  ggsave(out_file, p, width = w_in, height = h_in, dpi = 300, limitsize = FALSE)
  invisible(p)
}

plot_presence_heatmap_species <- function(presence_mat, genomes_subset, genome_order, gene_order,
                                          title_text, out_file, show_gene_labels = FALSE, map_df_ref) {
  genomes_subset <- intersect(genomes_subset, setdiff(colnames(presence_mat), "gene"))
  if (length(genomes_subset) == 0) return(invisible(NULL))
  
  # enforce global group/species ordering
  genome_order <- genome_order[genome_order %in% genomes_subset]
  if (length(genome_order) == 0) return(invisible(NULL))
  
  gene_order <- intersect(gene_order, presence_mat$gene)
  if (length(gene_order) == 0) return(invisible(NULL))
  
  dt <- as.data.table(presence_mat[, c("gene", genome_order), drop = FALSE])
  for (g in genome_order) dt[, (g) := as.logical(get(g))]
  long_dt <- melt(dt, id.vars = "gene", variable.name = "accession", value.name = "present")
  
  long_dt[, gene := factor(gene, levels = gene_order)]
  long_dt[, accession := factor(accession, levels = rev(genome_order))]
  
  plot_levels <- rev(genome_order)
  acc_labels <- make_accession_labels(plot_levels, map_df_ref)
  
  # accession breaks
  breaks_plot <- accession_break_positions(plot_levels)
  
  p <- ggplot(long_dt, aes(x = gene, y = accession, fill = present)) +
    geom_raster() +
    scale_fill_manual(values = c(`TRUE` = "#08306B", `FALSE` = "white"), na.value = "white") +
    scale_y_discrete(labels = acc_labels) +
    labs(
      title = title_text,
      subtitle = paste0("Genomes: ", length(genome_order), " | Genes: ", length(gene_order)),
      x = "Genes (sorted by prevalence within species)",
      y = "",
      fill = "Present"
    ) +
    theme_minimal(base_size = 11) +
    theme(
      panel.grid = element_blank(),
      axis.text.y = element_text(size = 6),
      axis.title.x = element_text(margin = margin(t = 6)),
      plot.subtitle = element_text(size = 9),
      plot.title = element_text(face = "bold")
    )
  
  if (!show_gene_labels) {
    p <- p + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
  } else {
    p <- p + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 3))
  }
  
  n_genes <- length(gene_order)
  n_genomes <- length(genome_order)
  w_in <- max(10, min(70, n_genes / 150))
  h_in <- max(6,  min(45, n_genomes / 25))
  
  ggsave(out_file, p, width = w_in, height = h_in, dpi = 300, limitsize = FALSE)
  invisible(p)
}

# output
dir.create("heatmaps", showWarnings = FALSE)
dir.create(file.path("heatmaps", "groups"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path("heatmaps", "species"), showWarnings = FALSE, recursive = TRUE)

# Plot 1: ALL genomes (ordered by GROUP -> SPECIES from mapping)

plot_presence_heatmap2(
  presence_mat = presence_mat,
  genomes_subset = genomes_all,
  genome_order = genome_order_all,
  gene_order = gene_order_global,
  title_text = "Plot 1: Presence/absence (ordered by GROUP -> SPECIES from mapping CSV)",
  out_file = file.path("heatmaps", "plot1_group_species_order_all_genomes.pdf"),
  show_gene_labels = FALSE,
  map_df_ref = map_df2
)

# Plot 2: 
plot_presence_heatmap2(
  presence_mat = presence_mat,
  genomes_subset = genomes_all,
  genome_order = genome_order_all,
  gene_order = gene_order_global,
  title_text = "Plot 2: Presence/absence (ordered by GROUP -> SPECIES from mapping CSV)",
  out_file = file.path("heatmaps", "plot2_group_species_order_all_genomes.pdf"),
  show_gene_labels = FALSE,
  map_df_ref = map_df2
)


# Group plots: one plot per group

groups_plot_order <- group_levels
for (g in groups_plot_order) {
  acc_g <- map_df2$accession[map_df2$group == g]
  acc_g <- intersect(acc_g, genomes_all)
  if (length(acc_g) < 1) next
  
  genome_order_g <- genome_order_all[genome_order_all %in% acc_g]
  gene_order_g <- get_gene_order(presence_mat, acc_g)
  if (length(gene_order_g) == 0) next
  
  safe_g <- gsub("[^A-Za-z0-9._-]+", "_", g)
  
  plot_presence_heatmap2(
    presence_mat = presence_mat,
    genomes_subset = acc_g,
    genome_order = genome_order_g,
    gene_order = gene_order_g,
    title_text = paste0("Group: ", g, " — genes present in group (N=", length(gene_order_g), ")"),
    out_file = file.path("heatmaps", "groups", paste0("group_", safe_g, "_presence_absence.pdf")),
    show_gene_labels = FALSE,
    map_df_ref = map_df2
  )
}


# Species plots: one plot per species (ordered by group/species order)

species_plot_order <- map_df2 %>%
  arrange(group_rank, species_rank) %>%
  distinct(species) %>%
  pull(species)

for (sp in species_plot_order) {
  acc_sp <- map_df2$accession[map_df2$species == sp]
  acc_sp <- intersect(acc_sp, genomes_all)
  if (length(acc_sp) < 1) next
  
  genome_order_sp <- genome_order_all[genome_order_all %in% acc_sp]
  gene_order_sp <- get_gene_order(presence_mat, acc_sp)
  if (length(gene_order_sp) == 0) next
  
  safe_sp <- gsub("[^A-Za-z0-9._-]+", "_", sp)
  
  plot_presence_heatmap_species(
    presence_mat = presence_mat,
    genomes_subset = acc_sp,
    genome_order = genome_order_sp,
    gene_order = gene_order_sp,
    title_text = paste0("Species: ", sp, " — accessions (N=", length(acc_sp),
                        "); genes present in species (N=", length(gene_order_sp), ")"),
    out_file = file.path("heatmaps", "species", paste0("species_", safe_sp, "_presence_absence_by_accession.pdf")),
    show_gene_labels = FALSE,
    map_df_ref = map_df2
  )
}

message("Done. PDFs written into 'heatmaps/' (groups & species directories).")
