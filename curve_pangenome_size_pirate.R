setwd("/Users/waodedwidaningrat/viridans/pirate_31July2024/results_pirate_rmduplicate_4/summary_pirate_4/analysis_in_R/")

library(data.table)

dt <- fread("PIRATE_gene_families_forR.csv", na.strings = c("", "NA"))

gene_names   <- dt[[1]]
genome_names <- names(dt)[-1]

# Pull matrix of genome columns
m <- as.matrix(dt[, -1, with = FALSE])

# Presence/absence: TRUE if non-empty, FALSE otherwise
present <- !is.na(m) & m != ""
storage.mode(present) <- "logical"

n_genomes <- ncol(present)

seen <- rep(FALSE, nrow(present))
new_genes <- integer(n_genomes)
cum_pan <- integer(n_genomes)

for (i in seq_len(n_genomes)) {
  newly <- present[, i] & !seen
  new_genes[i] <- sum(newly)
  seen <- seen | present[, i]
  cum_pan[i] <- sum(seen)
}


library(ggplot2)

df <- data.frame(
  step = seq_len(n_genomes),
  genome = genome_names,
  new_genes = new_genes,
  cum_pangenome = cum_pan
)

# New genes per added genome
ggplot(df, aes(x = step, y = new_genes)) +
  geom_line() +
  labs(x = "Genomes added", y = "New genes added",
       title = "Pangenome gene discovery curve") +
  theme_minimal()

# cumulative pangenome size
ggplot(df, aes(x = step, y = cum_pangenome)) +
  geom_line() +
  labs(x = "Genomes added", y = "Cumulative pangenome size (genes)",
       title = "Pangenome accumulation") +
  theme_minimal()


pangenome_permute <- function(present, n_perm = 100, seed = 1) {
  set.seed(seed)
  G <- ncol(present)
  R <- nrow(present)
  
  new_mat <- matrix(0L, nrow = n_perm, ncol = G)
  cum_mat <- matrix(0L, nrow = n_perm, ncol = G)
  
  for (p in seq_len(n_perm)) {
    ord <- sample.int(G)
    seen <- rep(FALSE, R)
    
    for (i in seq_len(G)) {
      g <- ord[i]
      newly <- present[, g] & !seen
      new_mat[p, i] <- sum(newly)
      seen <- seen | present[, g]
      cum_mat[p, i] <- sum(seen)
    }
  }
  
  data.frame(
    step = seq_len(G),
    mean_new = colMeans(new_mat),
    sd_new   = apply(new_mat, 2, sd),
    mean_cum = colMeans(cum_mat),
    sd_cum   = apply(cum_mat, 2, sd)
  )
}

avg <- pangenome_permute(present, n_perm = 200, seed = 42)

# Plot mean new genes ± SD
ggplot(avg, aes(x = step, y = mean_new)) +
  geom_line() +
  geom_ribbon(aes(ymin = mean_new - sd_new, ymax = mean_new + sd_new), alpha = 0.2) +
  labs(x = "Genomes added (random order)", y = "New genes added (mean ± SD)",
       title = "Gene discovery curve (permuted)") +
  theme_minimal()

# Plot mean cumulative pangenome ± SD
ggplot(avg, aes(x = step, y = mean_cum)) +
  geom_line() +
  geom_ribbon(aes(ymin = mean_cum - sd_cum, ymax = mean_cum + sd_cum), alpha = 0.2) +
  labs(x = "Genomes added (random order)", y = "Cumulative pangenome size (mean ± SD)",
       title = "Pangenome accumulation (permuted)") +
  theme_minimal()





library(data.table)

meta <- fread("accession_species_group.csv") 

meta2 <- meta[accession %chin% genome_names]

meta2 <- unique(meta2, by = "accession")

# anity checks
cat("Genomes in matrix:", length(genome_names), "\n")
cat("Genomes with metadata:", nrow(meta2), "\n")
cat("Missing metadata for (first 10):\n")
print(head(setdiff(genome_names, meta2$accession), 10))


pangenome_permute <- function(present, n_perm = 200, seed = 42) {
  set.seed(seed)
  G <- ncol(present)
  R <- nrow(present)
  
  new_mat <- matrix(0L, nrow = n_perm, ncol = G)
  cum_mat <- matrix(0L, nrow = n_perm, ncol = G)
  
  for (p in seq_len(n_perm)) {
    ord <- sample.int(G)
    seen <- rep(FALSE, R)
    
    for (i in seq_len(G)) {
      g <- ord[i]
      newly <- present[, g] & !seen
      new_mat[p, i] <- sum(newly)
      seen <- seen | present[, g]
      cum_mat[p, i] <- sum(seen)
    }
  }
  
  data.table(
    step = 1:G,
    mean_new = colMeans(new_mat),
    sd_new   = apply(new_mat, 2, sd),
    mean_cum = colMeans(cum_mat),
    sd_cum   = apply(cum_mat, 2, sd)
  )
}


col_idx <- match(meta2$accession, genome_names)
meta2[, col_idx := col_idx]

# Split metadata by group
meta_by_group <- split(meta2, meta2$vgs)

group_curves <- rbindlist(lapply(names(meta_by_group), function(g) {
  mg <- meta_by_group[[g]]
  mg <- mg[!is.na(col_idx)]  # safety
  
  # subset
  present_g <- present[, mg$col_idx, drop = FALSE]
  
  if (ncol(present_g) < 2) return(NULL)
  
  avg_g <- pangenome_permute(present_g, n_perm = 200, seed = 42)
  avg_g[, group := g]
  avg_g
}), fill = TRUE)


#COLOUR SCHEME:

hex_dt <- fread("group_hex.csv")

# Adjust column names
setnames(hex_dt, tolower(names(hex_dt)))

if (!("group" %in% names(hex_dt))) {
  stop("No 'group' column found in group_hex.csv. Columns are: ",
       paste(names(hex_dt), collapse = ", "))
}

hex_col <- intersect(names(hex_dt), c("hex", "group_hex", "color", "colour"))[1]
if (is.na(hex_col)) {
  stop("No hex/color column found in group_hex.csv. Columns are: ",
       paste(names(hex_dt), collapse = ", "))
}

hex_dt[, hex := get(hex_col)]
hex_dt[, hex := ifelse(startsWith(hex, "#"), hex, paste0("#", hex))]

# Named vector: names = group, values = hex
col_map <- setNames(hex_dt$hex, hex_dt$group)

# Check
missing_cols <- setdiff(unique(group_curves$group), names(col_map))
if (length(missing_cols)) {
  message("These groups are missing colors in group_hex.csv: ",
          paste(missing_cols, collapse = ", "))
}

#PLOTTING
library(ggplot2)

#New genes
ggplot(group_curves, aes(x = step, y = mean_new, color = group, fill = group)) +
  geom_line() +
  geom_ribbon(aes(ymin = mean_new - sd_new, ymax = mean_new + sd_new), alpha = 0.2) +
  facet_wrap(~ group, scales = "free_x") +
  scale_color_manual(values = col_map) +
  scale_fill_manual(values = col_map) +
  labs(x = "Genomes added (within group, random order)",
       y = "New genes added (mean ± SD)",
       title = "Gene discovery curves by group") +
  theme_minimal()

#pangenome accummulation
ggplot(group_curves, aes(x = step, y = mean_cum, color = group, fill = group)) +
  geom_ribbon(aes(ymin = mean_cum - sd_cum, ymax = mean_cum + sd_cum),
              alpha = 0.2, colour = NA) +
  geom_line(linewidth = 1) +
  facet_wrap(~ group, scales = "free_x") +
  scale_color_manual(values = col_map) +
  scale_fill_manual(values = col_map) +
  labs(x = "Genomes added (within group, random order)",
       y = "Cumulative pangenome size (mean ± SD)",
       title = "Pangenome accumulation by group") +
  theme_minimal() +
  theme(legend.position = "none")



plots_new <- lapply(split(group_curves, group_curves$group), function(d) {
  ggplot(d, aes(step, mean_new, color = group, fill = group)) +
    geom_line() +
    scale_color_manual(values = col_map) +
    scale_fill_manual(values = col_map) +
    geom_ribbon(aes(ymin = mean_new - sd_new, ymax = mean_new + sd_new), alpha = 0.2) +
    labs(x = "Genomes added", y = "New genes (mean ± SD)", title = unique(d$group)) +
    theme_minimal()
})

for (p in plots_new) print(p)

pdf("gene_discovery_by_group.pdf", width = 10, height = 7)
for (p in plots_new) print(p)
dev.off()





library(data.table)
library(ggplot2)

# Combined plot
#new genes added
ggplot(group_curves, aes(x = step, y = mean_new, color = group)) +
  geom_line(linewidth = 1) +
  scale_color_manual(values = col_map) +
  labs(x = "Genomes added (within group, random order)",
       y = "New genes added (mean)",
       title = "Gene discovery curves by group") +
  theme_minimal()


ggplot(group_curves, aes(x = step, y = mean_new, color = group, fill = group)) +
  geom_ribbon(aes(ymin = mean_new - sd_new, ymax = mean_new + sd_new),
              alpha = 0.15, colour = NA) +
  geom_line(linewidth = 1) +
  scale_color_manual(values = col_map) +
  scale_fill_manual(values = col_map) +
  labs(x = "Genomes added (within group, random order)",
       y = "New genes added (mean ± SD)",
       title = "Gene discovery curves by group") +
  theme_minimal()

#pangenome accumulation
ggplot(group_curves, aes(x = step, y = mean_cum, color = group)) +
  geom_line(linewidth = 1) +
  scale_color_manual(values = col_map) +
  labs(x = "Genomes added (within group, random order)",
       y = "Cumulative pangenome size (mean)",
       title = "Pangenome accumulation by group") +
  theme_minimal()

ggplot(group_curves, aes(x = step, y = mean_cum, color = group, fill = group)) +
  geom_ribbon(aes(ymin = mean_cum - sd_cum, ymax = mean_cum + sd_cum),
              alpha = 0.15, colour = NA) +
  geom_line(linewidth = 1) +
  scale_color_manual(values = col_map) +
  scale_fill_manual(values = col_map) +
  labs(x = "Genomes added (within group, random order)",
       y = "Cumulative pangenome size (mean ± SD)",
       title = "Pangenome accumulation by group") +
  theme_minimal()

