#Required packages
packages <- c("tidyverse","cluster","factoextra","ggplot2")
install.packages(setdiff(packages, rownames(installed.packages())), repos="https://cloud.r-project.org")
library(tidyverse); library(cluster); library(factoextra); library(ggplot2); library(ggsci); library(dplyr)

#Input files
setwd("/Users/waodedwidaningrat/viridans/fastani_output/")
ani_path <- "ani_matrix_original.csv"      #ANI matrix
ani_distpath <- "ani_distance.csv"      #ANI distance
meta_path <- "accession_speciescoregenome_forANI.csv"       #accession, species

#Parameters
k_chosen <- 29          # number of clusters (≈ number of species)
mds_dims <- 3           # number of MDS dimensions
random_reps <- 4        # random reps per species
output_prefix <- "ani_mds_kmeans"

#Load ANI distance matrix
anidist <- read.csv(ani_distpath, row.names=1, check.names=FALSE)
anidist <- as.matrix(anidist)

#detect scale and convert to distance if needed (this is if used the matrix ani csv, if use distance matrix no need to do this)
#rng <- range(ani, na.rm=TRUE)
#if (rng[2] > 1) ani <- ani/100
#dist_mat <- 1 - ani
#dobj <- as.dist(dist_mat)

#if use distance matrix no need to do this
dist_mat <- anidist

#Replace NA with maximum observed distance  and run MDS
#Assumes: dist_mat is a square numeric matrix with rownames/colnames (may contain NA)
#mds_dims is set (e.g. 2 or 3)
# sanity checks
if (!is.matrix(dist_mat) || nrow(dist_mat) != ncol(dist_mat)) stop("dist_mat must be a square matrix")
if (is.null(rownames(dist_mat)) || is.null(colnames(dist_mat))) stop("dist_mat must have row and column names")

#Ensure diagonal is zero (distance to self)
diag(dist_mat) <- 0

#Symmetrize: where one side exists and the other is NA, copy across
#where both present but slightly different, average them
for (i in seq_len(nrow(dist_mat))) {
  for (j in seq_len(ncol(dist_mat))) {
    if (i == j) next
    a <- dist_mat[i, j]
    b <- dist_mat[j, i]
    if (is.na(a) && !is.na(b)) {
      dist_mat[i, j] <- b
    } else if (!is.na(a) && is.na(b)) {
      dist_mat[j, i] <- a
    } else if (!is.na(a) && !is.na(b) && abs(a - b) > 1e-8) {
      # small asymmetry: average them to be safe
      avgv <- mean(c(a, b), na.rm = TRUE)
      dist_mat[i, j] <- dist_mat[j, i] <- avgv
    }
  }
}

#Compute maximum observed distance (excluding NA and diagonal zeros if any)
max_d <- max(dist_mat[upper.tri(dist_mat)], na.rm = TRUE)
if (!is.finite(max_d)) stop("No finite distances found to compute max distance — check your dist_mat")

#Replace remaining NAs with max_d
na_idx <- which(is.na(dist_mat), arr.ind = TRUE)
if (nrow(na_idx) > 0) {
  for (r in seq_len(nrow(na_idx))) {
    i <- na_idx[r, 1]; j <- na_idx[r, 2]
    dist_mat[i, j] <- max_d
    dist_mat[j, i] <- max_d
  }
  message("Replaced ", nrow(na_idx)/2, " NA pairs with max distance = ", format(max_d, digits = 6))
} else {
  message("No NA entries found (nothing to impute).")
}

#sanity: diagonal zero, symmetric
diag(dist_mat) <- 0
if (!all.equal(dist_mat, t(dist_mat), tolerance = 1e-8)) warning("Matrix not perfectly symmetric after imputation — check")

#Convert to dist and run classical MDS
dobj <- as.dist(dist_mat)
mds <- cmdscale(dobj, k = mds_dims)
mds_df <- as.data.frame(mds)
colnames(mds_df) <- paste0("MDS", seq_len(ncol(mds_df)))
mds_df$accession <- rownames(dist_mat)

message("MDS completed: ", nrow(mds_df), " points with ", ncol(mds_df), " dimensions.")
#mds_df now contains MDS1, MDS2, ... and accession
write.csv(dist_mat, "dist_matrix_imputed_maxna.csv", row.names = TRUE)

#MDS
#mds <- cmdscale(dobj, k=mds_dims)
#mds_df <- as.data.frame(mds)
#colnames(mds_df) <- paste0("MDS",1:mds_dims)
#mds_df$accession <- rownames(ani)

#K-means
set.seed(42)
km <- kmeans(mds_df[,1:mds_dims], centers=k_chosen, nstart=25)
mds_df$cluster <- factor(km$cluster)

#Choose medoids per cluster
dist_full <- as.matrix(dobj)
genomes <- rownames(anidist)
medoids <- sapply(levels(mds_df$cluster), function(cl){
  members <- mds_df$accession[mds_df$cluster==cl]
  sub <- dist_full[members,members,drop=FALSE]
  members[which.min(rowSums(sub))]
})
medoid_df <- data.frame(cluster=levels(mds_df$cluster), medoid=medoids)
write.csv(medoid_df, paste0(output_prefix,"_medoids.csv"), row.names=FALSE)

#add random per species
meta <- read.csv(meta_path, stringsAsFactors=FALSE)
if (!"accession" %in% colnames(meta) || !"species" %in% colnames(meta)) {
  stop("metadata must contain columns 'accession' and 'species'.")
}

sel <- medoid_df$medoid
set.seed(1)
for(sp in unique(meta$species)){
  ids <- meta$accession[meta$species==sp]
  add <- setdiff(ids, sel)
  sel <- c(sel, sample(add, min(random_reps, length(add))))
}
write.csv(data.frame(accession=sel), paste0(output_prefix,"_selected.csv"), row.names=FALSE)

#Within/between species ANI stats
ani <- read.csv(ani_path, row.names = 1, check.names = FALSE, stringsAsFactors = FALSE)
ani <- as.matrix(ani)   # numeric matrix with NAs
# Confirm row/col names match
if (!all(rownames(ani) == colnames(ani))) {
  warning("Row and column names differ; using rownames for columns where necessary.")
  colnames(ani) <- rownames(ani)
}

#compute within/between summaries
species <- meta$species; names(species) <- meta$accession
species_list <- unique(species)
within_summary <- list()
between_summary <- list()

for (sp in species_list) {
  ids <- names(species)[species == sp]
  ids <- intersect(ids, colnames(ani))  # only those present in ANI matrix
  others <- setdiff(colnames(ani), ids)
  
  #Within-species pairwise values
  if (length(ids) >= 2) {
    sub_w <- ani[ids, ids, drop = FALSE]
    vals_w <- sub_w[lower.tri(sub_w)]
    n_present <- sum(!is.na(vals_w))
    n_possible <- length(vals_w)
    mean_w <- ifelse(n_present>0, mean(vals_w, na.rm=TRUE), NA_real_)
    sd_w   <- ifelse(n_present>0, sd(vals_w, na.rm=TRUE), NA_real_)
    min_w  <- ifelse(n_present>0, min(vals_w, na.rm=TRUE), NA_real_)
    max_w  <- ifelse(n_present>0, max(vals_w, na.rm=TRUE), NA_real_)
  } else {
    n_present <- 0; n_possible <- 0
    mean_w <- sd_w <- min_w <- max_w <- NA_real_
  }
  within_summary[[sp]] <- tibble(species = sp, n_genomes = length(ids),
                                 n_within_pairs_present = n_present,
                                 n_within_pairs_possible = n_possible,
                                 mean_within_pct = mean_w,
                                 sd_within_pct = sd_w,
                                 min_within_pct = min_w,
                                 max_within_pct = max_w)
  
  #Between-species (members vs all others)
  if (length(ids) >= 1 && length(others) >= 1) {
    sub_b <- ani[ids, others, drop = FALSE]
    vals_b <- as.vector(sub_b)
    n_present_b <- sum(!is.na(vals_b))
    n_possible_b <- length(vals_b)
    mean_b <- ifelse(n_present_b>0, mean(vals_b, na.rm=TRUE), NA_real_)
    sd_b   <- ifelse(n_present_b>0, sd(vals_b, na.rm=TRUE), NA_real_)
    min_b  <- ifelse(n_present_b>0, min(vals_b, na.rm=TRUE), NA_real_)
    max_b  <- ifelse(n_present_b>0, max(vals_b, na.rm=TRUE), NA_real_)
  } else {
    n_present_b <- 0; n_possible_b <- 0
    mean_b <- sd_b <- min_b <- max_b <- NA_real_
  }
  between_summary[[sp]] <- tibble(species = sp, n_genomes = length(ids),
                                  n_between_pairs_present = n_present_b,
                                  n_between_pairs_possible = n_possible_b,
                                  mean_between_pct = mean_b,
                                  sd_between_pct = sd_b,
                                  min_between_pct = min_b,
                                  max_between_pct = max_b)
}

within_df <- bind_rows(within_summary)
between_df <- bind_rows(between_summary)
summary_tbl <- left_join(within_df, between_df, by = c("species","n_genomes"))

# Save summary
write.csv(summary_tbl, paste0(output_prefix, "_within_between_species.csv"), row.names = FALSE)
message("Wrote species summary: ", paste0(output_prefix, "_within_between_species.csv"))



######PLOTTING STARTING HERE
#PLOT MDS BY CLUSTER
library(tidyverse); library(cluster); library(factoextra); library(ggplot2); library(ggsci); library(dplyr)
library(showtext)

# common font folders on macOS
list.files("/Library/Fonts", pattern = "Times", ignore.case = TRUE, full.names = TRUE)
list.files("/System/Library/Fonts", pattern = "Times", ignore.case = TRUE, full.names = TRUE)
# user fonts (if installed per-user)
list.files("~/Library/Fonts", pattern = "Times", ignore.case = TRUE, full.names = TRUE)

base_font <- "Times New Roman"


showtext_auto()

ggplot(mds_df, aes(MDS1, MDS2, color=cluster)) +
  geom_point(size=0.8, alpha=0.8) +
  theme_minimal(base_size = 12) +
  scale_color_igv() +  # distinct palette for many species
  guides(color = guide_legend(override.aes = list(size = 4))) +
  labs(title = "K-means Clustering of MDS - ANI distances - colored by cluster") +
  theme(
    legend.position = "right",
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10, face="italic")
  )

mds_df <- merge(mds_df, meta[, c("accession","species")], by="accession", all.x=TRUE)

#Load species colours
colour_path <- "species_colours_hex_final.csv"
col_df <- read.csv(colour_path, stringsAsFactors = FALSE)

#Ensure all species in MDS data have a colour assigned
missing_col <- setdiff(unique(mds_df$species), col_df$species)
if (length(missing_col) > 0) {
  warning("These species have no assigned colour: ", paste(missing_col, collapse=", "))
}

#Build a named colour vector for ggplot
col_vec <- setNames(col_df$colour, col_df$species)

#MDS for each species using automatic colours for species
ggplot(mds_df, aes(MDS1, MDS2, color=species)) +
  geom_point(size=0.8, alpha=0.8) +
  theme_minimal(base_size = 12) +
  scale_color_igv() +
  guides(color = guide_legend(override.aes = list(size = 4))) +  # make legend dots larger
  labs(title = "K-means Clustering of MDS - ANI distances") +
  theme(
    legend.position = "right",
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10, face="italic")
  )

# Compute convex hulls for each species 
hulls <- mds_df %>%
  group_by(species) %>%
  slice(chull(MDS1, MDS2))

# mds with convex hulls using automatic colours for species
ggplot(mds_df, aes(MDS1, MDS2, color=species, fill=species)) +
  geom_polygon(data = hulls, aes(group=species), alpha = 0.2, color="black", linetype="dashed") +
  geom_point(size=0.8, alpha=0.8) +
  theme_minimal(base_size = 12) +
  scale_color_igv() +
  scale_fill_igv() +
  guides(color = guide_legend(override.aes = list(size = 4))) +
  labs(title="K-means Clustering of MDS - ANI distances") +
  theme(
    legend.position="right",
    legend.title = element_text(size=12),
    legend.text = element_text(size=10, face="italic")
  )

#Plot using custom colours from CSV
ggplot(mds_df, aes(MDS1, MDS2, color=species, fill=species)) +
  geom_point(size=0.9, alpha=0.8) +
  theme_minimal(base_size=12) +
  scale_color_manual(values=col_vec) +
  scale_fill_manual(values=col_vec) +
  guides(color = guide_legend(override.aes = list(size=5))) +  # larger legend dots
  labs(title="K-means Clustering of MDS - ANI distances") +
  theme(
    legend.position="right",
    legend.title = element_text(size=12),
    legend.text = element_text(size=10, face="italic")
  )

# Compute convex hulls for cluster outlines (optional)
hulls <- mds_df %>%
  group_by(species) %>%
  slice(chull(MDS1, MDS2))

#linetypes:
#solid	"solid" or 1
#dashed	"dashed" or 2
#dotted	"dotted" or 3
#dotdash	"dotdash" or 4
#longdash	"longdash" or 5
#twodash	"twodash" or 6

# --- Plot using custom colours from CSV ---
ggplot(mds_df, aes(MDS1, MDS2, color=species, fill=species)) +
  geom_polygon(data=hulls, aes(group=species), alpha=0.2,
               color="black", linetype="solid") +
  geom_point(size=0.9, alpha=0.8) +
  theme_minimal(base_size=12) +
  scale_color_manual(values=col_vec) +
  scale_fill_manual(values=col_vec) +
  guides(color = guide_legend(override.aes = list(size=5))) +  # larger legend dots
  labs(title="K-means Clustering of MDS - ANI distances") +
  theme(
    legend.position="right",
    legend.title = element_text(size=12),
    legend.text = element_text(size=10, face="italic")
  )

ggplot(mds_df, aes(MDS1, MDS2, color=species, fill=species)) +
  geom_polygon(data=hulls, aes(group=species), alpha=0.2, linetype="solid") +
  geom_point(size=0.9, alpha=0.8) +
  theme_minimal(base_size=12) +
  scale_color_manual(values=col_vec) +
  scale_fill_manual(values=col_vec) +
  guides(color = guide_legend(override.aes = list(size=5))) +  # larger legend dots
  labs(title="K-means Clustering of MDS - ANI distances") +
  theme(
    legend.position="right",
    legend.title = element_text(size=12),
    legend.text = element_text(size=10, face="italic")
  )

ggplot(mds_df, aes(MDS1, MDS2, color=species, fill=species)) +
  geom_polygon(data=hulls, aes(group=species), alpha=0.2, linetype="solid") +
  geom_point(size=0.9, alpha=0.8) +
  theme_minimal(base_size=12) +
  scale_color_manual(values=col_vec) +
  scale_fill_manual(values=col_vec) +
  guides(color = guide_legend(override.aes = list(size=5))) +  # larger legend dots
  labs(title="K-means Clustering of MDS - ANI distances") +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 14,colour = "black", face = "bold"),
    axis.text.y = element_text(size = 14, colour = "black", face = "bold"),
    plot.title  = element_text(size = 14, face = "bold"),
    axis.title.x.bottom = element_text(size = 14, face = "bold"),
    axis.title.y.left = element_text(size = 14, face = "bold"),
    legend.position="right",
    legend.title = element_text(size=12),
    legend.text = element_text(size=12, face="italic")
  )


ggplot(mds_df, aes(MDS1, MDS2, color=species)) +
  stat_ellipse(aes(group=species), level=0.95, alpha=0.3, linetype="solid") +
  geom_point(size=0.8, alpha=0.8) +
  theme_minimal(base_size = 12) +
  scale_color_manual(values=col_vec) +
  scale_fill_manual(values=col_vec) +
  guides(color = guide_legend(override.aes = list(size = 4))) +
  labs(title="K-means Clustering of MDS - ANI distances with 95% ellipses") +
  theme(
    legend.position="right",
    legend.title = element_text(size=12),
    legend.text = element_text(size=10, face="italic")
  )





# Full pipeline:
#  - read ANI matrix (0..100, NAs allowed) and metadata (accession,species)
#  - compute within-/between-species ANI summaries
#  - choose medoid per species (minimize mean distance to other members)
#  - compute medoid-to-medoid distance and ANI matrices
#  - classical MDS (imputing missing distances)
#  - convex-hull overlap (intersection/union proportion) in MDS 2D
#  - species×species mean ANI matrix + counts + long table
#  - save CSVs and a single MDS hull plot


# User parameters
ani_csv <- "ani_matrix_original.csv"
meta_csv <- "accession_speciescoregenome_forANI.csv"
output_prefix <- "species_dist"
impute_method <- "max"   # "max", "mean", or numeric
mds_dims <- 2

# Packages
pkgs <- c("tidyverse", "sf", "ggplot2")
install.packages(setdiff(pkgs, rownames(installed.packages())), repos = "https://cloud.r-project.org")
library(tidyverse)
library(sf)
library(ggplot2)

# Load inputs
message("Loading inputs...")
ani <- as.matrix(read.csv(ani_csv, row.names = 1, check.names = FALSE, stringsAsFactors = FALSE))
meta <- read.csv(meta_csv, stringsAsFactors = FALSE)

# Basic checks
if (!("accession" %in% colnames(meta)) || !("species" %in% colnames(meta))) {
  stop("metadata must have columns 'accession' and 'species'")
}
if (!all(rownames(ani) == colnames(ani))) {
  warning("Row and column names of ANI matrix differ; forcing columns to match rownames.")
  colnames(ani) <- rownames(ani)
}

# Restrict to common IDs between ANI matrix and metadata
common_ids <- intersect(rownames(ani), meta$accession)
if (length(common_ids) == 0) stop("No overlap between ANI matrix IDs and metadata$accession")
if (length(common_ids) < nrow(ani)) {
  message(sprintf("Restricting to %d genomes present in both ANI matrix and metadata (out of %d).",
                  length(common_ids), nrow(ani)))
  ani <- ani[common_ids, common_ids, drop = FALSE]
}
meta <- meta %>% filter(accession %in% common_ids)

# Convert ANI to similarity (0..1) and distance (0..1)
sim <- ani / 100   # similarity 0..1; NAs allowed
dist_mat <- 1 - sim  # distance 0..1 (small = similar)

# Within- and between-species stats
message("Computing within- and between-species stats...")
species_vec <- meta$species; names(species_vec) <- meta$accession
species_list <- unique(species_vec)

# Check species with <3 genomes
species_counts <- table(species_vec) %>% as.data.frame()
colnames(species_counts) <- c("species", "n_genomes")

# identify species with small number of genomes
small_species <- species_counts %>% filter(n_genomes < 3)

# print summary
message("Species genome counts:")
print(species_counts)

if (nrow(small_species) > 0) {
  message("\n⚠️  The following species have fewer than 3 genomes:")
  print(small_species)
} else {
  message("\n✅ All species have 3 or more genomes.")
}

# save outputs
write.csv(species_counts, paste0(output_prefix, "_species_counts.csv"), row.names = FALSE)
message("Wrote species genome counts to: ", paste0(output_prefix, "_species_counts.csv"))

within_rows <- list()
between_rows <- list()
for (sp in species_list) {
  ids <- names(species_vec)[species_vec == sp] %>% intersect(colnames(sim))
  others <- setdiff(colnames(sim), ids)
  
  # within
  if (length(ids) >= 2) {
    sub_sim <- sim[ids, ids, drop = FALSE]
    vals_within <- sub_sim[lower.tri(sub_sim)]
    n_present <- sum(!is.na(vals_within))
    n_possible <- length(vals_within)
    within_rows[[sp]] <- tibble(species = sp,
                                n_genomes = length(ids),
                                n_within_present = n_present,
                                n_within_possible = n_possible,
                                mean_within_pct = ifelse(n_present>0, mean(vals_within, na.rm = TRUE)*100, NA_real_),
                                sd_within_pct = ifelse(n_present>0, sd(vals_within, na.rm = TRUE)*100, NA_real_),
                                min_within_pct = ifelse(n_present>0, min(vals_within, na.rm = TRUE)*100, NA_real_),
                                max_within_pct = ifelse(n_present>0, max(vals_within, na.rm = TRUE)*100, NA_real_))
  } else {
    within_rows[[sp]] <- tibble(species = sp,
                                n_genomes = length(ids),
                                n_within_present = 0,
                                n_within_possible = 0,
                                mean_within_pct = NA_real_,
                                sd_within_pct = NA_real_,
                                min_within_pct = NA_real_,
                                max_within_pct = NA_real_)
  }
  
  # between
  if (length(ids) >= 1 && length(others) >= 1) {
    sub_sim_b <- sim[ids, others, drop = FALSE]
    vals_between <- as.vector(sub_sim_b)
    n_present_b <- sum(!is.na(vals_between))
    n_possible_b <- length(vals_between)
    between_rows[[sp]] <- tibble(species = sp,
                                 n_genomes = length(ids),
                                 n_between_present = n_present_b,
                                 n_between_possible = n_possible_b,
                                 mean_between_pct = ifelse(n_present_b>0, mean(vals_between, na.rm=TRUE)*100, NA_real_),
                                 sd_between_pct = ifelse(n_present_b>0, sd(vals_between, na.rm=TRUE)*100, NA_real_),
                                 min_between_pct = ifelse(n_present_b>0, min(vals_between, na.rm=TRUE)*100, NA_real_),
                                 max_between_pct = ifelse(n_present_b>0, max(vals_between, na.rm=TRUE)*100, NA_real_))
  } else {
    between_rows[[sp]] <- tibble(species = sp,
                                 n_genomes = length(ids),
                                 n_between_present = 0,
                                 n_between_possible = 0,
                                 mean_between_pct = NA_real_,
                                 sd_between_pct = NA_real_,
                                 min_between_pct = NA_real_,
                                 max_between_pct = NA_real_)
  }
}

within_df <- bind_rows(within_rows)
between_df <- bind_rows(between_rows)
species_summary <- left_join(within_df, between_df, by = c("species","n_genomes"))
write.csv(species_summary, paste0(output_prefix, "_species_within_between.csv"), row.names = FALSE)
message("Wrote: ", paste0(output_prefix, "_species_within_between.csv"))

# Choose medoid per species (minimize mean distance to other members)
message("Selecting medoids per species...")
medoid_rows <- list()
for (sp in species_list) {
  members <- names(species_vec)[species_vec == sp] %>% intersect(colnames(dist_mat))
  if (length(members) == 0) next
  if (length(members) == 1) {
    medoid_rows[[sp]] <- tibble(species = sp, medoid = members, medoid_mean_distance = 0, n_members = 1, fraction_missing = 0)
    next
  }
  sub_dist <- dist_mat[members, members, drop = FALSE]
  diag(sub_dist) <- 0
  mean_d <- sapply(seq_len(nrow(sub_dist)), function(i) {
    vals <- sub_dist[i, -i]
    if (all(is.na(vals))) return(Inf)
    mean(vals, na.rm = TRUE)
  })
  best <- which.min(mean_d)
  medoid_id <- members[best]
  total_pairs <- length(members)*(length(members)-1)/2
  present_pairs <- sum(!is.na(sub_dist[lower.tri(sub_dist)]))
  fraction_missing <- ifelse(total_pairs>0, 1 - (present_pairs / total_pairs), NA_real_)
  medoid_rows[[sp]] <- tibble(species = sp, medoid = medoid_id, medoid_mean_distance = mean_d[best],
                              n_members = length(members), fraction_missing = fraction_missing)
}
medoid_df <- bind_rows(medoid_rows)
write.csv(medoid_df, paste0(output_prefix, "_species_medoids.csv"), row.names = FALSE)
message("Wrote: ", paste0(output_prefix, "_species_medoids.csv"))

# Inter-species (medoid-to-medoid) matrices
message("Computing medoid-to-medoid matrices...")
medoids <- medoid_df$medoid
medoid_dist <- dist_mat[medoids, medoids, drop = FALSE]
write.csv(medoid_dist, paste0(output_prefix, "_medoid_to_medoid_distance.csv"), row.names = TRUE)

medoid_ani <- (1 - medoid_dist) * 100
write.csv(medoid_ani, paste0(output_prefix, "_medoid_to_medoid_ani.csv"), row.names = TRUE)
message("Wrote medoid-to-medoid distance and ANI matrices.")

# MDS for visualization (impute missing distances first)
message("Preparing MDS (imputing missing distances)...")
d_obs <- dist_mat
na_count <- sum(is.na(d_obs[lower.tri(d_obs)]))
message("Missing pairwise distances (lower-tri): ", na_count)

if (identical(impute_method, "max")) {
  fill_val <- max(d_obs, na.rm = TRUE)
} else if (identical(impute_method, "mean")) {
  fill_val <- mean(d_obs, na.rm = TRUE)
} else if (is.numeric(impute_method)) {
  fill_val <- as.numeric(impute_method)
} else {
  stop("impute_method must be 'max', 'mean', or numeric value")
}
message("Imputing missing distances with: ", round(fill_val, 6))
d_imputed <- d_obs
d_imputed[is.na(d_imputed)] <- fill_val

# classical MDS
dobj <- as.dist(d_imputed)
mds_res <- cmdscale(dobj, k = mds_dims, eig = TRUE)
mds_coords <- as.data.frame(mds_res$points)
colnames(mds_coords) <- paste0("MDS", seq_len(ncol(mds_coords)))
# attach genome id column (accession)
mds_coords$accession <- rownames(d_imputed)
write.csv(mds_coords, paste0(output_prefix, "_mds_coords.csv"), row.names = FALSE)
message("Wrote: ", paste0(output_prefix, "_mds_coords.csv"))

# Convex hull overlap (intersection / union proportion)
message("Building convex hulls and computing overlap proportions...")
# Helper to get coords DF robustly
get_coords_df <- function(mds_df) {
  mds_df <- as.data.frame(mds_df)
  if ("accession" %in% colnames(mds_df)) {
    coords <- mds_df
    rownames(coords) <- coords$accession
    coords$accession <- NULL
  } else if ("genome_id" %in% colnames(mds_df)) {
    coords <- mds_df
    rownames(coords) <- coords$genome_id
    coords$genome_id <- NULL
  } else if (!is.null(rownames(mds_df)) && any(rownames(mds_df) != "")) {
    coords <- mds_df
    if (!("accession" %in% colnames(coords))) coords$accession <- rownames(coords)
  } else {
    stop("mds_coords has no accession/genome_id column and no meaningful rownames.")
  }
  if (!all(c("MDS1","MDS2") %in% colnames(coords))) stop("mds_coords must contain columns 'MDS1' and 'MDS2'.")
  return(coords)
}

coords <- get_coords_df(mds_coords)

# Build hull polygons per species
hull_polys <- list()
for (sp in species_list) {
  ids <- names(species_vec)[species_vec == sp] %>% intersect(rownames(coords))
  if (length(ids) < 3) { hull_polys[[sp]] <- NULL; next }
  pts <- coords[ids, c("MDS1","MDS2"), drop = FALSE]
  pts <- pts[complete.cases(pts), , drop = FALSE]
  if (nrow(pts) < 3) { hull_polys[[sp]] <- NULL; next }
  unique_pts <- unique(pts)
  if (nrow(unique_pts) < 3) { hull_polys[[sp]] <- NULL; next }
  ch_idx <- chull(unique_pts$MDS1, unique_pts$MDS2)
  ch_coords <- unique_pts[c(ch_idx, ch_idx[1]), c("MDS1","MDS2")]
  poly_sfc <- tryCatch({
    s <- st_sfc(st_polygon(list(as.matrix(ch_coords))), crs = NA)
    s <- st_make_valid(s)
    s
  }, error = function(e) NULL)
  hull_polys[[sp]] <- poly_sfc
}

# Pairwise overlap matrix
nsp <- length(species_list)
overlap_mat <- matrix(NA_real_, nrow = nsp, ncol = nsp, dimnames = list(species_list, species_list))
for (i in seq_along(species_list)) {
  for (j in seq_along(species_list)) {
    sp1 <- species_list[i]; sp2 <- species_list[j]
    p1 <- hull_polys[[sp1]]; p2 <- hull_polys[[sp2]]
    if (is.null(p1) || is.null(p2)) { overlap_mat[i,j] <- NA_real_; next }
    inter <- tryCatch(st_intersection(p1, p2), error = function(e) NULL)
    union_ <- tryCatch(st_union(p1, p2), error = function(e) NULL)
    inter_area <- 0; union_area <- 0
    if (!is.null(inter) && length(inter) > 0) {
      inter_area <- tryCatch(sum(as.numeric(st_area(inter))), error = function(e) 0)
    }
    if (!is.null(union_) && length(union_) > 0) {
      union_area <- tryCatch(sum(as.numeric(st_area(union_))), error = function(e) 0)
    }
    if (is.na(union_area) || union_area <= 0) {
      overlap_mat[i,j] <- NA_real_
    } else {
      overlap_mat[i,j] <- as.numeric(inter_area / union_area)
    }
  }
}
write.csv(overlap_mat, paste0(output_prefix, "_hull_overlap_matrix.csv"), row.names = TRUE)
message("Wrote: ", paste0(output_prefix, "_hull_overlap_matrix.csv"))

# Prepare hull_plot_df for plotting
hull_plot_df <- bind_rows(lapply(names(hull_polys), function(sp) {
  poly <- hull_polys[[sp]]
  if (is.null(poly)) return(NULL)
  coords_mat <- as.data.frame(st_coordinates(poly)[,1:2])
  names(coords_mat) <- c("X","Y")
  coords_mat$species <- sp
  coords_mat$ord <- seq_len(nrow(coords_mat))
  coords_mat
}), .id = NULL)

# MDS with hulls and medoids
message("Preparing final plot...")
plot_df <- coords
if (!("accession" %in% colnames(plot_df))) plot_df$accession <- rownames(plot_df)
plot_df$species <- species_vec[plot_df$accession]

# medoid positions (join to MDS coords)
medoid_pos <- medoid_df %>% left_join(plot_df %>% dplyr::select(accession, MDS1, MDS2),
                                      by = c("medoid" = "accession"))

# Build plot
p <- ggplot(plot_df, aes(x = MDS1, y = MDS2, color = species)) +
  geom_point(alpha = 0.6, size = 0.9) +
  theme_minimal(base_size=12, base_family = base_font) +
  labs(title = "MDS with species convex hulls and medoids")

p <- ggplot(plot_df, aes(MDS1, MDS2, color=species, fill=species)) +
  geom_polygon(data=hulls, aes(group=species), alpha=0.2, linetype="solid") +
  geom_point(size=0.9, alpha=0.8) +
  theme_minimal(base_size=12, base_family = base_font) +
  scale_color_manual(values=col_vec) +
  scale_fill_manual(values=col_vec) +
  guides(color = guide_legend(override.aes = list(size=5))) +  # larger legend dots
  labs(title="MDS with species convex hulls and medoids") +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 14,colour = "black", face = "bold"),
    axis.text.y = element_text(size = 14, colour = "black", face = "bold"),
    plot.title  = element_text(size = 14, face = "bold"),
    axis.title.x.bottom = element_text(size = 14, face = "bold"),
    axis.title.y.left = element_text(size = 14, face = "bold"),
    legend.position="right",
    legend.title = element_text(size=12),
    legend.text = element_text(size=10, face="italic")
  )

p

if (exists("hull_plot_df") && nrow(hull_plot_df) > 0) {
  p <- p + geom_polygon(data = hull_plot_df, aes(x = X, y = Y, group = species, fill = species),
                        color = "black", alpha = 0.2, show.legend = FALSE)
}
if (nrow(medoid_pos) > 0) {
  medoid_pos_valid <- medoid_pos %>% filter(!is.na(MDS1) & !is.na(MDS2))
  if (nrow(medoid_pos_valid) > 0) {
    p <- p + geom_point(data = medoid_pos_valid, aes(x = MDS1, y = MDS2),
                        color = "black", shape = 8, size = 2)
  }
}

p

ggsave(paste0(output_prefix, "_mds_hulls.pdf"), p, width = 10, height = 8, dpi = 200)
message("Wrote: ", paste0(output_prefix, "_mds_hulls.pdf"))

# Species × species mean ANI matrices and long table 
message("Computing species × species mean ANI and counts...")
species_list_sorted <- sort(unique(species_vec))
S <- length(species_list_sorted)
mean_mat <- matrix(NA_real_, nrow = S, ncol = S, dimnames = list(species_list_sorted, species_list_sorted))
sd_mat   <- matrix(NA_real_, nrow = S, ncol = S, dimnames = list(species_list_sorted, species_list_sorted))
n_present_mat <- matrix(0L, nrow = S, ncol = S, dimnames = list(species_list_sorted, species_list_sorted))
n_possible_mat <- matrix(0L, nrow = S, ncol = S, dimnames = list(species_list_sorted, species_list_sorted))

genomes_by_species <- lapply(species_list_sorted, function(sp) {
  names(species_vec)[species_vec == sp]
})
names(genomes_by_species) <- species_list_sorted

for (i in seq_along(species_list_sorted)) {
  for (j in seq_along(species_list_sorted)) {
    sp1 <- species_list_sorted[i]; sp2 <- species_list_sorted[j]
    ids1 <- intersect(genomes_by_species[[sp1]], colnames(sim))
    ids2 <- intersect(genomes_by_species[[sp2]], colnames(sim))
    if (length(ids1) == 0 || length(ids2) == 0) {
      mean_mat[i,j] <- NA_real_; sd_mat[i,j] <- NA_real_
      n_present_mat[i,j] <- 0L; n_possible_mat[i,j] <- 0L
      next
    }
    sub <- sim[ids1, ids2, drop = FALSE]
    if (sp1 == sp2) {
      if (nrow(sub) >= 2) {
        vals <- sub[lower.tri(sub)]
      } else {
        vals <- numeric(0)
      }
      n_possible <- length(sub[lower.tri(sub)])
    } else {
      vals <- as.vector(sub)
      n_possible <- length(vals)
    }
    n_present <- sum(!is.na(vals))
    mean_val <- if (n_present > 0) mean(vals, na.rm = TRUE) else NA_real_
    sd_val <- if (n_present > 0) sd(vals, na.rm = TRUE) else NA_real_
    mean_mat[i,j] <- mean_val * 100
    sd_mat[i,j] <- sd_val * 100
    n_present_mat[i,j] <- n_present
    n_possible_mat[i,j] <- n_possible
  }
}

write.csv(as.data.frame(mean_mat), paste0(output_prefix, "_species_pairwise_mean_ani_matrix.csv"), row.names = TRUE)
write.csv(as.data.frame(n_present_mat), paste0(output_prefix, "_species_pairwise_counts_matrix.csv"), row.names = TRUE)

# Long form
long_rows <- list()
for (i in seq_along(species_list_sorted)) {
  for (j in seq_along(species_list_sorted)) {
    long_rows[[length(long_rows) + 1]] <- data.frame(
      species1 = species_list_sorted[i],
      species2 = species_list_sorted[j],
      mean_ani_pct = mean_mat[i,j],
      sd_ani_pct = sd_mat[i,j],
      n_present = n_present_mat[i,j],
      n_possible = n_possible_mat[i,j],
      stringsAsFactors = FALSE
    )
  }
}
long_df <- bind_rows(long_rows)
write.csv(long_df, paste0(output_prefix, "_species_pairwise_long.csv"), row.names = FALSE)
message("Wrote species pairwise ANI outputs with prefix: ", output_prefix)

message("All done.")


# Heatmap of species × species mean ANI

library(reshape2)
library(viridis)
if (!"viridis" %in% rownames(installed.packages())) install.packages("viridis")

mean_ani_path <- paste0(output_prefix, "_species_pairwise_mean_ani_matrix.csv")
mean_ani_mat <- read.csv(mean_ani_path, row.names = 1, check.names = FALSE, stringsAsFactors = FALSE)
mat <- as.matrix(mean_ani_mat)  # percent, may contain NAs

# Build a distance matrix for clustering:
# distance = 1 - mean_ani/100 ; set NA entries to 1 (max distance) so clustering can proceed.
dist_species <- 1 - (mat / 100)
# Replace any NaN / NA / Inf with 1 (max distance)
dist_species[!is.finite(dist_species)] <- NA_real_
# If entire row or column is NA, leave NAs and then set to 1
dist_species[is.na(dist_species)] <- 1.0

# Ensure diagonal is zero
diag(dist_species) <- 0

# hierarchical clustering
species_order <- rownames(dist_species)  # default order if clustering fails
try({
  hc <- hclust(as.dist(dist_species), method = "average")
  species_order <- hc$labels[hc$order]
}, silent = TRUE)

# Melt for ggplot and apply ordering
# Melt matrix to long form 
# mat must have row/col names that are the species labels
df_long <- reshape2::melt(mat, varnames = c("species1","species2"), value.name = "mean_ani")

# Read plotting order and apply it to both axes
ord <- readr::read_csv("species_order.csv", show_col_types = FALSE)
species_order <- ord$species

# Keep listed species in given order, then append any others that appear in the data
present <- unique(c(as.character(df_long$species1), as.character(df_long$species2)))
levels_order <- c(intersect(species_order, present), setdiff(present, species_order))

df_long <- df_long %>%
  mutate(
    species1 = factor(as.character(species1), levels = levels_order),
    species2 = factor(as.character(species2), levels = levels_order)
  )

# Colour scheme

pal_cols  <- colorRampPalette(c("#f7fbff", "#c6dbef", "#6baed6", "#2171b5", "#08306b"))(256)

# Plot heatmap
p_heat <- ggplot(df_long, aes(x = species2, y = species1, fill = mean_ani)) +
  geom_tile(color = "white", linewidth = 0.2) +
  scale_fill_viridis(name = "%ANI", na.value = "grey90", option = "C") +
  coord_fixed() +
  theme_minimal(base_size = 10, base_family = base_font) +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10,colour = "black", face = "italic"),
    axis.text.y = element_text(size = 10, colour = "black", face = "italic"),
    plot.title  = element_text(size = 10, face = "bold"),
    legend.title = element_text(size = 10, face = "bold"),
    axis.title.x.bottom = element_text(size = 10, face = "bold"),
    axis.title.y.left = element_text(size = 10, face = "bold"),
    legend.text  = element_text(size = 10)) +
  labs(title = "Species × species mean ANI (%)", x = "Species", y = "Species")

p_heat

p_heat <- ggplot(df_long, aes(x = species2, y = species1, fill = mean_ani)) +
  geom_tile(color = "grey80", linewidth = 0.2, width = 1, height = 1) +
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0)) +
  scale_fill_gradientn(
    colors   = pal_cols,
    na.value = "white",
    name     = "Mean ANI"
  ) +
  coord_fixed() +
  theme_minimal(base_size = 10, base_family = base_font) +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10,colour = "black", face = "italic"),
    axis.text.y = element_text(size = 10, colour = "black", face = "italic"),
    plot.title  = element_text(size = 10, face = "bold"),
    legend.title = element_text(size = 10, face = "bold"),
    axis.title.x.bottom = element_text(size = 10, face = "bold"),
    axis.title.y.left = element_text(size = 10, face = "bold"),
    legend.text  = element_text(size = 10)) +
  labs(title = "Species × species mean ANI (%)", x = "Species", y = "Species")

p_heat

thr <- median(df_long$mean_ani, na.rm = TRUE)
  
p_heat <- p_heat +
    geom_text(
      data = subset(df_long, !is.na(mean_ani)),
      aes(label = sprintf("%.1f", mean_ani), colour = mean_ani > thr),
      size = 2.5, show.legend = FALSE
    ) +
    scale_colour_manual(values = c("black", "white"))
  
out_pdf <- paste0(output_prefix, "_species_pairwise_heatmap.pdf")



# Packages
library(ggplot2)
library(dplyr)
library(readr)
library(tidyr)

#distance matrix
medoid_dist_mat <- as.matrix(medoid_dist)

med_long <- as.data.frame(as.table(medoid_dist_mat)) %>%
  rename(row_acc = Var1, col_acc = Var2, value = Freq) %>%
  mutate(row_acc = as.character(row_acc),
         col_acc = as.character(col_acc))

# Preserve the original order from the matrix
row_order <- rownames(medoid_dist_mat)
col_order <- colnames(medoid_dist_mat)

# Read accession _> species map
map <- read_csv("species_dist_species_medoids.csv", show_col_types = FALSE)

# Build lookup vectors (fallback to accession if a mapping is missing)
row_lookup <- setNames(map$species, map$medoid)[row_order]
row_lookup[is.na(row_lookup)] <- row_order[is.na(row_lookup)]

col_lookup <- setNames(map$species, map$medoid)[col_order]
col_lookup[is.na(col_lookup)] <- col_order[is.na(col_lookup)]

# If species names are duplicated, disambiguate by appending accession
if (any(duplicated(row_lookup))) {
  dups <- duplicated(row_lookup) | duplicated(row_lookup, fromLast = TRUE)
  row_lookup[dups] <- paste0(row_lookup[dups], " (", row_order[dups], ")")
}
if (any(duplicated(col_lookup))) {
  dups <- duplicated(col_lookup) | duplicated(col_lookup, fromLast = TRUE)
  col_lookup[dups] <- paste0(col_lookup[dups], " (", col_order[dups], ")")
}

# Map species to the long data frame and lock factor levels to original order
med_long <- med_long %>%
  mutate(
    row_species = factor(
      ifelse(!is.na(setNames(map$species, map$medoid)[row_acc]),
             setNames(map$species, map$medoid)[row_acc], row_acc),
      levels = unique(row_lookup)
    ),
    col_species = factor(
      ifelse(!is.na(setNames(map$species, map$medoid)[col_acc]),
             setNames(map$species, map$medoid)[col_acc], col_acc),
      levels = unique(col_lookup)
    )
  )

# Read  plotting order
ord <- readr::read_csv("species_order.csv", show_col_types = FALSE)

# Build a lookup: accession/medoid -> DISPLAY label actually used on the plot
label_by_acc <- setNames(row_lookup, row_order)

# Resolve order to the exact display labels used on the axes
if ("species" %in% names(ord)) {
  desired_labels <- ord$species
} else if ("medoid" %in% names(ord)) {
  desired_labels <- label_by_acc[ord$medoid]
} else if ("accession" %in% names(ord)) {
  desired_labels <- label_by_acc[ord$accession]
} else {
  stop("species_order.csv must contain a 'species', 'medoid', or 'accession' column.")
}
desired_labels <- unique(na.omit(as.character(desired_labels)))

# Labels actually present in the data (post-mapping)
present_labels <- unique(as.character(row_lookup))

# Final levels: first the desired order, then any leftovers not listed in the ordering csv
row_levels <- c(intersect(desired_labels, present_labels),
                setdiff(present_labels, desired_labels))

# Use the same order for columns so the heatmap reads consistently
col_levels <- row_levels

# Apply the ordering to the data used for plotting
med_long <- med_long %>%
  mutate(
    row_species = factor(as.character(row_species), levels = row_levels),
    col_species = factor(as.character(col_species), levels = col_levels)
  )

# Use ComplexHeatmap color function if available; otherwise a fallback palette.
val_range <- range(med_long$value, na.rm = TRUE)
pal_vals  <- seq(val_range[1], val_range[2], length.out = 256)
pal_cols  <- if (exists("col_fun")) {
  col_fun(pal_vals)                           # reuse your colorRamp2 mapping
} else {
  colorRampPalette(c("#f7fbff", "#c6dbef", "#6baed6", "#2171b5", "#08306b"))(256)
}

pal_cols_rev  <- colorRampPalette(c("#08306b", "#2171b5", "#6baed6", "#c6dbef", "#f7fbff"))(256)
pal_cols  <- colorRampPalette(c("#f7fbff", "#c6dbef", "#6baed6", "#2171b5", "#08306b"))(256)

# Plot
p_med <- ggplot(med_long, aes(x = col_species, y = row_species, fill = value)) +
  geom_tile(color = "grey80", linewidth = 0.2) +
  scale_fill_gradientn(
    colors   = pal_cols,
    na.value = "white",
    name     = "Distances"
  ) +
  labs(
    title = "Species medoid-to-medoid distances",
    x = "Species", y = "Species"
  ) +
  coord_equal(expand = FALSE) +
  theme_minimal(base_size = 16, base_family = base_font) +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 12,colour = "black", face = "italic"),
    axis.text.y = element_text(size = 10, colour = "black", face = "italic"),
    plot.title  = element_text(size = 10, face = "bold"),
    legend.title = element_text(size = 10, face = "bold"),
    axis.title.x.bottom = element_text(size = 10, face = "bold"),
    axis.title.y.left = element_text(size = 10, face = "bold"),
    legend.text  = element_text(size = 10)
  )

p_med

p_med <- ggplot(med_long, aes(x = col_species, y = row_species, fill = value)) +
  geom_tile(color = "grey80", linewidth = 0.2) +
  scale_fill_gradientn(
    colors   = pal_cols,
    na.value = "black",
    name     = "Distances"
  ) +
  labs(
    title = "Species medoid-to-medoid distances",
    x = "Species", y = "Species"
  ) +
  coord_equal(expand = FALSE) +
  theme_minimal(base_size = 16, base_family = base_font) +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10,colour = "black", face = "italic"),
    axis.text.y = element_text(size = 10, colour = "black", face = "italic"),
    plot.title  = element_text(size = 10, face = "bold"),
    legend.title = element_text(size = 10, face = "bold"),
    axis.title.x.bottom = element_text(size = 10, face = "bold"),
    axis.title.y.left = element_text(size = 10, face = "bold"),
    legend.text  = element_text(size = 10)
  )

p_med

p_med <- ggplot(med_long, aes(x = col_species, y = row_species, fill = value)) +
  geom_tile(color = "grey80", linewidth = 0.2, width = 1, height = 1) +
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0)) +
  scale_fill_gradientn(colors = pal_cols, na.value = "black", name = "Distances") +
  labs(title = "Species medoid-to-medoid distances", x = "Species", y = "Species") +
  coord_equal(expand = FALSE) +
  theme_minimal(base_size = 16, base_family = base_font) +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10, colour = "black", face = "italic"),
    axis.text.y = element_text(size = 10, colour = "black", face = "italic"),
    plot.title  = element_text(size = 10, face = "bold"),
    legend.title = element_text(size = 10, face = "bold"),
    axis.title.x.bottom = element_text(size = 10, face = "bold"),
    axis.title.y.left = element_text(size = 10, face = "bold"),
    legend.text  = element_text(size = 10)
  )

p_med

thr <- median(med_long$value, na.rm = TRUE)

p_med +
  geom_text(
    data = subset(med_long, !is.na(value)),
    aes(label = sprintf("%.2f", value), colour = value > thr),
    size = 2.2, fontface = "bold", show.legend = FALSE
  ) +
  scale_colour_manual(values = c("black", "white"))

p_med <- p_med +
  geom_text(
    data = subset(med_long, !is.na(value)),
    aes(label = sprintf("%.2f", value), colour = value > thr),
    size = 2.5, show.legend = FALSE
  ) +
  scale_colour_manual(values = c("black", "white"))

p_med

out_pdf <- paste0(output_prefix, "_species_medoid-to-medoid_dist_heatmap.pdf")



#ANI matrix
medoid_ani_mat <- as.matrix(medoid_ani)

med_long_ani <- as.data.frame(as.table(medoid_ani_mat)) %>%
  rename(row_acc = Var1, col_acc = Var2, value = Freq) %>%
  mutate(row_acc = as.character(row_acc),
         col_acc = as.character(col_acc))

# Preserve the original order from the matrix
row_order <- rownames(medoid_ani_mat)
col_order <- colnames(medoid_ani_mat)

# Read accession → species map
# CSV should have columns: accession, species
map <- read_csv("species_dist_species_medoids.csv", show_col_types = FALSE)

# Build lookup vectors (fallback to accession if a mapping is missing)
row_lookup <- setNames(map$species, map$medoid)[row_order]
row_lookup[is.na(row_lookup)] <- row_order[is.na(row_lookup)]

col_lookup <- setNames(map$species, map$medoid)[col_order]
col_lookup[is.na(col_lookup)] <- col_order[is.na(col_lookup)]

# If species names are duplicated, disambiguate by appending accession
if (any(duplicated(row_lookup))) {
  dups <- duplicated(row_lookup) | duplicated(row_lookup, fromLast = TRUE)
  row_lookup[dups] <- paste0(row_lookup[dups], " (", row_order[dups], ")")
}
if (any(duplicated(col_lookup))) {
  dups <- duplicated(col_lookup) | duplicated(col_lookup, fromLast = TRUE)
  col_lookup[dups] <- paste0(col_lookup[dups], " (", col_order[dups], ")")
}

# Map species to the long data frame and lock factor levels to original order
med_long_ani <- med_long_ani %>%
  mutate(
    row_species = factor(
      ifelse(!is.na(setNames(map$species, map$medoid)[row_acc]),
             setNames(map$species, map$medoid)[row_acc], row_acc),
      levels = unique(row_lookup)
    ),
    col_species = factor(
      ifelse(!is.na(setNames(map$species, map$medoid)[col_acc]),
             setNames(map$species, map$medoid)[col_acc], col_acc),
      levels = unique(col_lookup)
    )
  )


# Read desired plotting order
# species_order.csv should have a column named 'species' OR 'medoid' (accession).
ord <- readr::read_csv("species_order.csv", show_col_types = FALSE)

# Build a lookup: accession/medoid -> DISPLAY label actually used on the plot
label_by_acc <- setNames(row_lookup, row_order)

# Resolve the desired order to the exact display labels used on the axes
if ("species" %in% names(ord)) {
  desired_labels <- ord$species
} else if ("medoid" %in% names(ord)) {
  desired_labels <- label_by_acc[ord$medoid]
} else if ("accession" %in% names(ord)) {
  desired_labels <- label_by_acc[ord$accession]
} else {
  stop("species_order.csv must contain a 'species', 'medoid', or 'accession' column.")
}
desired_labels <- unique(na.omit(as.character(desired_labels)))

# Labels actually present in the data (post-mapping)
present_labels <- unique(as.character(row_lookup))

# Final levels: first the desired order, then any leftovers not listed in the CSV
row_levels <- c(intersect(desired_labels, present_labels),
                setdiff(present_labels, desired_labels))

# Use the same order for columns so the heatmap reads consistently
col_levels <- row_levels

# Apply the ordering to the data used for plotting
med_long_ani <- med_long_ani %>%
  mutate(
    row_species = factor(as.character(row_species), levels = row_levels),
    col_species = factor(as.character(col_species), levels = col_levels)
  )

#  Colors 
# Use ComplexHeatmap color function if available; otherwise a fallback palette.
val_range <- range(med_long_ani$value, na.rm = TRUE)
pal_vals  <- seq(val_range[1], val_range[2], length.out = 256)
pal_cols  <- if (exists("col_fun")) {
  col_fun(pal_vals)                           # reuse your colorRamp2 mapping
} else {
  colorRampPalette(c("#f7fbff", "#c6dbef", "#6baed6", "#2171b5", "#08306b"))(256)
}

pal_cols  <- colorRampPalette(c("#f7fbff", "#c6dbef", "#6baed6", "#2171b5", "#08306b"))(256)

p_med_ani <- ggplot(med_long_ani, aes(x = col_species, y = row_species, fill = value)) +
  geom_tile(color = "grey80", linewidth = 0.2, width = 1, height = 1) +
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0)) +
  scale_fill_gradientn(
    colors   = pal_cols,
    na.value = "white",
    name     = "ANI"
  ) +
  labs(
    title = "Species medoid-to-medoid ANI",
    x = "Species", y = "Species"
  ) +
  coord_equal(expand = FALSE) +
  theme_minimal(base_size = 16, base_family = base_font) +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10,colour = "black", face = "italic"),
    axis.text.y = element_text(size = 10, colour = "black", face = "italic"),
    plot.title  = element_text(size = 10, face = "bold"),
    legend.title = element_text(size = 10, face = "bold"),
    axis.title.x.bottom = element_text(size = 10, face = "bold"),
    axis.title.y.left = element_text(size = 10, face = "bold"),
    legend.text  = element_text(size = 10)
  )

p_med_ani

thr <- median(med_long_ani$value, na.rm = TRUE)

p_med_ani <- p_med_ani +
  geom_text(
    data = subset(med_long_ani, !is.na(value)),
    aes(label = sprintf("%.1f", value), colour = value > thr),
    size = 2.3, show.legend = FALSE
  ) +
  scale_colour_manual(values = c("black", "white"))

p_med_ani

out_pdf <- paste0(output_prefix, "_species_medoid-to-medoid_ani_heatmap.pdf")

#convex hull matrix
mds_df <- merge(mds_df, meta[, c("accession","species")], by="accession", all.x=TRUE)

library(dplyr)
library(purrr)

# ensure unique column names (prevents dplyr from erroring)
names(mds_df) <- make.unique(names(mds_df), sep = "_dup")

# grab all species-like columns
sp_cols <- grep("^species", names(mds_df), value = TRUE)
if (length(sp_cols) == 0) stop("No species* columns found in mds_df.")

# build a clean 'species' by taking the first non-NA/non-empty across the species* columns
tmp <- mds_df[sp_cols] %>%
  mutate(across(everything(), ~{
    x <- as.character(.); x[trimws(x) == ""] <- NA_character_; x
  }))

mds_df$species <- reduce(tmp, coalesce)

# make sure MDS cols are numeric
mds_df <- mds_df %>% mutate(MDS1 = as.numeric(MDS1), MDS2 = as.numeric(MDS2))

message("Assigned species to ", sum(!is.na(mds_df$species)), " / ", nrow(mds_df), " rows.")




library(sf)
# a robust "no CRS" object for this sf version
crs_none <- tryCatch(sf::NA_crs_, error = function(e) sf::st_crs(NA))

# ensure MDS are numeric
mds_ok <- mds_df %>% dplyr::filter(!is.na(species), is.finite(MDS1), is.finite(MDS2)) %>%
  dplyr::mutate(MDS1 = as.numeric(MDS1), MDS2 = as.numeric(MDS2))

# tiny buffer scale for degenerate hulls
rng_x <- diff(range(mds_ok$MDS1, na.rm = TRUE))
rng_y <- diff(range(mds_ok$MDS2, na.rm = TRUE))
eps   <- 1e-6 * max(rng_x, rng_y); if (!is.finite(eps) || eps <= 0) eps <- 1e-6

species_polys <- list(); status <- list()

for (sp in sort(unique(mds_ok$species))) {
  pts <- mds_ok %>% dplyr::filter(species == sp) %>% dplyr::select(MDS1, MDS2) %>%
    dplyr::distinct() %>% tidyr::drop_na()
  if (nrow(pts) == 0) { species_polys[[sp]] <- NULL; next }
  
  # build MULTIPOINT and convex hull
  p_sf <- sf::st_as_sf(pts, coords = c("MDS1","MDS2"), remove = FALSE, crs = crs_none)
  mp   <- sf::st_union(p_sf)  # MULTIPOINT
  hull <- tryCatch(sf::st_convex_hull(mp), error = function(e) NULL)
  
  # geometry type + area
  gt   <- if (is.null(hull)) NA_character_ else as.character(sf::st_geometry_type(hull))
  area <- if (is.null(hull)) NA_real_ else suppressWarnings(as.numeric(sf::st_area(hull)))
  
  # if not a polygon or area <= 0, buffer the points slightly
  if (is.null(hull) || !(gt %in% c("POLYGON","MULTIPOLYGON")) || !is.finite(area) || area <= 0) {
    hull <- tryCatch(sf::st_buffer(mp, eps), error = function(e) NULL)
    gt   <- if (is.null(hull)) NA_character_ else as.character(sf::st_geometry_type(hull))
    area <- if (is.null(hull)) NA_real_ else suppressWarnings(as.numeric(sf::st_area(hull)))
  }
  
  if (!is.null(hull)) hull <- tryCatch(sf::st_make_valid(hull), error=function(e) hull)
  
  species_polys[[sp]] <- hull
  status[[sp]] <- data.frame(species = sp, n_points = nrow(pts),
                             geom_type = gt, area = area, stringsAsFactors = FALSE)
}

hull_status <- dplyr::bind_rows(status)
write.csv(hull_status, paste0(output_prefix, "_hull_status.csv"), row.names = FALSE)
message("Hull types:"); print(table(hull_status$geom_type, useNA = "ifany"))

library(sf)

# Sum area of any sf geometry into a single numeric scalar (or NA if not computable)
area_sum <- function(g) {
  if (is.null(g) || length(g) == 0) return(NA_real_)
  # if it's a collection, extract polygons; if none, keep original
  g2 <- tryCatch(suppressWarnings(st_collection_extract(g, "POLYGON")), error = function(e) g)
  if (length(g2) == 0) g2 <- g
  if (all(suppressWarnings(st_is_empty(g2)))) return(0)  # empty -> zero area
  a <- tryCatch(sum(as.numeric(st_area(g2))), error = function(e) NA_real_)
  if (!is.finite(a)) NA_real_ else a
}

safe_overlap <- function(g1, g2) {
  inter <- tryCatch(st_intersection(g1, g2), error = function(e) NULL)
  union_ <- tryCatch(st_union(g1, g2),       error = function(e) NULL)
  ia <- area_sum(inter)
  ua <- area_sum(union_)
  if (!is.finite(ua) || ua <= 0) return(NA_real_)
  if (!is.finite(ia) || ia < 0)  return(NA_real_)
  ia / ua
}

spn <- names(species_polys)
overlap <- matrix(NA_real_, length(spn), length(spn), dimnames = list(spn, spn))
for (i in seq_along(spn)) {
  for (j in seq_along(spn)) {
    p1 <- species_polys[[spn[i]]]; p2 <- species_polys[[spn[j]]]
    if (is.null(p1) || is.null(p2)) { overlap[i, j] <- NA_real_; next }
    overlap[i, j] <- safe_overlap(p1, p2)
  }
}

write.csv(overlap, paste0(output_prefix, "_hull_overlap_matrix.csv"), row.names = TRUE)
message("Wrote polygon-overlap: ", paste0(output_prefix, "_hull_overlap_matrix.csv"))

# sanity checks
message("NA cells in overlap: ", sum(is.na(overlap)))
message("Diagonal (should be 1's or close): ", paste(round(diag(overlap), 3), collapse = ", "))


library(sf); crs_none <- tryCatch(sf::NA_crs_, error=function(e) sf::st_crs(NA))

# Build point-in-hull matrix
pts_sf <- st_as_sf(mds_ok %>% dplyr::select(accession, species, MDS1, MDS2) %>% dplyr::distinct(),
                   coords = c("MDS1","MDS2"), crs = crs_none)

sp <- names(species_polys)
A_in_B <- matrix(NA_real_, length(sp), length(sp), dimnames = list(sp, sp))
for (i in seq_along(sp)) {
  ptsA <- pts_sf %>% dplyr::filter(species == sp[i])
  if (nrow(ptsA)==0) next
  for (j in seq_along(sp)) {
    polyB <- species_polys[[sp[j]]]
    if (is.null(polyB) || length(polyB)==0) next
    inside <- lengths(st_within(ptsA, polyB)) > 0
    A_in_B[i,j] <- sum(inside) / nrow(ptsA)
  }
}
sym <- (A_in_B + t(A_in_B)) / 2

# Fill NA cells in polygon-overlap with the symmetric point-in-hull
overlap_filled <- overlap
overlap_filled[is.na(overlap_filled)] <- sym[rownames(overlap_filled), colnames(overlap_filled)][is.na(overlap_filled)]

areas <- setNames(hull_status$area, hull_status$species)
pair_has_polys <- outer(rownames(overlap_filled), colnames(overlap_filled),
                        Vectorize(function(a,b) is.finite(areas[a]) && areas[a]>0 && is.finite(areas[b]) && areas[b]>0))
overlap_filled[is.na(overlap_filled) & pair_has_polys] <- 0

write.csv(overlap_filled, paste0(output_prefix, "_hull_overlap_matrix_filled.csv"), row.names = TRUE)

library(writexl)

write_xlsx(med_long, "medoid-to-medoid-long-dist.xlsx")
write_xlsx(med_long_ani, "medoid-to-medoid-long-ani.xlsx")
