# =====================================================
# Species Co-occurrence Analysis Script (Mammals) — Baluran NP
# Location    : Baluran National Park (East Java, Indonesia)
# Data        : co_occur_bnp_56.csv (Presence–Absence Matrix)
# Software    : R (version 2025.06.13 ucrt)
# Purpose     : Compute pairwise co-occurrence among species 
# Method      : Probabilistic Co-occurrence (Griffith et al., 2016; R package `cooccur`)
# =====================================================

# -----------------------------
# Load required libraries
# -----------------------------
library(cooccur)
library(ggplot2)
library(reshape2)
library(RColorBrewer)

# (Optional) Cite the `cooccur` package
# citation("cooccur")

# -----------------------------------------------------
# Read the presence–absence .csv (rows = species; cols = camera stations)
# Notes:
# - The first column should contain species names (header: 'species')
# - Values must be 0/1 for absence/presence
# - File is assumed to be in the current working directory
# -----------------------------------------------------
data <- read.csv("co_occur_bnp_56.csv", sep = ";", header = TRUE, stringsAsFactors = FALSE)
species_names <- as.character(data$species)

# Convert to a binary matrix (0/1), excluding the species name column
species_data <- data[, -1]
species_data <- as.data.frame(lapply(species_data, as.numeric))

# >>> Change #1: attach species names as rownames so they propagate to results
rownames(species_data) <- species_names

# -----------------------------------------------------
# Run probabilistic co-occurrence analysis
# -----------------------------------------------------
# `cooccur` will compute observed & expected co-occurrences under a hypergeometric model.
cooccur_result <- cooccur(species_data)

# Extract the results table and build a symmetric O (observed joint stations) matrix
cooccur_table <- cooccur_result$results
num_species <- length(species_names)

# Initialize an empty symmetric matrix with species names as dimnames
co_occurrence_matrix <- matrix(0, nrow = num_species, ncol = num_species,
                               dimnames = list(species_names, species_names))

# Fill the matrix with observed joint-station counts (O_ij)
for (i in 1:nrow(cooccur_table)) {
  sp1_index <- cooccur_table$sp1[i]
  sp2_index <- cooccur_table$sp2[i]
  obs_cooccur <- cooccur_table$obs_cooccur[i]
  sp1 <- species_names[sp1_index]
  sp2 <- species_names[sp2_index]
  co_occurrence_matrix[sp1, sp2] <- obs_cooccur
  co_occurrence_matrix[sp2, sp1] <- obs_cooccur
}

# -----------------------------------------------------
# Heatmap visualization
# -----------------------------------------------------
if (sum(co_occurrence_matrix) > 0) {
  cooccur_matrix_melted <- melt(co_occurrence_matrix)
  colnames(cooccur_matrix_melted) <- c("Species1", "Species2", "CoDetection_Count")

  ggplot(cooccur_matrix_melted, aes(x = Species1, y = Species2, fill = CoDetection_Count)) +
    geom_tile() +
    scale_fill_gradientn(colors = RColorBrewer::brewer.pal(9, "Blues"), name = "Co-detection") +
    geom_text(aes(label = round(CoDetection_Count, 2)), color = "black", size = 3) +
    theme_minimal() +
    labs(title = "Species Co-detection Heatmap", x = "Species", y = "Species") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, face = "italic", size = 8),
          axis.text.y = element_text(angle = 0, hjust = 1, face = "italic", size = 8))
} else {
  message("The co-occurrence matrix is empty or does not contain sufficient data.")
}

# -----------------------------------------------------
# Export tables
# 1) Pairwise (probabilistic) table from `cooccur`
# -----------------------------------------------------
pairs <- cooccur_result$results

# >>> Change #2: ensure species names are exported along with indices
if (!("sp1_name" %in% names(pairs) && "sp2_name" %in% names(pairs))) {
  pairs$sp1_name <- species_names[pairs$sp1]
  pairs$sp2_name <- species_names[pairs$sp2]
}
pairs_out <- pairs[, c("sp1_name", "sp2_name",
                       setdiff(names(pairs), c("sp1", "sp2", "sp1_name", "sp2_name")))]
names(pairs_out)[1:2] <- c("Species1", "Species2")
write.csv(pairs_out, "table_S5_pairs_results_BNP.csv", row.names = FALSE)

# -----------------------------------------------------
# 2) Symmetric O (observed joint-station counts) matrix
# -----------------------------------------------------
write.csv(co_occurrence_matrix, "table_S5_obs_counts_matrix_BNP.csv", row.names = TRUE)

# -----------------------------------------------------
# 3) Summary table: C1 (% park-wide), Score (1–3), n_i, Mean/SD/Range (%)
# -----------------------------------------------------
N  <- ncol(species_data)               # total stations (park-wide)
M  <- co_occurrence_matrix             # O_ij counts
diag(M) <- NA                          # exclude diagonal for means

# C1_i (% park-wide) = mean_j!=i [100 * O_ij / N]
C1 <- rowMeans(M / N * 100, na.rm = TRUE)

# n_i = stations where species i detected
n_i <- rowSums(species_data)

# Coverage_i→j (%) = 100 * O_ij / n_j, summarized across j!=i
nj <- n_i
D  <- matrix(nj, nrow = length(nj), ncol = length(nj), byrow = TRUE); D[D==0] <- NA
covg <- (M / D) * 100
Mean <- round(rowMeans(covg, na.rm = TRUE), 1)
SD   <- round(apply(covg, 1, sd, na.rm = TRUE), 1)
rng  <- t(apply(covg, 1, function(x) if (all(is.na(x))) c(NA, NA) else range(x, na.rm = TRUE)))
Range <- paste0(round(rng[,1], 0), "–", round(rng[,2], 0))

# Tertiles (tie-friendly: exactly Q2 → Score 2)
Q1 <- as.numeric(quantile(C1, 1/3, na.rm = TRUE))
Q2 <- as.numeric(quantile(C1, 2/3, na.rm = TRUE))
Score <- ifelse(C1 <= Q1, 1, ifelse(C1 <= Q2, 2, 3))

# Assemble and export
tab <- data.frame(
  Species = rownames(co_occurrence_matrix),
  `C1 (% park-wide)` = round(C1, 2),
  `Score (1–3; tertiles, tie-friendly)` = Score,
  `#Stations detected (n_i)` = as.integer(n_i),
  `Mean (%)` = Mean, `SD (%)` = SD, `Range (%)` = Range,
  check.names = FALSE
)
tab <- tab[order(-tab$`C1 (% park-wide)`, tab$Species), ]
write.csv(tab, "table_S5_C1_parkwide_BNP.csv", row.names = FALSE)
cat(sprintf("Thresholds: Q1 = %.4f%%, Q2 = %.4f%% (exactly Q2 → Score 2)\n", Q1, Q2))

