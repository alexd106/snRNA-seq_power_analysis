# Simulation based power analysis for single nuclei RNA-seq
# For: Speakman, Mitchell (BBSRC application)
# By: Alex

# Date created: 2025-08-01 09:15:00
# Date last modified: 2025-08-12 10:37:00 

# GENERAL APPROACH & RATIONALE:
#
# 1) Goal:
#    Estimate statistical power to detect differential expression (DE) 
#    between two biological groups across a range of sample sizes and 
#    numbers of cells per sample, for realistic single-cell/nucleus data.
#
# 2) Rationale:
#    - Power depends jointly on biological replication (samples per group)
#      and technical replication (cells per sample).
#    - Simulations allow us to control these factors and quantify their 
#      effect on power in a controlled setting.
#    - Using parameters estimated from a real dataset ensures simulated data 
#      reflect realistic mean–variance relationships and cluster composition.
#
# 3) Strategy:
#    a) Load and trim an existing single-nucleus RNA-seq dataset to only 
#       include relevant groups, top clusters, and a compact sparse counts assay.
#    b) Filter to highly variable and sufficiently expressed genes to reduce 
#       computational cost without distorting mean–variance patterns.
#    c) Fit negative binomial parameters for each gene × cluster using muscat::prepSim().
#    d) For each (sample size, cells per sample) combination:
#         - Simulate data using muscat::simData() with a known proportion of 
#           truly DE genes and defined log2 fold-change.
#         - Aggregate counts per (cluster, sample) to perform pseudobulk DE.
#         - Compare DE test results to the known truth to compute power.
#    e) Summarise power across clusters and visualise as:
#         - Heatmap: mean power vs. sample size and cells per sample.
#         - Line plots: power by cluster over sample sizes.
#
# 4) Benefits:
#    - Produces empirically grounded power estimates without requiring new data.
#    - Informs study design by revealing trade-offs between sequencing depth 
#      (cells) and biological replication (samples).
#    - Scales better by trimming and sparsifying the reference dataset before 
#      parameter estimation.

# NOTE
# Data to parameterise simulation obtained from :
# https://www.repository.cam.ac.uk/items/8f9c3683-29fd-44f3-aad5-7acf5e963a75

# Lam, Y. H., Yeo, G., Steuernagel, L., Klemm, P., & Brüning, J. (2022). Research data 
# supporting “HypoMap – a unified single cell gene expression atlas of the murine hypothalamus”. 
# Apollo - University of Cambridge Repository. https://doi.org/10.17863/CAM.87955
#
# Downloaded hypoMap.rds (4.07 GB) for import into R in ./data/ on 2025-08-04

# ---- Dependencies ----
library(muscat)
library(Seurat)
library(SingleCellExperiment)
library(SummarizedExperiment)
library(scater)
library(edgeR)
library(Matrix)      # sparse matrices
library(ggplot2)
library(dplyr)
library(data.table)

suppressPackageStartupMessages({
  if (!requireNamespace("scran", quietly = TRUE)) {
    message("Package 'scran' not installed; falling back to simple HVG filter.")
  }
})

# =========================
# Import, trim, and prepare ref
# =========================
seurat_obj <- readRDS("data/hypoMap.rds")
sce_full   <- as.SingleCellExperiment(seurat_obj)
rm(seurat_obj); gc()

# Minimal colData
cd <- colData(sce_full)
cd$sample_id  <- as.factor(cd$Sample_ID)
cd$group_id   <- as.factor(cd$Diet)
cd$cluster_id <- as.factor(cd$Author_Class_Curated)
colData(sce_full) <- cd[, c("sample_id", "group_id", "cluster_id"), drop = FALSE]

# Keep comparison groups and top 5 clusters
keep_grp <- sce_full$group_id %in% c("Fasted", "Normal chow")
sce_diet <- sce_full[, keep_grp, drop = FALSE]
sce_diet$group_id <- droplevels(sce_diet$group_id)
rm(sce_full); gc()

top5 <- names(sort(table(sce_diet$cluster_id), decreasing = TRUE))[1:5]
sce  <- sce_diet[, sce_diet$cluster_id %in% top5, drop = FALSE]
sce$cluster_id <- droplevels(sce$cluster_id)
rm(sce_diet); gc()

# Keep only sparse counts; drop heavy baggage
stopifnot("counts" %in% assayNames(sce))
assay(sce, "counts") <- as(assay(sce, "counts"), "dgCMatrix")
assayNames(sce) <- "counts"
reducedDims(sce) <- NULL
altExpNames(sce) <- character(0)
metadata(sce)    <- list()
rowData(sce)     <- DataFrame()
gc()

# Gene filtering BEFORE prepSim
target_ng <- 1000L
expr_frac <- Matrix::rowMeans(assay(sce, "counts") > 0)
sce_f     <- sce[expr_frac >= 0.01, ]

if (requireNamespace("scran", quietly = TRUE)) {
  hvg <- scran::modelGeneVar(sce_f, block = sce_f$cluster_id)
  ord <- order(hvg$bio, decreasing = TRUE)
  keep_ng <- head(rownames(sce_f)[ord], n = min(target_ng, nrow(sce_f)))
} else {
  cpm <- edgeR::cpm(assay(sce_f, "counts"), log = TRUE, prior.count = 1)
  v   <- matrixStats::rowVars(cpm)
  keep_ng <- head(rownames(sce_f)[order(v, decreasing = TRUE)], n = min(target_ng, nrow(sce_f)))
  rm(cpm, v)
}
sce_f <- sce_f[keep_ng, ]
gc()

# Prepare muscat reference on trimmed object
ref <- prepSim(
  sce_f,
  group_keep = c("Fasted", "Normal chow"),
  min_genes  = 200
)

# =========================
# Grid & global params
# =========================
sample_sizes <- 4:10                  # samples per group
cell_counts  <- seq(1000, 5000, 200)  # cells per sample
ngenes       <- min(1000L, nrow(ref)) # <= nrow(ref)
nclusters    <- length(unique(colData(ref)$cluster_id))
lfc_val      <- 2
fdr_cutoff   <- 0.05
# p_dd must have length 6 in this muscat version: (ee, ep, de, dp, dm, db)
p_dd <- c(0.9, 0, 0.1, 0, 0, 0)

# =========================
# Helpers (kept minimal)
# =========================

# Enforce minimum cells per (sample, cluster) so pbDS won't drop everything
min_cells_sc <- 10
keep_min_cells <- function(sim, min_cells = 10) {
  cd  <- as.data.frame(colData(sim))
  key <- interaction(cd$sample_id, cd$cluster_id, drop = TRUE)
  n   <- ave(seq_along(key), key, FUN = length)
  sim[, n >= min_cells]
}

# Standardise gene_info to gene_id / cluster_id / category
get_gene_info_std <- function(sim) {
  gi <- metadata(sim)$gene_info
  nm <- colnames(gi)
  g_col <- if ("gene_id" %in% nm) "gene_id" else if ("gene" %in% nm) "gene" else stop("No gene column in gene_info")
  k_col <- if ("cluster_id" %in% nm) "cluster_id" else if ("cluster" %in% nm) "cluster" else stop("No cluster column in gene_info")
  c_col <- if ("category" %in% nm) "category" else if ("status" %in% nm) "status" else stop("No category/status column in gene_info")
  data.frame(
    gene_id    = as.character(gi[[g_col]]),
    cluster_id = as.character(gi[[k_col]]),
    category   = as.character(gi[[c_col]]),
    stringsAsFactors = FALSE
  )
}

# Detect FDR/p-adj column in pbDS table
detect_fdr_col <- function(df) {
  cand <- c("FDR","fdr","p_adj.loc","p_adj.glb","p_adj_global","padj","q_value","qvalue")
  hit  <- cand[cand %in% colnames(df)]
  if (length(hit)) hit[1] else NA_character_
}

# Detect gene ID column in pbDS table
detect_gene_col_tbl <- function(df) {
  cand <- c("gene_id","gene","feature_id","symbol","id")
  hit  <- cand[cand %in% colnames(df)]
  if (length(hit)) hit[1] else NA_character_
}

# Normalise cluster labels like "cluster3" -> "3"
norm_cl <- function(x) sub("^cluster", "", as.character(x))

# =========================
# One run with retries until each cluster has >=1 true DE gene
# =========================
one_run <- function(ns, nc, max_tries = 10) {
  message(sprintf("Simulating ns=%d, nc=%d", ns, nc))

  # Retry loop to ensure DE per cluster in muscat truth
  tried <- 0L
  repeat {
    tried <- tried + 1L
    sim <- simData(
      ref, ns = ns, nc = nc, ng = ngenes, nk = nclusters,
      p_dd = p_dd, lfc = lfc_val, paired = FALSE
    )
    gi_std <- get_gene_info_std(sim)
    de_counts <- table(norm_cl(gi_std$cluster_id)[gi_std$category == "de"])
    sim_cls <- unique(norm_cl(as.character(colData(sim)$cluster_id)))
    ok <- length(sim_cls) > 0 && all(sim_cls %in% names(de_counts)) && all(de_counts[sim_cls] > 0)
    if (ok || tried >= max_tries) {
      if (!ok) warning("Proceeding without DE in every cluster after max_tries.")
      break
    }
    set.seed(sample.int(.Machine$integer.max, 1))
  }

  # Avoid pbDS dropping clusters due to extreme sparsity
  sim <- keep_min_cells(sim, min_cells_sc)

  # Pseudobulk & DE
  pb  <- aggregateData(sim, assay = "counts", fun = "sum",
                       by = c("cluster_id", "sample_id"))
  de  <- pbDS(pb)

  tbl_list <- de$table[[1]]
  if (is.null(tbl_list) || length(tbl_list) == 0) {
    warning("pbDS returned no cluster tables; skipping this run.")
    return(data.table::data.table(ns=ns, nc=nc, cluster=character(), power=numeric()))
  }

  fdr_col  <- detect_fdr_col(tbl_list[[1]])
  if (is.na(fdr_col)) stop("No FDR/p-adj column in pbDS table: ", paste(colnames(tbl_list[[1]]), collapse=", "))
  gene_col <- detect_gene_col_tbl(tbl_list[[1]])
  if (is.na(gene_col)) stop("No gene ID column in pbDS table: ", paste(colnames(tbl_list[[1]]), collapse=", "))

  # Truth (standardised + normalised cluster labels)
  gi_std <- get_gene_info_std(sim)
  gi_std$cl_norm <- norm_cl(gi_std$cluster_id)

  # Compute power per cluster using the label INSIDE each table
  out <- lapply(tbl_list, function(tbl) {
    cl_lab   <- as.character(unique(tbl$cluster_id))[1]
    cl_normv <- norm_cl(cl_lab)

    degenes  <- gi_std$gene_id[gi_std$cl_norm == cl_normv & gi_std$category == "de"]
    detected <- tbl[[gene_col]][ !is.na(tbl[[fdr_col]]) & tbl[[fdr_col]] <= fdr_cutoff ]

    pow <- if (!length(degenes)) 0 else sum(degenes %in% detected) / length(degenes)
    data.frame(ns = ns, nc = nc, cluster = cl_lab, power = pow, stringsAsFactors = FALSE)
  })

  # FYI if pbDS dropped clusters
  sim_cls <- unique(norm_cl(as.character(colData(sim)$cluster_id)))
  present_cls <- unique(norm_cl(vapply(tbl_list, function(t) as.character(unique(t$cluster_id)[1]), character(1))))
  dropped <- setdiff(sim_cls, present_cls)
  if (length(dropped)) {
    message("pbDS dropped clusters (few cells per sample): ", paste(dropped, collapse = ", "))
  }

  rm(sim, pb, de, tbl_list, gi_std); gc(FALSE)
  data.table::rbindlist(out, use.names = TRUE, fill = TRUE)
}

# =========================
# Run grid
# =========================
results_list <- vector("list", length(sample_sizes) * length(cell_counts))
k <- 0L
for (ns in sample_sizes) {
  for (nc in cell_counts) {
    k <- k + 1L
    results_list[[k]] <- one_run(ns, nc)
    if (k %% 10 == 0) gc(verbose = FALSE)
  }
}

results <- data.table::rbindlist(results_list, use.names = TRUE, fill = TRUE)

# =========================
# Summaries & plots
# =========================
# Mean power by ns × nc (summary table)
summary_power <- results %>%
  group_by(ns, nc) %>%
  summarise(mean_power = mean(power, na.rm = TRUE), .groups = "drop") %>%
  arrange(ns, nc)

print(head(summary_power, 10))
# Save summary as CSV for the grant appendix
write.csv(summary_power, file = "power_summary_ns_nc.csv", row.names = FALSE)

# Heatmap of mean power
p1 <- ggplot(summary_power, aes(x = factor(ns), y = nc, fill = mean_power)) +
  geom_tile(color = "white", linewidth = 0.2) +
  scale_fill_gradient(low = "white", high = "steelblue", na.value = "grey90") +
  labs(x = "Samples per group", y = "Cells per sample", fill = "Mean power",
       title = "Power vs. sample size and cells") +
  theme_minimal(base_size = 12)

# Power by cluster (optional diagnostic)
p2 <- ggplot(results, aes(x = ns, y = power, color = cluster, group = cluster)) +
  geom_line() +
  facet_wrap(~ nc, labeller = label_both) +
  theme_minimal(base_size = 11) +
  labs(title = "Power by cluster", x = "Samples per group", y = "Power")

print(p1); print(p2)


# SUGGESTED TEXT

# We conducted a simulation-based power analysis for single-nucleus RNA-seq using parameters estimated from the murine hypothalamus atlas HypoMap (DOI: 10.17863/CAM.87955). The reference dataset was trimmed to the two study-relevant groups and five most abundant cell types, retaining only sparse counts for highly variable genes. Negative binomial parameters were estimated using the muscat R package, and simulations were run across a grid of biological replicates and cells per sample, introducing 10% truly DE genes with fixed log₂ fold-change. Pseudobulk DE testing with edgeR was applied, and power was calculated as the proportion of true DE genes detected at FDR ≤ 0.05. For a design with 8 mice per group and 3,600 cells per sample, the estimated mean power across clusters was [XX%], with cluster-specific power ranging from [YY%] to [ZZ%].

# GENERAL SUMMARY

# 1) Load and convert
#    Read the Seurat .rds file and convert it to a SingleCellExperiment (sce_full).
#    Immediately drop the Seurat object to free RAM.

# 2) Keep only essential metadata
#    In colData, keep just three factors: sample_id, group_id, and cluster_id.
#    This strips out lots of per-cell columns you don’t need for power simulations.

# 3) Subset to the comparison and most abundant cell types
#    Keep only the two groups to compare: “Fasted” vs “Normal chow”.
#    Identify the top 5 most frequent clusters and keep only those.
#    Dropping rarer clusters speeds things up and stabilizes pseudobulk DE.

# 4) Slim the assays to a single sparse count matrix
#    Ensure there’s a counts assay; convert it to a sparse dgCMatrix.
#    Remove other assays (logcounts, data, etc.), reducedDims, altExps, and heavy metadata.
#    Result: a tiny SCE that’s just counts + 3 columns of colData.

# 5) Pre-filter genes before parameter fitting
#    Keep genes expressed in at least ~1% of cells (tuneable).
#    From those, select up to target_ng highly variable genes (HVGs):
#       - Prefer scran::modelGeneVar(..., block = cluster_id) to avoid bias to one cell type.
#       - Fallback: rank by variance of log-CPM.
#    This ensures the next step fits NB parameters on a small, informative gene set.

# 6) Build the muscat reference on the trimmed object
#    Run prepSim() on the slimmed SCE, keeping only the two groups;
#    set min_genes = 200 to discard ultra-sparse cells quickly.
#    This ref stores the compact NB/mean–dispersion estimates muscat needs for fast simData() calls.

# Why this matters for power analysis:
#    The heavy RAM/time hit usually comes from running prepSim() on a fat SCE
#    (many assays, metadata, lowly expressed genes).
#    Trimming to (i) two groups, (ii) top clusters, (iii) sparse counts only,
#    (iv) HVGs only shrinks the parameter-learning step from hundreds of GB
#    to single-digit GB and makes simulations snappy.

# Quick sanity checks to run once:
#    - dim(sce) and dim(sce_f) to confirm cell/gene reductions look reasonable.
#    - table(sce$group_id, sce$cluster_id) to ensure both groups are present across clusters.
#    - nrow(ref) (or nrow(rowData(ref))) should be ≥ your planned ng in simData();
#      if not, drop ng or relax filters.

