# =========================================================================
# scRNA-seq Best-Practices (Seurat v5 + Harmony + SCT v2 + SingleR + scuttle)
# v5-layer-safe accessors, S4 checks, memory-safe SCT, Leiden clustering
# =========================================================================

suppressPackageStartupMessages({
  if (!requireNamespace("pacman", quietly = TRUE)) install.packages("pacman")
  library(pacman)
})

options(repos = c(CRAN = "https://cran.r-project.org"))
options(dplyr.summarise.inform = FALSE)

pacman::p_load(
  # CRAN
  tidyverse, data.table, glue, here, fs, stringr, future, patchwork,
  cowplot, ggplot2, RColorBrewer, ggrepel, Matrix, matrixStats, scales,
  # Bioconductor
  BiocManager, SingleCellExperiment, scater, scran, scDblFinder, scuttle,
  # Seurat
  Seurat, SeuratObject,
  # Integration
  harmony,
  # Annotation
  SingleR, celldex,
  # DE
  edgeR, DESeq2, limma,
  # Pathways
  clusterProfiler, fgsea, msigdbr, enrichplot, DOSE,
  # I/O
  qs
)

  BiocManager::install("glmGamPoi", ask = FALSE, update = FALSE)
if (!requireNamespace("glmGamPoi", quietly = TRUE)) {
}

set.seed(42)

# Parallel & BLAS
workers <- min(8L, parallel::detectCores())
options(future.globals.maxSize = 110 * 1024^3)
future::plan(future::multisession, workers = workers)
Sys.setenv(
  OMP_NUM_THREADS = "1",
  MKL_NUM_THREADS = "1",
  OPENBLAS_NUM_THREADS = "1",
  VECLIB_MAXIMUM_THREADS = "1",
  NUMEXPR_NUM_THREADS = "1"
)

# ---------------------- User parameters ------------------------------------
params <- list(
  h5_dir = "/home/jacky/scRNASeq_Steffi/data/matrices",
  project_dir = "/home/jacky/scRNASeq_Steffi/analysis_outputs",
  h5_pattern = "*.h5",
  # filename parser expects DONOR_CONDITION_*.h5 (first two '_' tokens)
  min_features = NA, max_features = NA, min_counts = NA, max_mt_pct = NA,
  run_doublets = TRUE,
  use_sctransform = TRUE,
  sct_vars_to_regress = c("percent.mt"),         # lean regression for SCT v2
  integration_method = "harmony",                # "harmony" | "seurat_sct" | "none"
  harmony_group_var = "donor",
  n_pcs = 50,
  resolutions = c(0.1, 0.2, 0.5, 0.8),
  de_contrast = c("condition", "il18", "wo"),
  msigdb_collection = "H",                       # Hallmark in msigdbr >=10
  save_intermediates = TRUE
)

# ---------------------- Paths & logging ------------------------------------
proj <- fs::path_norm(params$project_dir)
dir_create <- function(p) if (!fs::dir_exists(p)) fs::dir_create(p, recurse = TRUE)
dirs <- list(
  root = proj,
  rds = fs::path(proj, "rds"),
  plots = fs::path(proj, "plots"),
  plots_qc = fs::path(proj, "plots", "qc"),
  plots_dim = fs::path(proj, "plots", "dimensionality"),
  plots_cluster = fs::path(proj, "plots", "clustering"),
  plots_markers = fs::path(proj, "plots", "markers"),
  plots_annot = fs::path(proj, "plots", "annotation"),
  plots_de = fs::path(proj, "plots", "dge"),
  plots_pea = fs::path(proj, "plots", "pea"),
  tables = fs::path(proj, "tables"),
  logs = fs::path(proj, "logs"),
  text = fs::path(proj, "text")
)
purrr::walk(dirs, dir_create)

log_msg <- function(...) {
  cat(glue::glue(...), "\n")
  write(glue::glue("{Sys.time()} | ", ...), file = fs::path(dirs$logs, "pipeline.log"), append = TRUE)
}

# ---------------------- Load data (10x HDF5) --------------------------------
h5_files <- fs::dir_ls(params$h5_dir, glob = params$h5_pattern, type = "file")
stopifnot("No .h5 files found in h5_dir" = length(h5_files) > 0)
log_msg("Found {length(h5_files)} H5 files")

parse_info <- function(fn) {
  bn <- basename(fn)
  core <- tools::file_path_sans_ext(bn)
  parts <- strsplit(core, "_", fixed = TRUE)[[1]]
  donor <- dplyr::na_if(parts[1], "")
  condition <- if (length(parts) >= 2) parts[2] else NA_character_
  tibble::tibble(file = fn, donor = donor, condition = condition, sample = core)
}
info <- purrr::map_dfr(h5_files, parse_info)
stopifnot("Filename parser didn’t extract donor/condition; expected DONOR_CONDITION_*.h5" =
            all(!is.na(info$donor)) & all(!is.na(info$condition)))
data.table::fwrite(info, fs::path(dirs$tables, "samples_parsed.csv"))

objs <- vector("list", nrow(info))
for (i in seq_len(nrow(info))) {
  fn <- info$file[i]
  log_msg("Reading {basename(fn)}")
  mat <- Seurat::Read10X_h5(fn)
  so  <- CreateSeuratObject(counts = mat, project = info$donor[i], min.cells = 3, min.features = 200)
  so$sample    <- info$sample[i]
  so$donor     <- info$donor[i]
  so$condition <- info$condition[i]
  colnames(so) <- paste(info$sample[i], colnames(so), sep = "_")  # stable barcodes
  objs[[i]]    <- so
}
stopifnot(length(objs) > 0)
obj <- Reduce(function(a, b) merge(a, b, add.cell.ids = NULL), objs)

# --- S4/Seurat v5 sanity
stopifnot(isS4(obj), inherits(obj, "Seurat"))

# ---------------------- QC metrics & filtering ------------------------------
obj[["percent.mt"]]   <- PercentageFeatureSet(obj, pattern = "^MT-")
obj[["percent.ribo"]] <- PercentageFeatureSet(obj, pattern = "^RP[SL]")

qc_by_sample <- obj@meta.data %>%
  dplyr::select(nCount_RNA, nFeature_RNA, percent.mt, percent.ribo, sample) %>%
  dplyr::mutate(across(c(nCount_RNA, nFeature_RNA), ~log10(.x + 1))) %>%
  dplyr::group_by(sample) %>%
  dplyr::summarize(
    min_features_auto = { x <- .$nFeature_RNA; thr <- min(x[!scater::isOutlier(x, nmads = 3, type = "lower",  log = FALSE)], na.rm = TRUE); ifelse(is.finite(thr), 10^thr - 1, 200) },
    max_features_auto = { x <- .$nFeature_RNA; thr <- max(x[!scater::isOutlier(x, nmads = 3, type = "higher", log = FALSE)], na.rm = TRUE); ifelse(is.finite(thr), 10^thr - 1, 8000) },
    min_counts_auto   = { x <- .$nCount_RNA;   thr <- min(x[!scater::isOutlier(x, nmads = 3, type = "lower",  log = FALSE)], na.rm = TRUE); ifelse(is.finite(thr), 10^thr - 1, 500) },
    max_mt_auto       = { x <- .$percent.mt;   thr <- max(x[!scater::isOutlier(x, nmads = 3, type = "higher", log = FALSE)], na.rm = TRUE); ifelse(is.finite(thr), thr, 20) }
  )

thr <- list(
  min_features = ifelse(is.na(params$min_features), median(qc_by_sample$min_features_auto), params$min_features),
  max_features = ifelse(is.na(params$max_features), median(qc_by_sample$max_features_auto), params$max_features),
  min_counts   = ifelse(is.na(params$min_counts),   median(qc_by_sample$min_counts_auto),   params$min_counts),
  max_mt_pct   = ifelse(is.na(params$max_mt_pct),   median(qc_by_sample$max_mt_auto),       params$max_mt_pct)
)

writeLines(c(
  "QC thresholds:",
  paste("  min_features:", thr$min_features),
  paste("  max_features:", thr$max_features),
  paste("  min_counts:",   thr$min_counts),
  paste("  max_mt_pct:",   thr$max_mt_pct)
), con = fs::path(dirs$text, "qc_thresholds.txt"))

data.table::fwrite(qc_by_sample, fs::path(dirs$tables, "qc_by_sample_auto.csv"))

keep <- with(obj@meta.data,
             nFeature_RNA >= thr$min_features &
             nFeature_RNA <= thr$max_features &
             nCount_RNA   >= thr$min_counts   &
             percent.mt   <= thr$max_mt_pct)
obj <- subset(obj, cells = colnames(obj)[keep])
stopifnot("All cells filtered out by QC" = ncol(obj) > 0)

png(fs::path(dirs$plots_qc, "qc_violin.png"), width = 1500, height = 900, res = 150)
print(VlnPlot(obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.ribo"),
              group.by = "sample", pt.size = 0.1, ncol = 2))
dev.off()

png(fs::path(dirs$plots_qc, "qc_scatter.png"), width = 1200, height = 800, res = 150)
print(FeatureScatter(obj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") +
        FeatureScatter(obj, "nCount_RNA", "percent.mt"))
dev.off()

if (params$save_intermediates) saveRDS(obj, fs::path(dirs$rds, "01_qc_filtered.rds"))

# ---------------------- Doublet detection (v5-layer-safe) -------------------
if (isTRUE(params$run_doublets)) {
  DefaultAssay(obj) <- ifelse("RNA" %in% names(obj@assays), "RNA", DefaultAssay(obj))
  aa <- DefaultAssay(obj)
  counts_mat <- SeuratObject::LayerData(obj[[aa]], layer = "counts")   # v5 accessor

  sce <- SingleCellExperiment::SingleCellExperiment(
    assays  = list(counts = counts_mat),
    colData = obj@meta.data
  )

  set.seed(42)
  sce <- scDblFinder::scDblFinder(sce, samples = sce$sample)

  obj$doublet_score <- sce$scDblFinder.score
  obj$doublet_class <- sce$scDblFinder.class

  png(fs::path(dirs$plots_qc, "doublet_scores.png"), width = 1200, height = 800, res = 150)
  print(ggplot(as.data.frame(SummarizedExperiment::colData(sce)),
               aes(x = sample, y = scDblFinder.score)) +
          geom_boxplot() + theme_bw() + labs(y = "Doublet score"))
  dev.off()

  obj <- subset(obj, cells = colnames(obj)[obj$doublet_class == "singlet"])
  stopifnot("All cells removed as doublets" = ncol(obj) > 0)
  if (params$save_intermediates) saveRDS(obj, fs::path(dirs$rds, "02_singlets.rds"))
}

# ---------------------- SCTransform v2 (memory-safe) -----------------------
if (isTRUE(params$use_sctransform)) {
  prev_plan <- future::plan(); on.exit(future::plan(prev_plan), add = TRUE)
  future::plan(future::sequential)   # avoid huge globals export

  obj <- SCTransform(
    obj,
    vst.flavor = "v2",
    variable.features.n = 3000,
    vars.to.regress = params$sct_vars_to_regress,  # lean: percent.mt
    method = "glmGamPoi",
    verbose = FALSE
  )

  future::plan(prev_plan)
} else {
  obj <- NormalizeData(obj)
  obj <- FindVariableFeatures(obj, selection.method = "vst", nfeatures = 3000)
  obj <- ScaleData(obj)
}
if (params$save_intermediates) saveRDS(obj, fs::path(dirs$rds, "03_normalized.rds"))

# --- S4/v5 sanity (assay layers exist)
aa <- DefaultAssay(obj)
stopifnot(isS4(obj), inherits(obj, "Seurat"))
stopifnot(aa %in% names(obj@assays))
stopifnot(all(c("counts","data") %in% names(obj[[aa]]@layers)))

# ---------------------- Dimensionality reduction & integration -------------
obj <- RunPCA(obj, npcs = params$n_pcs, verbose = FALSE)
stopifnot("pca" %in% names(obj@reductions))

png(fs::path(dirs$plots_dim, "pca_elbow.png"), width = 1000, height = 700, res = 150)
print(ElbowPlot(obj, ndims = params$n_pcs))
dev.off()

if (identical(params$integration_method, "harmony")) {
  stopifnot("Batch column for Harmony not found in meta.data" =
              params$harmony_group_var %in% colnames(obj@meta.data))

  obj <- RunHarmony(
    object        = obj,
    group.by.vars = params$harmony_group_var,
    reduction     = "pca",
    dims.use      = 1:params$n_pcs,
    verbose       = TRUE
  )
  stopifnot("harmony" %in% names(obj@reductions))

  obj <- RunUMAP(obj, reduction = "harmony", dims = 1:30)
  obj <- FindNeighbors(obj, reduction = "harmony", dims = 1:30, k.param = 30)

} else if (identical(params$integration_method, "seurat_sct")) {
  objs_list <- SplitObject(obj, split.by = params$harmony_group_var)
  objs_list <- lapply(objs_list, function(x) SCTransform(x, vst.flavor = "v2", method = "glmGamPoi", verbose = FALSE))
  features  <- SelectIntegrationFeatures(objs_list, nfeatures = 3000)
  objs_list <- PrepSCTIntegration(objs_list, anchor.features = features)
  anchors   <- FindIntegrationAnchors(objs_list, normalization.method = "SCT", anchor.features = features)
  obj       <- IntegrateData(anchorset = anchors, normalization.method = "SCT")
  obj <- RunPCA(obj, npcs = params$n_pcs, verbose = FALSE)
  obj <- RunUMAP(obj, dims = 1:30)
  obj <- FindNeighbors(obj, dims = 1:30, k.param = 30)

} else {
  obj <- RunUMAP(obj, dims = 1:30)
  obj <- FindNeighbors(obj, dims = 1:30, k.param = 30)
}
stopifnot("umap" %in% names(obj@reductions))

# Leiden clustering (algorithm = 4)
for (res in params$resolutions) {
  obj <- FindClusters(obj, resolution = res, algorithm = 4)
}
stopifnot(any(grepl("^SCT_snn_res\\.", colnames(obj@meta.data))))
if (params$save_intermediates) saveRDS(obj, fs::path(dirs$rds, "04_integrated_clustered.rds"))

png(fs::path(dirs$plots_cluster, "umap_by_cluster.png"), width = 1400, height = 900, res = 150)
print(DimPlot(obj, reduction = "umap", group.by = paste0("SCT_snn_res.", tail(params$resolutions, 1))) +
        ggtitle("UMAP (highest resolution)"))
dev.off()

png(fs::path(dirs$plots_cluster, "umap_by_sample.png"), width = 1400, height = 900, res = 150)
print(DimPlot(obj, reduction = "umap", group.by = "sample"))
dev.off()

# ---------------------- Annotation (SingleR, v5-layer-safe) ----------------
hpca <- celldex::HumanPrimaryCellAtlasData()
active_assay <- DefaultAssay(obj)
sce <- as.SingleCellExperiment(obj, assay = active_assay)

# Provide counts/logcounts explicitly using LayerData()
SummarizedExperiment::assay(sce, "counts")    <- SeuratObject::LayerData(obj[[active_assay]], layer = "counts")
SummarizedExperiment::assay(sce, "logcounts") <- SeuratObject::LayerData(obj[[active_assay]], layer = "data")
stopifnot(all(c("counts","logcounts") %in% SummarizedExperiment::assayNames(sce)))

pred <- SingleR(test = sce, ref = hpca, labels = hpca$label.main)
obj$SingleR_label <- pred$labels

clus <- Idents(obj)
tab  <- table(Cluster = clus, Label = obj$SingleR_label)
data.table::fwrite(as.data.frame(tab), fs::path(dirs$tables, "annotation_confusion.csv"))

png(fs::path(dirs$plots_annot, "umap_by_SingleR.png"), width = 1400, height = 900, res = 150)
print(DimPlot(obj, reduction = "umap", group.by = "SingleR_label", label = TRUE, repel = TRUE))
dev.off()

if (params$save_intermediates) saveRDS(obj, fs::path(dirs$rds, "05_annotated.rds"))

# ---------------------- Markers per cluster --------------------------------
Idents(obj) <- paste0("SCT_snn_res.", tail(params$resolutions, 1))
markers <- FindAllMarkers(obj, only.pos = TRUE, test.use = "wilcox", min.pct = 0.1, logfc.threshold = 0.25)
data.table::fwrite(markers, fs::path(dirs$tables, "markers_wilcox.csv"))

if (nrow(markers) > 0) {
  top_markers <- markers %>% group_by(cluster) %>% slice_max(order_by = avg_log2FC, n = 10)
  writeLines(capture.output(print(top_markers)), fs::path(dirs$text, "top10_markers.txt"))
}

# ---------------------- Pseudobulk with scuttle (sparse-safe) --------------
sce_pb <- as.SingleCellExperiment(obj, assay = DefaultAssay(obj))
SummarizedExperiment::assay(sce_pb, "counts") <- SeuratObject::LayerData(obj[[DefaultAssay(obj)]], layer = "counts")

meta  <- as.data.frame(SummarizedExperiment::colData(sce_pb)) %>%
  dplyr::mutate(cluster = as.character(Seurat::Idents(obj)))

ids <- interaction(meta$sample, meta$cluster, drop = TRUE)
pb <- scuttle::sumCountsAcrossCells(sce_pb, ids = ids)
pb_counts <- SummarizedExperiment::assay(pb, "counts")
colnames(pb_counts) <- levels(ids)

samples <- as.data.frame(SummarizedExperiment::colData(pb)) %>%
  tibble::rownames_to_column("sample_cluster") %>%
  dplyr::rename(sample = sample, cluster = cluster)

meta_simple <- meta %>% distinct(sample, donor, condition)
samples <- dplyr::left_join(samples, meta_simple, by = "sample")

dge <- edgeR::DGEList(pb_counts)
keep_genes <- edgeR::filterByExpr(dge, group = samples$condition)
dge <- dge[keep_genes, , keep.lib.sizes = FALSE]
dge <- edgeR::calcNormFactors(dge, method = "TMM")

design <- model.matrix(~ 0 + condition + donor, data = samples)
colnames(design) <- gsub("^condition", "", colnames(design))

v    <- limma::voom(dge, design, plot = FALSE)
fit  <- limma::lmFit(v, design)
cn   <- glue::glue("{params$de_contrast[2]}-{params$de_contrast[3]}")
contr <- limma::makeContrasts(contrasts = cn, levels = colnames(design))
fit2 <- limma::contrasts.fit(fit, contr)
fit2 <- limma::eBayes(fit2)
res  <- limma::topTable(fit2, number = Inf, sort.by = "P")
res$gene_symbol <- rownames(res)
data.table::fwrite(res, fs::path(dirs$tables, "pseudobulk_edgeR_allclusters.csv"))

# ---------------------- Pathway enrichment ---------------------------------
gene_sets_df <- msigdbr(species = "Homo sapiens", collection = params$msigdb_collection)
TERM2GENE   <- gene_sets_df %>% dplyr::select(gs_name, gene_symbol)

res$diffexpressed <- dplyr::case_when(
  res$logFC > 0 & res$adj.P.Val < 0.05 ~ "UP",
  res$logFC < 0 & res$adj.P.Val < 0.05 ~ "DOWN",
  TRUE ~ "NO"
)
res_sig <- res %>% dplyr::filter(diffexpressed != "NO")

deg_list <- split(res_sig, res_sig$diffexpressed)
ora <- lapply(deg_list, function(df) clusterProfiler::enricher(
  gene = df$gene_symbol, TERM2GENE = TERM2GENE, qvalueCutoff = 0.2))
saveRDS(ora, fs::path(dirs$rds, "ora_results.rds"))

rankings <- sign(res$logFC) * (-log10(res$P.Value))
names(rankings) <- res$gene_symbol
rankings[!is.finite(rankings)] <- 0
rankings <- sort(rankings, decreasing = TRUE)
gene_sets <- split(gene_sets_df$gene_symbol, gene_sets_df$gs_name)

fg <- fgsea::fgsea(pathways = gene_sets, stats = rankings, scoreType = "std",
                   minSize = 10, maxSize = 500, nproc = 1)
data.table::fwrite(as.data.frame(fg), fs::path(dirs$tables, "fgsea_results.csv"))

if (!is.null(ora$UP) && length(ora$UP) > 0) {
  png(fs::path(dirs$plots_pea, "ora_barplot_up.png"), width = 1200, height = 900, res = 150)
  print(enrichplot::barplot(ora$UP, showCategory = 15, title = "ORA – Upregulated"))
  dev.off()
}
if (!is.null(ora$DOWN) && length(ora$DOWN) > 0) {
  png(fs::path(dirs$plots_pea, "ora_barplot_down.png"), width = 1200, height = 900, res = 150)
  print(enrichplot::barplot(ora$DOWN, showCategory = 15, title = "ORA – Downregulated"))
  dev.off()
}

fg_df <- as.data.frame(fg) %>% arrange(padj) %>% head(30)
png(fs::path(dirs$plots_pea, "fgsea_top30.png"), width = 1500, height = 800, res = 150)
print(ggplot(fg_df, aes(x = NES, y = reorder(pathway, NES), size = size, color = -log10(padj))) +
        geom_point() + theme_minimal() + labs(x = "NES", y = "Pathway"))
dev.off()

# ---------------------- Save & run summary ---------------------------------
saveRDS(obj, fs::path(dirs$rds, "final_seurat_object.rds"))
writeLines(capture.output(sessionInfo()), fs::path(dirs$text, "sessionInfo.txt"))

summary_lines <- c(
  "Accessible run summary",
  paste0("Cells after QC+singlets: ", ncol(obj)),
  paste0("Genes: ", nrow(obj)),
  paste0("Clusters (highest res ", tail(params$resolutions, 1), "): ",
         length(unique(as.character(Idents(obj))))),
  paste0("Annotation labels: ",
         paste(head(sort(table(obj$SingleR_label), decreasing = TRUE), 10) %>% names(), collapse = ", "))
)
writeLines(summary_lines, fs::path(dirs$text, "run_summary.txt"))

log_msg("Pipeline completed successfully.")
# =========================================================================
