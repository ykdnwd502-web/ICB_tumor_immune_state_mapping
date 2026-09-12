############################################################
## 09A-06. Generate hg38 inferCNV gene order
##
## Source: TxDb.Hsapiens.UCSC.hg38.knownGene + org.Hs.eg.db.
## No package installation and no downstream script patching are performed.
############################################################

options(stringsAsFactors = FALSE)
SCRIPT_ID <- "09A_06_generate_inferCNV_gene_order"

project_dir <- Sys.getenv("ICB_PROJECT_DIR")
if (!nzchar(project_dir)) project_dir <- getwd()
project_dir <- normalizePath(project_dir, winslash = "/", mustWork = TRUE)

required_pkgs <- c(
  "Seurat", "SeuratObject", "readr", "dplyr", "tibble",
  "AnnotationDbi", "org.Hs.eg.db", "TxDb.Hsapiens.UCSC.hg38.knownGene",
  "GenomicFeatures", "GenomeInfoDb", "BiocGenerics"
)
missing_pkgs <- required_pkgs[!vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_pkgs) > 0L) {
  stop(
    "Missing required annotation package(s): ", paste(missing_pkgs, collapse = ", "),
    ". Install/restore them before the full clean run; this script never auto-installs.",
    call. = FALSE
  )
}

out_table_dir <- file.path(project_dir, "results", "tables", "GSE244983", "doublet_CNV_sensitivity")
out_obj_dir <- file.path(project_dir, "results", "objects", "GSE244983", "doublet_CNV_sensitivity")
input_rds <- file.path(out_obj_dir, "GSE244983_singlets_after_pANN_doublet_filter.rds")
gene_order_header <- file.path(out_table_dir, "D06_inferCNV_gene_order_hg38_header.tsv")
gene_order_noheader <- file.path(out_table_dir, "D06_inferCNV_gene_order_hg38_noheader.txt")
write_csv <- function(x, name) readr::write_csv(x, file.path(out_table_dir, name), na = "")

if (!file.exists(input_rds)) stop("09A-03 singlet object missing: ", input_rds, call. = FALSE)
obj <- readRDS(input_rds)
Seurat::DefaultAssay(obj) <- if ("RNA" %in% names(obj@assays)) "RNA" else Seurat::DefaultAssay(obj)
seurat_genes <- unique(as.character(rownames(obj)))

map_symbol <- AnnotationDbi::select(
  org.Hs.eg.db::org.Hs.eg.db,
  keys = seurat_genes,
  keytype = "SYMBOL",
  columns = c("SYMBOL", "ENTREZID")
)
map_symbol <- as.data.frame(map_symbol, stringsAsFactors = FALSE)
map_symbol <- map_symbol[!is.na(map_symbol$SYMBOL) & !is.na(map_symbol$ENTREZID), , drop = FALSE]
map_symbol$input_gene <- map_symbol$SYMBOL
map_symbol$mapping_source <- "SYMBOL"

remaining <- setdiff(seurat_genes, unique(map_symbol$input_gene))
map_alias <- data.frame()
if (length(remaining) > 0L) {
  map_alias <- AnnotationDbi::select(
    org.Hs.eg.db::org.Hs.eg.db,
    keys = remaining,
    keytype = "ALIAS",
    columns = c("ALIAS", "SYMBOL", "ENTREZID")
  )
  map_alias <- as.data.frame(map_alias, stringsAsFactors = FALSE)
  map_alias <- map_alias[!is.na(map_alias$ALIAS) & !is.na(map_alias$ENTREZID), , drop = FALSE]
  map_alias$input_gene <- map_alias$ALIAS
  map_alias$mapping_source <- "ALIAS"
}

map_all <- rbind(
  map_symbol[, c("input_gene", "SYMBOL", "ENTREZID", "mapping_source"), drop = FALSE],
  if (nrow(map_alias)) map_alias[, c("input_gene", "SYMBOL", "ENTREZID", "mapping_source"), drop = FALSE] else map_symbol[0, c("input_gene", "SYMBOL", "ENTREZID", "mapping_source"), drop = FALSE]
)
map_all <- dplyr::distinct(map_all, input_gene, SYMBOL, ENTREZID, .keep_all = TRUE)

write_csv(map_all, "D06_gene_symbol_to_entrez_mapping_details.csv")
write_csv(
  data.frame(unmapped_gene = setdiff(seurat_genes, unique(map_all$input_gene))),
  "D06_unmapped_genes.csv"
)

txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene::TxDb.Hsapiens.UCSC.hg38.knownGene
gr <- GenomicFeatures::genes(txdb)
gene_pos <- data.frame(
  ENTREZID = names(gr),
  chr = as.character(GenomeInfoDb::seqnames(gr)),
  start = as.integer(BiocGenerics::start(gr)),
  stop = as.integer(BiocGenerics::end(gr)),
  stringsAsFactors = FALSE
)
gene_pos <- gene_pos[!is.na(gene_pos$ENTREZID) & !is.na(gene_pos$chr), , drop = FALSE]

standard_chr <- paste0("chr", c(1:22, "X", "Y"))
gene_order <- dplyr::inner_join(map_all, gene_pos, by = "ENTREZID")
gene_order <- gene_order[gene_order$chr %in% standard_chr, , drop = FALSE]
gene_order$chr_order <- match(gene_order$chr, standard_chr)
gene_order$gene <- gene_order$input_gene
gene_order <- dplyr::arrange(gene_order, chr_order, start, stop, gene)
gene_order <- dplyr::group_by(gene_order, gene)
gene_order <- dplyr::slice_head(gene_order, n = 1L)
gene_order <- dplyr::ungroup(gene_order)
gene_order <- dplyr::select(gene_order, gene, chr, start, stop)
gene_order <- dplyr::distinct(gene_order, gene, .keep_all = TRUE)

if (nrow(gene_order) < 3000L) {
  stop("Generated inferCNV gene-order table has fewer than 3,000 genes.", call. = FALSE)
}

readr::write_tsv(gene_order, gene_order_header)
readr::write_tsv(gene_order, gene_order_noheader, col_names = FALSE)

audit <- data.frame(
  input_rds = input_rds,
  n_singlet_cells = ncol(obj),
  n_seurat_genes = length(seurat_genes),
  n_mapped_input_genes = length(unique(map_all$input_gene)),
  n_unmapped_input_genes = length(setdiff(seurat_genes, unique(map_all$input_gene))),
  n_gene_order_genes = nrow(gene_order),
  fraction_seurat_genes_in_gene_order = nrow(gene_order) / length(seurat_genes),
  genome_build = "hg38",
  source = "TxDb.Hsapiens.UCSC.hg38.knownGene + org.Hs.eg.db",
  gene_order_noheader = gene_order_noheader,
  stringsAsFactors = FALSE
)
write_csv(audit, "D06_inferCNV_gene_order_audit.csv")
writeLines(capture.output(sessionInfo()), file.path(out_table_dir, paste0(SCRIPT_ID, "_sessionInfo.txt")))
message("09A-06 complete: ", nrow(gene_order), " genes in hg38 inferCNV gene order.")
