install.packages("ggplot2")
install.packages("BiocManager")
BiocManager::install("limma")
BiocManager::install("edgeR")

library(limma)
library(ggplot2)
library(edgeR)
# ==============================================================================
# 1. Load expression data sets


general_expr <- read.delim("C:/R_Data/RNASeq_Control_Gene_Expression_Data.tsv",
                           header = TRUE,
                           row.names = 1,
                           check.names = FALSE)

tumor_expr <- read.delim("C:/R_Data/RNASeq_Tumor_Gene_Expression_Data.tsv",
                        header = TRUE,
                        row.names = 1,
                        check.names = FALSE)

# ==============================================================================
# 2. Make sure all expression values are numeric

general_expr <- data.frame(lapply(general_expr, as.numeric),
                           row.names = rownames(general_expr),
                           check.names = FALSE)

tumor_expr <- data.frame(lapply(tumor_expr, as.numeric),
                        row.names = rownames(tumor_expr),
                        check.names = FALSE)

# ==============================================================================
# 3. Making unique groups before combining to run differential gene expression analysis

colnames(general_expr) <- paste0("Control_", colnames(general_expr))
colnames(tumor_expr) <- paste0("Tumor_", colnames(tumor_expr))

# ==============================================================================
# 4. Combine the two labelled groups

stopifnot(all(rownames(general_expr) == rownames(tumor_expr)))
expr_combined <- cbind(general_expr, tumor_expr)

# ==============================================================================
# 5. Create group information

group <- c(
  rep("Control", ncol(general_expr)),
  rep("Tumor", ncol(tumor_expr))
)

group <- factor(group, levels = c("Control", "Tumor"))

meta_combined <- data.frame(
  sample_id = colnames(expr_combined),
  group = group
)

table(meta_combined$group)

# ==============================================================================
# 6. Log-transform normalized expression

log_expr <- log2(expr_combined + 1)

# ==============================================================================
# 7. Build design matrix

design <- model.matrix(~ 0 + group, data = meta_combined)
colnames(design) <- levels(group)

design[1:5, ]

# ==============================================================================
# 8. Run limma-trend

fit <- lmFit(log_expr, design)

contrast.matrix <- makeContrasts(
  Tumor_vs_General = Tumor - Control,
  levels = design
)

fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2, trend = TRUE)

# ==============================================================================
# 9. Extract results

res <- topTable(fit2,
                coef = "Tumor_vs_General",
                number = Inf,
                adjust.method = "BH",
                sort.by = "P")

res$gene_id <- rownames(res)

# ==============================================================================
# 10. Classify UP / DOWN / NS genes

logFC_cut <- 1
adjP_cut <- 0.05

res$regulation <- "NS"

res$regulation[res$adj.P.Val < adjP_cut & res$logFC >  logFC_cut] <- "UP"
res$regulation[res$adj.P.Val < adjP_cut & res$logFC < -logFC_cut] <- "DOWN"

table(res$regulation)

up_genes <- res[res$regulation == "UP", ]
down_genes <- res[res$regulation == "DOWN", ]

# ==============================================================================
# 11. Volcano plot

volcano_df <- res

# Incase p-values = 0 are present
volcano_df$negLogAdjP <- -log10(pmax(volcano_df$adj.P.Val, 1e-300))

ggplot(volcano_df, aes(x = logFC, y = negLogAdjP, color = regulation)) +
  geom_point(alpha = 0.7, size = 1.5) +
  scale_color_manual(values = c("UP" = "red", "DOWN" = "blue", "NS" = "grey")) +
  geom_vline(xintercept = c(-logFC_cut, logFC_cut),
             linetype = "dashed", color = "black") +
  geom_hline(yintercept = -log10(adjP_cut),
             linetype = "dashed", color = "black") +
  theme_minimal(base_size = 14) +
  labs(
    title = "Volcano Plot: Tumor versus Control",
    x = "log2 Fold Change",
    y = "-log10 Adjusted P-value",
    color = "Regulation"
  )

# Right side red genes = upregulated in Tumor (wrt Control)
# Left side blue genes = downregulated in Tumor (wrt Control)

# ==============================================================================
# 12. Save output files

outdir <- "Tumor_vs_Control_limma_trend_results"

if (!dir.exists(outdir)) {
  dir.create(outdir, recursive = TRUE)
}

write.csv(res,
          file.path(outdir, "Tumor_vs_Control_full.csv"),
          row.names = FALSE)

write.csv(up_genes,
          file.path(outdir, "Tumor_vs_Control_UP.csv"),
          row.names = FALSE)

write.csv(down_genes,
          file.path(outdir, "Tumor_vs_Control_DOWN.csv"),
          row.names = FALSE)
# ==============================================================================
# DAVID input files
write.table(up_genes$gene_id,
            file.path(outdir, "DAVID_Tumor_vs_Control_UP.txt"),
            quote = FALSE,
            row.names = FALSE,
            col.names = FALSE)

write.table(down_genes$gene_id,
            file.path(outdir, "DAVID_Tumor_vs_Control_DOWN.txt"),
            quote = FALSE,
            row.names = FALSE,
            col.names = FALSE)

write.table(rownames(log_expr),
            file.path(outdir, "DAVID_background_Tumor_vs_Control.txt"),
            quote = FALSE,
            row.names = FALSE,
            col.names = FALSE)
# ==============================================================================

