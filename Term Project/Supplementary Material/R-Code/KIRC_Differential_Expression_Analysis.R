library(limma)
library(ggplot2)

# ==============================================================================
# 1. Load expression datasets


general_expr <- read.delim("D:/Semester_Files_(Currently-V)/AIBD_Final_Project/Central_cluster_data/Central_cluster_tumor_Gene_Expression_Data.tsv",
                           header = TRUE,
                           row.names = 1,
                           check.names = FALSE)

kirc_expr <- read.delim("D:/Semester_Files_(Currently-V)/AIBD_Final_Project/KIRC_Cluster_Analysis/KIRC_Gene_Expression_Data.tsv",
                        header = TRUE,
                        row.names = 1,
                        check.names = FALSE)

# ==============================================================================
# 2. Make sure all expression values are numeric

general_expr <- data.frame(lapply(general_expr, as.numeric),
                           row.names = rownames(general_expr),
                           check.names = FALSE)

kirc_expr <- data.frame(lapply(kirc_expr, as.numeric),
                        row.names = rownames(kirc_expr),
                        check.names = FALSE)

# ==============================================================================
# 3. Making unique groups before combining to run differential gene expression analysis

colnames(general_expr) <- paste0("General_", colnames(general_expr))
colnames(kirc_expr) <- paste0("KIRC_", colnames(kirc_expr))

# ==============================================================================
# 4. Combine the two labelled groups

stopifnot(all(rownames(general_expr) == rownames(kirc_expr)))
expr_combined <- cbind(general_expr, kirc_expr)

# ==============================================================================
# 5. Create group information

group <- c(
  rep("GeneralCluster", ncol(general_expr)),
  rep("KIRCCluster", ncol(kirc_expr))
)

group <- factor(group, levels = c("GeneralCluster", "KIRCCluster"))

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
  KIRC_vs_General = KIRCCluster - GeneralCluster,
  levels = design
)

fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2, trend = TRUE)

# ==============================================================================
# 9. Extract results

res <- topTable(fit2,
                coef = "KIRC_vs_General",
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
    title = "Volcano Plot: KIRC vs Central Cancer Cluster",
    x = "log2 Fold Change",
    y = "-log10 Adjusted P-value",
    color = "Regulation"
  )

# Right side red genes = upregulated in KIRC
# Left side blue genes = downregulated in KIRC

# ==============================================================================
# 12. Save output files

outdir <- "D:/Semester_Files_(Currently-V)/AIBD_Final_Project/KIRC_vs_General_limma_trend_results"

# ==============================================================================
if (!dir.exists(outdir)) {
  dir.create(outdir, recursive = TRUE)
}

write.csv(res,
          file.path(outdir, "KIRC_vs_General_full.csv"),
          row.names = FALSE)

write.csv(up_genes,
          file.path(outdir, "KIRC_vs_General_UP.csv"),
          row.names = FALSE)

write.csv(down_genes,
          file.path(outdir, "KIRC_vs_General_DOWN.csv"),
          row.names = FALSE)
# ==============================================================================
# DAVID input files
write.table(up_genes$gene_id,
            file.path(outdir, "DAVID_KIRC_vs_General_UP.txt"),
            quote = FALSE,
            row.names = FALSE,
            col.names = FALSE)

write.table(down_genes$gene_id,
            file.path(outdir, "DAVID_KIRC_vs_General_DOWN.txt"),
            quote = FALSE,
            row.names = FALSE,
            col.names = FALSE)

write.table(rownames(log_expr),
            file.path(outdir, "DAVID_background_KIRC_vs_General.txt"),
            quote = FALSE,
            row.names = FALSE,
            col.names = FALSE)
# ==============================================================================