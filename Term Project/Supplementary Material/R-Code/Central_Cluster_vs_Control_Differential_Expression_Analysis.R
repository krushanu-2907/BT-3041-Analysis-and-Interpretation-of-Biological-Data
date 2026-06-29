library(limma)
library(ggplot2)

# ==============================================================================
# 1. Load expression datasets

central_control_expr <- read.delim("D:/Semester_Files_(Currently-V)/AIBD_Final_Project/Central_cluster_data/Central_Cluster_Control_Gene_Expression_Data.tsv",
                                   header = TRUE,
                                   row.names = 1,
                                   check.names = FALSE)

central_tumor_expr <- read.delim("D:/Semester_Files_(Currently-V)/AIBD_Final_Project/Central_cluster_data/Central_cluster_tumor_Gene_Expression_Data.tsv",
                                 header = TRUE,
                                 row.names = 1,
                                 check.names = FALSE)

# ==============================================================================
# 2. Make sure all expression values are numeric

central_control_expr <- data.frame(lapply(central_control_expr, as.numeric),
                                   row.names = rownames(central_control_expr),
                                   check.names = FALSE)

central_tumor_expr <- data.frame(lapply(central_tumor_expr, as.numeric),
                                 row.names = rownames(central_tumor_expr),
                                 check.names = FALSE)

# ==============================================================================
# 3. Making unique groups before combining to run differential gene expression analysis

colnames(central_control_expr) <- paste0("Central_Control_", colnames(central_control_expr))
colnames(central_tumor_expr) <- paste0("Central_Tumor_", colnames(central_tumor_expr))

# ==============================================================================
# 4. Combine the two labelled groups

stopifnot(all(rownames(central_control_expr) == rownames(central_tumor_expr)))
expr_combined <- cbind(central_control_expr, central_tumor_expr)

# ==============================================================================
# 5. Create group information

group <- c(
  rep("CentralControlCluster", ncol(central_control_expr)),
  rep("CentralTumorCluster", ncol(central_tumor_expr))
)

group <- factor(group, levels = c("CentralControlCluster", "CentralTumorCluster"))

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
  Central_Tumor_vs_Central_Control = CentralTumorCluster - CentralControlCluster,
  levels = design
)

fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2, trend = TRUE)

# ==============================================================================
# 9. Extract results

res <- topTable(fit2,
                coef = "Central_Tumor_vs_Central_Control",
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

# In case p-values = 0 are present
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
    title = "Volcano Plot: Central Tumor Cluster vs Central Control Cluster",
    x = "log2 Fold Change",
    y = "-log10 Adjusted P-value",
    color = "Regulation"
  )

# Right side red genes = upregulated in Central Tumor Cluster
# Left side blue genes = downregulated in Central Tumor Cluster

# ==============================================================================
# 12. Save output files

outdir <- "D:/Semester_Files_(Currently-V)/AIBD_Final_Project/Central_cluster_data/Central_Tumor_vs_Central_Control_limma_trend_results"

# ==============================================================================
if (!dir.exists(outdir)) {
  dir.create(outdir, recursive = TRUE)
}

write.csv(res,
          file.path(outdir, "Central_Tumor_vs_Central_Control_full.csv"),
          row.names = FALSE)

write.csv(up_genes,
          file.path(outdir, "Central_Tumor_vs_Central_Control_UP.csv"),
          row.names = FALSE)

write.csv(down_genes,
          file.path(outdir, "Central_Tumor_vs_Central_Control_DOWN.csv"),
          row.names = FALSE)

# ==============================================================================
# DAVID input files

write.table(up_genes$gene_id,
            file.path(outdir, "DAVID_Central_Tumor_vs_Central_Control_UP.txt"),
            quote = FALSE,
            row.names = FALSE,
            col.names = FALSE)

write.table(down_genes$gene_id,
            file.path(outdir, "DAVID_Central_Tumor_vs_Central_Control_DOWN.txt"),
            quote = FALSE,
            row.names = FALSE,
            col.names = FALSE)

write.table(rownames(log_expr),
            file.path(outdir, "DAVID_background_Central_Tumor_vs_Central_Control.txt"),
            quote = FALSE,
            row.names = FALSE,
            col.names = FALSE)

# ==============================================================================