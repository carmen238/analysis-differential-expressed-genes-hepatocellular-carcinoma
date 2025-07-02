# ===========================================
# Differential Gene Expression Analysis - GSE22058 (HCC)
# ===========================================

# INSTALL AND LOAD PACKAGES (only first time)
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(c("GEOquery", "clusterProfiler", "org.Hs.eg.db"))

# LOAD LIBRARIES
library(GEOquery)
library(clusterProfiler)
library(org.Hs.eg.db)
library(pheatmap)

# CREATE RESULTS FOLDERS

# Main folder
dirRes <- "Results/"
dataset <- "GEO_HCC"
dirDataset <- paste0(dirRes, dataset, "/")
if (!dir.exists(dirRes)) dir.create(dirRes)
if (!dir.exists(dirDataset)) dir.create(dirDataset)

# ===========================================
# STEP 1: Downloading data
# ===========================================
series <- "GSE22058"
tmp <- getGEO(GEO = series)
set <- tmp[["GSE22058-GPL6793_series_matrix.txt.gz"]]

pData <- phenoData(set)
metadata <- pData(set)
aData <- assayData(set)
matrix <- data.frame(aData$exprs)

# ===========================================
# STEP 2: Preparing data
# ===========================================
annotation <- fData(set)
geneSymbol <- annotation$GeneSymbol
matrix <- matrix[rownames(annotation), ]
matrix <- aggregate(matrix, list(geneSymbol), mean)
ind <- which(matrix$Group.1 == "")
matrix <- matrix[-ind, ]
rownames(matrix) <- matrix$Group.1
matrix <- matrix[, -1]

# ===========================================
# STEP 3: Extract tumor and normal samples
# ===========================================
metadata <- metadata[, c("geo_accession", "individual:ch1", "tissue:ch1")]
list <- split(metadata, metadata$`tissue:ch1`)
normal <- list$`adjacent liver non-tumor`
tumor  <- list$`liver tumor`

normal <- normal[!duplicated(normal$`individual:ch1`), ]
tumor  <- tumor[!duplicated(tumor$`individual:ch1`), ]

normal <- normal[order(normal$`individual:ch1`), "geo_accession"]
tumor  <- tumor[order(tumor$`individual:ch1`), "geo_accession"]

data <- matrix[, c(normal, tumor)]
dataN <- matrix[, normal]
dataC <- matrix[, tumor]

# ===========================================
# STEP 4: Export raw data
# ===========================================
write.table(data, paste0(dirDataset, "matrix.txt"), sep = "\t", col.names = NA)
write.table(normal, paste0(dirDataset, "normal.txt"), sep = "\t", col.names = FALSE, row.names = FALSE)
write.table(tumor, paste0(dirDataset, "tumor.txt"), sep = "\t", col.names = FALSE, row.names = FALSE)
write.table(metadata, paste0(dirDataset, "metadata.txt"), sep = "\t", col.names = TRUE, row.names = FALSE)

# ===========================================
# STEP 5: Preprocessing (log2, IQR, filtering)
# ===========================================
dataN <- log2(dataN + 1)
dataC <- log2(dataC + 1)
dataAll <- cbind(dataN, dataC)

# Remove genes with mean = 0 in both groups
mean_N <- rowMeans(dataN)
mean_C <- rowMeans(dataC)
keep_mean <- !(mean_N == 0 & mean_C == 0)
dataN <- dataN[keep_mean, ]
dataC <- dataC[keep_mean, ]
dataAll <- cbind(dataN, dataC)
threshold <- quantile(iqr, 0.10)
filtered_data <- dataAll[iqr > threshold, ]

# IQR filtering
iqr <- apply(dataAll, 1, IQR)
png("istogramma_iqr.png", width = 800, height = 600)
hist(iqr, breaks = 50, main = "Distribuzione IQR", xlab = "IQR", ylab="Frequency", col = "skyblue")
abline(v= threshold, lty = 2 , lwd = 4 , col="darkgoldenrod1")
dev.off()


# DEG analysis
N <- ncol(dataN)
M <- ncol(dataC)
filtered_dataN <- filtered_data[, 1:N]
filtered_dataC <- filtered_data[, (N + 1):(N + M)]

logFC <- rowMeans(filtered_dataC) - rowMeans(filtered_dataN)
pval <- apply(filtered_data, 1, function(x) t.test(x[1:N], x[(N+1):(N+M)], paired = TRUE)$p.value)
pval_adj <- p.adjust(pval, method = "fdr")

results <- data.frame(gene = rownames(filtered_data), logFC = logFC, pval = pval, adj_pval = pval_adj)
DEGs <- results[abs(results$logFC) > 1 & results$adj_pval < 0.05, ]

# Volcano plot
# ======== Volcano Plot 1: tutti i geni filtrati per IQR =========
png(paste0(dirDataset, "volcano_plot_1.png"), width = 800, height = 600)

plot(results$logFC, -log10(results$adj_pval),
     pch = 20, col = "grey",
     main = "Volcano Plot - Tutti i geni (post-IQR)",
     xlab = "log2 Fold Change",
     ylab = "-log10 Adjusted p-value")

abline(h = -log10(0.05), col = "red", lty = 2)  # soglia FDR
abline(v = c(-1, 1), col = "blue", lty = 2)     # soglia logFC

dev.off()

# ======== Volcano Plot 2: geni significativi evidenziati =========

# Definizione colori: nero = non significativo, arancio = up, verde = down
color <- ifelse(results$adj_pval > 0.05 | abs(results$logFC) <= 1, "black",
                ifelse(results$logFC > 0, "orange", "seagreen"))

png(paste0(dirDataset, "volcano_plot_2.png"), width = 800, height = 600)

plot(results$logFC, -log10(results$adj_pval),
     pch = 20, col = color,
     main = "Volcano Plot - Geni significativi",
     xlab = "log2 Fold Change",
     ylab = "-log10 Adjusted p-value")

abline(h = -log10(0.05), col = "red", lty = 2)  # soglia FDR
abline(v = c(-1, 1), col = "blue", lty = 2)     # soglia logFC

dev.off()

# ===========================================
# STEP 6: Boxplot of most significant DEG
# ===========================================
#Box-plot with expression levels of the most down-regulated
gene <- DEGs$gene[which.min(DEGs$adj_pval)]
x <- as.numeric(dataC[gene, ])
y <- as.numeric(dataN[gene, ])
pval_gene <- t.test(x, y, paired = TRUE)$p.value

png("boxplot_min_pval_gene.png", width = 800, height = 600)
boxplot(c(x, y) ~ rep(c("Tumor", "Normal"), each = length(x)), col = c("tomato", "skyblue"), main = paste("Boxplot –", gene), ylab = "Expression (log2)", xlab = "Group")
legend("topright", legend = paste("p =", signif(pval_gene, 3)), bty = "n")
dev.off()

#Box-plot with expression levels of the most up-regulated
gene <- DEGs$gene[which.max(DEGs$adj_pval)]
x <- as.numeric(dataC[gene, ])
y <- as.numeric(dataN[gene, ])
pval_gene <- t.test(x, y, paired = TRUE)$p.value

png("boxplot_max_pval_gene.png", width = 800, height = 600)
boxplot(c(x, y) ~ rep(c("Tumor", "Normal"), each = length(x)), col = c("tomato", "skyblue"), main = paste("Boxplot –", gene), ylab = "Expression (log2)", xlab = "Group")
legend("topright", legend = paste("p =", signif(pval_gene, 3)), bty = "n")
dev.off()

# ===========================================
# STEP 7: Heatmap of top 50 DEGs
# ===========================================
top_genes <- rownames(DEGs[order(DEGs$adj_pval), ])[1:50]
ann <- data.frame(Group = rep(c("Normal", "Tumor"), c(N, M)))
rownames(ann) <- colnames(filtered_data)

png("heatmap_top50_DEGs.png", width = 1000, height = 1200)
pheatmap(filtered_data[top_genes, ], scale = "row", show_rownames = TRUE, show_colnames = FALSE, color = colorRampPalette(c("navy", "white", "skyblue"))(50), main = "Top 50 DEGs - Heatmap", annotation_col = ann)
dev.off()

# ===========================================
# STEP 6: PCA on DEGs
# ===========================================
data_pca <- t(filtered_data[rownames(DEGs), ])
group_pca <- rep(c("Normal", "Tumor"), c(N, M))
pca <- prcomp(data_pca, scale. = TRUE)

png("pca_plot.png", width = 800, height = 600)
plot(pca$x[, 1:2], col = ifelse(group_pca == "Tumor", "tomato", "skyblue"), pch = 19, xlab = "PC1", ylab = "PC2", main = "PCA - DEGs")
legend("topright", legend = c("Tumor", "Normal"), col = c("tomato", "skyblue"), pch = 19)
dev.off()

# ===========================================
# STEP 7: Export final DEGs
# ===========================================
write.table(DEGs, file = paste0(dirDataset, "DEGs_filtered.txt"), sep = "\t", quote = FALSE, row.names = FALSE)

# ===========================================
# STEP 8: Functional Enrichment (GO)
# ===========================================
genes <- DEGs$gene
entrez <- bitr(genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
ego <- enrichGO(gene = entrez$ENTREZID, OrgDb = org.Hs.eg.db, ont = "BP", pAdjustMethod = "fdr", qvalueCutoff = 0.05, readable = TRUE)

#enrichment_barplot
png("GO_enrichment_barplot.png", width = 1000, height = 600)
barplot(ego, showCategory = 10, title = "GO Enrichment - Biological Process")
dev.off()

#KEGG pathway enrichment
png("KEGG_barplot.png", width = 1000, height = 600)
ekegg <- enrichKEGG(gene = entrez$ENTREZID, organism = 'hsa', pAdjustMethod = 'fdr', qvalueCutoff = 0.05)
barplot(ekegg, showCategory = 10, title = "KEGG Pathway Enrichment")

# ===========================================
# Save gene names for external tools
# ===========================================
write.table(DEGs$gene, file = paste0(dirDataset, "top_genes.txt"), quote = FALSE, row.names = FALSE, col.names = FALSE)
