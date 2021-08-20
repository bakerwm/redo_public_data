### Analysis of auxin treated embryos

library("ggplot2")
library("DESeq2")
library("pheatmap")
library(ggrepel)
library(stringr)
# library(svglite)
# library("ggpubr")

setwd("~/work/public/2021_elife_maternal_piwi_greg/results/redo_RNAseq/DESeq2_output")

##----------------------------------------------------------------------------##
## BEGIN: by WM
## TE - degron: by htseq-count, sense
directory <- "count"
f_list <- list.files(directory, ".gene_te.count.htseq$", recursive = TRUE, full.names = TRUE)
f_list <- purrr::discard(f_list, grepl("gap90|old.files|w1118", f_list))

## construct table
sampleTable <- data.frame(
  sampleName = f_list,
  fileName   = basename(f_list),
  condition  = stringr::str_extract(basename(f_list), "(auxin|pbs)"),
  batch      = rep(c("a", "b", "c"), 2)
) %>%
  dplyr::mutate(condition = factor(condition , c("pbs", "auxin")),
                batch     = as.factor(batch))

ddsHTSeq <- DESeqDataSetFromHTSeqCount(sampleTable = sampleTable,
                                       directory   = directory,
                                       design      = ~ batch + condition)
## load te names
df_te <- readr::read_delim(
  "count/RNA-Seq_degron_embryo_2_2.5h_auxin_treat_rna_rep1.te.count.htseq",
  comment = "_", col_names = c("gene", "count"), col_types = readr::cols())

dds <- DESeq(ddsHTSeq)
res <- results(dds)
saveRDS(dds, "degron.DESeq2.dds.rds")

resOrdered <- res[order(res$pvalue),]
plotMA(res, ylim=c(-4,4))

# vsd <- vst(dds, blind=FALSE)
# rld <- rlog(dds, blind=FALSE)
# head(assay(vsd), 3)
#
# plotPCA(vsd, intgroup=c("condition"))

results <- as.data.frame(res) %>%
  tibble::rownames_to_column("gene") %>%
  dplyr::mutate(
    te       = gene %in% df_te$gene,
    baseMean = log10(baseMean),
    sig      = ifelse(is.na(padj) | padj >= 0.05, "not",
                      ifelse(padj < 0.05 & log2FoldChange <= -1, "down",
                             ifelse(padj < 0.05 & log2FoldChange >= 1, "up", "not"))),
    fc       = ifelse(log2FoldChange < 1, "fc<2",
                      ifelse(log2FoldChange > 1 & padj < 0.05, "fc<2 & p<0.05", "not"))
  )

## subset
results_gene <- dplyr::filter(results, ! te)
results_te   <- dplyr::filter(results, te)
results_roo  <- dplyr::filter(results, gene %in% c("297", "roo"))

p1 <- ggplot(results, aes(baseMean, log2FoldChange)) +
  geom_point(data = results_gene,
             shape = 1, size = 1, color = "grey70") +
  geom_point(data = results_te,
             shape = 16, size = 2, color = "darkgreen") +
  geom_point(data = results_roo,
             shape = 16, size = 2, color = "red") +
  ggrepel::geom_text_repel(aes(baseMean, log2FoldChange, label = gene),
                           data = results_roo) +
  geom_hline(yintercept = c(-1, 1), linetype = 2, color = "grey50") +
  scale_x_continuous(limits = c(1, 5), name = "log10 mean expression") +
  scale_y_continuous(limits = c(-1.6, 1.6), name = "log2 Fold Change auxin vs ctrl") +
  ggtitle("degron") +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    axis.text = element_text(size = 15),
    axis.title = element_text(size = 15)
  )
##----------------------------------------------------------------------------##



##----------------------------------------------------------------------------##
## for not lfcShrink
resLFC <- lfcShrink(dds, coef="condition_auxin_vs_pbs", type="apeglm")
plotMA(resLFC, ylim=c(-2,2))

resultsLFC <- as.data.frame(resLFC) %>%
  tibble::rownames_to_column("gene") %>%
  dplyr::mutate(
    te       = gene %in% df_te$gene,
    baseMean = log10(baseMean),
    sig      = ifelse(is.na(padj) | padj >= 0.05, "not",
                      ifelse(padj < 0.05 & log2FoldChange <= -1, "down",
                             ifelse(padj < 0.05 & log2FoldChange >= 1, "up", "not"))),
    fc       = ifelse(log2FoldChange < 1, "fc<2",
                      ifelse(log2FoldChange > 1 & padj < 0.05, "fc<2 & p<0.05", "not"))
  )

resultsLFC_gene <- dplyr::filter(resultsLFC, ! te)
resultsLFC_te   <- dplyr::filter(resultsLFC, te)
resultsLFC_roo  <- dplyr::filter(resultsLFC, gene %in% c("297", "roo"))

p2 <- ggplot(resultsLFC, aes(baseMean, log2FoldChange)) +
  geom_point(data = resultsLFC_gene,
             shape = 1, size = 1, color = "grey70") +
  geom_point(data = resultsLFC_te,
             shape = 16, size = 2, color = "darkgreen") +
  geom_point(data = resultsLFC_roo,
             shape = 16, size = 2, color = "red") +
  ggrepel::geom_text_repel(aes(baseMean, log2FoldChange, label = gene),
                           data = resultsLFC_roo) +
  geom_hline(yintercept = c(-1, 1), linetype = 2, color = "grey50") +
  scale_x_continuous(limits = c(1, 5), name = "log10 mean expression") +
  scale_y_continuous(limits = c(-1.6, 1.6), name = "log2 Fold Change auxin vs ctrl") +
  ggtitle("degron") +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    axis.text = element_text(size = 15),
    axis.title = element_text(size = 15)
  )


## add annotation
p1a <- p1 + ggtitle("not estimate: batch + condition")
p1b <- p1a +
  scale_x_continuous(limits = c(0, 5), name = "log10 mean expression") +
  ggtitle("not estimate: batch + condition")

p2a <- p2 + ggtitle("lfcShrink: batch + condition")
p2b <- p2a +
  scale_x_continuous(limits = c(0, 5), name = "log10 mean expression") +
  ggtitle("lfcShrink: batch + condition")


##----------------------------------------------------------------------------##
## merge
df3a <- results %>%
  dplyr::filter(gene %in% c("roo", "297"))
df3b <- resultsLFC %>%
  dplyr::filter(! gene %in% c("roo", "297"))
df3 <- dplyr::bind_rows(df3a, df3b)
df3_gene <- dplyr::filter(df3, ! te)
df3_te   <- dplyr::filter(df3, te)
df3_roo  <- dplyr::filter(df3, gene %in% c("297", "roo"))


p3a <- ggplot(df3, aes(baseMean, log2FoldChange)) +
  geom_point(data = df3_gene,
             shape = 1, size = 1, color = "grey70") +
  geom_point(data = df3_te,
             shape = 16, size = 2, color = "darkgreen") +
  geom_point(data = df3_roo,
             shape = 16, size = 2, color = "red") +
  ggrepel::geom_text_repel(aes(baseMean, log2FoldChange, label = gene),
                           data = df3_roo) +
  geom_hline(yintercept = c(-1, 1), linetype = 2, color = "grey50") +
  scale_x_continuous(limits = c(1, 5), name = "log10 mean expression") +
  scale_y_continuous(limits = c(-1.6, 1.6), name = "log2 Fold Change auxin vs ctrl") +
  ggtitle("merge: batch + condition") +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    axis.text = element_text(size = 15),
    axis.title = element_text(size = 15)
  )

p3b <- p3a +
  scale_x_continuous(limits = c(0, 5), name = "log10 mean expression")

## combine png files
p1x <- wrap_plots(p1a, p1b, p2a, p2b, p3a, p3b, ncol = 2) +
  patchwork::plot_annotation(title = "Degron")

ggsave("DESeq2.degron.lfc.ma.png", p1x, width = 8, height = 11, dpi = 300)







##----------------------------------------------------------------------------##
## BEGIN: gap90
##
## TE, by htseq-count, sense
directory <- "count"
f_list <- list.files(directory, ".gene_te.count.htseq$", recursive = TRUE, full.names = TRUE)
f_list <- purrr::discard(f_list, grepl("degron|old.files|w1118", f_list))

## construct table
sampleTable <- data.frame(
  sampleName = f_list,
  fileName   = basename(f_list),
  condition  = stringr::str_extract(basename(f_list), "(auxin|pbs)"),
  batch      = rep(c("a", "b", "c"), 2)
) %>%
  dplyr::mutate(condition = factor(condition , c("pbs", "auxin")))

ddsHTSeq <- DESeqDataSetFromHTSeqCount(sampleTable = sampleTable,
                                       directory   = directory,
                                       design      = ~ batch + condition)
## load te names
df_te <- readr::read_delim(
  "count/RNA-Seq_gap90_embryo_2_2.5h_auxin_treat_rna_rep1.te.count.htseq",
  comment = "_", col_names = c("gene", "count"), col_types = readr::cols())

dds <- DESeq(ddsHTSeq)
res <- results(dds)
saveRDS(dds, "gap90.DESeq2.dds.rds")

resOrdered <- res[order(res$pvalue),]
plotMA(res, ylim=c(-4,4))

# vsd <- vst(dds, blind=FALSE)
# rld <- rlog(dds, blind=FALSE)
# head(assay(vsd), 3)
# plotPCA(vsd, intgroup=c("condition"))

results <- as.data.frame(res) %>%
  tibble::rownames_to_column("gene") %>%
  dplyr::mutate(
    te       = gene %in% df_te$gene,
    baseMean = log10(baseMean),
    sig      = ifelse(is.na(padj) | padj >= 0.05, "not",
                      ifelse(padj < 0.05 & log2FoldChange <= -1, "down",
                             ifelse(padj < 0.05 & log2FoldChange >= 1, "up", "not"))),
    fc       = ifelse(log2FoldChange < 1, "fc<2",
                      ifelse(log2FoldChange > 1 & padj < 0.05, "fc<2 & p<0.05", "not"))
  )

## subset
results_gene <- dplyr::filter(results, ! te)
results_te   <- dplyr::filter(results, te)
results_roo  <- dplyr::filter(results, gene %in% c("297", "roo"))

p1 <- ggplot(results, aes(baseMean, log2FoldChange)) +
  geom_point(data = results_gene,
             shape = 1, size = 1, color = "grey70") +
  geom_point(data = results_te,
             shape = 16, size = 2, color = "darkgreen") +
  geom_point(data = results_roo,
             shape = 16, size = 2, color = "red") +
  ggrepel::geom_text_repel(aes(baseMean, log2FoldChange, label = gene),
                           data = results_roo) +
  geom_hline(yintercept = c(-1, 1), linetype = 2, color = "grey50") +
  scale_x_continuous(limits = c(1, 5), name = "log10 mean expression") +
  scale_y_continuous(limits = c(-1.6, 1.6), name = "log2 Fold Change auxin vs ctrl") +
  ggtitle("degron") +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    axis.text = element_text(size = 15),
    axis.title = element_text(size = 15)
  )


##----------------------------------------------------------------------------##
## for not lfcShrink
resLFC <- lfcShrink(dds, coef="condition_auxin_vs_pbs", type="apeglm")
plotMA(resLFC, ylim=c(-2,2))

resultsLFC <- as.data.frame(resLFC) %>%
  tibble::rownames_to_column("gene") %>%
  dplyr::mutate(
    te       = gene %in% df_te$gene,
    baseMean = log10(baseMean),
    sig      = ifelse(is.na(padj) | padj >= 0.05, "not",
                      ifelse(padj < 0.05 & log2FoldChange <= -1, "down",
                             ifelse(padj < 0.05 & log2FoldChange >= 1, "up", "not"))),
    fc       = ifelse(log2FoldChange < 1, "fc<2",
                      ifelse(log2FoldChange > 1 & padj < 0.05, "fc<2 & p<0.05", "not"))
  )

resultsLFC_gene <- dplyr::filter(resultsLFC, ! te)
resultsLFC_te   <- dplyr::filter(resultsLFC, te)
resultsLFC_roo  <- dplyr::filter(resultsLFC, gene %in% c("297", "roo"))

p2 <- ggplot(resultsLFC, aes(baseMean, log2FoldChange)) +
  geom_point(data = resultsLFC_gene,
             shape = 1, size = 1, color = "grey70") +
  geom_point(data = resultsLFC_te,
             shape = 16, size = 2, color = "darkgreen") +
  geom_point(data = resultsLFC_roo,
             shape = 16, size = 2, color = "red") +
  ggrepel::geom_text_repel(aes(baseMean, log2FoldChange, label = gene),
                           data = resultsLFC_roo) +
  geom_hline(yintercept = c(-1, 1), linetype = 2, color = "grey50") +
  scale_x_continuous(limits = c(1, 5), name = "log10 mean expression") +
  scale_y_continuous(limits = c(-1.6, 1.6), name = "log2 Fold Change auxin vs ctrl") +
  ggtitle("degron") +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    axis.text = element_text(size = 15),
    axis.title = element_text(size = 15)
  )


## add annotation
p1a <- p1 + ggtitle("not estimate: batch + condition")
p1b <- p1a +
  scale_x_continuous(limits = c(0, 5), name = "log10 mean expression") +
  ggtitle("not estimate: batch + condition")

p2a <- p2 + ggtitle("lfcShrink: batch + condition")
p2b <- p2a +
  scale_x_continuous(limits = c(0, 5), name = "log10 mean expression") +
  ggtitle("lfcShrink: batch + condition")


##----------------------------------------------------------------------------##
## merge
df3a <- results %>%
  dplyr::filter(gene %in% c("roo", "297"))
df3b <- resultsLFC %>%
  dplyr::filter(! gene %in% c("roo", "297"))
df3 <- dplyr::bind_rows(df3a, df3b)
df3_gene <- dplyr::filter(df3, ! te)
df3_te   <- dplyr::filter(df3, te)
df3_roo  <- dplyr::filter(df3, gene %in% c("297", "roo"))


p3a <- ggplot(df3, aes(baseMean, log2FoldChange)) +
  geom_point(data = df3_gene,
             shape = 1, size = 1, color = "grey70") +
  geom_point(data = df3_te,
             shape = 16, size = 2, color = "darkgreen") +
  geom_point(data = df3_roo,
             shape = 16, size = 2, color = "red") +
  ggrepel::geom_text_repel(aes(baseMean, log2FoldChange, label = gene),
                           data = df3_roo) +
  geom_hline(yintercept = c(-1, 1), linetype = 2, color = "grey50") +
  scale_x_continuous(limits = c(1, 5), name = "log10 mean expression") +
  scale_y_continuous(limits = c(-1.6, 1.6), name = "log2 Fold Change auxin vs ctrl") +
  ggtitle("merge: batch + condition") +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    axis.text = element_text(size = 15),
    axis.title = element_text(size = 15)
  )

p3b <- p3a +
  scale_x_continuous(limits = c(0, 5), name = "log10 mean expression")

## combine png files
p1x <- wrap_plots(p1a, p1b, p2a, p2b, p3a, p3b, ncol = 2) +
  patchwork::plot_annotation(title = "Degron")

ggsave("DESeq2.gap90.lfc.ma.png", p1x, width = 8, height = 11, dpi = 300)



