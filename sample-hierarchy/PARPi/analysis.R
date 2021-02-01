# load packages----
setRepositories(ind = c(1:8))
for (package in c('biomaRt', 'tidyverse', 'tximport', 'ensembldb', 'EnsDb.Hsapiens.v86', 
	'edgeR', 'matrixStats', 'cowplot',
	'DT', 'plotly', 'gt'
	'limma', 'ggrepel',
	'GSEABase', 'Biobase', 'GSVA', 'gprofiler2', 'clusterProfiler', 'msigdbr', 'enrichplot',
	'ggthemes')) {
    if (!require(package, character.only=T, quietly=T)) {
        install.packages(package)
        library(package, character.only=T)
    }
}

mart <- useMart("ENSEMBL_MART_ENSEMBL", host = "useast.ensembl.org") # coding genes for cleaning
mart <- useDataset("hsapiens_gene_ensembl", mart)
coding_genes <- getBM(attributes = c( "hgnc_symbol"), 
                      filters = c("biotype"), 
                      values = list(biotype="protein_coding"), 
                      mart = mart)$hgnc_symbol

# import transcripts----
targets <- read_tsv("studydesign.txt")
paths <- file.path(targets$dir, "abundance.tsv") # set file paths to your mapped data
Tx <- transcripts(EnsDb.Hsapiens.v86, columns=c("tx_id", "gene_name"))
Tx <- as_tibble(Tx)
Tx <- dplyr::rename(Tx, target_id = tx_id)
Tx <- dplyr::select(Tx, "target_id", "gene_name")
Txi_gene <- tximport(paths, 
                     type = "kallisto", 
                     tx2gene = Tx, 
                     txOut = FALSE,
                     countsFromAbundance = "lengthScaledTPM",
                     ignoreTxVersion = TRUE)

## define important global variables
sampleLabels <- targets$sample
group <- factor(targets$group)

# log2(CountsPerMillion)----
# unfiltered, non-normalized
myDGEList <- DGEList(Txi_gene$counts)
log2.cpm <- cpm(myDGEList, log=TRUE)
log2.cpm.df <- as_tibble(log2.cpm, rownames = "geneID") #converts to .df
colnames(log2.cpm.df) <- c("geneID", sampleLabels) #renames columns

# filtered, non-normalized
cpm <- cpm(myDGEList)
keepers <- rowSums(cpm>1)>=(dim(targets)[1]/length(levels(group))) #filters wherever the CPM < 1 in BOTH samples
myDGEList.filtered <- myDGEList[keepers,]
log2.cpm.filtered <- cpm(myDGEList.filtered, log=TRUE)
log2.cpm.filtered.df <- as_tibble(log2.cpm.filtered, rownames = "geneID")
colnames(log2.cpm.filtered.df) <- c("geneID", sampleLabels)

# filtered, normalized
myDGEList.filtered.norm <- calcNormFactors(myDGEList.filtered, method = "TMM") #normalize using TMM
log2.cpm.filtered.norm <- cpm(myDGEList.filtered.norm, log=TRUE)
log2.cpm.filtered.norm.df <- as_tibble(log2.cpm.filtered.norm, rownames = "geneID")
colnames(log2.cpm.filtered.norm.df) <- c("geneID", sampleLabels)

# multivariate analysis----
## distance mapping
distance <- dist(t(log2.cpm.filtered.norm), method = "maximum")
clusters <- hclust(distance, method = "average")
dendrogram <- plot(clusters, labels=sampleLabels, cex=1.2)

## pca
pca.res <- prcomp(t(log2.cpm.filtered.norm), scale.=F, retx=T)
pc.var <- pca.res$sdev^2
pc.per <- round(pc.var/sum(pc.var)*100, 1) 
pca.res.df <- as_tibble(pca.res$x)
pca.plot <- ggplot(pca.res.df) +
  aes(x=PC1, y=PC2, label=sampleLabels, color = group) +
  geom_point(size=3) +
  #stat_ellipse() +
  xlab(paste0("PC1 (",pc.per[1],"%",")")) + 
  ylab(paste0("PC2 (",pc.per[2],"%",")")) +
  labs(title="PCA plot (filtered, normalized)") +
  #     caption=paste0("produced on ", Sys.time())) +
  #coord_fixed() +
  #theme_bw() + 
  theme_fivethirtyeight()
ggsave(path = "./plots/dim", filename = "pca.png", plot = pca.plot)

# model matrix generation----
mydata.df <- mutate(log2.cpm.filtered.norm.df) %>% mutate_if(is.numeric, round, 2)
design <- model.matrix(~0 + group)
colnames(design) <- levels(group)
v.DEGList.filtered.norm <- voom(myDGEList.filtered.norm, design, plot = FALSE)
fit <- lmFit(v.DEGList.filtered.norm, design)
mydata.df$Control.AVG <- rowMeans(subset(mydata.df, select = c("Control_sample1", "Control_sample2")), na.rm = TRUE)
mydata.df$KO5.AVG <- rowMeans(subset(mydata.df, select = c("KO5_sample1", "KO5_sample2")), na.rm = TRUE)
mydata.df$E6.AVG <- rowMeans(subset(mydata.df, select = c("E6_sample1", "E6_sample2")), na.rm = TRUE)
mydata.df$LacZ.AVG <- rowMeans(subset(mydata.df, select = c("LacZ_sample1", "LacZ_sample2")), na.rm = TRUE)
mydata.df$WT.AVG <- rowMeans(subset(mydata.df, select = c("WT_sample1", "WT_sample2")), na.rm = TRUE)

## KO5 vs E6----
mydata.df$KO5.E6 <- mydata.df$E6.AVG - mydata.df$KO5.AVG
KO5.E6.contrast.matrix <- makeContrasts(infection = E6 - KO5,
                                 levels=design)
KO5.E6.fits <- contrasts.fit(fit, KO5.E6.contrast.matrix)
KO5.E6.ebFit <- eBayes(KO5.E6.fits)
KO5.E6.hits <- topTable(KO5.E6.ebFit, adjust ="BH", coef=1, number=40000, sort.by="logFC")
KO5.E6.hits.df <- KO5.E6.hits %>%
  as_tibble(rownames = "geneID")
KO5.E6.hits.df$significant <- ((KO5.E6.hits.df$adj.P.Val < 0.05) & (abs(KO5.E6.hits.df$logFC) > 1))
KO5.E6.hits.df$sorter <- -1*abs(-log10(KO5.E6.hits.df$adj.P.Val)*KO5.E6.hits.df$logFC)
KO5.E6.hits.df <- KO5.E6.hits.df[order(KO5.E6.hits.df$sorter),]
KO5.E6.hits.df <- dplyr::filter(KO5.E6.hits.df, geneID %in% coding_genes)
KO5_E6_vplot <- ggplot(KO5.E6.hits.df) +
  aes(y=-log10(adj.P.Val), x=logFC, colour=significant) +
  scale_color_manual(values = c("TRUE" = "black", "FALSE"="grey")) +
  geom_point(size=1, alpha=.2) +
  geom_text_repel(size=3, data=head(KO5.E6.hits.df, 20), aes(label=geneID)) +
  geom_hline(yintercept = -log10(0.05), linetype="longdash", colour="black", size=.5) +
  geom_vline(xintercept = 1, linetype="longdash", colour="#BE684D", size=.5) +
  geom_vline(xintercept = -1, linetype="longdash", colour="#2C467A", size=.5) +
  #annotate("rect", xmin = 1, xmax = 10, ymin = -log2(0.05), ymax = 6, alpha=.2, fill="#BE684D") +
  #annotate("rect", xmin = -1, xmax = -10, ymin = -log2(0.05), ymax = 6, alpha=.2, fill="#2C467A") +
  labs(title="Volcano plot",
       subtitle = "KO5 to E6") +
  theme_fivethirtyeight()
ggsave(path = "./plots/vplots", filename = "KO5_E6_vplot.png", plot = KO5_E6_vplot)

write_csv(KO5.E6.hits.df, "./dge/KO5_E6.csv")

## control vs KO5----
mydata.df$Control.KO5 <- mydata.df$KO5.AVG - mydata.df$Control.AVG
Control.KO5.contrast.matrix <- makeContrasts(infection = KO5 - Control,
                                 levels=design)
Control.KO5.fits <- contrasts.fit(fit, Control.KO5.contrast.matrix)
Control.KO5.ebFit <- eBayes(Control.KO5.fits)
Control.KO5.hits <- topTable(Control.KO5.ebFit, adjust ="BH", coef=1, number=40000, sort.by="logFC")
Control.KO5.hits.df <- Control.KO5.hits %>%
  as_tibble(rownames = "geneID")
Control.KO5.hits.df$significant <- ((Control.KO5.hits.df$adj.P.Val < 0.05) & (abs(Control.KO5.hits.df$logFC) > 1))
Control.KO5.hits.df$sorter <- -1*abs(-log10(Control.KO5.hits.df$adj.P.Val)*Control.KO5.hits.df$logFC)
Control.KO5.hits.df <- Control.KO5.hits.df[order(Control.KO5.hits.df$sorter),]
Control.KO5.hits.df <- dplyr::filter(Control.KO5.hits.df, geneID %in% coding_genes)
Control_KO5_vplot <- ggplot(Control.KO5.hits.df) +
  aes(y=-log10(adj.P.Val), x=logFC, colour=significant) +
  scale_color_manual(values = c("TRUE" = "black", "FALSE"="grey")) +
  geom_point(size=1, alpha=.2) +
  geom_text_repel(size=3, data=head(Control.KO5.hits.df, 20), aes(label=geneID)) +
  geom_hline(yintercept = -log10(0.05), linetype="longdash", colour="grey", size=.5) +
  geom_vline(xintercept = 1, linetype="longdash", colour="#BE684D", size=.5) +
  geom_vline(xintercept = -1, linetype="longdash", colour="#2C467A", size=.5) +
  #annotate("rect", xmin = 1, xmax = 10, ymin = -log2(0.05), ymax = 6, alpha=.2, fill="#BE684D") +
  #annotate("rect", xmin = -1, xmax = -10, ymin = -log2(0.05), ymax = 6, alpha=.2, fill="#2C467A") +
  labs(title="Volcano plot",
       subtitle = "Control to KO5") +
#       caption=paste0("produced on ", Sys.time())) +
  theme_fivethirtyeight()
ggsave(path = "./plots/vplots", filename = "Control_KO5_vplot.png", plot = Control_KO5_vplot)

write_csv(Control.KO5.hits.df, "./dge/Control_KO5.csv")

## control vs E6----
mydata.df$Control.E6 <- mydata.df$E6.AVG - mydata.df$Control.AVG
Control.E6.contrast.matrix <- makeContrasts(infection = E6 - Control,
                                             levels=design)
Control.E6.fits <- contrasts.fit(fit, Control.E6.contrast.matrix)
Control.E6.ebFit <- eBayes(Control.E6.fits)
Control.E6.hits <- topTable(Control.E6.ebFit, adjust ="BH", coef=1, number=40000, sort.by="logFC")
Control.E6.hits.df <- Control.E6.hits %>%
  as_tibble(rownames = "geneID")
Control.E6.hits.df$significant <- ((Control.E6.hits.df$adj.P.Val < 0.05) & (abs(Control.E6.hits.df$logFC) > 1))
Control.E6.hits.df$sorter <- -1*abs(-log10(Control.E6.hits.df$adj.P.Val)*Control.E6.hits.df$logFC)
Control.E6.hits.df <- Control.E6.hits.df[order(Control.E6.hits.df$sorter),]
Control.E6.hits.df <- dplyr::filter(Control.E6.hits.df, geneID %in% coding_genes)
Control_E6_vplot <- ggplot(Control.E6.hits.df) +
  aes(y=-log10(adj.P.Val), x=logFC, colour=significant) +
  scale_color_manual(values = c("TRUE" = "black", "FALSE"="grey")) +
  geom_point(size=1, alpha=.2) +
  geom_text_repel(size=3, data=head(Control.E6.hits.df, 20), aes(label=geneID)) +
  geom_hline(yintercept = -log10(0.05), linetype="longdash", colour="grey", size=.5) +
  geom_vline(xintercept = 1, linetype="longdash", colour="#BE684D", size=.5) +
  geom_vline(xintercept = -1, linetype="longdash", colour="#2C467A", size=.5) +
  #annotate("rect", xmin = 1, xmax = 10, ymin = -log2(0.05), ymax = 6, alpha=.2, fill="#BE684D") +
  #annotate("rect", xmin = -1, xmax = -10, ymin = -log2(0.05), ymax = 6, alpha=.2, fill="#2C467A") +
  labs(title="Volcano plot",
       subtitle = "Control to E6") +
  #       caption=paste0("produced on ", Sys.time())) +
  theme_fivethirtyeight()
ggsave(path = "./plots/vplots", filename = "Control_E6_vplot.png", plot = Control_E6_vplot)

write_csv(Control.E6.hits.df, "./dge/Control_E6.csv")

## control vs LacZ----
mydata.df$Control.LacZ <- mydata.df$LacZ.AVG - mydata.df$Control.AVG
Control.LacZ.contrast.matrix <- makeContrasts(infection = LacZ - Control,
                                              levels=design)
Control.LacZ.fits <- contrasts.fit(fit, Control.LacZ.contrast.matrix)
Control.LacZ.ebFit <- eBayes(Control.LacZ.fits)
Control.LacZ.hits <- topTable(Control.LacZ.ebFit, adjust ="BH", coef=1, number=40000, sort.by="logFC")
Control.LacZ.hits.df <- Control.LacZ.hits %>%
  as_tibble(rownames = "geneID")
Control.LacZ.hits.df$significant <- ((Control.LacZ.hits.df$adj.P.Val < 0.05) & (abs(Control.LacZ.hits.df$logFC) > 1))
Control.LacZ.hits.df$sorter <- -1*abs(-log10(Control.LacZ.hits.df$adj.P.Val)*Control.LacZ.hits.df$logFC)
Control.LacZ.hits.df <- Control.LacZ.hits.df[order(Control.LacZ.hits.df$sorter),]
Control.LacZ.hits.df <- dplyr::filter(Control.LacZ.hits.df, geneID %in% coding_genes)
Control_LacZ_vplot <- ggplot(Control.LacZ.hits.df) +
  aes(y=-log10(adj.P.Val), x=logFC, colour=significant) +
  scale_color_manual(values = c("TRUE" = "black", "FALSE"="grey")) +
  geom_point(size=1, alpha=.2) +
  geom_text_repel(size=3, data=head(Control.LacZ.hits.df, 20), aes(label=geneID)) +
  geom_hline(yintercept = -log10(0.05), linetype="longdash", colour="grey", size=.5) +
  geom_vline(xintercept = 1, linetype="longdash", colour="#BE684D", size=.5) +
  geom_vline(xintercept = -1, linetype="longdash", colour="#2C467A", size=.5) +
  #annotate("rect", xmin = 1, xmax = 10, ymin = -log2(0.05), ymax = 6, alpha=.2, fill="#BLacZ84D") +
  #annotate("rect", xmin = -1, xmax = -10, ymin = -log2(0.05), ymax = 6, alpha=.2, fill="#2C467A") +
  labs(title="Volcano plot",
       subtitle = "Control to LacZ") +
  #       caption=paste0("produced on ", Sys.time())) +
  theme_fivethirtyeight()
ggsave(path = "./plots/vplots", filename = "Control_LacZ_vplot.png", plot = Control_LacZ_vplot)

write_csv(Control.LacZ.hits.df, "./dge/Control_LacZ.csv")

## KO5 vs LacZ----
mydata.df$KO5.LacZ <- mydata.df$LacZ.AVG - mydata.df$KO5.AVG
KO5.LacZ.contrast.matrix <- makeContrasts(infection = LacZ - KO5,
                                        levels=design)
KO5.LacZ.fits <- contrasts.fit(fit, KO5.LacZ.contrast.matrix)
KO5.LacZ.ebFit <- eBayes(KO5.LacZ.fits)
KO5.LacZ.hits <- topTable(KO5.LacZ.ebFit, adjust ="BH", coef=1, number=40000, sort.by="logFC")
KO5.LacZ.hits.df <- KO5.LacZ.hits %>%
  as_tibble(rownames = "geneID")
KO5.LacZ.hits.df$significant <- ((KO5.LacZ.hits.df$adj.P.Val < 0.05) & (abs(KO5.LacZ.hits.df$logFC) > 1))
KO5.LacZ.hits.df$sorter <- -1*abs(-log10(KO5.LacZ.hits.df$adj.P.Val)*KO5.LacZ.hits.df$logFC)
KO5.LacZ.hits.df <- KO5.LacZ.hits.df[order(KO5.LacZ.hits.df$sorter),]
KO5.LacZ.hits.df <- dplyr::filter(KO5.LacZ.hits.df, geneID %in% coding_genes)
KO5_LacZ_vplot <- ggplot(KO5.LacZ.hits.df) +
  aes(y=-log10(adj.P.Val), x=logFC, colour=significant) +
  scale_color_manual(values = c("TRUE" = "black", "FALSE"="grey")) +
  geom_point(size=1, alpha=.2) +
  geom_text_repel(size=3, data=head(KO5.LacZ.hits.df, 20), aes(label=geneID)) +
  geom_hline(yintercept = -log10(0.05), linetype="longdash", colour="black", size=.5) +
  geom_vline(xintercept = 1, linetype="longdash", colour="#BE684D", size=.5) +
  geom_vline(xintercept = -1, linetype="longdash", colour="#2C467A", size=.5) +
  #annotate("rect", xmin = 1, xmax = 10, ymin = -log2(0.05), ymax = 6, alpha=.2, fill="#BE684D") +
  #annotate("rect", xmin = -1, xmax = -10, ymin = -log2(0.05), ymax = 6, alpha=.2, fill="#2C467A") +
  labs(title="Volcano plot",
       subtitle = "KO5 to LacZ") +
  theme_fivethirtyeight()
ggsave(path = "./plots/vplots", filename = "KO5_LacZ_vplot.png", plot = KO5_LacZ_vplot)

write_csv(KO5.LacZ.hits.df, "./dge/KO5_LacZ.csv")

## KO5 vs WT----
mydata.df$KO5.WT <- mydata.df$WT.AVG - mydata.df$KO5.AVG
KO5.WT.contrast.matrix <- makeContrasts(infection = WT - KO5,
                                        levels=design)
KO5.WT.fits <- contrasts.fit(fit, KO5.WT.contrast.matrix)
KO5.WT.ebFit <- eBayes(KO5.WT.fits)
KO5.WT.hits <- topTable(KO5.WT.ebFit, adjust ="BH", coef=1, number=40000, sort.by="logFC")
KO5.WT.hits.df <- KO5.WT.hits %>%
  as_tibble(rownames = "geneID")
KO5.WT.hits.df$significant <- ((KO5.WT.hits.df$adj.P.Val < 0.05) & (abs(KO5.WT.hits.df$logFC) > 1))
KO5.WT.hits.df$sorter <- -1*abs(-log10(KO5.WT.hits.df$adj.P.Val)*KO5.WT.hits.df$logFC)
KO5.WT.hits.df <- KO5.WT.hits.df[order(KO5.WT.hits.df$sorter),]
KO5.WT.hits.df <- dplyr::filter(KO5.WT.hits.df, geneID %in% coding_genes)
KO5_WT_vplot <- ggplot(KO5.WT.hits.df) +
  aes(y=-log10(adj.P.Val), x=logFC, colour=significant) +
  scale_color_manual(values = c("TRUE" = "black", "FALSE"="grey")) +
  geom_point(size=1, alpha=.2) +
  geom_text_repel(size=3, data=head(KO5.WT.hits.df, 20), aes(label=geneID)) +
  geom_hline(yintercept = -log10(0.05), linetype="longdash", colour="black", size=.5) +
  geom_vline(xintercept = 1, linetype="longdash", colour="#BE684D", size=.5) +
  geom_vline(xintercept = -1, linetype="longdash", colour="#2C467A", size=.5) +
  #annotate("rect", xmin = 1, xmax = 10, ymin = -log2(0.05), ymax = 6, alpha=.2, fill="#BE684D") +
  #annotate("rect", xmin = -1, xmax = -10, ymin = -log2(0.05), ymax = 6, alpha=.2, fill="#2C467A") +
  labs(title="Volcano plot",
       subtitle = "KO5 to WT") +
  theme_fivethirtyeight()
ggsave(path = "./plots/vplots", filename = "KO5_WT_vplot.png", plot = KO5_WT_vplot)

write_csv(KO5.WT.hits.df, "./dge/KO5_WT.csv")

## LacZ vs WT----
mydata.df$LacZ.WT <- mydata.df$WT.AVG - mydata.df$LacZ.AVG
LacZ.WT.contrast.matrix <- makeContrasts(infection = WT - LacZ,
                                         levels=design)
LacZ.WT.fits <- contrasts.fit(fit, LacZ.WT.contrast.matrix)
LacZ.WT.ebFit <- eBayes(LacZ.WT.fits)
LacZ.WT.hits <- topTable(LacZ.WT.ebFit, adjust ="BH", coef=1, number=40000, sort.by="logFC")
LacZ.WT.hits.df <- LacZ.WT.hits %>%
  as_tibble(rownames = "geneID")
LacZ.WT.hits.df$significant <- ((LacZ.WT.hits.df$adj.P.Val < 0.05) & (abs(LacZ.WT.hits.df$logFC) > 1))
LacZ.WT.hits.df$sorter <- -1*abs(-log10(LacZ.WT.hits.df$adj.P.Val)*LacZ.WT.hits.df$logFC)
LacZ.WT.hits.df <- LacZ.WT.hits.df[order(LacZ.WT.hits.df$sorter),]
LacZ.WT.hits.df <- dplyr::filter(LacZ.WT.hits.df, geneID %in% coding_genes)
LacZ_WT_vplot <- ggplot(LacZ.WT.hits.df) +
  aes(y=-log10(adj.P.Val), x=logFC, colour=significant) +
  scale_color_manual(values = c("TRUE" = "black", "FALSE"="grey")) +
  geom_point(size=1, alpha=.2) +
  geom_text_repel(size=3, data=head(LacZ.WT.hits.df, 20), aes(label=geneID)) +
  geom_hline(yintercept = -log10(0.05), linetype="longdash", colour="black", size=.5) +
  geom_vline(xintercept = 1, linetype="longdash", colour="#BE684D", size=.5) +
  geom_vline(xintercept = -1, linetype="longdash", colour="#2C467A", size=.5) +
  #annotate("rect", xmin = 1, xmax = 10, ymin = -log2(0.05), ymax = 6, alpha=.2, fill="#BE684D") +
  #annotate("rect", xmin = -1, xmax = -10, ymin = -log2(0.05), ymax = 6, alpha=.2, fill="#2C467A") +
  labs(title="Volcano plot",
       subtitle = "LacZ to WT") +
  theme_fivethirtyeight()
ggsave(path = "./plots/vplots", filename = "LacZ_WT_vplot.png", plot = LacZ_WT_vplot)

write_csv(LacZ.WT.hits.df, "./dge/LacZ_WT.csv")

# overlap analysis----
## load and prep control --> KO5
Control.KO5.df <- read.csv("./dge/Control_KO5.csv")
colnames(Control.KO5.df) <- paste("Control_KO5", colnames(Control.KO5.df), sep = "_")
colnames(Control.KO5.df)[colnames(Control.KO5.df) == 'Control_KO5_geneID'] <- 'geneID'

## load and prep KO5 --> WT
KO5.WT.df <- read.csv("./dge/KO5_WT.csv")
colnames(KO5.WT.df) <- paste("KO5_WT", colnames(KO5.WT.df), sep = "_")
colnames(KO5.WT.df)[colnames(KO5.WT.df) == 'KO5_WT_geneID'] <- 'geneID'

## load and prep LacZ --> WT
LacZ.WT.df <- read.csv("./dge/LacZ_WT.csv")
colnames(LacZ.WT.df) <- paste("LacZ_WT", colnames(LacZ.WT.df), sep = "_")
colnames(LacZ.WT.df)[colnames(LacZ.WT.df) == 'LacZ_WT_geneID'] <- 'geneID'

## load and prep Control --> E6
Control.E6.df <- read.csv("./dge/Control_E6.csv")
colnames(Control.E6.df) <- paste("Control_E6", colnames(Control.E6.df), sep = "_")
colnames(Control.E6.df)[colnames(Control.E6.df) == 'Control_E6_geneID'] <- 'geneID'

## merge and sort Control.KO5.WT
Control.KO5.WT.df <- merge(Control.KO5.df, KO5.WT.df, by=c("geneID"))
Control.KO5.WT.df$sorter <- -1*abs(rowMeans(Control.KO5.WT.df[,c("Control_KO5_logFC", 
                                                                 "KO5_WT_logFC")]))
Control.KO5.WT.df <- Control.KO5.WT.df[order(Control.KO5.WT.df$sorter),]

write.csv(Control.KO5.WT.df, "./overlap/Control.KO5.WT.csv")

## merge and sort Control.LacZ.WT
Control.LacZ.WT.df <- merge(Control.KO5.df, LacZ.WT.df, by=c("geneID"))
Control.LacZ.WT.df$sorter <- -1*abs(rowMeans(Control.LacZ.WT.df[,c("Control_KO5_logFC", 
                                                                   "LacZ_WT_logFC")]))
Control.LacZ.WT.df <- Control.LacZ.WT.df[order(Control.LacZ.WT.df$sorter),]

write.csv(Control.LacZ.WT.df, "./overlap/Control.LacZ.WT.csv")

## merge and sort Control.KO5.E6
Control.KO5.E6.df <- merge(Control.KO5.df, Control.E6.df, by=c("geneID"))
Control.KO5.E6.df$sorter <- -1*abs(rowMeans(Control.KO5.E6.df[,c("Control_KO5_logFC", 
                                                                 "Control_E6_logFC")]))
Control.KO5.E6.df <- Control.KO5.E6.df[order(Control.KO5.E6.df$sorter),]

write.csv(Control.KO5.E6.df, "./overlap/Control.KO5.E6.csv")

## Control.KO5.WT analysis
Control.KO5.WT.significant <- Control.KO5.WT.df %>% 
  filter(Control.KO5.WT.df$Control_KO5_significant == "TRUE" & 
           Control.KO5.WT.df$KO5_WT_significant == "TRUE")
write.csv(Control.KO5.WT.significant, "./overlap/Control.KO5.WT_significant.csv")

Control.KO5.WT.reexpressed <- Control.KO5.WT.df %>% 
  filter(Control.KO5.WT.df$Control_KO5_logFC < 0 & 
           Control.KO5.WT.df$KO5_WT_logFC > 0)
write.csv(Control.KO5.WT.reexpressed, "./overlap/Control.KO5.WT_reexpressed.csv")

Control.KO5.WT.significant.reexpressed <- Control.KO5.WT.df %>% 
  filter(Control.KO5.WT.df$Control_KO5_significant == "TRUE" & 
           Control.KO5.WT.df$KO5_WT_significant == "TRUE" &
           Control.KO5.WT.df$Control_KO5_logFC < 0 & 
           Control.KO5.WT.df$KO5_WT_logFC > 0)
Control.KO5.WT.significant.reexpressed$sorter <- -1*abs(Control.KO5.WT.significant.reexpressed$KO5_WT_logFC-Control.KO5.WT.significant.reexpressed$Control_KO5_logFC)
Control.KO5.WT.significant.reexpressed <- Control.KO5.WT.significant.reexpressed[order(Control.KO5.WT.significant.reexpressed$sorter),]
write.csv(Control.KO5.WT.significant.reexpressed, "./overlap/Control.KO5.WT_significant_reexpressed.csv")

## Control.LacZ.WT analysis
Control.LacZ.WT.significant <- Control.LacZ.WT.df %>% 
  filter(Control.LacZ.WT.df$Control_KO5_significant == "TRUE" & 
           Control.LacZ.WT.df$LacZ_WT_significant == "TRUE")
write.csv(Control.LacZ.WT.significant, "./overlap/Control.LacZ.WT_significant.csv")

Control.LacZ.WT.reexpressed <- Control.LacZ.WT.df %>% 
  filter(Control.LacZ.WT.df$Control_KO5_logFC < 0 & 
           Control.LacZ.WT.df$LacZ_WT_logFC > 0)
write.csv(Control.LacZ.WT.reexpressed, "./overlap/Control.LacZ.WT_reexpressed.csv")

Control.LacZ.WT.significant.reexpressed <- Control.LacZ.WT.df %>% 
  filter(Control.LacZ.WT.df$Control_KO5_significant == "TRUE" & 
           Control.LacZ.WT.df$LacZ_WT_significant == "TRUE" &
           Control.LacZ.WT.df$Control_KO5_logFC < 0 & 
           Control.LacZ.WT.df$LacZ_WT_logFC > 0)
Control.LacZ.WT.significant.reexpressed$sorter <- -1*abs(Control.LacZ.WT.significant.reexpressed$LacZ_WT_logFC-Control.LacZ.WT.significant.reexpressed$Control_KO5_logFC)
Control.LacZ.WT.significant.reexpressed <- Control.LacZ.WT.significant.reexpressed[order(Control.LacZ.WT.significant.reexpressed$sorter),]
write.csv(Control.LacZ.WT.significant.reexpressed, "./overlap/Control.LacZ.WT_significant_reexpressed.csv")

## consistents
Control.KO5.WT <- read.csv("./overlap/Control.KO5.WT_significant_reexpressed.csv")
Control.LacZ.WT <- read.csv("./overlap/Control.LacZ.WT_significant_reexpressed.csv")
Control.KO5.LacZ.WT <- merge(Control.KO5.WT, Control.LacZ.WT, by=c("geneID"))
Control.KO5.LacZ.WT$sorter <- -1*log(Control.KO5.LacZ.WT$sorter.x*Control.KO5.LacZ.WT$sorter.y)
Control.KO5.LacZ.WT <- Control.KO5.LacZ.WT[order(Control.KO5.LacZ.WT$sorter),]
write.csv(Control.KO5.LacZ.WT, "./overlap/Control.KO5.LacZ.WT.csv")

# GSEA----
hs_gsea <- msigdbr(species = "Homo sapiens")
hs_gsea %>% 
  dplyr::distinct(gs_cat, gs_subcat) %>% 
  dplyr::arrange(gs_cat, gs_subcat)
hs_gsea_h <- msigdbr(species = "Homo sapiens",
                     category = "H") %>% # H for Hallmark dataset; C5 for gene ontology; C6 cor oncogenic
  dplyr::select(gs_name, gene_symbol)

## GSEA - Control to KO5----
mydata.df.Control.KO5 <- dplyr::select(mydata.df, geneID, Control.KO5)
Control.KO5.gsea <- mydata.df.Control.KO5$Control.KO5
names(Control.KO5.gsea) <- as.character(mydata.df.Control.KO5$geneID)
Control.KO5.gsea <- sort(Control.KO5.gsea, decreasing = TRUE)

Control.KO5.res <- GSEA(Control.KO5.gsea, pvalueCutoff = 1, TERM2GENE=hs_gsea_h, verbose=FALSE)
ControlKO5.GSEA.df <- as_tibble(Control.KO5.res@result)
write_csv(ControlKO5.GSEA.df, "./gsea/Control.KO5.GSEA.csv")
mydata.df.Control.KO5 <- as_tibble(Control.KO5.res@result)

mydata.df.Control.KO5 <- mydata.df.Control.KO5 %>%
  mutate(phenotype = case_when(
    NES > 0 ~ "KO",
    NES < 0 ~ "Control"))
#mydata.df.Control.KO5 <- mydata.df.Control.KO5[order(abs(mydata.df.Control.KO5$NES)),]

Control.KO5.GSEA <- ggplot(mydata.df.Control.KO5[1:20,], aes(x=NES, y=ID)) + 
  geom_point(aes(size=setSize, color = p.adjust)) +
  scale_color_gradient(low="black", high="grey", limits=c(0,.05), na.value="white") +
  theme_fivethirtyeight()
ggsave(path = "./plots/gsea", filename = "Control.KO5.GSEA.png", plot = Control.KO5.GSEA)

## GSEA - LacZ to WT----
mydata.df.LacZ.WT <- dplyr::select(mydata.df, geneID, LacZ.WT)
LacZ.WT.gsea <- mydata.df.LacZ.WT$LacZ.WT
names(LacZ.WT.gsea) <- as.character(mydata.df.LacZ.WT$geneID)
LacZ.WT.gsea <- sort(LacZ.WT.gsea, decreasing = TRUE)

LacZ.WT.res <- GSEA(LacZ.WT.gsea, pvalueCutoff = 1, TERM2GENE=hs_gsea_h, verbose=FALSE)
LacZ.WT.GSEA.df <- as_tibble(LacZ.WT.res@result)
write_csv(LacZ.WT.GSEA.df, "./gsea/LacZ.WT.GSEA.csv")
mydata.df.LacZ.WT <- as_tibble(LacZ.WT.res@result)

mydata.df.LacZ.WT <- mydata.df.LacZ.WT %>%
  mutate(phenotype = case_when(
    NES > 0 ~ "WT",
    NES < 0 ~ "LacZ"))
#mydata.df.LacZ.WT <- mydata.df.LacZ.WT[order(abs(mydata.df.LacZ.WT$NES)),]

LacZ.WT.GSEA <- ggplot(mydata.df.LacZ.WT[1:20,], aes(x=NES, y=ID)) + 
  geom_point(aes(size=setSize, color = p.adjust)) +
  scale_color_gradient(low="black", high="grey", limits=c(0,.05), na.value="white") +
  theme_fivethirtyeight()
ggsave(path = "./plots/gsea", filename = "LacZ.WT.GSEA.png", plot = LacZ.WT.GSEA)

## GSEA - KO5 to WT----
mydata.df.KO5.WT <- dplyr::select(mydata.df, geneID, KO5.WT)
KO5.WT.gsea <- mydata.df.KO5.WT$KO5.WT
names(KO5.WT.gsea) <- as.character(mydata.df.KO5.WT$geneID)
KO5.WT.gsea <- sort(KO5.WT.gsea, decreasing = TRUE)

KO5.WT.res <- GSEA(KO5.WT.gsea, pvalueCutoff = 1, TERM2GENE=hs_gsea_h, verbose=FALSE)
ControlKO5.GSEA.df <- as_tibble(KO5.WT.res@result)
write_csv(ControlKO5.GSEA.df, "./gsea/KO5.WT.GSEA.csv")
mydata.df.KO5.WT <- as_tibble(KO5.WT.res@result)

mydata.df.KO5.WT <- mydata.df.KO5.WT %>%
  mutate(phenotype = case_when(
    NES > 0 ~ "WT",
    NES < 0 ~ "KO5"))
#mydata.df.KO5.WT <- mydata.df.KO5.WT[order(abs(mydata.df.KO5.WT$NES)),]

KO5.WT.GSEA <- ggplot(mydata.df.KO5.WT[1:20,], aes(x=NES, y=ID)) + 
  geom_point(aes(size=setSize, color = p.adjust)) +
  scale_color_gradient(low="black", high="grey", limits=c(0,.05), na.value="white") +
  theme_fivethirtyeight()
ggsave(path = "./plots/gsea", filename = "KO5.WT.GSEA.png", plot = KO5.WT.GSEA)

## GSEA - Control to E6----
mydata.df.Control.E6 <- dplyr::select(mydata.df, geneID, Control.E6)
Control.E6.gsea <- mydata.df.Control.E6$Control.E6
names(Control.E6.gsea) <- as.character(mydata.df.Control.E6$geneID)
Control.E6.gsea <- sort(Control.E6.gsea, decreasing = TRUE)

Control.E6.res <- GSEA(Control.E6.gsea, pvalueCutoff = 1, TERM2GENE=hs_gsea_h, verbose=FALSE)
ControlE6.GSEA.df <- as_tibble(Control.E6.res@result)
write_csv(ControlE6.GSEA.df, "./gsea/Control.E6.GSEA.csv")
mydata.df.Control.E6 <- as_tibble(Control.E6.res@result)

mydata.df.Control.E6 <- mydata.df.Control.E6 %>%
  mutate(phenotype = case_when(
    NES > 0 ~ "E6",
    NES < 0 ~ "Control"))
#mydata.df.Control.E6 <- mydata.df.Control.E6[order(abs(mydata.df.Control.E6$NES)),]

Control.E6.GSEA <- ggplot(mydata.df.Control.E6[1:20,], aes(x=NES, y=ID)) + 
  geom_point(aes(size=setSize, color = p.adjust)) +
  scale_color_gradient(low="black", high="grey", limits=c(0,.05), na.value="white") +
  theme_fivethirtyeight()
ggsave(path = "./plots/gsea", filename = "Control.E6.GSEA.png", plot = Control.E6.GSEA)
