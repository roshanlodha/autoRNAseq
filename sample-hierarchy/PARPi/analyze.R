# load packages----
setRepositories(ind = c(1:8))

library(tidyverse) # provides access to Hadley Wickham's collection of R packages for data science, which we will use throughout the course
library(tximport) # package for getting Kallisto results into R
library(ensembldb) #helps deal with ensembl
library(EnsDb.Hsapiens.v86) #replace with your organism-specific database package

library(edgeR) # well known package for differential expression analysis, but we only use for the DGEList object and for normalization methods
library(matrixStats) # let's us easily calculate stats on rows or columns of a data matrix
library(cowplot) # allows you to combine multiple plots in one figure

library(DT) # for making interactive tables
library(plotly) # for making interactive plots
library(gt) # A layered 'grammar of tables' - think ggplot, but for tables

library(limma) # venerable package for differential gene expression using linear modeling

# transcripts import----
targets <- read_tsv("studydesign.txt")# read in your study design
paths <- file.path(targets$dir, "abundance.tsv") # set file paths to your mapped data
Tx <- transcripts(EnsDb.Hsapiens.v86, columns=c("tx_id", "gene_name"))
Tx <- as_tibble(Tx)
Tx <- dplyr::rename(Tx, target_id = tx_id)
Tx <- dplyr::select(Tx, "target_id", "gene_name")
Txi_gene <- tximport(paths, 
                     type = "kallisto", 
                     tx2gene = Tx, 
                     txOut = FALSE, #determines whether your data represented at transcript or gene level
                     countsFromAbundance = "lengthScaledTPM",
                     ignoreTxVersion = TRUE)

# log2 CPM unfiltered, non-normalized----
sampleLabels <- targets$sample
myDGEList <- DGEList(Txi_gene$counts)
log2.cpm <- cpm(myDGEList, log=TRUE)

log2.cpm.df <- as_tibble(log2.cpm, rownames = "geneID")
colnames(log2.cpm.df) <- c("geneID", sampleLabels)
log2.cpm.df.pivot <- pivot_longer(log2.cpm.df, # dataframe to be pivoted
                                  cols = sampleLabels[1]:sampleLabels[length(sampleLabels)], # column names to be stored as a SINGLE variable
                                  names_to = "samples", # name of that new variable (column)
                                  values_to = "expression") # name of new variable (column) storing all the values (data)

p1 <- ggplot(log2.cpm.df.pivot) +
  aes(x=samples, y=expression, fill=samples) +
  geom_violin(trim = FALSE, show.legend = FALSE) +
  stat_summary(fun = "median", 
               geom = "point", 
               shape = 95, 
               size = 10, 
               color = "black", 
               show.legend = FALSE) +
  labs(y="log2 expression", x = "sample",
       title="Log2 Counts per Million (CPM)",
       subtitle="unfiltered, non-normalized",
       caption=paste0("produced on ", Sys.time())) +
  theme_bw()

# log2 CPM filtered, non-normalized----
cpm <- cpm(myDGEList)
keepers <- rowSums(cpm>1)>=(dim(targets)[1]/2) #user defined
myDGEList.filtered <- myDGEList[keepers,]

log2.cpm.filtered <- cpm(myDGEList.filtered, log=TRUE)
log2.cpm.filtered.df <- as_tibble(log2.cpm.filtered, rownames = "geneID")
colnames(log2.cpm.filtered.df) <- c("geneID", sampleLabels)
log2.cpm.filtered.df.pivot <- pivot_longer(log2.cpm.filtered.df, # dataframe to be pivoted
                                           cols = sampleLabels[1]:sampleLabels[length(sampleLabels)], # column names to be stored as a SINGLE variable
                                           names_to = "samples", # name of that new variable (column)
                                           values_to = "expression") # name of new variable (column) storing all the values (data)

p2 <- ggplot(log2.cpm.filtered.df.pivot) +
  aes(x=samples, y=expression, fill=samples) +
  geom_violin(trim = FALSE, show.legend = FALSE) +
  stat_summary(fun = "median", 
               geom = "point", 
               shape = 95, 
               size = 10, 
               color = "black", 
               show.legend = FALSE) +
  labs(y="log2 expression", x = "sample",
       title="Log2 Counts per Million (CPM)",
       subtitle="filtered, non-normalized",
       caption=paste0("produced on ", Sys.time())) +
  theme_bw()

# log2 CPM filtered, normalized----
myDGEList.filtered.norm <- calcNormFactors(myDGEList.filtered, method = "TMM")
log2.cpm.filtered.norm <- cpm(myDGEList.filtered.norm, log=TRUE)
log2.cpm.filtered.norm.df <- as_tibble(log2.cpm.filtered.norm, rownames = "geneID")
colnames(log2.cpm.filtered.norm.df) <- c("geneID", sampleLabels)
log2.cpm.filtered.norm.df.pivot <- pivot_longer(log2.cpm.filtered.norm.df, # dataframe to be pivoted
                                                cols = sampleLabels[1]:sampleLabels[length(sampleLabels)], # column names to be stored as a SINGLE variable
                                                names_to = "samples", # name of that new variable (column)
                                                values_to = "expression") # name of new variable (column) storing all the values (data)

p3 <- ggplot(log2.cpm.filtered.norm.df.pivot) +
  aes(x=samples, y=expression, fill=samples) +
  geom_violin(trim = FALSE, show.legend = FALSE) +
  stat_summary(fun = "median", 
               geom = "point", 
               shape = 95, 
               size = 10, 
               color = "black", 
               show.legend = FALSE) +
  labs(y="log2 expression", x = "sample",
       title="Log2 Counts per Million (CPM)",
       subtitle="filtered, TMM normalized",
       caption=paste0("produced on ", Sys.time())) +
  theme_bw()

plot_grid(p1, p2, p3, labels = c('A', 'B', 'C'), label_size = 12)

# multivariate analysis----
group <- targets$group
group <- factor(group)

## distance mapping----
#distance <- dist(t(log2.cpm.filtered.norm), method = "maximum") #other distance methods are "euclidean", maximum", "manhattan", "canberra", "binary" or "minkowski"
#clusters <- hclust(distance, method = "average") #other agglomeration methods are "ward.D", "ward.D2", "single", "complete", "average", "mcquitty", "median", or "centroid"
#plot(clusters, labels=sampleLabels)

## pca analysis----
pca.res <- prcomp(t(log2.cpm.filtered.norm), scale.=F, retx=T)
pc.var <- pca.res$sdev^2 # sdev^2 captures these eigenvalues from the PCA result
pc.per <- round(pc.var/sum(pc.var)*100, 1) 
pca.res.df <- as_tibble(pca.res$x)
pca.plot <- ggplot(pca.res.df) +
  aes(x=PC1, y=PC2, label=sampleLabels, color = group) +
  geom_point(size=4) +
  geom_label() +
  stat_ellipse() +
  xlab(paste0("PC1 (",pc.per[1],"%",")")) + 
  ylab(paste0("PC2 (",pc.per[2],"%",")")) +
  labs(title="PCA plot",
       caption=paste0("produced on ", Sys.time())) +
  coord_fixed() +
  theme_bw()

## average analysis----
mydata.df <- mutate(log2.cpm.filtered.norm.df,
                    naive.AVG = (LTL331sample1 + LTL331sample2)/2, 
                    resistant.AVG = (LTL331Rsample3 + LTL331Rsample4)/2,
                    #now make columns comparing each of the averages above that you're interested in
                    LogFC = (resistant.AVG - naive.AVG)) %>% #note that this is the first time you've seen the 'pipe' operator
  mutate_if(is.numeric, round, 2)

# differential gene analysis----
group <- factor(targets$group)
design <- model.matrix(~0 + group)
colnames(design) <- levels(group)

v.DEGList.filtered.norm <- voom(myDGEList.filtered.norm, design, plot = FALSE)
fit <- lmFit(v.DEGList.filtered.norm, design)
contrast.matrix <- makeContrasts(infection = resistant - naive,
                                 levels=design)
fits <- contrasts.fit(fit, contrast.matrix)
ebFit <- eBayes(fits)
myTopHits <- topTable(ebFit, adjust ="BH", coef=1, number=40000, sort.by="logFC")
myTopHits.df <- myTopHits %>%
  as_tibble(rownames = "geneID")
vplot <- ggplot(myTopHits) +
  aes(y=-log2(adj.P.Val), x=logFC) +
  geom_point(size=2) +
  #geom_hline(yintercept = -log2(0.05), linetype="longdash", colour="grey", size=1) +
  #geom_vline(xintercept = 1, linetype="longdash", colour="#BE684D", size=1) +
  #geom_vline(xintercept = -1, linetype="longdash", colour="#2C467A", size=1) +
  annotate("rect", xmin = 1, xmax = 10, ymin = -log2(0.05), ymax = 6, alpha=.2, fill="#BE684D") +
  annotate("rect", xmin = -1, xmax = -10, ymin = -log2(0.05), ymax = 6, alpha=.2, fill="#2C467A") +
  labs(title="Volcano plot",
       subtitle = "naive to resistant",
       caption=paste0("produced on ", Sys.time())) +
  theme_bw()
vplot

write_csv(myTopHits.df, "./logFC.csv")

results <- decideTests(ebFit, method="global", adjust.method="BH", p.value=0.05, lfc=1)
colnames(v.DEGList.filtered.norm$E) <- sampleLabels
diffGenes <- v.DEGList.filtered.norm$E[results[,1] !=0,]
diffGenes.df <- as_tibble(diffGenes, rownames = "geneID")
