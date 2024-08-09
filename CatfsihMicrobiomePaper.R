## Install and load needed packages

install.packages("BiocManager")
BiocManager::install("phyloseq", force = TRUE)
BiocManager::install("phyloseq")
BiocManager::install("Biostrings")
BiocManager::install("metagenomeSeq" = TRUE)
install.packages("tidyverse")
installed.packages("ggplot2")
installed.packages("ggpubr")
BiocManager::install("decontam", force = TRUE)
BiocManager::install("indicspecies", force = TRUE)
BiocManager::install("ggpubr", force = TRUE)
BiocManager::install("DESeq2" = TRUE)
install.packages("pairwiseAdonis", force = TRUE)


library(phyloseq)
library(vegan)
library(tidyverse)
library(ggplot2)
library(Biostrings)
library(ggpubr)
library(decontam)
library(metagenomeSeq)
library(indicspecies)
library(dada2); packageVersion("dada2")
library(DESeq2)
library(ggrepel)
library(pairwiseAdonis)

# color blind pallet
cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")


#load forward and reverse reads
#select file directory
getwd()
path <- "/Users/abdulmalikoladipupo/Downloads/qimme2 data/ccf trial 4/dada2/zr15092.rawdata.231213"
list.files(path)

# Forward and reverse fastq filenames have format: ID_SAMPLENAME_R1_001.fastq.gz and ID_SAMPLENAME_R2_001.fastq.gz
fnFs <- sort(list.files(path, pattern="_R1.fastq.gz", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_R2.fastq.gz", full.names = TRUE))

# Extract sample names, assuming filenames have format: ID_SAMPLENAME_XXX.fastq
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 2)

#plot quality of foward and reverse
plotQualityProfile(fnFs[1:4])
plotQualityProfile(fnRs[1:4])

# Place filtered files in filtered/ subdirectory
filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))
names(filtFs) <- sample.names
names(filtRs) <- sample.names

# trim forward and reverse reads based on quality and trim primer length
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(320,180), trimLeft = c(16,24),
                     maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
                     compress=TRUE, multithread=TRUE, matchIDs = TRUE) # On Windows set multithread=FALSE

head(out)

# learning error rate
errF <- learnErrors(filtFs, multithread=TRUE)
errR <- learnErrors(filtRs, multithread=TRUE)
plotErrors(errF, nominalQ=TRUE)

#dereplication
derepFs <- derepFastq(filtFs, verbose=TRUE)
derepRs <- derepFastq(filtRs, verbose=TRUE)
# Name the derep-class objects by the sample names
names(derepFs) <- sample.names
names(derepRs) <- sample.names

dadaFs <- dada(derepFs, err=errF, multithread=TRUE)
dadaRs <- dada(derepRs, err=errR, multithread=TRUE)
dadaFs[[2]]

#merge paired reads
mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose=TRUE)
# Inspect the merger data.frame from the first sample
head(mergers[[1]])

# construct asv table
seqtab <- makeSequenceTable(mergers)
dim(seqtab)


# Inspect distribution of sequence lengths
table(nchar(getSequences(seqtab)))

#remove chimera
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim)
sum(seqtab.nochim)/sum(seqtab)

#track the pipeline
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
# If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
head(track)

#load reference database classifier
getwd()
#assign taxonomy
taxa <- assignTaxonomy(seqtab.nochim, "/Users/abdulmalikoladipupo/Downloads/qimme2 data/ccf trial 4/dada2/zr15092.rawdata.231213/silva_nr99_v138.1_train_set.fa.gz", multithread=TRUE)
taxa <- addSpecies(taxa, "/Users/abdulmalikoladipupo/Downloads/qimme2 data/ccf trial 4/dada2/zr15092.rawdata.231213/silva_species_assignment_v138.1.fa.gz")

#inspect taxonomic assignment
taxa.print <- taxa # Removing sequence rownames for display only
rownames(taxa.print) <- NULL
head(taxa.print)

# Convert taxa.print to a data frame
taxa_df <- data.frame(taxa.print)
head(taxa_df)


#constructing phylogenetic tree
install.packages("BiocStyle") 
## load packages
library("knitr")
library("BiocStyle")
.cran_packages <- c("ggplot2", "gridExtra")
.bioc_packages <- c("dada2", "phyloseq", "DECIPHER", "phangorn")
.inst <- .cran_packages %in% installed.packages()
if(any(!.inst)) {
  install.packages(.cran_packages[!.inst])
}
.inst <- .bioc_packages %in% installed.packages()
if(any(!.inst)) {
  source("http://bioconductor.org/biocLite.R")
  biocLite(.bioc_packages[!.inst], ask = F)
}
# Load packages into session, and print package version
sapply(c(.cran_packages, .bioc_packages), require, character.only = TRUE)

seqs <- getSequences(seqtab.nochim)
names(seqs) <- seqs # This propagates to the tip labels of the tree
alignment <- AlignSeqs(DNAStringSet(seqs), anchor=NA,verbose=FALSE)

## using the package phangorn to build the tree using maximum likelihood
phangAlign <- phyDat(as(alignment, "matrix"), type="DNA")
dm <- dist.ml(phangAlign)
treeNJ <- NJ(dm) # Note, tip order != sequence order
fit = pml(treeNJ, data=phangAlign)
fitGTR <- update(fit, k=4, inv=0.2)
fitGTR <- optim.pml(fitGTR, model="GTR", optInv=TRUE, optGamma=TRUE,
                    rearrangement = "stochastic", control = pml.control(trace = 0))
detach("package:phangorn", unload=TRUE)



### Handing off to phyloseq
## load meta data
samdf <- read.csv("/Users/abdulmalikoladipupo/Downloads/qimme2 data/ccf trial 4/dada2/zr15092.rawdata.231213/metadata.csv",header=TRUE)
samples.out <- rownames(seqtab.nochim)
rownames(samdf) <- samples.out

## fix all files into phyloseq (metadata, otu table, phy tree, taxa table)
ps <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows=FALSE), 
               sample_data(samdf), 
               tax_table(taxa),phy_tree(fitGTR$tree))
ps

# removing chloroplast or taxa not assigned at the domain level
physeq.clean <- ps %>% subset_taxa(Family!= "Mitochondria" & Order!= "Chloroplast")
physeq.clean

# Remove samples with less than 1000 reads
physeq.no.chloro <- prune_samples(sample_sums(physeq.clean) >= 1000, physeq.clean)
physeq.no.chloro

# Now lets look at the read distribution per sample and decide if we need to get rid of some samples because of low sequence depth
sort(sample_sums(physeq.no.chloro), decreasing = T) # read distribution

# how many total reads are we now working with?
sum(sample_sums(physeq.no.chloro))

# What is our mean and median read depth per sample? 
mean(sample_sums(physeq.no.chloro))
median(sample_sums(physeq.no.chloro))

#Lets make a histogram of the read distribution and put the median read depth
read.depths <- data.frame(sample_sums(physeq.no.chloro))
colnames(read.depths) <- "read.depth"
read.depth.plot <- ggplot(read.depths, aes(read.depth)) +
  geom_histogram(fill = cbbPalette[[3]], color = "black") + 
  geom_vline(xintercept = median(sample_sums(physeq.no.chloro)), linetype = "dashed") + 
  theme_classic() + 
  xlab("Read Depth")


# Rarefaction analysis 
sam.data <- data.frame(physeq.no.chloro@sam_data)
physeq.no.chloro@otu_table
bOTU.table <- otu_table(physeq.no.chloro) %>%
  as.data.frame() %>%
  as.matrix()

raremax <- min(rowSums(bOTU.table)) ## without the "t" tranform
rare.fun <- rarecurve((bOTU.table), step = 1000, sample = raremax, tidy = T) ## without the "t" tranform

bac.rare.curve.extract2 <- left_join(sam.data, rare.fun, by = c("SampleID" = "Site"))
bac.rare <- ggplot(bac.rare.curve.extract2, aes(x = Sample, y = Species, group = SampleID, color = treatment)) + 
  #geom_point() +
  geom_line() + 
  xlab("Reads") + 
  ylab("Number of OTUs") +
  ggtitle("") +
  theme_classic() + 
  geom_vline(xintercept = median(sample_sums(physeq.no.chloro)), linetype = "dashed")
#scale_color_manual(values = cbbPalette)

# publication ready figure
SuplementalFig1 <- ggarrange(bac.rare, nrow = 1, labels = c(""))
ggsave("SupplementalFig1.tiff", dpi=300, width = 15, height = 10, units = "cm")


# Normalize Sampling reads based on cumulative sum scaling (CSS normalization)
MGS <- phyloseq_to_metagenomeSeq(physeq.no.chloro)
p <- metagenomeSeq::cumNormStatFast(MGS)
MGS <- metagenomeSeq::cumNorm(MGS, p =p)
metagenomeSeq::normFactors(MGS) # exports the normalized factors for each sample
norm.bacteria <- metagenomeSeq::MRcounts(MGS, norm = T)
norm.bacteria.OTU <- phyloseq::otu_table(norm.bacteria, taxa_are_rows = TRUE)

# fix back the normalized OTU table into phyloseq and save as new object
physeq.css <- phyloseq(otu_table(norm.bacteria.OTU, taxa_are_rows=FALSE), 
                       sample_data(samdf), 
                       tax_table(taxa),phy_tree(fitGTR$tree))


## Analysis

# ALPHA DIVERSITY
physeq.no.chloro@sam_data$shannon <- estimate_richness(physeq.no.chloro, measures=c("Shannon"))$Shannon
physeq.no.chloro@sam_data$invsimpson <- estimate_richness(physeq.no.chloro, measures=c("InvSimpson"))$InvSimpson
physeq.no.chloro@sam_data$richness <- estimate_richness(physeq.no.chloro, measures=c("Observed"))$Observed
physeq.no.chloro@sam_data$even <- physeq.no.chloro@sam_data$shannon/log(physeq.no.chloro@sam_data$richness)

# convert to a data frame
sample.data.bac <- data.frame(physeq.no.chloro@sam_data)

# Richness by treatment*time interaction plot
richness.treatment.time <- ggplot(sample.data.bac, aes(x = treatment, y = richness, fill = treatment)) + # fill = treatment
  geom_boxplot() +
  #geom_jitter() + 
  ylab("Richness") + 
  stat_compare_means(method = "") + 
  xlab("")+
  theme_bw() +
  facet_wrap(~time) +
  theme(axis.text.x = element_text(size = 11, vjust = 0.5, hjust=0.5),
        axis.text.y = element_text(size = 11),
        axis.title.x = element_text(size = 14, face = "bold"),
        axis.title.y = element_text(size = 14, face = "bold"))

richness.treatment.time

# richness statistics
rich_trt_time <- aov(richness ~ treatment+time+treatment:time, data = sample.data.bac)
summary(rich_trt_time)
TukeyHSD(rich_trt_time)
rich_trt_time$residuals
shapiro.test(rich_trt_time$residuals)

anova_table <- anova(rich_trt_time)
print(anova_table)
tukey_results <- TukeyHSD(rich_trt_time)
print(tukey_results)

# Shannon diversity by time*treatment interaction plot
shan.treatment.time <- ggplot(sample.data.bac2, aes(x = treatment, y = shannon, fill = treatment)) + 
  geom_boxplot() +
  #geom_jitter() + 
  ylab("Shannon") + 
  stat_compare_means(method = "") + 
  xlab("")+
  theme_bw() +
  facet_wrap(~time) +
  theme(axis.text.x = element_text(size = 11, vjust = 0.5, hjust=0.5),
        axis.text.y = element_text(size = 11),
        axis.title.y = element_text(size = 14, face = "bold"))

shan.treatment.time

# shannon statistics
shan_trt_time <- aov(shannon ~ treatment+time+treatment:time, data = sample.data.bac)
anova_table <- anova(shan_trt_time)
print(anova_table)
tukey_results <- TukeyHSD(shan_trt_time)
print(tukey_results)


# even diversity by time*treatment interaction plot
even.treatment.time <- ggplot(sample.data.bac, aes(x = treatment, y = even, fill = treatment)) + 
  geom_boxplot() +
  #geom_jitter() + 
  ylab("Eveness") + 
  stat_compare_means(method = "") + 
  xlab("")+
  theme_bw() +
  facet_wrap(~time) +
  theme(axis.text.x = element_text(size = 11, vjust = 0.5, hjust=0.5),
        axis.text.y = element_text(size = 11),
        axis.title.y = element_text(size = 14, face = "bold"))

even.treatment.time

# eveness statistics
even_trt_time <- aov(even ~ treatment+time+treatment:time, data = sample.data.bac)
anova_table <- anova(even_trt_time)
print(anova_table)
tukey_results <- TukeyHSD(even_trt_time)
print(tukey_results)

#export figures
SuplementalFig2 <- ggarrange(richness.treatment.time, nrow = 1, labels = c(""))
ggsave("SupplementalFig2.tiff", dpi=300, width = 15, height = 10, units = "cm")

SuplementalFig3 <- ggarrange(even.treatment.time, nrow = 1, labels = c(""))
ggsave("SupplementalFig3.tiff", dpi=300, width = 15, height = 10, units = "cm")

SuplementalFig4 <- ggarrange(shan.treatment.time, nrow = 1, labels = c(""))
ggsave("SupplementalFig4.tiff", dpi=300, width = 15, height = 10, units = "cm")



## BETA DIVERSITY
# Principle coordinates analysis with Bray-Curtis distances ## 30 days data
ordination.pcoa <- ordinate(phy30, "PCoA", "bray") # calculate the resemblance and ordinate using PCoA
ordination.pcoa$vectors # positions of your points on the pcoa graph
ordination.pcoa$values #values to calculate the variance explained on each axis (dimension)

pcoa <- plot_ordination(phy30, ordination = ordination.pcoa, type = "samples", color = "treatment") +
  theme_classic() + 
  scale_color_manual(values = cbbPalette)
pcoa

pcoa.data <- pcoa$data # taking the data to make a fancy plot
pcoa30.plot.2 <- ggplot(pcoa.data, aes(x = Axis.1, y = Axis.2, fill = treatment)) +
  geom_point(alpha = 0.8, size = 2, aes(color = treatment)) +
  theme_bw() +
  ylab("PCoA2 (19.4%)") + 
  xlab("PCoA1 (33.2%)") +
  stat_ellipse(geom = "polygon", alpha = 0.05, linewidth = 0.5, aes(color = treatment))+
  scale_fill_manual(values = cbbPalette) +
  scale_shape_manual(values = c(21, 24, 23)) +
  facet_wrap(~time) +
  theme(axis.text.x = element_text(size = 11, vjust = 0.5, hjust=0.5),
        axis.text.y = element_text(size = 11),
        axis.title.y = element_text(size = 14, face = "bold"),
        axis.title.x = element_text(size = 14, face = "bold"))

pcoa30.plot.2

# Principle coordinates analysis with Bray-Curtis distances ## 60 days data
ordination.pcoa <- ordinate(phy60, "PCoA", "bray") # calculate the resemblance and ordinate using PCoA
ordination.pcoa$vectors # positions of your points on the pcoa graph
ordination.pcoa$values #values to calculate the variance explained on each axis (dimension)

pcoa <- plot_ordination(phy60, ordination = ordination.pcoa, type = "samples", color = "treatment") +
  theme_classic() + 
  scale_color_manual(values = cbbPalette)
pcoa

pcoa.data <- pcoa$data # taking the data to make a fancy plot
pcoa60.plot.2 <- ggplot(pcoa.data, aes(x = Axis.1, y = Axis.2, fill = treatment)) +
  geom_point(alpha = 0.8, size = 2, aes(color = treatment)) +
  theme_bw() +
  ylab("PCoA2 (20.5%)") + 
  xlab("PCoA1 (40.5%)") +
  stat_ellipse(geom = "polygon", alpha = 0.05, linewidth = 0.5, aes(color = treatment))+
  scale_fill_manual(values = cbbPalette) +
  scale_shape_manual(values = c(21, 24, 23)) +
  facet_wrap(~time) +
  #guides(fill = guide_legend(override.aes = list(shape = 21)))
  theme(axis.text.x = element_text(size = 11, vjust = 0.5, hjust=0.5),
        axis.text.y = element_text(size = 11),
        axis.title.y = element_text(size = 14, face = "bold"),
        axis.title.x = element_text(size = 14, face = "bold"))

pcoa60.plot.2

SuplementalFig6c <- ggarrange(pcoa30.plot.2, nrow = 1, labels = c())
ggsave("SupplementalFig6c.tiff", dpi=300, width = 15, height = 10, units = "cm")

SuplementalFig6a <- ggarrange(pcoa60.plot.2, nrow = 1, labels = c())
ggsave("SupplementalFig6a.tiff", dpi=300, width = 15, height = 10, units = "cm")


# Principle coordinates analysis with Weighted-unifrac distances 30 days (plot)
phy30@phy_tree
unifrac <- UniFrac(physeq.css, weighted = FALSE)
ordination.unifrac.pcoa <- ordinate(physeq.css, "PCoA", distance = unifrac) # calculate the resemblance and ordinate using PCoA
ordination.unifrac.pcoa$vectors # positions of your points on the pcoa graph
ordination.unifrac.pcoa$values #values to calculate the variance explained on each axis (dimension)

pcoa.unifrac <- plot_ordination(physeq.css, ordination = ordination.unifrac.pcoa, type = "samples", color = "treatment", shape = "time") +
  theme_bw() +
  theme(axis.text.x = element_text(size = 11, vjust = 0.5, hjust=0.5),
        axis.text.y = element_text(size = 11),
        axis.title.y = element_text(size = 14, face = "bold"),
        axis.title.x = element_text(size = 14, face = "bold"))
#scale_color_manual(values = cbbPalette) + 
#stat_ellipse()
pcoa.unifrac

unifrac.data <- pcoa.unifrac$data # taking the data to make a fancy plot
unifrac.plot <- ggplot(unifrac.data, aes(x = Axis.1, y = Axis.2, fill = treatment)) +
  geom_point(alpha = 0.8, size = 2, aes(color = treatment)) +
  theme_bw() +
  ylab("PCoA2 (10.4%)") + 
  xlab("PCoA1 (26.5%)") +
  stat_ellipse(geom = "polygon", alpha = 0.05, linewidth = 0.5, aes(color = treatment))+
  scale_fill_manual(values = cbbPalette) +
  scale_shape_manual(values = c(21, 24, 23)) +
  facet_wrap(~time) +
  #guides(fill = guide_legend(override.aes = list(shape = 21)))
  theme(axis.text.x = element_text(size = 11, vjust = 0.5, hjust=0.5),
        axis.text.y = element_text(size = 11),
        axis.title.y = element_text(size = 14, face = "bold"),
        axis.title.x = element_text(size = 14, face = "bold"))   


SuplementalFig7a <- ggarrange(unifrac.plot, nrow = 1, labels = c())
ggsave("SupplementalFig7a.tiff", dpi=300, width = 20, height = 10, units = "cm")


# PERMANOVA - testing for differences in centroids
# bray-curtis
set.seed(2000)
prok.dist.bray = phyloseq::distance(physeq.css, "bray") # create bray-curtis distance matrix
prok.dist.bray30 = phyloseq::distance(phy30, "bray") # create bray-curtis distance matrix
prok.dist.bray60 = phyloseq::distance(phy60, "bray") # create bray-curtis distance matrix

result.bray <- adonis2(prok.dist.bray~treatment*time, as(sample_data(physeq.css), "data.frame")) #Are there significant changes? Significant intercation? Yes
result.bray30 <- adonis2(prok.dist.bray30~treatment, as(sample_data(phy30), "data.frame")) #Are there significant changes?
result.bray60 <- adonis2(prok.dist.bray60~treatment, as(sample_data(phy60), "data.frame")) #Are there significant changes?

#pairwise comparison
pairwise.adonis2(prok.dist.bray30~treatment, as(sample_data(phy30), "data.frame"))
pairwise.adonis2(prok.dist.bray60~treatment, as(sample_data(phy60), "data.frame"))


unifrac <- UniFrac(physeq.css, weighted = FALSE)
unifrac2 <- UniFrac(phy30, weighted = FALSE)
unifrac3 <- UniFrac(phy60, weighted = FALSE)

# weight unifrac
set.seed(2000)
result.unifrac <- adonis2(unifrac~treatment*time, as(sample_data(physeq.css), "data.frame")) #Are there significant changes?
result.unifrac30 <- adonis2(unifrac2~treatment, as(sample_data(phy30), "data.frame")) #Are there significant changes?
result.unifrac60 <- adonis2(unifrac3~treatment, as(sample_data(phy60), "data.frame")) #Are there significant changes?

result.unifrac
result.unifrac30
result.unifrac60

pairwise.adonis2(unifrac~treatment*time, as(sample_data(physeq.css), "data.frame"))
pairwise.adonis2(unifrac2~treatment, as(sample_data(phy30), "data.frame"))
pairwise.adonis2(unifrac3~treatment, as(sample_data(phy60), "data.frame"))

## END


unifrac2 <- UniFrac(phy30, weighted = FALSE)
unifrac3 <- UniFrac(phy60, weighted = FALSE)
# Unweight unifrac
set.seed(2000)
result.unifrac <- adonis2(unifrac~treatment*time, as(sample_data(physeq.css), "data.frame")) #Are there significant changes?
result.unifrac30 <- adonis2(unifrac2~treatment, as(sample_data(phy30), "data.frame")) #Are there significant changes?
result.unifrac60 <- adonis2(unifrac3~treatment, as(sample_data(phy60), "data.frame")) #Are there significant changes?

result.unifrac
result.unifrac30
result.unifrac60

pairwise.adonis2(unifrac~treatment + time, as(sample_data(physeq.css), "data.frame"))
pairwise.adonis2(unifrac2~treatment, as(sample_data(phy30), "data.frame"))
pairwise.adonis2(unifrac3~treatment, as(sample_data(phy60), "data.frame"))


### Relative Abundance plot ###

# fecal samples relative amount/treatment*time
top20 <- names(sort(taxa_sums(physeq.no.chloro), decreasing=TRUE))[1:30]
ps.top20 <- transform_sample_counts(physeq.no.chloro, function(OTU) OTU/sum(OTU))
ps.top20 <- prune_taxa(top20, ps.top20)
dat.dataframe <- psmelt(ps.top20)
dat.agr <- aggregate(Abundance ~ treatment*time + Phylum, data = dat.dataframe, FUN = mean)
rel_abun_treatment <-ggplot(dat.agr, aes(x = treatment, y = Abundance, fill = Phylum)) +
  facet_wrap(~time) +
  geom_bar(stat = "identity", position = "fill") +
  labs(title = "") +
  xlab("Dietary treatments") +
  ylab("Relative abundance (%)") +
  theme(axis.text.x = element_text(size = 11, vjust = 0.5, hjust=0.5),
        axis.text.y = element_text(size = 11),
        axis.title.y = element_text(size = 11, face = "bold"),
        axis.title.x = element_text(size = 11, face = "bold"))   


dat.agr.GEN <- aggregate(Abundance ~ treatment*time + Genus, data = dat.dataframe, FUN = mean)
rel_abun_treatment.Gen <-ggplot(dat.agr.GEN, aes(x = treatment, y = Abundance, fill = Genus)) +
  facet_wrap(~time) +
  geom_bar(stat = "identity", position = "fill") +
  labs(title = "") +
  xlab("Dietary treatments") +
  ylab("Relative abundance (%)") +
  theme(axis.text.x = element_text(size = 11, vjust = 0.5, hjust=0.5),
        axis.text.y = element_text(size = 11),
        axis.title.y = element_text(size = 11, face = "bold"),
        axis.title.x = element_text(size = 11, face = "bold"))   


SuplementalFig8 <- ggarrange(rel_abun_treatment, rel_abun_treatment.Gen, nrow = 1, labels = c(""))
ggsave("SupplementalFig8.tiff", dpi=300, width = 37, height = 10, units = "cm")


# per sample
dat.agr2 <- aggregate(Abundance ~ SampleID*samples + Family, data = dat.dataframe, FUN = mean)
rel_abun_treatment2 <-ggplot(dat.agr2, aes(x = SampleID, y = Abundance, fill = Family)) +
  geom_bar(stat = "identity", position = "fill") +
  theme(axis.text.x = element_text(size = 10, angle = 90, vjust = 0.5, hjust=1)) +
  labs(title = "")

# Plot each bacteria Abundance by phylum using ggplot2
dat.agr3 <- aggregate(Abundance ~ SampleID*time + Phylum + treatment, data = dat.dataframe, FUN = mean)
#dat.agr3 <- aggregate(Abundance ~ treatment*time + Phylum, data = dat.dataframe, FUN = mean)
abun.fam <- ggplot(dat.agr3, aes(x = Phylum, y = Abundance, fill = treatment)) +
  geom_bar(stat = "identity", position = "dodge") +
  facet_grid(~time) +
  labs(title = "",
       x = "",
       y = "Mean Abundance") +
  theme(axis.text.x = element_text(size = 10, angle = 90, vjust = 0.5, hjust=1)) +
  ylim(0, 0.5) +  
  coord_flip()




SuplementalFig8a <- ggarrange(rel_abun_treatment2, nrow = 1, labels = c(""))
ggsave("SupplementalFig8a.tiff", dpi=300, width = 22, height = 10, units = "cm")

SuplementalFig8b <- ggarrange(abun.fam, nrow = 1, labels = c(""))
ggsave("SupplementalFig8b.tiff", dpi=300, width = 15, height = 10, units = "cm")



# Plot using ggplot2
ggplot(abundance_table, aes(x = Family, y = mean_proportion, fill = treatment)) +
  geom_bar(stat = "identity", position = "dodge") +
  facet_grid(~time) +
  labs(title = "Top 6 Bacterial Families by Abundance",
       x = "Bacterial Family",
       y = "Mean Proportion (%)") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  scale_y_continuous(breaks = seq(0, max(abundance_table$mean_proportion), by = 5))

# Plot using ggplot2
ggplot(taxa_count, aes(x = treatment, y = number_of_taxa, fill = time)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(title = "Number of Taxa by Treatment and Time",
       x = "Treatment",
       y = "Number of Taxa") +
  scale_fill_discrete(name = "Time")


# fecal samples relative amount/sample*time
top20 <- names(sort(taxa_sums(physeq.no.chloro), decreasing=TRUE))[1:30]
ps.top20 <- transform_sample_counts(physeq.no.chloro, function(OTU) OTU/sum(OTU))
ps.top20 <- prune_taxa(top20, ps.top20)
dat.dataframe <- psmelt(ps.top20)
dat.agr <- aggregate(Abundance ~ samples*time + Family, data = dat.dataframe, FUN = mean)
rel_abun_sample <- ggplot(dat.agr, aes(x = samples, y = Abundance, fill = Family)) +
  facet_wrap(~time) +
  geom_bar(stat = "identity", position = "fill") +
  labs(title = "Relative Abundance/Sample") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

##End
