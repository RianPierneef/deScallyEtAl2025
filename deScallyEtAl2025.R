library(phyloseq)
library(MetaLonDA)
library(FSA)
library(ggstatsplot)
library(ggplot2)
library(ggpubr)
library(ggsci)
library(vegan)
library(pairwiseAdonis)
library(readxl)
library(reshape2)
library(phylosmith)

#16S
#Create asvOTUs
ps16s <- readRDS("16S.RDS")
rank_names(ps16s)
OTU <- paste(ps16s@tax_table[,1], ps16s@tax_table[,2], ps16s@tax_table[,3],
             ps16s@tax_table[,4], ps16s@tax_table[,5], ps16s@tax_table[,6], sep = "__")
taxOTU <- cbind(ps16s@tax_table[,1], ps16s@tax_table[,2], ps16s@tax_table[,3],
                ps16s@tax_table[,4], ps16s@tax_table[,5], ps16s@tax_table[,6], OTU=OTU)
tax_table(ps16s) <- taxOTU
ntaxa(ps16s)
ps16sOTU <- tax_glom(ps16s, taxrank="OTU")
ntaxa(ps16sOTU)
ps16sOTU <- subset_taxa(ps16sOTU, !is.na(Phylum) & !Phylum %in% c("", "uncharacterized"))
ntaxa(ps16sOTU)
saveRDS(ps16sOTU, "16S.OTU.RDS")
rm(list=ls())

#Normalize
ps16sOTU <- readRDS("16S.OTU.RDS")
x <- as(otu_table(ps16sOTU), "matrix")
x <- x+1
deseqOTU <- MetaLonDA::normalize(x, method = "median_ratio")
deseqOTU <- floor(deseqOTU)
ps16sOTU.deseq <- phyloseq(otu_table(as.matrix(deseqOTU), taxa_are_rows=FALSE), 
                           ps16sOTU@sam_data, 
                           ps16sOTU@tax_table)
saveRDS(ps16sOTU.deseq, "16S.OTU.DESeq.RDS")
rm(list=ls())

#ITS
#Create asvOTUs
psITS <- readRDS("ITS.RDS")
rank_names(psITS)
OTU <- paste(gsub("k__", "", psITS@tax_table[,1]), sub("p__", "", psITS@tax_table[,2]),
             sub("c__", "", psITS@tax_table[,3]), sub("o__", "", psITS@tax_table[,4]),
             sub("f__", "", psITS@tax_table[,5]), sub("g__", "", psITS@tax_table[,6]), sep = "__")
taxOTU <- cbind(gsub("k__", "", psITS@tax_table[,1]), sub("p__", "", psITS@tax_table[,2]),
                sub("c__", "", psITS@tax_table[,3]), sub("o__", "", psITS@tax_table[,4]),
                sub("f__", "", psITS@tax_table[,5]), sub("g__", "", psITS@tax_table[,6]), OTU=OTU)

tax_table(psITS) <- taxOTU
ntaxa(psITS)
psITSOTU <- tax_glom(psITS, taxrank="OTU")
ntaxa(psITSOTU)
psITSOTU <- subset_taxa(psITSOTU, !is.na(Phylum) & !Phylum %in% c("", "uncharacterized"))
ntaxa(psITSOTU)
saveRDS(psITSOTU, "ITS.OTU.RDS")
rm(list=ls())

#Normalize
psITSOTU <- readRDS("ITS.OTU.RDS")
x <- as(otu_table(psITSOTU), "matrix")
x <- x+1
deseqOTU <- normalize(x, method = "median_ratio")
deseqOTU <- floor(deseqOTU)
psITSOTU.deseq <- phyloseq(otu_table(as.matrix(deseqOTU), taxa_are_rows=FALSE), 
                           psITSOTU@sam_data, 
                           psITSOTU@tax_table)
saveRDS(psITSOTU.deseq, "ITS.OTU.DESeq.RDS")
rm(list=ls())

#16S Alpha Diversity
ps16s.deseq <- readRDS("16S.OTU.DESeq.RDS")
alpha16s <- estimate_richness(ps16s.deseq, measures = c("Observed", "ACE", "Shannon", 
                                                        "Simpson"))
metalpha16s <- cbind(ps16s.deseq@sam_data, alpha16s)

richness.test <- kruskal.test(metalpha16s$Observed~metalpha16s$TreatmentAll)
summary(aov(metalpha16s$Observed~metalpha16s$TreatmentAll))
shannon.test <- kruskal.test(metalpha16s$Shannon~metalpha16s$TreatmentAll)
summary(aov(metalpha16s$Shannon~metalpha16s$TreatmentAll))
simpson.test <- kruskal.test(metalpha16s$Simpson~metalpha16s$TreatmentAll)
summary(aov(metalpha16s$Simpson~metalpha16s$TreatmentAll))
ace1.test <- kruskal.test(metalpha16s$ACE~metalpha16s$TreatmentAll)
summary(aov(metalpha16s$ACE~metalpha16s$TreatmentAll))
pvalues <- matrix(NA, nrow=4, ncol=2)
colnames(pvalues) <- c("Indice", "pvalue")
pvalues[1,1] <- "richness"
pvalues[1,2] <- richness.test$p.value
pvalues[2,1] <- "shannon"
pvalues[2,2] <- shannon.test$p.value
pvalues[3,1] <- "simpson"
pvalues[3,2] <- simpson.test$p.value
pvalues[4,1] <- "ace1"
pvalues[4,2] <- ace1.test$p.value
padjust <- p.adjust(as.numeric(pvalues[,"pvalue"]), method="BH")
pvalues <- cbind(pvalues, padjust)
pvalues

dunnTest(metalpha16s$Observed~metalpha16s$TreatmentAll, method = "bh")
boxplot(metalpha16s$Observed~metalpha16s$TreatmentAl)
dunnTest(metalpha16s$ACE~metalpha16s$TreatmentAll, method = "bh")
boxplot(metalpha16s$ACE~metalpha16s$TreatmentAl)
dunnTest(metalpha16s$Shannon~metalpha16s$TreatmentAll, method = "bh")
boxplot(metalpha16s$Shannon~metalpha16s$TreatmentAl)
dunnTest(metalpha16s$Simpson~metalpha16s$TreatmentAll, method = "bh")
boxplot(metalpha16s$Simpson~metalpha16s$TreatmentAl)

a <- ggbetweenstats(
  data = metalpha16s,
  x = TreatmentAll,
  y = Observed,
  type = "nonparametric",
  p.adjust.method = "BH",
  xlab = "Treatment",
  palette = "default_aaas",
  package = "ggsci",
  ggtheme = theme_classic(),
  k = 4,
  point.args = list(alpha = 1, size = 3),
  violin.args = list(width = 0.5, alpha = 0)
)
b <- ggbetweenstats(
  data = metalpha16s,
  x = TreatmentAll,
  y = ACE,
  type = "nonparametric",
  p.adjust.method = "BH",
  xlab = "Treatment",
  palette = "default_aaas",
  package = "ggsci",
  ggtheme = theme_classic(),
  k = 4,
  point.args = list(alpha = 1, size = 3),
  violin.args = list(width = 0.5, alpha = 0)
)
c <- ggbetweenstats(
  data = metalpha16s,
  x = TreatmentAll,
  y = Shannon,
  type = "nonparametric",
  p.adjust.method = "BH",
  xlab = "Treatment",
  palette = "default_aaas",
  package = "ggsci",
  ggtheme = theme_classic(),
  k = 4,
  point.args = list(alpha = 1, size = 3),
  violin.args = list(width = 0.5, alpha = 0)
)
d <- ggbetweenstats(
  data = metalpha16s,
  x = TreatmentAll,
  y = Simpson,
  type = "nonparametric",
  p.adjust.method = "BH",
  xlab = "Treatment",
  palette = "default_aaas",
  package = "ggsci",
  ggtheme = theme_classic(),
  k = 4,
  point.args = list(alpha = 1, size = 3),
  violin.args = list(width = 0.5, alpha = 0)
)
ggarrange(a, b, c, d, nrow=2, ncol=2, labels=c("a)", "b)", "c)", "d)"))

cor.test(metalpha16s$Observed, metalpha16s$pH, method = "spearman")
cor.test(metalpha16s$ACE, metalpha16s$pH, method = "spearman")
cor.test(metalpha16s$Shannon, metalpha16s$pH, method = "spearman")
cor.test(metalpha16s$Simpson, metalpha16s$pH, method = "spearman")

a <- ggscatterstats(
  data = metalpha16s,
  x = pH,
  y = Observed,
  type = "nonparametric",
  p.adjust.method = "BH",
  xlab = "pH",
  palette = "default_aaas",
  package = "ggsci",
  k = 4,
  ggtheme = theme_classic()
)
b <- ggscatterstats(
  data = metalpha16s,
  x = pH,
  y = ACE,
  type = "nonparametric",
  p.adjust.method = "BH",
  xlab = "pH",
  palette = "default_aaas",
  package = "ggsci",
  k = 4,
  ggtheme = theme_classic()
)
c <- ggscatterstats(
  data = metalpha16s,
  x = pH,
  y = Shannon,
  type = "nonparametric",
  p.adjust.method = "BH",
  xlab = "pH",
  palette = "default_aaas",
  package = "ggsci",
  k = 4,
  ggtheme = theme_classic()
)
d <- ggscatterstats(
  data = metalpha16s,
  x = pH,
  y = Simpson,
  type = "nonparametric",
  p.adjust.method = "BH",
  xlab = "pH",
  palette = "default_aaas",
  package = "ggsci",
  k = 4,
  ggtheme = theme_classic()
)
ggarrange(a, b, c, d, nrow=2, ncol=2, labels=c("a)", "b)", "c)", "d)"))

rm(list=ls())
dev.off()

#ITS Alpha Diversity
psITS.deseq <- readRDS("ITS.OTU.DESeq.RDS")
alphaITS <- estimate_richness(psITS.deseq, measures = c("Observed", "ACE", "Shannon", 
                                                        "Simpson"))
metalphaITS <- cbind(psITS.deseq@sam_data, alphaITS)

richness.test <- kruskal.test(metalphaITS$Observed~metalphaITS$TreatmentAll)
summary(aov(metalphaITS$Observed~metalphaITS$TreatmentAll))
shannon.test <- kruskal.test(metalphaITS$Shannon~metalphaITS$TreatmentAll)
summary(aov(metalphaITS$Shannon~metalphaITS$TreatmentAll))
simpson.test <- kruskal.test(metalphaITS$Simpson~metalphaITS$TreatmentAll)
summary(aov(metalphaITS$Simpson~metalphaITS$TreatmentAll))
ace1.test <- kruskal.test(metalphaITS$ACE~metalphaITS$TreatmentAll)
summary(aov(metalphaITS$ACE~metalphaITS$TreatmentAll))
pvalues <- matrix(NA, nrow=4, ncol=2)
colnames(pvalues) <- c("Indice", "pvalue")
pvalues[1,1] <- "richness"
pvalues[1,2] <- richness.test$p.value
pvalues[2,1] <- "shannon"
pvalues[2,2] <- shannon.test$p.value
pvalues[3,1] <- "simpson"
pvalues[3,2] <- simpson.test$p.value
pvalues[4,1] <- "ace1"
pvalues[4,2] <- ace1.test$p.value
padjust <- p.adjust(as.numeric(pvalues[,"pvalue"]), method="BH")
pvalues <- cbind(pvalues, padjust)
pvalues

dunnTest(metalphaITS$Observed~metalphaITS$TreatmentAll, method = "bh")
boxplot(metalphaITS$Observed~metalphaITS$TreatmentAl)
dunnTest(metalphaITS$ACE~metalphaITS$TreatmentAll, method = "bh")
boxplot(metalphaITS$ACE~metalphaITS$TreatmentAl)
dunnTest(metalphaITS$Shannon~metalphaITS$TreatmentAll, method = "bh")
boxplot(metalphaITS$Shannon~metalphaITS$TreatmentAl)
dunnTest(metalphaITS$Simpson~metalphaITS$TreatmentAll, method = "bh")
boxplot(metalphaITS$Simpson~metalphaITS$TreatmentAl)

a <- ggbetweenstats(
  data = metalphaITS,
  x = TreatmentAll,
  y = Observed,
  type = "nonparametric",
  p.adjust.method = "BH",
  xlab = "Treatment",
  palette = "default_aaas",
  package = "ggsci",
  ggtheme = theme_classic(),
  k = 4,
  point.args = list(alpha = 1, size = 3),
  violin.args = list(width = 0.5, alpha = 0)
)
b <- ggbetweenstats(
  data = metalphaITS,
  x = TreatmentAll,
  y = ACE,
  type = "nonparametric",
  p.adjust.method = "BH",
  xlab = "Treatment",
  palette = "default_aaas",
  package = "ggsci",
  ggtheme = theme_classic(),
  k = 4,
  point.args = list(alpha = 1, size = 3),
  violin.args = list(width = 0.5, alpha = 0)
)
c <- ggbetweenstats(
  data = metalphaITS,
  x = TreatmentAll,
  y = Shannon,
  type = "nonparametric",
  p.adjust.method = "BH",
  xlab = "Treatment",
  palette = "default_aaas",
  package = "ggsci",
  ggtheme = theme_classic(),
  k = 4,
  point.args = list(alpha = 1, size = 3),
  violin.args = list(width = 0.5, alpha = 0)
)
d <- ggbetweenstats(
  data = metalphaITS,
  x = TreatmentAll,
  y = Simpson,
  type = "nonparametric",
  p.adjust.method = "BH",
  xlab = "Treatment",
  palette = "default_aaas",
  package = "ggsci",
  ggtheme = theme_classic(),
  k = 4,
  point.args = list(alpha = 1, size = 3),
  violin.args = list(width = 0.5, alpha = 0)
)
ggarrange(a, b, c, d, nrow=2, ncol=2, labels=c("a)", "b)", "c)", "d)"))

cor.test(metalphaITS$Observed, metalphaITS$pH, method = "spearman")
cor.test(metalphaITS$ACE, metalphaITS$pH, method = "spearman")
cor.test(metalphaITS$Shannon, metalphaITS$pH, method = "spearman")
cor.test(metalphaITS$Simpson, metalphaITS$pH, method = "spearman")

a <- ggscatterstats(
  data = metalphaITS,
  x = pH,
  y = Observed,
  type = "nonparametric",
  p.adjust.method = "BH",
  xlab = "Treatment",
  palette = "default_aaas",
  package = "ggsci",
  k = 4,
  ggtheme = theme_classic()
)
b <- ggscatterstats(
  data = metalphaITS,
  x = pH,
  y = ACE,
  type = "nonparametric",
  p.adjust.method = "BH",
  xlab = "Treatment",
  palette = "default_aaas",
  package = "ggsci",
  k = 4,
  ggtheme = theme_classic()
)
c <- ggscatterstats(
  data = metalphaITS,
  x = pH,
  y = Shannon,
  type = "nonparametric",
  p.adjust.method = "BH",
  xlab = "Treatment",
  palette = "default_aaas",
  package = "ggsci",
  k = 4,
  ggtheme = theme_classic()
)
d <- ggscatterstats(
  data = metalphaITS,
  x = pH,
  y = Simpson,
  type = "nonparametric",
  p.adjust.method = "BH",
  xlab = "Treatment",
  palette = "default_aaas",
  package = "ggsci",
  k = 4,
  ggtheme = theme_classic()
)
ggarrange(a, b, c, d, nrow=2, ncol=2, labels=c("a)", "b)", "c)", "d)"))

rm(list=ls())
dev.off()

#16S Beta Diversity
ps16s.deseq <- readRDS("16S.OTU.deseq.RDS")
ps16s.deseq@sam_data$Day <- as.factor(ps16s.deseq@sam_data$Day)

ps16s.bdiver <- betadiver(ps16s.deseq@otu_table, "z")
ps16s.bdispr <- betadisper(ps16s.bdiver, ps16s.deseq@sam_data$TreatmentAll)
anova(ps16s.bdispr)
adonis2(ps16s.bdiver ~ TreatmentAll, data.frame(ps16s.deseq@sam_data), perm=999)
pairwise.adonis(ps16s.bdiver, ps16s.deseq@sam_data$TreatmentAll, p.adjust.m = "BH")
x <- "MDS"
p1 <- plot_ordination(ps16s.deseq, ordinate(ps16s.deseq, method = x), color = "TreatmentAll", shape = "Day") +
  scale_color_aaas() +
  geom_point(size=3) +
  theme_classic() +
  theme(aspect.ratio=1) +
  guides(color=guide_legend(title="Treatment"))
p1
x <- "PCoA"
p1.2 <- plot_ordination(ps16s.deseq, ordinate(ps16s.deseq, method = x), 
                      color = "TreatmentAll") +
  scale_color_aaas() +
  geom_point(size=3) +
  theme_classic(base_size = 16) +
  theme(aspect.ratio=1) +
  guides(color=guide_legend(title="Treatment"))
p1.2 <- p1.2  + 
  stat_ellipse(type = "norm", linetype = 2) +
  stat_ellipse(type = "t")
p1.2

#ITS Beta Diversity
psITS.deseq <- readRDS("ITS.OTU.DESeq.RDS")
psITS.deseq@sam_data$Day <- as.factor(c("0", "7", "14", "21", "14", "0", "7", "14",
                                        "21", "0", "7", "14", "21", "0", "14", "21",
                                        "14", "21", "0", "7", "14", "21", "0", "7",
                                        "14", "21", "0", "7", "14", "21"))
psITS.deseq@sam_data$Day <- factor(psITS.deseq@sam_data$Day, levels = c("0", "7", "14", "21"))

psITS.bdiver <- betadiver(psITS.deseq@otu_table, "z")
psITS.bdispr <- betadisper(psITS.bdiver, psITS.deseq@sam_data$TreatmentAll)
anova(psITS.bdispr)
adonis2(psITS.bdiver ~ TreatmentAll, data.frame(psITS.deseq@sam_data), perm=999)
pairwise.adonis(psITS.bdiver, psITS.deseq@sam_data$TreatmentAll, p.adjust.m = "BH")
x <- "MDS"
p2 <- plot_ordination(psITS.deseq, ordinate(psITS.deseq, method = x), color = "TreatmentAll", shape = "Day") +
  scale_color_aaas() +
  geom_point(size=3) +
  theme_classic() +
  theme(aspect.ratio=1) +
  guides(color=guide_legend(title="Treatment"))
p2
x <- "PCoA"
p2.2 <- plot_ordination(psITS.deseq, ordinate(psITS.deseq, method = x), color = "TreatmentAll") +
  scale_color_aaas() +
  geom_point(size=3) +
  theme_classic(base_size = 16) +
  theme(aspect.ratio=1) +
  guides(color=guide_legend(title="Treatment"))
p2.2 <- p2.2  + 
  stat_ellipse(type = "norm", linetype = 2) +
  stat_ellipse(type = "t")
p2.2

ggarrange(p1, p2, ncol=2, labels=c("a)", "b)"), vjust = 11)
ggarrange(p1.2, p2.2, ncol=2, labels=c("a)", "b)"), vjust = 11)

rm(list=ls())
dev.off()

#16S Phylum Relative Abundance
ps16s.deseq <- readRDS("16S.OTU.DESeq.RDS")
ps16s.deseq@sam_data$Day <- as.factor(ps16s.deseq@sam_data$Day)

x <- psmelt(ps16s.deseq)
control <- subset(x, x$TreatmentAll=="Control")
low <- subset(x, x$TreatmentAll=="Low")
medium <- subset(x, x$TreatmentAll=="Medium")

control.phyla <- cbind(aggregate(Abundance ~ Phylum, control, sum),
                       rep("Control",nrow(aggregate(Abundance ~ Phylum, control, sum))))
control.phyla <- cbind(control.phyla, control.phyla$Abundance/sum(control.phyla$Abundance)*100)
colnames(control.phyla) <- c("Phylum", "Abundance", "Treatment", "Relative Abundance (%)")
low.phyla <- cbind(aggregate(Abundance ~ Phylum, low, sum),
                   rep("Low",nrow(aggregate(Abundance ~ Phylum, low, sum))))
low.phyla <- cbind(low.phyla, low.phyla$Abundance/sum(low.phyla$Abundance)*100)
colnames(low.phyla) <- c("Phylum", "Abundance", "Treatment", "Relative Abundance (%)")
medium.phyla <- cbind(aggregate(Abundance ~ Phylum, medium, sum),
                      rep("Medium",nrow(aggregate(Abundance ~ Phylum, medium, sum))))
medium.phyla <- cbind(medium.phyla, medium.phyla$Abundance/sum(medium.phyla$Abundance)*100)
colnames(medium.phyla) <- c("Phylum", "Abundance", "Treatment", "Relative Abundance (%)")
control.phyla <- control.phyla[order(-control.phyla$`Relative Abundance (%)`),]
low.phyla <- low.phyla[order(-low.phyla$`Relative Abundance (%)`),]
medium.phyla <- medium.phyla[order(-medium.phyla$`Relative Abundance (%)`),]
res.phyla <- rbind(control.phyla, low.phyla, medium.phyla)

p1 <- ggplot(res.phyla, aes(fill=Phylum, y=`Relative Abundance (%)`, x=Treatment)) + 
  geom_bar(position="stack", stat="identity") +
  theme_classic()+
  scale_fill_igv()
p1

#ITS Phylum Relative Abundance
psITS.deseq <- readRDS("ITS.OTU.DESeq.RDS")
psITS.deseq@sam_data$Day <- as.factor(psITS.deseq@sam_data$Day)

x <- psmelt(psITS.deseq)
control <- subset(x, x$TreatmentAll=="Control")
low <- subset(x, x$TreatmentAll=="Low")
medium <- subset(x, x$TreatmentAll=="Medium")

control.phyla <- cbind(aggregate(Abundance ~ Phylum, control, sum),
                       rep("Control",nrow(aggregate(Abundance ~ Phylum, control, sum))))
control.phyla <- cbind(control.phyla, control.phyla$Abundance/sum(control.phyla$Abundance)*100)
colnames(control.phyla) <- c("Phylum", "Abundance", "Treatment", "Relative Abundance (%)")
low.phyla <- cbind(aggregate(Abundance ~ Phylum, low, sum),
                   rep("Low",nrow(aggregate(Abundance ~ Phylum, low, sum))))
low.phyla <- cbind(low.phyla, low.phyla$Abundance/sum(low.phyla$Abundance)*100)
colnames(low.phyla) <- c("Phylum", "Abundance", "Treatment", "Relative Abundance (%)")
medium.phyla <- cbind(aggregate(Abundance ~ Phylum, medium, sum),
                      rep("Medium",nrow(aggregate(Abundance ~ Phylum, medium, sum))))
medium.phyla <- cbind(medium.phyla, medium.phyla$Abundance/sum(medium.phyla$Abundance)*100)
colnames(medium.phyla) <- c("Phylum", "Abundance", "Treatment", "Relative Abundance (%)")
control.phyla <- control.phyla[order(-control.phyla$`Relative Abundance (%)`),]
low.phyla <- low.phyla[order(-low.phyla$`Relative Abundance (%)`),]
medium.phyla <- medium.phyla[order(-medium.phyla$`Relative Abundance (%)`),]
res.phyla <- rbind(control.phyla, low.phyla, medium.phyla)

p2 <- ggplot(res.phyla, aes(fill=Phylum, y=`Relative Abundance (%)`, x=Treatment)) + 
  geom_bar(position="stack", stat="identity") +
  theme_classic()+
  scale_fill_igv()

ggarrange(p1, p2, ncol=2, labels=c("a)", "b)"))

rm(list=ls())
dev.off()

#SEM
SEM_Datafile_OA <- read_excel("SEM_Datafile_OA.xlsx")
EnzymeTreat <- data.frame(cbind(SampleID=SEM_Datafile_OA$Sample, Treatment=SEM_Datafile_OA$Treatment, 
                                BG=SEM_Datafile_OA$BG, AP=SEM_Datafile_OA$AP, NAG=SEM_Datafile_OA$NAG, 
                                LAP=SEM_Datafile_OA$LAP, BX=SEM_Datafile_OA$BX))
EnzymeTreat.long <- melt(EnzymeTreat, id.vars = c("SampleID", "Treatment"))

#Difference
a <- ggbetweenstats(
  data = SEM_Datafile_OA,
  x = Treatment,
  y = AP,
  type = "nonparametric",
  p.adjust.method = "BH",
  xlab = "Treatment",
  palette = "default_aaas",
  package = "ggsci",
  ggtheme = theme_classic(),
  k = 4,
  point.args = list(alpha = 1, size = 3),
  violin.args = list(width = 0.5, alpha = 0),
  ylab = "Alkaline phosphatase"
) + theme(plot.subtitle = element_text(size = 9))
b <- ggbetweenstats(
  data = SEM_Datafile_OA,
  x = Treatment,
  y = BG,
  type = "nonparametric",
  p.adjust.method = "BH",
  xlab = "Treatment",
  palette = "default_aaas",
  package = "ggsci",
  ggtheme = theme_classic(),
  k = 4,
  point.args = list(alpha = 1, size = 3),
  violin.args = list(width = 0.5, alpha = 0),
  ylab = expression(paste(Beta, "-glucosidase"))
) + theme(plot.subtitle = element_text(size = 9))
c <- ggbetweenstats(
  data = SEM_Datafile_OA,
  x = Treatment,
  y = BX,
  type = "nonparametric",
  p.adjust.method = "BH",
  xlab = "Treatment",
  palette = "default_aaas",
  package = "ggsci",
  ggtheme = theme_classic(),
  k = 4,
  point.args = list(alpha = 1, size = 3),
  violin.args = list(width = 0.5, alpha = 0),
  ylab = expression(paste(Beta, "-xylosidase"))
) + theme(plot.subtitle = element_text(size = 9))
d <- ggbetweenstats(
  data = SEM_Datafile_OA,
  x = Treatment,
  y = LAP,
  type = "nonparametric",
  p.adjust.method = "BH",
  xlab = "Treatment",
  palette = "default_aaas",
  package = "ggsci",
  ggtheme = theme_classic(),
  k = 4,
  point.args = list(alpha = 1, size = 3),
  violin.args = list(width = 0.5, alpha = 0),
  ylab = "Leucine aminopeptidase"
) + theme(plot.subtitle = element_text(size = 9))
e <- ggbetweenstats(
  data = SEM_Datafile_OA,
  x = Treatment,
  y = NAG,
  type = "nonparametric",
  p.adjust.method = "BH",
  xlab = "Treatment",
  palette = "default_aaas",
  package = "ggsci",
  ggtheme = theme_classic(),
  k = 4,
  point.args = list(alpha = 1, size = 3),
  violin.args = list(width = 0.5, alpha = 0),
  ylab = expression(paste("N-acetyl-", Beta, "-D-glucosaminide"))
) + theme(plot.subtitle = element_text(size = 9))

ggarrange(a, b, c,
          d, e, ncol=3, nrow=2, align = "hv")

#Correlation
a <- ggscatterstats(
  data = SEM_Datafile_OA,
  x = pH,
  y = AP,
  type = "nonparametric",
  p.adjust.method = "BH",
  palette = "default_aaas",
  package = "ggsci",
  k = 4,
  ggtheme = theme_classic(),
  ylab = "Alkaline phosphatase",
  ggstatsplot.layer = FALSE,
  ggplot.component = list(theme(text = element_text(size = 9)))
)
b <- ggscatterstats(
  data = SEM_Datafile_OA,
  x = pH,
  y = BG,
  type = "nonparametric",
  p.adjust.method = "BH",
  palette = "default_aaas",
  package = "ggsci",
  k = 4,
  ggtheme = theme_classic(),
  ylab = expression(paste(Beta, "-glucosidase")),
  ggstatsplot.layer = FALSE,
  ggplot.component = list(theme(text = element_text(size = 9)))
)
c <- ggscatterstats(
  data = SEM_Datafile_OA,
  x = pH,
  y = BX,
  type = "nonparametric",
  p.adjust.method = "BH",
  palette = "default_aaas",
  package = "ggsci",
  k = 4,
  ggtheme = theme_classic(),
  ylab = expression(paste(Beta, "-xylosidase")),
  ggstatsplot.layer = FALSE,
  ggplot.component = list(theme(text = element_text(size = 9)))
)
d <- ggscatterstats(
  data = SEM_Datafile_OA,
  x = pH,
  y = LAP,
  type = "nonparametric",
  p.adjust.method = "BH",
  palette = "default_aaas",
  package = "ggsci",
  k = 4,
  ggtheme = theme_classic(),
  ylab = "Leucine aminopeptidase",
  ggstatsplot.layer = FALSE,
  ggplot.component = list(theme(text = element_text(size = 9)))
)
e <- ggscatterstats(
  data = SEM_Datafile_OA,
  x = pH,
  y = NAG,
  type = "nonparametric",
  p.adjust.method = "BH",
  palette = "default_aaas",
  package = "ggsci",
  k = 4,
  ggtheme = theme_classic(),
  ylab = expression(paste("N-acetyl-", Beta, "-D-glucosaminide")),
  ggstatsplot.layer = FALSE,
  ggplot.component = list(theme(text = element_text(size = 9)))
)

ggarrange(a, b, c, d, e, ncol=3, nrow=2)

rm(list=ls())
dev.off()

#ByDay
SEM_Datafile_OA <- read_excel("SEM_Datafile_OA.xlsx")
unique(SEM_Datafile_OA$Day)
#Change SEM_Datafile_OA$Day== to view other days
day <- subset(SEM_Datafile_OA, SEM_Datafile_OA$Day=="21")

a <- ggbetweenstats(
  data = day,
  x = Treatment,
  y = AP,
  type = "nonparametric",
  p.adjust.method = "BH",
  xlab = "Treatment",
  palette = "default_aaas",
  package = "ggsci",
  ggtheme = theme_classic(),
  k = 4,
  point.args = list(alpha = 1, size = 3),
  violin.args = list(width = 0.5, alpha = 0),
  ylab = "Alkaline phosphatase"
) + theme(plot.subtitle = element_text(size = 9))
b <- ggbetweenstats(
  data = day,
  x = Treatment,
  y = BG,
  type = "nonparametric",
  p.adjust.method = "BH",
  xlab = "Treatment",
  palette = "default_aaas",
  package = "ggsci",
  ggtheme = theme_classic(),
  k = 4,
  point.args = list(alpha = 1, size = 3),
  violin.args = list(width = 0.5, alpha = 0),
  ylab = expression(paste(Beta, "-glucosidase"))
) + theme(plot.subtitle = element_text(size = 9))
c <- ggbetweenstats(
  data = day,
  x = Treatment,
  y = BX,
  type = "nonparametric",
  p.adjust.method = "BH",
  xlab = "Treatment",
  palette = "default_aaas",
  package = "ggsci",
  ggtheme = theme_classic(),
  k = 4,
  point.args = list(alpha = 1, size = 3),
  violin.args = list(width = 0.5, alpha = 0),
  ylab = expression(paste(Beta, "-xylosidase"))
) + theme(plot.subtitle = element_text(size = 9))
d <- ggbetweenstats(
  data = day,
  x = Treatment,
  y = LAP,
  type = "nonparametric",
  p.adjust.method = "BH",
  xlab = "Treatment",
  palette = "default_aaas",
  package = "ggsci",
  ggtheme = theme_classic(),
  k = 4,
  point.args = list(alpha = 1, size = 3),
  violin.args = list(width = 0.5, alpha = 0),
  ylab = "Leucine aminopeptidase"
) + theme(plot.subtitle = element_text(size = 9))
e <- ggbetweenstats(
  data = day,
  x = Treatment,
  y = NAG,
  type = "nonparametric",
  p.adjust.method = "BH",
  xlab = "Treatment",
  palette = "default_aaas",
  package = "ggsci",
  ggtheme = theme_classic(),
  k = 4,
  point.args = list(alpha = 1, size = 3),
  violin.args = list(width = 0.5, alpha = 0),
  ylab = expression(paste("N-acetyl-", Beta, "-D-glucosaminide"))
) + theme(plot.subtitle = element_text(size = 9))

ggarrange(a, b, c,
          d, e, ncol=3, nrow=2, align = "hv")

#Correlation
a <- ggscatterstats(
  data = day,
  x = pH,
  y = AP,
  type = "nonparametric",
  p.adjust.method = "BH",
  palette = "default_aaas",
  package = "ggsci",
  k = 4,
  ggtheme = theme_classic(),
  ylab = "Alkaline phosphatase",
  ggstatsplot.layer = FALSE,
  ggplot.component = list(theme(text = element_text(size = 9)))
)
b <- ggscatterstats(
  data = day,
  x = pH,
  y = BG,
  type = "nonparametric",
  p.adjust.method = "BH",
  palette = "default_aaas",
  package = "ggsci",
  k = 4,
  ggtheme = theme_classic(),
  ylab = expression(paste(Beta, "-glucosidase")),
  ggstatsplot.layer = FALSE,
  ggplot.component = list(theme(text = element_text(size = 9)))
)
c <- ggscatterstats(
  data = day,
  x = pH,
  y = BX,
  type = "nonparametric",
  p.adjust.method = "BH",
  palette = "default_aaas",
  package = "ggsci",
  k = 4,
  ggtheme = theme_classic(),
  ylab = expression(paste(Beta, "-xylosidase")),
  ggstatsplot.layer = FALSE,
  ggplot.component = list(theme(text = element_text(size = 9)))
)
d <- ggscatterstats(
  data = day,
  x = pH,
  y = LAP,
  type = "nonparametric",
  p.adjust.method = "BH",
  palette = "default_aaas",
  package = "ggsci",
  k = 4,
  ggtheme = theme_classic(),
  ylab = "Leucine aminopeptidase",
  ggstatsplot.layer = FALSE,
  ggplot.component = list(theme(text = element_text(size = 9)))
)
e <- ggscatterstats(
  data = day,
  x = pH,
  y = NAG,
  type = "nonparametric",
  p.adjust.method = "BH",
  palette = "default_aaas",
  package = "ggsci",
  k = 4,
  ggtheme = theme_classic(),
  ylab = expression(paste("N-acetyl-", Beta, "-D-glucosaminide")),
  ggstatsplot.layer = FALSE,
  ggplot.component = list(theme(text = element_text(size = 9)))
)

ggarrange(a, b, c, d, e, ncol=3, nrow=2)

rm(list=ls())
dev.off()

#Network Robustness

library(igraph)
library(brainGraph)
library(DescTools)

#16S
ps16s.deseq <- readRDS("16S.OTU.DESeq.RDS")
ps16s.deseq@sam_data$Day <- as.factor(ps16s.deseq@sam_data$Day)

#Control 16s
control16s <- subset_samples(ps16s.deseq, TreatmentAll=="Control")
controlCoOc16s <- co_occurrence(control16s, p = 0.05, rho = 0.85)
controlCoOc16sPos <- subset(controlCoOc16s, controlCoOc16s$rho>0)

controlCoOc16sPos.links <- data.frame(
  source=controlCoOc16sPos$X,
  target=controlCoOc16sPos$Y,
  weight=controlCoOc16sPos$rho
)
controlCoOc16sPos.igraph <- graph_from_data_frame(controlCoOc16sPos.links, directed=F)
is_directed(controlCoOc16sPos.igraph)

#Medium 16s
medium16s <- subset_samples(ps16s.deseq, TreatmentAll=="Medium")
mediumCoOc16s <- co_occurrence(medium16s, p = 0.05, rho = 0.85)
mediumCoOc16sPos <- subset(mediumCoOc16s, mediumCoOc16s$rho>0)

mediumCoOc16sPos.links <- data.frame(
  source=mediumCoOc16sPos$X,
  target=mediumCoOc16sPos$Y,
  weight=mediumCoOc16sPos$rho
)
mediumCoOc16sPos.igraph <- graph_from_data_frame(d=mediumCoOc16sPos.links, directed=F)
is_directed(mediumCoOc16sPos.igraph)

#Low 16s
low16s <- subset_samples(ps16s.deseq, TreatmentAll=="Low")
lowCoOc16s <- co_occurrence(low16s, p = 0.05, rho = 0.85)
lowCoOc16sPos <- subset(lowCoOc16s, lowCoOc16s$rho>0)

lowCoOc16sPos.links <- data.frame(
  source=lowCoOc16sPos$X,
  target=lowCoOc16sPos$Y,
  weight=lowCoOc16sPos$rho
)
lowCoOc16sPos.igraph <- graph_from_data_frame(d=lowCoOc16sPos.links, directed=F)
is_directed(lowCoOc16sPos.igraph)

S16aucStats <- data.frame(matrix(ncol = 2, nrow = 3))
x <- c("Treatment", "AUC")
colnames(S16aucStats) <- x

S16aucStats[1,1] <- "Low"
low.rob <- robustness(lowCoOc16sPos.igraph)
AUC(x = low.rob$removed.pct, y = low.rob$comp.pct)
S16aucStats[1,2] <- AUC(x = low.rob$removed.pct, y = low.rob$comp.pct)
low.rob$Treatment <- c(rep("Low", nrow(low.rob)))

S16aucStats[2,1] <- "Medium"
med.rob <- robustness(mediumCoOc16sPos.igraph)
AUC(x = med.rob$removed.pct, y = med.rob$comp.pct)
S16aucStats[2,2] <- AUC(x = med.rob$removed.pct, y = med.rob$comp.pct)
med.rob$Treatment <- c(rep("Medium", nrow(med.rob)))

S16aucStats[3,1] <- "Control"
con.rob <- robustness(controlCoOc16sPos.igraph)
AUC(x = con.rob$removed.pct, y = con.rob$comp.pct)
S16aucStats[3,2] <- AUC(x = con.rob$removed.pct, y = con.rob$comp.pct)
con.rob$Treatment <- c(rep("Control", nrow(con.rob)))

rob <- rbind(low.rob, med.rob, con.rob)

S16aucGraph <- ggplot(data=rob,
                      aes(x=removed.pct, y=comp.pct, colour = Treatment)) +
  geom_line(linewidth = 1.5) +
  geom_abline(slope=-1, intercept=1, col='gray', lty=2) +
  theme_classic() +
  scale_color_aaas() +
  theme_classic(base_size = 14) +
  xlab("Vertices Removed") +
  ylab("Maximal Component Size")

S16aucGraph

S16aucStatsGraph <- ggplot(S16aucStats, aes(y=AUC, x=Treatment, fill=Treatment)) + 
  geom_bar(stat='identity', position='dodge') +
  theme_classic(base_size = 14) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  scale_fill_aaas()

S16aucStatsGraph

#ITS
psITS.deseq <- readRDS("ITS.OTU.DESeq.RDS")
psITS.deseq@sam_data$Day <- as.factor(c("0", "7", "14", "21", "14", "0", "7", "14",
                                        "21", "0", "7", "14", "21", "0", "14", "21",
                                        "14", "21", "0", "7", "14", "21", "0", "7",
                                        "14", "21", "0", "7", "14", "21"))

#Control ITS
controlITS <- subset_samples(psITS.deseq, TreatmentAll=="Control")
controlCoOcITS <- co_occurrence(controlITS, p = 0.05, rho = 0.85)
controlCoOcITSPos <- subset(controlCoOcITS, controlCoOcITS$rho>0)

controlCoOcITSPos.links <- data.frame(
  source=controlCoOcITSPos$X,
  target=controlCoOcITSPos$Y,
  weight=controlCoOcITSPos$rho
)
controlCoOcITSPos.igraph <- graph_from_data_frame(d=controlCoOcITSPos.links, directed=F)
is_directed(controlCoOcITSPos.igraph)

#Medium ITS
mediumITS <- subset_samples(psITS.deseq, TreatmentAll=="Medium")
mediumCoOcITS <- co_occurrence(mediumITS, p = 0.05, rho = 0.85)
mediumCoOcITSPos <- subset(mediumCoOcITS, mediumCoOcITS$rho>0)

mediumCoOcITSPos.links <- data.frame(
  source=mediumCoOcITSPos$X,
  target=mediumCoOcITSPos$Y,
  weight=mediumCoOcITSPos$rho
)
mediumCoOcITSPos.igraph <- graph_from_data_frame(d=mediumCoOcITSPos.links, directed=F)
is_directed(mediumCoOcITSPos.igraph)

#Low ITS
lowITS <- subset_samples(psITS.deseq, TreatmentAll=="Low")
lowCoOcITS <- co_occurrence(lowITS, p = 0.05, rho = 0.85)
lowCoOcITSPos <- subset(lowCoOcITS, lowCoOcITS$rho>0)

lowCoOcITSPos.links <- data.frame(
  source=lowCoOcITSPos$X,
  target=lowCoOcITSPos$Y,
  weight=lowCoOcITSPos$rho
)
lowCoOcITSPos.igraph <- graph_from_data_frame(d=lowCoOcITSPos.links, directed=F)
is_directed(lowCoOcITSPos.igraph)

ITSaucStats <- data.frame(matrix(ncol = 2, nrow = 3))
x <- c("Treatment", "AUC")
colnames(ITSaucStats) <- x

ITSaucStats[1,1] <- "Low"
low.rob <- robustness(lowCoOcITSPos.igraph)
AUC(x = low.rob$removed.pct, y = low.rob$comp.pct)
ITSaucStats[1,2] <- AUC(x = low.rob$removed.pct, y = low.rob$comp.pct)
low.rob$Treatment <- c(rep("Low", nrow(low.rob)))

ITSaucStats[2,1] <- "Medium"
med.rob <- robustness(mediumCoOcITSPos.igraph)
AUC(x = med.rob$removed.pct, y = med.rob$comp.pct)
ITSaucStats[2,2] <- AUC(x = med.rob$removed.pct, y = med.rob$comp.pct)
med.rob$Treatment <- c(rep("Medium", nrow(med.rob)))

ITSaucStats[3,1] <- "Control"
con.rob <- robustness(controlCoOcITSPos.igraph)
AUC(x = con.rob$removed.pct, y = con.rob$comp.pct)
ITSaucStats[3,2] <- AUC(x = con.rob$removed.pct, y = con.rob$comp.pct)
con.rob$Treatment <- c(rep("Control", nrow(con.rob)))

rob <- rbind(low.rob, med.rob, con.rob)

ITSaucGraph <- ggplot(data=rob,
                      aes(x=removed.pct, y=comp.pct, colour = Treatment)) +
  geom_line(linewidth = 1.5) +
  geom_abline(slope=-1, intercept=1, col='gray', lty=2) +
  theme_classic() +
  scale_color_aaas() +
  theme_classic(base_size = 14) +
  xlab("Vertices Removed") +
  ylab("Maximal Component Size")

ITSaucGraph

ITSaucStatsGraph <- ggplot(ITSaucStats, aes(y=AUC, x=Treatment, fill=Treatment)) + 
  geom_bar(stat='identity', position='dodge') +
  theme_classic(base_size = 14) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  scale_fill_aaas()

ITSaucStatsGraph

ggarrange(S16aucGraph, S16aucStatsGraph,
          ITSaucGraph, ITSaucStatsGraph,
          nrow = 2, ncol = 2, align = "hv", labels = c("a)", "b)", "c)", "d)"))

rm(list=ls())
dev.off()

sessionInfo()
writeLines(capture.output(sessionInfo()), "sessionInfo.txt")
