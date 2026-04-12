# Gut fungal overgrowth, anaerobic bacteria and parasitome

This is the data used for the main analyses in "Assessment of Relationships Between Faecal Yeast Overgrowth and Yeast Species, Protozoan Parasites, and Bacterial Population Dynamics: a Retrospective, Single-Centre, Cohort Study" (submitted). For questions: Jason Biemond at j.j.biemond@amsterdamumc.nl or Bob Kullberg at r.f.j.kullberg@amsterdamumc.nl


The following files are available:
- Fungal quantification and presence of parasites per sample
- Fungal species identification
- Clinical metadata (detailed data are available upon request)
- Raw sequencing data are availble at the European Nucleotide Archive (accession number PRJEB82678)


## Code for bacterial microbiota analyses
```
library(tidyverse)
library(yingtools2)
library(phyloseq)
library(vegan)
library(microbiome)
library(RColorBrewer) 
library(ggpubr) 
library(data.table)
library(cowplot)
library(DESeq2)
library(splitstackshape)
library(DirichletMultinomial)
library(readxl)
library(scales)
library(decontam)
library(ggpmisc)
library(tableone)
```
```
df <- read_excel("~/Documents/Fungal overgrowth/Data/metadata.xlsx")
```

Microbiota sequence data (phyloseq file) are integrated with the taxonomy and a phylogenetic tree using the phyloseq package (details described in the manuscript). Contaminants were identified using the package decontam and removed from the dataset. 
```
ps <- readRDS("~/Documents/Fungal overgrowth/Data/ps.2021_22_Fungal_overgrowth.2021-08-27.RDS")

## Delete contaminants
# correct Nucl_Acid_AMP_Conc
ps@sam_data$Nucl_Acid_AMP_Conc <- as.numeric(gsub(",",".", as.character(ps@sam_data$Nucl_Acid_AMP_Conc)))
ps@sam_data$Nucl_Acid_AMP_Conc_RS <- ps@sam_data$Nucl_Acid_AMP_Conc
ps@sam_data$Nucl_Acid_AMP_Conc_RS[is.na(ps@sam_data$Nucl_Acid_AMP_Conc_RS)] <- min(ps@sam_data$Nucl_Acid_AMP_Conc, na.rm = T)
ps@sam_data$Nucl_Acid_AMP_Conc_RS <- ps@sam_data$Nucl_Acid_AMP_Conc_RS - min(ps@sam_data$Nucl_Acid_AMP_Conc_RS) + 0.01
#primer id
ps@sam_data$i5 <- substr(ps@sam_data$Index_Name,1,5)
ps@sam_data$i7 <- substr(ps@sam_data$Index_Name,6,10)
dco <- isContaminant(ps, conc = "Nucl_Acid_AMP_Conc_RS", method = "frequency")

noncontaminants <- get.tax(ps) %>% #list of non-contaminans
  filter(dco$contaminant == F)
contaminants <- get.tax(ps) %>% #list of contaminans
  filter(dco$contaminant == T)

tax_table(ps) <- set.tax(noncontaminants) # delete contaminants

# Add metadata to phyloseq file
s <- df %>%
  filter(bacteriome == "y") %>%
  filter(subject_id != "NEG_CON_1") %>%
  filter(subject_id != "NEG_CON_2") %>%
  filter(match != "1")
sample_data(ps) <- set.samp(s) 
rm(contaminants, noncontaminants, dco)
```
```
# Remove negative controls from the dataset for subsequent analyses
df <- df %>%
  filter(subject_id != "NEG_CON_1" | is.na(subject_id)) %>%
  filter(subject_id != "NEG_CON_2" | is.na(subject_id))
```

```
# Bacterial diversity
alpha <- estimate_richness(ps) 
alpha$sample <- row.names(alpha)
alpha$sample <- gsub("X","",as.character(alpha$sample))

alpha <- alpha %>%
  left_join(s)
alpha <- alpha %>%
  arrange(match)
lev <- levels(alpha$group) # get the variables
L.pairs <- combn(seq_along(lev), 2, simplify = FALSE, FUN = function(i)lev[i])

alpha %>%
  ggplot(aes(x = yeast_cfu, y = Shannon))+
  geom_point(size=2.5, colour = "#000000") +
  geom_smooth(method="lm", se=F, fullrange=FALSE, span = 0.99,  level=0.95, colour = "black") +
  scale_x_continuous(trans=log_epsilon_trans(epsilon=10000))+
  theme_bw()+
  xlab("Fecal yeast CFU per gram stool")+
  ylab("Shannon Bacterial Diversity")+
  theme(legend.position = "none")
cor.test((alpha$log_yeast_cfu), alpha$Shannon, method="spearman") 
```
```
alpha %>%
  ggplot(aes(x = group, y = Shannon, fill = group))+
  geom_boxplot(alpha = 0.5, outlier.shape = NA, show.legend = FALSE) +
  geom_point(color = "black", pch = 21, alpha =.75, size = 2, show.legend = FALSE)+
  geom_line(aes(group=match), alpha = 0.3) +
  theme_cowplot(11)+
  ggtitle("") +
  ylab("Shannon Bacterial Diversity")+
  xlab("")+
  stat_compare_means(method = "wilcox.test", comparisons = L.pairs, paired = TRUE,exact = FALSE,
                     size=4, label = "p.value")
rm(alpha, L.pairs, lev)
```

We used a DESeq2 model to identify bacterial taxa that were correlated with faecal fungal growth
```
# compare overgrowth vs absent
ps.overgrowthabsent <- ps
overgrowthabsent <- s %>%
  filter(group == "overgrowth" | group == "absent")
sample_data(ps.overgrowthabsent) <- set.samp(overgrowthabsent) 

ps.deseq <- tax_glom(ps.overgrowthabsent, "Genus")
ps.deseq <- core(ps.deseq, detection = 1, prevalence = 10/100, include.lowest = T) # prevalence of 10%

gm_mean <- function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))}

dsq <- phyloseq_to_deseq2(ps.deseq,~group ) 
geoMeans <- apply(counts(dsq), 1, gm_mean)
dsq <- estimateSizeFactors(dsq, geoMeans = geoMeans) 
dsq <- DESeq(dsq,  fitType="local")    
res <- results(dsq, cooksCutoff = FALSE, pAdjustMethod = "BH" )
deseq <- res[which(res$padj < 0.05), ]  #adjusted p-value <0.05
deseq <- cbind(as(deseq, "data.frame"), as(tax_table(ps.deseq)[rownames(deseq), ], "matrix"))
deseq <- deseq %>% 
  select(Genus,log2FoldChange) %>%
  group_by(Genus) %>%
  summarise_at(c("log2FoldChange"), sum, na.rm=T)

ggplot(deseq, aes(x=reorder(Genus,log2FoldChange), y=log2FoldChange, fill=group), 
                     stat="identity", color= "black")+
  geom_bar(stat = "identity") + 
  coord_flip() +
  ylab("log 2-Fold Change")+
  xlab("") +
  theme_cowplot(11) +
  scale_fill_manual(values = c("#752936","#44aa99"))+
  theme(legend.position = "none") 
```
Similar code was used to compare low with overgrowth. 
