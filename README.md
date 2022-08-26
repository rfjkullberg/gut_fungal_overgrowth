# Gut fungal overgrowth, anaerobic bacteria and parasitome

This is the code used for the main analyses in "Gut fungal overgrowth is associated with reduced prevalence of Blastocystis species in a large human cohort" (submitted). For questions: Jason Biemond at j.j.biemond@amsterdamumc.nl or Bob Kullberg at r.f.j.kullberg@amsterdamumc.nl


## Step 1 - Load libraries
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

## Step 2 - Load data
```
df <- read_excel("~/Documents/PhD/Fungal overgrowth/Data/Github_metadata.xlsx", 
                 col_types = c("text",  "text", "numeric", "numeric", "numeric", 
                               "text", "text", "text", "text", "text", "text", "numeric", 
                               "text", "text", "numeric", "text", "numeric")) %>%
  mutate(group = fct_relevel(group, "absent", "low", "high"))
```

Microbiota sequence data (phyloseq file) are integrated with the taxonomy and a phylogenetic tree using the phyloseq package (details described in the manuscript). Contaminants were identified using the package decontam and removed from the dataset. 
```
ps <- readRDS("~/Documents/PhD/Fungal overgrowth/Data/ps.2021_22_Fungal_overgrowth.2021-08-27.RDS")

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

Remove negative controls from the dataset for subsequent analyses
```
df <- df %>%
  filter(subject_id != "NEG_CON_1" | is.na(subject_id)) %>%
  filter(subject_id != "NEG_CON_2" | is.na(subject_id))
```

## Step 3 - Fungal growth
Gut fungal growth ranges from 0 to 191318760 CFUs per gram stool. We divided our cohort into three groups: no fungal growth (n=987), low-to-normal fungal growth (n=1745) and fungal overgrowth (n=115). 
```
summary(df$yeast_cfu)
table(df$group) 
```
```
ggplot(df, aes(x=reorder(sample_id,yeast_cfu), y=yeast_cfu, color = group)) +
  geom_bar(stat = "identity")+
  ylab("Yeast CFUs per gram stool")+ 
  xlab("Participants (n=2847)")+ 
  theme_cowplot(18)+
  scale_y_continuous(trans=log_epsilon_trans(epsilon=100),expand =c(0,0))+
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        legend.title = element_blank(),
        legend.position = "none")
```

## Step 4 - Fungal growth and parasitome
Colonization with the anaerobic gut parasites Blastocystis species, Dientamoeba fragilis, and Giardia lamblia 
```
ggplot(df, aes(x=group, fill = blasto)) +
  geom_bar(position = "fill")+
  ylab("% of participants")+ 
  xlab("")+ 
  ggtitle("Blastocystis spp") +
  theme_cowplot(18)+
  scale_y_continuous(expand =c(0,0))+
  scale_fill_manual(values=c("#9f514d", "#49a258")) +
  theme(plot.title = element_text(size=14, face="plain"), 
        legend.position = "none") +
  coord_flip()
ggplot(df, aes(x=group, fill = dienta)) +
  geom_bar(position = "fill")+
  ylab("% of participants")+   
  xlab("")+ 
  ggtitle("Dientamoeba fragilis") +
  theme_cowplot(18)+
  scale_y_continuous(expand =c(0,0))+
  scale_fill_manual(values=c('#9f514d', "#49a258")) +
  theme(plot.title = element_text(size=14, face="plain"), 
        legend.position = "none") +
  coord_flip()
ggplot(df, aes(x=group, fill = giardia)) +
  geom_bar(position = "fill")+
  ylab("% of participants")+   
  xlab("")+ 
  ggtitle("Giardia lamblia") +
  theme_cowplot(18)+
  scale_y_continuous(expand =c(0,0))+
  scale_fill_manual(values=c('#9f514d', "#49a258")) +
  theme(plot.title = element_text(size=14, face="plain"), 
        legend.position = "none") +
  coord_flip()
```

Colonization with Blastocystis spp was less frequent in participants with fungal overgrowth compared to those with no or low-to-normal fungal growth. 
```
highabsent <- df %>%
  filter(group != "low") 
chisq.test(highabsent$group, highabsent$blasto, correct = T)
#chisq.test(highabsent$group, highabsent$dienta, correct = T) 
#chisq.test(highabsent$group, highabsent$giardia, correct = T) 

highlow <- df %>%
  filter(group != "absent") 
chisq.test(highlow$group, highlow$blasto, correct = T) 
#chisq.test(highlow$group, highlow$dienta, correct = T) 
#chisq.test(highlow$group, highlow$giardia, correct = T) 

#lowabsent <- df %>%
  #filter(group != "high") 
#chisq.test(lowabsent$group, lowabsent$blasto, correct = T)
#chisq.test(lowabsent$group, lowabsent$dienta, correct = T)
#chisq.test(lowabsent$group, lowabsent$giardia, correct = T)

rm(highabsent, highlow)
```

Participants colonized by Blastocystis spp. had a higher faecal fungal burden compared to non-colonized participants. 
```
comparisons <- list(c("pos", "neg"))
ggplot(df, aes(x = blasto, y = yeast_cfu, fill = blasto))+
  geom_boxplot(alpha = 0.5, outlier.shape = NA) +
  geom_jitter(color = "black", pch = 21, alpha =.75, size = 2, width=0.15) +
  scale_y_continuous(trans=log_epsilon_trans(epsilon=100),expand = expansion(mult = c(0, .05)))+
  stat_compare_means(method = "wilcox.test", comparisons = comparisons, size=6, label = "p.value") +
  scale_color_manual(values=c('#9f514d', "#49a258")) +
  scale_fill_manual(values=c('#9f514d', "#49a258")) +
  theme_bw(base_size=14) +
  ylab("Fungal CFUs per gram stool") +
  theme(legend.position = "none")
```

## Step 5 - Fungal growth and bacteriome
Faecal fungal growth was inversely correlated with gut bacterial diversity

```
alpha <- estimate_richness(ps) 
alpha$sample <- row.names(alpha)
alpha$sample <- gsub("X","",as.character(alpha$sample))

alpha <- alpha %>%
  left_join(s)
alpha <- alpha %>%
  arrange(match)
alpha$group <- factor(alpha$group, levels = c("absent", "low", "high"))
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

Accordingly, bacterial diversity was lower in patients with fungal overgrowth compared to patients without faecal fungal growth
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
# compare high vs absent
ps.highabsent <- ps
highabsent <- s %>%
  filter(group != "low")
sample_data(ps.highabsent) <- set.samp(highabsent) 

ps.deseq <- tax_glom(ps.highabsent, "Genus")
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
  select(Family,log2FoldChange) %>%
  group_by(Family) %>%
  summarise_at(c("log2FoldChange"), sum, na.rm=T)%>%
  mutate(group=ifelse(log2FoldChange<0, "Deceased or Intubated", "Extubated")) 

ggplot(deseq, aes(x=reorder(Family,log2FoldChange), y=log2FoldChange, fill=group), 
                     stat="identity", color= "black")+
  geom_bar(stat = "identity") + 
  coord_flip() +
  ylab("log 2-Fold Change")+
  xlab("") +
  theme_cowplot(11) +
  scale_fill_manual(values = c("#752936","#44aa99"))+
  theme(legend.position = "none") 

rm(highabsent, ps.highabsent, deseq, dsq, geoMeans, res, ps.deseq, gm_mean)
```
```
# compare high vs low
ps.highlow <- ps
highlow <- s %>%
  filter(group != "absent")
sample_data(ps.highlow) <- set.samp(highlow) 

ps.deseq <- tax_glom(ps.highlow, "Genus")
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
  select(Family,log2FoldChange) %>%
  group_by(Family) %>%
  summarise_at(c("log2FoldChange"), sum, na.rm=T)%>%
  mutate(group=ifelse(log2FoldChange<0, "Deceased or Intubated", "Extubated")) 

ggplot(deseq, aes(x=reorder(Family,log2FoldChange), y=log2FoldChange, fill=group), 
       stat="identity", color= "black")+
  geom_bar(stat = "identity") + 
  coord_flip() +
  ylab("log 2-Fold Change")+
  xlab("") +
  theme_cowplot(11) +
  scale_fill_manual(values = c("#752936","#44aa99"))+
  theme(legend.position = "none") 
rm(highlow, ps.highlow, deseq, dsq, geoMeans, res, ps.deseq)
```



  
