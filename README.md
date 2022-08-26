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
  
