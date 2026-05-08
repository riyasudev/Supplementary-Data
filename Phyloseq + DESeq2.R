######Install packages###########

if(!requireNamespace("BiocManager")){
  install.packages("BiocManager")
}
BiocManager::install("phyloseq", version = "3.20") 

library("phyloseq")

BiocManager::install("DESeq2")

library(DESeq2)

install.packages("remotes")
remotes::install_github("vmikk/metagMisc")

library("ggplot2")
library("scales")
library("grid")
library("magrittr")
library("ggpubr")
library("vegan")
library("dplyr")
library("ggsci")
library("viridis")
library("data.table")
library("metagMisc")
library("iNEXT")

#STEP 2. Set WORKING DIRECTORY on Rstudio main menu -> Session -> set working directory##############

#############Load data tables#######################
ASVTable<-read.csv("ASV_Count_Table.csv", header = TRUE, stringsAsFactors=FALSE,row.names = 1)

MetadataTable<-read.csv("Metadata_Table.csv", header = TRUE, stringsAsFactors=FALSE, row.names = 1)

TaxonomyTable<-read.csv("Taxonomy_Table_backfilled.csv", header = TRUE, stringsAsFactors=FALSE, row.names = 1)

sam1 <- sample_data(MetadataTable) 
tax1 <- tax_table(TaxonomyTable)
asv1 <- otu_table(ASVTable, taxa_are_rows=TRUE)

colnames(tax1) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus")
rownames(tax1)<- rownames(asv1)


Microbiome<- phyloseq(asv1,tax1,sam1)

#check data
most_abundant_taxa <- sort(taxa_sums(Microbiome), TRUE)[1:25]
ex2 <- prune_taxa(names(most_abundant_taxa), Microbiome)
topFamilies <- tax_table(ex2)[, "Family"]
as(topFamilies, "vector")

#######################################################################################
#Distribution of read counts
sample_sum_df <- data.frame(sum = sample_sums(Microbiome))

ggplot(sample_sum_df, aes(x = sum)) + 
  geom_histogram(color = "black", fill = "indianred", binwidth = 2500) +
  ggtitle("Distribution of sample sequencing depth") + 
  xlab("Read counts") +
  theme(axis.title.y = element_blank())

# sequencing stats
smin <- min(sample_sums(Microbiome))
smean <- mean(sample_sums(Microbiome))
ssd<- sd(sample_sums(Microbiome))
smax <- max(sample_sums(Microbiome))

smin
smean
ssd
smax

#coverage
coverage_all<-phyloseq_coverage(Microbiome,correct_singletons = FALSE, add_attr = T)
coverage_all

write.csv(coverage_all, file = "phyloseq_coverage_all.csv")

####################prune taxa#################

Microbiome_prune = prune_taxa(taxa_sums(Microbiome) > 25, Microbiome)

#################################### GROUPING VARIABLE #########################################

sample_data(Microbiome_prune)$group <- factor(
  sample_data(Microbiome_prune)$Disease,
  levels = c("Control","ASD")
)

#################################### ALPHA DIVERSITY #########################################

group_color <- c("#009E73","#D55E00")

alpha_all<-plot_richness(
  Microbiome_prune,
  x="group",
  color="group",
  measures=c("InvSimpson", "Shannon","Observed")
) +
  scale_color_manual(values=group_color) +
  geom_boxplot(lwd=1,width=0.5) +
  geom_point()

alpha_all

#################################### STATISTICS #########################################

alpha_stats <- alpha_all +
  stat_compare_means(method = "wilcox.test", label = "p.format") +
  theme_bw()

alpha_stats

#################################### SAVE TABLE #########################################

richness.rare<- cbind(
  estimate_richness(
    Microbiome_prune,
    measures = c('Observed','Shannon',"Simpson","InvSimpson")
  ),
  sample_data(Microbiome_prune)$group
)

write.csv(richness.rare, file = "alphaDiv_ASD_vs_Control.csv")

library(dplyr)

alpha_df <- estimate_richness(
  Microbiome_prune,
  measures = c("Observed","Shannon","Simpson","InvSimpson")
)

alpha_df$group <- sample_data(Microbiome_prune)$group


################# DESCRIPTIVE STATISTICS ########################

alpha_descriptives <- alpha_df %>%
  group_by(group) %>%
  summarise(
    n = n(),
    Observed_median = median(Observed),
    Observed_IQR = IQR(Observed),
    Shannon_median = median(Shannon),
    Shannon_IQR = IQR(Shannon),
    InvSimpson_median = median(InvSimpson),
    InvSimpson_IQR = IQR(InvSimpson)
  )
alpha_descriptives


################### WILCOXON TESTS (W, p, n) #########################

wilcox_observed <- wilcox.test(Observed ~ group, data = alpha_df)
wilcox_shannon <- wilcox.test(Shannon ~ group, data = alpha_df)
wilcox_invsimpson <- wilcox.test(InvSimpson ~ group, data = alpha_df)

alpha_stats_table <- data.frame(
  Metric = c("Observed","Shannon","InvSimpson"),
  W = c(wilcox_observed$statistic,
        wilcox_shannon$statistic,
        wilcox_invsimpson$statistic),
  p_value = c(wilcox_observed$p.value,
              wilcox_shannon$p.value,
              wilcox_invsimpson$p.value),
  n_group1 = sum(alpha_df$group == "Control"),
  n_group2 = sum(alpha_df$group == "ASD")
)

alpha_stats_table

write.csv(alpha_descriptives, "alpha_descriptives_ASD_vs_Control.csv")
write.csv(alpha_stats_table, "alpha_stats_ASD_vs_Control.csv")

#################################### RELATIVE ABUNDANCE #########################################

Microbiome_Rel  = transform_sample_counts(Microbiome_prune, function(x) x / sum(x) )

write.csv(as.data.frame(otu_table(Microbiome_Rel)) , file = "Microbiome_Rel_abundance.csv")

#################### COMPOSITIONAL ANALYSIS #################

##FAMILIES

library(RColorBrewer)
Rel_family <- tax_glom(Microbiome_Rel, taxrank = 'Family', NArm=FALSE) 
FamilyRel <- psmelt(Rel_family) 
FamilyRel$Family <- as.character(FamilyRel$Family) 
FamilyRel$Family[FamilyRel$Abundance < 0.05] <- "Other"

high_contrast_palette <- c(
  "#E63946", "#2196F3", "#FF9800", "#4CAF50", "#9C27B0",
  "#00BCD4", "#FF5722", "#8BC34A", "#673AB7", "#FFC107",
  "#F06292", "#26C6DA", "#D4E157", "#AB47BC", "#42A5F5",
  "#FF7043", "#66BB6A", "#EC407A", "#29B6F6", "#FFCA28",
  "#EF5350", "#26A69A", "#7E57C2", "#9CCC65", "#FFA726"
)

n_families <- length(unique(FamilyRel$Family))
family_colors <- setNames(high_contrast_palette[1:n_families], unique(FamilyRel$Family))

ggplot(FamilyRel, aes(x = Sample, y = Abundance/3, fill = Family)) +
  geom_bar(stat="identity") +
  ylab("Relative abundance") +
  xlab("") +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust=1, size=7)
  ) +
  facet_wrap(~ group, scales="free_x", nrow=1) +
  scale_fill_manual(values=family_colors)

##GENUS

library(RColorBrewer)
Rel_genus <- tax_glom(Microbiome_Rel, taxrank = 'Genus', NArm=FALSE) 
GenusRel <- psmelt(Rel_genus) 
GenusRel$Genus <- as.character(GenusRel$Genus) 
GenusRel$Genus[GenusRel$Abundance < 0.05] <- "Other"

high_contrast_palette <- c(
  "#E63946", "#2196F3", "#FF9800", "#4CAF50", "#9C27B0",
  "#00BCD4", "#FF5722", "#8BC34A", "#673AB7", "#FFC107",
  "#F06292", "#26C6DA", "#D4E157", "#AB47BC", "#42A5F5",
  "#FF7043", "#66BB6A", "#EC407A", "#29B6F6", "#FFCA28",
  "#EF5350", "#26A69A", "#7E57C2", "#9CCC65", "#FFA726",
  "#B71C1C", "#1A237E", "#E65100", "#1B5E20", "#4A148C",
  "#006064", "#BF360C", "#33691E", "#311B92", "#F57F17",
  "#880E4F", "#0D47A1", "#E84393", "#558B2F", "#6A1B9A",
  "#00838F", "#D84315", "#558B2F", "#4527A0", "#F9A825",
  "#C62828", "#283593", "#FF6F00", "#2E7D32", "#6A1B9A"
)

n_genera <- length(unique(GenusRel$Genus))
genus_colors <- setNames(high_contrast_palette[1:n_genera], unique(GenusRel$Genus))

ggplot(GenusRel, aes(x = Sample, y = Abundance/3, fill = Genus)) +
  geom_bar(stat="identity") +
  ylab("Relative abundance") +
  xlab("") +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust=1, size=7)
  ) +
  facet_wrap(~ group, scales="free_x", nrow=1) +
  scale_fill_manual(values=genus_colors)

#################################### BETA DIVERSITY #########################################
# Bray Curtis distance
bc_dist = phyloseq::distance(Microbiome_Rel, method="bray", weighted=FALSE)

################ PERMANOVA ################
adonis2(bc_dist ~ group, data = data.frame(sample_data(Microbiome_Rel)))

################ CHECK DISPERSION ################
dispersion <- betadisper(bc_dist, sample_data(Microbiome_Rel)$group)
anova(dispersion)

################ PCoA ################
ordination = ordinate(Microbiome_Rel, method="PCoA", distance=bc_dist)
pCoA_plot <- plot_ordination(Microbiome_Rel, ordination, color="group") +
  scale_color_manual(values=group_color) +
  theme_bw()
pCoA_plot

################ PCoA ################
ordination <- ordinate(Microbiome_Rel, method="PCoA", distance=bc_dist)

pCoA_plot <- plot_ordination(
  Microbiome_Rel,
  ordination,
  color="group"
) +
  geom_point(size=4, alpha=0.9) +
  stat_ellipse(type="t", linetype=2, linewidth=0.8) +
  scale_color_manual(values=group_color) +
  theme_classic(base_size = 14) +
  theme(
    legend.title = element_blank(),
    legend.position = "right"
  ) +
  coord_equal()

# ADD IT HERE
pCoA_plot <- pCoA_plot +
  labs(
    title = "PCoA (Bray-Curtis)",
    x = paste0("PCoA1 (", round(ordination$values$Relative_eig[1]*100,1), "%)"),
    y = paste0("PCoA2 (", round(ordination$values$Relative_eig[2]*100,1), "%)")
  )

pCoA_plot

################ NMDS ################
ordination_nmds <- ordinate(Microbiome_Rel, method="NMDS", distance="bray")
NMDS_plot <- plot_ordination(Microbiome_Rel, ordination_nmds, color="group") +
  scale_color_manual(values=group_color) +
  theme_bw()
NMDS_plot

ordination_nmds <- ordinate(Microbiome_Rel, method="NMDS", distance="bray")

NMDS_plot <- plot_ordination(
  Microbiome_Rel,
  ordination_nmds,
  color="group"
) +
  geom_point(size=4, alpha=0.9) +
  stat_ellipse(type="t", linetype=2, linewidth=0.8) +
  scale_color_manual(values=group_color) +
  theme_classic(base_size = 14) +
  theme(
    legend.title = element_blank(),
    legend.position = "right"
  ) +
  coord_equal()

NMDS_plot

ggsave("PCoA_plot.tiff", pCoA_plot, width=7, height=5, dpi=600)
ggsave("NMDS_plot.tiff", NMDS_plot, width=7, height=5, dpi=600)

###########TABLE OF COMPOSITIONS#############

family_summary <- FamilyRel %>%
  group_by(Family, Condition) %>%
  summarise(Mean_Abundance = mean(Abundance), .groups="drop") %>%
  arrange(Condition, desc(Mean_Abundance)) %>%
  filter(Family != "Other") %>%
  rename(Taxon = Family) %>%
  mutate(Level = "Family")

genus_summary <- GenusRel %>%
  group_by(Genus, Condition) %>%
  summarise(Mean_Abundance = mean(Abundance), .groups="drop") %>%
  arrange(Condition, desc(Mean_Abundance)) %>%
  filter(Genus != "Other") %>%
  rename(Taxon = Genus) %>%
  mutate(Level = "Genus")

combined_table <- bind_rows(family_summary, genus_summary) %>%
  select(Level, Taxon, Condition, Mean_Abundance) %>%
  arrange(Level, Condition, desc(Mean_Abundance))

print(combined_table)
write.csv(combined_table, "taxonomy_abundance_table.csv", row.names=FALSE)

############COLLAPSED COMPOSITION PLOTS#####################

FamilyRel_collapsed <- FamilyRel %>%
  group_by(Disease, Family) %>%
  summarise(Mean_Abundance = mean(Abundance, na.rm = TRUE), .groups = "drop") %>%
  group_by(Disease) %>%
  mutate(Mean_Abundance = Mean_Abundance / sum(Mean_Abundance)) # This line fixes the bar height

ggplot(FamilyRel_collapsed, aes(x = Disease, y = Mean_Abundance, fill = Family)) +
  geom_col() + # Use geom_col() instead of geom_bar(stat="identity")
  scale_y_continuous(labels = scales::percent) + # Optional: change 1.0 to 100%
  ylab("Relative Abundance (%)") +
  xlab("") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_fill_manual(values = family_colors)


#################################### DESEQ2 + TOP TAXA PLOT #########################################

Microbiome_genus <- tax_glom(Microbiome_prune, taxrank = "Genus")

dds_asd <- phyloseq_to_deseq2(Microbiome_genus, ~ Disease)

dds_asd <- estimateSizeFactors(dds_asd, type = "poscounts")

dds_asd <- DESeq(dds_asd, fitType = "local")

res_asd <- results(dds_asd, contrast = c("Disease","ASD","Control"))
res_asd <- res_asd[order(res_asd$padj, na.last=NA), ]

#################################### SIGNIFICANT TAXA #########################################

alpha <- 0.05

sigtab <- res_asd[res_asd$padj < alpha, ]


sigtab_annot <- cbind(
  as(sigtab, "data.frame"),
  as(tax_table(Microbiome_genus)[rownames(sigtab), ], "matrix")
)


sigtabgen <- subset(sigtab_annot, log2FoldChange > 1 | log2FoldChange < -1)


sigtabgen <- subset(sigtabgen, !is.na(Genus))

#################################### ORDERING #########################################


x <- tapply(sigtabgen$log2FoldChange, sigtabgen$Phylum, max)
x <- sort(x, TRUE)
sigtabgen$Phylum <- factor(as.character(sigtabgen$Phylum), levels = names(x))


x <- tapply(sigtabgen$log2FoldChange, sigtabgen$Genus, max)
x <- sort(x, TRUE)
sigtabgen$Genus <- factor(as.character(sigtabgen$Genus), levels = names(x))

#################################### PLOT #########################################

library(ggplot2)

deseq_plot <- ggplot(sigtabgen, aes(y = Genus, x = log2FoldChange, colour = Phylum)) + 
  geom_vline(xintercept = 0, colour = "grey50", linewidth = 0.5) +
  geom_point(size = 4) + 
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  ) +
  xlab("Log2 Fold Change (ASD vs Control)") +
  ylab("Genus")

deseq_plot + coord_flip()



#################################### ASD SEVERITY ANALYSIS #########################################

Microbiome_prune <- subset_samples(
  Microbiome_prune,
  Condition %in% c("ASD_Mild","ASD_Severe")
)


sample_data(Microbiome_prune)$severity <- factor(
  sample_data(Microbiome_prune)$Condition,
  levels = c("ASD_Mild","ASD_Severe")
)

Microbiome_Rel <- transform_sample_counts(Microbiome_prune, function(x) x / sum(x))

#################################### ALPHA DIVERSITY #########################################

group_color <- c("red","blue")

alpha_all <- plot_richness(
  Microbiome_prune,
  x="severity",
  color="severity",
  measures=c("InvSimpson","Shannon","Observed")
) +
  scale_color_manual(values=group_color) +
  geom_boxplot(lwd=1,width=0.5) +
  geom_point()

alpha_all

# Wilcoxon test (2 groups)
alpha_stats <- alpha_all +
  stat_compare_means(method="wilcox.test", label="p.format") +
  theme_bw()

alpha_stats

#################################### SAVE ALPHA TABLE #########################################

richness.rare <- cbind(
  estimate_richness(
    Microbiome_prune,
    measures=c('Observed','Shannon',"Simpson","InvSimpson")
  ),
  sample_data(Microbiome_prune)$severity
)

write.csv(richness.rare, file="alphaDiv_ASD_severity.csv")


library(dplyr)

alpha_df <- estimate_richness(
  Microbiome_prune,
  measures = c("Observed","Shannon","Simpson","InvSimpson")
)

alpha_df$severity <- sample_data(Microbiome_prune)$severity


#################### DESCRIPTIVE STATISTICS ########################

alpha_descriptives <- alpha_df %>%
  group_by(severity) %>%
  summarise(
    n = n(),
    Observed_median = median(Observed),
    Observed_IQR = IQR(Observed),
    Shannon_median = median(Shannon),
    Shannon_IQR = IQR(Shannon),
    InvSimpson_median = median(InvSimpson),
    InvSimpson_IQR = IQR(InvSimpson)
  )


alpha_descriptives


########################## WILCOXON TESTS ################################

wilcox_observed <- wilcox.test(Observed ~ severity, data = alpha_df)
wilcox_shannon <- wilcox.test(Shannon ~ severity, data = alpha_df)
wilcox_invsimpson <- wilcox.test(InvSimpson ~ severity, data = alpha_df)

alpha_stats_table <- data.frame(
  Metric = c("Observed","Shannon","InvSimpson"),
  W = c(wilcox_observed$statistic,
        wilcox_shannon$statistic,
        wilcox_invsimpson$statistic),
  p_value = c(wilcox_observed$p.value,
              wilcox_shannon$p.value,
              wilcox_invsimpson$p.value),
  n_mild = sum(alpha_df$severity == "ASD_Mild"),
  n_severe = sum(alpha_df$severity == "ASD_Severe")
)

alpha_stats_table

write.csv(alpha_descriptives, "alpha_descriptives_ASD_severity.csv")
write.csv(alpha_stats_table, "alpha_stats_ASD_severity.csv")

#################################### BETA DIVERSITY #########################################

# Bray Curtis distance
bc_dist <- phyloseq::distance(Microbiome_Rel, method="bray", weighted=FALSE)

# PERMANOVA
adonis2(bc_dist ~ Condition, data=data.frame(sample_data(Microbiome_Rel)))

# Check dispersion (important assumption for PERMANOVA)
dispersion <- betadisper(bc_dist, sample_data(Microbiome_Rel)$Condition)
anova(dispersion)

#################################### PCoA #########################################

ordination <- ordinate(Microbiome_Rel, method="PCoA", distance=bc_dist)

pCoA_plot <- plot_ordination(
  Microbiome_Rel,
  ordination,
  color="severity"
) +
  scale_color_manual(values=group_color) +
  theme_bw()

pCoA_plot

#################################### PCoA #########################################

ordination <- ordinate(Microbiome_Rel, method="PCoA", distance=bc_dist)

condition_color <- c(
  "ASD_Mild" = "red",
  "ASD_Severe" = "blue"
)

pCoA_plot <- plot_ordination(
  Microbiome_Rel,
  ordination,
  color="severity"
) +
  geom_point(size=4, alpha=0.9) +
  stat_ellipse(type="t", linetype=2, linewidth=0.8) +
  scale_color_manual(values=condition_color) +
  theme_classic(base_size = 14) +
  theme(
    legend.title = element_blank(),
    legend.position = "right"
  ) +
  coord_equal() +
  labs(
    title = "PCoA (Bray-Curtis)",
    x = paste0("PCoA1 (", round(ordination$values$Relative_eig[1]*100,1), "%)"),
    y = paste0("PCoA2 (", round(ordination$values$Relative_eig[2]*100,1), "%)")
  )

pCoA_plot

ggsave(
  "PCoA_ASD_severity.tiff",
  pCoA_plot,
  width = 7,
  height = 5,
  dpi = 600
)

#################################### NMDS #########################################

condition_color <- c(
  "ASD_Mild" = "red",
  "ASD_Severe" = "blue"
)

pCoA_plot <- plot_ordination(
  Microbiome_Rel,
  ordination,
  color = "Condition"
) +
  scale_color_manual(values = condition_color) +
  theme_bw()
pCoA_plot

ordination_nmds <- ordinate(Microbiome_Rel, method="NMDS", distance="bray")

NMDS_plot <- plot_ordination(
  Microbiome_Rel,
  ordination_nmds,
  color = "Condition"
) +
  scale_color_manual(values = condition_color) +
  theme_bw()
NMDS_plot

#################################### NMDS #########################################

ordination_nmds <- ordinate(Microbiome_Rel, method="NMDS", distance="bray")

NMDS_plot <- plot_ordination(
  Microbiome_Rel,
  ordination_nmds,
  color="severity"
) +
  geom_point(size=4, alpha=0.9) +
  stat_ellipse(type="t", linetype=2, linewidth=0.8) +
  scale_color_manual(values=condition_color) +
  theme_classic(base_size = 14) +
  theme(
    legend.title = element_blank(),
    legend.position = "right"
  ) +
  coord_equal()

NMDS_plot

ggsave(
  "NMDS_ASD_severity.tiff",
  NMDS_plot,
  width = 7,
  height = 5,
  dpi = 600
)

################################ FAMILY LEVEL BAR PLOT ################################

Rel_family <- tax_glom(Microbiome_Rel, taxrank = 'Family', NArm=FALSE) 
FamilyRel <- psmelt(Rel_family) 
FamilyRel$Family <- as.character(FamilyRel$Family) 
FamilyRel$Family[FamilyRel$Abundance < 0.05] <- "Other"

high_contrast_palette <- c(
  "#E63946", "#2196F3", "#FF9800", "#4CAF50", "#9C27B0",
  "#00BCD4", "#FF5722", "#8BC34A", "#673AB7", "#FFC107",
  "#F06292", "#26C6DA", "#D4E157", "#AB47BC", "#42A5F5",
  "#FF7043", "#66BB6A", "#EC407A", "#29B6F6", "#FFCA28",
  "#EF5350", "#26A69A", "#7E57C2", "#9CCC65", "#FFA726",
  "#B71C1C", "#1A237E", "#E65100", "#1B5E20", "#4A148C",
  "#006064", "#BF360C", "#33691E", "#311B92", "#F57F17",
  "#880E4F", "#0D47A1", "#E84393", "#558B2F", "#6A1B9A",
  "#00838F", "#D84315", "#558B2F", "#4527A0", "#F9A825",
  "#C62828", "#283593", "#FF6F00", "#2E7D32", "#6A1B9A"
)

n_families <- length(unique(FamilyRel$Family))
family_colors <- setNames(high_contrast_palette[1:n_families], unique(FamilyRel$Family))

ggplot(FamilyRel, aes(x = Sample, y = Abundance/3, fill = Family)) +
  geom_bar(stat="identity") +
  ylab("Relative abundance") +
  xlab("") +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust=1, size=7)
  ) +
  facet_wrap(~ Condition, scales="free_x", nrow=1) +
  scale_fill_manual(values=family_colors)

library(dplyr)
library(ggplot2)

library(scales) # Needed for the percent labels

# 1. Filter, calculate means, and re-normalize for 100% stack
FamilyRel_collapsed <- FamilyRel %>%
  filter(Condition %in% c("ASD_Severe", "ASD_Mild")) %>% 
  group_by(Condition, Family) %>%
  summarise(Mean_Abundance = mean(Abundance, na.rm = TRUE), .groups = "drop") %>%
  group_by(Condition) %>%
  mutate(Mean_Abundance = Mean_Abundance / sum(Mean_Abundance)) %>%
  mutate(Family = ifelse(Mean_Abundance < 0.05, "Other", Family)) %>%
  group_by(Condition, Family) %>%
  summarise(Mean_Abundance = sum(Mean_Abundance), .groups = "drop") %>%
  group_by(Condition) %>%
  mutate(Mean_Abundance = Mean_Abundance / sum(Mean_Abundance))
# 2. Plotting code
ggplot(FamilyRel_collapsed, aes(x = Condition, y = Mean_Abundance, fill = Family)) +
  geom_col() + 
  scale_y_continuous(labels = scales::percent) + 
  ylab("Mean Relative Abundance (%)") +
  xlab("Severity Group") +
  theme_bw() +
  theme(axis.text.x = element_text(size = 10, face = "bold")) +
  scale_fill_manual(values = family_colors)


Family_Percentages <- FamilyRel_collapsed %>%
 
  mutate(Percentage = Mean_Abundance * 100) %>%

  select(Condition, Family, Percentage)


write.csv(Family_Percentages, "Severity_Family_Relative_Abundance_Averages.csv", row.names = FALSE)


print(Family_Percentages)



################################ Genus LEVEL BAR PLOT ################################

Rel_genus <- tax_glom(Microbiome_Rel, taxrank = 'Genus', NArm=FALSE) 
GenusRel <- psmelt(Rel_genus) 
GenusRel$Genus <- as.character(GenusRel$Genus) 
GenusRel$Genus[GenusRel$Abundance < 0.05] <- "Other"

n_genera <- length(unique(GenusRel$Genus))
genus_colors <- setNames(high_contrast_palette[1:n_genera], unique(GenusRel$Genus))

ggplot(GenusRel, aes(x = Sample, y = Abundance/3, fill = Genus)) +
  geom_bar(stat="identity") +
  ylab("Relative abundance") +
  xlab("") +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust=1, size=7)
  ) +
  facet_wrap(~ Condition, scales="free_x", nrow=1) +
  scale_fill_manual(values=genus_colors)

#################################### DESEQ2 + TOP TAXA (ASD SEVERITY) #########################################


dds_severity <- phyloseq_to_deseq2(Microbiome_genus, ~ Condition)


dds_severity <- estimateSizeFactors(dds_severity, type = "poscounts")


dds_severity <- DESeq(dds_severity, fitType = "local")


res_severity <- results(dds_severity, contrast = c("Condition","ASD_Severe","ASD_Mild"))
res_severity <- res_severity[order(res_severity$padj, na.last=NA), ]

#################################### SIGNIFICANT TAXA #########################################

alpha <- 0.05

sigtab <- res_severity[res_severity$padj < alpha, ]


sigtab_annot <- cbind(
  as(sigtab, "data.frame"),
  as(tax_table(Microbiome_genus)[rownames(sigtab), ], "matrix")
)


sigtabgen <- subset(sigtab_annot, log2FoldChange > 1 | log2FoldChange < -1)


sigtabgen <- subset(sigtabgen, !is.na(Genus))

#################################### ORDERING #########################################


x <- tapply(sigtabgen$log2FoldChange, sigtabgen$Phylum, max)
x <- sort(x, TRUE)
sigtabgen$Phylum <- factor(as.character(sigtabgen$Phylum), levels = names(x))


x <- tapply(sigtabgen$log2FoldChange, sigtabgen$Genus, max)
x <- sort(x, TRUE)
sigtabgen$Genus <- factor(as.character(sigtabgen$Genus), levels = names(x))

#################################### PLOT #########################################

library(ggplot2)

deseq_plot_severity <- ggplot(sigtabgen, aes(y = Genus, x = log2FoldChange, colour = Phylum)) + 
  geom_vline(xintercept = 0, colour = "grey50", linewidth = 0.5) +
  geom_point(size = 4) + 
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  ) +
  xlab("Log2 Fold Change (Severe vs Mild)") +
  ylab("Genus")

deseq_plot_severity + coord_flip()
