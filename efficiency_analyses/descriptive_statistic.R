#efficiency may have two measures, sensitivity, or percentage of exons covered by at least one read
#and specificity, or the percentage of mapped reads that map to target regions

# Load data
#This data includes the name of the suborder and superfamily of each species in the taxon sampling, as well as the number of genes extracted indicated in the hybpiper starts results hybpiper_stats.tsv
data <- read.csv("D:/Sofi/Desktop/sofia/AMNH - postdoc/jessica's Lab/Postdoc project/4. target capture valdes et al/results/nucleotides_new_baitset/exons.csv")

# Set working directory
setwd("D:/Sofi/Desktop/sofia/AMNH - postdoc/jessica's Lab/Postdoc project/4. target capture valdes et al/results/nucleotides_new_baitset/R analyses/")

# Order Species by Suborder and Superfamily
data2 <- data[order(data$Suborder, data$Superfamily), ]
data2$Suborder <- factor(data2$Suborder, levels = unique(data2$Suborder))
data2$Superfamily <- factor(data2$Superfamily, levels = unique(data2$Superfamily))
data2$Species <- factor(data2$Species, levels = unique(data2$Species))

# Load libraries
library(ggplot2)
library(scales)
library(RColorBrewer)
library(colorspace)

# --- Step 1: Clean the data ---
data2$Superfamily <- as.character(trimws(data2$Superfamily))
data2$Suborder <- as.character(trimws(data2$Suborder))

# --- Step 2: Build lookup table ---
lookup_df <- unique(data2[, c("Superfamily", "Suborder")])

# Create the named lookup
suborder_lookup <- setNames(lookup_df$Suborder, lookup_df$Superfamily)

# --- Step 3: Define base colors ---
superfamilies <- sort(unique(data2$Superfamily))
base_colors <- hue_pal()(length(superfamilies))
names(base_colors) <- superfamilies

# --- Step 4: Build color mapping ---
color_mapping <- sapply(names(base_colors), function(sf) {
  suborder <- suborder_lookup[[sf]]
  
  # Debug print
  cat("Superfamily:", sf, "→ Suborder:", suborder, "\n")
  
  if (suborder == "Doridina") {
    lighten(base_colors[sf], amount = 0.5)
  } else if (suborder == "Cladobranchia") {
    darken(base_colors[sf], amount = 0.3)
  } else {
    "grey60"
  }
})

names(color_mapping) <- names(base_colors)

print("FINAL COLOR MAPPING:")
print(color_mapping)

# --- Step 5: Factor levels for plotting ---
data2$Superfamily <- factor(data2$Superfamily, levels = names(color_mapping))
data2$Species <- factor(data2$Species, levels = unique(data2$Species))

# Set Superfamily factor order based on Superfamily's appearance in data2 (for legend order)
sf_order <- unique(data2$Superfamily[order(data2$Species)])

# Reorder Superfamily factor to match bar order
data2$Superfamily <- factor(data2$Superfamily, levels = sf_order)

##Lets do it in percentage, to show the sensitivity of the dataset. So all numbers /304 original dataset

data2$sensitivity <- (data2$genesExtracted / 304) * 100 #304 is the number of genes in my baitset

bar_plot = ggplot(data2, aes(x = Species, y = sensitivity, fill = Superfamily)) +
  geom_bar(stat = "identity", width = 0.7) +
  labs(title = "Capture Sensitivity",
       x = "Species",
       y = "% of Genes Extracted") +
  scale_fill_manual(values = color_mapping, guide = guide_legend(reverse = TRUE)) +  # Reverse legend order
  theme_minimal() +
  coord_flip() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.text.y = element_text(size = 10),
        legend.text = element_text(size = 14),         # Adjust legend text size
        legend.title = element_text(size = 13),
        plot.title = element_text(hjust = 0.5))        # Adjust legend title size

ggsave("barplot_exons_sensitivity_colors.pdf", bar_plot, width = 10, height = 13, bg = "white")
ggsave("barplot_exons_sensitivity_colors.png", bar_plot, width = 10, height = 13, bg = "white")


###Now lets see the specificity, species vs hybpiper stat PctOnTarget

specif <- read.csv("D:/Sofi/Desktop/sofia/AMNH - postdoc/jessica's Lab/Postdoc project/4. target capture valdes et al/results/nucleotides_new_baitset/hybpiper_stats_new_baitset.csv", header = TRUE, sep = ",")

data2$PctOnTarget <- specif$PctOnTarget[match(data2$Species, specif$Species)]

###Layton 2020 are exceptions with so many reads, so exclude outliers in my analyses
species_to_remove <- c("Ardeadoris_egretta", "Chromodoris_magnifica", "Chromodoris_westraliensis", "Doriprismatica_atromarginata", "Goniobranchus_fidelis")

data2 <- data2[!data2$Species %in% species_to_remove, ]

bar_plot = ggplot(data2, aes(x = Species, y = PctOnTarget, fill = Superfamily)) +
  geom_bar(stat = "identity", width = 0.7) +
  labs(title = "Capture Specificity",
       x = "Species",
       y = "% reads mapped to targets") +
  scale_fill_manual(values = color_mapping, guide = guide_legend(reverse = TRUE)) +
  coord_flip() +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 0),
        axis.text.y = element_text(size = 10),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 13),
        plot.title = element_text(hjust = 0.5))

ggsave("barplot_exons_specificity_colors.pdf", bar_plot, width = 10, height = 13, bg = "white")
ggsave("barplot_exons_specificity_colors.png", bar_plot, width = 10, height = 13, bg = "white")

#to loop for figure
bar_plot = ggplot(data2, aes(x = Species, y = PctOnTarget, fill = Superfamily)) +
  geom_bar(stat = "identity", width = 0.7) +
  labs(title = "Capture Specificity",
       x = "Species",
       y = "% reads mapped to targets") +
  scale_fill_manual(values = color_mapping, guide = guide_legend(reverse = TRUE)) +
  coord_flip() +
  scale_x_discrete(position = "top") +   # Move species labels to the right
  scale_y_reverse() +                   # Reverse the data direction (right to left)
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 0),
        axis.text.y = element_text(size = 10),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 13),
        plot.title = element_text(hjust = 0.5))
ggsave("barplot_exons_specificity_colorsv2.pdf", bar_plot, width = 10, height = 13, bg = "white")
ggsave("barplot_exons_specificity_colorsv2.png", bar_plot, width = 10, height = 13, bg = "white")
