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
    "grey20"
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


#we could see a difference in the percentate of genes recovered. To be more specific, plot the species and the specific genes extracted with the aim to see which genes could be missing in Cladobranchia vs Doridina.

# Load required library
library(pheatmap)
library(dplyr)

# Read the full CSV (no header, preserve everything)
full_data <- read.csv(
  "D:/Sofi/Desktop/sofia/AMNH - postdoc/jessica's Lab/Postdoc project/4. target capture valdes et al/results/nucleotides_new_baitset/gene_presence_matrix.csv",
  header = FALSE,
  stringsAsFactors = FALSE
)

# Extract suborder (first row) and species names (second row), skipping first column
suborder_vec <- as.character(full_data[1, -1])
species_names <- as.character(full_data[2, -1])

# Create annotation dataframe
annotation_col <- data.frame(Suborder = suborder_vec)
rownames(annotation_col) <- species_names

# Define suborder colors for annotations
suborder_colors <- c(
  "Doridina" = "pink",
  "Cladobranchia" = "purple",
  "Outgroup" = "gray20"
)
annotation_colors <- list(Suborder = suborder_colors)

# Extract gene presence/absence data
gene_data <- full_data[-c(1, 2), ]  # Remove first 2 rows
gene_names <- gene_data[, 1]       # First column = gene names
gene_data <- gene_data[, -1]       # Remove gene names column

# Convert to matrix and to numeric values
presence_matrix <- as.matrix(gene_data)
presence_matrix <- apply(presence_matrix, 2, as.numeric)

# Restore row and column names
rownames(presence_matrix) <- gene_names
colnames(presence_matrix) <- species_names

# Sort species by suborder: Doridina → Cladobranchia → Outgroup
sorted_species <- annotation_col %>%
  mutate(Species = rownames(.)) %>%
  arrange(factor(Suborder, levels = c("Doridina", "Cladobranchia", "Outgroup"))) %>%
  pull(Species)

# Reorder matrix and annotations accordingly
presence_matrix <- presence_matrix[, sorted_species]
annotation_col <- annotation_col[sorted_species, , drop = FALSE]

# Sort genes by number of species present (descending)
gene_presence_counts <- rowSums(presence_matrix, na.rm = TRUE)
presence_matrix <- presence_matrix[order(-gene_presence_counts), ]

svg("heatmap_presence-absence.svg", width = 13, height = 25)
pheatmap(presence_matrix,
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         annotation_col = annotation_col,
         annotation_colors = annotation_colors,
         color = c("white", "black"),
         legend_breaks = c(0, 1),
         legend_labels = c("Absent", "Present"),
         fontsize_row = 5,
         fontsize_col = 8,
         show_rownames = TRUE,
         border_color = NA)
dev.off()

#Mark teasdale genes
library(ComplexHeatmap)
library(circlize)  # for color ramp

# Read your Teasdale gene list
teasdale_genes <- readLines("D:/Sofi/Desktop/sofia/AMNH - postdoc/jessica's Lab/Postdoc project/4. target capture valdes et al/results/nucleotides_new_baitset/Teasdale IDs.txt")

# Ensure gene names are character
gene_names <- rownames(presence_matrix)

# Create color vector for gene labels
gene_label_colors <- ifelse(gene_names %in% teasdale_genes, "red", "black")


# 2. Update presence matrix: 0 = absent, 1 = present non-Teasdale, 2 = present Teasdale
modified_matrix <- presence_matrix  # Copy original matrix
gene_names <- rownames(presence_matrix)
for (gene in teasdale_genes) {
  if (gene %in% gene_names) {
    gene_index <- which(rownames(modified_matrix) == gene)
    modified_matrix[gene_index, modified_matrix[gene_index, ] == 1] <- 2
  }
}

# 3. Create color mapping: 0 = white, 1 = darkblue, 2 = red
col_fun <- c("0" = "white", "1" = "lightskyblue1", "2" = "steelblue")

# 4. Set row label colors
gene_label_colors <- ifelse(gene_names %in% teasdale_genes, "seagreen", "black")
annotation_col$Suborder <- factor(annotation_col$Suborder, levels = c("Doridina", "Cladobranchia", "Outgroup"))
# 5. Annotation for species
ha_column <- HeatmapAnnotation(
  Suborder = annotation_col$Suborder,
  col = list(Suborder = c("Doridina" = "pink", "Cladobranchia" = "purple", "Outgroup" = "gray")),
  annotation_name_side = "left"
)

# 6. Heatmap plotting (with split by suborder for a visible dividing line)
draw_heatmap <- function(file, device_fun) {
  device_fun(file, width = 13, height = 30)
  
  ht <- Heatmap(modified_matrix,
                name = "Presence",
                col = col_fun,
                cluster_rows = FALSE,
                cluster_columns = FALSE,
                show_column_names = TRUE,
                column_names_gp = gpar(fontsize = 8),
                row_names_gp = gpar(fontsize = 7, col = gene_label_colors),
                top_annotation = ha_column,
                show_row_names = TRUE,
                column_split = annotation_col$Suborder)  # Add vertical split lines
  
  draw(ht)  # This actually renders the plot
  dev.off()
}

# 7. Save all formats
draw_heatmap("heatmap_with_teasdale2.pdf", pdf)
draw_heatmap("heatmap_with_teasdale2.png", function(...) png(..., units = "in", res = 300))
draw_heatmap("heatmap_with_teasdale2.svg", svg)

###lets try to clean it.

draw_heatmap2 <- function(file, device_fun) {
  device_fun(file, width = 10, height = 15)
ht2 <- Heatmap(modified_matrix,
              name = "Presence",
              col = col_fun,
              cluster_rows = FALSE,
              cluster_columns = FALSE,
              show_column_names = FALSE,
              column_names_gp = gpar(fontsize = 2),
              row_names_gp = gpar(fontsize = 2, col = gene_label_colors),
              top_annotation = ha_column,
              show_row_names = FALSE,
              column_split = annotation_col$Suborder)  # Add vertical split lines

draw(ht2)  # This actually renders the plot
dev.off()
}

# 7. Save all formats
draw_heatmap2("heatmap_with_teasdale3.pdf", pdf)
draw_heatmap2("heatmap_with_teasdale3.png", function(...) png(..., units = "in", res = 300))
draw_heatmap2("heatmap_with_teasdale3.svg", svg)


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
