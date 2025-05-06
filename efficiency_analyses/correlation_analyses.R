####analyses using distances of species with target dataset, query vs target
#and then, correlation between distances and number of exons

#Reading data
data <- read.csv("D:/Sofi/Desktop/sofia/AMNH - postdoc/jessica's Lab/Postdoc project/4. target capture valdes et al/results/nucleotides_new_baitset/exons2.csv")

#reading distances calculated
uncorr_distance <- read.csv("D:/Sofi/Desktop/sofia/AMNH - postdoc/jessica's Lab/Postdoc project/4. target capture valdes et al/results/nucleotides_new_baitset/p_distance_matrix.csv", header = TRUE, sep = ",", row.names = 1)
uncorr_distance <- as.data.frame(sapply(uncorr_distance, as.numeric))

mean_uncorrected <- colMeans(uncorr_distance, na.rm = TRUE)
# Add the mean values as the last row
uncorr_distance <- rbind(uncorr_distance, mean_uncorrected)

# Assign a name to the new row
rownames(uncorr_distance)[nrow(uncorr_distance)] <- "Mean_PDistance"

# Now, we want to extract the Mean_PDistance row
mean_pdistance_vector <- mean_uncorrected

# Create a data frame with species and mean p-distance
mean_pdistance_df <- data.frame(Species = colnames(uncorr_distance),
                                Mean_PDistance = mean_pdistance_vector)

# Merge the two data frames by species
merged_data <- merge(data, mean_pdistance_df, by = "Species")

write.csv(merged_data, "merged_data.csv", row.names = FALSE)

####test normality

merged_data$sensitivity <- (merged_data$genesExtracted / 304) * 100 #304 is the number of genes in my baitset

# Shapiro-Wilk test for normality
shapiro.test(merged_data$sensitivity)
shapiro.test(merged_data$Mean_PDistance)

#data are normal in Shapiro-wilk test when W (messure how close your data distribution is to nnormal) is close to 1 
#and p-value greater than 0.05 (probability that the observed data could come from normal distribution by chance)

## not normal, so no parametric analyses are fine. like Kendall correlation

# Perform Kendall rank correlation test
kendall_correlation_result <- cor.test(merged_data$Mean_PDistance, merged_data$sensitivity, method = "kendall")

# Print the Kendall correlation result
print(kendall_correlation_result)

#tau, kendall's tau coefficient, +1 possitive association, -1 negatice association 
#p-value threshold 0.05, mine 2.2e-16, correlation is statistically significant.
#z-score, further from 0 evidence against the null hypothesis, -8.9 strong evidence against no correlation

library(ggplot2)
library(scales)
merged_data$Superfamily <- as.character(trimws(merged_data$Superfamily))
merged_data$Suborder <- as.character(trimws(merged_data$Suborder))

# Create a lookup table for Superfamily to Suborder
lookup_df <- unique(merged_data[, c("Superfamily", "Suborder")])
suborder_lookup <- setNames(lookup_df$Suborder, lookup_df$Superfamily)

superfamilies <- sort(unique(merged_data$Superfamily))
base_colors <- hue_pal()(length(superfamilies))
names(base_colors) <- superfamilies

# Manually define colors based on Superfamily and Suborder
# Assign colors based on the Suborder classification for each Superfamily
color_mapping <- sapply(names(base_colors), function(sf) {
  suborder <- suborder_lookup[[sf]]
  
  if (suborder == "Doridina") {
    lighten(base_colors[sf], amount = 0.5)
  } else if (suborder == "Cladobranchia") {
    darken(base_colors[sf], amount = 0.3)
  } else {
    "grey60"
  }
})
names(color_mapping) <- names(base_colors)

merged_data$Superfamily <- factor(merged_data$Superfamily, levels = names(color_mapping))
merged_data$Species <- factor(merged_data$Species, levels = unique(merged_data$Species))

# Set Superfamily factor order based on Superfamily's appearance in merged_data (for legend order)
sf_order <- unique(data2$Superfamily[order(data2$Species)])

# Reorder Superfamily factor to match the desired order
merged_data$Superfamily <- factor(merged_data$Superfamily, levels = sf_order)

# Create the ggplot object with the custom color palette
gg_plot <- ggplot(merged_data, aes(x = Mean_PDistance, y = genesExtracted, color = Superfamily)) +
  geom_point(size = 4, alpha = 0.7) +
  labs(title = paste("Kendall Rank Correlation between genesExtracted and Mean_PDistance (τ =", round(kendall_correlation_result$estimate, 2), ")"),
       x = "Mean P-Distance", y = "Genes Extracted") +
  theme_minimal() +  # Minimal theme with white background
  scale_color_manual(values = color_mapping, guide = guide_legend(reverse = TRUE)) +  # Apply the custom color mapping
  theme(panel.background = element_rect(fill = "white"),  # Ensure the plot background is white
        plot.background = element_rect(fill = "white"),    # Ensure the overall background is white
        legend.background = element_rect(fill = "white"),  # Optional: make legend background white
        axis.text = element_text(size = 10),  # Adjust axis text size
        axis.title = element_text(size = 12),  # Adjust axis title size
        plot.title = element_text(size = 14, face = "bold"))  # Adjust plot title size

# Save the plot with a white background
ggsave("ggplot_distancevsgenes_color.png", plot = gg_plot, width = 8.5, height = 6, dpi = 300)
ggsave("ggplot_distancevsgenes_color.pdf", plot = gg_plot, width = 8.5, height = 6, dpi = 300)

###My genes decrease, but when? Is there a specific genetic distance after which exon numbers start to drop off more quickly?

##more clean, only one line with two clear linear trends, one before and one after the break. with a single breakpoint
#Piecewise (Segmented) Regression

# Load packages
library(segmented)
library(ggplot2)
library(scales)

### %of genes vs distance

merged_data$sensitivity <- (merged_data$genesExtracted / 304) * 100 #304 is the number of genes in my baitset

# Fit initial linear model
lm_fit <- lm(sensitivity ~ Mean_PDistance, data = merged_data)

# Fit segmented regression with one breakpoint
seg_fit <- segmented(lm_fit, seg.Z = ~Mean_PDistance)

# View estimated breakpoint
summary(seg_fit)
breakpoint <- seg_fit$psi[2]  # the estimated distance value

# Create predicted values for plotting
merged_data$seg_pred <- predict(seg_fit)

gg_plot <- ggplot(merged_data, aes(x = Mean_PDistance, y = sensitivity, color = Superfamily)) +
  geom_point(alpha = 0.6, size = 4) +
  geom_line(aes(y = seg_pred), color = "grey", size = 1.2, alpha = 0.5) +
  geom_vline(xintercept = breakpoint, linetype = "dashed", color = "red", size = 1) +
  labs(
    title = paste("Piecewise Regression: % Genes Count vs Genetic Distance"),
    subtitle = paste("Breakpoint at Mean P-Distance ≈", round(breakpoint, 3)),
    x = "Mean P-Distance", y = "% Genes Extracted"
  ) +
  scale_color_manual(values = color_mapping, guide = guide_legend(reverse = TRUE)) +
  theme_minimal() +
  theme(
    panel.background = element_rect(fill = "white"),
    plot.background = element_rect(fill = "white"),
    legend.background = element_rect(fill = "white"),
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 13))

# Save the plot
ggsave("ggplot_piecewise_percentage_color.pdf", plot = gg_plot, width = 8, height = 6, dpi = 300)
ggsave("ggplot_piecewise_percentage_color.png", plot = gg_plot, width = 8, height = 6, dpi = 300)


############what about the correlation between PcTOnTarget and distances

merged_data$PctOnTarget <- specif$PctOnTarget[match(merged_data$Species, specif$Species)]

###Layton 2020 are exceptions with so many reads, so exclude outliers in my analyses
species_to_remove <- c("Ardeadoris_egretta", "Chromodoris_magnifica", "Chromodoris_westraliensis", "Doriprismatica_atromarginata", "Goniobranchus_fidelis")

merged_data <- merged_data[!merged_data$Species %in% species_to_remove, ]

####test normality

# Shapiro-Wilk test for normality
shapiro.test(merged_data$PctOnTarget)
shapiro.test(merged_data$Mean_PDistance)

#data are normal in Shapiro-wilk test when W (messure how close your data distribution is to nnormal) is close to 1 
#and p-value greater than 0.05 (probability that the observed data could come from normal distribution by chance)

##not normal, so no parametric analyses are fine. like Kendall correlation

# Perform Kendall rank correlation test
kendall_correlation_result <- cor.test(merged_data$Mean_PDistance, merged_data$PctOnTarget, method = "kendall")

# Print the Kendall correlation result
print(kendall_correlation_result)

#tau, kendall's tau coefficient, +1 possitive association, -1 negatice association 
#p-value threshold 0.05, mine 2.2e-16, correlation is statistically significant.
#z-score, further from 0 evidence against the null hypothesis, -8.9 strong evidence against no correlation

library(ggplot2)
library(scales)
merged_data$Superfamily <- as.character(trimws(merged_data$Superfamily))
merged_data$Suborder <- as.character(trimws(merged_data$Suborder))

# Create a lookup table for Superfamily to Suborder
lookup_df <- unique(merged_data[, c("Superfamily", "Suborder")])
suborder_lookup <- setNames(lookup_df$Suborder, lookup_df$Superfamily)

superfamilies <- sort(unique(merged_data$Superfamily))
base_colors <- hue_pal()(length(superfamilies))
names(base_colors) <- superfamilies

# Manually define colors based on Superfamily and Suborder
# Assign colors based on the Suborder classification for each Superfamily
color_mapping <- sapply(names(base_colors), function(sf) {
  suborder <- suborder_lookup[[sf]]
  
  if (suborder == "Doridina") {
    lighten(base_colors[sf], amount = 0.5)
  } else if (suborder == "Cladobranchia") {
    darken(base_colors[sf], amount = 0.3)
  } else {
    "grey60"
  }
})
names(color_mapping) <- names(base_colors)

merged_data$Superfamily <- factor(merged_data$Superfamily, levels = names(color_mapping))
merged_data$Species <- factor(merged_data$Species, levels = unique(merged_data$Species))

# Set Superfamily factor order based on Superfamily's appearance in merged_data (for legend order)
sf_order <- unique(data2$Superfamily[order(data2$Species)])

# Reorder Superfamily factor to match the desired order
merged_data$Superfamily <- factor(merged_data$Superfamily, levels = sf_order)

# Create the ggplot object with the custom color palette
gg_plot <- ggplot(merged_data, aes(x = Mean_PDistance, y = PctOnTarget, color = Superfamily)) +
  geom_point(size = 4, alpha = 0.7) +
  geom_smooth(method = "lm", se = TRUE, color = "black", linetype = "dashed", size = 1) +  # Add linear model
  labs(title = paste("Kendall Rank Correlation between PctOnTarget and Mean_PDistance (τ =", 
                     round(kendall_correlation_result$estimate, 2), ")"),
       x = "Mean P-Distance", 
       y = "PctOnTarget") +
  theme_minimal() +
  scale_color_manual(values = color_mapping, guide = guide_legend(reverse = TRUE)) +
  theme(panel.background = element_rect(fill = "white"),
        plot.background = element_rect(fill = "white"),
        legend.background = element_rect(fill = "white"),
        axis.text = element_text(size = 10),
        axis.title = element_text(size = 12),
        plot.title = element_text(size = 14, face = "bold"),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 13))

# Save the plot with a white background
ggsave("ggplot_distancevsPctOnTarget_color.png", plot = gg_plot, width = 8.5, height = 6, dpi = 300)
ggsave("ggplot_distancevsPctOnTarget_color.pdf", plot = gg_plot, width = 8.5, height = 6, dpi = 300)

###I saw also differences in the reads, let's see if there is a correlation between number of exons and initial reads

#Read data
stats <- read.csv("D:/Sofi/Desktop/sofia/AMNH - postdoc/jessica's Lab/Postdoc project/4. target capture valdes et al/results/nucleotides_new_baitset/hybpiper_stats_new_baitset.csv", header = TRUE, sep = ",", row.names = 1)

stats$Species <- rownames(stats)
Number_reads <- stats[, c("Species", "NumReads")]  # extract Species + NumReads
merged_data <- merge(data, Number_reads, by = "Species")

write.csv(merged_data, "merged_data.csv", row.names = FALSE)
merged_data$sensitivity <- (merged_data$genesExtracted / 304) * 100 #304 is the number of genes in my baitset

#####excluding Layton et al 2020 species:

species_to_remove <- c("Ardeadoris_egretta", "Chromodoris_magnifica", "Chromodoris_westraliensis", "Doriprismatica_atromarginata", "Goniobranchus_fidelis")

merged_data <- merged_data[!merged_data$Species %in% species_to_remove, ]

####test normality

# Shapiro-Wilk test for normality
shapiro.test(merged_data$NumReads)

#no normal

# Perform Kendall rank correlation test
kendall_correlation_result <- cor.test(merged_data$NumReads, merged_data$sensitivity, method = "kendall")

# Print the Kendall correlation result
print(kendall_correlation_result)

# Set Superfamily factor order based on Superfamily's appearance in merged_data (for legend order)
sf_order <- unique(data2$Superfamily[order(data2$Species)])

# Reorder Superfamily factor to match the desired order
merged_data$Superfamily <- factor(merged_data$Superfamily, levels = sf_order)

# Create the ggplot object with a white background
gg_plot <- ggplot(merged_data, aes(x = NumReads, y = sensitivity, color = Superfamily)) +
  geom_point(size = 4, alpha = 0.5) +
  #geom_smooth(method = "loess", color = "black", se = FALSE) +
  labs(title = paste("Kendall Rank Correlation between %genesExtracted and NumberReads (τ =", round(kendall_correlation_result$estimate, 2), ")"),
       x = "Number of Reads", y = "% Genes Extracted") +
  theme_minimal() +  # Minimal theme with white background
  scale_color_manual(values = color_mapping, guide = guide_legend(reverse = TRUE)) +
  theme(panel.background = element_rect(fill = "white"),  # Ensure the plot background is white
        plot.background = element_rect(fill = "white"),    # Ensure the overall background is white
        legend.background = element_rect(fill = "white"),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 13))  # Optional: make legend background white

# Save the plot with a white background
ggsave("ggplot_readsvsexonsperct_colors.png", plot = gg_plot, width = 8, height = 6, dpi = 300)
ggsave("ggplot_readsvsexonsperct_colors.pdf", plot = gg_plot, width = 8, height = 6, dpi = 300)
