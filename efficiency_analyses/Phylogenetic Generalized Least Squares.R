#Phylogenetic Generalized Least Squares (PGLS)

library(ape)
library(geiger)
library(nlme)
library(phytools)
library(segmented)
library(ggplot2)

setwd("D:/Sofi/Desktop/sofia/AMNH - postdoc/jessica's Lab/Postdoc project/4. target capture valdes et al/Mol-Phylo-Evol/review1/PGLS")

data <- read.csv("D:/Sofi/Desktop/sofia/AMNH - postdoc/jessica's Lab/Postdoc project/4. target capture valdes et al/results/nucleotides_new_baitset/exons2_R1.csv")


#reading distances calculated
uncorr_distance <- read.csv("D:/Sofi/Desktop/sofia/AMNH - postdoc/jessica's Lab/Postdoc project/4. target capture valdes et al/results/nucleotides_new_baitset/p_distance_matrix-nogap_R1.csv", header = TRUE, sep = ",", row.names = 1)
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

write.csv(merged_data, "PGLS_data_R1.csv", row.names = FALSE)

merged_data$sensitivity <- (merged_data$genesExtracted / 304) * 100 #304 is the number of genes in my baitset

# Set row names, To have the species names in the rows
rownames(merged_data) <- merged_data$Species

Tree <- read.tree("D:/Sofi/Desktop/sofia/AMNH - postdoc/jessica's Lab/Postdoc project/4. target capture valdes et al/Mol-Phylo-Evol/review1/PGLS/alignment_75_R1.contree")

plot(Tree)

#check that the names match between the tree and the data frame
name.check(Tree, merged_data)

#Che(ck if there is a correlation between Mean_PDistance and sensitivity
plot(merged_data[, c("Mean_PDistance","sensitivity")])

##Because my correlation is segmented

# Estimate breakpoint using non-phylogenetic model
lm_fit <- lm(sensitivity ~ Mean_PDistance, data = merged_data)

seg_fit <- segmented(
  lm_fit,
  seg.Z = ~ Mean_PDistance
)

# Extract breakpoint
bp <- seg_fit$psi[2]

# Create hinge (piecewise) term
merged_data$hinge <- pmax(0, merged_data$Mean_PDistance - bp)

###PGLS analyses
pglsModel <- gls(sensitivity ~ Mean_PDistance + hinge, correlation = corBrownian(phy = Tree),
                 data = merged_data, method = "ML")
summary(pglsModel)

# Slopes before and after breakpoint
slope_before <- coef(pglsModel)["Mean_PDistance"]
slope_after  <- slope_before + coef(pglsModel)["hinge"]

slope_before
slope_after

#Add fitted values for plotting
merged_data$fit_piecewise <- predict(pglsModel)

# Plot

p <- ggplot(merged_data, aes(x = Mean_PDistance, y = sensitivity)) +
  geom_point(alpha = 0.5) +
  geom_line(aes(y = fit_piecewise), color = "black", linewidth = 1) +
  geom_vline(xintercept = bp, linetype = "dashed", color = "red") +
  labs(
    x = "Mean P-distance",
    y = "% Genes Extracted",
    title = "Piecewise PGLS: Sensitivity vs Genetic Distance"
  ) +
  theme_minimal()

ggsave(filename = "piecewise_PGLS_sensitivity.png", plot = p, width = 8, height = 6, dpi = 300)


###SPECIFICITY

specif <- read.csv("D:/Sofi/Desktop/sofia/AMNH - postdoc/jessica's Lab/Postdoc project/4. target capture valdes et al/results/nucleotides_new_baitset/hybpiper_stats_new_baitset.csv", header = TRUE, sep = ",")

merged_data$PctOnTarget <- specif$PctOnTarget[match(merged_data$Species, specif$Species)]

###Layton 2020 are exceptions with so many reads, so exclude outliers in my analyses
species_to_remove <- c("Ardeadoris_egretta", "Chromodoris_magnifica", "Chromodoris_westraliensis", "Doriprismatica_atromarginata", "Goniobranchus_fidelis")

merged_data2 <- merged_data[!merged_data$Species %in% species_to_remove, ]

name.check(Tree, merged_data2)


pglsspecif <- gls(PctOnTarget ~ Mean_PDistance, correlation = corBrownian(phy = Tree),
                 data = merged_data2, method = "ML")
summary(pglsspecif)

coef(pglsspecif)

merged_data2$fit_pgls <- predict(pglsspecif)

p <- ggplot(merged_data2, aes(x = Mean_PDistance, y = PctOnTarget)) +
  geom_point(alpha = 0.6, size = 3) +
  geom_line(aes(y = fit_pgls), color = "black", linewidth = 1) +
  labs(
    x = "Mean P-distance",
    y = "Specifiticy",
    title = "PGLS: Specificity vs P-Distance"
  ) +
  theme_minimal()

print(p)

ggsave(filename = "PGLS_specificity.png", plot = p, width = 8, height = 6, dpi = 300)
                                                                              
