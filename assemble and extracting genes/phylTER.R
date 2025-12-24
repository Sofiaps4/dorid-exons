setwd("D:/Sofi/Desktop/sofia/AMNH - postdoc/jessica's Lab/Postdoc project/4. target capture valdes et al/Mol-Phylo-Evol/review1/")


install.packages("phylter")
library("phylter")

tree_files <- list.files(
  "D:/Sofi/Desktop/sofia/AMNH - postdoc/jessica's Lab/Postdoc project/4. target capture valdes et al/Mol-Phylo-Evol/review1/phylter_input",
  pattern = "\\.treefile$",
  full.names = TRUE
)

trees <- lapply(tree_files, ape::read.tree)
class(trees) <- "multiPhylo"

names <- tools::file_path_sans_ext(basename(tree_files))

results <- phylter(trees, gene.names = names)

#To get the list of outliers detected by phylter, simply type:
results$Final$Outliers

# Assuming results is your Phylter output
outliers <- results$Final$Outliers

# save results in a table
write.table(outliers, "phylter_outliers_table.txt",
            row.names = FALSE, col.names = TRUE, quote = FALSE)


# Get a summary: nb of outliers, gain in concordance, etc.
summary(results)

# Show the number of species in each gene, and how many per gene are outliers
plot(results, "genes") 

png("plot_genes.png", width = 2400, height = 2000, res = 300)  # larger image for better resolution

# Set plotting parameters to increase text size
par(cex = 4,       # overall character expansion
    cex.axis = 1.8, # axis labels
    cex.lab = 15,    # axis titles
    cex.main = 20) # main titles

# Create the plot
plot(results, "genes")

# Close device
dev.off()


# Show the number of genes where each species is found, and how many are outliers
plot(results, "species") 

png("plot_species2.png", width = 2400, height = 2000, res = 300)

plot(results, "species",
     cex.main = 30,    # main title
     cex.lab = 30,     # axis titles
     cex.axis = 30,    # axis numbers
     cexNames = 25)    # Y-axis species names

dev.off()

# Compare before and after genes x species matrices, highlighting missing data and outliers 
# identified (not efficient for large datasets)
plot2WR(results) 

png("plot2WR.png", width = 2400, height = 2000, res = 300)
plot2WR(results)
dev.off()

# Plot the dispersion of data before and after outlier removal. One dot represents one 
# gene x species association
plotDispersion(results) 

png("plot_dispersion.png", width = 2400, height = 2000, res = 300)
plotDispersion(results)
dev.off()

# Plot the genes x genes matrix showing pairwise correlation between genes
plotRV(results) 

png("plot_RV.png", width = 2400, height = 2000, res = 300)
plotRV(results)
dev.off()

# Plot optimization scores during optimization
plotopti(results)

png("plot_optimization.png", width = 2400, height = 2000, res = 300)
plotopti(results)
dev.off()
