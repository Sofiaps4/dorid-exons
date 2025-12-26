# Load required packages
library(tidyverse)

# Cargar el archivo sin nombres de columna
raw <- read_tsv("D:/Sofi/Desktop/sofia/AMNH - postdoc/jessica's Lab/Postdoc project/4. target capture valdes et al/results/nucleotides_new_baitset/seq_lengths.tsv", col_names = FALSE)

# Extraer los nombres de los genes (desde la columna 3 en adelante de la primera fila)
gene_names <- as.character(unlist(raw[1, -c(1,2)]))

# Extraer longitudes promedio (desde columna 3 en adelante de la segunda fila)
mean_lengths <- as.numeric(unlist(raw[2, -c(1,2)]))

# Crear un named vector para luego emparejar por nombre de gen
mean_length_df <- tibble(Gene = gene_names, MeanLength = mean_lengths)

# Extraer datos reales (desde la fila 3 en adelante)
data <- raw[-c(1,2), ]
colnames(data) <- c("Suborder", "Species", gene_names)

# Convertir columnas de genes a numérico
data[gene_names] <- lapply(data[gene_names], as.numeric)

# Convertir a formato largo
long_data <- pivot_longer(data, cols = all_of(gene_names),
                          names_to = "Gene", values_to = "GeneLength")

# Agregar la columna de longitud promedio a cada gen
long_data <- left_join(long_data, mean_length_df, by = "Gene")

# Calcular el porcentaje de longitud
long_data <- long_data %>%
  mutate(PercentLength = (GeneLength / MeanLength) * 100)

# Bin percentages
long_data <- long_data %>%
  mutate(Bin = cut(PercentLength,
                   breaks = c(0, seq(10, 100, by = 10), Inf),
                   labels = c("0-10%", "10-20%", "20-30%", "30-40%", "40-50%",
                              "50-60%", "60-70%", "70-80%", "80-90%", "90-100%", ">100%"),
                   include.lowest = TRUE, right = FALSE))

# Count number of genes per bin per suborder
bin_counts <- long_data %>%
  group_by(Suborder, Bin) %>%
  summarise(GeneCount = n(), .groups = "drop")

# Calculate total genes per Suborder to compute %
total_per_suborder <- long_data %>%
  group_by(Suborder) %>%
  summarise(TotalGenes = n(), .groups = "drop")

# Merge and calculate percent
plot_data <- left_join(bin_counts, total_per_suborder, by = "Suborder") %>%
  mutate(PercentGenes = (GeneCount / TotalGenes) * 100)

# Before plotting, set factor levels to include all Suborders
plot_data$Suborder <- factor(plot_data$Suborder, levels = c("Doridina", "Cladobranchia", "Outgroup"))

# Plot
gg_plot = ggplot(plot_data, aes(x = Bin, y = PercentGenes, fill = Suborder)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(x = "Gene Length (% of Bait set)", y = "Percentage of Genes",
       title = "Gene length recovery by Suborder") +
  theme_minimal() +
  scale_fill_manual(values = c("Doridina" = "#F17CB0",
                               "Cladobranchia" = "#984EA3",
                               "Outgroup" = "grey85")) +
  theme(panel.background = element_rect(fill = "white"),  # Ensure the plot background is white
        plot.background = element_rect(fill = "white"),    # Ensure the overall background is white
        legend.background = element_rect(fill = "white"),
        legend.text = element_text(size = 12),
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.title = element_blank())


ggsave("ggplot_percentageLength_R1.png", plot = gg_plot, width = 8, height = 6, dpi = 300)
ggsave("ggplot_percentageLength_R1.pdf", plot = gg_plot, width = 8, height = 6, dpi = 300)
ggsave("ggplot_percentageLength_R1.svg", plot = gg_plot, width = 8, height = 6, dpi = 300)
