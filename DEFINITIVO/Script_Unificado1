rm(list=ls()) # Limpiar entorno

setwd("C:/Users/yanni/Desktop/MASTER/CUATRI 1/Algoritmos/Actividad grupal")

# Carga de datos
# Usamos check.names = FALSE para que los nombres de los genes se mantengan exactos
gene_exp <- read.csv("gene_expression.csv", sep = ";", header = FALSE)
genes <- read.table("column_names.txt", stringsAsFactors = FALSE)$V1
metadata <- read.csv("classes.csv", sep = ";", header = FALSE, col.names = c("SampleID", "Clase"))

# Asignación de nombres de genes a las columnas y IDs de muestras a las filas
colnames(gene_exp) <- genes
rownames(gene_exp) <- metadata$SampleID

# Se añade la columna de clase al final
df_final <- gene_exp
df_final$Clase <- metadata$Clase

# Verificamos el resultado
dim(df_final)      # Debería ser [801 filas x 501 columnas]
head(df_final[, 1:5]) # Ver las primeras 5 columnas (genes)

View(df_final)

# Carga de librerias
library(factoextra) # WARD
library(cluster) # DIANA
library(RDRToolbox) # LLE
library(Rtsne) # t-SNE
library(uwot) #UMAP
library(ggplot2)
library(dplyr)
library(stats)
library(ggdendro)
library(gridExtra)
library(pheatmap)

################################################################################ Clustering jerarquico divisivo: DIANA

# Separar datos numéricos
df_mat <- as.matrix(df_final[, -ncol(df_final)])  # Todas las columnas excepto la última (Clase)

# Implementación del clustering divisivo - euclidean
diana_euclidean <- diana(df_mat, metric = "euclidean", stand = FALSE)

# Visualización
colors <- rainbow(5)
clust_diana_euclidean <- fviz_dend(diana_euclidean, 
                                   cex = 0.5, # tamaño general del texto
                                   k = 5, # numero de grupos
                                   palette = colors, 
                                   lwd = 0.5, # hacer mas finas las lineas del dendrograma
                                   rect = 5,# pone cuadrados punteados al rededor de cada grupo (k)
                                   main = 'Euclidean',
                                   xlab = "Índice de Observaciones",
                                   ylab = "Distancia") + 
  theme_classic()

clust_diana_euclidean 

ggsave(
  filename = "Diana_Euc.png", # nombre del archivo
  plot = clust_diana_euclidean,# objeto ggplot
  width = 12, height = 5, dpi = 300
)

################################################################################ Clustering jerarquico aglomerativo: WARD

# Separar datos numéricos
df_mat <- as.matrix(df_final[, -ncol(df_final)])  # Todas las columnas excepto la última (Clase)

# Imputamos los NA por la media de cada gen
for (i in 1:ncol(df_mat)) {
  df_mat[is.na(df_mat[, i]), i] <- mean(df_mat[, i], na.rm = TRUE)
}

# Escalamos los datos
df_scaled <- scale(df_mat)

# Eliminamos columnas que hayan quedado con NA tras el escalado
df_scaled <- df_scaled[, colSums(is.na(df_scaled)) == 0]

# Calculamos la matriz de distancias (euclídea)
dist_matrix <- dist(df_scaled, method = "euclidean")

# Clustering jerárquico con Ward
hclust_ward <- hclust(dist_matrix, method = "ward.D2")

# Visualización
colors <- rainbow(5)
clust_ward <- fviz_dend(hclust_ward,
                        cex = 0.5, # tamaño general del texto
                        k = 5, # numero de grupos
                        palette = colors, 
                        lwd = 0.5, # hacer mas finas las lineas del dendrograma
                        rect = 5, # pone cuadrados punteados al rededor de cada grupo (k)
                        main = 'Ward',
                        xlab = "Índice de Observaciones",
                        ylab = "Distancia") + 
  theme_classic()

clust_ward

ggsave(
  filename = "Ward.png", # nombre del archivo
  plot = clust_ward, # objeto ggplot
  width = 12, height = 5, dpi = 300
)

################################################################################ t-SNE


# 1. Preparamos los datos: eliminamos la columna 'Clase' para el cálculo
# y nos aseguramos de que no haya filas duplicadas (t-SNE falla si las hay)
data_matrix <- as.matrix(df_final[, -ncol(df_final)]) 
labels <- df_final$Clase

# Ejecutamos t-SNE
# perplexity: suele ser entre 5 y 50 (depende del tamaño de tu muestra, 801 filas)
set.seed(42) # Para que el resultado sea reproducible
tsne_out <- Rtsne(X = data_matrix, dims = 2, perplexity = 30, verbose = TRUE, check_duplicates = FALSE)

# Creamos el dataframe para graficar
tsne_result <- data.frame(X1 = tsne_out$Y[,1], X2 = tsne_out$Y[,2], Grupo = labels)

# Graficamos con ggplot2
TSNE_plot <- ggplot(tsne_result, aes(x = X1, y = X2, color = Grupo)) +
  geom_point(size = 2, alpha = 0.7) +
  scale_color_brewer(palette = "Set1") + # Paleta automática para las clases de cáncer
  labs(title = "t-SNE - Tipos de Cáncer", 
       x = "Dimensión 1", y = "Dimensión 2", color = "Clase") +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5))

TSNE_plot

ggsave(
  filename = "TSNE.png", # nombre del archivo
  plot = TSNE_plot, # objeto ggplot
  width = 8, height = 6, dpi = 300
)


################################################################################ UMAP

# Procesamiento de los datos
## Imputamos los NA por la media de cada gen
data_imputed <- df_final
for (i in 1:(ncol(data_imputed)-1)) {
  data_imputed[is.na(data_imputed[, i]), i] <- 
    mean(data_imputed[, i], na.rm = TRUE)
}

## Escalamos los datos
data_scaled <- scale(data_imputed[, -ncol(data_imputed)])

## Eliminamos columnas que hayan quedado con NA tras el escalado
data_scaled <- data_scaled[, colSums(is.na(data_scaled)) == 0]

# Método no supervisado de reducción de la dimensionalidad: UMAP 
set.seed(123)
umap_results <- umap(
  data_scaled,
  n_neighbors = round(0.2 * nrow(data_scaled)),
  n_components = 2,
  min_dist = 0.1,
  local_connectivity = 1,
  ret_model = TRUE
)

umap_df <- data.frame(
  UMAP1 = umap_results$embedding[,1],
  UMAP2 = umap_results$embedding[,2],
  Clase = df_final$Clase
)

## Visualización
UMAP_plot <- ggplot(umap_df, aes(x = UMAP1, y = UMAP2, color = Clase)) +
  geom_point(size = 3) +
  theme_classic() +
  labs(title = "UMAP de expresión génica",
       x = "UMAP1", y = "UMAP2")

UMAP_plot

ggsave(
  filename = "UMAP.png", # nombre del archivo
  plot = UMAP_plot, # objeto ggplot
  width = 8, height = 6, dpi = 300
)
