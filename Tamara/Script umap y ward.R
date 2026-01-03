
#ACTIVIDAD 3. ANÁLISIS DE UN CONJUNTO DE DATOS DE ORIGEN BIOLÓGICO MEDIANTE TÉCNICAS DE MACHINE LEARNING SUPERVISADAS Y NO SUPERVISADAS 


# Preparación del entorno

## Primero limpiamos el entorno
rm(list = ls())

## Y cargamos las librerías básicas
library(dplyr)
library(ggplot2)
library(cluster)
library(factoextra)
library(uwot)
library(pheatmap)

## Cargamos los datos
gene_expression <- read.csv("gene_expression.csv", header = FALSE, sep = ";")
column_names <- read.csv("column_names.txt", header = FALSE)
classes <- read.csv("classes.csv", header = FALSE, sep = ";", stringsAsFactors = FALSE)

## Asignamos los nombres de los genes
colnames(gene_expression) <- column_names$V1

## Añadimos la clase
data <- gene_expression
data$Class <- as.factor(classes[, 2])

## Comprobamos la estructura
str(data)


# Procesamiento de los datos

## Imputamos los NA por la media de cada gen
data_imputed <- data
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
  Class = data$Class
)

## Visualización

ggplot(umap_df, aes(x = UMAP1, y = UMAP2, color = Class)) +
  geom_point(size = 3) +
  theme_classic() +
  labs(title = "UMAP de expresión génica",
       x = "UMAP1", y = "UMAP2")


## Método no supervisado de clusterización: CLUSTERING JERÁRQUICO MÉTODO WARD 

# Calculamos la matriz de distancias (euclídea)
dist_matrix <- dist(data_scaled, method = "euclidean")

# Clustering jerárquico con Ward
hclust_ward <- hclust(dist_matrix, method = "ward.D2")

fviz_dend(
  hclust_ward,
  k = length(unique(data$Class)),     # número de clusters
  rect = TRUE,                        # colorea rectángulos de clusters
  cex = 0.5,                          # tamaño de texto de etiquetas
  main = "Clustering jerárquico - Método Ward",
  xlab = "Muestras",
  ylab = "Distancia"
)


