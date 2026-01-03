rm(list=ls()) # Limpiar entorno

setwd("C:/Users/ritap/Downloads")

# Carga de datos
library(dplyr)


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


# Imputación de los datos NA, Verificar NAs antes de escalar
if(anyNA(gene_exp)) { 
  print("Hay valores NA. Considera imputar con la media o 0.")
}

# Quitar genes constantes (Varianza cero)
df_num <- gene_exp[, apply(gene_exp, 2, function(x) sd(x, na.rm = TRUE) != 0)]

# 3. Escalar los datos (Z-score)
df_scaled <- data.frame(scale(df_num))

# Verificar valores finitos (Inf, -Inf, NaN)
if(any(!is.finite(as.matrix(df_scaled)))) {
  print("Aún existen valores no finitos.")
} else {
  print("Datos escalados y listos.")
}



######################
# LE # Rita #
######################

# Librerías
library(Rdimtools)
library(ggplot2)

# Aseguramos que los datos sean numéricos
X <- as.matrix(gene_exp)
X <- apply(X, 2, as.numeric)

# Laplacian Eigenmaps con proporción
le.results <- do.lapeig(
  X = X,
  Ndim = 3,
  type = c("proportion", 0.20),  # 20% de vecinos
  weighted = FALSE
)

# Laplacian Eigenmaps con k numérico
#le.results <- do.lapeig(
#  X = X,
#  Ndim = 3,
#  type = c("knn", 150),  # 20 puntos vecinos
#  weighted = FALSE
#)

# Resultados en data frame
le.df <- data.frame(le.results$Y)
le.df$Clase <- metadata$Clase

# Gráfico
ggplot(le.df, aes(x = X1, y = X2, color = Clase)) +
  geom_point(size = 2, alpha = 0.8) +
  labs(
    title = "Laplacian Eigenmaps (k=20%) - Expresión Génica",
    x = "Dimensión 1",
    y = "Dimensión 2",
    color = "Clase"
  ) +
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5)
  )

grafico_le <- ggplot(le.df, aes(x = X1, y = X2, color = Clase)) +
  geom_point(size = 1, alpha = 0.8) +
  labs(
    title = "Laplacian Eigenmaps (k=20%) - Expresión Génica",
    x = "Dimensión 1",
    y = "Dimensión 2",
    color = "Clase"
  ) +
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5)
  )

ggsave("Laplacian_Eigenmaps.png", plot = grafico_le, width = 8, height = 6, dpi = 300)

###################################
# Clustering Single Linkage # Rita #
###################################

# Librerías
library(factoextra)
library(pheatmap)

# Preparación de datos
# Usamos la matriz escalada que ya tenemos
X_clust <- as.matrix(df_scaled)

# Cálculo de Distancia y Modelo Single
dist_matrix <- dist(X_clust, method = "manhattan") # Probado con manhattan y euclídea
hclust_single <- hclust(dist_matrix, method = "single")

# Visualización del Dendrograma
# 'k = 5' marca los 5 grupos principales con colores
clust_plot <- fviz_dend(hclust_single, 
                        k = 5, 
                        cex = 0.6, 
                        palette = "jco", 
                        rect = TRUE, 
                        main = "Clustering Jerárquico: Método Single",
                        xlab = "Muestras", 
                        ylab = "Distancia (Proximidad)") + 
  theme_classic()

print(clust_plot)

ggsave("Dendrograma_Single.png", plot = clust_plot, width = 10, height = 7)

#######################################################################
# Clustering Jerárquico (Single, Average, Complete y Ward.2D) # Rita #
#######################################################################

# Librerías necesarias
library(dplyr)
library(factoextra)
library(stats)
library(cluster)
library(ggplot2) 
library(gridExtra)
library(pheatmap)

# Preparación de datos (usando lo anterior)
# Usamos df_scaled que ya tiene las columnas constantes eliminadas y valores normalizados
X_clust <- as.matrix(df_scaled)

# Mapa de Calor (Heatmap)
# Añadimos anotación de clase para comparar visualmente
annotation_col <- data.frame(Clase = metadata$Clase)
rownames(annotation_col) <- rownames(X_clust)


# Cálculo de la Matriz de Distancia
# En expresión génica, la distancia euclídea es el estándar
dist_matrix <- dist(X_clust)

# Modelos de Clustering Jerárquico
hclust_model_single   <- hclust(dist_matrix, method = "single") 
hclust_model_complete <- hclust(dist_matrix, method = "complete") 
hclust_model_average  <- hclust(dist_matrix, method = "average") 
hclust_model_ward     <- hclust(dist_matrix, method = "ward.D2") # ward.D2 es más preciso

# Visualización de Dendrogramas
colors <- rainbow(5)

# Función rápida para generar los gráficos de dendrogramas
crear_dend <- function(modelo, titulo) {
  fviz_dend(modelo, k = 5, cex = 0.5, palette = colors, 
            main = titulo, xlab = "Muestras", ylab = "Distancia",
            rect = TRUE) + theme_classic()
}

clust_single   <- crear_dend(hclust_model_single, "Single (Encadenamiento)")
clust_complete <- crear_dend(hclust_model_complete, "Complete (Máxima)")
clust_average  <- crear_dend(hclust_model_average, "Average (Promedio)")
clust_ward     <- crear_dend(hclust_model_ward, "Ward.D2 (Mínima Varianza)")

# Organizar en una cuadrícula
grid.arrange(clust_single, clust_complete, clust_average, clust_ward, nrow = 2)

# Asignación de Clusters al DataFrame Original
# Usamos df_final para mantener la Clase original junto a los clusters
df_final$cluster_single   <- as.factor(cutree(hclust_model_single, k = 5))
df_final$cluster_complete <- as.factor(cutree(hclust_model_complete, k = 5))
df_final$cluster_average  <- as.factor(cutree(hclust_model_average, k = 5))
df_final$cluster_ward     <- as.factor(cutree(hclust_model_ward, k = 5))

# Verificamos la concordancia entre Clase real y Cluster de Ward (ejemplo)
table(Real = df_final$Clase, Cluster_Ward = df_final$cluster_ward)



###########################################
# Heatmaps: Complete vs Ward # Rita #
###########################################

# Preparar la anotación de las columnas (las muestras)
# Esto nos permite ver si los clusters coinciden con la "Clase" real
annotation_col <- data.frame(Clase = metadata$Clase)
rownames(annotation_col) <- rownames(df_scaled)

# Mapa de Calor con método COMPLETE
pheatmap(t(as.matrix(df_scaled)), 
         clustering_method = "complete", 
         annotation_col = annotation_col,
         show_colnames = FALSE,
         main = "Heatmap: Clustering Jerárquico (Complete)",
         fontsize_row = 6,
         color = colorRampPalette(c("blue", "white", "red"))(50))

# Mapa de Calor con método WARD.D2
pheatmap(t(as.matrix(df_scaled)), 
         clustering_method = "ward.D2", 
         annotation_col = annotation_col,
         show_colnames = FALSE,
         main = "Heatmap: Clustering Jerárquico (Ward.D2)",
         fontsize_row = 6,
         color = colorRampPalette(c("blue", "white", "red"))(50))

