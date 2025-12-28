rm(list=ls()) # Limpiar entorno


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


###################
######T-SNE########
###################



library(ggplot2)
library(Rtsne)

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
ggplot(tsne_result, aes(x = X1, y = X2, color = Grupo)) +
  geom_point(size = 2, alpha = 0.7) +
  scale_color_brewer(palette = "Set1") + # Paleta automática para las clases de cáncer
  labs(title = "t-SNE - Tipos de Cáncer", 
       x = "Dimensión 1", y = "Dimensión 2", color = "Clase") +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5))



#######################
####cluster-average####
#######################




library(factoextra)

# 1. Calculamos la matriz de distancias (Euclídea por defecto)
dist_matrix <- dist(data_matrix)

# 2. Aplicamos el clustering con el método 'average' (UPGMA)
hclust_model_average <- hclust(dist_matrix, method = "average")

# 3. Visualizamos el dendrograma
# 'k = 5' asume que tienes 5 tipos de cáncer según tu código de soporte
clust_average <- fviz_dend(hclust_model_average, 
                           k = 5, 
                           cex = 0.5, 
                           palette = "jco", 
                           rect = TRUE, 
                           rect_fill = TRUE,
                           main = "Clustering Jerárquico - Método Average",
                           xlab = "Muestras", 
                           ylab = "Distancia (Disimilitud)") + 
  theme_classic()

print(clust_average)
