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
library(RDRToolbox) # LLE
library(ggplot2)
library(dplyr)
library(factoextra)
library(stats)
library(ggdendro)
library(cluster) # DIANA
library(gridExtra)
library(pheatmap)

################################################################################ LLE
# Algoritmo
# Funcion LLE()
#   x: matriz sobre la cual se va a reducir la dimensionalidad
#   dim: numero de dimensiones de salida
#   k: numero de vecinos cercanos. Puede aproximarse con la funcion calc_k()

#   la salida será un dataframe con dimensión = dim


# Separar datos numéricos y clase
df_mat <- as.matrix(df_final[, -ncol(df_final)])  # Todas las columnas excepto la última (Clase), que es character
clase <- df_final$Clase

# Verificar dimensiones
dim(df_mat)  # Debería ser [801 x 500]

# Ahora aplicar LLE. Se podria construir un optimizador del k (entre 5%-15% de la muestra)
lle.results <- LLE(df_mat, dim = 2, k = 50)

# Crear dataframe para visualización
lle_df <- data.frame(
  LLE1 = lle.results[, 1],
  LLE2 = lle.results[, 2],
  Clase = clase
)

# Graficamos
LLE_Yannis <- ggplot(lle_df, aes(x = LLE1, y = LLE2, color = df_final$Clase)) +
                geom_point(size = 3) +
                scale_color_manual(values = c("red", "blue", "green", "orange", "purple")) +
                labs(title = "Método LLE Types of Cancer", x = "X1", y = "X2", color = "Grupo") +
                theme_classic() +
                theme(panel.grid.major = element_line(color = "gray90"), panel.grid.minor = element_blank(),
                      panel.background = element_rect(fill = "gray95"), plot.title=element_text(hjust=0.5))

ggsave(
  filename = "LLE_Yannis.png", # nombre del archivo. ES IMPORTANTE PONER LA EXTENSION .PNG...
  plot = LLE_Yannis,# objeto ggplot
  width = 8, height = 6, dpi = 300
)



################################################################################ Clustering jerarquico divisivo: DIANA

# Implementación del clustering divisivo
df_scaled <- scale(df_mat)


# ideal para datos donde las distancias más pequeñas
diana_euclidean <- diana(df_mat, metric = "euclidean", stand = FALSE) #uso DIANA, meto metrica euclidean (aunque esta ya es default) y stand = FALSE porque ya hemos estandarizado antes
# ideal para datos con diferentes escalas o datos categóricos
diana_manhattan <- diana(df_mat, metric = "manhattan", stand = FALSE) #esta vez lo uso con distancia manhattan

#euclidean
colors <- rainbow(5)
clust_diana_euclidean <- fviz_dend(diana_euclidean, 
                                   cex = 0.5, # tamño general del texto
                                   k = 5, # numero de grupos
                                   palette = colors, 
                                   lwd = 0.5, # hacer mas finas las lineas del dendrograma
                                   rect = 5,# pone cuadrados punteados al rededor de cada grupo (k)
                                   main = 'Euclidean',
                                   xlab = "Índice de Observaciones",
                                   ylab = "Distancia") + 
  theme_classic()

clust_diana_euclidean 

ggsave(
  filename = "Diana_Euc.png", # nombre del archivo
  plot = clust_diana_euclidean ,# objeto ggplot
  width = 8, height = 6, dpi = 300
)

#manhattan
colors <- rainbow(5)
clust_diana_manhattan <- fviz_dend(diana_manhattan, 
                                   cex = 0.3, 
                                   k = 5,
                                   palette = colors, 
                                   lwd = 0.5,
                                   rect = 5,
                                   main = 'Manhattan',
                                   xlab = "Índice de Observaciones",
                                   ylab = "Distancia") + 
  theme_classic()

clust_diana_manhattan


ggsave(
  filename = "Diana_Man.png", # nombre del archivo
  plot = clust_diana_manhattan ,# objeto ggplot
  width = 8, height = 6, dpi = 300
)

MAn_Euc_Diana <- grid.arrange(clust_diana_euclidean, clust_diana_manhattan, nrow = 2) #PARA REPRESENTAR TODOS LOS DENDROGRAMAS JUNTOS
#vemos que la distancia de manhattan es mejor en este caso

ggsave(
  filename = "Man_Euc_Diana.png", # nombre del archivo
  plot = MAn_Euc_Diana ,# objeto ggplot
  width = 8, height = 6, dpi = 300
)

################################################################################ LDA

# Teoría:
# Para contruir el modelo, hay que relacionar varibles como el radius, compactness... con B (benigno) o M (metastasico), es decir,
# se va a intentar entrenar un modelo que pueda predecir la malignidad de un tumor a partir de variables observadas.
# Los predicores (X) van a ser las medidas de los tumores, y el vector respuesta Y (la malignidad)
# X = MEDIDAS Y = DIAGNOSIS para el training testing 



library(MASS)      # Para LDA
library(caret)     # Para partición de datos y matriz de confusión


## 1. PREPARACIÓN DE DATOS


# Convertir Clase a factor
df_final$Clase <- as.factor(df_final$Clase)

# Partición de datos (70% entrenamiento, 30% prueba)
set.seed(42)
train_index <- createDataPartition(df_final$Clase, p = 0.7, list = FALSE)

training_genes <- df_final[train_index, ]
testing_genes <- df_final[-train_index, ]

# Verificar distribución de clases
print(table(training_genes$Clase))
print(table(testing_genes$Clase))


## 2. TRATAMIENTO DE VARIABLES PROBLEMÁTICAS


# Identificar variables con varianza cero (constantes)
# Estas variables no aportan información discriminante
var_no_cero <- apply(training_genes[, -ncol(training_genes)], 2, var) != 0

# Filtrar solo variables con varianza
train_filtered <- training_genes[, var_no_cero]
test_filtered <- testing_genes[, var_no_cero]

# Limpiar nombres de columnas (eliminar guiones y caracteres problemáticos)
colnames(train_filtered) <- make.names(colnames(train_filtered))
colnames(test_filtered) <- make.names(colnames(test_filtered))

# Añadir la clase de nuevo
train_filtered$Clase <- training_genes$Clase
test_filtered$Clase <- testing_genes$Clase


## 3. ESCALADO DE VARIABLES


# Extraer solo columnas numéricas
numerical_cols_train <- train_filtered[, sapply(train_filtered, is.numeric)]
numerical_cols_test <- test_filtered[, sapply(test_filtered, is.numeric)]

# Escalar (centrar y reducir)
# IMPORTANTE: usar los parámetros del training para escalar el test
train_means <- colMeans(numerical_cols_train)
train_sds <- apply(numerical_cols_train, 2, sd)

scaled_train <- scale(numerical_cols_train, center = train_means, scale = train_sds)
scaled_test <- scale(numerical_cols_test, center = train_means, scale = train_sds)

# Reconstruir datasets
train_lda <- as.data.frame(scaled_train)
train_lda$Clase <- train_filtered$Clase

test_lda <- as.data.frame(scaled_test)
test_lda$Clase <- test_filtered$Clase


## 4. CONSTRUCCIÓN DE LA FÓRMULA


genes <- colnames(train_lda)[1:(ncol(train_lda)-1)]
formula_lda <- as.formula(paste("Clase ~", paste(genes, collapse = "+")))


## 5. MODELO LDA 

set.seed(42) # Semilla 42 en todo el grupo para reproducibilidad
lda_model <- lda(formula_lda, data = train_lda)


## 6. PREDICCIONES Y EVALUACIÓN


# Predicciones en conjunto de prueba
lda_predictions <- predict(lda_model, newdata = test_lda)

# Clases predichas vs reales
predicted_lda <- lda_predictions$class
true_classes <- test_lda$Clase

# Matriz de confusión
conf_lda <- confusionMatrix(predicted_lda, true_classes)
print(conf_lda)

# Guardar métricas
accuracy_lda <- conf_lda$overall['Accuracy']
kappa_lda <- conf_lda$overall['Kappa']


## 7. VISUALIZACIÓN LDA


# Proyección LDA del conjunto de entrenamiento
lda_proj <- predict(lda_model)$x
lda_data <- cbind(train_lda, lda_proj)

# Determinar número de funciones discriminantes
num_lds <- ncol(lda_proj)

# Gráfico LDA (bidimensional si hay LD1 y LD2)
if(num_lds >= 2) {
  
  plot_lda <- ggplot(lda_data, aes(x = LD1, y = LD2, color = Clase)) +
    geom_point(size = 2.5, alpha = 0.7) +
    stat_ellipse(level = 0.95, linetype = 2) + # Elipses de confianza
    scale_color_manual(values = c("red", "blue", "green", "orange", "purple")) +
    labs(title = "LDA - Separación de Tipos de Cáncer",
         subtitle = sprintf("Accuracy: %.2f%% | Kappa: %.3f", 
                            accuracy_lda*100, kappa_lda),
         x = "LD1 (Primera función discriminante)",
         y = "LD2 (Segunda función discriminante)",
         color = "Tipo de Cáncer") +
    theme_classic() +
    theme(
      panel.grid.major = element_line(color = "gray90"),
      panel.background = element_rect(fill = "gray95"),
      plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
      plot.subtitle = element_text(hjust = 0.5, size = 10)
    )
  
  print(plot_lda)
  
  ggsave(
    filename = "LDA_Cancer_Types.png",
    plot = plot_lda,
    width = 10, height = 7, dpi = 300
  )
  
} else {
  
  # Gráfico unidimensional si solo hay LD1
  plot_lda <- ggplot(lda_data, aes(x = LD1, y = 0, color = Clase)) +
    geom_jitter(height = 0.05, size = 2.5, alpha = 0.7) +
    scale_color_manual(values = c("red", "blue", "green", "orange", "purple")) +
    labs(title = "LDA - Separación de Tipos de Cáncer (1D)",
         subtitle = sprintf("Accuracy: %.2f%% | Kappa: %.3f", 
                            accuracy_lda*100, kappa_lda),
         x = "LD1 (Primera función discriminante)",
         y = "",
         color = "Tipo de Cáncer") +
    theme_classic() +
    theme(
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
      plot.subtitle = element_text(hjust = 0.5, size = 10)
    )
  
  print(plot_lda)
  
  ggsave(
    filename = "LDA_Cancer_Types_1D.png",
    plot = plot_lda,
    width = 10, height = 4, dpi = 300
  )
}


## 8. ANÁLISIS DE GENES MÁS IMPORTANTES


# Top 10 genes en LD1
ld1_weights <- abs(lda_model$scaling[, 1])
top_ld1 <- sort(ld1_weights, decreasing = TRUE)[1:10]

print(top_ld1)

if(num_lds >= 2) {
  # Top 10 genes en LD2
  ld2_weights <- abs(lda_model$scaling[, 2])
  top_ld2 <- sort(ld2_weights, decreasing = TRUE)[1:10]
  
  print(top_ld2)
}


## 9. CURVA DE SEPARACIÓN ENTRE GRUPOS


# Gráfico de densidad para cada función discriminante por clase
if(num_lds >= 1) {
  
  plot_density_ld1 <- ggplot(lda_data, aes(x = LD1, fill = Clase)) +
    geom_density(alpha = 0.5) +
    scale_fill_manual(values = c("red", "blue", "green", "orange", "purple")) +
    labs(title = "Distribución de LD1 por Tipo de Cáncer",
         x = "LD1",
         y = "Densidad",
         fill = "Tipo de Cáncer") +
    theme_classic() +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold")
    )
  
  print(plot_density_ld1)
  
  ggsave(
    filename = "LDA_LD1_Density.png",
    plot = plot_density_ld1,
    width = 10, height = 6, dpi = 300
  )
}

