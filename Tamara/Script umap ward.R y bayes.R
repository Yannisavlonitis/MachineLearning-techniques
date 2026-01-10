
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
)


umap_df <- data.frame(
  UMAP1 = umap_results[, 1],
  UMAP2 = umap_results[, 2],
  Class = data$Class
)


## Visualización

ggplot(umap_df, aes(x = UMAP1, y = UMAP2, color = Class)) +
  geom_point(size = 3) +
  theme_classic() +
  labs(title = "UMAP de expresión génica",
       x = "UMAP1", y = "UMAP2")


# Método no supervisado de clusterización: CLUSTERING JERÁRQUICO MÉTODO WARD 

## Calculamos la matriz de distancias (euclídea)
dist_matrix <- dist(data_scaled, method = "euclidean")

## Clustering jerárquico con Ward
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




# Método de apredizaje supervisado: NAIVE BAYES


## Librerías necesarias
library(caret)


## Separación entrenamiento/prueba
set.seed(123)

trainIndex <- createDataPartition(
  data_imputed$Class,
  p = 0.8,          # 80% entrenamiento, 20% prueba
  list = FALSE
)

trainData <- data_imputed[trainIndex, ]
testData  <- data_imputed[-trainIndex, ]

trainData$Class <- as.factor(trainData$Class)
testData$Class  <- as.factor(testData$Class)


## Identificar variables con varianza cero en TRAIN
var_no_cero <- apply(trainData[, -ncol(trainData)], 2, var) != 0

## Filtrar variables
train_filtered <- trainData[, c(var_no_cero, TRUE)]
test_filtered  <- testData[,  c(var_no_cero, TRUE)]

## Limpiar nombres de columnas
colnames(train_filtered) <- make.names(colnames(train_filtered))
colnames(test_filtered)  <- make.names(colnames(test_filtered))

## Entrenamiento Naive Bayes
nb_model <- train(
  Class ~ .,
  data = train_filtered,
  method = "nb",
  trControl = trainControl(
    method = "cv",
    number = 10
  )
)


## Predicción en el conjunto de test
nb_predictions <- predict(
  nb_model,
  newdata = test_filtered
)

nb_predictions


## Matriz de confusión
conf_nb <- confusionMatrix(
  nb_predictions,
  test_filtered$Class
)

conf_nb


## Imprimirlas específicamente
conf_nb$table
conf_nb$overall["Accuracy"]
conf_nb$byClass[c("Sensitivity", "Specificity", "F1")]


## Probabilidades
probabilities_nb <- predict(
  nb_model,
  newdata = test_filtered,
  type = "prob"
)

head(probabilities_nb)


## Visualización: usamos un UMAP coloreado por clase predicha para la visualización

# Añadimos predicciones al dataframe UMAP
umap_nb_df <- umap_df
umap_nb_df$Predicted <- nb_predictions[match(rownames(umap_df), rownames(test_filtered))]

# Filtramos solo muestras de test
umap_nb_df_test <- umap_nb_df[rownames(umap_nb_df) %in% rownames(test_filtered), ]

# Gráfico
ggplot(umap_nb_df_test, aes(x = UMAP1, y = UMAP2)) +
  geom_point(aes(color = Predicted, shape = Class), size = 3) +
  theme_classic() +
  labs(
    title = "UMAP de expresión génica\nColoreado por clase predicha (Naive Bayes)",
    color = "Clase predicha",
    shape = "Clase real"
  )



