
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


## Entrenamos el modelo con validación cruzada de 10 folds
nb_model <- train(
  Class ~ .,
  data = trainData,
  method = "nb",
  trControl = trainControl(
    method = "cv",
    number = 10
  )
)



## Durante el entrenamiento del modelo se obtiene el siguiente error: 

## Something is wrong; all the Accuracy metric values are missing:
##  Accuracy       Kappa    
##  Min.   : NA   Min.   : NA  
##  1st Qu.: NA   1st Qu.: NA  
##  Median : NA   Median : NA  
##  Mean   :NaN   Mean   :NaN  
##  3rd Qu.: NA   3rd Qu.: NA  
##  Max.   : NA   Max.   : NA  
##  NA's   :2     NA's   :2    
##  Error: Stopping

## Esto se debe a que el método Naive Bayes utilizado (modelo gaussiano) asume que cada variable sigue una
## distribución normal dentro de cada clase, por lo que necesita estimar correctamente la media y la varianza, 
## además de asumir independencia entre genes.

## Sin embargo, en datos de expresión génica es frecuente encontrar genes con expresión muy baja o nula en determinadas 
## clases. En nuestro conjunto de datos, algunos genes presentan varianza cero dentro de alguna clase, como 
## ocurre con el gen DBX2 en la clase HPB, donde todas las muestras tienen el mismo valor de expresión:

boxplot(
  DBX2 ~ Class,
  data = data_imputed,
  main = "Expresión de DBX2 por clase",
  ylab = "Expresión",
  las = 2
)

## Esto hace que el modelo no pueda calcular probabilidades y falle. Otros métodos como Random Forest, SVM o k-NN,
## probablemente funcionen mejor en este caso porque no necesitan estimar varianzas.


