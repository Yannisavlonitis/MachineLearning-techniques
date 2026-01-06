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

# -------------- Metodos de aprendizaje supervisado: --------------------------

#                               
# _______________________________ kNN __________________________________________

library(caret)

set.seed(42) #fijamos una semilla para reproducibilidad del codigo

trainIndex <- createDataPartition(df_final$Clase, p = 0.8, list = FALSE)
df_final$Clase <- as.factor(df_final$Clase)
trainData <- df_final[trainIndex, ]
testData <- df_final[-trainIndex, ]

set.seed(42)
knnModel <- train(Clase ~ .,
                  data = trainData,
                  method = "knn",
                  trControl = trainControl(method = "cv", number = 10),
                  preProcess = c("zv", "center", "scale"),
                  tuneLength = 30)
plot(knnModel)

knnModel
# Utiliza k=11 como k optimo


# Realizar predicciones en el conjunto de prueba utilizando el modelo entrenado

# El modelo kNN calcula la distancia entre el punto test y todos los puntos de
# entrenamiento, ordenando las distancias y seleccionando los k vecinos mas
# cercanos para cada punto. Finalmente, asigna la etiqueta al punto basandose
# en la mayoria de etiquetas de entre sus vecinos

set.seed(42)
predictions_knn <- predict(knnModel, newdata = testData )


# Evaluar la precisión del modelo gracias a la matriz de confusión
conf_knn <- confusionMatrix(predictions_knn, testData$Clase)
conf_knn

#Accuracy : 0.9937, muy buen modelo

# Obtener probabilidades
probabilities_knn <- predict(knnModel, newdata = testData, type = "prob")
probabilities_knn

