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


# Algoritmo de Machine Learning: Gradient Boosting Machine (GBM)

library(caret)

# División de datos para el Entrainamiento (70%) y Prueba (30%)

set.seed(123) # Para que los resultados sean reproducibles

train_index <- createDataPartition(df_final$Clase, p = 0.7, list = FALSE)
train <- df_final[train_index, ]
test  <- df_final[-train_index, ]

cat("Pacientes en el grupo de ENtrainAMIENTO (Training):")
print(cbind(Frecuencia = table(train$Clase), 
            Porcentaje = round(prop.table(table(train$Clase)) * 100, 2)))

cat("Pacientes en el grupo de PRUEBA (Test):")
print(cbind(Frecuencia = table(test$Clase), 
            Porcentaje = round(prop.table(table(test$Clase)) * 100, 2)))


# Entrainamiento del modelo GBM (trada unos 3 minutos!! se ha reducido el tiempo)
gbm_model <- train(
  Clase ~ ., 
  data = train,
  method = "gbm",
  trControl = trainControl(method = "cv", number = 5), # Reducido a 5 para ganar velocidad
  preProcess = c("center", "scale"), # Normalización
  tuneLength = 5,  # Prueba 5 valores impares diferentes de 'k' (pocos para ganar velocidad)
  verbose = FALSE
)


# Resumen del modelo entrainado
print(gbm_model)
plot(gbm_model)

mejor_k <- gbm_model$bestTune # Extraer y mostrar el mejor k automáticamente
print(paste("El algoritmo ha seleccionado automáticamente el valor óptimo de k =", mejor_k))


# Importancia de variables (Genes más relevantes)
importancia <- varImp(gbm_model, scale = TRUE)

# Representación de los 10 genes con mayor poder predictivo
plot(importancia, top = 10, main = "Top 10 Genes - Importancia GBM")



# Predicciones y evaluación en el conjunto de prueba

# Predicción de etiquetas
predictions <- predict(gbm_model, newdata = test)

# Matriz de confusión para ver sensibilidad y especificidad
cm <- confusionMatrix(predictions, as.factor(test$Clase))
print(cm)

# Probabilidades de clase (para curvas ROC)
probabilities_gbm <- predict(
  gbm_model, 
  newdata = test, 
  type = "prob"
)

head(probabilities_gbm)

# para guardar las probabilidades para el ROC
saveRDS(probabilities_gbm, file = "gbm_probabilities.rds")

# para cargar las probabilidades
probabilidades_recuperadas <- readRDS("gbm_probabilities.rds")

head(probabilidades_recuperadas)
