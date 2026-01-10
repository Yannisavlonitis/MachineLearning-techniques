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


# -------------- Metodos de aprendizaje supervisado: --------------------------

#                               
# _________________________ DECCISION TREE ____________________________________

# librerías 
library(rpart)
library(rpart.plot)
library(caret)

# la variable objetivo (Clase) sea un factor
df_final$Clase <- as.factor(df_final$Clase)

# Visualización de las clases y porcentaje de las dos categorías
print(cbind(Frecuencia = table(df_final$Clase), 
            Porcentaje = round(prop.table(table(df_final$Clase)) * 100, 2)))



# División de datos para el Entrenamiento (70%) y Prueba (30%)

set.seed(123) # Para que los resultados sean reproducibles

train_index <- createDataPartition(df_final$Clase, p = 0.7, list = FALSE)
tren <- df_final[train_index, ]
test  <- df_final[-train_index, ]

cat("Pacientes en el grupo de ENTRENAMIENTO (Training):")
print(cbind(Frecuencia = table(tren$Clase), 
            Porcentaje = round(prop.table(table(tren$Clase)) * 100, 2)))

cat("Pacientes en el grupo de PRUEBA (Test):")
print(cbind(Frecuencia = table(test$Clase), 
            Porcentaje = round(prop.table(table(test$Clase)) * 100, 2)))



# Entrenamiento del Árbol de Decisión
# con a fórmula Clase ~ . queremos predecir la Clase usando todos los genes
modelo_arbol <- rpart(Clase ~ ., data = tren, method = "class")

# Visualización del Árbol
# Esto te permitirá ver qué genes son los más importantes para clasificar las muestras
prp(modelo_arbol, type = 2, extra = 104, 
    main = "Árbol de Decisión: Expresión Génica", 
    under = TRUE, faclen = 0, box.palette = "RdYlGn")

#  Predicciones con el conjunto de prueba
predicciones <- predict(modelo_arbol, test, type = "class")

# Evaluación del Modelo
matriz_confusion <- confusionMatrix(predicciones, test$Clase)
print(matriz_confusion)

# Ver importancia de las variables (los genes que más influyen)
importancia <- modelo_arbol$variable.importance
head(importancia, 10)




library(pROC)

# Obtenemos probabilidades del árbol
prob_arbol <- predict(modelo_arbol, test, type = "prob")
# Suponiendo que la clase positiva es la segunda columna
roc_arbol <- roc(test$Clase, prob_arbol[, 2])


# para guardar las probabilidades para el ROC
saveRDS(prob_arbol, file = "arbol_probabilities.rds")

# para cargar las probabilidades de gbm
prob_gbm <- readRDS("gbm_probabilities.rds")
head(prob_gbm)

roc_gbm <- roc(test$Clase, prob_gbm[, 2])


# Graficar ambas para comparar
plot(roc_arbol, col = "blue", lwd = 2, main = "Comparativa Curvas ROC")
plot(roc_gbm, col = "red", lwd = 2, add = TRUE)
legend("bottomright", legend = c(
  paste("Árbol (AUC =", round(auc(roc_arbol), 3), ")"),
  paste("GBM (AUC =", round(auc(roc_gbm), 3), ")")), 
  col = c("blue", "red"), lwd = 2)
abline(a = 0, b = 1, lty = 2, col = "grey")

