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
# _________________ Gradient Boosting Machine (GBM) ____________________________

# librerías 
library(gbm)
library(caret)

# Asegurar que la clase sea factor
df_final$Clase <- as.factor(df_final$Clase)

# División entrenamiento/prueba
set.seed(123)
train_index <- createDataPartition(df_final$Clase, p = 0.7, list = FALSE)
tren <- df_final[train_index, ]
test  <- df_final[-train_index, ]


# Entrenar el modelo gbm
modelo_gbm <- gbm(
  formula = Clase ~ ., 
  data = tren, 
  distribution = "multinomial", 
  n.trees = 500,               
  interaction.depth = 3,       
  shrinkage = 0.1,             
  cv.folds = 5
)

# Encontrar el número óptimo de árboles
best_iter <- gbm.perf(modelo_gbm, method = "cv")



# Visualización
png("importancia_genes_gbm.png", width = 900, height = 1100, res = 120)

par(mar = c(5, 12, 4, 2)) 
# Graficamos los 20 genes más importantes
summary(modelo_gbm, n.trees = best_iter, cBars = 20, las = 2)

dev.off()



# Predicciones
pred_probs <- predict(modelo_gbm, test, n.trees = best_iter, type = "response")
pred_clase <- colnames(pred_probs)[apply(pred_probs, 1, which.max)]
pred_clase <- factor(pred_clase, levels = levels(df_final$Clase))

# Matriz de confusión
matriz_confusion <- confusionMatrix(pred_clase, test$Clase)
print(matriz_confusion)



library(pROC)

# --- ROC PARA GBM ---
# Usamos las probabilidades que calculamos anteriormente (pred_probs)
# pred_probs suele ser una matriz; tomamos la columna de la clase de interés
roc_gbm <- roc(test$Clase, pred_probs[, , 1][, 2]) # Ajustar índices según niveles de clase

# --- VISUALIZACIÓN ---
# Graficar ambas para comparar
plot(roc_arbol, col = "blue", lwd = 2, main = "Comparativa Curvas ROC")
plot(roc_gbm, col = "red", lwd = 2, add = TRUE)
legend("bottomright", legend = c(
  paste("Árbol (AUC =", round(auc(roc_arbol), 3), ")"),
  paste("GBM (AUC =", round(auc(roc_gbm), 3), ")")), 
  col = c("blue", "red"), lwd = 2)
abline(a = 0, b = 1, lty = 2, col = "grey")