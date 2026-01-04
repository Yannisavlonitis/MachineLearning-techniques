# 1. PREPARACIÓN Y LIBRERÍAS
library(pROC)
library(PRROC)
library(mda)      
library(gbm)      
library(e1071)    

df_final$Clase <- as.factor(df_final$Clase)

# 2. LISTA DE PROBABILIDADES
# Corrección 1: Aseguramos que 'nb' y otros extraigan la columna correcta (prob. de la clase positiva)
# Nota: Verifica si tu clase positiva es la 1 o la 2 en levels(testData$Clase)
modelos_probs <- list(
  knn      = probabilities_knn[,2],
  svm_lin  = probabilities_svm_linear[,2],
  svm_rad  = probabilities_svm_kernel[,2],
  svm_pol  = probabilities_svm_kernelpol[,2],
  dt       = probabilities_dt[,2],
  lda      = probabilities_lda$posterior[,2],
  qda      = probabilities_qda$posterior[,2],
  rda      = probabilities_rda$posterior[,2],
  # Corrección 2: predict de mda (FDA) devuelve matriz; necesitamos la columna de la clase positiva
  fda      = predict(model_fda, testData, type = "posterior")[,2],
  nb       = probabilities_nb[,2],
  boosting = predict(model_gbm, testData, n.trees = 100, type = "response")
)

# 3. CÁLCULO DE CURVAS ROC
# Corrección 3: En pROC, a veces hay que especificar los niveles para evitar confusión de dirección
calcular_roc <- function(prob, labels) {
  valid <- !is.na(prob)
  # Usamos levels(labels) para que pROC sepa exactamente el orden de las clases
  roc(labels[valid], prob[valid], quiet = TRUE, levels = levels(labels))
}

rocs <- lapply(modelos_probs, calcular_roc, labels = testData$Clase)

# 4. VISUALIZACIÓN ROC UNIFICADA
colores <- c("blue", "red", "darkgreen", "orange", "purple", "deeppink", "gold", "brown", "cyan", "darkmagenta", "darkgray")

# Aseguramos que el gráfico se limpie antes de dibujar
plot.new() 
plot(rocs[[1]], col = colores[1], main = "Comparativa ROC: Diagnóstico Genético", lwd = 2)

for(i in 2:length(rocs)) {
  # Usamos la función genérica lines() o plot(..., add=TRUE)
  plot(rocs[[i]], col = colores[i], add = TRUE, lwd = 2)
}

# Leyenda con AUC redondeado
legend_labels <- sapply(names(rocs), function(n) {
  paste(toupper(n), ":", round(as.numeric(auc(rocs[[n]])), 2))
})

legend("bottomright", legend = legend_labels, col = colores, lwd = 2, cex = 0.6, ncol = 2, bty = "n")
