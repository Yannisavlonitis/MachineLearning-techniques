# 1. PREPARACIÓN PREVIA (Asegurar formato)
# Es vital que la Clase sea factor para que los modelos y las curvas funcionen
df_final$Clase <- as.factor(df_final$Clase)

# 2. CÁLCULO DE AUC Y CURVAS ROC
# Nota: Asumo que testData es una partición de df_final
# Sustituimos 'primaryormetastasis' por 'Clase'

# Lista de modelos y sus probabilidades para automatizar un poco
modelos_list <- list(
  knn = probabilities_knn[,2],
  svm_lin = probabilities_svm_linear[,2],
  svm_ker = probabilities_svm_kernel[,2],
  svm_pol = probabilities_svm_kernelpol[,2],
  dt = probabilities_dt[,2],
  lda = probabilities_lda$posterior[,2],
  rda = probabilities_rda$posterior[,2]
)

# Calcular ROCs
roc_knn      <- roc(testData$Clase, modelos_list$knn)
roc_svm_lin  <- roc(testData$Clase, modelos_list$svm_lin)
roc_svm_ker  <- roc(testData$Clase, modelos_list$svm_ker)
roc_svm_pol  <- roc(testData$Clase, modelos_list$svm_pol)
roc_dt       <- roc(testData$Clase, modelos_list$dt)
roc_lda      <- roc(testData$Clase, modelos_list$lda)

# Manejo de NA para RDA (común en datasets genómicos con alta colinealidad)
valid_rda    <- !is.na(probabilities_rda$posterior[, 2])
roc_rda      <- roc(testData$Clase[valid_rda], modelos_list$rda[valid_rda])

# --- VISUALIZACIÓN ROC ---
plot(roc_knn, col = "blue", main = "Curvas ROC: Clasificación de Tumores", lwd = 2)
plot(roc_svm_lin, col = "red", add = TRUE, lwd = 2)
plot(roc_svm_ker, col = "green", add = TRUE, lwd = 2)
plot(roc_svm_pol, col = "orange", add = TRUE, lwd = 2)
plot(roc_dt, col = "purple", add = TRUE, lwd = 2)
plot(roc_lda, col = "pink", add = TRUE, lwd = 2)
plot(roc_rda, col = "yellow", add = TRUE, lwd = 2)

legend("bottomright", 
       legend = c(paste("k-NN:", round(auc(roc_knn), 2)),
                  paste("SVM Lin:", round(auc(roc_svm_lin), 2)),
                  paste("SVM Ker:", round(auc(roc_svm_ker), 2)),
                  paste("SVM Pol:", round(auc(roc_svm_pol), 2)),
                  paste("Tree:", round(auc(roc_dt), 2)),
                  paste("LDA:", round(auc(roc_lda), 2)),
                  paste("RDA:", round(auc(roc_rda), 2))),
       col = c("blue", "red", "green", "orange", "purple", "pink", "yellow"), 
       lwd = 2, cex = 0.8)

