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










###############################################################################
### METODOS SUPERVISADOS ESCOGIDOS
###############################################################################

# librerias a usar
library(caret) # para knn, gbm, LDA, decision tree y partición de datos y matriz de confusion
library(MASS)  # Para LDA
library(gbm)   # Para GBM
library(rpart) # Para decision tree
library(rpart.plot) # Para decision tree



# Preparación de datos
df_final$Clase <- as.factor(df_final$Clase) # Convertir Clase a factor

set.seed(42) #fijamos una semilla para reproducibilidad del codigo

# Partición de datos (80% entrenamiento, 20% prueba)
train_index <- createDataPartition(df_final$Clase, p = 0.8, list = FALSE)

train_data <- df_final[train_index, ]
test_data <- df_final[-train_index, ]


# Verificación de la partición de datos
# Frecuencias y proporciones para training
train_tab <- table(train_data$Clase)
train_prop <- prop.table(train_tab)

# Frecuencias y proporciones para testing
test_tab <- table(test_data$Clase)
test_prop <- prop.table(test_tab)

# Crear tabla combinada para comparar si ambas partes tienen la misma proporcion de datos
summary_table <- data.frame(
  Clase = names(train_tab),  # asumimos que ambas tablas tienen las mismas clases
  Train_Frecuencia = as.vector(train_tab),
  Train_Porcentaje = round(100 * as.vector(train_prop), 2),
  Test_Frecuencia = as.vector(test_tab),
  Test_Porcentaje = round(100 * as.vector(test_prop), 2)
)

print(summary_table)










# _______________________________ kNN __________________________________________

# modelo del knn (tarda aproximadamente 30 segundos!!)
set.seed(42)
knnModel <- train(Clase ~ .,
                  data = train_data,
                  method = "knn",
                  trControl = trainControl(method = "cv", number = 10),
                  preProcess = c("zv", "center", "scale"),
                  tuneLength = 30)
plot(knnModel)


mejor_k_knn <- knnModel$bestTune # Extraer y mostrar el mejor k automáticamente
print(paste("El modelo KNN ha seleccionado automáticamente el valor óptimo de k =", mejor_k_knn))


# Importancia de variables (Genes más relevantes)
importancia_knn <- varImp(knnModel, scale = TRUE)

# Representación de los 10 genes con mayor poder predictivo
plot(importancia_knn, top = 10, main = "Top 10 Genes - Importancia KNN")


# Realizar predicciones en el conjunto de prueba utilizando el modelo entrenado

# El modelo kNN calcula la distancia entre el punto test y todos los puntos de
# entrenamiento, ordenando las distancias y seleccionando los k vecinos mas
# cercanos para cada punto. Finalmente, asigna la etiqueta al punto basandose
# en la mayoria de etiquetas de entre sus vecinos

set.seed(42)
predictions_knn <- predict(knnModel, newdata = test_data )


# Evaluar la precisión del modelo gracias a la matriz de confusión
conf_knn <- confusionMatrix(predictions_knn, test_data$Clase)
conf_knn

#Accuracy : 0.9937, muy buen modelo

# Obtener probabilidades, util para curvas ROC
prob_knn <- predict(knnModel, newdata = test_data, type = "prob")
prob_knn










# _______________________________ LDA __________________________________________

# Teoría:
# Para contruir el modelo, hay que relacionar varibles como el radius, compactness... 
# con B (benigno) o M (metastasico), es decir, se va a intentar entrenar un modelo 
# que pueda predecir la malignidad de un tumor a partir de variables observadas.
# Los predicores (X) van a ser las medidas de los tumores, y el vector respuesta Y (la malignidad)
# X = MEDIDAS Y = DIAGNOSIS para el training testing 


#----------------TRATAMIENTO DE VARIABLES PROBLEMÁTICAS ------------------------
# Identificar variables con varianza cero (constantes), no aportan información discriminante
var_no_cero <- apply(train_data[, -ncol(train_data)], 2, var) != 0

# Mostrar qué variables tienen varianza cero
variables_cero <- names(var_no_cero)[!var_no_cero]  # !var_no_cero = varianza = 0
print("Variables con varianza cero (constantes):")
print(variables_cero)


# Data filtrada (solo variables con varianza)
train_filtered <- train_data[, var_no_cero]
test_filtered <- test_data[, var_no_cero]

# Limpiar nombres de columnas (eliminar guiones y caracteres problemáticos)
colnames(train_filtered) <- make.names(colnames(train_filtered))
colnames(test_filtered) <- make.names(colnames(test_filtered))

# Añadir la clase de nuevo
train_filtered$Clase <- train_data$Clase
test_filtered$Clase <- test_data$Clase


#----------------------ESCALADO DE VARIABLES------------------------------------
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

#----------------------CONSTRUCCIÓN DE LA FÓRMULA-------------------------------

genes <- colnames(train_lda)[1:(ncol(train_lda)-1)]
formula_lda <- as.formula(paste("Clase ~", paste(genes, collapse = "+")))



#--------------------------MODELO LDA-------------------------------------------

# modelo del LDA (tarda 0 segundos!!)
set.seed(42) 
lda_model <- lda(formula_lda, data = train_lda)


# Predicciones en conjunto de prueba
lda_predictions <- predict(lda_model, newdata = test_lda)

# Clases predichas vs reales
predicted_lda <- lda_predictions$class
true_classes <- test_lda$Clase

# Matriz de confusión
conf_lda <- confusionMatrix(predicted_lda, true_classes)
conf_lda

# Guardar métricas
accuracy_lda <- conf_lda$overall['Accuracy']
kappa_lda <- conf_lda$overall['Kappa']

# Obtener probabilidades
prob_LDA <- predict(lda_model, newdata = test_lda)$posterior
head(prob_LDA)




#---------------------------VISUALIZACIÓN LDA-----------------------------------
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

# Reemplazar "-" por "." en todos los nombres de test_data
colnames(test_data) <- gsub("-", ".", colnames(test_data))

# Verificamos
colnames(test_data)










# _______________________________ gbm __________________________________________
# Algoritmo de Machine Learning: Gradient Boosting Machine (GBM)

# Entrainamiento del modelo GBM (trada unos 3 minutos!! se ha reducido el tiempo)
gbm_model <- train(
  Clase ~ ., 
  data = train_data,
  method = "gbm",
  trControl = trainControl(method = "cv", number = 5), # Reducido a 5 para ganar velocidad
  preProcess = c("center", "scale"), # Normalización
  tuneLength = 5,  # Prueba 5 valores impares diferentes de 'k' (pocos para ganar velocidad)
  verbose = FALSE
)


# Resumen del modelo entrainado
print(gbm_model)
plot(gbm_model)

mejor_k_gbm <- gbm_model$bestTune # Extraer y mostrar el mejor k automáticamente
print(paste("El Modelo GBM ha seleccionado automáticamente el valor óptimo de k =", mejor_k_gbm))


# Importancia de variables (Genes más relevantes)
importancia_gbm <- varImp(gbm_model, scale = TRUE)

# Representación de los 10 genes con mayor poder predictivo
plot(importancia_gbm, top = 10, main = "Top 10 Genes - Importancia GBM")



# Predicciones y evaluación en el conjunto de prueba
# Predicción de etiquetas
predicted_gbm <- predict(gbm_model, newdata = test_data)

# Matriz de confusión para ver sensibilidad y especificidad
conf_gbm <- confusionMatrix(predicted_gbm, as.factor(test_data$Clase))
conf_gbm

# Probabilidades de clase (para curvas ROC)
prob_gbm <- predict(gbm_model, newdata = test_data, type = "prob")
head(prob_gbm)










# _________________________CURVAS ROC____________________________________
library(pROC)
# multiclass.roc calcula la curva ROC promedio para clasificación multiclase
# no puede ser roc normal, debido hay más de una clase (no binomial)
roc_knn <- multiclass.roc(test_data$Clase, prob_knn)
roc_LDA <- multiclass.roc(test_lda$Clase, prob_LDA)
roc_gbm <- multiclass.roc(test_data$Clase, prob_gbm)

# Tabla de AUCs metodos seleccionados
tabla_auc_seleccionados <- data.frame(
  Modelo = c("KNN", "LDA", "GBM"),
  AUC = c(as.numeric(roc_knn$auc),
          as.numeric(roc_LDA$auc),
          as.numeric(roc_gbm$auc))
)
tabla_auc_seleccionados




# ----------------------- ROC One-vs-Rest por clase ----------------------------
# Para problemas multiclase es útil ver ROC para cada clase individual

# Obtener nombres de clases
clases <- levels(test_data$Clase)

# Función para calcular ROC One-vs-Rest
roc_por_clase <- function(prob_matrix, true_labels, clases) {
  lapply(clases, function(clase) { 
    roc(response = true_labels == clase, # Cada clase se toma como "positiva" y el resto como "negativa"
        predictor = prob_matrix[, clase])
  }) -> roc_list
  names(roc_list) <- clases
  return(roc_list)
}

# Calcular ROC One-vs-Rest para cada modelo
roc_knn_list <- roc_por_clase(prob_knn, test_data$Clase, clases)
roc_LDA_list <- roc_por_clase(prob_LDA, test_lda$Clase, clases)
roc_gbm_list <- roc_por_clase(prob_gbm, test_data$Clase, clases)


# -------------------- Graficar ROC One-vs-Rest por clase ----------------------
# Definir colores para cada modelo
colores <- c(KNN="#377EB8", LDA="#4DAF4A", GBM="#E41A1C")

# Graficar una figura por clase
for(clase in clases){
  
  # Graficar ROC del modelo KNN
  plot(roc_knn_list[[clase]], col=colores["KNN"], lwd=2,
       main=paste("Curvas ROC One-vs-Rest - Clase:", clase))
  plot(roc_LDA_list[[clase]], col=colores["LDA"], lwd=2, add=TRUE)
  plot(roc_gbm_list[[clase]], col=colores["GBM"], lwd=2, add=TRUE)
  
  legend("bottomright",
         legend=c(
           paste0("KNN (AUC=", round(auc(roc_knn_list[[clase]]),3), ")"),
           paste0("LDA (AUC=", round(auc(roc_LDA_list[[clase]]),3), ")"),
           paste0("GBM (AUC=", round(auc(roc_gbm_list[[clase]]),3), ")")
         ),
         col=colores, lwd=2)
  abline(0,1, lty=2, col="darkgrey") # Línea diagonal de referencia (clasificador aleatorio)
}


# ------------------------- tabla con AUC por clase ----------------------------
# Extraer AUC de cada ROC One-vs-Rest
auc_knn <- sapply(roc_knn_list, auc)
auc_LDA <- sapply(roc_LDA_list, auc)
auc_gbm <- sapply(roc_gbm_list, auc)

# Crear tabla comparativa
tabla_auc_binaria <- data.frame(
  Clase = clases,
  KNN = auc_knn,
  LDA = auc_LDA,
  GBM = auc_gbm
)

# Mostrar tabla
print(tabla_auc_binaria)















###############################################################################
### METODOS SUPERVISADOS NO ESCOGIDOS
###############################################################################

# _________________________ DECCISION TREE ____________________________________

# Entrenamiento del Árbol de Decisión
# con a fórmula Clase ~ . queremos predecir la Clase usando todos los genes
modelo_arbol <- rpart(Clase ~ ., data = train_data, method = "class")

# Visualización del Árbol
# Esto te permitirá ver qué genes son los más importantes para clasificar las muestras
prp(modelo_arbol, type = 2, extra = 104, 
    main = "Árbol de Decisión: Expresión Génica", 
    under = TRUE, faclen = 0, box.palette = "RdYlGn")

#  Predicciones con el conjunto de prueba
predicted_arbol <- predict(modelo_arbol, test_data, type = "class")

# Evaluación del Modelo
conf_arbol <- confusionMatrix(predicted_arbol, test_data$Clase)
conf_arbol

# Ver importancia de las variables (los genes que más influyen)
importancia <- modelo_arbol$variable.importance
head(importancia, 10)


# Obtenemos probabilidades del árbol
prob_arbol <- predict(modelo_arbol, test_data, type = "prob")
head(prob_arbol)










# ________________________________ SVM _________________________________________

#cargar librerias
library(caret)
library(tidyverse)
library(glmnet)
library(dplyr)
library(gridExtra)
library(e1071)

# -----------------------------SVM lineal---------------------------------------
# Entrainamiento del modelo SVM lineal (trada unos 3 minutos!!)
svmModelLineal <- train(Clase ~.,
                        data = train_data,
                        method = "svmLinear",
                        trControl = trainControl(method = "cv", number = 10),
                        preProcess = c("center", "scale"),
                        tuneGrid = expand.grid(C = seq(0, 2, length = 30)), #C grande lleva al sobreajuste, C pequeño al infraajuste. Prueba valores de "C" entre 0 y 20 cada 2
                        prob.model = TRUE)

#The final value used for the model was C = 0.2758621.
svmModelLineal
plot(svmModelLineal)

# Predicciones en el conjunto de prueba utilizando el modelo entrenado
predicted_svm_lineal <- predict(svmModelLineal, newdata = test_data )
head(predicted_svm_lineal)

# Evaluar la precisión del modelo utilizando la matriz de confusión
conf_SVM_lineal <- confusionMatrix(predicted_svm_lineal, test_data$Clase)
conf_SVM_lineal

prob_svm_Lineal <- predict(svmModelLineal, newdata = test_data, type = "prob")
head(prob_svm_Lineal)






# --------------------------SVM Kernel radial-----------------------------------
# Entrainamiento del modelo SVM Kernel radial (trada 1 minuto!!)
svmModelKernel <- train(Clase ~.,
                        data = train_data,
                        method = "svmRadial", #metodo radial
                        trControl = trainControl(method = "cv", number = 10),
                        preProcess = c("center", "scale"),
                        tuneLength = 10,
                        prob.model = TRUE) 

#The final values used for the model were sigma = 0.001241186 and C = 16.
svmModelKernel 
plot(svmModelKernel)

# Predicciones en el conjunto de prueba utilizando el modelo entrenado
predicted_SVM_k_radial <- predict(svmModelKernel, newdata = test_data )
predicted_SVM_k_radial

# Evaluar la precisión del modelo utilizando la matriz de confusión
conf_SVM_k_radial <- confusionMatrix(predicted_SVM_k_radial, test_data$Clase)
conf_SVM_k_radial

# Análisis: 
#Modelo muy bueno para detectar AGH, CFB, CGC (alta sensibilidad)
#Muy bueno para descartar: AGH, CFB, CHC, HPB (alta especificidad)
#El modelo está bien entrenado para la clasificacion de: AGH y CFB (prevalencia real y detectada, similares)
# En el caso de CGC tiene un elevado porcentaje de falsos positivos (56% de aciertos al id +)
#En el caso de CHC sucede al revés, elevado porcentaje de falsos negativos (85% de aciertos al id-)

prob_svm_k_radial <- predict(svmModelKernel, newdata = test_data, type = "prob")
head(prob_svm_k_radial)






# --------------------------SVM Kernel Polinomial-----------------------------------
# Entrainamiento del modelo SVM Kernel Polinomial (trada aprox 10 minutos!!)
svmModelKernelPolynomial <- train(Clase ~.,
                                  data = train_data,
                                  method = "svmPoly",
                                  trControl = trainControl(method = "cv", number = 10),
                                  preProcess = c("center", "scale"),
                                  tuneLength = 5,
                                  prob.model = TRUE) 
# The final values used for the model were degree = 1, scale = 0.001 and C = 0.25.
svmModelKernelPolynomial
plot(svmModelKernelPolynomial)


# Realizar predicciones en el conjunto de prueba utilizando el modelo entrenado
predicted_svm_k_polinomial <- predict(svmModelKernelPolynomial, newdata = test_data )
predicted_svm_k_polinomial

# Evaluar la precisión del modelo utilizando la matriz de confusión
conf_svm_k_polinomial <- confusionMatrix(predicted_svm_k_polinomial, test_data$Clase)
conf_svm_k_polinomial

# SVM kernelpol
prob_svm_k_polinomial <- predict(svmModelKernelPolynomial, newdata = test_data, type = "prob")
head(prob_svm_k_polinomial)










# ________________________________ NAIVE BAYES _________________________________________
## Entrenamiento Naive Bayes (trada unos 3 minutos!!)
nb_model <- train(
  Clase ~ .,
  data = train_filtered,
  method = "nb",
  trControl = trainControl(
    method = "cv",
    number = 10
  )
)


## Predicción en el conjunto de test
predicted_nb <- predict(nb_model, newdata = test_filtered)
predicted_nb

## Matriz de confusión
conf_nb <- confusionMatrix(predicted_nb, test_filtered$Clase)
print(conf_nb)

# Obtenemos probabilidades para curvas ROC
prob_nb <- predict(nb_model, newdata = test_filtered, type = "raw") 
head(prob_nb)










# _________________________CURVAS ROC (TODOS)____________________________________
library(pROC)

# Extraer AUC para cada modelo
roc_knn <- multiclass.roc(test_data$Clase, prob_knn)
roc_LDA <- multiclass.roc(test_lda$Clase, prob_LDA)
roc_gbm <- multiclass.roc(test_data$Clase, prob_gbm)
roc_arbol <- multiclass.roc(test_data$Clase, prob_arbol)
roc_svm_lineal <- multiclass.roc(test_data$Clase, prob_svm_Lineal)
roc_svm_k_radial <- multiclass.roc(test_data$Clase, prob_svm_k_radial)
roc_svm_k_polinomial <- multiclass.roc(test_data$Clase, prob_svm_k_polinomial)



tabla_auc <- data.frame(
  Modelo = c(
    "KNN",
    "LDA",
    "GBM",
    "Árbol de decisión",
    "SVM lineal",
    "SVM kernel radial",
    "SVM kernel polinomial"
  ),
  AUC = c(
    roc_knn$auc,
    roc_LDA$auc,
    roc_gbm$auc,
    roc_arbol$auc,
    roc_svm_lineal$auc,
    roc_svm_k_radial$auc,
    roc_svm_k_polinomial$auc
  )
)

tabla_auc


clases <- levels(test_data$Clase)  # Extraer niveles de clase

# ---------------------- Función ROC One-vs-Rest -------------------------------
roc_por_clase <- function(prob_matrix, true_labels, clases) {
  # Calcula ROC binario One-vs-Rest para cada clase
  roc_list <- lapply(clases, function(clase) {
    roc(response = true_labels == clase, predictor = prob_matrix[, clase])
  })
  names(roc_list) <- clases
  return(roc_list)
}

# Calcular ROC por clase para cada modelo 

roc_knn_list <- roc_por_clase(prob_knn, test_data$Clase, clases)
roc_LDA_list <- roc_por_clase(prob_LDA, test_lda$Clase, clases)
roc_gbm_list <- roc_por_clase(prob_gbm, test_data$Clase, clases)
roc_arbol_list <- roc_por_clase(prob_arbol, test_data$Clase, clases)
roc_svm_lineal_list <- roc_por_clase(prob_svm_Lineal, test_data$Clase, clases)
roc_svm_k_radial_list <- roc_por_clase(prob_svm_k_radial, test_data$Clase, clases)
roc_svm_k_polinomial_list <- roc_por_clase(prob_svm_k_polinomial, test_data$Clase, clases)



# Crear tabla AUC por clase y modelo

# Extraer AUC para cada modelo y clase
auc_knn <- sapply(roc_knn_list, auc)
auc_LDA <- sapply(roc_LDA_list, auc)
auc_gbm <- sapply(roc_gbm_list, auc)
auc_arbol <- sapply(roc_arbol_list, auc)
auc_svm_lineal <- sapply(roc_svm_lineal_list, auc)
auc_svm_k_radial <- sapply(roc_svm_k_radial_list, auc)
auc_svm_k_polinomial <- sapply(roc_svm_k_polinomial_list, auc)

tabla_auc_binaria <- data.frame(
  Clase = clases,
  KNN = round(auc_knn, 3),
  LDA = round(auc_LDA, 3),
  GBM = round(auc_gbm, 3),
  Arbol = round(auc_arbol, 3),
  SVM_Lineal = round(auc_svm_lineal, 3),
  SVM_Radial = round(auc_svm_k_radial, 3),
  SVM_Polinomial = round(auc_svm_k_polinomial, 3)
)

tabla_auc_binaria




# ---------------------- Graficar ROC por clase ---------------------------------

colores <- c("KNN"="#377EB8", "LDA"="gold", "GBM"="#E41A1C",
             "Arbol"="#4DAF4A", "SVM_Lineal"="hotpink",
             "SVM_Radial"="deeppink", "SVM_Polinomial"="darkmagenta")

for(clase in clases){
  plot(roc_knn_list[[clase]], col=colores["KNN"], lwd=2,
       main=paste("Curvas ROC One-vs-Rest - Clase:", clase))
  plot(roc_LDA_list[[clase]], col=colores["LDA"], lwd=2, add=TRUE)
  plot(roc_gbm_list[[clase]], col=colores["GBM"], lwd=2, add=TRUE)
  plot(roc_arbol_list[[clase]], col=colores["Arbol"], lwd=2, add=TRUE)
  plot(roc_svm_lineal_list[[clase]], col=colores["SVM_Lineal"], lwd=2, add=TRUE)
  plot(roc_svm_k_radial_list[[clase]], col=colores["SVM_Radial"], lwd=2, add=TRUE)
  plot(roc_svm_k_polinomial_list[[clase]], col=colores["SVM_Polinomial"], lwd=2, add=TRUE)
  
  legend("bottomright",
         legend = c(
           paste("KNN (AUC=", round(auc(roc_knn_list[[clase]]),3), ")"),
           paste("LDA (AUC=", round(auc(roc_LDA_list[[clase]]),3), ")"),
           paste("GBM (AUC=", round(auc(roc_gbm_list[[clase]]),3), ")"),
           paste("Árbol (AUC=", round(auc(roc_arbol_list[[clase]]),3), ")"),
           paste("SVM Lineal (AUC=", round(auc(roc_svm_lineal_list[[clase]]),3), ")"),
           paste("SVM Radial (AUC=", round(auc(roc_svm_k_radial_list[[clase]]),3), ")"),
           paste("SVM Polinomial (AUC=", round(auc(roc_svm_k_polinomial_list[[clase]]),3), ")")
         ),
         col=colores,
         lwd=2
  )
  
  abline(0,1, lty=2, col="darkgrey")   # Línea diagonal de referencia
}
