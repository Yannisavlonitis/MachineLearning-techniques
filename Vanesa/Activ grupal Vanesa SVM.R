rm(list=ls()) # Limpiar entorno

setwd("~/AIA/Activ grupal")

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

####SVM  Kernel (Lineal, polinomial y gaussiano)
anyNA(df_final)

#cargar librerias
library(caret)
library(tidyverse)
library(glmnet)
library(dplyr)
library(gridExtra)
library(e1071)

# Dividir el conjunto de datos en conjuntos de entrenamiento y prueba
set.seed(1995)
trainIndex <- createDataPartition(df_final$Clase, p = 0.8, list = FALSE) #80% trainig, 20% teste
df_final$Clase <- as.factor(df_final$Clase)
trainData <- df_final[trainIndex,]
testData <- df_final[-trainIndex,]

###SVM lineal###
svmModelLineal <- train(Clase ~.,
                        data = trainData,
                        method = "svmLinear",
                        trControl = trainControl(method = "cv", number = 10),
                        preProcess = c("center", "scale"),
                        tuneGrid = expand.grid(C = seq(0, 2, length = 30)), #C grande lleva al sobreajuste, C pequeño al infraajuste. Prueba valores de "C" entre 0 y 20 cada 2
                        prob.model = TRUE)

svmModelLineal    #The final value used for the model was C = 0.2758621.
plot(svmModelLineal)

# Predicciones en el conjunto de prueba utilizando el modelo entrenado
predictionsL <- predict(svmModelLineal, newdata = testData )
predictionsL

# Evaluar la precisión del modelo utilizando la matriz de confusión
confusionMatrix(predictionsL, testData$Clase)

#          Reference
#Prediction AGH CFB CGC CHC HPB
#      AGH  29   0   0   0   0
#      CFB   0  60   0   3   0
#      CGC   0   0  28  20   2   poco fiable al predecir CGC (falsos +)
#      CHC   0   0   0   4   0   mal modelo para CHC(falsos-)
#      HPB   0   0   0   0  13

#Overall Statistics

#Accuracy : 0.8428          
#95% CI : (0.7767, 0.8956)
#No Information Rate : 0.3774          
#P-Value [Acc > NIR] : < 2.2e-16       

#Kappa : 0.7903          

#Mcnemar's Test P-Value : NA              

#Statistics by Class:

 #                    Class: AGH Class: CFB Class: CGC Class: CHC Class: HPB
#Sensitivity              1.0000     1.0000     1.0000    0.14815    0.86667
#Specificity              1.0000     0.9697     0.8321    1.00000    1.00000
#Pos Pred Value           1.0000     0.9524     0.5600    1.00000    1.00000
#Neg Pred Value           1.0000     1.0000     1.0000    0.85161    0.98630
#Prevalence               0.1824     0.3774     0.1761    0.16981    0.09434
#Detection Rate           0.1824     0.3774     0.1761    0.02516    0.08176
#Detection Prevalence     0.1824     0.3962     0.3145    0.02516    0.08176
#Balanced Accuracy        1.0000     0.9848     0.9160    0.57407    0.93333

probabilities_svm_Lineal <- predict(svmModelLineal, newdata = testData, type = "prob")
probabilities_svm_Lineal

### SVM Kernel radial###
svmModelKernel <- train(Clase ~.,
                        data = trainData,
                        method = "svmRadial", #metodo radial
                        trControl = trainControl(method = "cv", number = 10),
                        preProcess = c("center", "scale"),
                        tuneLength = 10,
                        prob.model = TRUE) 
svmModelKernel #The final values used for the model were sigma = 0.001241186 and C = 16.

plot(svmModelKernel)

# Predicciones en el conjunto de prueba utilizando el modelo entrenado
predictions <- predict(svmModelKernel, newdata = testData )
predictions

# Evaluar la precisión del modelo utilizando la matriz de confusión
confusionMatrix(predictions, testData$Clase)
#               Reference
#Prediction AGH CFB CGC CHC HPB
      #AGH  29   0   0   0   0
      #CFB   0  60   0   3   0
      #CGC   0   0  28  20   2  #poco fiable al predecir CGC (falsos +)  
      #CHC   0   0   0   4   0  #mal modelo para CHC(falsos-)
      #HPB   0   0   0   0  13

#Overall Statistics

#Accuracy : 0.8428          
#95% CI : (0.7767, 0.8956) # ¿demasiado amplio?
#No Information Rate : 0.3774          
#P-Value [Acc > NIR] : < 2.2e-16    Significativo   

#Kappa : 0.7903          

#Mcnemar's Test P-Value : NA              

#Statistics by Class:

#                     Class: AGH Class: CFB Class: CGC Class: CHC Class: HPB
#Sensitivity              1.0000     1.0000     1.0000    0.14815    0.86667
#Specificity              1.0000     0.9697     0.8321    1.00000    1.00000
#Pos Pred Value           1.0000     0.9524     0.5600    1.00000    1.00000
#Neg Pred Value           1.0000     1.0000     1.0000    0.85161    0.98630
#Prevalence               0.1824     0.3774     0.1761    0.16981    0.09434
#Detection Rate           0.1824     0.3774     0.1761    0.02516    0.08176
#Detection Prevalence     0.1824     0.3962     0.3145    0.02516    0.08176
#Balanced Accuracy        1.0000     0.9848     0.9160    0.57407    0.93333

# Análisis: 
#Modelo muy bueno para detectar AGH, CFB, CGC (alta sensibilidad)
#Muy bueno para descartar: AGH, CFB, CHC, HPB (alta especificidad)
#El modelo está bien entrenado para la clasificacion de: AGH y CFB (prevalencia real y detectada, similares)
# En el caso de CGC tiene un elevado porcentaje de falsos positivos (56% de aciertos al id +)
#En el caso de CHC sucede al revés, elevado porcentaje de falsos negativos (85% de aciertos al id-)
probabilities_svm_kernel <- predict(svmModelKernel, newdata = testData, type = "prob")
probabilities_svm_kernel

###SVM KernelPolinomial###

svmModelKernelPolynomial <- train(Clase ~.,
                                  data = trainData,
                                  method = "svmPoly",
                                  trControl = trainControl(method = "cv", number = 10),
                                  preProcess = c("center", "scale"),
                                  tuneLength = 5,
                                  prob.model = TRUE) 
svmModelKernelPolynomial

plot(svmModelKernelPolynomial)


# Realizar predicciones en el conjunto de prueba utilizando el modelo entrenado
predictionsP <- predict(svmModelKernelPolynomial, newdata = testData )
predictionsP

# Evaluar la precisión del modelo utilizando la matriz de confusión
confusionMatrix(predictionsP, testData$Clase)
#          Reference
#Prediction AGH CFB CGC CHC HPB
      #AGH 29   0   0   0   0
    #   CFB 0  60   1   2   0
    #  CGC  0   0  27  18   4   #poco fiable al predecir CGC (falsos +)
#      CHC  0   0   0   7   0   #mal predictor de CHC (falsos-)
#      HPB  0   0   0   0  11

#Overall Statistics

#Accuracy : 0.8428          
#95% CI : (0.7767, 0.8956)
#No Information Rate : 0.3774          
#P-Value [Acc > NIR] : < 2.2e-16       

#Kappa : 0.7901          

#Mcnemar's Test P-Value : NA              

#Statistics by Class:

 #                    Class: AGH Class: CFB Class: CGC Class: CHC Class: HPB
#Sensitivity              1.0000     1.0000     0.9643    0.25926    0.73333
#Specificity              1.0000     0.9697     0.8321    1.00000    1.00000
#Pos Pred Value           1.0000     0.9524     0.5510    1.00000    1.00000
#Neg Pred Value           1.0000     1.0000     0.9909    0.86842    0.97297
#Prevalence               0.1824     0.3774     0.1761    0.16981    0.09434
#Detection Rate           0.1824     0.3774     0.1698    0.04403    0.06918
#Detection Prevalence     0.1824     0.3962     0.3082    0.04403    0.06918
#Balanced Accuracy        1.0000     0.9848     0.8982    0.62963    0.86667

# SVM kernelpol
probabilities_svm_kernelpol <- predict(svmModelKernelPolynomial, newdata = testData, type = "prob")
probabilities_svm_kernelpol


