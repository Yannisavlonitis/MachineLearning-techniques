# -------------- Metodos de aprendizaje supervisado: --------------------------

#                               
# _________________________________ LDA ________________________________________
# En el dataset hay dos problemas para poder hacer las formulas del lda y el modelo
#   1. porque hay variables que son constantes para todos los grupos y da problemas
#   2. porque los nombres de las muestras contienen guiones y seran problematicas
#   a la hora de crear la formula que evalua el lda. 
# Por ello buscamos las variables que tienen varianza != 0 y son las unicas que
# incluimos en el analisis (497 variables frente a 500).

con_varianza <- apply(trainData[, -501], 2, var) != 0 #buscamos var != 0
train_lda <- trainData[, con_varianza] # nos quedamos solo con esas variables
colnames(train_lda) <- make.names(colnames(train_lda)) #para quitar los guiones de los nombres
train_lda$Clase <- trainData$Clase #metemos las clases, tenemos que mantener el mismo orden que en trainData
# Generamos la formula para el lda
formula <- as.formula(paste("Clase ~", paste(colnames(train_lda[, -498]), collapse = "+")))
formula # Var dependiente ~ vars independientes

library(MASS)

# Ajustar el modelo LDA
set.seed(42)

lda_model <- lda(formula, data = train_lda)
lda_model$scaling # contribuciones al modelo de cada variable independiente

# Realizar predicciones
colnames(testData) <- make.names(colnames(testData))
lda_predictions <- predict(lda_model, newdata = testData)
lda_predictions$x

# Obtener la predicción (predicciones de la clase)
predicted_classes <- lda_predictions$class 
# Obtener las verdaderas etiquetas (las clases reales en el conjunto de prueba)
true_classes <- as.factor(testData$Clase)

# Crear la matriz de confusión
conf_lda <- confusionMatrix(predicted_classes, true_classes)
conf_lda
# Obtener probabilidades
probabilities_lda <- predict(lda_model, newdata = testData, type = "prob")
probabilities_lda



## ----------------------------------------------------------
## -------------------------- QDA ---------------------------
## ----------------------------------------------------------
set.seed(42)
qda_model <- qda(formula, data = train_lda)
qda_model



qda_predictions <- predict(qda_model, newdata = testing_extension)

## Clases predichas por QDA
predicted_classes <- qda_predictions$class
predicted_classes
length(predicted_classes)

## Clases reales del conjunto de prueba
true_classes <- as.factor(testing_extension$extension)
true_classes
length(true_classes)


## ----------------------------------------------------------
## --------------------- MATRIZ DE CONFUSIÓN ----------------
## ----------------------------------------------------------

confusion <- confusionMatrix(predicted_classes, true_classes)

#### ---- RDA (regularizado) ----

set.seed(42)
library(klaR)
rda_model <- rda(formula, data = train_lda)
rda_model$scaling # contribuciones al modelo de cada variable independiente

# Realizar predicciones
rda_predictions <- predict(rda_model, newdata = testData)
rda_predictions$posterior

# Obtener la predicción (predicciones de la clase)
predicted_classes <- rda_predictions$class 
# Obtener las verdaderas etiquetas (las clases reales en el conjunto de prueba)
true_classes <- as.factor(testData$Clase)



# Crear la matriz de confusión
conf_rda <- confusionMatrix(predicted_classes, true_classes)
conf_rda
# Obtener probabilidades
probabilities_rda <- predict(rda_model, newdata = testData, type = "prob")
probabilities_rda