############################## Vir ############################## 
rm(list=ls()) # Limpiar entorno

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




#------ K means -------

# Kmeans es un algoritmo de aprendizaje no supervisado que agrupa observaciones
# similares intentando minimizar la suma de distancias cuadradas de cada punto a
# su centroide. 

# funcion kmeans (argumentos)
#   x: matriz numerica o dataframe
#   centers: numero de clusteres que se desean formar
#   iter.max: numero de iteraciones para el algoritmo
#   nstart: veces que se ejecuta el algoritmo con diferentes centroides (10 - 25)
#   algorithm: varias opciones (Hartingan-wong por defecto)

# Primero hay que ver si existen valores NA, constantes o infinitos porque Kmeans si no,
# se rompe:

anyNA(gene_exp, recursive = TRUE)

any(!is.finite(as.matrix(df_scaled)))

df_num <- gene_exp[, apply(gene_exp, 2, sd) != 0] # quitamos los constantes
df_scaled <- data.frame(scale(df_num)) # escalamos para que todas las variables
# tengan la misma contribucion


library(factoextra) #para visualizar los graficos del kmeans

# k=3
set.seed(42)
kmeans.result <- kmeans(df_scaled, centers = 6, iter.max = 100, nstart = 25)
# Visualizacion
fviz_cluster(kmeans.result, df_scaled, xlab = '', ylab = '', labelsize = 0) +
  ggtitle("Cluster plot", subtitle = "Centers = 6") +
  theme_classic()


# n optimo de clusters
fviz_nbclust(df_scaled, kmeans, method = "wss") +
  ggtitle("Optimal number of clusters", subtitle = "") +
  theme_classic()

# Comprobaciones de calidad del kmeans
kmeans.result$size
kmeans.result$withinss
kmeans.result$betweenss / kmeans.result$totss

# Se puede ver como el grafico de sedimentacion no indica claramente
# un numero de clusters optimo por la regla del codo. Pues no hay un punto de
# inflexion claro, quedandose entre 4-6 el descenso de suma de cuadrados ya algo 
# menos pronunciado que antes, pero no de manera clara...
# Ademas, viendo el grafico de clusteres, se ve como 4 de los 6 clusteres estan
# casi completamente solapados en dos dimemsiones... y que si cambiamos de k, la
# varianza explicada (kmeans.result$betweenss / kmeans.result$totss )
# no sube de 0.35.
# Esto indica que el algoritmo de Kmeans no es el mejor para este dataset.


# ------------------------------ Isomap -----------------------------

library(RDRToolbox)
library(ggplot2)

# Funcion Isomap()
#   data -> datos (matriz) sobre los que haremos reduccion de dimensionalidad
#   dim -> dimensiones de las columnas del espacio reducido
#   k -> numero de vecinos cercanos a cada punto. A mayor k mayor computacion
#   potResiduals -> devuelve la varianza explicada por las diferentes dimensiones

#   Si se ha elegido una unica dimension devuelve una matriz
#   Si se ha elegido un vector de dimensiones devolvera una matriz por cada elemento del vector

# En el caso de Isomap tambien hay que usar las variables escaladas 
# pero deben ser en forma de matriz:
scaled_matrix <- as.matrix(df_scaled)

# Calculamos isomap de 1 a 10 dimensiones para aproximar
isomap.results <- Isomap(data=scaled_matrix, dims=1:10, k=15, plotResiduals=TRUE)

# Dataframe con los puntos que queremos dibujar en el plano 2D
#     (elegiriamos otro si queremos otra dimension)
isomap.df <- data.frame(isomap.results$dim3) 

# Graficamos en 2d
ggplot(isomap.df, aes(x = X1, y = X2, color = metadata$Clase)) +
  geom_point(size = 3) +
  scale_color_manual(values = c("red", "deepskyblue", "lightgreen", "#ffde21", "purple")) +
  labs(title = "Isomap - Types of Cancer", x = 'Dim 1', y = 'Dim 2', color = "Grupo") +
  theme_classic() +
  theme(panel.grid.major = element_line(color = "gray90"), panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "gray95"), plot.title=element_text(hjust=0.5))

# Graficamos en 3d
coords <- isomap.df

library(plotly)
plot_ly(
  x = coords[,1],
  y = coords[,2],
  z = coords[,3],
  color = as.factor(metadata$Clase),
  colors =c("red", "deepskyblue", "lightgreen", "#ffde21", "purple"),
  type = "scatter3d",
  mode = "markers"
)

# Interpretacion isomap
# Al ver el grafico de sedimentacion de isomap indica que las dimensiones donde la
# varianza disminuye menos es a partir de la dimension 3. Por ello se tienen en
# cuenta las 10 primeras dimensiones simplemente para tener un amplio margen de
# calculo. A la hora de graficar grafico en 2d y 3d para ver si se consiguen separar
# mas los grupos de uno a otro pero no consigo separar el verde del morado.
# Ante el cambio de k, he probado con 10 y 20 y numeros por encima y debajo y este k
# es el que mejor consigue separar los grupos, pero podria mejorarse quizas con otro
# algoritmo.