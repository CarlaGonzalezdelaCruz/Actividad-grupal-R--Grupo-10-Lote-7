---
title: "Resolución Actividad 3 máster Bioinformática UNIR"
author: "Valeria Pesa"
date: "2025-01-27"
output: html_document
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
```

# Actividad 3. Análisis de un caso práctico en R

## Introducción

En primer lugar, importamos el dataset e instalamos y cargamos los paquetes apropiados:
* tidyverse para asistir en el análisis de los datos
* factoextra para el análisis de componentes principales
* gtsummary, gt y broom para la creacion de tablas
* car para realizar el test de homogeneidad de varianzes de Levene

```{r, include=FALSE}
library(readr)
Dataset_expresión_genes <- read_csv(file.choose())

library(tidyverse)
library(factoextra)
library(gtsummary)
library(gt)
library(car) 
library(broom)
library(plotly)
```

Antes de comenzar con el análisis, confirmamos la ausencia de valores "NA" en el dataset.

```{r}
any(is.na(Dataset_expresión_genes))
```

Al realizar un paneo de los datos, observamos que el paciente con id 14 presenta valores anómalos en la expresión génica: para los genes AQ_ADIPOQ y AQ_NOX5 presenta una expresión anormalmente alta de 1, mientras que en el resto de los genes la expresión es 0. Consideramos que se trata de un error o datos faltantes, y eliminamos al paciente del análisis. 

```{r}
Dataset_expresión_genes <- Dataset_expresión_genes %>%  filter(id != 14)
```

## Análisis de Componentes Principales (PCA)

Para poder proceder con el análisis de componentes principales, creamos primero un dataframe que incluya solo los datos a analizar: aquellos de expresión génica.

```{r expresión génica} 
datos_expresiongenica <- Dataset_expresión_genes %>%
  select(starts_with("AQ_"))
```

Ya que estaremos realizando análisis estadísticos de los datos de expresión génica, es importante determinar si su distribución es normal o no. para esto, utilizamos el test de Shapiro.

```{r shapiro} 
resultados_saphiro <- apply(datos_expresiongenica, 2, function(x) shapiro.test(x)$p.value)
resultados_saphiro

#creo una tabla con los datos de normalidad
tabla_resultados <- data.frame(
  Variantes = colnames(datos_expresiongenica),
  p_valor = format(resultados_saphiro, digits = 3, scientific = TRUE), 
  Test = rep("Shapiro-Wilk", length(resultados_saphiro)), 
  Interpretacion = ifelse(resultados_saphiro < 0.05, "No normal", "Normal")
)
tabla_resultados
```

Como puede observarse en la tabla, el p-valor para todas las expresiones génicas es menor a 0.05, por lo que se concluye que ningún gen mostró distribución normal según el test de Shapiro. Por lo tanto, los siguientes análisis se basarán en medidas no paramétricas, utilizando la mediana y el rango intercuartílico como estadísticos descriptivos.

Finalmente, aplicamos la función pca() para generar componentes a partir de las variables originales (la expresión de 46 genes). Se observa a continuación un resumen de la data generada: 

```{r pca} 
pca<- prcomp(datos_expresiongenica, scale=TRUE)
summary(pca)
```
Para determinar el número de componentes a seleccionar, calculamos y analizamos los eigenvalues y la varianza explicada acumulativa.

```{r eigenvalues}
eigenvalues <- get_eigenvalue(pca)
eigenvalues_formatted <- eigenvalues
eigenvalues_formatted[] <- lapply(eigenvalues_formatted, function(x) sprintf("%.3f", as.numeric(x)))
print(eigenvalues_formatted) 
```

Decidimos proceder con las primeras 5 dimensiones, ya que en su conjunto explican el 70% de la varianza de los datos.

A continuación, pueden observarse los scores o cargas de cada gen y componente:

```{r} 
# Obtener las cargas de las variables en los componentes principales
cargas <- pca$rotation

# Seleccionar solo los primeros 5 componentes
cargas_seleccionadas <- as.data.frame(cargas[, 1:5]) 

# Añadir los nombres de las variables como una columna
cargas_seleccionadas <- tibble::rownames_to_column(cargas_seleccionadas, var = "Variable") 

# Renombrar las columnas para mayor clardad
colnames(cargas_seleccionadas) <- c("Variable", "Componente 1", "Componente 2", 
                                    "Componente 3", "Componente 4", "Componente 5")

#Formatear las cargas con 3 decimales
cargas_seleccionadas <- cargas_seleccionadas %>%
  mutate(across(-Variable, ~ sprintf("%.3f", .)))

cargas_seleccionadas
```

## Gráficos descriptivos de los componentes principales
Realizamos un gráfico de correlación de las variables en las dos primeras dimensiones.
```{r}
# Gráfico de correlación de las varibles en las dimensiones 1-2
fviz_pca_var(pca, 
             geom = "point",  # Usamos puntos (círculos)
             col.var = "cos2",  
             gradient.col = c("blue", "yellow", "red"),  
             pointsize = 5,  
             scale = "contrib",  
             repel = TRUE)

# Representación de las 10 variables (genes) con más contribución
fviz_pca_var(pca,
             col.var="black",
             select.var =list(contrib =10),
             repel = TRUE)
```
En el gráfico de correlación se puede observar la distribución de las variables en las dos primeras dimensiones del PCA. Las variables en rojo presentan una alta contribución a los dos primeros componentes, lo que indica que estas variables están fuertemente representadas por las dimensiones. Por otro lado, las variables en azul tienen una menor contribución, sugiriendo que su influencia en los componentes es más débil. Es importante destacar que en el primer componente todas las variables presentan una relación positiva, indicando que están alineadas de forma similar al componente. En cambio, el segundo componente presenta una relación más compleja con las variables con mayor dispersión y diversidad. Los genes presentan asociaciones tanto positivas como negativas al segundo componente.

Gráfico de correlación de las variables en 3 dimensiones(1-2-3)
```{r}
# Extraer las coordenadas de los genes (variables) en los 3 primeros componentes principales
pca_genes <- as.data.frame(pca$rotation[, 1:3])
colnames(pca_genes) <- c("PC1", "PC2", "PC3")

# Extraemos los nombres de los genes (filas de la matriz de rotación)
genes <- rownames(pca_genes)

# Generamos el gráfico 3D
grafico3d <- plot_ly(pca_genes, 
               x = ~PC1, y = ~PC2, z = ~PC3, 
               type = 'scatter3d', 
               mode = 'markers+text',
               text = genes, 
               marker = list(size = 5, color = ~PC1, colorscale = 'blues', opacity = 0.7),
               textposition = 'top center')

grafico3d <- grafico3d %>% layout(title = "Distribución de Genes en PCA (3D)",
                      scene = list(xaxis = list(title = "Componente 1"),
                                   yaxis = list(title = "Componente 2"),
                                   zaxis = list(title = "Componente 3")))
grafico3d
```
La distribución de las variables en los tres primeros componentes principales muestra un grupo de genes como AQ_JAK1, AQ_CCL5 y AQ_NFE2L2 que se agrupan en el centro del gráfico, indicando una similitud en su contribución a los componentes. No obstante, también se diferencian variables más alejadas como AQ_SLC2A4, AQ_ARG1 y AQ_FOXO3, que presentan un comportamiento diferente a su contribución a los componentes respecto a los demás genes. Al igual que en el gráfico anterior, el primer componente mantiene una contribución positiva para todas las variables. En cambio, el segundo y tercer componente muestran una mayor variabilidad y dispersión con asociaciones tanto positivas como negativas entre los genes.

Variables con mayor contribución en las dimensiones. Se ha analizado tanto las dimensiones por separado como la combinación de estas con la dimensión 1, que es la que contribuye más (52'49%). 
```{r}
fviz_contrib(pca, choice = "var", axes = 1, top = 50)
fviz_contrib(pca, choice = "var", axes = 2, top = 50)
fviz_contrib(pca, choice = "var", axes = 3, top = 50)
fviz_contrib(pca, choice = "var", axes = 4, top = 50)
fviz_contrib(pca, choice = "var", axes = 5, top = 50)
```
En el análisis de la Dimensión 1, los genes más influyentes incluyen AQ_JAK1, AQ_CCL5, AQ_SREBF1 y AQ_NFE2L2.   
En la Dimensión 2, las contribuciones están más dispersas, destacándose genes como AQ_ARG1, AQ_CCL1 y AQ_LIF, los cuales presentan valores superiores en comparación con los demás genes.  
En la Dimensión 3, se identifican algunas variables con contribuciones elevadas como  son AQ_SLC2A4 y AQ_JAK3.  
Por otro lado, la Dimensión 4 tienen una gran contribución AQ_NOX5, AQ_ADIPOQ y AQ_NOS2, superando el 30%.  
Finalmente, en la Dimensión 5, los genes con mayor contribución son AQ_IL6, AQ_CSF1, AQ_NOS2, AQ_BMP2 y AQ_SLC2A4.

```{r}
fviz_contrib(pca, choice = "var", axes = c(1, 2), top = 50)
fviz_contrib(pca, choice = "var", axes = c(1, 3), top = 50)
fviz_contrib(pca, choice = "var", axes = c(1, 4), top = 50)
fviz_contrib(pca, choice = "var", axes = c(1, 5), top = 50)
```
En la combinación de las Dimensiones 1 y 2, se observa una distribución más uniforme de las contribuciones, con valores más equilibrados en comparación con las gráficas por individual. Entre los genes con mayores contribuciones se encuentran AQ_JAK1, AQ_PTAR, AQ_SREBF1 y AQ_NFE2L2.  
Al considerar las Dimensiones 1 y 3, se evidencia que algunos genes que eran dominantes en la Dimensión 3, como AQ_SLC2A4, AQ_JAK3 y AQ_FOXO3, ya no presentan un aporte tan marcado. En contraste, genes como AQ_JAK1, AQ_IL5 y AQ_TGFBI adquieren mayor importancia en esta gráfica.   
En la combinación de las Dimensiones 1 y 4, la contribución se encuentra más dispersa y distribuida entre varios genes, con AQ_JAK1, AQ_CCL5, AQ_SREBF1 y AQ_NFE2L2 mostrando valores similares, en torno al 3%.  
Finalmente, en la combinación con la Dimensión 5, los genes con mayor contribución son AQ_JAK1, AQ_CCL5, AQ_SREBF1, AQ_NFE2L2 y AQ_PTARF.

Teniendo en cuenta los valores de los scores y la contribución de cada gen a cada dimensión, así como la función de dichos genes, decidimos nombrar a los componentes de la siguiente manera:

* Componente 1: Inflamación Sistémica y Señalización Celular

* Componente 2: Regulación Inmune y Estrés Oxidativo

* Componente 3: Metabolismo Celular y Resistencia al Estrés

* Componente 4: Homeostasis Metabólica y Respuesta Antioxidante

* Componente 5: Inflamación y Diferenciación Inmune


Gráfico de correlación de los pacientes en clústers en función de su expresión génica.
```{r}
kmeans <- kmeans(datos_expresiongenica, centers = 3)
grupo <- as.factor(kmeans$cluster)

fviz_pca_ind(pca,
             col.ind = grupo,     
             palette = c("lightgreen", "lightcoral", "skyblue"),
             addEllipses = TRUE,
             legend.title = "Cluster")
```
En función de los dos primeros componentes principales que representan la inflamación sistemática y señalización celular, y la regulación inmune y estrés oxidativo, se han agrupado los pacientes en 3 clústers en base a sus perfiles de expresión génica. El análisis de K-means ha identificado las 3 arupaciones del gráfico:

Clúster 1 (verde, círculos): Representa pacientes con un perfil genético caracterizado por una mayor expresión de genes relacionados con la regulación inmune y el estrés oxidativo, esto queda reflejado por el desplazamiento que se oberva hacia valores positivos en la dimensión 1. Por otro lado, presenta una mayor dispersión lo que puede indicar una mayor heterogeneidad en sus perfiles de expresión génica.

Clúster 2 (rojo, triángulos): Representa a un grupo intermedio de pacientes que combina genes relacionados con la inflamación y la regulación inmune. Los individuos de este grupo se superponen en parte con los otros clústeres, indicando una posible transición entre los perfiles genéticos.

Clúster 3 (azul, cuadrados): Representa pacientes con menor activación de procesos inflamatorios y de señalización celular, localizados en la zona negativa de la dimensión 1 y más concentrados en el gráfico en comparación con los dos clústeres anteriores.

## Expresión génica en función de las cargas de los pacientes

Para proceder con el análisis, determinaremos las estadísticas descriptivas de la expresión génica (mediana y rango intercuartílico, ya que su distribución no es normal) luego de dividir a los individuos en tres terciles según su score de CPA. Por ende, en primer lugar calculamos dichos terciles. Generamos una nueva tabla que indica a qué tercil pertenece cada paciente dentro de cada componente.

```{r} 
pca_ind <- as.data.frame(pca$x) # Extraer los valores de los scores para los individuos

# Dividir en terciles para cada componente
terciles_componente_1 <- quantile(pca_ind$PC1, probs = c(0, 1/3, 2/3, 1))
terciles_componente_2 <- quantile(pca_ind$PC2, probs = c(0, 1/3, 2/3, 1))
terciles_componente_3 <- quantile(pca_ind$PC3, probs = c(0, 1/3, 2/3, 1))
terciles_componente_4 <- quantile(pca_ind$PC4, probs = c(0, 1/3, 2/3, 1))
terciles_componente_5 <- quantile(pca_ind$PC5, probs = c(0, 1/3, 2/3, 1))

# Asignar a cada muestra su tercil correspondiente dentro de una nueva tabla 
pca_terciles <- data.frame(
  row = rownames(pca_ind), # Requeriremos esta columna para unir esta tabla a otra
  Componente_1 = cut(pca_ind$PC1, breaks = terciles_componente_1, labels = c("T1", "T2", "T3"), include.lowest = TRUE),
  Componente_2 = cut(pca_ind$PC2, breaks = terciles_componente_2, labels = c("T1", "T2", "T3"), include.lowest = TRUE),
  Componente_3 = cut(pca_ind$PC3, breaks = terciles_componente_3, labels = c("T1", "T2", "T3"), include.lowest = TRUE),
  Componente_4 = cut(pca_ind$PC4, breaks = terciles_componente_4, labels = c("T1", "T2", "T3"), include.lowest = TRUE),
  Componente_5 = cut(pca_ind$PC5, breaks = terciles_componente_5, labels = c("T1", "T2", "T3"), include.lowest = TRUE)
)

pca_terciles
```

Para proceder con el análisis, requeriremos una tabla que resuma tanto la expresión génica por paciente como a qué tercil pertenecen. Se observa a continuación el código para crear dicha tabla, así como un resumen de las columnas presentes en esta:

```{r}
Dataset_expresión_genes <- Dataset_expresión_genes %>%
  mutate(row = row_number()) # Tal como en el código anterior, creamos una columna row para unir esta tabla con la anterior (ya que id no coincide con el número de fila)

expgenica_terciles <- Dataset_expresión_genes %>%
  select(row, starts_with("AQ_")) %>% # Seleccionamos las columnas de expresión génica
  mutate(row = as.numeric(row)) %>%
  left_join(pca_terciles %>% mutate(row = as.numeric(row)), by = "row") # Nos aseguramos de que "row" esté en el mismo formato en ambas tablas y las unimos por dicha columna

colnames(expgenica_terciles)
```

Cómo último paso antes de comenzar a generar las tablas que resuman los estadísticos de los genes, creamos un vector que contiene el nombre de los genes del dataset:

```{r genes}
genes <- names(Dataset_expresión_genes)[startsWith(names(Dataset_expresión_genes), "AQ")]
```

Con el objetivo de poder aplicar todos los tests estadísticos pertinentes, realizamos tablas separadas para cada componente. Este resumen, para cada tercil, la mediana y el rango intercuartílico de la expresión de cada gen, así como el valor p que indica si la diferencia entre los cuartiles es significativa.

Respetando el orden de los componentes, analizamos primero la expresión génica en función de los scores individuales calculados para el componente 1 (Inflamación Sistémica y Señalización Celular). Utilizamos el test de Levene para determinar si, al agrupar a los pacientes por tercil, la expresión génica es paramétrica u homogénea entre los tres grupos. 

```{r}
expgenica_CP1 <- select(expgenica_terciles, starts_with("AQ_"), Componente_1) # Extraer la expresión génica y solo la data sobre los terciles del componente 1

levene_CP1 <- data.frame(
  Gen = character(46),                     
  pvalue = numeric(46)    
)

for (i in 1:46) {
  levene_CP1[i,1] <- genes[i]
  levene <- leveneTest(expgenica_CP1[[genes[i]]], expgenica_CP1$Componente_1)
  levene_CP1[i,2] <- levene$`Pr(>F)`[1]
}

levene_CP1$Homogeneidad <- ifelse(levene_CP1$pvalue > 0.05, "Sí", "No")

print(arrange(levene_CP1, pvalue))
```

Como puede observarse en la tabla, solo 5 genes evidencian una distribución homogénea. Guardamos estos genes en un vector para poder luego aplicar ANOVA, el cual es indicado cuando hay 2 o más grupos independientes y presentan homogeneidad de varianzas. El resto de los genes son guardados en un vector separado, para poder luego aplicar el test de Kruskal-Wallis, que es utilizado cuando hay 2 o más grupos independientes, pero no presentan homogeneidad de varianzas.

```{r}
CP1_hom <- levene_CP1 %>% # Usamos ANOVA
  filter(Homogeneidad == "Sí") %>%
  pull(Gen)

CP1_no_hom <- levene_CP1 %>% # Usamos K-W
  filter(Homogeneidad == "No") %>%
  pull(Gen)  
```

Finalmente, creamos una tabla resumen que contiene las medianas y rangos intercuartílicos de la expresión de cada gen agrupado por tercil, así como el valor p que indica si la diferencia entre dichos grupos es o no significativa.

```{r, message=FALSE, warning=FALSE}
resumen_CP1 <- expgenica_CP1 %>%
  tbl_summary(by = Componente_1,
              statistic = all_continuous() ~ "{median} ({p25} - {p75})",
              digits = all_continuous() ~ function(x) format(x, digits = 2, scientific = TRUE)) %>%
  add_p(test = list(all_of(CP1_hom) ~ "aov",
                    all_of(CP1_no_hom) ~ "kruskal.test"),
        pvalue_fun = ~ style_pvalue(.x, digits = 3)) %>%
  modify_caption("**Componente 1. Inflamación Sistémica y Señalización Celular**")


resumen_CP1
```

Para la gran mayoría de los genes (exceptuando AQ_ADIPOQ, AQ_LIF, AQ_NOX5 y AQ_SLC2A4), se observan p-valores menores a 0.05. Estos indican que se rechaza la hipótesis nula, lo que sugiere que al menos un grupo de pacientes, clasificado según el tercil del score del PCA, tiene una mediana de expresión génica significativamente diferente a las otras. No es de sorprender que, para el componente 1—que explica la mayor varianza de los datos y está relacionado con la inflamación sistémica y señalización celular—se observe una diferencia estadísticamente significativa de la expresión de la gran mayoría de los genes según la magnitud del score del paciente.

Procedemos de la misma manera con el segundo componente (Regulación Inmune y Estrés Oxidativo):

```{r}
expgenica_CP2 <- select(expgenica_terciles, starts_with("AQ_"), Componente_2)

levene_CP2 <- data.frame(
  Gen = character(46),                     
  pvalue = numeric(46)    
)

for (i in 1:46) {
  levene_CP2[i,1] <- genes[i]
  levene <- leveneTest(expgenica_CP2[[genes[i]]], expgenica_CP2$Componente_2)
  levene_CP2[i,2] <- levene$`Pr(>F)`[1]
}

levene_CP2$Homogeneidad <- ifelse(levene_CP2$pvalue > 0.05, "Sí", "No")

print(arrange(levene_CP2, pvalue))
```

Nuevamente, consideramos los resultados del test de Levene para agrupar los genes en función de la homogeneidad de varianzas y aplicar ANOVA o Kruskal-Wallis según corresponda:

```{r, message=FALSE, warning=FALSE}
CP2_hom <- levene_CP2 %>% # Usamos ANOVA
  filter(Homogeneidad == "Sí") %>%
  pull(Gen)

CP2_no_hom <- levene_CP2 %>% # Usamos K-W
  filter(Homogeneidad == "No") %>%
  pull(Gen)  

resumen_CP2 <- expgenica_CP2 %>%
  tbl_summary(by = Componente_2,
              statistic = all_continuous() ~ "{median} ({p25} - {p75})",
              digits = all_continuous() ~ function(x) format(x, digits = 2, scientific = TRUE)) %>%
  add_p(test = list(all_of(CP2_hom) ~ "aov",
                    all_of(CP2_no_hom) ~ "kruskal.test"),
        pvalue_fun = ~ style_pvalue(.x, digits = 3))  %>%
  modify_caption("**Componente 2. Regulación Inmune y Estrés Oxidativo**")

resumen_CP2
```
Como se puede obsevar en el resumen del Componente 2 (Regulación Inmune y Estrés Oxidativo), se ha encontrado significancia en 36 de 46 genes entre ellos AQ_ALOX5, AQ_ARG1, AQ_CCL2, AQ_CCL5 y AQ_CCR5. Al tener estos un p-value inferior a 0'05, aceptaríamos la hipotesis nula que indica que existen diferencias entre las expresiónes génicas entre los terciles del componente.


A continuación, pueden observarse el código y las tablas correspondientes al análisis de los componentes 3, 4 y 5:
```{r}
expgenica_CP3 <- select(expgenica_terciles, starts_with("AQ_"), Componente_3)

levene_CP3 <- data.frame(
  Gen = character(46),                     
  pvalue = numeric(46)    
)

for (i in 1:46) {
  levene_CP3[i,1] <- genes[i]
  levene <- leveneTest(expgenica_CP3[[genes[i]]], expgenica_CP3$Componente_3)
  levene_CP3[i,2] <- levene$`Pr(>F)`[1]
}

levene_CP3$Homogeneidad <- ifelse(levene_CP3$pvalue > 0.05, "Sí", "No")

print(arrange(levene_CP3, pvalue))
```

```{r, message=FALSE, warning=FALSE}
CP3_hom <- levene_CP3 %>% # Usamos ANOVA
  filter(Homogeneidad == "Sí") %>%
  pull(Gen)

CP3_no_hom <- levene_CP3 %>% # Usamos K-W
  filter(Homogeneidad == "No") %>%
  pull(Gen)  

resumen_CP3 <- expgenica_CP3 %>%
  tbl_summary(by = Componente_3,
              statistic = all_continuous() ~ "{median} ({p25} - {p75})",
              digits = all_continuous() ~ function(x) format(x, digits = 2, scientific = TRUE)) %>%
  add_p(test = list(all_of(CP3_hom) ~ "aov",
                    all_of(CP3_no_hom) ~ "kruskal.test"),
        pvalue_fun = ~ style_pvalue(.x, digits = 3)) %>%
  modify_caption("**Componente 3. Metabolismo Celular y Resistencia al Estrés**")

resumen_CP3
```
En referencia al componente 3 (Metabolismo Celular y Resistencia al Estrés), de los 46 genes analizados, únicamente 16 de ellos presenta un valor de p-value >0.05, lo que implica la aceptación de la hipótesis nula. En otras palabras, la mayoría de las medias de expresión génica muestran diferencias estadísticamente significativas entre algunos de los terciles.

```{r}
expgenica_CP4 <- select(expgenica_terciles, starts_with("AQ_"), Componente_4)

levene_CP4 <- data.frame(
  Gen = character(46),                     
  pvalue = numeric(46)    
)

for (i in 1:46) {
  levene_CP4[i,1] <- genes[i]
  levene <- leveneTest(expgenica_CP4[[genes[i]]], expgenica_CP4$Componente_4)
  levene_CP4[i,2] <- levene$`Pr(>F)`[1]
}

levene_CP4$Homogeneidad <- ifelse(levene_CP4$pvalue > 0.05, "Sí", "No")

print(arrange(levene_CP4, pvalue))
```

```{r, message=FALSE, warning=FALSE}
CP4_hom <- levene_CP4 %>% # Usamos ANOVA
  filter(Homogeneidad == "Sí") %>%
  pull(Gen)

CP4_no_hom <- levene_CP4 %>% # Usamos K-W
  filter(Homogeneidad == "No") %>%
  pull(Gen)  

resumen_CP4 <- expgenica_CP4 %>%
  tbl_summary(by = Componente_4,
              statistic = all_continuous() ~ "{median} ({p25} - {p75})",
              digits = all_continuous() ~ function(x) format(x, digits = 2, scientific = TRUE)) %>%
  add_p(test = list(all_of(CP4_hom) ~ "aov",
                    all_of(CP4_no_hom) ~ "kruskal.test"),
        pvalue_fun = ~ style_pvalue(.x, digits = 3)) %>%
  modify_caption("**Componente 4. Homeostasis Metabólica y Respuesta Antioxidante**")

resumen_CP4
```
Como se observa, la mayoría de los genes presentan un valor de p-value <0.05 aceptando la hipótesis alternativa, indicando la presencia de diferencias significativas entre los terciles de expresión génica.

```{r}
expgenica_CP5 <- select(expgenica_terciles, starts_with("AQ_"), Componente_5)

levene_CP5 <- data.frame(
  Gen = character(46),                     
  pvalue = numeric(46)    
)

for (i in 1:46) {
  levene_CP5[i,1] <- genes[i]
  levene <- leveneTest(expgenica_CP5[[genes[i]]], expgenica_CP5$Componente_5)
  levene_CP5[i,2] <- levene$`Pr(>F)`[1]
}

levene_CP5$Homogeneidad <- ifelse(levene_CP5$pvalue > 0.05, "Sí", "No")

print(arrange(levene_CP5, pvalue))

```

```{r, message=FALSE, warning=FALSE}
CP5_hom <- levene_CP5 %>% # Usamos ANOVA
  filter(Homogeneidad == "Sí") %>%
  pull(Gen)

CP5_no_hom <- levene_CP5 %>% # Usamos K-W
  filter(Homogeneidad == "No") %>%
  pull(Gen)  

resumen_CP5 <- expgenica_CP5 %>%
  tbl_summary(by = Componente_5,
              statistic = all_continuous() ~ "{median} ({p25} - {p75})",
              digits = all_continuous() ~ function(x) format(x, digits = 2, scientific = TRUE)) %>%
  add_p(test = list(all_of(CP5_hom) ~ "aov",
                    all_of(CP5_no_hom) ~ "kruskal.test"),
        pvalue_fun = ~ style_pvalue(.x, digits = 3)) %>%
        modify_caption("**Componente 5. Inflamación y Diferenciación Inmune**")

resumen_CP5
```
En los terciles basados en el componente 5 (Inflamación y Diferenciación Inmune), a pesar de que este componente explica solo el 3.82% de la varianza total, la mayoría de genes presentan diferencias estadísticamente significativas. No obstante, cinco genes (AQ_ADIPOQ, AQ_IL10, AQ_LIF, AQ_NOX5 y AQ_TLR3) no muestran diferencias significativas, indicando que su expresión no varía de forma relevante en función de los terciles.


Finalmente, el siguiente código reúne las 5 tablas en una única tabla que resume los estadísticos obtenidos:

```{r}
tabla_combinada <- tbl_merge(
  list(resumen_CP1, resumen_CP2, resumen_CP3, resumen_CP4, resumen_CP5),
  tab_spanner = c("**Componente 1. Inflamación Sistémica y Señalización Celular**", 
                  "**Componente 2. Regulación Inmune y Estrés Oxidativo**", 
                  "**Componente 3. Metabolismo Celular y Resistencia al Estrés**", 
                  "**Componente 4. Homeostasis Metabólica y Respuesta Antioxidante**", 
                  "**Componente 5. Inflamación y Diferenciación Inmune**")  
  ) %>%
  modify_caption("**Expresión génica por componente y tercil**")%>% 
  modify_header(label ~ "**Gen**")

tabla_combinada
```
### Resumen de las tablas de los 5 componentes

## Implementar un modelo de regresión logística

Finalmente, utilizamos un modelo de regresión logística para analizar el riesgo (OR del inglés *odds ratio*) de que un paciente sufra metástasis en funcion de otras variables. 

Para comenzar el análisis, observamos los posibles grados de extensión del tumor en el dataset:

```{r} 
unique(Dataset_expresión_genes$extension)
```

Debemos reducir los tres factores a dos para tener una respuesta binaria y poder aplicar el modelo deseado. Luego, generamos una nueva columna llamada metástasis, que presenta como posibles valores "sí" y "no". Dicha columna está contenida en un nuevo dataframe que contiene las variables del dataframe original, así como la clasificación de los pacientes en terciles según su score por cada componente estudiado.

```{r}
dataset_expresión_terciles <- Dataset_expresión_genes %>% # Creación del un nuevo dataframe
  mutate(row = as.numeric(row)) %>%
  left_join(pca_terciles %>% mutate(row = as.numeric(row)), by = "row") %>%
  mutate(metastasis = factor(if_else(extension == "metastasico", "sí", "no"))) %>% # Creación de la columna metastasis para ralizar la regresion logistica
  select(-c(1, 2, 57, 105)) # Eliminar las columnas innecesarias: ...1, id, row y extension (reemplazada por metastasis)
```
Utilizamos este para aplicar el modelo de regresión logística, seleccionando la columna metastasis como variable dependiente, y a los componentes como variables independientes. Se observan a continuación los resultados del modelo:

```{r}
modelo_logistica_componentes <- glm(metastasis ~ Componente_1 + Componente_2 + Componente_3 + Componente_4 + Componente_5, 
                        data= dataset_expresión_terciles, family = "binomial")

# Convertir coeficientes en OR (exponenciar)
resultado_logistica_componentes <- tidy(modelo_logistica_componentes, exponentiate = TRUE, conf.int = TRUE) %>%
  mutate(term = case_when( # Modificación de las etiquetas
    term == "Componente_1T2" ~ "Componente 1 - T2",
    term == "Componente_1T3" ~ "Componente 1 - T3",
    term == "Componente_2T2" ~ "Componente 2 - T2",
    term == "Componente_2T3" ~ "Componente 2 - T3",
    term == "Componente_3T2" ~ "Componente 3 - T2",
    term == "Componente_3T3" ~ "Componente 3 - T3",
    term == "Componente_4T2" ~ "Componente 4 - T2",
    term == "Componente_4T3" ~ "Componente 4 - T3",
    term == "Componente_5T2" ~ "Componente 5 - T2",
    term == "Componente_5T3" ~ "Componente 5 - T3",
    TRUE ~ term  # Mantiene los términos que no necesitan cambios
  ))

resultado_logistica_componentes
```

Para simplificar su visualización, interpetación y análisis, generamos una nueva tabla que extrae la información pertinente:

```{r, message=FALSE, warning=FALSE}
gt_componentes <- resultado_logistica_componentes %>%
  select(term, estimate, conf.low, conf.high, p.value) %>%
  mutate(OR_IC = sprintf("%.1f (%.1f - %.1f)", estimate, conf.low, conf.high),
         p.value = format.pval(p.value, digits = 3, eps = 0.001)) %>%
  select(term, OR_IC, p.value) %>%
  gt() %>%
  tab_header(title = "Riesgo de metástasis según componente y tercil")

gt_componentes
```

Para cada componente observamos cual es la disminución (OR<1) o el aumento (OR>1) del riesgo de presentar metastasis según si el paciente pertenece al tercil 2 o 3, comparado con el 1. Sin embargo, este riesgo no es estadisticamente significativo. Esto se concluye tanto a partir de los valores p mayores o iguales a 0.05, como el hecho de que los intervalos de confianza (excepto el último) contienen al valor 1 de riesgo nulo.

Sin embargo, dentro del componente 5 (Inflamación y Diferenciación Inmune) podemos afirmar que hay una clara tedencia para los pacientes dentro del tercil 3 a tener menos riesgo de presentar metástasis.

Los hayazgos del modelo de regresión logística para los componentes se encuentran resumidos en el siguiete gráfico:

```{r}
resultado_logistica_componentes_resumen <- resultado_logistica_componentes %>%
  filter(term != "(Intercept)") %>%
  mutate(term = as.character(term)) %>%
  mutate(term = case_when( # Modificación de las etiquetas
    term == "dolor_abdo" ~ "Dolor abdominal",
    term == "disnea" ~ "Disnea",
    TRUE ~ term  # Mantener los valores no modificados
  )) %>%
  mutate(term = factor(term, levels = rev(unique(term))))  # Volver a factor con los nuevos nombres

grafico_OR <- ggplot(resultado_logistica_componentes_resumen, aes(x = estimate, y = term)) +
  geom_point(size = 3, color = "blue") +  # Punto para cada OR
  geom_errorbarh(aes(xmin = conf.low, xmax = conf.high), height = 0.2, color = "black") + # Intervalo de confianza
  geom_vline(xintercept = 1, linetype = "dashed", color = "red") +  # Línea de referencia en OR = 1
  labs(title = "Regresión logística: OR e Intervalos de confianza",
       x = "Odds Ratio (OR)",
       y = "Variable",
       caption = "Línea roja discontinua representa OR = 1") +
  theme_minimal()

grafico_OR
```

Repetimos esta análisis considerando como variables independientes a los parametro bioquímicos, las variables sociodemográficas, las comorbilidades de los pacientes, y los síntomas presentados. Para todas estas variables, exceptuando la presencia de dolor abdominal y disnea, los resultados obtenidos no fueron estadisticamente significativos.

A continuación, se observa el código, tabla, y análisis de los sintomas. En primer lugar creamos un dataset que resuma la presencia o ausencia de metástasis, así como de los síntomas. Para evitar problemas de formato, convertimos los valores de los síntomas "no" y "si" en 0 y 1 respectivamente.

```{r}
sintomas <- c('tos', 'disnea', 'expect', 'secrecion', 'dolor_garg', 'escalofrios', 'fiebre', 
              'cansancio', 'cefalea', 'mareo', 'nauseas', 'vomitos', 'diarrea', 'dolor_hueso', 
              'dolor_abdo', 'perd_ape', 'disgueusia')

dataset_expresión_terciles_sintomas <- dataset_expresión_terciles %>%
  select (metastasis, all_of(sintomas)) %>%
  mutate(across(all_of(sintomas), ~ if_else(. == "si", 1, if_else(. == "no", 0, NA_real_))))

head (dataset_expresión_terciles_sintomas)
```

Luego, aplicamos el modelo de regresión logística y creamos una tabla que resuma los resultados.

```{r, message=FALSE, warning=FALSE}
modelo_logistica_sintomas <- glm(metastasis ~ ., 
                                   data = dataset_expresión_terciles_sintomas %>% select(metastasis, all_of(sintomas)), 
                                   family = "binomial")

resultado_logistica_sintomas <- tidy(modelo_logistica_sintomas, exponentiate = TRUE, conf.int = TRUE) 

gt_sintomas <- resultado_logistica_sintomas %>%
  select(term, estimate, conf.low, conf.high, p.value) %>%
  mutate(OR_IC = sprintf("%.1f (%.1f - %.1f)", estimate, conf.low, conf.high),
         p.value = format.pval(p.value, digits = 2, eps = 0.001)) %>%
  select(term, OR_IC, p.value) %>%
  gt() %>%
  tab_header(title = "Riesgo de metástasis según síntoma")

gt_sintomas
```

Como se adelantó, el p valor e IC del OR superior a 1 para las variables disnea y dolor abdominal sugieren un aumento de riesgo estadisticamente significativo. El OR de la variable disnea indica que en un paciente con disnea el riesgo de presentar metástasis sería unas 300 veces mayor. Sin embargo, el amplio IC hace cuestionable este valor. Dicho esto, sí puede afirmarse con confianza que el riesgo aumenta.

De forma similar, hay un aumento estadisticamente significativo de que un paciente con dolor abdominal presente metástasis. Aunque el OR indica que este aumento es de más de 1700 veces, este valor vuelve a ser cuestionable teniendo en cuenta el IC.

## Conclusión
