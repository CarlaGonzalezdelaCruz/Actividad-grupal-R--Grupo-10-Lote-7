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
Dataset_expresión_genes <- read_csv("Master en Bioinformática/Estadística y R para Ciencias de la Salud/Actividad grupal/Dataset expresión genes.csv")

library(tidyverse)
library(factoextra)
library(gtsummary)
library(gt)
library(car) 
library(broom)
```

Antes de comenzar con el análisis, confirmamos la ausencia de valores "NA" en el dataset.

```{r}
any(is.na(Dataset_expresión_genes))
```

#### explicación del paciente eliminado

```{r}
Dataset_expresión_genes <- Dataset_expresión_genes %>%  filter(id != 14)
```


## Análisis de Componentes Principales (PCA)

```{r expresión génica} 
datos_expresiongenica <- Dataset_expresión_genes %>%
  select(starts_with("AQ_"))
```


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

```{r pca} 
pca<- prcomp(datos_expresiongenica, scale=TRUE)
summary(pca)
```

```{r eigenvalues}
eigenvalues <- get_eigenvalue(pca)
eigenvalues_formatted <- eigenvalues
eigenvalues_formatted[] <- lapply(eigenvalues_formatted, function(x) sprintf("%.3f", as.numeric(x)))
print(eigenvalues_formatted) 
```

#### elegimos 5 componentes

```{r} 
##Tabla PCA cargas
# Obtener las cargas de las variables en los componentes principales
cargas <- pca$rotation  # Las cargas están en la matriz 'rotation'

# Seleccionar solo los primeros n componentes (en este caso, 5 componentes)
cargas_seleccionadas <- as.data.frame(cargas[, 1:5])

# Añadir los nombres de las variables como una columna
cargas_seleccionadas <- tibble::rownames_to_column(cargas_seleccionadas, var = "Variable")

# Renombrar las columnas para mayor claridad
colnames(cargas_seleccionadas) <- c("Variable", "Componente 1", "Componente 2", 
                                    "Componente 3", "Componente 4", "Componente 5")

# Formatear las cargas con 3 decimales
cargas_seleccionadas <- cargas_seleccionadas %>%
  mutate(across(-Variable, ~ sprintf("%.3f", .)))

# Mostrar la tabla
cargas_seleccionadas
```

## Gráficos descriptivos de los componentes principales


#### Tenendo en cuenta los graficos de barra y las funciones de los genes vamos nombrando a los componentes:

Componente 1:Inflamación Sistémica y Señalización Celular

Componente 2: Regulación Inmune y Estrés Oxidativo

Componente 3: Metabolismo Celular y Resistencia al Estrés

Componente 4: Homeostasis Metabólica y Respuesta Antioxidante

Componente 5: Inflamación y Diferenciación Inmune


## Expresión génica en función de las cargas de los pacientes

```{r} 
pca_ind <- as.data.frame(pca$x) # Extraer los valores de los componentes principales (scores) para los individuos

# Dividir en terciles para cada componente
terciles_componente_1 <- quantile(pca_ind$PC1, probs = c(0, 1/3, 2/3, 1))
terciles_componente_2 <- quantile(pca_ind$PC2, probs = c(0, 1/3, 2/3, 1))
terciles_componente_3 <- quantile(pca_ind$PC3, probs = c(0, 1/3, 2/3, 1))
terciles_componente_4 <- quantile(pca_ind$PC4, probs = c(0, 1/3, 2/3, 1))
terciles_componente_5 <- quantile(pca_ind$PC5, probs = c(0, 1/3, 2/3, 1))

# Asignar a cada muestra su tercil correspondiente
# Creo una nueva tabla que dice para cada paciente y componente, a cual tercil corresponde el valor del score
pca_terciles <- data.frame(
  row = rownames(pca_ind),
  Componente_1 = cut(pca_ind$PC1, breaks = terciles_componente_1, labels = c("T1", "T2", "T3"), include.lowest = TRUE),
  Componente_2 = cut(pca_ind$PC2, breaks = terciles_componente_2, labels = c("T1", "T2", "T3"), include.lowest = TRUE),
  Componente_3 = cut(pca_ind$PC3, breaks = terciles_componente_3, labels = c("T1", "T2", "T3"), include.lowest = TRUE),
  Componente_4 = cut(pca_ind$PC4, breaks = terciles_componente_4, labels = c("T1", "T2", "T3"), include.lowest = TRUE),
  Componente_5 = cut(pca_ind$PC5, breaks = terciles_componente_5, labels = c("T1", "T2", "T3"), include.lowest = TRUE)
)

pca_terciles
```

```{r}
# Creo primero una tabla que combina solo los datos que necesito, osea genes y tercil (e id)
# La convierto a long para poder aplicar luego strata

Dataset_expresión_genes <- Dataset_expresión_genes %>%
  mutate(row = row_number())

expgenica_terciles <- Dataset_expresión_genes %>%
  select(row, starts_with("AQ_")) %>%
  mutate(row = as.numeric(row)) %>%
  left_join(pca_terciles %>% mutate(row = as.numeric(row)), by = "row")

colnames(expgenica_terciles)
 ```

```{r genes}
genes <- names(Dataset_expresión_genes)[startsWith(names(Dataset_expresión_genes), "AQ")]
```

```{r, message=FALSE, warning=FALSE}
expgenica_CP1 <- select(expgenica_terciles, starts_with("AQ_"), Componente_1)

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

CP1_hom <- levene_CP1 %>% # Usamos ANOVA
  filter(Homogeneidad == "Sí") %>%
  pull(Gen)

CP1_no_hom <- levene_CP1 %>% # Usamos K-W
  filter(Homogeneidad == "No") %>%
  pull(Gen)  

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

```{r, message=FALSE, warning=FALSE}
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
        pvalue_fun = ~ style_pvalue(.x, digits = 3)) %>%
  modify_caption("**Componente 2. Regulación Inmune y Estrés Oxidativo**")

resumen_CP2
```

```{r, message=FALSE, warning=FALSE}
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

```{r, message=FALSE, warning=FALSE}
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

```{r, message=FALSE, warning=FALSE}
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
        modify_caption("**Componente 5**") %>%
        modify_caption("**Componente 5. Inflamación y Diferenciación Inmune**")

resumen_CP5
```

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
```

## Implementar un modelo de regresión logística

```{r} 
unique(Dataset_expresión_genes$extension)
```

```{r}
dataset_expresión_terciles <- Dataset_expresión_genes %>% #creo un nuevo dataframe con toda la info que puedo llegar a necestiar
  mutate(row = as.numeric(row)) %>%
  left_join(pca_terciles %>% mutate(row = as.numeric(row)), by = "row") %>% #Uno la df original con los terciles
  mutate(metastasis = factor(if_else(extension == "metastasico", "sí", "no"))) %>% #creo la columna metastaris para ralizar la regresion logistica
  select(-c(1, 2, 57, 105)) #elimino las columnas ...1, id, row (que cree en el codigo del 3 para unir las df) y extension (ya que ahora tengo la de metastasis)

modelo_logistica_componentes <- glm(metastasis ~ Componente_1 + Componente_2 + Componente_3 + Componente_4 + Componente_5, 
                        data= dataset_expresión_terciles, family = "binomial")

# Convertir coeficientes en OR (exponenciar)
resultado_logistica_componentes <- tidy(modelo_logistica_componentes, exponentiate = TRUE, conf.int = TRUE) 

modelo_logistica_componentes
```

```{r}
resultado_logistica_componentes %>%
  select(term, estimate, conf.low, conf.high, p.value) %>%
  mutate(OR_IC = sprintf("%.1f (%.1f - %.1f)", estimate, conf.low, conf.high),
         p.value = format.pval(p.value, digits = 3, eps = 0.001)) %>%
  select(term, OR_IC, p.value) %>%
  gt()
```

```{r}
sintomas <- c('tos', 'disnea', 'expect', 'secrecion', 'dolor_garg', 'escalofrios', 'fiebre', 
              'cansancio', 'cefalea', 'mareo', 'nauseas', 'vomitos', 'diarrea', 'dolor_hueso', 
              'dolor_abdo', 'perd_ape', 'disgueusia')

dataset_expresión_terciles_sintomas <- dataset_expresión_terciles %>%
  select (all_of(sintomas)) %>%
  mutate(across(all_of(sintomas), ~ if_else(. == "si", 1, if_else(. == "no", 0, NA_real_)))) #convierto a 0 y 1 para que no salgan raros los nombres de la primera fila

modelo_logistica_sintomas <- glm(metastasis ~ ., 
                                   data = dataset_expresión_terciles_sintomas %>% select(metastasis, all_of(sintomas)), 
                                   family = "binomial")

resultado_logistica_sintomas <- tidy(modelo_logistica_sintomas, exponentiate = TRUE, conf.int = TRUE) 

resultado_logistica_sintomas %>%
  select(term, estimate, conf.low, conf.high, p.value) %>%
  mutate(OR_IC = sprintf("%.1f (%.1f - %.1f)", estimate, conf.low, conf.high),
         p.value = format.pval(p.value, digits = 2, eps = 0.001)) %>%
  select(term, OR_IC, p.value) %>%
  gt()

# la presencia de dolor abdominar y disnea son un riesgo estadisticamente significativo
```

```{r}
variables_incluidas <- c('metastasis', 'Componente_1', 'Componente_2', 'Componente_3', 'Componente_4', 
                         'Componente_5', 'dolor_abdo', 'disnea', 'sexo', 'exfumador')

modelo_logistica_resumen <- glm(metastasis ~ ., 
                                 data = dataset_expresión_terciles %>% select(all_of(variables_incluidas)), 
                                 family = "binomial")

resultado_logistica_resumen <- tidy(modelo_logistica_resumen, exponentiate = TRUE, conf.int = TRUE) 

resultado_logistica_resumen <- resultado_logistica_resumen %>%
  filter(term != "(Intercept)") %>%
  mutate(term = factor(term, levels = rev(term))) # Ordenar en el gráfico

# Crear Dot Plot con intervalos de confianza
grafico <- ggplot(resultado_logistica_resumen, aes(x = estimate, y = term)) +
  geom_point(size = 3, color = "blue") +  # Punto para cada OR
  geom_errorbarh(aes(xmin = conf.low, xmax = conf.high), height = 0.2, color = "black") + # Intervalo de confianza
  geom_vline(xintercept = 1, linetype = "dashed", color = "red") +  # Línea de referencia en OR = 1
  labs(title = "Regresión logística: OR e Intervalos de confianza",
       x = "Odds Ratio (OR)",
       y = "Componentes principales",
       caption = "Línea roja discontinua representa OR = 1") +
  theme_minimal()

grafico
```



## Conclusión
