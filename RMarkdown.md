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
