# Parches aplicados a iMKT (v0.2 → v0.2.1)

Este documento resume los cambios aplicados sobre el código fuente descargado
de https://github.com/BGD-UAB/iMKT (rama master), a partir de la sesión de
depuración del 25-26 de septiembre de 2026. Todos los cambios están
confirmados contra el comportamiento ya validado en producción, salvo donde
se indica lo contrario.

## Bugs corregidos

### 1. `R/PopFlyAnalysis.R` — validación de `genes`/`pops` incompleta
Las comprobaciones `genes == ''` y `pops == ''` no usaban `all()`, por lo que
con más de un gen/población fallaban con un error de coerción de R
("'length = N' in coercion to logical(1)") en vez de validar correctamente.
Corregido a `all(genes == '')` / `all(pops == '')`.

### 2. `R/PopHumanAnalysis.R` — mismo bug que el punto 1
Idéntico fix aplicado.

### 3. `R/PopHumanAnalysis.R` — columna de emparejamiento de genes incorrecta
El código original filtraba y agrupaba genes usando `data$globalID`
(identificadores tipo Ensembl, `ENSG...`). El frontend de producción
(formularios web, ficheros de listas de genes) siempre ha enviado
**símbolos de gen** (p. ej. `TSPAN6`, `SCYL3`), no IDs de Ensembl. Se ha
sustituido `globalID` por `symbol` en todo el fichero (filtro inicial,
`subsetGenes`, y los bucles de agrupamiento en ambos bloques `recomb`).

### 4. `R/PopFlyAnalysis.R` — variable inexistente `cutoffs` en el bloque `recomb==TRUE`
En las llamadas a `FWW()`/`imputedMKT()` dentro del bloque `recomb==TRUE`, el
código usaba `listCutoffs=cutoffs` (plural) — pero el parámetro real de la
función es `cutoff` (singular). Esto habría producido "object 'cutoffs' not
found" en cualquier análisis con `recomb=TRUE` y `test="FWW"` o
`test="imputedMKT"`. Corregido a `cutoff` en las 6 líneas afectadas.

### 5. `R/completeMKT.R` (`completimputedMKT()`) — llamaba a funciones inexistentes
El wrapper llamaba a `DGRP()` y `iMKT()`, dos nombres de función que ya no
existen en el paquete desde el renombrado de 2021 (`DGRP`→`imputedMKT`,
`iMKT`→`aMKT`/`asymptoticMKT`). Además llamaba a `FWW(daf, divergence)` sin
el argumento obligatorio `listCutoffs` (sin valor por defecto en `FWW()`),
lo que también habría fallado. Reescrito para:
- Aceptar un nuevo parámetro `listCutoffs` (por defecto `c(0, 0.05, 0.1)`)
- Llamar a `standardMKT`, `FWW`, `eMKT`, `imputedMKT` y `aMKT` con argumentos
  válidos

## Función añadida

### `R/eMKT.R` — no existía en el repositorio oficial
Confirmado en esta sesión que `eMKT()` es una función que no forma parte del
código fuente publicado en GitHub — es una adición local, ya validada
numéricamente en producción (α=0.373 con cutoff=0.15 sobre los datos de
prueba estándar). Se ha incorporado como función exportada del paquete,
usando exactamente la misma fórmula ya validada (basada en Mackay et al.
2012, con el término `f_neutral`).

Se ha conectado como opción `test="eMKT"` en `PopFlyAnalysis()` y
`PopHumanAnalysis()`, en ambos bloques (`recomb=TRUE`/`FALSE`), junto a
`standardMKT`, `FWW`, `imputedMKT` y `aMKT`.

## Confirmado que NO tenía el bug de otra sesión de depuración

**Corrección respecto a la versión 0.2.1 de este documento:** en un primer
análisis se concluyó erróneamente que `R/aMKT.R` y `R/asymptoticMKT.R` ya
eran correctos. Una prueba en sesión de R completamente limpia (`rm(list=ls())`)
demostró que **sí tenían el mismo bug de scoping** ya documentado hace
tiempo por el equipo: `fitMKmodel()` en ambos ficheros llamaba a
`nls2::nls2(formula, start=st, ...)` **sin pasar `data=` explícitamente**.
Sin ese argumento, `nls2` puede depender de forma no fiable de qué objetos
existan en el entorno de llamada para resolver `alpha_trimmed`/`f_trimmed`
— exactamente el mismo mecanismo que causaba que `aMKT()` funcionara o no
según qué variables residuales hubiera en `.GlobalEnv` de una sesión de
depuración anterior.

## Bug adicional corregido (v0.2.2)

### 6 y 7. `R/aMKT.R` y `R/asymptoticMKT.R` — `nls2()` sin `data=` explícito
En ambos ficheros, `fitMKmodel()` construye ahora un `data.frame(alpha_trimmed=..., f_trimmed=...)`
explícito y lo pasa como `data=` en las dos llamadas a `nls2()` (ajuste
inicial y refinamiento). Esto elimina la dependencia del entorno de
llamada, haciendo el ajuste determinista y reproducible en sesión limpia.
Confirmado con `myDafData`/`myDivergenceData`: `aMKT()` ahora converge de
forma consistente (α asintótico ≈ 0.657, IC [0.633, 0.673]) tanto en
sesión limpia como con objetos residuales en el entorno.

### 8. `R/loadPopFly.R` — dataset de PopFly desactualizado
La URL original (`GenesData_recomb_comeron.tab`) solo contiene 2
poblaciones de las 16 que el resto del código espera poder usar (la lista
completa está hardcodeada en `PopFlyAnalysis.R`: AM, AUS, CHB, EA, EF, EG,
ENA, EQA, FR, RAL, SA, SD, SP, USI, USW, ZI). Se ha actualizado a
`GenesData_recomb_comeron_new.tab`, que sí trae las 16 poblaciones
completas (mismos 13753 genes, mismas 12 columnas — solo más filas por
población).

## Pendiente / no verificable sin R

Este parche se ha aplicado únicamente editando texto (sin `R CMD build` ni
`R CMD check`, por no disponer de un entorno R en el momento de editar).
**Antes de usarlo en producción, instálalo y ejecuta al menos una prueba
manual de cada test (`standardMKT`, `FWW`, `eMKT`, `imputedMKT`, `aMKT`) en
los tres pipelines (daf/div directo, PopFly, PopHuman)**, comparando los
resultados con los ya validados hoy en el servidor.

## Bug grave corregido (v0.2.4) — pérdida silenciosa de genes al binear por recombinación

### 9 y 10. `R/PopFlyAnalysis.R` y `R/PopHumanAnalysis.R` — algoritmo de binning defectuoso

Con `recomb=TRUE`, ambos ficheros repartían los genes en `bins` grupos
usando un bucle manual (`for (i in 0:nrow(x)) if (i%%binsize==0) ...`)
con `binsize = round(nrow(x)/bins)`. **Si `nrow(x)` no era exactamente
divisible por `bins`, el último bin se perdía por completo**, no solo
unos pocos genes sobrantes.

Confirmado empíricamente: con 2822 genes y `bins=3`,
`binsize = round(2822/3) = 941`, y `941 × 3 = 2823 ≠ 2822` — el índice
del tercer bin (`i1 = 2823`) supera `nrow(x)`, así que la condición
`i1 <= nrow(x)` nunca se cumple para ese tramo. Resultado: **940 genes
(un tercio del dataset) excluidos silenciosamente**, con solo un aviso
genérico de "excluidos para igualar tamaños de bin" que no refleja la
magnitud real del problema. Con `bins=2` sobre un total par (2822/2=1411
exacto) el bug no se manifestaba, lo que hizo parecer que el código
funcionaba bien — solo aparecía con divisiones no exactas.

**Corrección:** sustituido el bucle manual por
`cut(seq_len(nrow(x)), breaks=bins, labels=FALSE)`, que reparte
**todos** los genes en exactamente `bins` grupos sin excepción,
independientemente de si la división es exacta. El aviso de "genes
excluidos" (que ahora solo dispara por genes no encontrados en los
datos, su propósito original) se mantiene sin cambios.

## Cosmético (v0.2.3) — warnings de ggplot2

`aes_string()` está deprecada desde ggplot2 3.0.0 y `scale_y_discrete()`
no es apropiada para un eje con límites continuos (0 a 1). Ninguno de los
dos afectaba a los resultados numéricos, solo generaban warnings en pantalla.

- **`R/eMKT.R`, `R/imputedMKT.R`, `R/aMKT.R`** (`mutFractionsPlot`,
  `plotDaf`): sustituido `aes_string(x='col', y='col2', ...)` por
  `aes(x=col, y=col2, ...)` sin comillas, ya que en todos estos casos los
  nombres eran literales (columnas fijas del `data.frame` melteado:
  `test`, `value`, `variable`, `daf`), no variables dinámicas.
- **`R/aMKT.R`** (`plotAsymptotic`): el caso `aes_string(x='daf', y=alpha)`
  sí usa un nombre de columna dinámico (`alpha` es un parámetro de la
  función). Sustituido por la forma moderna `aes(x=daf, y=.data[[alpha]])`.
- **`scale_y_discrete(limit=seq(0,1,0.25), expand=c(0,0))`** → 
  `scale_y_continuous(limits=c(0,1), breaks=seq(0,1,0.25), expand=c(0,0))`
  en los tres ficheros — el eje de fracciones (0 a 1) es continuo, no
  discreto; esto es lo que generaba el aviso "Continuous limits supplied
  to discrete scale".

Quedan sin tocar dos apariciones de `aes_string`/`scale_y_discrete` dentro
de bloques de código ya comentado (`#`) en `aMKT.R` — no se ejecutan.

## Viñetas actualizadas (v0.2.5)

### `vignettes/iMKTPipeline.Rmd`
- Corregido typo en la cabecera YAML (`/---` en vez de `---` en la
  primera línea — podía romper el renderizado)
- Renombrados los últimos restos de nomenclatura antigua: la sección final
  (que describe el método asintótico de Messer & Petrov) pasó de
  `### iMKT` a `### aMKT`, coherente con el nombre exportado actual
- Añadida una sección `eMKT correction` completa (paralela a `imputedMKT`),
  con nota explícita señalando que `imputedMKT` y `FWW` dan el mismo α
  (identidad algebraica) mientras que `eMKT` da un α genuinamente distinto
- Corregidos los enlaces de instalación (`devtools::install_github(...)`
  apuntaba a un fork personal desactualizado, `sergihervas/iMKT`; ahora
  apunta a `BGD-UAB/iMKT`)
- **Genealogía de nombres aclarada** (confirmado por Marta Coronado-Zamora,
  coautora de Murga-Moreno et al. 2022 y con conocimiento directo del
  desarrollo del paquete): la corrección de Mackay et al. 2012 se
  implementó originalmente como `DGRP()`, cuya fórmula (`f_neutral`, sobre
  el rango completo de frecuencias) es la que hoy se exporta como
  `eMKT()`. Una mejora posterior y distinta sobre FWW se llamó
  originalmente `iMKT()`, cuya fórmula (`P0Minus/P0Greater`, equivalente
  algebraicamente a FWW) es la que hoy se exporta como `imputedMKT()`. En
  algún punto el nombre "iMKT" se reutilizó también para el método
  asintótico de Messer & Petrov (hoy `aMKT()`) — esa colisión de nombres
  entre dos métodos no relacionados fue un error de documentación de
  Jesús Murga (mantenedor original), no un cambio en las fórmulas de
  ninguno de los dos métodos.

### `vignettes/PopDataPipeline.Rmd`
- Renombrados los ejemplos y prosa que aún decían "DGRP"/"iMKT" para
  que coincidan con el código de los chunks (que ya usaban
  `test="imputedMKT"`/`test="aMKT"` correctamente)
- Corregido el recuento de genes de PopFly (13,745 → 13,753) y de
  PopHuman (18,145 → 20,643, tomado de la propia web de iMKT)
- Añadida mención de `eMKT` como alternativa a `imputedMKT` en el
  ejemplo 1

## Corregido (v0.2.6) — `devtools::document()` no podía regenerar `eMKT.Rd`/`completimputedMKT.Rd`

Ambos `.Rd` los escribí a mano sin la marca `% Generated by roxygen2: do
not edit by hand`, así que roxygen2 los trataba como ficheros "de autor" y
los saltaba en vez de regenerarlos. Se han eliminado del paquete — al
correr `devtools::document()`, roxygen2 los recrea limpiamente a partir
de los comentarios `#'` ya presentes en `R/eMKT.R` y `R/completeMKT.R`,
igual que el resto de `.Rd` del paquete.

## Bug grave corregido (v0.2.7) — `aMKT()` solo devolvía uno de los 3 gráficos

### 11. `R/aMKT.R` — `output$graphs` guardaba `plotAlpha` en vez de `plotsiMKT`

Con `plot=TRUE`, la función construye correctamente los tres paneles
(`dafPoints`=A, `plotAlpha`=B, `plotFraction`=C) y los combina con
`plot_grid(...)` en el objeto `plotsiMKT` — pero justo después, la línea
que guarda el resultado final asignaba `plotAlpha` (solo el panel B, la
curva de ajuste exponencial) al campo `$graphs` del output, en vez de
`plotsiMKT` (el grid combinado). El resultado: quien llamaba a `aMKT()`
solo veía **un** gráfico (el central), nunca la distribución de
frecuencias alélicas (panel A) ni las fracciones de selección negativa
(panel C) — a pesar de que la documentación de la función (y la propia
viñeta) siempre prometió los tres. `plotsiMKT` se calculaba y se
descartaba sin usar. Corregido para que `$graphs` contenga `plotsiMKT`.

Diversos typos corregidos en los manuales de ayuda de las funciones


## Instalación

```r
install.packages("iMKT_0.2.6.tar.gz", repos = NULL, type = "source")
library(iMKT)
```

Dependencias necesarias (ya declaradas en `DESCRIPTION`, se instalan solas
si tienes acceso a CRAN):
```r
install.packages(c("ggplot2", "cowplot", "reshape2", "nls2", "MASS",
                    "ggthemes", "knitr"))
```
