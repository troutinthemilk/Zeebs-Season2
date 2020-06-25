---
title: "Module 7: Generalized linear models"
output:
  xaringan::moon_reader:
  css: [default, metropolis, metropolis-fonts, "Title.css"]
lib_dir: libs
includes:
  after_body: logo.html
nature:
  highlightStyle: github
highlightLines: true
countIncrementalSlides: false
editor_options: 
  chunk_output_type: console
header-includes :
  - \usepackage {amsmath}
- \usepackage{mathrsfs}
---
class: inverse
```{r, echo=F, message=F, warning=F, eval=T}
#setting up my ggplot defaults. Update with your own preferences
library(ggplot2)#plotting functions
library(ggthemes) #more themes!
library(wesanderson)
library(RColorBrewer)

theme_set(theme_tufte()) # a theme I like.
theme_update(plot.title = element_text(hjust = 0.5), 
             axis.line.x = element_line(color="black", size = 1),
             axis.line.y = element_line(color="black", size = 1),
             text=element_text(size=20),
             axis.text=element_text(size=15)) #center all titles and and axis lines
```

