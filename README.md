# Habitat selection modeling guidance (in practice)

## A repository for:

Gerber, BD, Setash, C, Ivan, JS and Northrup, JM. A plain language review and guidance for modeling animal habitat-selection. 

---

## What's in this repository?

This repository comprises two main R Markdown files (`TraditionalHSF.Rmd`, `MovementHSF.Rmd`) and associated files to compile them into HTML and markdown files. 
You can download all necessary files in this [Zip](HSF.Guide.Files.zip). If you are interested in only the R code from each Markdown file, these are available in the .R files with the same name. 

## Viewing

The compiled HTMLs can be viewed here:

[Traditional HSF](https://bgerber123.github.io/hsfguide/TraditionalHSF.html)

[Movement-based HSF](https://bgerber123.github.io/hsfguide/MovementHSF.html)

---

## The working directory

### Main files

The two fundamental scripts are `TraditionalHSF.Rmd` and `MovementHSF.Rmd`. Both files are organized 
similarly, setting up the environment and data, and then follow the subsections of the manuscript. 

The first file focuses on addressing the points in the manuscript regarding traditional broad-scale
habitat selection, for example, selection within an individual's home range. The second R markdown file
does most of the same operations but regarding movement-based habitat selection (i.e., step-selection functions).

### Additional Files
- folder `data` has one file, `Covs`, which contain rasters of spatial covariates
- folder `functions` contains three files:
  - `sim.ind.movement.hsf.r`: simulation of individual movement based habibtat selection
  - `mean_ci_boot.r`: bootstrapping function to get population (across individual) inference
  - `sample.size.used.locs.r`: implements the power analysis methods of [Street et al. 2021](https://doi.org/10.1111/2041-210X.13701)
  