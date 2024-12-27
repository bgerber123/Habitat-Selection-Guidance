# Habitat selection modeling guidance (in practice)

## A repository for:

Gerber, BD, Setash, C, Ivan, J and Northrup, J. A plain language review and guidance for modeling animal habitat-selection. 

---

## What's in this repository?

This repository comprises two R markdown files and associated files.

---

## The working directory

### Main files

The two fundamental scripts are `TraditionalHSF.Rmd` and `MovementHSF.Rmd`. Both files are organized 
similarly, setting up the environment and data, and then follow the subsections of the manuscript. 
The first file focuses on addressing the manuscripts points regarding traditional broad-scale
habitat selection, for example, selection within an individual's home range. The second file
focues on addressing the manuscripts points regarding movement-based habitat selection or
step-selection functions.

<a href="https://bgerber123.github.io/Habitat-Selection-Guidance
/TraditionalHSF.html" title="">Traditional HSF (Link)</a>

<a href="MovementHSF.html" title="">Traditional HSF (Link)</a>

### Additional Files

- `sim.ind.movement.hsf.r`: simulation of individual movement based habibtat selection
- `mean_ci_boot.r`: bootsrapping function to get population (across individual) infernece
- `Covs`: rasters of spatial covariates
  