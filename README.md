# BrainSTEM

BrainSTEM (brain Single-cell Two tiEr Mapping) is a biologically-driven 
multi-resolution mapping strategy to assess the transcriptional fidelity of 
midbrain cultures via scRNA-seq. In the first tier, query datasets are projected 
onto the fetal whole brain atlas to assess the brain region specificity of the 
various neural populations. Following this, midbrain-specific cells are further 
projected onto a more intricately annotated midbrain subatlas. 

Here is some example code to run BrainSTEM on the La Manno 2016 in vivo midbrain 
and our in-house midbrain differentiation dataset: 

``` r
# Read in ref and query seurat objects
seuRefWB <- readRDS("path-to-whole-brain-atlas.rds")
seuRefMB <- readRDS("path-to-midbrain-subatlas.rds")


# Run BrainSTEM on La Manno 2016
seuQueryManno <- readRDS("path-to-lamanno2016.rds")
seuQueryManno <- mapToWB(seuQueryManno)
seuQueryManno <- mapToMB(seuQueryManno, min.nCell = 100)
seuQueryManno <- getWBstat(seuQueryManno, group.by = "Timepoint")
seuQueryManno <- getMBstat(seuQueryManno, group.by = "Timepoint")
plotWBstat(seuQueryManno)
plotMBstat(seuQueryManno)


# Run BrainSTEM on in-house midbrain differentiation dataset
seuQueryToh <- readRDS("path-to-inhousedata.rds")
seuQueryToh <- mapToWB(seuQueryToh)
seuQueryToh <- mapToMB(seuQueryToh)
seuQueryToh <- getWBstat(seuQueryToh, group.by = "group")
seuQueryToh <- getMBstat(seuQueryToh, group.by = "group")
plotWBstat(seuQueryToh)
plotMBstat(seuQueryToh)
```


