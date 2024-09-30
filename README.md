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
seuQueryManno <- readRDS("manno.rds")
seuQueryManno <- mapToWB(seuQueryManno)
seuQueryManno <- mapToMB(seuQueryManno, min.nCell = 0)
seuQueryMannoWBstats <- getWBstats(seuQueryManno, group.by = "timepoint")
seuQueryMannoMBstats <- getMBstats(seuQueryManno, group.by = "timepoint")
plotWBstats(seuQueryMannoWBstats)
plotMBstats(seuQueryMannoMBstats)

# Save plots
pdf(width = 6, height = 8.2, bg = "transparent", file = "plotMannoWBstats.pdf")
plotWBstats(seuQueryMannoWBstats)
dev.off()
pdf(width = 6, height = 12, bg = "transparent", file = "plotMannoMBstats.pdf")
plotMBstats(seuQueryMannoMBstats)
dev.off()


# Run BrainSTEM on in-house midbrain differentiation dataset
seuQueryToh <- readRDS("toh.inhouse.rds")
seuQueryToh <- mapToWB(seuQueryToh)
seuQueryToh <- mapToMB(seuQueryToh)
seuQueryTohWBstats <- getWBstats(seuQueryToh, group.by = "timepoint")
seuQueryTohMBstats <- getMBstats(seuQueryToh, group.by = "timepoint")
plotWBstat(seuQueryTohWBstats)
plotMBstat(seuQueryTohMBstats)
```


