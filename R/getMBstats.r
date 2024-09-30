#' Calculate statistics for the second-tier projection
#' 
#' Calculate statistics for the second-tier projection
#' 
#' @param seu Modified query Seurat object from \code{mapToMB()}.
#' @param group.by Name of one metadata column to group cells by.
#' 
#' @return A matrix of 154 rows. Each column corresponds to one group specified by group.by. 
#' Row 1-17 detail midbrain celltype proportion
#' Row 18-23 detail midbrain progenitor proportion
#' Row 24-31 detail midbrain neuron proportion
#' Row 32-131 detail ventral score
#' Row 132-137 detail the correlation between midbrain celltype proportion of current query dataset and that reported by La Manno 2016
#' Row 138-152 detail midbrain celltype identity NES
#' 
#' @author Lisheng Xu
#' 
#' @import Seurat
#' @import Matrix
#' @import dplyr
#' @import data.table
#' @import clusterProfiler
#' @export

getMBstats <- function(seu, group.by){ 
  
  seu@meta.data[["group"]] <- factor(seu@meta.data[[group.by]])
  seuMeta <- seu@meta.data
  seuMeta <- data.table(cellID = rownames(seuMeta), seuMeta)
  # Part 1: get MB celltype
  oup1 <- as.matrix(table(seuMeta$predicted.MB.celltype, seuMeta$group))
  oup1 <- t(t(oup1) / colSums(oup1))   #  Normalise to colsum=1
  # Part 2: get MB celltype NProg only
  subCT <- c("hNProg.DA","hNProg.DARN","hNProg.vGaba","hNProg.dGaba",
             "hNProg.Gaba","hNProg.Glu")
  oup2 <- seuMeta[predicted.MB.celltype %in% subCT]
  oup2$predicted.MB.celltype <- factor(oup2$predicted.MB.celltype, levels = subCT)
  oup2 <- as.matrix(table(oup2$predicted.MB.celltype, oup2$group))
  oup2 <- t(t(oup2) / colSums(oup2))   #  Normalise to colsum=1
  # Part 3: get celltype Neuron only
  subCT <- c("hDA","hDA.STN","hRN","hvGaba","hdGaba","hGlu","hOMTN","hSert")
  oup3 <- seuMeta[predicted.MB.celltype %in% subCT]
  oup3$predicted.MB.celltype <- factor(oup3$predicted.MB.celltype, levels = subCT)
  oup3 <- as.matrix(table(oup3$predicted.MB.celltype, oup3$group))
  oup3 <- t(t(oup3) / colSums(oup3))   #  Normalise to colsum=1
  # Part 4: get ventral score
  nSample = 100
  oup4 <- matrix(nrow = nSample, ncol = ncol(oup1))
  rownames(oup4) <- paste0("sample", 1:nSample)
  colnames(oup4) <- colnames(oup1)
  for(iGrp in colnames(oup4)){
    set.seed(42)
    oup4[,iGrp] <- sample(seuMeta[group == iGrp]$predicted.ventralScore, size = nSample,
                          replace = TRUE)
  }
  # Part 5: calculate corr with mannoProp
  mannoProp <- read.csv(system.file("extdata", "mannoProp.csv", package = "BrainSTEM"), row.name = 1)
  mannoProp <- as.matrix(mannoProp)
  oup5 <- cor(mannoProp, oup1) # predicted.MB.celltype needs to be factorised
  # Part 6: calculate NES with marker genes
  inpFC = data.table()
  for(iGrp in levels(seu$group)){
    seuSUB <- subset(seu, group == iGrp)
    Idents(seuSUB) <- factor(seuSUB$predicted.MB.celltype)
    for(iCT in names(which(table(Idents(seuSUB)) >= 5))){
      tmp = FoldChange(seuSUB, ident.1 = iCT)
      tmp = data.table(group = iGrp, celltype = iCT, 
                       gene = rownames(tmp), tmp)
      inpFC = rbindlist(list(inpFC, tmp))
    }
  }
  inpFC$group = factor(inpFC$group, levels = levels(seu$group))
  inpFC$celltype = factor(inpFC$celltype, levels = levels(seu$predicted.MB.celltype))
  
  markers <- fread(system.file("extdata", "celltypeMarkers.csv", package = "BrainSTEM"))
  markers$cluster <- factor(markers$cluster, levels = unique(markers$cluster))
  markers$V1 <- NULL
  oup6 <- matrix(nrow = uniqueN(markers$cluster), ncol = ncol(oup1))
  rownames(oup6) <- levels(markers$cluster)
  colnames(oup6) <- colnames(oup1)
  for(iGrp in levels(inpFC$group)){
    for(iCT in levels(inpFC$celltype)){
      inp <- inpFC[group == iGrp][celltype == iCT]
      if(nrow(inp) > 0){
        inpinp <- inp$avg_log2FC
        names(inpinp) <- inp$gene
        inpinp <- sort(inpinp, decreasing = T)
        tmpOut <- GSEA(inpinp, TERM2GENE = markers, minGSSize = 5,
                       verbose = FALSE, pvalueCutoff = 1)
        tmpOut <- data.table(data.frame(tmpOut))
        tmpOut <- tmpOut[ID == iCT]
        if(nrow(tmpOut) == 1){oup6[iCT, iGrp] = tmpOut$NES}
      }
    }
  }
  # Combine
  oupMat <- rbind(oup1, oup2, oup3, oup4, oup5, oup6)
  return(oupMat)
}