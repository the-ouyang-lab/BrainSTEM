#' Calculate statistics for the first-tier mapping
#' 
#' Calculate statistics for the first-tier mapping
#' 
#' @param seu Modified query Seurat object from \code{mapToMB()}.
#' @param group.by Name of one metadata column to group cells by.
#' 
#' @return A matrix of 37 rows. Each column corresponds to one group specified by group.by. 
#' Row 1-2 detail the proportion of confidently assigned and unassigned cells from the first-tier mapping
#' Row 3-11 detail the region.lineage proportions for confidently assigned cells
#' Row 12-34 detail the cell type proportions for confidently assigned cells
#' Row 35-37 detail the brain region proportions for confidently assigned DA neurons
#' 
#' @author Lisheng Xu
#' 
#' @import Seurat
#' @import Matrix
#' @import dplyr
#' @import data.table
#' @export

getWBstats <- function(seu, group.by){
  seuMeta <- seu@meta.data
  seuMeta[["group"]] <- factor(seuMeta[[group.by]])
  seuMeta <- data.table(cellID = rownames(seuMeta), seuMeta)
  # Part 1: get unassigned assigned
  oup1 <- as.matrix(table(seuMeta$is.assigned, seuMeta$group))
  oup1 <- t(t(oup1) / colSums(oup1))   #  Normalise to colsum=1
  # Part 2: get region_MajorCelltype
  oup2 <- as.matrix(table(seuMeta$predicted.reg.celltype.major, seuMeta$group))
  if(uniqueN(seuMeta$group) == 1){
    oup2 <- oup2[-grep("unassigned", rownames(oup2)),] %>% as.data.frame()
    colnames(oup2) <- unique(seuMeta$group)
  } else {
    oup2 <- oup2[-grep("unassigned", rownames(oup2)),]
  }
  oup2 <- t(t(oup2) / colSums(oup2))   #  Normalise to colsum=1
  # Part 3: get celltype
  oup3 <- seuMeta[is.assigned == "assigned"]
  oup3 <- as.matrix(table(oup3$predicted.celltype, oup3$group))
  oup3 <- t(t(oup3) / colSums(oup3))   #  Normalise to colsum=1
  # Part 4: get DA region
  oup4 <- seuMeta[is.assigned == "assigned"]
  oup4 <- oup4[predicted.celltype == "DA N"]
  oup4[["predicted.reg.celltype.major"]] <- factor(oup4[["predicted.reg.celltype.major"]], levels = c(
    "Midbrain.Neuron", "Forebrain.Neuron", "Hindbrain.Neuron"))
  oup4 <- as.matrix(table(oup4$predicted.reg.celltype.major, oup4$group))
  rownames(oup4) <- c("Midbrain.DA N", "Forebrain.DA N", "Hindbrain.DA N")
  oup4 <- t(t(oup4) / colSums(oup4))   #  Normalise to colsum=1
  # Combine
  oupMat <- rbind(oup1, oup2, oup3, oup4)
  return(oupMat)
}
