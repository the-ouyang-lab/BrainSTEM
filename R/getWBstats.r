getWBstats <- function(seuMeta, group.by){
  if (length(seuMeta$cellID != 0)){seuMeta$cellID.backup <- seuMeta$cellID;seuMeta$cellID <- NULL}
  seuMeta[["group"]] <- factor(seuMeta[[group.by]])
  seuMeta <- data.table(cellID = rownames(seuMeta), seuMeta)
  # Part 1: get unassigned assigned
  seuMeta$assign <- as.character(seuMeta$predicted.reg.celltype.major)
  seuMeta[assign != "unassigned"]$assign <- "assigned"
  oup1 <- as.matrix(table(seuMeta$assign, seuMeta$group))
  oup1 <- t(t(oup1) / colSums(oup1))   #  Normalise to colsum=1
  # Part 2: get region_MajorCelltype
  oup2 <- as.matrix(table(seuMeta$predicted.reg.celltype.major, seuMeta$group))
  oup2 <- oup2[-grep("unassigned", rownames(oup2)),]
  oup2 <- t(t(oup2) / colSums(oup2))   #  Normalise to colsum=1
  # Part 3: get celltype
  oup3 <- seuMeta[predicted.reg.celltype.major != "unassigned"]
  oup3 <- as.matrix(table(oup3$predicted.celltype, oup3$group))
  oup3 <- t(t(oup3) / colSums(oup3))   #  Normalise to colsum=1
  # Part 4: get DA region
  oup4 <- seuMeta[predicted.reg.celltype.major != "unassigned"]
  oup4 <- oup4[predicted.celltype == "DA N"]
  oup4 <- as.matrix(table(oup4$predicted.reg.celltype.major, oup4$group))
  oup4 <- oup4[c("Midbrain.Neuron","Forebrain.Neuron","Hindbrain.Neuron"), ]
  rownames(oup4) <- c("Midbrain.DA N", "Forebrain.DA N", "Hindbrain.DA N")
  oup4 <- t(t(oup4) / colSums(oup4))   #  Normalise to colsum=1
  # Combine
  oupMat <- rbind(oup1, oup2, oup3, oup4)
  return(oupMat)
}
