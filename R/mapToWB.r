#' Map query dataset to the fetal whole brain atlas
#' 
#' Map query dataset to the fetal whole brain atlas
#' 
#' @param seuWBRef The Seurat object of the fetal whole brain atlas, usually downloaded by function \code{downloadWBref()}
#' @param seuQuery The Seurat object of the query dataset, with log-normalised data and PC calculated
#' @param returnAnchor Logical. Setting it to TRUE allows the transfer anchors to be returned
#' 
#' @return A modified query Seurat object containing:
#' new assay corresponding to predicted.celltype
#' an integrated reduction
#' a projected UMAP reduction 
#' new columns in the metadata: is.assigned, predicted.celltype, predicted.reg.celltype, predicted.reg.celltype.major and multiple scores. 
#' If returnAnchor is set to TRUE, a list containing such Seurat object and the transfer anchors is returned.
#' 
#' @author Lisheng Xu
#' 
#' @import Seurat
#' @import dplyr
#' @export



mapToWB <- function(seuWBRef, seuQuery, returnAnchor = FALSE){
  markerList <- seuWBRef@misc$markerList # celltype specific markers
  markerList.2 <- seuWBRef@misc$markerList.2 # region specific markers for each celltype
  regspe <- seuWBRef@misc$regionSpecificityTb #conca terms to predicted terms mapping table
  
  # - Step1: Calculates module score for celltype and region (for all cells)
  # --- Loads marker genes from seuWBRef@misc
  seuQuery <- AddModuleScore(seuQuery, assay = "RNA", features = markerList, name = "celltypeScore")
  colnames(seuQuery@meta.data)[grep(pattern = "^celltypeScore", x = colnames(seuQuery@meta.data))] <- paste0("score_", names(markerList)) # need to deal with the case where none of the markers for the specific celltype exists in seuQuery
  seuQuery <- AddModuleScore(seuQuery, assay = "RNA", features = markerList.2, name = "regScore")
  colnames(seuQuery@meta.data)[grep(pattern = "^regScore", x = colnames(seuQuery@meta.data))] <- paste0("regScore_", names(markerList.2))
  
  # - Step2: Run Seurat FindTransferAnchors & MapQuery
  anchors <- FindTransferAnchors(reference = seuWBRef, query = seuQuery,
                                 reference.assay = "integrated", reference.reduction = "pca", dims =  1:50,
                                 reduction = "rpca",
                                 k.filter = NA)
  seuQuery <- MapQuery(anchorset = anchors, reference = seuWBRef, query = seuQuery, 
                       refdata = list(celltype = "celltype"), 
                       reference.reduction = "pca", reference.dims = 1:50, reduction.model = "umap")
  
  # - Step3: Set assign/unassign
  tmpMeta <- seuQuery@meta.data
  tmpMeta[["cellID"]] <- rownames(tmpMeta)
  # how to deal with existing cellID
  tmpMeta[["is.assigned"]] <- "unassigned"
  for(i in 1:nrow(tmpMeta)){
    predicted.celltype <- tmpMeta[["predicted.celltype"]][i]
    if(tmpMeta[[paste0("score_", predicted.celltype)]][i] > 0){
      tmpMeta[["is.assigned"]][i] <- "assigned"
    }
  }
  
  # - Step4: For each RS-celltype (among assigned cells), assign region
  # - Step5: label clean up e.g. NRS
  ctyNonNeural <- c("Epithelium", "Mesoderm", "LPM", "Endothelium", "Blood", "Immune", "Notochord", "Endoderm")
  tmpMeta.n <- tmpMeta[!tmpMeta[["predicted.celltype"]] %in% ctyNonNeural,]
  tmpMeta.nn <- tmpMeta[tmpMeta[["predicted.celltype"]] %in% ctyNonNeural,]
  
  tmpMeta.n[["conca"]] <- 0
  for(i in 1:nrow(tmpMeta.n)){
    predicted.celltype <- tmpMeta.n[["predicted.celltype"]][i]
    scoreName <- colnames(tmpMeta.n)[grepl("^regScore_",colnames(tmpMeta.n)) & grepl(predicted.celltype,colnames(tmpMeta.n), fixed = TRUE)]
    tmpMeta.n[["conca"]][i] <- apply(tmpMeta.n[i, scoreName], MARGIN = 1,
                                     FUN = function(x) {
                                       conca <- scoreName[which.max(x)]
                                       conca <- gsub("regScore_", "", conca)
                                       return(conca)
                                     })
  }
  
  tmpMeta.n <- merge(tmpMeta.n, regspe, by = "conca")
  tmpMeta.n[["conca"]] <- NULL
  colnames(tmpMeta.n)[grepl("^reg.anno$",colnames(tmpMeta.n))] <- "predicted.reg.celltype"
  colnames(tmpMeta.n)[grepl("^reg.anno.major$",colnames(tmpMeta.n))] <- "predicted.reg.celltype.major"
  rownames(tmpMeta.n) <- tmpMeta.n[["cellID"]]
  
  
  
  tmpMeta.nn[["predicted.reg.celltype"]] <- tmpMeta.nn[["predicted.celltype"]]
  tmpMeta.nn[["predicted.reg.celltype.major"]] <- "nonNeural"
  
  tmpMeta <- rbind(tmpMeta.n,tmpMeta.nn)
  tmpMeta[tmpMeta[["is.assigned"]] == "unassigned", "predicted.reg.celltype"] <- "unassigned"
  tmpMeta[tmpMeta[["is.assigned"]] == "unassigned", "predicted.reg.celltype.major"] <- "unassigned"
  
  
  
  # return value
  seuQuery <- AddMetaData(seuQuery, metadata = tmpMeta[,c("is.assigned", "predicted.reg.celltype", "predicted.reg.celltype.major")])
  seuQuery@meta.data[["is.assigned"]] <- factor(seuQuery@meta.data[["is.assigned"]], levels = c("assigned", "unassigned"))
  seuQuery@meta.data[["predicted.celltype"]] <- factor(seuQuery@meta.data[["predicted.celltype"]], levels = c("DA N", "Cholin N", "GABA N", "Glu N", "Ser N", 
                                                                                                              "CA10+TFAP2B+ N", "FSTL4+RELN+ N",  "LHX2+LHX9+ N", "LHX5+OTP+ N", "SIM1+PITX2+ N", 
                                                                                                              "Radial Glia", "Neuroblast", "Oligo", "Neural Crest", "Glioblast", 
                                                                                                              "Epithelium", "Mesoderm", "LPM", "Endothelium", "Blood", "Immune", "Notochord", "Endoderm"))
  seuQuery@meta.data[["predicted.reg.celltype.major"]] <- factor(seuQuery@meta.data[["predicted.reg.celltype.major"]], levels = c("Midbrain.Neuron", "Midbrain.Prog", 
                                                                                                                                  "Forebrain.Neuron", "Forebrain.Prog", 
                                                                                                                                  "Hindbrain.Neuron", "Hindbrain.Prog", 
                                                                                                                                  "NRS.Neuron", "NRS.Prog", "nonNeural", "unassigned"))
  if (isTRUE(returnAnchor)){
    return(list(seuQuery, anchors))
  } else{
    return(seuQuery)
  }
}