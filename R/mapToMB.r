#' Map query dataset to the fetal midbrain subatlas
#' 
#' Map query dataset to the fetal midbrain subatlas
#' 
#' @param seuMBRef The Seurat object of the fetal midbrain subatlas, usually downloaded by function downloadMBref()
#' @param seuQuery The Seurat object of the query dataset, with log-normalised data and PC calculated, 
#' which has been mapped to the fetal whole brain atlas and contains predicted.reg.celltype.major in the metadata
#' @param min.nCell Minimal count of midbrain cells according to the first-tier mapping in seuQuery to allow it 
#' to be mapped by mapToMB()
#' @param returnAnchor Logical. Setting it to TRUE allows the transfer anchors to be returned
#' 
#' @return A modified query Seurat object containing:
#' new assay corresponding to predicted.MB.celltype
#' an integrated reduction
#' a projected UMAP reduction 
#' new columns in the metadata: predicted.MB.celltype, predicted.MB.celltype.score and predicted.ventralScore. 
#' If returnAnchor is set to TRUE, a list containing such Seurat object and the transfer anchors is returned.
#' 
#' @author Lisheng Xu
#' 
#' @import Seurat
#' @import dplyr
#' @import ranger
#' @export



mapToMB <- function(seuMBRef, seuQuery, min.nCell = 1000,returnAnchor = FALSE){
  model <- seuMBRef@misc$ventralScoreModel
  
  if(nrow(seuQuery@meta.data[seuQuery@meta.data[["predicted.reg.celltype.major"]] %in% c("Midbrain.Neuron", "Midbrain.Prog"),]) <min.nCell){
    stop(paste0("Aborted: The query must contain not fewer than ", min.nCell, " midbrain cells."))
  }
  # - Step1: Filter out midbrain neural cells
  # --- Uses predicted.reg.celltype.major
  seuQuery <- subset(seuQuery, subset = predicted.reg.celltype.major %in% c("Midbrain.Neuron", "Midbrain.Prog"))
  seuQuery <- seuQuery %>% ScaleData() %>%
    RunPCA(npcs = 50) %>% RunUMAP(reduction = "pca", dims = 1:50)
  
  # - Step2: Run Seurat FindTransferAnchors & MapQuery
  anchors <- FindTransferAnchors(reference = seuMBRef, query = seuQuery,
                                 reference.assay = "integrated", reference.reduction = "pca", dims =  1:50,
                                 reduction = "rpca",
                                 k.filter = NA)
  seuQuery <- MapQuery(anchorset = anchors, reference = seuMBRef, query = seuQuery, 
                       refdata = list(MB.celltype = "celltype"), 
                       reference.reduction = "pca", reference.dims = 1:50, reduction.model = "umap")
  
  # - Step3: Calculate ventral score
  # --- Loads ranger model from seuMBRef@misc
  inp <- Embeddings(seuQuery, reduction = "ref.pca"); colnames(inp) <- paste0("PC_", 1:50)
  tmp <- predict(model, data = inp)
  seuQuery@meta.data[["predicted.ventralScore"]] <- tmp[["predictions"]][, 2]
  
  
  # return value
  if (returnAnchor){
    return(list(seuQuery, anchors))
  } else{
    return(seuQuery)
  }
}
