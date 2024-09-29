#' Map query dataset to the fetal midbrain subatlas
#'
#' Wrapper function to map a query dataset onto the fetal midbrain subatlas.
#' This function uses the Seurat algorithm to project the query dataset onto
#' the reference. Upon projection, the midbrain detailed celltype annotation is
#' predicted. Furthermore, a ventral score is predicted for each single cell.
#'
#' @param seuMBRef The Seurat object of the fetal midbrain subatlas, usually
#'   downloaded by function \code{downloadMBref()}
#' @param seuQuery The Seurat object of the query dataset, with log-normalised
#'   data and PC calculated. This Seurat object should already have been mapped
#'   to the fetal whole brain atlas using \code{mapToWB()} and contains
#'   predicted.reg.celltype.major in the metadata
#' @param min.nCell Minimal count of midbrain cells according to the first-tier
#'   mapping in seuQuery to allow for second tier mapping by \code{mapToMB()}
#' @param returnAnchor Logical. Setting it to TRUE allows the transfer anchors
#'   to be returned
#'
#' @return A modified query Seurat object containing (i) an projected PCA
#'   reduction, (ii) projected UMAP coordinates, (iii) new assay corresponding
#'   to predicted.MB.celltype probabilities and (iv) new columns in the
#'   metadata: predicted.MB.celltype, predicted.MB.celltype.score and
#'   predicted.ventralScore.
#'   If returnAnchor is set to TRUE, a list containing such Seurat object and
#'   the transfer anchors is returned.
#'
#' @author Lisheng Xu
#'
#' @import Seurat
#' @import dplyr
#' @import ranger
#' @export



mapToMB <- function(seuMBRef, seuQuery, min.nCell = 1000, returnAnchor = FALSE){
  model <- seuMBRef@misc$ventralScoreModel

  if(nrow(seuQuery@meta.data[seuQuery@meta.data[["predicted.reg.celltype.major"]] %in% c(
    "Midbrain.Neuron", "Midbrain.Prog"),]) < min.nCell){
    stop(paste0("Aborted: The query must contain not fewer than ", min.nCell, " midbrain cells."))
  }
  # - Step1: Filter out midbrain neural cells
  # --- Uses predicted.reg.celltype.major
  seuQuery <- subset(seuQuery, subset =
                       predicted.reg.celltype.major %in% c("Midbrain.Neuron", "Midbrain.Prog"))
  seuQuery <- seuQuery %>% ScaleData() %>%
    RunPCA(npcs = 50) %>% RunUMAP(reduction = "pca", dims = 1:50)

  # - Step2: Run Seurat FindTransferAnchors & MapQuery
  anchors <- FindTransferAnchors(reference = seuMBRef, query = seuQuery,
                                 reference.assay = "integrated", reference.reduction = "pca",
                                 dims =  1:50, reduction = "rpca",
                                 k.filter = NA)
  seuQuery <- MapQuery(anchorset = anchors, reference = seuMBRef, query = seuQuery,
                       refdata = list(MB.celltype = "celltype"),
                       reference.reduction = "pca", reference.dims = 1:50,
                       reduction.model = "umap")

  # - Step3: Calculate ventral score
  # --- Loads ranger model from seuMBRef@misc
  inp <- Embeddings(seuQuery, reduction = "ref.pca"); colnames(inp) <- paste0("PC_", 1:50)
  tmp <- predict(model, data = inp)
  seuQuery@meta.data[["predicted.ventralScore"]] <- tmp[["predictions"]][, 2]


  # return value
  seuQuery@meta.data[["predicted.MB.celltype"]] <-
    factor(seuQuery@meta.data[["predicted.MB.celltype"]], levels = c(
      "hDA", "hDA.STN", "hNProg.DA", "hRN", "hNProg.DARN",
      "hvGaba", "hNProg.vGaba", "hdGaba", "hNProg.dGaba",  "hNProg.Gaba",
      "hGlu", "hNProg.Glu", "hOMTN", "hSert",
      "hOPC", "hProg", "hRgl"))
  if (returnAnchor){
    return(list(seuQuery, anchors))
  } else{
    return(seuQuery)
  }
}
