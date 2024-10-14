#' Plot statistics for the first-tier mapping using \code{ComplexHeatmap}
#' 
#' Plot statistics for the first-tier mapping using \code{ComplexHeatmap} (Fig 3)
#' 
#' @param inpMat A matrix from \code{getWBstats()}
#' @param plot blocks of the plots to be printed, i.e. c(1,2,3,4)
#' 1: the proportion of confidently assigned and unassigned cells from the first-tier mapping
#' 2: the region.lineage proportions for confidently assigned cells
#' 3: the cell type proportions for confidently assigned cells
#' 4: detail the brain region proportions for confidently assigned DA neurons
#' @param no.legend If true, a heatmap without legend is returned
#' 
#' @return A heatmap constructed using \code{ComplexHeatmap}
#' 
#' @author Lisheng Xu
#' 
#' @import Seurat
#' @import ComplexHeatmap
#' @import data.table
#' @export

plotWBstats <- function(inpMat, plot = c(1,2,3,4), no.legend = FALSE){
  colWB1 <- c("#043259", "#BBBBBB")
  names(colWB1) <- c("assigned", "unassigned")
  colWB2 <- c("#d42f2f", "#d48c8c", "#16b87d","#b5e8d5", "#0b5394", "#6fa8dc", 
              "#f8d6b3", "#EDDFD6", "#dedede")
  names(colWB2) <- c("Midbrain.Neuron", "Midbrain.Prog", 
                     "Forebrain.Neuron", "Forebrain.Prog", 
                     "Hindbrain.Neuron", "Hindbrain.Prog", 
                     "NRS.Neuron", "NRS.Prog", "nonNeural")
  colWB3 <- c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00", "#FFFF33", "#A65628", "#F781BF", "#1B9E77", "#D95F02", "#7570B3", "#E7298A", "#66A61E", "#E6AB02", "#A6761D",
              "#FBB4AE", "#B3CDE3", "#CCEBC5", "#DECBE4", "#FED9A6", "#FFADAC", "#E5D8BD", "#FBD7EC")
  names(colWB3) <- c("DA N", "Cholin N", "GABA N", "Glu N", "Ser N", "CA10+TFAP2B+ N", "FSTL4+RELN+ N", "LHX2+LHX9+ N",  
                     "LHX5+OTP+ N", "SIM1+PITX2+ N", "Radial Glia", "Neuroblast", "Oligo", "Neural Crest", "Glioblast", "Epithelium",    
                     "Mesoderm", "LPM", "Endothelium", "Blood", "Immune", "Notochord", "Endoderm"  )
  colWB4 <- c("#d42f2f", "#16b87d", "#0b5394")
  names(colWB4) <- c("Midbrain.DA N", "Forebrain.DA N", "Hindbrain.DA N")
  
  if (nrow(inpMat) != 37){
    stop("The input matrix should contain 37 rows.")
  }
  if (!all(plot %in% c(1, 2, 3, 4))) {
    stop("\"plot\" must contain only values 1, 2, 3, and/or 4.")
  }
  oupAll <- inpMat
  oupAll[is.nan(oupAll)] <- 0
  
  # fontsize
  fs1 <- 14
  fs2 <- 12
  fs3 <- 8
  # Main plot
  ht_list = NULL
  if (1 %in% plot){
    ht_list = ht_list %v%
      HeatmapAnnotation("proportion of cells\nassigned" = anno_barplot(
        t(oupAll[1:2,]), height = unit(3, "cm"), bar_width = 0.8, gp = gpar(fill = colWB1, col = NA), 
        axis_param = list(gp = gpar(fontsize = fs2),
                          at = c(0.2, 0.4, 0.6, 0.8, 1.0), 
                          labels = c("20%","40%","60%","80%","100%"))),
        annotation_name_rot = 90, annotation_name_gp = gpar(fontsize=fs2)) 
  }
  if (2 %in% plot){
    ht_list = ht_list %v%
      HeatmapAnnotation("region.lineage\nproportion" = anno_barplot(
        t(oupAll[3:11,]), height = unit(6, "cm"), bar_width = 0.8, gp = gpar(fill = colWB2, col = NA), 
        axis_param = list(gp = gpar(fontsize = fs2),
                          at = c(0.2, 0.4, 0.6, 0.8, 1.0), 
                          labels = c("20%","40%","60%","80%","100%"))),
        annotation_name_rot = 90, annotation_name_gp = gpar(fontsize=fs2)) 
  }
  if (3 %in% plot){
    ht_list = ht_list %v%
      HeatmapAnnotation("celltype\nproportion" = anno_barplot(
        t(oupAll[12:34,]), height = unit(6, "cm"), bar_width = 0.8, gp = gpar(fill = colWB3, col = NA), 
        axis_param = list(gp = gpar(fontsize = fs2),
                          at = c(0.2, 0.4, 0.6, 0.8, 1.0), 
                          labels = c("20%","40%","60%","80%","100%"))),
        annotation_name_rot = 90, annotation_name_gp = gpar(fontsize=fs2))
  }
  if (4 %in% plot){
    ht_list = ht_list %v%
      HeatmapAnnotation("DA N region\nproportion" = anno_barplot(
        t(oupAll[35:37,]), height = unit(3, "cm"), bar_width = 0.8, gp = gpar(fill = colWB4, col = NA), ylim = c(0,1),
        axis_param = list(gp = gpar(fontsize = fs2),
                          at = c(0, 0.2, 0.4, 0.6, 0.8, 1.0),
                          labels = c("0%","20%","40%","60%","80%","100%"))),
        annotation_name_rot = 90, annotation_name_gp = gpar(fontsize=fs2)) 
  } 
  ht_list = ht_list %v%
    columnAnnotation(text = anno_text(colnames(oupAll), rot = 45, gp = gpar(fontsize = fs2))) %v% 
    Heatmap(matrix(ncol = ncol(oupAll), nrow = 1), name = "total cells", 
            height = unit(0, "npc"),
            column_title = NULL, show_row_names = FALSE, show_heatmap_legend = FALSE,
            cluster_rows = FALSE, cluster_columns = FALSE)
  # legends
  lgd_list = list()
  if (1 %in% plot){
    lgd_list = c(lgd_list, list(
      Legend(labels = rev(names(colWB1)), legend_gp = gpar(fill = rev(colWB1)), title = "assignment", labels_gp = gpar(fontsize = fs2), title_gp = gpar(fontsize = fs1, fontface = "bold"))
    ))
  }
  if (2 %in% plot){
    lgd_list = c(lgd_list, list(
      Legend(labels = rev(names(colWB2)), legend_gp = gpar(fill = rev(colWB2)), ncol = 2, title = "region.lineage", labels_gp = gpar(fontsize = fs2), title_gp = gpar(fontsize = fs1, fontface = "bold"))
    ))
  }
  if (3 %in% plot){
    lgd_list = c(lgd_list, list(
      Legend(labels = rev(names(colWB3)), legend_gp = gpar(fill = rev(colWB3)), ncol = 2, title = "celltype", labels_gp = gpar(fontsize = fs2), title_gp = gpar(fontsize = fs1, fontface = "bold"))
    ))
  }
  if (4 %in% plot){
    lgd_list = c(lgd_list, list(
      Legend(labels = rev(names(colWB4)), legend_gp = gpar(fill = rev(colWB4)), title = "DA N region", labels_gp = gpar(fontsize = fs2), title_gp = gpar(fontsize = fs1, fontface = "bold"))
    ))
  }
  # Draw
  if (isTRUE(no.legend)){
    return(draw(ht_list, background = "transparent"))
  } else{
    return(draw(ht_list, heatmap_legend_list = lgd_list, heatmap_legend_side = "right", background = "transparent"))
  }
}
