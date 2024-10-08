#' Plot statistics for the second-tier projection using \code{ComplexHeatmap}
#' 
#' Plot statistics for the second-tier projection using \code{ComplexHeatmap} (Fig 4)
#' 
#' @param inpMat A matrix from \code{getMBstats()}.
#' @param no.legend If true, a heatmap without legend is returned
#' 
#' @return A heatmap constructed using \code{ComplexHeatmap}
#' 
#' @author Lisheng Xu
#' 
#' @import Seurat
#' @import circlize
#' @import ComplexHeatmap
#' @import data.table
#' @import viridis
#' @export

plotMBstats <- function(inpMat, no.legend = FALSE){
  
  colMB1 <- c("#ea3737", "#b6523f", "#ecb1b1", "#ff8348", "#f2cece",
              "#58bd81", "#77d5bb", "#a8bd8a", "#cce6a6", "#b9e07d", 
              "#AB9EE8", "#e4cff8", "#a9b4cf", "#53a7f2", 
              "#d5f5e3", "#f8f4df", "#e8dd99")
  names(colMB1) <- c("hDA", "hDA.STN", "hNProg.DA", "hRN", "hNProg.DARN", 
                     "hvGaba", "hNProg.vGaba", "hdGaba", "hNProg.dGaba",  "hNProg.Gaba", 
                     "hGlu", "hNProg.Glu", "hOMTN", "hSert",   
                     "hOPC", "hProg", "hRgl")
  colMB2 <- colMB1[c("hNProg.DA","hNProg.DARN","hNProg.vGaba","hNProg.dGaba",
                     "hNProg.Gaba","hNProg.Glu")]
  colMB3 <- colMB1[c("hDA","hDA.STN","hRN","hvGaba","hdGaba","hGlu","hOMTN","hSert")]
  
  oupAll <- inpMat
  oupAll[is.nan(oupAll)] <- 0
  
  col_fun1 = colorRamp2(c(0, 0.5, 1), c("blue", "white", "red"))
  col_fun2 = colorRamp2(c(0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0), viridis(7))
  
  fs1 <- 14
  fs2 <- 12
  fs3 <- 8
  
  ht_list = 
    
    HeatmapAnnotation("midbrain celltype proportion" = anno_barplot(
      t(oupAll[1:17,]), height = unit(6, "cm"), bar_width = 0.8, gp = gpar(fill = colMB1, col = NA),
      axis_param = list(gp = gpar(fontsize = fs2),
                        at = c(0.2, 0.4, 0.6, 0.8, 1.0),
                        labels = c("20%","40%","60%","80%","100%"))),
      annotation_name_rot = 90, annotation_name_gp = gpar(fontsize=fs2)) %v%
    HeatmapAnnotation("midbrain progenitor\nproportion" = anno_barplot(
      t(oupAll[18:23,]), height = unit(4, "cm"), bar_width = 0.8, gp = gpar(fill = colMB2, col = NA),
      axis_param = list(gp = gpar(fontsize = fs2),
                        at = c(0.2, 0.4, 0.6, 0.8, 1.0),
                        labels = c("20%","40%","60%","80%","100%"))),
      annotation_name_rot = 90, annotation_name_gp = gpar(fontsize=fs2)) %v%
    HeatmapAnnotation("midbrain neuron\nproportion" = anno_barplot(
      t(oupAll[24:31,]), height = unit(4, "cm"), bar_width = 0.8, gp = gpar(fill = colMB3, col = NA), ylim = c(0,1),
      axis_param = list(gp = gpar(fontsize = fs2),
                        at = c(0, 0.2, 0.4, 0.6, 0.8, 1.0),
                        labels = c("0%","20%","40%","60%","80%","100%"))),
      annotation_name_rot = 90, annotation_name_gp = gpar(fontsize=fs2)) %v%
    HeatmapAnnotation("ventral\nscore" = anno_density(
      oupAll[32:131,], type = "lines", joyplot_scale = 1.2, height = unit(4, "cm"),  xlim = c(-0.5,1.5),
      axis_param = list(gp = gpar(fontsize = fs2),
                        at = c(0, 0.5, 1.0),
                        labels = c(0, 0.5, 1.0))),
      annotation_name_rot = 90, annotation_name_gp = gpar(fontsize=fs2)) %v% 
    Heatmap(oupAll[132:137,], name = "midbrain celltype\nproportion correlation", col = col_fun1, 
            
            column_title = NULL, 
            row_names_side = "left", row_names_gp = gpar(fontsize = fs2), row_title = "midbrain celltype\nproportion correlation", 
            row_title_side = "right", row_title_gp = gpar(fontsize = fs2), 
            row_title_rot = 90, show_heatmap_legend = FALSE,
            cluster_rows = FALSE, cluster_columns = FALSE, column_title_rot = 45, border = "black") %v% 
    Heatmap(oupAll[138:154,], name = "midbrain celltype\nidentity NES", col = col_fun2, na_col = "snow2",
            
            column_title = NULL, 
            row_names_side = "left", row_names_gp = gpar(fontsize = fs2), row_title = "midbrain celltype\nidentity NES", 
            row_title_side = "right", row_title_gp = gpar(fontsize = fs2), 
            row_title_rot = 90, show_heatmap_legend = FALSE,
            cluster_rows = FALSE, cluster_columns = FALSE, column_title_rot = 45, border = "black") %v% 
    columnAnnotation(text = anno_text(colnames(oupAll), rot = 45, gp = gpar(fontsize = fs2)))
  
  # legends
  lgd_list = list(
    Legend(labels = rev(names(colMB1)), legend_gp = gpar(fill = rev(colMB1)), ncol = 1, title = "midbrain celltype", labels_gp = gpar(fontsize = fs2), title_gp = gpar(fontsize = fs1, fontface = "bold")),
    Legend(col_fun = col_fun1, title = "midbrain celltype\nproportion correlation", at = c(0, 0.25, 0.5, 0.75, 1), labels_gp = gpar(fontsize = fs2), title_gp = gpar(fontsize = fs1, fontface = "bold")),
    Legend(col_fun = col_fun2, title = "midbrain celltype\nidentity NES", at = c(0, 1, 2, 3), labels_gp = gpar(fontsize = fs2), title_gp = gpar(fontsize = fs1, fontface = "bold"))
  )
  
  # Draw
  if (isTRUE(no.legend)){
    return(draw(ht_list, background = "transparent"))
  } else{
    return(draw(ht_list, heatmap_legend_list = lgd_list, heatmap_legend_side = "right", background = "transparent"))
  }
  
  
}
