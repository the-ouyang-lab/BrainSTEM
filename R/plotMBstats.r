#' Plot statistics for the second-tier projection using \code{ComplexHeatmap}
#' 
#' Plot statistics for the second-tier projection using \code{ComplexHeatmap} (Fig 4)
#' 
#' @param inpMat A matrix from \code{getMBstats()}, or a list of multiple matrices. If names detected in the list, they are shown in the legend as study names.
#' @param plot blocks of the plots to be printed, i.e. c(1,2,3,4,5,6)
#' 1: midbrain celltype proportion
#' 2: midbrain progenitor proportion
#' 3: midbrain neuron proportion
#' 4: ventral score
#' 5: the correlation between midbrain celltype proportion of current query dataset and that reported by La Manno 2016
#' 6: midbrain celltype identity NES
#' @param show.manno Display an in vivo fetal midbrain dataset published by La Manno 2016 for comparison.
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

plotMBstats <- function(inpMat, plot = c(1,2,3,4,5,6), show.manno = FALSE, no.legend = FALSE){
  colDat <- c("#e8b571", "#cc1e32", "#f2e291", "#71d4e8", "#e88fd7", 
              "#db53c2", "#ca75f4", "#047a60", "#cc262c", "#3c67ad",
              "#05af2a", "#66e8e3", "#d4ed84", "#898dff", "#ace03c",
              "#1a737f", "#dd4042", "#36c497", "#dd6aaf", "#7f71dd")
  
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
  
  if (is.list(inpMat)){
    
    if (unique(unlist(lapply(inpMat, nrow))) != 154){
      stop("The input matrix should contain 154 rows.")
    }
    
    oupAll <- do.call(cbind, inpMat)
    oupAll[is.nan(oupAll)] <- 0
    show.study <- TRUE
    
  } else {
    
    if (nrow(inpMat) != 154){
      stop("The input matrix should contain 154 rows.")
    }
    
    oupAll <- inpMat
    oupAll[is.nan(oupAll)] <- 0
    show.study <- FALSE
  }
  
  if (isTRUE(show.manno)){
    mannoMBstats <- read.csv(system.file("extdata", "mannoMBstats.csv", package = "BrainSTEM")) 
    mannoMBstats.rowname <- mannoMBstats[,1]
    mannoMBstats <- as.matrix(mannoMBstats[,-1])
    rownames(mannoMBstats) <- mannoMBstats.rowname
    show.study <- TRUE
    if (is.list(inpMat)){
      inpMat <- c(list("La Manno 2016" = mannoMBstats),
                  inpMat)
      
    } else {
      inpMat <- list("La Manno 2016" = mannoMBstats,
                     inpMat)
    }
    oupAll <- cbind(mannoMBstats, oupAll)
  }
  
  if (isTRUE(show.study)){
    
    if (!is.null(names(inpMat))){
      studyname <- names(inpMat)
      for (i in seq_along(studyname)){
        if (studyname[i] == ""){
          studyname[i] <- paste0("Dataset ", i)
        }
      }
    } else {
      studyname <- paste0("Dataset ", seq_len(length(inpMat)))
    }
    colDat <- colDat[seq_along(studyname)]
    names(colDat) <- studyname
    
    grpByDataset <- sapply(seq_along(studyname), function(x) return(rep(studyname[x], ncol(inpMat[[x]]))), USE.NAMES = FALSE) %>% do.call(c,.)
    grpByDataset <- factor(grpByDataset, levels = unique(grpByDataset))
    listDataset  = list(); for(iD in unique(grpByDataset)){listDataset[[iD]] = which(grpByDataset == iD)}
    
  }
  
  
  if (!all(plot %in% c(1, 2, 3, 4, 5, 6))) {
    stop("\"plot\" must contain only values 1, 2, 3, 4, 5 and/or 6.")
  }
  
  
  # col and fontsize
  col_fun1 = colorRamp2(c(0, 0.5, 1), c("blue", "white", "red"))
  col_fun2 = colorRamp2(c(0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0), viridis(7))
  
  fs1 <- 14
  fs2 <- 12
  
  
  # Main plot
  ht_list = NULL
  if (isTRUE(show.study)){
    ht_list = ht_list %v% HeatmapAnnotation(
      empty = anno_empty(border = FALSE),
      dataset = anno_block(align_to = listDataset, panel_fun = function(index, levels) {
        grid.rect(gp = gpar(fill = colDat[levels], col = "black"))
      }), show_legend = TRUE)
  }
  if (1 %in% plot){
    ht_list = ht_list %v%
      HeatmapAnnotation("midbrain celltype proportion" = anno_barplot(
        t(oupAll[1:17,]), height = unit(6, "cm"), bar_width = 0.8, gp = gpar(fill = colMB1, col = NA),
        axis_param = list(gp = gpar(fontsize = fs2),
                          at = c(0.2, 0.4, 0.6, 0.8, 1.0),
                          labels = c("20%","40%","60%","80%","100%"))),
        annotation_name_rot = 90, annotation_name_gp = gpar(fontsize=fs2))
  }
  if (2 %in% plot){
    ht_list = ht_list %v%
      HeatmapAnnotation("midbrain progenitor\nproportion" = anno_barplot(
        t(oupAll[18:23,]), height = unit(4, "cm"), bar_width = 0.8, gp = gpar(fill = colMB2, col = NA),
        axis_param = list(gp = gpar(fontsize = fs2),
                          at = c(0.2, 0.4, 0.6, 0.8, 1.0),
                          labels = c("20%","40%","60%","80%","100%"))),
        annotation_name_rot = 90, annotation_name_gp = gpar(fontsize=fs2)) 
  }
  if (3 %in% plot){
    ht_list = ht_list %v%
      HeatmapAnnotation("midbrain neuron\nproportion" = anno_barplot(
        t(oupAll[24:31,]), height = unit(4, "cm"), bar_width = 0.8, gp = gpar(fill = colMB3, col = NA), ylim = c(0,1),
        axis_param = list(gp = gpar(fontsize = fs2),
                          at = c(0, 0.2, 0.4, 0.6, 0.8, 1.0),
                          labels = c("0%","20%","40%","60%","80%","100%"))),
        annotation_name_rot = 90, annotation_name_gp = gpar(fontsize=fs2))
  }
  if (4 %in% plot){
    ht_list = ht_list %v%
      HeatmapAnnotation("ventral\nscore" = anno_density(
        oupAll[32:131,], type = "lines", joyplot_scale = 1.2, height = unit(4, "cm"),  xlim = c(-0.25,1.25),
        axis_param = list(gp = gpar(fontsize = fs2),
                          at = c(0, 0.5, 1.0),
                          labels = c(0, 0.5, 1.0))),
        annotation_name_rot = 90, annotation_name_gp = gpar(fontsize=fs2)) 
  }
  if (5 %in% plot){
    ht_list = ht_list %v%
      Heatmap(oupAll[132:137,], name = "midbrain celltype\nproportion correlation", col = col_fun1, 
              column_split = if(isTRUE(show.study)) grpByDataset else NULL, 
              column_gap = if(isTRUE(show.study)) unit(2, "mm") else NULL, 
              column_title = NULL, 
              row_names_side = "left", row_names_gp = gpar(fontsize = fs2), row_title = "midbrain celltype\nproportion correlation", 
              row_title_side = "right", row_title_gp = gpar(fontsize = fs2), 
              row_title_rot = 90, show_heatmap_legend = FALSE,
              cluster_rows = FALSE, cluster_columns = FALSE, column_title_rot = 45, border = "black") 
  }
  if (6 %in% plot){
    ht_list = ht_list %v%
      Heatmap(oupAll[138:154,], name = "midbrain celltype\nidentity NES", col = col_fun2, na_col = "snow2",
              column_split = if(isTRUE(show.study)) grpByDataset else NULL, 
              column_gap = if(isTRUE(show.study)) unit(2, "mm") else NULL, 
              column_title = NULL, 
              row_names_side = "left", row_names_gp = gpar(fontsize = fs2), row_title = "midbrain celltype\nidentity NES", 
              row_title_side = "right", row_title_gp = gpar(fontsize = fs2), 
              row_title_rot = 90, show_heatmap_legend = FALSE,
              cluster_rows = FALSE, cluster_columns = FALSE, column_title_rot = 45, border = "black") 
  }
  ht_list = ht_list %v%
    columnAnnotation(text = anno_text(colnames(oupAll), rot = 45, gp = gpar(fontsize = fs2))) %v% 
    Heatmap(matrix(ncol = ncol(oupAll), nrow = 1), name = "total cells", 
            height = unit(0, "npc"),
            column_title = NULL, show_row_names = FALSE, show_heatmap_legend = FALSE,
            cluster_rows = FALSE, cluster_columns = FALSE)
  
  # legends
  lgd_list = list()
  if (isTRUE(show.study)){
    lgd_list = c(lgd_list, list(
      Legend(labels = studyname, legend_gp = gpar(fill = colDat), title = "study", labels_gp = gpar(fontsize = fs2), ncol = ifelse(length(studyname) <= 6, 1, 2), title_gp = gpar(fontsize = fs1, fontface = "bold"))
    ))
  }
  if (any(c(1,2,3) %in% plot)){
    lgd_list = c(lgd_list, list(
      Legend(labels = rev(names(colMB1)), legend_gp = gpar(fill = rev(colMB1)), ncol = 1, title = "midbrain celltype", labels_gp = gpar(fontsize = fs2), title_gp = gpar(fontsize = fs1, fontface = "bold"))
    ))
  }
  if (5 %in% plot){
    lgd_list = c(lgd_list, list(
      Legend(col_fun = col_fun1, title = "midbrain celltype\nproportion correlation", at = c(0, 0.25, 0.5, 0.75, 1), labels_gp = gpar(fontsize = fs2), title_gp = gpar(fontsize = fs1, fontface = "bold"))
    ))
  }  
  if (6 %in% plot){
    lgd_list = c(lgd_list, list(
      Legend(col_fun = col_fun2, title = "midbrain celltype\nidentity NES", at = c(0, 1, 2, 3), labels_gp = gpar(fontsize = fs2), title_gp = gpar(fontsize = fs1, fontface = "bold"))
    ))
  }
  # Draw
  if (isTRUE(no.legend)){
    return(draw(ht_list, background = "transparent"))
  } else{
    return(draw(ht_list, heatmap_legend_list = lgd_list, heatmap_legend_side = "right", background = "transparent"))
  }
  
  
}
