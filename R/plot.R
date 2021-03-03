#' Generate barplot of identity group composition
#'
#' Generate a percentage barplot that shows the composition of each identity
#' (e.g. sample)in terms of groups (e.g. cell types)
#'
#' @param seurat_obj Seurat object (Seurat ^3.0)
#' @param var_ident the identify variable, character
#' @param var_group the group variable, character
#' @param vec_group_colors a vector of colors, named by corresponding group. Length must match number of groups. Character
#' @param f_color if vec_group_colors is not provided, the user may instead provide a function f_color() that takes as its only argument the number of colors to generate
#' @param do_plot Whether to plot, logical
#' @param title NULL to leave out
#' @param fontsize_title NULL to leave out
#' @param fontsize_axistitle_x NULL to leave out
#' @param fontsize_axistitle_y NULL to leave out
#' @param fontsize_axistext_x NULL to leave out
#' @param fontsize_axistext_y NULL to leave out
#' @param fontsize_legendtitle NULL to leave out
#' @param fontsize_legendtext NULL to leave out
#'
#' @return a ggplot object
#' @export TODO
#'
#' @examples TODO
f_plot_bar_ident_group_composition = function(seurat_obj,
                                          var_ident,
                                          var_group,
                                          vec_group_colors=NULL,
                                          f_color=colorRampPalette(brewer.pal(n=11, name="RdYlBu")),
                                          do_plot = F,
                                          title = NULL,
                                          fontsize_title = 24,
                                          fontsize_axistitle_x = 18,
                                          fontsize_axistitle_y = 18,
                                          fontsize_axistext_x = 12,
                                          fontsize_axistext_y = 12,
                                          fontsize_legendtitle = 12,
                                          fontsize_legendtext = 10,
                                          aspect.ratio=1.2) {

  stopifnot(require("data.table"))
  stopifnot(require("ggplot2"))

  #===============data.table with sums==================
  dt = data.table("ident" = as.character(seurat_obj@meta.data[[var_ident]]),
                  "group" = as.character(seurat_obj@meta.data[[var_group]]))
  dt[,n_ident := paste0(ident," (n=",.N, ")"), by=ident]
  vec_factorLevels <- dt$n_ident[gsub("\\ .*","",dt$n_ident) %>% as.numeric %>% order] %>% unique
  dt[,n_ident := factor(n_ident, levels = vec_factorLevels, ordered=T),]
  dt_sum <- dt[,.N, by=.(n_ident,group)]

  #===============ggplot==================
  # colors
  if (is.null(vec_group_colors)) {
    n_group <- length(unique(dt$group))
    vec_group_colors <- f_color(n_group)
    names(vec_group_colors) <- unique(dt$group)
  }

  p <- ggplot(dt_sum,
              aes(x = n_ident, y=N, fill = factor(group))) +

    geom_bar(
      position="fill",
      stat="identity",
      width=0.6,
      show.legend = if (!is.null(fontsize_legendtext)) TRUE else FALSE
      #position=position_dodge()
    ) +

    scale_y_continuous(labels = scales::percent) +
    scale_fill_manual(values=vec_group_colors) +

    theme(
      axis.title.x = if (is.null(fontsize_axistitle_x)) element_blank() else element_text(size=fontsize_axistitle_x, vjust=0),
      axis.text.x = if (is.null(fontsize_axistext_x)) element_blank() else element_text(angle = 90, size=fontsize_axistext_x,vjust=0.5),
      axis.title.y = if (is.null(fontsize_axistitle_y)) element_blank() else element_text(size=fontsize_axistitle_y),
      axis.text.y = if (is.null(fontsize_axistext_y)) element_blank() else  element_text(size=fontsize_axistext_y),
      legend.title = if (is.null(fontsize_legendtext)) element_blank() else element_text(size=fontsize_legendtitle),
      legend.text = if (is.null(fontsize_legendtext)) element_blank() else element_text(size = fontsize_legendtext),
      legend.background = element_blank(),
      legend.box.background = element_blank(),
      plot.background=element_blank(),
      aspect.ratio = aspect.ratio) +

    labs(x=var_ident, y="proportion", fill = var_group)

  if (do_plot) p

  return(p)
}


#' violin plot grid
#'
#' produce a n_cluster * n_feature grid of violin plots
#'
#' @param seurat_obj Seurat object (Seurat ^3.0)
#' @param assay  seurat_obj assay to use
#' @param var_group  the group variable, character
#' @param slot  seurat_obj slot to use
#' @param vec_features  a vector of features to plot in the violin plot
#' @param vec_group_colors  a vector of colors, named by corresponding group. Length must match number of groups. Character
#' @param f_color  if vec_group_colors is not provided, the user may instead provide a function f_color() that takes as its only argument the number of colors to generate
#' @param flip  if TRUE (default), groups are rows and features are columns, and vice-verso for FALSE
#' @param do_plot  Whether to plot, logical
#' @param pt.size size of jitter in the violin plots. Set to 0 (default) to omit
#'
#' @return ggplot object
#' @export
#'
#' @examples
f_plot_vln_grid = function(seurat_obj,
                            assay,
                            slot,
                            var_group,
                            vec_features,
                            vec_group_colors,
                            f_color = colorRampPalette(brewer.pal(n=11, name="RdYlBu")),
                            flip = T,
                            do_plot = F,
                            pt.size = 0,
                            feature_fontface = "bold.italic",
                            fontsize_axistext_x=12,
                            fontsize_axistext_y=12,
                            aspect.ratio =NULL
                           ) {

  stopifnot(require("data.table"))
  stopifnot(require("Seurat"))
  stopifnot(require("patchwork"))


  #=============prepare group and colors==================
  seurat_obj_tmp = seurat_obj
  Idents(seurat_obj_tmp) <- var_group
  levels(x = seurat_obj_tmp) = sort(unique(seurat_obj_tmp@meta.data[[var_group]]), decreasing = if (flip) T else F)

  if (is.null(vec_group_colors)) {
    n_group <- length(levels(x = seurat_obj_tmp))
    vec_group_colors <- f_color(n_group)
    names(vec_group_colors) <- levels(x = seurat_obj_tmp)
  }

  #=============generate plot list==================
  # produces a list of rows of violin plots, one per feature
  list_plot <- VlnPlot(object=seurat_obj_tmp,
                       assay=assay,
                       features = vec_features,
                       pt.size = pt.size,
                       cols=vec_group_colors,
                       sort = F,
                       #group.by = var_group,
                       same.y.lims = F,
                       slot=slot,
                       log = F,
                       combine = F,
                       flip=F)

  names(list_plot) <- vec_features

  if (is.null(aspect.ratio)) {
    aspect.ratio = 1.5*length(vec_group_colors)/length(vec_features)
    message(paste0("Using aspect ratio ", aspect.ratio))
  }

  list_plot_flip <- lapply(1:length(list_plot), function(i) {
    plot_tmp=list_plot[[i]]
    if (flip) {
      plot_tmp = plot_tmp +
        coord_flip()
    }
    plot_tmp = plot_tmp  +
      theme(
        plot.title = element_text(face = feature_fontface, size=fontsize_axistext_y),
        axis.text.y=element_blank(),
        plot.margin = margin(b=1, unit="cm"))
    #margin around entire plot (unit with the sizes of the top, right, bottom, and left margins)

    if (i==1) {
      plot_tmp <- plot_tmp +
        theme(
          axis.text.y=element_text(
            hjust=1,
            vjust=0,
            size=fontsize_axistext_x,
            angle=30
          )
        )
    }


    plot_tmp <- plot_tmp +
      theme(
        axis.title.y=element_blank(),
        axis.title.x=element_blank(),
        axis.text.x = element_text(
          size= fontsize_axistext_x,
          angle=30),
        aspect.ratio=aspect.ratio,
        legend.position="none")
  })

  p <- patchwork::wrap_plots(... = list_plot_flip,
                               ncol = length(list_plot_flip)
  )

  if (do_plot) p

  return(p)
}
