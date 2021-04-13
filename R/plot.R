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
#' @return a ggplot2 object
#' @export
#'
#' @examples p <- plot_barIdentGroup(seurat_obj=seu, var_ident="sample",var_group="cluster")
plot_barIdentGroup = function(seurat_obj,
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


#' Make grid of violin plots
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
#' @return ggplot2 object
#' @export
#'
#' @examples p <- plot_vlnGrid(seurat_obj=seu, assay="RNA", slot="data", var_group="cluster", vec_features=head(VariableFeatures(seu)))
plot_vlnGrid = function(seurat_obj,
                        assay,
                        slot,
                        var_group,
                        vec_features,
                        vec_group_colors=NULL,
                        f_color = colorRampPalette(brewer.pal(n=11, name="RdYlBu")),
                        flip = T,
                        do_plot = F,
                        pt.size = 0,
                        feature_fontface = "bold.italic",
                        fontsize_axistext_x=12,
                        fontsize_axistext_y=12,
                        aspect.ratio =NULL
                       ) {

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


plot_network = function(
  mat_datExpr,
  vec_geneImportance,
  vec_genes_highlight=c(),
  n_max_genes=50,
  igraph_algorithm = "drl",
  fontface_labels="bold.italic",
  color_edge = "grey70",
  edge_thickness = 1
) {

  # sort
  vec_geneImportance = sort(vec_geneImportance, decreasing = T)

  # take take hub genes and genes to highlight (if any)
  vec_logical = names(vec_geneImportance) %in% c(names(vec_geneImportance)[1:min(n_max_genes, length(vec_geneImportance))],  vec_genes_highlight)
  vec_geneImportance = vec_geneImportance[vec_logical]

  mat_datExpr <- mat_datExpr[names(vec_geneImportance),]

  # Compute network adjacency
  mat_adj <- WGCNA::adjacency(t(mat_datExpr),
                              power = 1,
                              corFnc = "cor",
                              corOptions = list(use = "p"),
                              type = "signed hybrid")

  df_hub.data <- cbind.data.frame(gene = names(vec_geneImportance), importance = vec_geneImportance)

  # hub.data$logfc <- DEGs$log2FoldChange[match(hub.data$gene, rownames(DEGs))]
  # hub.data$p.adj <- DEGs$padj[match(hub.data$gene, rownames(DEGs))]
  mat_adj %>%
    graph.adjacency(mode = "undirected", weighted = T, diag = FALSE) %>%
    tidygraph::as_tbl_graph() %>% upgrade_graph() %>% activate(nodes) %>%
    dplyr::mutate(importance = df_hub.data$importance*edge_thickness) %>%
    activate(edges) %>%
    activate(nodes) %>%
    filter(!node_is_isolated()) -> hub.plot

  module.plot <- ggraph(hub.plot,
                        layout = "igraph",
                        algorithm = igraph_algorithm,
                        #layout = "stress"
  ) +
    geom_edge_link(color=col_edge, show.legend=F, aes(alpha=weight)) +
    geom_node_point(aes(size = weight), fill = "#73BCC9", shape=21, alpha=0.8) +
    # scale_size(breaks = c(1,2,3),
    #            limits = c(0,4),
    #            range = c(0, 8),
    #            labels = c(" \U00B1 1", " \U00B1 2", " \U00B1 3")
    #            ) +
    geom_node_text(
      aes(label = as.character(name)),
      fontface="bold.italic",
      size=fontSize_label_lg,
      repel = T,
      # incredibly, there is a mismatch in how the hyphon "-" is encoded in the genes..
      color = ifelse(gsub("\\W","_",rownames(mat_adj)) %in% gsub("\\W","_",c(vec_genes_NAFLD_hep1,vec_genes_NAFLD_mac1)), "red", "black") ) +

    guides(
      size = guide_legend(override.aes = list(
        #size=c(2,4,6)),
      ),
      keywidth = 0.8,
      keyheight = 0.8, order = 1)) +
    #title = expression(bold(paste("Log"[2],~
    #"fold-change"))))) +
    theme_graph(base_family = 'Helvetica') +
    theme(
      legend.title.align=0.5,
      legend.position = "top",
      legend.margin = margin(0,0,0,0, unit="cm"),
      legend.title = element_text(size=fontSize_legend_xlg, face="bold"),
      legend.text = element_text(size=fontSize_legend_lg, face="bold"),
      legend.spacing.x = unit(0, "cm"),

      axis.text =element_blank(),
      axis.ticks = element_blank(),
      axis.line=element_blank(),
      axis.title = element_blank(),
      # margin: top, right, bottom, and left
      plot.margin = unit(c(0, 0, 0, 0), "cm"),

      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_blank(),

      plot.background = element_blank(),
      panel.grid = element_blank()
    ) +

    coord_cartesian(clip="off")

  if (do_plot) p
  return(p)

}
