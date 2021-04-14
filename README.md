
<!-- README.md is generated from README.Rmd. Please edit that file -->

# SCOPfunctions

<!-- badges: start -->

<!-- badges: end -->

An R package of functions for single cell -omics analysis. Mostly
wrappers around Seurat and Bioconductor packages.

## Install

Install using devtools:

    devtools::install_github("CBMR-Single-Cell-Omics-Platform/SCOPfunctions")

Manually:

    git clone https://www.github.com/CBMR-Single-Cell-Omics-Platform/SCOPfunctions.git 

then from R:

    install.packages("./SCOPfunctions", type="source", repos=NULL)

## Usage

    library("SCOPfunctions")

### Preprocess

generate a seurat object from kallisto output

    seurat_obj <- SCOPfunctions::kallisto_to_seu(
                    spliced_dir,
                    spliced_name,
                    unspliced_dir,
                    unspliced_name)

Identify optimal cut-off quantile for demultiplexing hash tag oligos
with Seurat::HTOdemux()

    HTO_positive.threshold(
      seurat_obj,
      range=seq(from=0.7,to=1, by=0.01),
      assay = "hto",
      do_plot=TRUE
      ) 

identify and remove cell-cycle effects

    TODO

### Reduce dimensions

estimate the optimal number of dimensions to use for PCA and downstream

    estimate_dims(x=seurat_obj,
                  max_dim=75,
                  seed.use=12345)

### Plot

plot the distribution of cell clusters in different samples

    plot_barIdentGroup(seurat_obj,
                        var_ident="sample_ID",
                        var_group="cluster",
                        vec_group_colors=NULL,
                        f_color=colorRampPalette(brewer.pal(n=11, name="RdYlBu")),
                        do_plot = F)

![networkplot](./assets/barplot.png)

plot a cluster \* feature grid of gene expression violin plots

    # Here we just use the top variable genes, but normally we would use cluster marker genes
    plot_vlnGrid(seurat_obj,
                  slot="data",
                  var_group="cluster",
                  vec_features=head(VariableFeatures(seurat_obj),n=15),
                  vec_group_colors=NULL,
                  f_color = colorRampPalette(brewer.pal(n=11, name="RdYlBu"))

![networkplot](./assets/vlnplot.png)

make a network plot of a set of co-expressed features

    plot_network(
      mat_datExpr=as.matrix(GetAssayData(seurat_obj, slot="data")),
      vec_geneImportance=vec_geneImportance,
      vec_genes_highlight=c(),
      n_max_genes=50,
      igraph_algorithm = "drl",
      fontface_labels="bold.italic",
      color_edge = "grey70",
      edge_thickness = 1)

## Contribute

Issues and pull requests are welcome\! All contributions should be in
line with the [usethis code of
conduct](https://usethis.r-lib.org/CODE_OF_CONDUCT.html). This package
uses the methods and R tools set out in [R
packages](https://r-pkgs.org/intro.html). All Pull Requests should
follow the [tidyverse style
guide](https://style.tidyverse.org/documentation.html).
