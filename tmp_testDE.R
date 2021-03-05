# test 1
df_test1 <- MAST_DE_rand_effect_seurat(
  object=,
  random_effect.var = "HTO_classification",
  ident.1 = 1,
  group.by = "seurat_clusters"
)

ident.1="322-KO"
ident.2 = NULL
assay=NULL
slot="data"
group.by = "HTO_classification"
features = NULL
min.pct = 0.1
max.cells.per.ident = 200
random.seed = 1
latent.vars = NULL
random_effect.var = "seurat_clusters"
# test 1
ident.1="322-KO"
ident.2 = "328-KO"
assay="SCT"
slot="data"
group.by = "HTO_classification"
features = NULL
min.pct = 0.1
max.cells.per.ident = Inf
random.seed = 1
latent.vars = NULL
latent.vars = NULL
random_effect.var = "seurat_clusters"



ident.2 = NULL
logfc.threshold = 0.25
base = exp(1)
assay=NULL
slot="data"
features = NULL
min.pct = 0.1
max.cells.per.ident = NULL
random.seed = 1
latent.vars = NULL
n_cores=NULL
verbose=TRUE
p.adjust.method="fdr"

destfile = "pbmc3k_filtered_gene_bc_matrices.tar.gz"
download.file(url =  "https://cf.10xgenomics.com/samples/cell/pbmc3k/pbmc3k_filtered_gene_bc_matrices.tar.gz", destfile = destfile)
system2(command = "tar", args = paste0("-xzvf ",destfile))
# create Seurat object
pbmc_data <- Seurat::Read10X(data.dir = "filtered_gene_bc_matrices/hg19/")
pbmc <- Seurat::CreateSeuratObject(counts = pbmc_data)
# delete downloaded files
system2(command = "rm", args = c("-r", destfile,"filtered_gene_bc_matrices"))
# add artifical metadata. KO is knock-out, WT is wildtype
# simulate that we have five samples of each condition
pbmc$condition = factor(rep(c("KO","WT"),c(ncol(pbmc)/2,ncol(pbmc)/2)))
pbmc$sample = factor(rep(paste0("sample_",1:10), rep(ncol(pbmc)/10,10)))




## Examples

```{r example}

install.
library(SCOPfunctions)
library(SeuratData)

```

Use human immune cells (PBMC) dataset from the Seurat Introduction to scRNA-seq integration Vignette


```{r example}
# https://satijalab.org/seurat/articles/integration_introduction.html

# install dataset
Seurat::InstallData("ifnb")

# load dataset
LoadData("ifnb")

# split the dataset into a list of two seurat objects (stim and CTRL)
ifnb.list <- SplitObject(ifnb, split.by = "stim")

# download data to current working directory

```

Run differential expression MAST test with sample as random effect

```{r example}
pbmc = Seurat::NormalizeData(object = pbmc) # using NormalizeData for speed. Better to use Seurat::SCTransform.
# Note how arguments are almost identical to those of Seurat::FindMarkers
df_DE <- SCOPfunctions::DE_MAST_RE_seurat(
  object = pbmc,
  random_effect.var = "sample", # NB: required!
  ident.1 = "KO",
  group.by = "condition",
  logfc.threshold = 0.1,
  max.cells.per.ident = 2000,)
```
