library(Seurat)
library(cowplot)
library(ggplot2)
library(dplyr)
library(tidyverse)
library(viridis)
library(RColorBrewer)
library(pheatmap)
library(scales)

####################################
# makedirs
####################################
makedirs <- function(dir_list){
  for (i_dir in dir_list){
    if (!dir.exists(i_dir) ){
      dir.create(i_dir, recursive = TRUE)
    }
  }
}

project_work_dir <- "/home/project/Single_cells_project/scrna_project/Guanjian_GJRE_GJLE/20210607_analysis"
makedirs(project_work_dir)
setwd(project_work_dir)

## load data
project_work_data_dir <- "/home/project/Single_cells_project/scrna_project/Guanjian_GJRE_GJLE/Guanjian_GJRE_GJLE_rna_seq_sct_report/Rdata"
data.filt <- get(load(file.path(project_work_data_dir, 'Cluster_all.Rdata')))
head(data.filt@meta.data)

Idents(data.filt) <- "cluster"
table(Idents(data.filt))
# 0    1    2    3    4    5    6    7    8    9   10   11   12   13   14   15 
# 1603 1153 1096  963  951  901  800  795  776  675  655  604  548  482  316  287 
# 16   17   18   19   20   21   22   23   24   25   26 
# 263  151  142  117  115  113  113  105   80   56   37


#########################
## 01_Tcells_TSNE_plot
#########################
makedirs(file.path(project_work_dir, "01_Tcells_TSNE_plot"))
setwd(file.path(project_work_dir, "01_Tcells_TSNE_plot"))

alpha.use <- 0.8
p1 <- DimPlot(object = data.filt, reduction = "tsne", label = TRUE, label.size = 5, pt.size=0.1, raster=FALSE) 
p1$layers[[1]]$mapping$alpha <- alpha.use
p2 <- DimPlot(object = data.filt, reduction = "umap", label = TRUE, label.size = 5, pt.size=0.1, raster=FALSE)
p2$layers[[1]]$mapping$alpha <- alpha.use
p1<- p1 + scale_alpha_continuous(range = alpha.use, guide = F)
p2<- p2 + scale_alpha_continuous(range = alpha.use, guide = F)
p <- plot_grid(p1, p2)
save_plot("Compare_TSNE_and_UMAP.png", p, base_height = 8, base_aspect_ratio = 2.5, base_width = NULL, dpi=600)
save_plot("Compare_TSNE_and_UMAP.pdf", p, base_height = 8, base_aspect_ratio = 2.5, base_width = NULL)
save_plot("TSNE_clustering.png", p1, base_height = 8, base_aspect_ratio = 1.3, base_width = NULL, dpi=600)
save_plot("TSNE_clustering.pdf", p1, base_height = 8, base_aspect_ratio = 1.3, base_width = NULL)
save_plot("UMAP_clustering.png", p2, base_height = 8, base_aspect_ratio = 1.3, base_width = NULL, dpi=600)
save_plot("UMAP_clustering.pdf", p2, base_height = 8, base_aspect_ratio = 1.3, base_width = NULL)


#######################################
## 02_gene_exprssion
## 做出以下marker gene的feature plot 和小提琴图
#######################################
makedirs(file.path(project_work_dir, "02_gene_exprssion"))
setwd(file.path(project_work_dir, "02_gene_exprssion"))


gene_list <- c("CD79A","CD79B","MS4A1","MZB1","CD27","CD38","IGHD","IGHM","IGHA1","IGHA2", "IGHE", "NEIL1","MKI67")

gene_list <- c("KRT14")

# DotPlot
DefaultAssay(data.filt) <- "RNA"
Idents(data.filt) <- "cluster"
plot1 <- DotPlot(data.filt,  assay = "RNA", features = rev(gene_list), cols = "RdBu") + RotatedAxis()+ coord_flip()+theme(panel.border = element_blank()) +
  ylab("Cluster") +
  theme(axis.text.x = element_text(
    angle = 0,
    hjust = 0.5,
    vjust = 0.5
  ))
save_plot("dotplot_marker_genes_expression.png", plot1, base_height = 8, base_aspect_ratio = 1.5, base_width = NULL, dpi=600)
save_plot("dotplot_marker_genes_expression.pdf", plot1, base_height = 8, base_aspect_ratio = 1.5, base_width = NULL)

  DefaultAssay(data.filt) <- "RNA"
  p2 <- DotPlot(data.filt, features = top_N_gene_names, cols = "RdYlBu") + RotatedAxis()
  
  
# VlnPlot
vlnplot_p1 <- VlnPlot(object = data.filt, features = gene_list, pt.size=0.1, ncol = 2)
save_plot("VlnPlot_marker_genes_expression.png", vlnplot_p1, base_height = 16, base_width = 18, dpi=600, limitsize = FALSE)
save_plot("VlnPlot_marker_genes_expression.pdf", vlnplot_p1, base_height = 16, base_width = 18, limitsize = FALSE)

vlnplot_p1 <- VlnPlot(object = data.filt, features = gene_list, split.by = 'orig.ident', pt.size=0.1, ncol = 2)
save_plot("VlnPlot_marker_genes_expression_split_sample.png", vlnplot_p1, base_height = 16, base_width = 30, dpi=600, limitsize = FALSE)
save_plot("VlnPlot_marker_genes_expression_split_sample.pdf", vlnplot_p1, base_height = 16, base_width = 30, limitsize = FALSE)

DefaultAssay(data.filt) <- "RNA"
Idents(data.filt) <- 'orig.ident'
stackVlnPlot <- VlnPlot(data.filt, gene_list, stack = TRUE, sort = F, flip = T, group.by='orig.ident') +
  theme(legend.position = "none")
save_plot("stackVlnPlot_marker_genes_expression.png", stackVlnPlot, base_height = 8, base_width = 8, dpi=600)
save_plot("stackVlnPlot_marker_genes_expression.pdf", stackVlnPlot, base_height = 8, base_width = 8)

stackVlnPlot <- VlnPlot(data.filt, gene_list, stack = TRUE, sort = F, flip = T, group.by='orig.ident') 
save_plot("stackVlnPlot_marker_genes_expression_among_Bcells1.png", stackVlnPlot, base_height = 8, base_width = 8, dpi=600)
save_plot("stackVlnPlot_marker_genes_expression_among_Bcells1.pdf", stackVlnPlot, base_height = 8, base_width = 8, dpi=600)



# FeaturePlot
DefaultAssay(data.data.filt) <- "RNA"
p5 <- FeaturePlot(data.filt, gene_list, reduction = "tsne", ncol = 5, cols = viridis(100),  order = T)
save_plot('FeaturePlot_genes_expression_among_Bcells.png', p5, base_height = 12, base_width = 20, dpi=600, limitsize = FALSE)
save_plot('FeaturePlot_genes_expression_among_Bcells.pdf', p5, base_height = 12, base_width = 20, limitsize = FALSE)

DefaultAssay(data.data.filt) <- "RNA"
p5 <- FeaturePlot(data.filt, gene_list, reduction = "tsne", ncol = 5, cols = c("lightgrey", "red"),  order = T)
save_plot('FeaturePlot_genes_expression_among_Bcells_1.png', p5, base_height = 12, base_width = 20, dpi=600, limitsize = FALSE)
save_plot('FeaturePlot_genes_expression_among_Bcells_1.pdf', p5, base_height = 12, base_width = 20, limitsize = FALSE)






