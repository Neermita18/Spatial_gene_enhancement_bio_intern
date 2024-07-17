library("ENGEP")

library("propr")

library(Seurat)

spatial_df<- read.csv('/home/pylab/Desktop/swainasish/Spatial_enhance/datasets/spage_datasets/Spatial/Starmap/STARmap_data_scvi.csv', row.names=1)




spatial_df=t(spatial_df)

spatial_df

saveRDS(spatial_df, file = "/home/pylab/Desktop/swainasish/Spatial_enhance/datasets/spage_datasets/Spatial/Starmap/vispspatial.rds")

sp=readRDS("/home/pylab/Desktop/swainasish/Spatial_enhance/datasets/spage_datasets/Spatial/Starmap/vispspatial.rds")

sp


dim(sp)


scrna_df=read.csv("/home/pylab/Desktop/swainasish/Spatial_enhance/datasets/spage_datasets/scRNAseq/Allen_VISp/mouse_VISp_2018-06-14_exon-matrix.csv")

scrna_df

genes=read.csv("/home/pylab/Desktop/swainasish/Spatial_enhance/datasets/spage_datasets/scRNAseq/Allen_VISp/mouse_VISp_2018-06-14_genes-rows.csv")

genes

rownames(scrna_df) <- genes$gene_symbol

scrna_df <- scrna_df[, -1]

scrna_df

saveRDS(scrna_df, file = "/home/pylab/Desktop/swainasish/Spatial_enhance/datasets/spage_datasets/scRNAseq/Allen_VISp/visprds.rds")

sc=readRDS("/home/pylab/Desktop/swainasish/Spatial_enhance/datasets/spage_datasets/scRNAseq/Allen_VISp/visprds.rds")

sc

scs <- CreateSeuratObject(sc,project = "sc")
scs <- NormalizeData(object = scs,normalization.method = "LogNormalize")



scs

new <- FindVariableFeatures(scs, nfeatures = 2000)
hvg_visp <- VariableFeatures(new)

str(hvg_visp)

install.packages("Seurat")

library(Seurat)
