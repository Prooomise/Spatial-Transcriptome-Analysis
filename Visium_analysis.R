###Library all packages
library(Seurat)
library(SeuratData)
library(ggplot2)
library(patchwork)
library(dplyr)
library(harmony)
library(purrr)
library(ggsci)


###Create Seurat object
p1 <- CreateSeuratObject(counts = Read10X("./GSM6833484/"), assay = "Spatial", project = "Patient1")
p1_spatial <- Read10X_Image(image.dir = file.path("./GSM6833484_spatial/"), 
                            image.name = "tissue_hires_image.png", slice = "slice1",
                            assay = "Spatial", filter.matrix = TRUE)
p1_spatial <- p1_spatial[Cells(x = p1)]
DefaultAssay(p1 = p1_spatial) <- "Spatial"
p1[["slice1"]] <- p1_spatial
#SpatialFeaturePlot(p1, features = "nFeature_Spatial", image.scale = "hires")+ theme(legend.position = "right")

p2 <- CreateSeuratObject(counts = Read10X("./GSM6833485/"), assay = "Spatial", project = "Patient2")
p2_spatial <- Read10X_Image(image.dir = file.path("./GSM6833485_spatial/"), 
                            image.name = "tissue_hires_image.png", slice = "slice2",
                            assay = "Spatial", filter.matrix = TRUE)
p2_spatial <- p2_spatial[Cells(x = p2)]
DefaultAssay(p2 = p2_spatial) <- "Spatial"
p2[["slice1"]] <- p2_spatial
#SpatialFeaturePlot(p2, features = "nFeature_Spatial")+ theme(legend.position = "right")

p3 <- CreateSeuratObject(counts = Read10X("./GSM6833486/"), assay = "Spatial", project = "Patient3")
p3_spatial <- Read10X_Image(image.dir = file.path("./GSM6833486_spatial/"), 
                            image.name = "tissue_hires_image.png", slice="slice3",
                            assay = "Spatial", filter.matrix = TRUE)
p3_spatial <- p3_spatial[Cells(x = p3)]
DefaultAssay(p3 = p3_spatial) <- "Spatial"
p3[["slice1"]] <- p3_spatial
#SpatialFeaturePlot(p3, features = "nFeature_Spatial")+ theme(legend.position = "right")

p4 <- CreateSeuratObject(counts = Read10X("./GSM6833487/"), assay = "Spatial", project = "Patient4")
p4_spatial <- Read10X_Image(image.dir = file.path("./GSM6833487_spatial/"), 
                            image.name = "tissue_hires_image.png", slice="slice4",
                            assay = "Spatial", filter.matrix = TRUE)
p4_spatial <- p4_spatial[Cells(x = p4)]
DefaultAssay(p4 = p4_spatial) <- "Spatial"
p4[["slice1"]] <- p4_spatial
#SpatialFeaturePlot(p4, features = "nFeature_Spatial")+ theme(legend.position = "right")

#Data quality control
p1[["percent.mt"]] <- PercentageFeatureSet(p1, pattern = "^MT-")
VlnPlot(p1, features = c("nFeature_Spatial", "nCount_Spatial", "percent.mt"), 
        ncol = 3, group.by = "orig.ident", pt.size = 0)
p1 <- SCTransform(p1, assay = "Spatial", verbose = FALSE)
#SpatialFeaturePlot(p1, features = c("KRT17"), alpha = c(0.1, 1))
p1 <- FindSpatiallyVariableFeatures(p1, assay = "SCT", features = VariableFeatures(p1)[1:1000],
                                    selection.method = "moransi")

p2[["percent.mt"]] <- PercentageFeatureSet(p2, pattern = "^MT-")
VlnPlot(p2, features = c("nFeature_Spatial", "nCount_Spatial", "percent.mt"), 
        ncol = 3, group.by = "orig.ident", pt.size = 0)
p2 <- SCTransform(p2, assay = "Spatial", verbose = FALSE)
#SpatialFeaturePlot(p2, features = c("KRT17"), alpha = c(0.1, 1))
p2 <- FindSpatiallyVariableFeatures(p2, assay = "SCT", features = VariableFeatures(p2)[1:1000],
                                    selection.method = "moransi")

p3[["percent.mt"]] <- PercentageFeatureSet(p3, pattern = "^MT-")
VlnPlot(p3, features = c("nFeature_Spatial", "nCount_Spatial", "percent.mt"), 
        ncol = 3, group.by = "orig.ident", pt.size = 0)
p3 <- SCTransform(p3, assay = "Spatial", verbose = FALSE)
#SpatialFeaturePlot(p3, features = c("KRT17"), alpha = c(0.1, 1))
p3 <- FindSpatiallyVariableFeatures(p3, assay = "SCT", features = VariableFeatures(p3)[1:1000],
                                    selection.method = "moransi")

p4[["percent.mt"]] <- PercentageFeatureSet(p4, pattern = "^MT-")
VlnPlot(p4, features = c("nFeature_Spatial", "nCount_Spatial", "percent.mt"), 
        ncol = 3, group.by = "orig.ident", pt.size = 0)
p4 <- SCTransform(p4, assay = "Spatial", verbose = FALSE)
#SpatialFeaturePlot(p4, features = c("KRT17"), alpha = c(0.1, 1))
p4 <- FindSpatiallyVariableFeatures(p4, assay = "SCT", features = VariableFeatures(p4)[1:1000],
                                    selection.method = "moransi")


'''
saveRDS(p1, "./p1.rds")
saveRDS(p2, "./p2.rds")
saveRDS(p3, "./p3.rds")
saveRDS(p4, "./p4.rds")
p1 <- readRDS("./p1.rds")
p2 <- readRDS("./p2.rds")
p3 <- readRDS("./p3.rds")
p4 <- readRDS("./p4.rds")
'''


### Merge Seurat object
OSCC_merge <- merge(p1, p2)
OSCC_merge <- merge(OSCC_merge, p3)
OSCC_merge <- merge(OSCC_merge, p4)
DefaultAssay(OSCC_merge) <- "SCT"
VariableFeatures(OSCC_merge) <- c(VariableFeatures(p1), VariableFeatures(p2),
                                  VariableFeatures(p3), VariableFeatures(p4))
#saveRDS(OSCC_merge, "./OSCC_merge.rds")


###PCA and Clustering
OSCC_merge <- RunPCA(OSCC_merge, verbose = FALSE)
OSCC_merge <- RunHarmony(OSCC_merge, dims = 1:50, "orig.ident", verbose = FALSE)

#ElbowPlot(OSCC_merge, ndims = 50)
#DimHeatmap(OSCC_merge, dims = 1:10, cells = 500, balanced = TRUE)
OSCC_merge <- FindNeighbors(OSCC_merge, dims = 1:50, reduction = "harmony")
OSCC_merge <- FindClusters(OSCC_merge, verbose = FALSE, resolution = 0.35) #0.35
OSCC_merge <- RunUMAP(OSCC_merge, reduction = "harmony", dims = 1:50)
DimPlot(OSCC_merge, reduction = "umap", group.by = c("ident", "orig.ident"))
SpatialDimPlot(OSCC_merge, image.alpha = 0.3, pt.size.factor = 2)
SpatialFeaturePlot(OSCC_merge, features = c("SPP1", "SLURP2"), image.alpha = 0.5, pt.size.factor = 2)

###Visualization cluster state in separate sample
OSCC_sample1 <- subset(OSCC_merge, subset = orig.ident == "Patient1")
OSCC_sample1 <- RunUMAP(OSCC_sample1, reduction = "harmony", dims = 1:30)
DimPlot(OSCC_sample1, reduction = "umap", group.by = "seurat_clusters", 
        label = TRUE, pt.size = 0.5) + ggtitle("Patient 1")


###Visualization particular gene spatial expression pattern
genelist <- c("JAG1", "JAG2", "NOTCH1", "NOTCH2", "CD44", "FN1")
pdf("./SpatialFeaturePlot_else.pdf")
map(genelist, ~{
  SpatialFeaturePlot(OSCC_merge, features = .x, image.alpha = 0, pt.size.factor = 2, stroke=NA)
})
dev.off()


###Calculate cell group ratio
Cellratio <- as.data.frame(prop.table(table(Idents(OSCC_merge), OSCC_merge$orig.ident), margin = 2))
ggplot(Cellratio) + 
  geom_bar(aes(x =Var2, y= Freq, fill = Var1),stat = "identity",width = 0.7,size = 0.5,colour = '#222222')+ 
  theme_classic() +
  labs(x='Sample',y = 'Ratio')+
  scale_fill_npg()+
  theme(panel.border = element_rect(fill=NA,color="black", size=0.5, linetype="solid"))

###Cell-type identification
OSCC_merge <- PrepSCTFindMarkers(OSCC_merge, assay = "SCT", verbose = TRUE)
all.markers.SCT <- FindAllMarkers(OSCC_merge, only.pos =T, min.pct = 0.5, logfc.threshold = 0.5)
table(all.markers.SCT$cluster)
#VlnPlot(subset(OSCC_merge, idents = "4"), features = c("percent.mt"), ncol = 1, group.by = "orig.ident")
#VlnPlot(OSCC_merge, features = c("percent.mt"), ncol = 1, group.by="seurat_clusters", pt.size = 0)

Genelist <- c("Flt1", "Pecam1", "Emcn", ### Vascular Endothelial cell
              "Ccl21a", "Prox1", "Lyve1", ### Lymphatic Endothelial cell
              "Wfdc2", "Tspan8", "Krt8", "Krt18", "Muc5b", ### Epithelial cell
              "Krt14", "Krt5", "Krt6a", "Cdh1", "Epcam", ### Epithelial cell
              "Cd14", "Fcer1g", "Csf3r", "Clec4e", "S100a8", ## Monocytes
              "Cd3g", "Cd3e", ## T cell
              "Cd74", ### Macrophage
              "Cd79a", "Igkc", "Ms4a1", ## B cell
              "Mki67", "Top2a", ## Cycling cell
              "Plp1", "Mpz", "Cdh19", "Csmd1", ### Neuron
              "Col1a2", "Col1a1", "Dcn", ## Fibroblast
              "Rgs5", "Acta2", "Myh11",
              "Hbb-bs", "Hba-a1", "Hbb-bt",
              "Pax7", "Myf5", "Ror1")
DotPlot(scRNA.merge.filter, features = Genelist, group.by = "seurat_clusters")+
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        legend.position = "top")+
  scale_color_gradient(low = "white", high = "red")


anno <- c("0"="Fibroblast", "1"="Epithelial Cell", "2"="Epithelial Cell", 
          "3"="Epithelial Cell", "4"="Glandular Epithelial Cell", "5" = "Smooth Muscle Cell", 
          "6" = "Myeloid Cell", "7" = "Epithelial Cell", "8" = "Vascular Endothelial Cell", 
          "9" = "Neuron", "10"="Myeloid Cell", "11"="T Cell", "12"="Muscle Satellite Cell",
          "13"="Neuron", "14"="Lymphatic Endothelial Cell", "15"="Fibroblast", 
          "16"="B Cell", "17"="Fibroblast", "18"="Erythrocyte","19"="T Cell", 
          "20"="Glandular Epithelial Cell", "21"="Epithelial Cell", "22"="Neuron", 
          "23"="Glandular Epithelial Cell", "24"="Epithelial Cell", "25"="Epithelial Cell")

scRNA.merge.filter[['cell_type']] <- unname(anno[scRNA.merge.filter@meta.data$seurat_clusters])








