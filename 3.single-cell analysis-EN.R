#####脑组织#####
setwd("D:/AD_GWAS/AD脑组织/GSE157827")

dir.create("~/SeuratV4")
# 然后安装的时候，指定安装目录
install.packages('Seurat', repos = c('https://satijalab.r-universe.dev'), lib = "~/SeuratV4")
.libPaths(c("~/SeuratV4", .libPaths()))

library(limma)
library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)
#BiocManager::install("SingleR")
library(SingleR)
library(CCA)
library(clustree)
library(cowplot)
#if (!require("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install("HSMMSingleCell")
library(monocle)
library(tidyverse)


dir_name <- list.files("D:/AD_GWAS/AD脑组织/GSE157827")
dir_name
####Read10X(data.dir绝对路径读取####
dir_name=list.files("D:/AD_GWAS/AD脑组织/GSE157827/")
scRNAlist <- list()
for(i in 1:length(dir_name)){
  counts <- Read10X(data.dir = paste("D:/AD_GWAS/AD脑组织/GSE157827/", dir_name[i], sep = ""))
  scRNAlist[[i]] <- CreateSeuratObject(counts, project = dir_name[i],min.cells = 3, min.features = 300)
}

####批量计算线粒体和红细胞比例####
for(i in 1:length(scRNAlist)){
  sc <- scRNAlist[[i]]#获取scRNAlist中的第i个Seurat对象
  # 计算线粒体比例
  sc[["mt_percent"]] <- PercentageFeatureSet(sc, pattern = "^MT-")#计算以"MT-"开头的基因比例
  # 计算红细胞比例
  HB_genes <- c("HBA1","HBA2","HBB","HBD","HBE1","HBG1","HBG2","HBM","HBQ1","HBZ") #定义红细胞基因列表
  HB_m <- match(HB_genes, rownames(sc@assays$RNA))#在Seurat对象的RNA数据中查找红细胞基因的索引位置
  HB_genes <- rownames(sc@assays$RNA)[HB_m]  # 获取匹配到的红细胞基因的行名
  HB_genes <- HB_genes[!is.na(HB_genes)]  # 删除NA值（未匹配到的基因）
  sc[["HB_percent"]] <- PercentageFeatureSet(sc, features=HB_genes)   #计算红细胞基因比例，将结果存储在名为"HB_percent"的新列中
  # 将sc赋值给scRNAlist[[i]]
  scRNAlist[[i]] <- sc
  # 删除sc
  rm(sc)
}

####批量绘制质控前小提琴图####
violin_before <- list()
for(i in 1:length(scRNAlist)){
  violin_before[[i]] <- VlnPlot(scRNAlist[[i]],
                                features = c("nFeature_RNA", "nCount_RNA", "mt_percent","HB_percent"),
                                pt.size = 0.01,
                                ncol = 4)
}
# 合并图片
violin_before_merge <- CombinePlots(plots = violin_before,nrow=length(scRNAlist),legend='none')
# 将图片输出到画板上
violin_before_merge
# 保存图片
ggsave("violin_before_merge.pdf", plot = violin_before_merge, width = 10, height = 50, limitsize = FALSE)

violin_before



####批量过滤细胞、MT、HB基因####
scRNAlist <- lapply(X = scRNAlist, FUN = function(x){
  x <- subset(x,
              subset = nFeature_RNA > 200 & nFeature_RNA < 5000 &
                mt_percent < 10 &
                HB_percent < 3 &
                nCount_RNA < quantile(nCount_RNA,0.97) &
                nCount_RNA > 1000)})


####merge合并样本####
scRNAlist <- merge(x=scRNAlist[[1]],y=scRNAlist[-1])

violin_after <- VlnPlot(scRNAlist,
                        features = c("nFeature_RNA", "nCount_RNA", "mt_percent","HB_percent"),
                        pt.size = 0.01,
                        ncol = 4)
# 将图片输出到画板上
violin_after
# 保存图片
ggsave("vlnplot_after_qc.pdf", plot = violin_after, width = 50, height =15, limitsize = FALSE)

####数据归一化、筛选高变基因与PCA降维####
# harmony整合是基于PCA降维结果进行的。
scRNAlist <- NormalizeData(scRNAlist) %>% # 数据归一化处理
  FindVariableFeatures(selection.method = "vst",nfeatures = 3000) %>% #筛选高变基因
  ScaleData() %>% #数据标准化
  RunPCA(npcs = 30, verbose = T)#npcs：计算和存储的PC数（默认为 50）
a=DimPlot(scRNAlist,reduction = "pca",group.by = "orig.ident")

top15 <- head(VariableFeatures(scRNAlist), 15)
plot1 <- VariableFeaturePlot(scRNAlist)
plot2 <- LabelPoints(plot = plot1, points = top15, repel = TRUE, size=3)
# 合并图片
feat_15 <- CombinePlots(plots = list(plot1,plot2),legend = "bottom")
feat_15
# 保存图片
ggsave(file = "feat_15.pdf",plot = feat_15,he = 10,wi = 15 )


library(harmony)
####RunHarmony去批次####
# 整合需要指定Seurat对象和metadata中需要整合的变量名。
scRNA_harmony <- RunHarmony(scRNAlist, group.by.vars = "orig.ident")
scRNA_harmony@reductions[["harmony"]][[1:5,1:5]]
b=DimPlot(scRNA_harmony,reduction = "harmony",group.by = "orig.ident")
#PCA图看到还有一点的批次效应（融合得比较好批次就弱）
b

pca_harmony_integrated <- CombinePlots(list(a,b),ncol=1)
pca_harmony_integrated

#去批次前后对比
ggsave(file = "RunHarmony去批次.pdf",plot = pca_harmony_integrated,he = 15,wi = 12 )
table(scRNA_harmony@meta.data$orig.ident)

scRNA_harmony$group <- NA  # Initialize group column in metadata

# Assign groups based on suffix
scRNA_harmony$group[grep("_AD", scRNA_harmony$orig.ident)] <- "AD"
scRNA_harmony$group[grep("_NC", scRNA_harmony$orig.ident)] <- "NC"

# Check the results
table(scRNA_harmony$group)


# 保存结果（如需保存映射后的 Seurat 对象）
saveRDS(scRNA_harmony, file = "scRNA_harmony.rds")

scRNA_harmony=readRDS('./scRNA_harmony.rds')##这个是T细胞类型的T_celltype,tissue_type,celltype
setwd("D:/AD_GWAS/AD脑组织/GSE157827/rest")
# 质控图按新分组绘制
pdf("分组质控.pdf")
VlnPlot(scRNA_harmony, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.rb"),
        group.by = "group", pt.size = 0.1, ncol = 2) +
  theme(legend.position = "right")
dev.off()


# 后续都是基于Harmony矫正之后的数据，不是基因表达数据和直接的PCA降维数据。
# 设置reduction = 'harmony'，后续分析是基于Harmony矫正之后的数据。
####聚类、umap/tsne降维降维####
ElbowPlot=ElbowPlot(scRNA_harmony, ndims=50, reduction="pca")
ElbowPlot(scRNA_harmony, ndims=50, reduction="harmony")

ElbowPlot
ggsave( ElbowPlot, filename = "ElbowPlot.pdf", width = 7, height = 5 )

scRNA_harmony <- FindNeighbors(scRNA_harmony, reduction = "harmony", dims = 1:20) %>% FindClusters(resolution = 0.8)
table(scRNA_harmony@meta.data$seurat_clusters)
table(scRNA_harmony@meta.data$seurat_clusters)
total_cells <- sum(table(scRNA_harmony@meta.data$seurat_clusters))
table(total_cells)

scRNA_harmony <- RunTSNE(scRNA_harmony, reduction = "harmony", dims = 1:20)
scRNA_harmony <- RunUMAP(scRNA_harmony, reduction = "harmony", dims = 1:20)
#scRNA_harmony=sce
# 绘图
umap_integrated1 <- DimPlot(scRNA_harmony, reduction = "umap", group.by = "group")
umap_integrated2 <- DimPlot(scRNA_harmony, reduction = "umap", label = TRUE)
tsne_integrated1 <- DimPlot(scRNA_harmony, reduction = "tsne", group.by = "group")
tsne_integrated2 <- DimPlot(scRNA_harmony, reduction = "tsne", label = TRUE)
# 合并图片
umap_tsne_integrated <- CombinePlots(list(tsne_integrated1,tsne_integrated2,umap_integrated1,umap_integrated2),ncol=2)
# 将图片输出到画板
umap_tsne_integrated
# 保存图片
ggsave("36-umap_tsne_integrated2.pdf",umap_tsne_integrated,wi=20,he=15)
#保存数据
save(scRNA_harmony,file = "35-cluster-scdata.Rdata")

# 差异分析
markers <- FindAllMarkers(object = scRNA_harmony, test.use="wilcox" ,
                          only.pos = TRUE,
                          logfc.threshold = 0.25)


#2.对计算好的每cluster的marker基因进行筛选
all.markers =markers %>% dplyr::select(gene, everything()) %>% subset(p_val_adj<0.05)
#筛选出P<0.05的marker基因
top10 = all.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC) #将每个cluster lgFC排在前10的marker基因挑选出来
top15 = all.markers %>% group_by(cluster) %>% top_n(n = 15, wt = avg_log2FC)
top20 = all.markers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_log2FC)
top30 = all.markers %>% group_by(cluster) %>% top_n(n = 30, wt = avg_log2FC)


View(top10)
write.csv(top10,"12cluster_top10.csv",row.names = T)
write.csv(top15,"12cluster_top15.csv",row.names = T)
write.csv(top20,"12cluster_top20.csv",row.names = T)
write.csv(top30,"12cluster_top30.csv",row.names = T)
write.csv(all.markers,"12cluster_all.markers.csv",row.names = T)


markers = c('PECAM1','VWF',#endothelial
            'DCN','COL1A2',#Fibroblasts
            "TYROBP","CTSS","C1QA","HEXB","AIF1",# Micro/Macro
            'C1QA','C1QB','APOE',#Macrophage
            'ALDH1L1','AQP4','GFAP','SLC1A2',#Astrocyte
            'MBP','MOBP','MOG','PLP1',#Oligo
            'RBFOX3','MAP2','GAD1',#Neurons
            'GAD1','GAD2','PVALB','CALB2','VIP',#Inhibitory Neurons
            'SLC17A7','CAMK2A','SATB2',#Excitatory Neurons
            'CSPG4','PDGFRA','GPR17'#opc
)

missing_genes <- markers[!markers %in% rownames(scRNA)]
print(missing_genes)
valid_markers <- markers[markers %in% rownames(scRNA)]
unique_markers <- unique(valid_markers)
DotPlot(scRNA, features = unique_markers) +
  coord_flip() +
  theme_bw() +
  theme(panel.grid = element_blank(), axis.text.x = element_text(hjust = 1, vjust = 0.5)) +
  labs(x = NULL, y = NULL) + guides(size = guide_legend(order = 3)) +
  theme(axis.text.x = element_text(angle = 90)) +
  scale_color_gradientn(values = seq(0, 1, 0.2), colours = c('#008D83', '#B0A8B9', '#FF8066', '#C34A36'))


ggsave(filename="7-cluster-marker4.pdf", width = 110, height = 170, units = "mm")
####按照自己想的细胞顺序进行排序####
# 将 seurat_clusters 转换为因子并按照指定顺序排序
scRNA_harmony$seurat_clusters <- factor(scRNA_harmony$seurat_clusters,
                                        levels = c("24","17","12","29","32","34","22","1","3","6","16","21",
                                                   "26","27","31","30","18","13","10","5","0","19","9",
                                                   "11","14","2","4","8","15","20","23","25","28","7",
                                                   "35","33"
                                        ))

# 确保 markers 都存在于数据中
missing_genes <- markers[!markers %in% rownames(scRNA_harmony)]
print(missing_genes)  # 检查缺失基因
valid_markers <- markers[markers %in% rownames(scRNA_harmony)]
unique_markers <- unique(valid_markers)  # 去除重复基因

# 绘制 DotPlot
dotplot=DotPlot(scRNA_harmony, features = unique_markers, group.by = "seurat_clusters") +
  coord_flip() +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(hjust = 1, vjust = 0.5, angle = 90)) +
  labs(x = NULL, y = NULL) +
  guides(size = guide_legend(order = 3)) +
  scale_color_gradientn(values = seq(0, 1, 0.2), colours = c('#008D83', '#B0A8B9', '#FF8066', '#C34A36'))

# Display the plot
dotplot

ggsave(filename="36-cluster-marker24-12-31_updated-4.pdf", width = 240, height = 170, units = "mm")
ggsave(filename = "36-cluster-marker_updated-4.pdf", plot = dotplot, width = 240, height = 140, units = "mm")


#####cluster命名#####
celltype=data.frame(ClusterID=0:35,
                    celltype='unkown')
celltype[celltype$ClusterID %in% c(9,11,14),2]='Inhibitory Neurons'
celltype[celltype$ClusterID %in% c(19,2,4,8,15,20,23,25),2]='Excitatory Neurons'
celltype[celltype$ClusterID %in% c(22,1,3,6,16,21),2]='Astrocyte'
celltype[celltype$ClusterID %in% c(33),2]='NA'
celltype[celltype$ClusterID %in% c(26,27,31,18,13,10,5,0),2]='Oligo'
celltype[celltype$ClusterID %in% c(28,7,35,30,34),2]='Opc'
celltype[celltype$ClusterID %in% c(12,32,29),2]='Macro/Micro'
celltype[celltype$ClusterID %in% c(24),2]='Endothelial cells'
celltype[celltype$ClusterID %in% c(17),2]='Fibroblasts'

celltype
table(celltype$celltype)

sce.in=scRNA_harmony
#先加一列celltype所有值为空
sce.in@meta.data$celltype = "NA"
###注释
for(i in 1:nrow(celltype)){
  sce.in@meta.data[which(sce.in@meta.data$seurat_clusters == celltype$ClusterID[i]),'celltype'] <- celltype$celltype[i]}
table(sce.in@meta.data$celltype)


p <- DotPlot(sce.in, features = unique(markers),
             assay='RNA' ,group.by = 'celltype' )  + coord_flip()

p
# ggsave(plot=p, filename="check_marker_by_celltype.pdf")
ggsave(filename="check_marker_by_celltype.pdf",width = 210,height = 297,units = "mm",plot = p)



table(sce.in@meta.data$celltype,sce.in@meta.data$seurat_clusters)
sce=sce.in
p.dim.cell=DimPlot(sce, reduction = "tsne", group.by = "celltype",label = T,pt.size = 1)
p.dim.cell
ggsave(plot=p.dim.cell,filename="7DimPlot_tsne_celltype.pdf",width=9, height=7)
p.dim.cell=DimPlot(sce, reduction = "umap", group.by = "celltype",label = T,pt.size = 1,repel = T)
p.dim.cell
library(ggplot2)
ggsave(plot=p.dim.cell,filename="7DimPlot_umap_celltype.pdf",width=9, height=7)

DefaultAssay(sce) <- "RNA"
# Defaultassay设置为RNA",意味着接下来的分析将基于原始值
save(sce,file='7cluster细胞注释-改Rdata_分16cluster.Rdata')
# load("2.基因注释/12cluster细胞注释-改Rdata_分16cluster")

Idents(scRNA)=scRNA$celltype
scRNA=sce
save(scRNA,file ='scRNA_anno.Rdata')
###注释
markers = c('PECAM1','VWF',#endothelial
            'DCN','COL1A2',#Fibroblasts
            "TYROBP","CTSS","C1QA","HEXB","AIF1",# Micro/Macro
            'C1QA','C1QB','APOE',#Macrophage
            'ALDH1L1','AQP4','GFAP','SLC1A2',#Astrocyte
            'MBP','MOBP','MOG','PLP1',#Oligo
            'RBFOX3','MAP2','GAD1',#Neurons
            'GAD1','GAD2','PVALB','CALB2','VIP',#Inhibitory Neurons
            'SLC17A7','CAMK2A','SATB2',#Excitatory Neurons
            'CSPG4','PDGFRA','GPR17'#opc
)

# Generate the DotPlot with the new order

levels_in_celltype <- levels(as.factor(scRNA@meta.data$celltype))
print(levels_in_celltype)

# 检查是否有重复水平
duplicated_levels <- levels_in_celltype[duplicated(levels_in_celltype)]
print(duplicated_levels)
scRNA@meta.data$celltype <- factor(scRNA@meta.data$celltype, levels = unique(scRNA@meta.data$celltype))


new_order <- c("Endothelial cells", "Fibroblasts", "Macro/Micro", "Astrocyte", "Neurons","Opc","Oligo","NA")
scRNA@meta.data$celltype <- factor(scRNA@meta.data$celltype, levels = new_order)

dotplot <- DotPlot(scRNA, features = unique(markers), group.by = "celltype") +
  coord_flip() +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(hjust = 1, vjust = 0.5, angle = 90)) +
  labs(x = NULL, y = NULL) +
  guides(size = guide_legend(order = 3)) +
  scale_color_gradientn(values = seq(0, 1, 0.2), colours = c('#008D83', '#B0A8B9', '#FF8066', '#C34A36'))

print(dotplot)

# Save the plot
ggsave(filename = "7-cluster-marker24-12-31-3.pdf", plot = dotplot, width = 110, height = 160, units = "mm")


load("scRNA_anno.Rdata")
###########差异分析热图#########
library(Seurat)
library(dplyr)
library(reticulate)
library(cowplot)
library(ggplot2)
library(viridis)
library(tidyr)
library(magrittr)
library(reshape2)
library(readr)
library(stringr)
# 输入数据预处理
scRNA=scRNA
# 查看细胞数量
table(scRNA@meta.data$celltype)
unique(scRNA@meta.data$celltype)
# 把annotations重新排个顺序，这里按照细胞数量排序
scRNA@meta.data$celltype <- factor(scRNA@meta.data$celltype,
                                   levels = c("Endothelial cells", "Fibroblasts", "Macro/Micro", "Astrocyte", "Neurons","Opc","Oligo"))

#提取设置好细胞类型的颜色
colour=c("#f8766d","#c49a00","#53b400", "#00c094","#00b6eb","#a58aff","#fb61d7")


Idents(scRNA) <- "celltype"
####FindAllMarkers差异分析/热图####
markers <- FindAllMarkers(object = scRNA, test.use="wilcox" ,
                          only.pos = TRUE,
                          logfc.threshold = 0.5,
                          min.pct = 0.5)
all.markers =markers %>% dplyr::select(gene, everything()) %>% subset(p_val<0.05)

#保存/读取
write.csv(all.markers,"分群之后all.markers.csv",row.names = T)
all.markers <- read.csv("分群之后all.markers.csv",header = T,row.names = 1)

#取top marker展示
top10 = all.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC) #将每个cluster lgFC排在前10的marker基因挑选出来
write.csv(top10,"分群之后top10.markers.csv",row.names = T)

top1 = all.markers %>% group_by(cluster) %>% top_n(n = 1, wt = avg_log2FC) #将每个cluster lgFC排在前10的marker基因挑选出来
markers <- top1$gene


####DoHeatmap热图####
top10 <- as.data.frame(top10)  # 将markers对象转换为数据框
markerdata <- ScaleData(scRNA, features = as.character(unique(top10$gene)), assay = "RNA")  # 对特定基因进行标准化
# 常规热图
DoHeatmap=DoHeatmap(markerdata,  # 使用DoHeatmap函数绘制热图
                    features = as.character(unique(top10$gene)),  # 指定要显示的特征
                    group.by = "celltype",  # 指定按照免疫细胞类型分组
                    assay = 'RNA')  # 指定分析的分析类型为RNA
ggsave(plot=DoHeatmap,filename="7DoHeatmap热图_celltype.pdf",width=18, height=15)

#修改颜色
DoHeatmap=DoHeatmap(markerdata,  # 再次绘制热图，附加自定义颜色和渐变色
                    features = as.character(unique(top10$gene)),
                    group.by = "celltype",
                    assay = 'RNA',
                    group.colors = c("#f8766d","#c49a00","#53b400", "#00c094","#00b6eb","#a58aff","#fb61d7")) +
  scale_fill_gradientn(colors = c("white","grey","firebrick3"))
ggsave(plot=DoHeatmap,filename="7DoHeatmap热图_celltype1.pdf",width=18, height=15)

####Umap图####
# 标准umap图
colour=c("#1E90FF","#7CFC00", "#808000","#FF00FF","#FA8072","#7B68EE","#800080","#A0522D","#D2B48C","#D2691E","#87CEEB","#40E0D0","#5F9EA0",
         
         "#FF1493","#0000CD","#008B8B")

p1 <- DimPlot(scRNA, reduction = "umap", group.by = "celltype",cols = colour,label = TRUE,repel = TRUE)
# group.by 指定分组umap图
p2 <- DimPlot(scRNA, reduction = "umap", group.by = "group")

p1 + p2

pdf("scRNA.UMAP.pdf", height = 5, width = 12)
p1 + p2
dev.off()

pdf("scRNA.UMAP.pdf", height = 5,width = 12)

#split.by = "stim"拆分分组umap图
p3 <- DimPlot(scRNA, reduction = "umap", group.by = "celltype",
              split.by = "group",cols = colour,repel = TRUE)
pdf("scRNA.UMAP1.pdf", height = 5, width = 12)
p3
dev.off()

####细胞比例图orig.ident####


table(scRNA$group)
prop.table(table(Idents(scRNA)))  # 计算细胞类型的比例
table(Idents(scRNA), scRNA$group)  # 计算分组细胞类型的数量
Cellratio <- prop.table(table(Idents(scRNA), sce.in$group), margin = 2)  # 计算细胞类型相对于原始身份的比例，并将结果转换为数据框
Cellratio <- as.data.frame(Cellratio)
colnames(Cellratio) <- c("celltype","sample","ratio")



colourCount = length(unique(Cellratio$celltype))  # 计算唯一的细胞类型数量
ggplot(Cellratio) +  # 创建 ggplot 对象并传入数据框 Cellratio
  geom_bar(aes(x = sample, y = ratio, fill = celltype), stat = "identity", width = 0.7, size = 0.5, colour = '#222222') +  # 添加条形图层，根据细胞类型填充颜色
  theme_classic() +  # 应用经典主题
  labs(x = 'Sample', y = 'Ratio') +  # 设置 x 和 y 轴标签
  coord_flip() +  # 翻转坐标轴
  theme(panel.border = element_rect(fill = NA, color = "black", size = 0.5, linetype = "solid"))+   # 设置面板边框属性
  scale_fill_manual(values = colour)  # 手动设置填充颜色

#直立柱状图
p3=ggplot(Cellratio) +
  geom_bar(aes(x =sample, y= ratio, fill = celltype),stat = "identity",width = 0.7,size = 0.5,colour = '#222222')+
  theme_classic() +
  labs(x='Sample',y = 'Ratio')+
  scale_fill_manual(values = colour)+
  theme(panel.border = element_rect(fill=NA,color="black", size=0.5, linetype="solid"))
p3
pdf("scRNA.细胞比例.pdf", height = 6, width = 5)
p3
dev.off()

#######危险因素评分##########
####all细胞AA代谢比例######
####先区分AD与nc####
library(dplyr)
library(hdf5r)
library(Seurat)
library(data.table, lib.loc = "C:/Users/Administrator/AppData/Local/R/win-library/4.4")
library(Seurat)
library(SeuratDisk)
library(ggplot2)
library(Seurat)
library(dplyr)
library(ggplot2)
library(tidyr)
library(patchwork)
library(clustree)
library(clusterProfiler)
scRNA_harmony=scRNA
####先定比较的类型celltype####
table(scRNA@meta.data$celltype)
head(scRNA@meta.data)

## ## ## ## ## ## ## ## ## ## ## ## ## ## ## up
SLC_mediated=read.gmt('./危险因素/SLC_mediated_transmembrane_transport.gmt')
DotPlot(scRNA,features = SLC_mediated$gene)

SLC_mediated=list(SLC_mediated$gene)
names(SLC_mediated)='SLC_mediated'

scRNA=AddModuleScore(scRNA,features = SLC_mediated,name = 'SLC_mediated')

####
Transport_of=read.gmt('./危险因素/Transport_of_small_molecules.gmt')
DotPlot(scRNA,features = Transport_of$gene)

Transport_of=list(Transport_of$gene)
names(Transport_of)='Transport_of'

scRNA=AddModuleScore(scRNA,features = Transport_of,name = 'Transport_of')

####


Disorders_of=read.gmt('./危险因素/Disorders_of_transmembrane_transporters.gmt')
DotPlot(scRNA,features = Disorders_of$gene)

Disorders_of=list(Disorders_of$gene)
names(Disorders_of)='Disorders_of'

scRNA=AddModuleScore(scRNA,features = Disorders_of,name = 'Disorders_of')
####


Astrocytic_Glutamate=read.gmt('./危险因素/Astrocytic_Glutamate_Glutamine_Uptake_And_Metabolism.gmt')
DotPlot(scRNA,features = Astrocytic_Glutamate$gene)

Astrocytic_Glutamate=list(Astrocytic_Glutamate$gene)
names(Astrocytic_Glutamate)='Astrocytic_Glutamate'

scRNA=AddModuleScore(scRNA,features = Astrocytic_Glutamate,name = 'Astrocytic_Glutamate')


####

Neurotransmitter_uptake=read.gmt('./危险因素/Neurotransmitter_uptake_and_metabolism_In_glial_cells.gmt')
DotPlot(scRNA,features = Neurotransmitter_uptake$gene)

Neurotransmitter_uptake=list(Neurotransmitter_uptake$gene)
names(Neurotransmitter_uptake)='Neurotransmitter_uptake'

scRNA=AddModuleScore(scRNA,features = Neurotransmitter_uptake,name = 'Neurotransmitter_uptake')


#####
SLC_transporter=read.gmt('./危险因素/SLC_transporter_disorders.gmt')
DotPlot(scRNA,features = SLC_transporter$gene)

SLC_transporter=list(SLC_transporter$gene)
names(SLC_transporter)='SLC_transporter'

scRNA=AddModuleScore(scRNA,features = SLC_transporter,name = 'SLC_transporter')


#####
Glycogen_breakdown=read.gmt('./危险因素/Glycogen_breakdown.gmt')
DotPlot(scRNA,features = Glycogen_breakdown$gene)

Glycogen_breakdown=list(Glycogen_breakdown$gene)
names(Glycogen_breakdown)='Glycogen_breakdown'

scRNA=AddModuleScore(scRNA,features = Glycogen_breakdown,name = 'Glycogen_breakdown')


#####
Glutathione_synthesis=read.gmt('./危险因素/Glutathione_synthesis_and_recycling.gmt')
DotPlot(scRNA,features = Glutathione_synthesis$gene)

Glutathione_synthesis=list(Glutathione_synthesis$gene)
names(Glutathione_synthesis)='Glutathione_synthesis'

scRNA=AddModuleScore(scRNA,features = Glutathione_synthesis,name = 'Glutathione_synthesis')


######
mRNA_protein=read.gmt('./危险因素/mRNA_protein_and_metabolite_inducation_pathway_by_cyclosporin.gmt')
DotPlot(scRNA,features = mRNA_protein$gene)

mRNA_protein=list(mRNA_protein$gene)
names(mRNA_protein)='mRNA_protein'

scRNA=AddModuleScore(scRNA,features = mRNA_protein,name = 'mRNA_protein')


##############
Carnitine_metabolism=read.gmt('./危险因素/GOBP_CARNITINE_METABOLIC_PROCESS.v2025.1.Hs.gmt')
DotPlot(scRNA,features = Carnitine_metabolism$gene)

Carnitine_metabolism=list(Carnitine_metabolism$gene)
names(Carnitine_metabolism)='Carnitine_metabolism'

scRNA=AddModuleScore(scRNA,features = Carnitine_metabolism,name = 'Carnitine_metabolism')


p_up_markers <- DotPlot(scRNA, features = c('SLC_mediated1','Transport_of1','Disorders_of1','Astrocytic_Glutamate1','Neurotransmitter_uptake1','Carnitine_metabolism1','mRNA_protein1','Glutathione_synthesis1','Glycogen_breakdown1','SLC_transporter1'),
                        assay='RNA' ,group.by = 'celltype' )

p_up_markers


####提取数据####
data<-p_up_markers$data

colnames(data)

colnames(data)<-c("AverageExpression_unscaled","Precent Expressed","Features","celltype","Average Expression")

unique(data$`Precent Expressed`)
library(ggplot2)
####用ggplot画图####
p = ggplot(data,aes(celltype,Features,size = `Precent Expressed` ))+
  geom_point(shape=21,aes(fill= `Average Expression`),position =position_dodge(0))+
  theme_minimal()+xlab(NULL)+ylab(NULL) +
  scale_size_continuous(range=c(1,10))+theme_bw()+
  scale_fill_gradient(low = "grey", high = "#C34A36")+
  theme(legend.position = "right",legend.box = "vertical", #图例位置
        legend.margin=margin(t= 0, unit='cm'),
        legend.spacing = unit(0,"in"),
        axis.text.x  = element_text(color="black",size=16,angle = 45,
                                    vjust = 0.5, hjust=0.5),#x轴
        axis.text.y  = element_text(color="black",size=12),#y轴
        legend.text = element_text(size =12,color="black"),#图例
        legend.title = element_text(size =12,color="black"),#图例
        axis.title.y=element_text(vjust=1,
                                  size=16)
  )+labs(x=" ",y = "Features")+coord_flip();p

ggsave("危险因素活性分析_plot1.pdf", plot = p, device = "pdf", width = 7, height = 6)

rm(list = ls())
load("7cluster细胞注释-改Rdata_分16cluster.Rdata")
scRNA=sce
rm(sce)
####保护因素活性分析down######
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## down
#######
Metabolism_of=read.gmt('./保护因素/Metabolism_of_amino_acids_and_derivativesv2025.1.Hs.gmt')
DotPlot(scRNA,features = Metabolism_of$gene)

Metabolism_of=list(Metabolism_of$gene)
names(Metabolism_of)='Metabolism_of'

scRNA=AddModuleScore(scRNA,features = Metabolism_of,name = 'Metabolism_of')


#####
SLC_mediated=read.gmt('./保护因素/SLC_mediated_transmembrane_transport.gmt')
DotPlot(scRNA,features = SLC_mediated$gene)

SLC_mediated=list(SLC_mediated$gene)
names(SLC_mediated)='SLC_mediated'

scRNA=AddModuleScore(scRNA,features = SLC_mediated,name = 'SLC_mediated')


####
amino_acid=read.gmt('./保护因素/amino_acid_metabolism.gmt')
DotPlot(scRNA,features = amino_acid$gene)

amino_acid=list(amino_acid$gene)
names(amino_acid)='amino_acid'

scRNA=AddModuleScore(scRNA,features = amino_acid,name = 'amino_acid')


####
Transport_of_small=read.gmt('./保护因素/Transport_of_small_molecules.gmt')
DotPlot(scRNA,features = Transport_of_small$gene)

Transport_of_small=list(Transport_of_small$gene)
names(Transport_of_small)='Transport_of_small'

scRNA=AddModuleScore(scRNA,features = Transport_of_small,name = 'Transport_of_small')


####
Transport_of_inorganic=read.gmt('./保护因素/Transport_of_inorganic_cations.gmt')
DotPlot(scRNA,features = Transport_of_inorganic$gene)

Transport_of_inorganic=list(Transport_of_inorganic$gene)
names(Transport_of_inorganic)='Transport_of_inorganic'

scRNA=AddModuleScore(scRNA,features = Transport_of_inorganic,name = 'Transport_of_inorganic')

####
Arginine_and_proline=read.gmt('./保护因素/Arginine_and_proline_metabolism.gmt')
DotPlot(scRNA,features = Arginine_and_proline$gene)

Arginine_and_proline=list(Arginine_and_proline$gene)
names(Arginine_and_proline)='Arginine_and_proline'

scRNA=AddModuleScore(scRNA,features = Arginine_and_proline,name = 'Arginine_and_proline')


####
Amino_acid_transport=read.gmt('./保护因素/Amino_acid_transport_across_the_plasma_membrane.gmt')
DotPlot(scRNA,features = Amino_acid_transport$gene)

Amino_acid_transport=list(Amino_acid_transport$gene)
names(Amino_acid_transport)='Amino_acid_transport'

scRNA=AddModuleScore(scRNA,features = Amino_acid_transport,name = 'Amino_acid_transport')


####
Trans_sulfuration=read.gmt('./保护因素/Trans_sulfuration_one_carbon_metabolism.gmt')
DotPlot(scRNA,features = Trans_sulfuration$gene)

Trans_sulfuration=list(Trans_sulfuration$gene)
names(Trans_sulfuration)='Trans_sulfuration'

scRNA=AddModuleScore(scRNA,features = Trans_sulfuration,name = 'Trans_sulfuration')

####
Argininemia=read.gmt('./保护因素/Argininemia.gmt')
DotPlot(scRNA,features = Argininemia$gene)

Argininemia=list(Argininemia$gene)
names(Argininemia)='Argininemia'

scRNA=AddModuleScore(scRNA,features = Argininemia,name = 'Argininemia')

p_down_markers <- DotPlot(scRNA, features = c('Metabolism_of1','SLC_mediated1','amino_acid1','Transport_of_small1','Transport_of_inorganic1','Arginine_and_proline1','Amino_acid_transport1','Trans_sulfuration1','Argininemia1'),
                          assay='RNA' ,group.by = 'celltype' )

p_down_markers


####提取数据####
data<-p_down_markers$data

colnames(data)

colnames(data)<-c("AverageExpression_unscaled","Precent Expressed","Features","celltype","Average Expression")

unique(data$`Precent Expressed`)
library(ggplot2)
####用ggplot画图####
p = ggplot(data,aes(celltype,Features,size = `Precent Expressed` ))+
  geom_point(shape=21,aes(fill= `Average Expression`),position =position_dodge(0))+
  theme_minimal()+xlab(NULL)+ylab(NULL) +
  scale_size_continuous(range=c(1,10))+theme_bw()+
  scale_fill_gradient(low = "grey", high = "#008B74")+
  theme(legend.position = "right",legend.box = "vertical", #图例位置
        legend.margin=margin(t= 0, unit='cm'),
        legend.spacing = unit(0,"in"),
        axis.text.x  = element_text(color="black",size=16,angle = 45,
                                    vjust = 0.5, hjust=0.5),#x轴
        axis.text.y  = element_text(color="black",size=12),#y轴
        legend.text = element_text(size =12,color="black"),#图例
        legend.title = element_text(size =12,color="black"),#图例
        axis.title.y=element_text(vjust=1,
                                  size=16)
  )+labs(x=" ",y = "Features")+coord_flip();p

ggsave("保护因素活性分析_plot1.pdf", plot = p, device = "pdf", width = 7, height = 6)




## 评估受体通路活性
receptor=read.gmt('./REACTOME_METABOLISM_OF_AMINO_ACIDS_AND_DERIVATIVES.v2024.1.Hs.gmt')
DotPlot(scRNA,features = receptor$gene)

receptor=list(receptor$gene)
names(receptor)='AAmetabolism'

scRNA_Ne_AA=AddModuleScore(scRNA_Ne,features = receptor,name = 'AAmetabolism')

median(scRNA_Ne_AA$AAmetabolism1)

scRNA_Ne_AA$group1=ifelse(scRNA_Ne_AA$AAmetabolism1> median(scRNA_Ne_AA$AAmetabolism1),'AAmetabolism_high','AAmetabolism_low')
Idents(scRNA_Ne_AA)=scRNA_Ne_AA$group1


###############细胞高低代谢分布#######################
scRNAsub=scRNA_Ma_i
median(scRNAsub$AAmetabolism1)
scRNAsub$group1=ifelse(scRNAsub$AAmetabolism1> median(scRNAsub$AAmetabolism1),'AAmetabolism_high','AAmetabolism_low')
Idents(scRNAsub)=scRNAsub$group1



library(ggplot2)
library(Seurat)

scRNAsub$group_AAmetabolism <- ifelse(scRNAsub$AAmetabolism1 > 0.7896365, "high", "low")

# 自定义颜色TRUE
custom_colors <- c("high" = "blue", "low" = "yellow")

# 按照 Nucleic_immunity1 的高低分组，基于 UMAP 图显示
DimPlot(scRNAsub, reduction = "umap", group.by = "group_Nucleic",
        cols = custom_colors, label = F) +
  ggtitle("UMAP - group_AAmetabolism Group Distribution") +
  theme_minimal()



#########GSEA富集分析###########

library(Seurat)

## findmarker寻找差异基因
df=FindMarkers(scRNAsub,ident.1 = 'NAMmetabolism_high',ident.2 = 'NAMmetabolism_low',logfc.threshold = 0)
df1 = df %>%
  filter( pct.1 > 0.1 & p_val_adj < 0.05 ) %>% # &：和；|：或。
  filter( abs( avg_log2FC ) > 0.25 ) #%>% # 按 logFC 绝对值筛选
# filter( avg_log2FC > 0 ) # 按 logFC 正负筛选



gc()
df1$mRNAs=rownames(df1)
setwd("D:\\AD-大脑单细胞\\GSE157827")
load('AAmetabolism_Ne_deg.Rdata')
load('scRNA_T_Nucleic_immunity.Rdata')
## 保存
save(df,file ='deg.Rdata')
load('deg.Rdata')
#4.制作genelist
gene <- df1$mRNAs
## 转换
library(clusterProfiler)
gene = bitr(gene, fromType="SYMBOL",
            toType="ENTREZID",
            OrgDb="org.Hs.eg.db")
## 去重
gene <- dplyr::distinct(gene,
                        SYMBOL,.keep_all=TRUE)

gene_df <- data.frame(logFC=df$avg_log2FC,
                      SYMBOL = df$mRNAs)
gene_df <- merge(gene_df,
                 gene,by="SYMBOL")

## geneList 三部曲
## 1.获取基因logFC
geneList <- gene_df$logFC
## 2.命名
names(geneList) = gene_df$ENTREZID
## 3.排序很重要
geneList = sort(geneList,
                decreasing = TRUE)

#5.运行GSEA分析
library(clusterProfiler)
## 读入hallmarks gene set，从哪来？
hallmarks <- read.gmt("c5.go.bp.v2023.1.Hs.entrez.gmt")
# 需要网络

y <- GSEA(geneList,
          TERM2GENE =hallmarks,pvalueCutoff = 0.25,seed = 123)


y=setReadable(y,OrgDb="org.Hs.eg.db",keyType = 'ENTREZID')


rt=y@result
write.csv(rt, file = "AAmetabolism low-vs-high_Ne_GO.csv")
## 显著的通路
rt=rt[rt$p.adjust<0.25,]

## NES>0代表正相关
rt1=rt[rt$NES>0,]
rt1=rt1[order(rt1$NES,decreasing = T),]
rt1=rt1[1:10,]
rt2=rt[rt$NES<0,]
rt2=rt2[order(rt2$NES),]
rt2=rt2[1:10,]

# 重新组合用来画图
rt=rbind(rt1,rt2)
library(ggplot2)
library(tidyverse)
rt$group <- ''
rt$group[which(rt$NES>0)]='up'
rt$group[which(rt$NES <0)]='down'
colnames(rt)
ggplot(rt,aes(reorder(Description,NES),NES,fill=group))+
  geom_col()+
  theme_bw()+
  theme(panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        panel.border = element_blank(),
        legend.title = element_blank(),
        axis.text = element_text(color="black",size=10),
        axis.line.x = element_line(color='black'),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),
        legend.position = 'none')+
  coord_flip()+
  geom_segment(aes(y=0, yend=0,x=0,xend=18.5))+
  geom_text(data = rt[which(rt$NES>0),],aes(x=Description, y=-0.01, label=Description),
            hjust=1, size=4)+
  geom_text(data = rt[which(rt$NES<0),],aes(x=Description, y=0.01, label=Description),
            hjust=0, size=4)+
  geom_text(data = rt[which(rt$NES>0),],aes(label=p.adjust),
            hjust=-0.1, size=4, color='#845EC2')+
  geom_text(data = rt[which(rt$NES<0),],aes(label=p.adjust),
            hjust=1.1, size=4, color="#845EC2")+
  scale_fill_manual(values = c("#4D8076",
                               "#C34A36"))+
  scale_x_discrete(expand = expansion(mult = c(0,0)))+
  labs(x='', y='NES')



setwd("D:/AD_GWAS/AD-血液单细胞/GSE226602")

dir.create("~/SeuratV4")
# 然后安装的时候，指定安装目录
install.packages('Seurat', repos = c('https://satijalab.r-universe.dev'), lib = "~/SeuratV4")
.libPaths(c("~/SeuratV4", .libPaths()))

library(limma)
library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)
#BiocManager::install("SingleR")
library(SingleR)
library(CCA)
library(clustree)
library(cowplot)
#if (!require("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install("HSMMSingleCell")
library(monocle)
library(tidyverse)
library(Seurat)

gc()
expression=readRDS('./GSE226602_rna_lognorm_expression1.rds')
counts=readRDS('./GSE226602_rna_raw_counts1.rds')

# 创建Seurat对象
sc <- CreateSeuratObject(counts = counts)

# 计算线粒体比例
sc[["mt_percent"]] <- PercentageFeatureSet(sc, pattern = "^MT-")

# 计算红细胞比例
HB_genes <- c("HBA1","HBA2","HBB","HBD","HBE1","HBG1","HBG2","HBM","HBQ1","HBZ")  # 定义红细胞基因列表
HB_m <- match(HB_genes, rownames(sc))  # 在Seurat对象的RNA数据中查找红细胞基因的索引位置
HB_genes <- HB_genes[!is.na(HB_m)]  # 删除NA值（未匹配到的基因）
sc[["HB_percent"]] <- PercentageFeatureSet(sc, features = HB_genes)  # 计算红细胞基因比例
# 查看 Seurat 对象中的总细胞数
total_cells <- ncol(sc)

# 打印总细胞数
cat("Total number of cells:", total_cells, "\n")
library(ggplot2)

# 绘制质量控制小提琴图
violin_before <- VlnPlot(sc, features = c("nFeature_RNA", "nCount_RNA", "mt_percent", "HB_percent"), pt.size = 0.01, ncol = 4)
ggsave("violin_before1.pdf", plot = violin_before, width = 25, height = 5)


# 过滤细胞
sc <- subset(sc, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & mt_percent < 10 & HB_percent < 3 & nCount_RNA < quantile(nCount_RNA, 0.97) & nCount_RNA > 1000)

# 您可以在这里打印过滤后的结果以确认
table(sc$orig.ident)
total_cells <- ncol(sc)

# 打印总细胞数
cat("Total number of cells:", total_cells, "\n")

# 绘图
violin_after <- VlnPlot(sc,
                        features = c("nFeature_RNA", "nCount_RNA", "mt_percent","HB_percent"),
                        pt.size = 0.01,
                        ncol = 4)
# 将图片输出到画板上
violin_after
# 保存图片
ggsave("vlnplot_after_qc.pdf", plot = violin_after, width = 50, height =15, limitsize = FALSE)




library(magrittr) # 只加载管道操作符
# 或者
library(tidyverse) # 加载一系列数据处理包，包括 %>%

####数据归一化、筛选高变基因与PCA降维####
# harmony整合是基于PCA降维结果进行的。
sc <- NormalizeData(sc) %>% # 数据归一化处理
  FindVariableFeatures(selection.method = "vst",nfeatures = 3000) %>% #筛选高变基因
  ScaleData() %>% #数据标准化
  RunPCA(npcs = 30, verbose = T)#npcs：计算和存储的PC数（默认为 50）
a=DimPlot(sc,reduction = "pca",group.by = "orig.ident")
#PCA图看到还有一点的批次效应（融合得比较好批次就弱）
a

##看下高变基因有哪些可视化
# 提取前15个高变基因ID
top15 <- head(VariableFeatures(sc), 15)
plot1 <- VariableFeaturePlot(sc)
plot2 <- LabelPoints(plot = plot1, points = top15, repel = TRUE, size=3)
# 合并图片
feat_15 <- CombinePlots(plots = list(plot1,plot2),legend = "bottom")
feat_15
# 保存图片
ggsave(file = "feat_15.pdf",plot = feat_15,he = 10,wi = 15 )


####细胞周期评分####
# 提取g2m特征向量
g2m_genes = cc.genes$g2m.genes
g2m_genes = CaseMatch(search = g2m_genes, match = rownames(sc))
# 提取s期特征向量
s_genes = cc.genes$s.genes
s_genes = CaseMatch(search = s_genes, match = rownames(sc))
# 对细胞周期阶段进行评分
sc <- CellCycleScoring(object=sc,  g2m.features=g2m_genes,  s.features=s_genes)
sc=CellCycleScoring(object = sc,
                    s.features = s_genes,
                    g2m.features = g2m_genes,
                    set.ident = TRUE)#set.ident 是否给每个细胞标注一个细胞周期标签
sc <- CellCycleScoring(object=sc,  g2m.features=g2m_genes,  s.features=s_genes)

p4=VlnPlot(sc, features = c("S.Score", "G2M.Score"), group.by = "orig.ident",
           ncol = 2, pt.size = 0.1)
p4

P5=sc@meta.data%>% ggplot(aes(S.Score,G2M.Score))+geom_point(aes(color=Phase))+
  theme_minimal()
ggsave(file = "s.G2M,score.pdf",plot = p4,he = 5,wi = 20 )
ggsave(file = "s.G2M,score2.pdf",plot = P5,he = 5,wi = 6 )



library(harmony)
####RunHarmony去批次####
# 整合需要指定Seurat对象和metadata中需要整合的变量名。
scRNA_harmony <- RunHarmony(sc, group.by.vars = "orig.ident")
scRNA_harmony@reductions[["harmony"]][[1:5,1:5]]
b=DimPlot(scRNA_harmony,reduction = "harmony",group.by = "orig.ident")
#PCA图看到还有一点的批次效应（融合得比较好批次就弱）
b
# 合并图片
pca_harmony_integrated <- CombinePlots(list(a,b),ncol=1)
pca_harmony_integrated

#去批次前后对比
ggsave(file = "RunHarmony去批次.pdf",plot = pca_harmony_integrated,he = 15,wi = 12 )
table(scRNA_harmony@meta.data$orig.ident)

#去批次前后对比
# View the existing metadata
# 1) 定义样本列表（注意：orig.ident 中是类似 "G1028" 的格式）

# 您给定的数字清单
control_nums <- c("1028","1273","911","1180","1279","598","780","1020","968","781",
                  "1111","912","978","1282","820","1162","989","1010","836","1200",
                  "970","905")

ad_nums <- c("1241","906","1034","1055","1092","1120","1160","254","1052","516",
             "696","738","773","863","917","932","802","942","1147","656","953",
             "965","921","230","1236","1237","947","70")

# 构造“带G”和“裸数字”两套ID，并合并去重
control_ids <- unique(c(control_nums, paste0("G", control_nums)))
ad_ids      <- unique(c(ad_nums,      paste0("G", ad_nums)))

# 取出当前对象中的样本ID
orig_ids <- scRNA_harmony@meta.data$orig.ident

# 分组（同时兼容两种格式）
group <- ifelse(orig_ids %in% control_ids, "Control",
                ifelse(orig_ids %in% ad_ids, "AD", NA))

scRNA_harmony@meta.data$Group <- factor(group, levels = c("Control", "AD"))

# 核对：每组细胞数
cat("\nCell counts by Group:\n")
print(table(scRNA_harmony@meta.data$Group, useNA = "ifany"))

# 核对：每个样本在各组内的细胞数
cat("\nPer-sample cell counts within each Group:\n")
per_sample <- with(scRNA_harmony@meta.data, table(Group, orig.ident, useNA = "no"))
print(per_sample)



# 保存结果（如需保存映射后的 Seurat 对象）
saveRDS(scRNA_harmony, file = "scRNA_harmony.rds")

scRNA_harmony=readRDS('./scRNA_harmony.rds')##这个是T细胞类型的T_celltype,tissue_type,celltype
setwd("D:/AD_GWAS/AD脑组织/GSE157827/rest")

colnames(scRNA_harmony@meta.data)[1:20]
# 质控图按新分组绘制
pdf("分组质控.pdf")
VlnPlot(scRNA_harmony, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.rb"),
        group.by = "Group", pt.size = 0.1, ncol = 2) +
  theme(legend.position = "right")
dev.off()


# 后续都是基于Harmony矫正之后的数据，不是基因表达数据和直接的PCA降维数据。
# 设置reduction = 'harmony'，后续分析是基于Harmony矫正之后的数据。
####聚类、umap/tsne降维降维####
ElbowPlot=ElbowPlot(scRNA_harmony, ndims=50, reduction="pca")
ElbowPlot(scRNA_harmony, ndims=50, reduction="harmony")

ElbowPlot
ggsave( ElbowPlot, filename = "ElbowPlot.pdf", width = 7, height = 5 )

#选择dim30，reso0.3
scRNA_harmony <- FindNeighbors(scRNA_harmony, reduction = "harmony", dims = 1:20) %>% FindClusters(resolution = 0.7)
table(scRNA_harmony@meta.data$seurat_clusters)
table(scRNA_harmony@meta.data$seurat_clusters)
total_cells <- sum(table(scRNA_harmony@meta.data$seurat_clusters))
table(total_cells)



Reductions(scRNA_harmony)
# or
names(scRNA_harmony@reductions)
scRNA_harmony <- FindNeighbors(scRNA_harmony, reduction = "pca", dims = 1:20)
scRNA_harmony <- RunTSNE(scRNA_harmony, reduction = "pca", dims = 1:20)
DimPlot(scRNA_harmony, reduction = "tsne", group.by = "Group")
##umap/tsne降维
scRNA_harmony <- RunTSNE(scRNA_harmony, reduction = "harmony", dims = 1:20)
scRNA_harmony <- RunUMAP(scRNA_harmony, reduction = "harmony", dims = 1:20)
#scRNA_harmony=sce
# 绘图
umap_integrated1 <- DimPlot(scRNA_harmony, reduction = "umap", group.by = "Group")
umap_integrated2 <- DimPlot(scRNA_harmony, reduction = "umap", label = TRUE)
tsne_integrated1 <- DimPlot(scRNA_harmony, reduction = "tsne", group.by = "Group")
tsne_integrated2 <- DimPlot(scRNA_harmony, reduction = "tsne", label = TRUE)
# 合并图片
umap_tsne_integrated <- CombinePlots(list(umap_integrated1,umap_integrated2),ncol=2)
# 将图片输出到画板
umap_tsne_integrated
# 保存图片
ggsave("22-umap_tsne_integrated1.pdf",umap_tsne_integrated,wi=20,he=7.5)
#保存数据
save(scRNA_harmony,file = "22-cluster-scdata.Rdata")


load("22-cluster-scdata.Rdata")

table(scRNA_harmony@meta.data$seurat_clusters)
###### 计算各组中每种细胞类型的数量######
table(Idents(scRNA_harmony), scRNA_harmony$Group)

####细胞注释####


# 差异分析
markers <- FindAllMarkers(object = scRNA_harmony, test.use="wilcox" ,
                          only.pos = TRUE,
                          logfc.threshold = 0.25)


#2.对计算好的每cluster的marker基因进行筛选
all.markers =markers %>% dplyr::select(gene, everything()) %>% subset(p_val_adj<0.05)
#筛选出P<0.05的marker基因
top10 = all.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC) #将每个cluster lgFC排在前10的marker基因挑选出来
top15 = all.markers %>% group_by(cluster) %>% top_n(n = 15, wt = avg_log2FC)
top20 = all.markers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_log2FC)
top30 = all.markers %>% group_by(cluster) %>% top_n(n = 30, wt = avg_log2FC)


View(top10)
write.csv(top10,"12cluster_top10.csv",row.names = T)
write.csv(top15,"12cluster_top15.csv",row.names = T)
write.csv(top20,"12cluster_top20.csv",row.names = T)
write.csv(top30,"12cluster_top30.csv",row.names = T)
write.csv(all.markers,"12cluster_all.markers.csv",row.names = T)



markers <- c(
  "CD3D","CD3E","CCR7","IL7R",                # Naive CD4+ T Cell
  "CD8A","CD8B","GZMK","GZMB",                # CD8+ T Cell
  "FCGR3A","LYZ","VCAN",                      # CD16+ Monocyte
  "MS4A1","CD79A","CD19",                     # B Cell
  "MZB1","IGKC","JCHAIN",                     # Plasma Cell
  "GNLY","NKG7","PRF1",                       # Natural Killer (NK) Cell
  "IL3RA","GZMB","CLEC4C",                    # Dendritic Cell (pDC)
  "S100A8","S100A9","CSF3R",                  # Neutrophil
  "PPBP","PF4","ITGA2B"                       # Platelet
)

scRNA=scRNA_harmony
missing_genes <- markers[!markers %in% rownames(scRNA)]
print(missing_genes)
valid_markers <- markers[markers %in% rownames(scRNA)]
unique_markers <- unique(valid_markers)
DotPlot(scRNA, features = unique_markers) +
  coord_flip() +
  theme_bw() +
  theme(panel.grid = element_blank(), axis.text.x = element_text(hjust = 1, vjust = 0.5)) +
  labs(x = NULL, y = NULL) + guides(size = guide_legend(order = 3)) +
  theme(axis.text.x = element_text(angle = 90)) +
  scale_color_gradientn(values = seq(0, 1, 0.2), colours = c('#008D83', '#B0A8B9', '#FF8066', '#C34A36'))


ggsave(filename="7-cluster-marker4.pdf", width = 180, height = 170, units = "mm")
####按照自己想的细胞顺序进行排序####
# 将 seurat_clusters 转换为因子并按照指定顺序排序
scRNA_harmony$seurat_clusters <- factor(scRNA_harmony$seurat_clusters,
                                        levels = c("0","1","4","7","12","5","9","13","8","16","6",
                                                   "11","18","21","2","3","10","14","17","19","20","15"
                                        ))

# 确保 markers 都存在于数据中
missing_genes <- markers[!markers %in% rownames(scRNA_harmony)]
print(missing_genes)  # 检查缺失基因
valid_markers <- markers[markers %in% rownames(scRNA_harmony)]
unique_markers <- unique(valid_markers)  # 去除重复基因

# 绘制 DotPlot
dotplot=DotPlot(scRNA_harmony, features = unique_markers, group.by = "seurat_clusters") +
  coord_flip() +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(hjust = 1, vjust = 0.5, angle = 90)) +
  labs(x = NULL, y = NULL) +
  guides(size = guide_legend(order = 3)) +
  scale_color_gradientn(values = seq(0, 1, 0.2), colours = c('#008D83', '#B0A8B9', '#FF8066', '#C34A36'))

# Display the plot
dotplot

ggsave(filename="22-cluster-marker24-12-31_updated-4.pdf", width = 240, height = 170, units = "mm")
ggsave(filename = "36-cluster-marker_updated-4.pdf", plot = dotplot, width = 240, height = 140, units = "mm")


#####cluster命名#####
celltype=data.frame(ClusterID=0:22,
                    celltype='unkown')
celltype[celltype$ClusterID %in% c(0,1,4,7,12),2]='Naive_CD4+_T_Cell'
celltype[celltype$ClusterID %in% c(5,9,13),2]='CD8+_T_Cell'
celltype[celltype$ClusterID %in% c(8,16),2]='CD16+Monocyte'
celltype[celltype$ClusterID %in% c(6,11,18),2]='B_Cell'
celltype[celltype$ClusterID %in% c(21),2]='Plasma_Cell'
celltype[celltype$ClusterID %in% c(2,3,10,14,17,19),2]='NK_Cell'
#celltype[celltype$ClusterID %in% c(12,32,29),2]='pDC'
#celltype[celltype$ClusterID %in% c(24),2]='Neutrophil'
celltype[celltype$ClusterID %in% c(20),2]='Platelet'
celltype[celltype$ClusterID %in% c(15),2]='NA'



celltype=data.frame(ClusterID=0:22,
                    celltype='unkown')
celltype[celltype$ClusterID %in% c(0,1,4,7,12),2]='T_Cell'
celltype[celltype$ClusterID %in% c(5,9,13),2]='T_Cell'
celltype[celltype$ClusterID %in% c(8,16),2]='CD16+Monocyte'
celltype[celltype$ClusterID %in% c(6,11,18),2]='B_Cell'
celltype[celltype$ClusterID %in% c(21),2]='Plasma_Cell'
celltype[celltype$ClusterID %in% c(2,3,10,14,17,19),2]='NK_Cell'
#celltype[celltype$ClusterID %in% c(12,32,29),2]='pDC'
#celltype[celltype$ClusterID %in% c(24),2]='Neutrophil'
celltype[celltype$ClusterID %in% c(20),2]='Platelet'
celltype[celltype$ClusterID %in% c(15),2]='T_Cell'
celltype
table(celltype$celltype)

sce.in=scRNA_harmony
#先加一列celltype所有值为空
sce.in@meta.data$celltype = "NA"
###注释
for(i in 1:nrow(celltype)){
  sce.in@meta.data[which(sce.in@meta.data$seurat_clusters == celltype$ClusterID[i]),'celltype'] <- celltype$celltype[i]}
table(sce.in@meta.data$celltype)


p <- DotPlot(sce.in, features = unique(markers),
             assay='RNA' ,group.by = 'celltype' )  + coord_flip()

p
# ggsave(plot=p, filename="check_marker_by_celltype.pdf")
ggsave(filename="check_marker_by_celltype.pdf",width = 210,height = 297,units = "mm",plot = p)



table(sce.in@meta.data$celltype,sce.in@meta.data$seurat_clusters)
sce=sce.in
p.dim.cell=DimPlot(sce, reduction = "tsne", group.by = "celltype",label = T,pt.size = 1)
p.dim.cell
ggsave(plot=p.dim.cell,filename="7DimPlot_tsne_celltype.pdf",width=9, height=7)
p.dim.cell=DimPlot(sce, reduction = "umap", group.by = "celltype",label = T,pt.size = 1,repel = T)
p.dim.cell
library(ggplot2)
ggsave(plot=p.dim.cell,filename="7DimPlot_umap_celltype.pdf",width=9, height=7)

DefaultAssay(sce) <- "RNA"
# Defaultassay设置为RNA",意味着接下来的分析将基于原始值
save(sce,file='6cluster细胞注释-改Rdata_分16cluster.Rdata')
# load("6cluster细胞注释-改Rdata_分16cluster.Rdata")


scRNA=sce
Idents(scRNA)=scRNA$celltype



# Generate the DotPlot with the new order

levels_in_celltype <- levels(as.factor(scRNA@meta.data$celltype))
print(levels_in_celltype)
markers <- c(
  "CD3D","CD3E","CCR7","IL7R",                # Naive CD4+ T Cell
  "CD8A","CD8B","GZMK","GZMB",                # CD8+ T Cell
  "FCGR3A","LYZ","VCAN",                      # CD16+ Monocyte
  "MS4A1","CD79A","CD19",                     # B Cell
  "MZB1","IGKC","JCHAIN",                     # Plasma Cell
  "GNLY","NKG7","PRF1",                       # Natural Killer (NK) Cell
  "IL3RA","GZMB","CLEC4C",                    # Dendritic Cell (pDC)
  "S100A8","S100A9","CSF3R",                  # Neutrophil
  "PPBP","PF4","ITGA2B"                       # Platelet
)

# 检查是否有重复水平
duplicated_levels <- levels_in_celltype[duplicated(levels_in_celltype)]
print(duplicated_levels)
scRNA@meta.data$celltype <- factor(scRNA@meta.data$celltype, levels = unique(scRNA@meta.data$celltype))


new_order <- c("Endothelial cells", "Fibroblasts", "Macro/Micro", "Astrocyte", "Neurons","Opc","Oligo","NA")
scRNA@meta.data$celltype <- factor(scRNA@meta.data$celltype, levels = new_order)

dotplot <- DotPlot(scRNA, features = unique(markers), group.by = "celltype") +
  coord_flip() +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(hjust = 1, vjust = 0.5, angle = 90)) +
  labs(x = NULL, y = NULL) +
  guides(size = guide_legend(order = 3)) +
  scale_color_gradientn(values = seq(0, 1, 0.2), colours = c('#008D83', '#B0A8B9', '#FF8066', '#C34A36'))

print(dotplot)

# Save the plot
ggsave(filename = "6-cluster-marker24-12-31-3.pdf", plot = dotplot, width = 110, height = 160, units = "mm")


load("scRNA_anno.Rdata")
###########差异分析热图#########
library(Seurat)
library(dplyr)
library(reticulate)
library(cowplot)
library(ggplot2)
library(viridis)
library(tidyr)
library(magrittr)
library(reshape2)
library(readr)
library(stringr)
# 输入数据预处理
scRNA=scRNA
# 查看细胞数量
table(scRNA@meta.data$celltype)
unique(scRNA@meta.data$celltype)
# 把annotations重新排个顺序，这里按照细胞数量排序
scRNA@meta.data$celltype <- factor(scRNA@meta.data$celltype,
                                   levels = c("Endothelial cells", "Fibroblasts", "Macro/Micro", "Astrocyte", "Neurons","Opc","Oligo"))

#提取设置好细胞类型的颜色
colour=c("#f8766d","#c49a00","#53b400", "#00c094","#00b6eb","#a58aff","#fb61d7")


Idents(scRNA) <- "celltype"
####FindAllMarkers差异分析/热图####
markers <- FindAllMarkers(object = scRNA, test.use="wilcox" ,
                          only.pos = TRUE,
                          logfc.threshold = 0.5,
                          min.pct = 0.5)
all.markers =markers %>% dplyr::select(gene, everything()) %>% subset(p_val<0.05)

#保存/读取
write.csv(all.markers,"分群之后all.markers.csv",row.names = T)
all.markers <- read.csv("分群之后all.markers.csv",header = T,row.names = 1)

#取top marker展示
top10 = all.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC) #将每个cluster lgFC排在前10的marker基因挑选出来
write.csv(top10,"分群之后top10.markers.csv",row.names = T)

top1 = all.markers %>% group_by(cluster) %>% top_n(n = 1, wt = avg_log2FC) #将每个cluster lgFC排在前10的marker基因挑选出来
markers <- top1$gene


####DoHeatmap热图####
top10 <- as.data.frame(top10)  # 将markers对象转换为数据框
markerdata <- ScaleData(scRNA, features = as.character(unique(top10$gene)), assay = "RNA")  # 对特定基因进行标准化
# 常规热图
DoHeatmap=DoHeatmap(markerdata,  # 使用DoHeatmap函数绘制热图
                    features = as.character(unique(top10$gene)),  # 指定要显示的特征
                    group.by = "celltype",  # 指定按照免疫细胞类型分组
                    assay = 'RNA')  # 指定分析的分析类型为RNA
ggsave(plot=DoHeatmap,filename="7DoHeatmap热图_celltype.pdf",width=18, height=15)

#修改颜色
DoHeatmap=DoHeatmap(markerdata,  # 再次绘制热图，附加自定义颜色和渐变色
                    features = as.character(unique(top10$gene)),
                    group.by = "celltype",
                    assay = 'RNA',
                    group.colors = c("#f8766d","#c49a00","#53b400", "#00c094","#00b6eb","#a58aff","#fb61d7")) +
  scale_fill_gradientn(colors = c("white","grey","firebrick3"))
ggsave(plot=DoHeatmap,filename="7DoHeatmap热图_celltype1.pdf",width=18, height=15)

####Umap图####
# 标准umap图
colour=c("#1E90FF","#7CFC00", "#808000","#FF00FF","#FA8072","#7B68EE","#800080","#A0522D","#D2B48C","#D2691E","#87CEEB","#40E0D0","#5F9EA0",
         
         "#FF1493","#0000CD","#008B8B")

p1 <- DimPlot(scRNA, reduction = "umap", group.by = "celltype",cols = colour,label = TRUE,repel = TRUE)
# group.by 指定分组umap图
p2 <- DimPlot(scRNA, reduction = "umap", group.by = "Group")

p1 + p2

pdf("scRNA.UMAP.pdf", height = 5, width = 12)
p1 + p2
dev.off()

pdf("scRNA.UMAP.pdf", height = 5,width = 12)

#split.by = "stim"拆分分组umap图
p3 <- DimPlot(scRNA, reduction = "umap", group.by = "celltype",
              split.by = "Group",cols = colour,repel = TRUE)
p3
pdf("scRNA.UMAP1.pdf", height = 5, width = 12)
p3
dev.off()

####细胞比例图orig.ident####


table(scRNA$Group)
prop.table(table(Idents(scRNA)))  # 计算细胞类型的比例
table(Idents(scRNA), scRNA$Group)  # 计算分组细胞类型的数量
Cellratio <- prop.table(table(Idents(scRNA), sce.in$Group), margin = 2)  # 计算细胞类型相对于原始身份的比例，并将结果转换为数据框
Cellratio <- as.data.frame(Cellratio)
colnames(Cellratio) <- c("celltype","sample","ratio")



colourCount = length(unique(Cellratio$celltype))  # 计算唯一的细胞类型数量
ggplot(Cellratio) +  # 创建 ggplot 对象并传入数据框 Cellratio
  geom_bar(aes(x = sample, y = ratio, fill = celltype), stat = "identity", width = 0.7, size = 0.5, colour = '#222222') +  # 添加条形图层，根据细胞类型填充颜色
  theme_classic() +  # 应用经典主题
  labs(x = 'Sample', y = 'Ratio') +  # 设置 x 和 y 轴标签
  coord_flip() +  # 翻转坐标轴
  theme(panel.border = element_rect(fill = NA, color = "black", size = 0.5, linetype = "solid"))+   # 设置面板边框属性
  scale_fill_manual(values = colour)  # 手动设置填充颜色

#直立柱状图
p3=ggplot(Cellratio) +
  geom_bar(aes(x =sample, y= ratio, fill = celltype),stat = "identity",width = 0.7,size = 0.5,colour = '#222222')+
  theme_classic() +
  labs(x='Sample',y = 'Ratio')+
  scale_fill_manual(values = colour)+
  theme(panel.border = element_rect(fill=NA,color="black", size=0.5, linetype="solid"))
p3
pdf("scRNA.细胞比例.pdf", height = 6, width = 5)
p3
dev.off()

#######危险因素评分##########
scRNA=sce
rm(sce.in,Cellratio)
gc()
#scRNA_Pla=readRDS('./scRNA_Pla1.RDS')##这个是T细胞类型的T_celltype,tissue_type,celltype
DimPlot(scRNA, label=TRUE,split.by = 'group')
T_scRNA=scRNA
seurat_clusters=scRNA1
load('scRNA_Nucleic_immunity.Rdata')

library(dplyr)
library(hdf5r)
library(Seurat)
library(data.table, lib.loc = "C:/Users/Administrator/AppData/Local/R/win-library/4.4")
library(Seurat)
library(SeuratDisk)
library(ggplot2)
library(Seurat)
library(dplyr)
library(ggplot2)
library(tidyr)
library(patchwork)
library(clustree)
library(clusterProfiler)

####先定比较的类型celltype####
table(scRNA@meta.data$celltype)
head(scRNA@meta.data)

## ## ## ## ## ## ## ## ## ## ## ## ## ## ## up
SLC_mediated=read.gmt('./危险因素/SLC_mediated_transmembrane_transport.gmt')
DotPlot(scRNA,features = SLC_mediated$gene)

SLC_mediated=list(SLC_mediated$gene)
names(SLC_mediated)='SLC_mediated'

scRNA=AddModuleScore(scRNA,features = SLC_mediated,name = 'SLC_mediated')

####
Transport_of=read.gmt('./危险因素/Transport_of_small_molecules.gmt')
DotPlot(scRNA,features = Transport_of$gene)

Transport_of=list(Transport_of$gene)
names(Transport_of)='Transport_of'

scRNA=AddModuleScore(scRNA,features = Transport_of,name = 'Transport_of')

####


Disorders_of=read.gmt('./危险因素/Disorders_of_transmembrane_transporters.gmt')
DotPlot(scRNA,features = Disorders_of$gene)

Disorders_of=list(Disorders_of$gene)
names(Disorders_of)='Disorders_of'

scRNA=AddModuleScore(scRNA,features = Disorders_of,name = 'Disorders_of')
####


Astrocytic_Glutamate=read.gmt('./危险因素/Astrocytic_Glutamate_Glutamine_Uptake_And_Metabolism.gmt')
DotPlot(scRNA,features = Astrocytic_Glutamate$gene)

Astrocytic_Glutamate=list(Astrocytic_Glutamate$gene)
names(Astrocytic_Glutamate)='Astrocytic_Glutamate'

scRNA=AddModuleScore(scRNA,features = Astrocytic_Glutamate,name = 'Astrocytic_Glutamate')


####

Neurotransmitter_uptake=read.gmt('./危险因素/Neurotransmitter_uptake_and_metabolism_In_glial_cells.gmt')
DotPlot(scRNA,features = Neurotransmitter_uptake$gene)

Neurotransmitter_uptake=list(Neurotransmitter_uptake$gene)
names(Neurotransmitter_uptake)='Neurotransmitter_uptake'

scRNA=AddModuleScore(scRNA,features = Neurotransmitter_uptake,name = 'Neurotransmitter_uptake')


#####
SLC_transporter=read.gmt('./危险因素/SLC_transporter_disorders.gmt')
DotPlot(scRNA,features = SLC_transporter$gene)

SLC_transporter=list(SLC_transporter$gene)
names(SLC_transporter)='SLC_transporter'

scRNA=AddModuleScore(scRNA,features = SLC_transporter,name = 'SLC_transporter')


#####
Glycogen_breakdown=read.gmt('./危险因素/Glycogen_breakdown.gmt')
DotPlot(scRNA,features = Glycogen_breakdown$gene)

Glycogen_breakdown=list(Glycogen_breakdown$gene)
names(Glycogen_breakdown)='Glycogen_breakdown'

scRNA=AddModuleScore(scRNA,features = Glycogen_breakdown,name = 'Glycogen_breakdown')


#####
Glutathione_synthesis=read.gmt('./危险因素/Glutathione_synthesis_and_recycling.gmt')
DotPlot(scRNA,features = Glutathione_synthesis$gene)

Glutathione_synthesis=list(Glutathione_synthesis$gene)
names(Glutathione_synthesis)='Glutathione_synthesis'

scRNA=AddModuleScore(scRNA,features = Glutathione_synthesis,name = 'Glutathione_synthesis')


######
mRNA_protein=read.gmt('./危险因素/mRNA_protein_and_metabolite_inducation_pathway_by_cyclosporin.gmt')
DotPlot(scRNA,features = mRNA_protein$gene)

mRNA_protein=list(mRNA_protein$gene)
names(mRNA_protein)='mRNA_protein'

scRNA=AddModuleScore(scRNA,features = mRNA_protein,name = 'mRNA_protein')


##############
Carnitine_metabolism=read.gmt('./危险因素/GOBP_CARNITINE_METABOLIC_PROCESS.v2025.1.Hs.gmt')
DotPlot(scRNA,features = Carnitine_metabolism$gene)

Carnitine_metabolism=list(Carnitine_metabolism$gene)
names(Carnitine_metabolism)='Carnitine_metabolism'

scRNA=AddModuleScore(scRNA,features = Carnitine_metabolism,name = 'Carnitine_metabolism')


p_up_markers <- DotPlot(scRNA, features = c('SLC_mediated1','Transport_of1','Disorders_of1','Astrocytic_Glutamate1','Neurotransmitter_uptake1','Carnitine_metabolism1','mRNA_protein1','Glutathione_synthesis1','Glycogen_breakdown1','SLC_transporter1'),
                        assay='RNA' ,group.by = 'celltype' )

p_up_markers


####提取数据####
data<-p_up_markers$data

colnames(data)

colnames(data)<-c("AverageExpression_unscaled","Precent Expressed","Features","celltype","Average Expression")

unique(data$`Precent Expressed`)
library(ggplot2)
####用ggplot画图####
p = ggplot(data,aes(celltype,Features,size = `Precent Expressed` ))+
  geom_point(shape=21,aes(fill= `Average Expression`),position =position_dodge(0))+
  theme_minimal()+xlab(NULL)+ylab(NULL) +
  scale_size_continuous(range=c(1,10))+theme_bw()+
  scale_fill_gradient(low = "grey", high = "#C34A36")+
  theme(legend.position = "right",legend.box = "vertical", #图例位置
        legend.margin=margin(t= 0, unit='cm'),
        legend.spacing = unit(0,"in"),
        axis.text.x  = element_text(color="black",size=16,angle = 45,
                                    vjust = 0.5, hjust=0.5),#x轴
        axis.text.y  = element_text(color="black",size=12),#y轴
        legend.text = element_text(size =12,color="black"),#图例
        legend.title = element_text(size =12,color="black"),#图例
        axis.title.y=element_text(vjust=1,
                                  size=16)
  )+labs(x=" ",y = "Features")+coord_flip();p

ggsave("血液危险因素活性分析_plot1.pdf", plot = p, device = "pdf", width = 7, height = 6)

rm(list = ls())
load("7cluster细胞注释-改Rdata_分16cluster.Rdata")
scRNA=sce
rm(sce)
####保护因素活性分析down######
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## down
#######
Metabolism_of=read.gmt('./保护因素/Metabolism_of_amino_acids_and_derivativesv2025.1.Hs.gmt')
DotPlot(scRNA,features = Metabolism_of$gene)

Metabolism_of=list(Metabolism_of$gene)
names(Metabolism_of)='Metabolism_of'

scRNA=AddModuleScore(scRNA,features = Metabolism_of,name = 'Metabolism_of')


#####
SLC_mediated=read.gmt('./保护因素/SLC_mediated_transmembrane_transport.gmt')
DotPlot(scRNA,features = SLC_mediated$gene)

SLC_mediated=list(SLC_mediated$gene)
names(SLC_mediated)='SLC_mediated'

scRNA=AddModuleScore(scRNA,features = SLC_mediated,name = 'SLC_mediated')


####
amino_acid=read.gmt('./保护因素/amino_acid_metabolism.gmt')
DotPlot(scRNA,features = amino_acid$gene)

amino_acid=list(amino_acid$gene)
names(amino_acid)='amino_acid'

scRNA=AddModuleScore(scRNA,features = amino_acid,name = 'amino_acid')


####
Transport_of_small=read.gmt('./保护因素/Transport_of_small_molecules.gmt')
DotPlot(scRNA,features = Transport_of_small$gene)

Transport_of_small=list(Transport_of_small$gene)
names(Transport_of_small)='Transport_of_small'

scRNA=AddModuleScore(scRNA,features = Transport_of_small,name = 'Transport_of_small')


####
Transport_of_inorganic=read.gmt('./保护因素/Transport_of_inorganic_cations.gmt')
DotPlot(scRNA,features = Transport_of_inorganic$gene)

Transport_of_inorganic=list(Transport_of_inorganic$gene)
names(Transport_of_inorganic)='Transport_of_inorganic'

scRNA=AddModuleScore(scRNA,features = Transport_of_inorganic,name = 'Transport_of_inorganic')

####
Arginine_and_proline=read.gmt('./保护因素/Arginine_and_proline_metabolism.gmt')
DotPlot(scRNA,features = Arginine_and_proline$gene)

Arginine_and_proline=list(Arginine_and_proline$gene)
names(Arginine_and_proline)='Arginine_and_proline'

scRNA=AddModuleScore(scRNA,features = Arginine_and_proline,name = 'Arginine_and_proline')


####
Amino_acid_transport=read.gmt('./保护因素/Amino_acid_transport_across_the_plasma_membrane.gmt')
DotPlot(scRNA,features = Amino_acid_transport$gene)

Amino_acid_transport=list(Amino_acid_transport$gene)
names(Amino_acid_transport)='Amino_acid_transport'

scRNA=AddModuleScore(scRNA,features = Amino_acid_transport,name = 'Amino_acid_transport')


####
Trans_sulfuration=read.gmt('./保护因素/Trans_sulfuration_one_carbon_metabolism.gmt')
DotPlot(scRNA,features = Trans_sulfuration$gene)

Trans_sulfuration=list(Trans_sulfuration$gene)
names(Trans_sulfuration)='Trans_sulfuration'

scRNA=AddModuleScore(scRNA,features = Trans_sulfuration,name = 'Trans_sulfuration')

####
Argininemia=read.gmt('./保护因素/Argininemia.gmt')
DotPlot(scRNA,features = Argininemia$gene)

Argininemia=list(Argininemia$gene)
names(Argininemia)='Argininemia'

scRNA=AddModuleScore(scRNA,features = Argininemia,name = 'Argininemia')

p_down_markers <- DotPlot(scRNA, features = c('Metabolism_of1','SLC_mediated1','amino_acid1','Transport_of_small1','Transport_of_inorganic1','Arginine_and_proline1','Amino_acid_transport1','Trans_sulfuration1','Argininemia1'),
                          assay='RNA' ,group.by = 'celltype' )

p_down_markers


####提取数据####
data<-p_down_markers$data

colnames(data)

colnames(data)<-c("AverageExpression_unscaled","Precent Expressed","Features","celltype","Average Expression")

unique(data$`Precent Expressed`)
library(ggplot2)
####用ggplot画图####
p = ggplot(data,aes(celltype,Features,size = `Precent Expressed` ))+
  geom_point(shape=21,aes(fill= `Average Expression`),position =position_dodge(0))+
  theme_minimal()+xlab(NULL)+ylab(NULL) +
  scale_size_continuous(range=c(1,10))+theme_bw()+
  scale_fill_gradient(low = "grey", high = "#008B74")+
  theme(legend.position = "right",legend.box = "vertical", #图例位置
        legend.margin=margin(t= 0, unit='cm'),
        legend.spacing = unit(0,"in"),
        axis.text.x  = element_text(color="black",size=16,angle = 45,
                                    vjust = 0.5, hjust=0.5),#x轴
        axis.text.y  = element_text(color="black",size=12),#y轴
        legend.text = element_text(size =12,color="black"),#图例
        legend.title = element_text(size =12,color="black"),#图例
        axis.title.y=element_text(vjust=1,
                                  size=16)
  )+labs(x=" ",y = "Features")+coord_flip();p

ggsave("血液保护因素活性分析_plot1.pdf", plot = p, device = "pdf", width = 7, height = 6)


library(ggplot2)
library(ggrepel)   # 用于避免标签重叠

## 1. 在 df1 里增加一列显著性标记
df1 <- df1 %>%
  mutate(
    significance = case_when(
      p_val_adj < 0.05 & avg_log2FC >=  0.25 ~ "Up",
      p_val_adj < 0.05 & avg_log2FC <= -0.25 ~ "Down",
      TRUE ~ "NotSig"
    )
  )

## 2. 准备 TOP 基因标签（可选）
top_genes <- df1 %>%
  filter(significance != "NotSig") %>%
  slice_max(abs(avg_log2FC), n = 20)   # 取 |logFC| 最大的 20 个

## 3. 画图
volcano <- ggplot(df1, aes(x = avg_log2FC, y = -log10(p_val_adj), colour = significance)) +
  geom_point(alpha = 0.8, size = 1.5) +
  scale_color_manual(values = c("Up" = "#E64B35", "Down" = "#4DBBD5", "NotSig" = "grey70")) +
  geom_vline(xintercept = c(-0.25, 0.25), linetype = "dashed", colour = "grey50") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", colour = "grey50") +
  geom_label_repel(data = top_genes,
                   aes(label = rownames(top_genes)),
                   colour = "black",
                   box.padding = 0.15,
                   segment.color = "grey50",
                   max.overlaps = 20) +
  labs(title = "Plasma cell: AAmetabolism high vs low",
       x = "log2 fold change",
       y = "-log10 adjusted p-value") +
  theme_minimal(base_size = 14) +
  theme(legend.title = element_blank())

print(volcano)

## 4. 保存
ggsave("Plasma cell-volcano_AAmetabolism.pdf", volcano, width = 7, height = 6)
write.csv(df1, file = "AAmetabolism low-vs-high_Plasma cell_gene.csv")




