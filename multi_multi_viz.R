################################################################################
## Sample multiplexing for joint RNA-ATAC multiomic profiling using MULTI-seq ##
## Nature Protocols (2026) ###### Chris McGinnis, PhD ###### May 2026 ##########
################################################################################

library(Seurat)
library(Signac)
library(deMULTIplex)
library(deMULTIplex2)
library(ggplot2)

##############
## Figure 4 ##
##############
load('seu_preQC.Robj')
load('seu_postQC.Robj')
load('seu_sub1.Robj')
load('seu_final.Robj')
load('multi_umap_good.Robj')
load('multi_umap_poor.Robj')
load('reclass_res.Robj')
load('dub_enrich.Robj')
load('neg_enrich.Robj')

## Fig. 4A: Multiome QC metric scatter plots
# load('seu_preQC.Robj')
DensityScatter(seu_preQC, x = 'nCount_ATAC', y = 'TSS.enrichment', log_x = TRUE, quantiles = TRUE)
DensityScatter(seu_preQC, x = 'nCount_ATAC', y = 'nCount_RNA', log_x = TRUE, log_y = TRUE, quantiles = TRUE)
DensityScatter(seu_preQC, x = 'pmito', y = 'nCount_RNA', log_y = TRUE, quantiles = TRUE)

## Fig. 4B: High-quality MULTI-seq barcode UMAP
# load('multi_umap_good.Robj')
ggplot(multi_umap_good, aes(x=UMAP_1, y=UMAP_2, color=multi)) + geom_point(size=0.1) + theme_classic() + scale_color_manual(values=c('black','grey','lightcoral','darkred','orchid4','cadetblue3','dodgerblue','seagreen','palegreen3','goldenrod')) + theme(legend.position = 'none')

## Fig. 4C: Poor-quality MULTI-seq barcode UMAP colored by pre-rescue sample classifications and MULTI nUMI ECDF
# load('multi_umap_poor.Robj')
ggplot(multi_umap_poor, aes(x=UMAP_1, y=UMAP_2, color=multi)) + geom_point(size=0.1) + theme_classic() + scale_color_manual(values=c('black','grey','goldenrod','lightcoral','steelblue4','seagreen')) + theme(legend.position = 'none')
ggplot(multi_umap_poor, aes(x=UMAP_1, y=UMAP_2, color=numi_ecdf)) + geom_point(size=0.1) + theme_classic() + scale_color_viridis_c(option = "B")+theme(legend.position = 'none')

## Fig. 4D: Negative cell rescue classification threshold inflection point scatter plot
# load('reclass_res.Robj')
ggplot(reclass_res[which(reclass_res$ClassStability <= 20), ], aes(x=ClassStability, y=MatchRate_mean)) + geom_point() + theme_classic() + geom_vline(xintercept = 6)

## Fig. 4E: Poor-quality MULTI-seq barcode UMAP colored by post-rescue sample classifications and scaled classification stablity
# load('multi_umap_poor.Robj')
ggplot(multi_umap_poor, aes(x=UMAP_1, y=UMAP_2, color=multi_res)) + geom_point(size=0.1) + theme_classic() + scale_color_manual(values=c('black','grey','goldenrod','lightcoral','steelblue4','seagreen')) + theme(legend.position = 'none')
ggplot(multi_umap_poor, aes(x=UMAP_1, y=UMAP_2)) + geom_point(size=0.1, color='grey70') + theme_classic() + 
  geom_point(data = multi_umap_poor[which(multi_umap_poor$multi == 'negative'), ], aes(color = cs), size=0.1) + scale_color_viridis_c(option = "B")

## Fig. 4G: Cluster-level doublet frequency plots: Joint RNA-ATAC UMAP and scatter plot 
# load('seu_preQC.Robj')
# load('dub_enrich.Robj)
FeaturePlot(seu_preQC, 'dubfreq', max.cutoff = 'q95', reduction = 'umap') + NoAxes() + scale_color_viridis_c(option = "B") + theme(plot.title = element_blank()) + NoLegend()
jitter_pos <- position_jitter(width=0.1, height=0, seed=1)
ggplot(dub_enrich, aes(x=multi, y=freq)) + geom_jitter(position = jitter_pos) + geom_hline(yintercept = mean(dub_enrich$freq)+sd(dub_enrich$freq), lty=2) + 
  geom_label_repel(data = dub_enrich, aes(label = cluster_dub), position=jitter_pos) + theme_classic() 

## Fig. 4H: Cluster-level negative frequency plots: Joint RNA-ATAC UMAP and scatter plot 
# load('seu_preQC.Robj')
# load('neg_enrich.Robj)
FeaturePlot(seu_preQC, 'negfreq', max.cutoff = 'q95', reduction = 'umap') + NoAxes() + scale_color_viridis_c(option = "B") + theme(plot.title = element_blank()) + NoLegend()
jitter_pos <- position_jitter(width=0.1, height=0, seed=1)
ggplot(neg_enrich, aes(x=multi, y=freq)) + geom_jitter(position = jitter_pos) + geom_hline(yintercept = mean(neg_enrich$freq)+sd(neg_enrich$freq), lty=2) + 
  geom_label_repel(data = neg_enrich, aes(label = cluster_neg), position=jitter_pos) + theme_classic() 

## Fig. 4I: Post-QC joint RNA-ATAC UMAP colored by MULTI-seq classification with marker gene violin plots
# load('seu_postQC.Robj)
# seu_postQC <- SetIdent(seu_postQC, value=seu_sub1$multi)
DimPlot(seu_postQC, cols=c('grey','black','goldenrod','palegreen3','seagreen','dodgerblue','cadetblue3','orchid4','darkred','lightcoral'))+NoLegend()+NoAxes()+theme(plot.title = element_blank())
VlnPlot(seu_postQC, pt.size = 0, idents=paste0('Tag',1:8), features=rev(c('CPED1','FOXP2','PSAT1','DPF3','FBXL7')), stack=T, fill.by='ident', cols=c('goldenrod','palegreen3','seagreen','dodgerblue','cadetblue3','orchid4','darkred','lightcoral')) + NoLegend() 

## Fig. 4J: Subset joint RNA-ATAC UMAP colored by MULTI-seq classification
# load('seu_sub1.Robj')
# seu_sub1 <- SetIdent(seu_sub1, value=seu_sub1$multi)
DimPlot(seu_sub1, cells=colnames(seu_sub1)[which(seu_sub1$multi %in% c('multiplet','negative',paste0('Tag',4:7)))], cols=c('grey','cadetblue3','dodgerblue','seagreen','palegreen3','black'), order='multiplet')+NoLegend()+NoAxes()+theme(plot.title = element_blank())
# seu_sub1 <- SetIdent(seu_sub1, value=seu_sub1$predicted_dub)
VlnPlot(seu_sub1, pt.size = 0, idents=c('dub',paste0('Tag',4:7)),features=rev(c('CPED1','FOXP2','PSAT1','DPF3','FBXL7')), stack=T, fill.by='ident', cols=rev(c('black','cadetblue3','dodgerblue','seagreen','palegreen3'))) + NoLegend() 

## Fig. 4K: Final cleaned UMAP
# load('seu_final.Robj')
DimPlot(seu_final, group.by = 'multi', cols=c('goldenrod','palegreen3','seagreen','dodgerblue','cadetblue3','orchid4','darkred','lightcoral'))+NoLegend()+NoAxes()+theme(plot.title = element_blank())
