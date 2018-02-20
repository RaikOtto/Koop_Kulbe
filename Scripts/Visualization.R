
# visualization
pdf("~/Koop_Kulbe/Results/Healthy_Benign/Volcano.pdf")
Plot.NanoStringNorm(eset,plot.type = "volcano")
dev.off()
subtype_top_bar = as.matrix( c(
  colorRampPalette(colors = c("blue"))( dim( traits )[1] ) 
))

subtype_top_bar[ case, 1  ] = "yellow"

colnames(subtype_top_bar) = c("Cohort")

logFC_side_bar = t(
  c(
    colorRampPalette(colors = c("green"))( length( dif[ dif < 0 ]) ),
    colorRampPalette(colors = c("yellow"))( length( dif[ dif == 0 ]) ),
    colorRampPalette(colors = c("red"))( length( dif[ dif > 0 ]) )
  )
)
rownames(logFC_side_bar) = c("FoldChange")
mcols = colorRampPalette(colors = c("green","black","red"))( 75 )
#library(devtools)
#install_github("ggbiplot", "vqv")

eset_pca = prcomp(
  t( pure_data ),
  center = TRUE,
  scale. = TRUE
) 
summary(eset_pca)

cohorts_bar = subtype_top_bar
cohorts_bar[cohorts_bar =="yellow"  ] = "Malignant"
cohorts_bar[cohorts_bar !="Malignant"  ] = "Benign"

library(ggbiplot)
g = ggbiplot( 
  eset_pca, 
  obs.scale = 1, 
  var.scale = 1,
  groups = as.character( cohorts_bar ),
  ellipse = TRUE, 
  circle = TRUE,
  labels = colnames(pure_data)
)
g <- g + scale_color_discrete(name = '')
g <- g + theme(legend.direction = 'horizontal', 
               legend.position = 'top')
print(g)

meta_data = traits
rownames(meta_data) = cohorts_t$New_Id[s_match]
sce <- scater::newSCESet( exprsData = pure_data)
pData(sce) = meta_data

scater::plotPCASCESet(
  sce,
  colour_by = "Group"
)
jpg("~/MAPTor_Net_RNA_data/Results/tSNE_classification_KNN4.jpeg")
scater::plotTSNE(
  sce,
  colour_by = "Group",
  #perplexity = 5,
  check_duplicates = FALSE
)
dev.off()

df = meta_data
#df$Group[is.na(df$Group)]= "healthy"
df$Group[ df$Group == "1"]= "benign"
df$Group[ df$Group == "2"]= "malignant"
rownames(df) = colnames(cor_mat)

pdf("~/Koop_Kulbe/Results/Malignant_Benign//PCA.pdf")
ggbiplot::ggbiplot(
  pcr_mat,
  obs.scale = 1, 
  var.scale = 1, 
  labels.size = 4,
  alpha = 1,
  groups = df$Group,
  ellipse = TRUE, 
  circle = TRUE,
  var.axes = F
)
dev.off()
### additional vis

cor_mat = cor(pure_data)
pcr_mat = prcomp(t(cor_mat))

#aka3 = list(Subtype = c(Alpha = "red", Beta = "blue", Gamma = "brown", Delta = "green", Epsilon = "black", Unknown = "gray", Whole = "purple"))

pdf("~/Koop_Kulbe/Results/Healthy_Malignant/Correlation_heatmap.pdf")
pheatmap::pheatmap(
  cor_mat,
  annotation_col = data.frame( df ),
  show_rownames = F,
  show_colnames = F
  #annotation_colors = aka3,
  #color =  colorRampPalette(rev(RColorBrewer::brewer.pal(n = 11, name = "RdYlBu")))(100)
)
dev.off()

pdf("~/Koop_Kulbe/Results/Healthy_Malignant/Heatmap_absolute.pdf")
pheatmap::pheatmap(
  pure_data,
  annotation_col = data.frame( df ),
  show_rownames = F,
  show_colnames = F
  #annotation_colors = aka3,
  #color =  colorRampPalette(rev(RColorBrewer::brewer.pal(n = 11, name = "RdYlBu")))(100)
)
dev.off()
