
#pure_data = read.table("~/Koop_Kulbe/Results/Malignant_Benign/mean_housekeeping_normalized_data_blood.tsv",sep ="\t",stringsAsFactors = F, header = T, row.names = 1)
pure_data = pure_data[,4:ncol(pure_data)]
#traits = str_detect(colnames(pure_data), pattern = "B")
#traits[ ! traits ] = 2
traits = eset$traits

# visualization
pdf("~/Koop_Kulbe/Results/Healthy_Benign/Volcano.pdf")
    Plot.NanoStringNorm(eset,plot.type = "volcano")
dev.off()

subtype_top_bar = as.matrix( c(
  colorRampPalette(colors = c("blue"))( length( traits ) ) 
))
subtype_top_bar[ traits == "1"  ] = "yellow"
colnames(subtype_top_bar) = c("Cohort")

cohorts_bar = subtype_top_bar
cohorts_bar[cohorts_bar =="yellow"  ] = "Malignant"
cohorts_bar[cohorts_bar !="Malignant"  ] = "Benign"


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

library(ggbiplot)
g = ggbiplot( 
  eset_pca,
  var.axes = F,
  #obs.scale = 1, 
  #var.scale = 1,
  groups = as.character( cohorts_bar ),
  ellipse = TRUE, 
  #circle = TRUE,
  labels = colnames(pure_data)
)
g <- g + scale_color_discrete(name = '')
g <- g + theme(legend.direction = 'horizontal', 
               legend.position = 'top')
print(g)

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
