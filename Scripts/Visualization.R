draw_colnames_45 <- function (coln, ...) {
  m = length(coln)
  x = (1:m)/m - 1/2/m
  grid.text(coln, x = x, y = unit(0.96, "npc"), vjust = .5, 
            hjust = 1, rot = 45, gp = gpar(...)) ## Was 'hjust=0' and 'rot=270'
}

## For pheatmap_1.0.8 and later:
draw_colnames_45 <- function (coln, gaps, ...) {
  coord = pheatmap:::find_coordinates(length(coln), gaps)
  x = coord$coord - 0.5 * coord$size
  res = textGrob(coln, x = x, y = unit(1, "npc") - unit(3,"bigpts"), vjust = 0.5, hjust = 1, rot = 45, gp = gpar(...))
  return(res)}
assignInNamespace(x="draw_colnames", value="draw_colnames_45",ns=asNamespace("pheatmap"))


#pure_data = read.table("~/Koop_Kulbe/Results/Malignant_Benign/mean_housekeeping_normalized_data_blood.tsv",sep ="\t",stringsAsFactors = F, header = T, row.names = 1)
#traits = str_detect(colnames(pure_data), pattern = "B")
#traits[ ! traits ] = 2
traits = eset$traits
df = as.data.frame(traits)
df$Group[df$Group == "1"] = "Benign"
df$Group[df$Group == "2"] = "Malignant"
df$Cohort = df$Group

# visualization
#pdf("~/Koop_Kulbe/Results/Healthy_Benign/Volcano.pdf")
    Plot.NanoStringNorm(eset,plot.type = "volcano")
#dev.off()

subtype_top_bar = as.matrix( c(
  colorRampPalette(colors = c("blue"))( nrow( traits ) ) 
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
  cor( pure_data ),
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
g <- g + theme(legend.direction = 'horizontal', legend.position = 'top')
g# + theme(plot.margin = unit(c(0,0,0,0), "cm"))

### additional vis

cor_mat = cor(pure_data)
pcr_mat = prcomp(t(cor_mat))

#aka3 = list(Subtype = c(Alpha = "red", Beta = "blue", Gamma = "brown", Delta = "green", Epsilon = "black", Unknown = "gray", Whole = "purple"))

pdf("~/Koop_Kulbe/Results/Healthy_Malignant/Correlation_heatmap.pdf")
    pheatmap::pheatmap(
      cor_mat,
      annotation_col = df[c("Cohort")],
      show_rownames = F,
      show_colnames = T
        
      #annotation_colors = aka3,
      #color =  colorRampPalette(rev(RColorBrewer::brewer.pal(n = 11, name = "RdYlBu")))(100)
    )
dev.off()

#pdf("~/Koop_Kulbe/Results/Healthy_Malignant/Heatmap_absolute.pdf")
diff_mat = pure_data- rowMeans(pure_data)
selector_variance = apply(diff_mat, FUN=var, MARGIN = 1)
diff_mat = diff_mat[
  rownames(diff_mat) %in% names(sort(selector_variance,decreasing = T)[1:10]),
]
    pheatmap::pheatmap(
      diff_mat,
      annotation_col = df[c("Cohort")],
      show_rownames = T,
      show_colnames = T,
      cluster_cols = F,
      cluster_rows = T,
      gaps_col = 12
      #annotation_colors = aka3,
      #color =  colorRampPalette(rev(RColorBrewer::brewer.pal(n = 11, name = "RdYlBu")))(100)
    )
#dev.off()

### boxplot
    
r_mat = reshape::melt(pure_data)
r_mat$Cohort = df$Cohort[match(r_mat$X2, rownames(df))]
colnames(r_mat) = c("Gene","Sample","Expression","Cohort")
r_mat = r_mat[r_mat$Gene %in% names(selector_variance)[1:10],]

gene_box = ggplot(data=r_mat, aes(x = Gene,y = Expression, fill = Cohort)) 
gene_box = gene_box + geom_boxplot() + theme(legend.direction = 'horizontal', legend.position = 'top')
gene_box + coord_flip()

### p_value_plot

library(stringr)
library(ggplot2)

p_tab = read.table("~/Koop_Kulbe/Misc/Boxplot_data.tsv",sep ="\t", header = T, stringsAsFactors = F)
p_tab$P_value = str_replace(p_tab$P_value, pattern = ",", ".")
p_tab$P_value = as.double(p_tab$P_value)

p_tab$P_value = log10(p_tab$P_value) * -1
new_lev = unique( factor(p_tab$Gene) )[ order(p_tab$P_value[p_tab$Group == "OSE_v_OVCA"], decreasing = T) ]
new_lev = new_lev[ ! is.na(new_lev)]
p_tab$Gene = factor(p_tab$Gene, levels = new_lev)

p = ggplot(
    data = p_tab,
    aes( x = Gene, y = P_value, fill = Group ) 
) + geom_bar( stat = "identity", position = position_dodge())
p = p + geom_hline(yintercept = 1.3)
p = p + theme( legend.position = "top")
p = p + annotate("text", x = 9.75, y = 1.25, vjust = -1, label = "Significant")
p = p + annotate("text", x = 9.75, y = 0.75, vjust = -1, label = "Not significant")
p
