require('NanoStringNorm')
library("stringr")

cohorts_t = read.table("~/Koop_Kulbe/Misc/cohorts.tsv",sep="\t",header =T, stringsAsFactors = F)
#raw_data  = read.markup.RCC( rcc.path = "~/Koop_Kulbe/Raw_data/Nanostring_final_samples/")
raw_data  = read.markup.RCC( rcc.path = "~/Koop_Kulbe/Raw_data/24_samples/Malignant_Benign//")
#raw_data  = read.markup.RCC( rcc.path = "~/Koop_Kulbe/Raw_data/24_samples/Healthy_Benign//")

sample_names = colnames(raw_data$x)[-seq(3)]
sample_names = str_replace(sample_names, pattern = "^X","")
colnames(raw_data$x) = str_replace(colnames(raw_data$x), pattern = "^X","")
sample_names = str_replace_all(sample_names, pattern = "\\.","_")

#s_match = match( sample_names, cohorts_t$New_Id, nomatch = 0)
s_match = match( sample_names, cohorts_t$Sample_Id, nomatch = 0)
cohort_names = cohorts_t$New_Id[s_match]

meta_data = data.frame( cohorts_t$Malignant[ s_match ] )
rownames( meta_data ) = colnames(raw_data$x)[-seq(3)]#cohorts_t$New_Id[ s_match]
colnames(meta_data) = "Group"
rownames(meta_data) = str_replace( rownames(meta_data), pattern = "^X", "" )

#meta_data$Group[ meta_data$Group == 1 ] = 2
#meta_data$Group[ is.na(meta_data$Group)] = 1

#norm.comp.results.test = norm.comp(raw_data, verbose = F)
eset = NanoStringNorm::NanoStringNorm( 
  raw_data,
  CodeCount.methods = "sum",
  Background.methods = "mean.2sd",
  SampleContent.methods = "none",
  OtherNorm.methods = "vsn",
  take.log = T,
  round.values = T,
  return.matrix.of.endogenous.probes = F,
  traits = meta_data,
  verbose = T
)

m    = matrix( as.character(unlist( eset$normalized.data)), nrow=  dim(eset$normalized.data)[1], ncol = dim(eset$normalized.data)[2])
info = m[,seq(3)]
data = matrix( as.double(m[,-seq(3)]), nrow=  dim(eset$normalized.data)[1], ncol = dim(eset$normalized.data)[2]-3)
data = round(data,1)

rownames(data) = rownames(eset$normalized.data)
col_labels = str_replace( colnames(eset$normalized.data)[-seq(3)], pattern = "^X", "") 

colnames(data) = col_labels
res  = cbind( info,data )
res = cbind(rownames(eset$normalized.data), res)

pure_data = eset$normalized.data[,c(-1,-2,-3)]
info = eset$normalized.data[,1]
pure_data = pure_data[ info == "Endogenous",]
#eset$normalized.data = eset$normalized.data[ eset$normalized.data$Code.Class == "Endogenous",]
#eset$gene.summary.stats.norm = eset$gene.summary.stats.norm[ eset$normalized.data$Code.Class == "Endogenous",]
#eset$gene.summary.stats.norm = eset$gene.summary.stats.norm[ !str_detect(rownames(eset$gene.summary.stats.norm), pattern = "POS_"),]
#eset$gene.summary.stats.norm = eset$gene.summary.stats.norm[ !str_detect(rownames(eset$gene.summary.stats.norm), pattern = "NEG_"),]
#eset$gene.summary.stats.norm = eset$gene.summary.stats.norm[ !str_detect(rownames(eset$gene.summary.stats.norm), pattern = "GAPDH"),]
#eset$gene.summary.stats.norm = eset$gene.summary.stats.norm[ !str_detect(rownames(eset$gene.summary.stats.norm), pattern = "B2M"),]
  
colnames( pure_data ) = cohorts_t$New_Id[s_match]

controls = which( meta_data$Group == 1 )
case     = which( meta_data$Group == 2 )

mean_controls = rowMeans( pure_data[, controls] )
mean_case     = rowMeans( pure_data[, case] )

dif = as.double(mean_case - mean_controls)

pure_data = pure_data[ order( dif ) , ]
controls  = controls[order(meta_data$Group)]
case      = case[order(meta_data$Group)]

export_pure_data = cbind( rownames(pure_data) ,round( dif[order(dif)], 1 ),round(mean_controls[order(dif)],1),round(mean_case[order(dif)],1),pure_data)
colnames( export_pure_data ) = c( c("Gene","Dif","Mean_Ctrl","Mean_Case")  ,colnames(pure_data))

export_pure_data = export_pure_data[order(abs(as.double(as.character(export_pure_data[,2]))), decreasing = T),]

#write.table("~/Koop_Kulbe/Results/Healthy_Benign//differential_expression_blood.tsv",x=cbind( info,eset$gene.summary.stats.norm ),sep="\t",row.names =F, quote =F)
