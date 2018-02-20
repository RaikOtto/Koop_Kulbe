require('NanoStringNorm')
library("stringr")

cohorts_t = read.table("~/Koop_Kulbe/Misc/cohorts.tsv",sep="\t",header =T, stringsAsFactors = F)
#raw_data  = read.markup.RCC( rcc.path = "~/Koop_Kulbe/Raw_data/24_samples/Malignant_Benign/")
#raw_data  = read.markup.RCC( rcc.path = "~/Koop_Kulbe/Raw_data/24_samples/Healthy_Malignant//")
raw_data  = read.markup.RCC( rcc.path = "~/Koop_Kulbe/Raw_data/24_samples/Healthy_Benign//")

sample_names = colnames(raw_data$x)[-seq(3)]
sample_names = str_replace(sample_names, pattern = "^X","")
sample_names = str_replace_all(sample_names, pattern = "\\.","_")
s_match = match(sample_names, cohorts_t$Sample_Id, nomatch = 0)
cohort_names = cohorts_t$New_Id[s_match]

meta_data = data.frame( cohorts_t$Malignant[ s_match ] )
rownames( meta_data ) = colnames(raw_data$x)[-seq(3)]#cohorts_t$New_Id[ s_match]
colnames(meta_data) = "Group"
meta_data$Group[ meta_data$Group == 1 ] = 2
meta_data$Group[ is.na(meta_data$Group)] = 1

norm.comp.results.test = norm.comp(raw_data, verbose = F)
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

pure_data = as.character(res)[-seq(dim(eset$normalized.data)[1]*4)]
pure_data = matrix( as.double( pure_data ), nrow = dim(eset$normalized.data)[1] )
rownames( pure_data ) = info[ ,2 ]
colnames( pure_data ) = cohorts_t$New_Id[s_match]

pure_data = pure_data[ info[,1] == "Endogenous"  ,]

controls = which( meta_data$Group == 1 )
case     = which( meta_data$Group == 2 )

mean_controls = rowMeans( pure_data[, controls] )
mean_case     = rowMeans( pure_data[, case] )

dif = as.double(mean_case - mean_controls)

pure_data = pure_data[ order( dif ) ,order(meta_data$Group)]
controls  = controls[order(meta_data$Group)]
case      = case[order(meta_data$Group)]

export_pure_data = cbind( rownames(pure_data) ,round( dif[order(dif)], 1 ),round(mean_controls[order(dif)],1),round(mean_case[order(dif)],1),pure_data)
colnames( export_pure_data ) = c( c("Gene","Dif","Mean_Ctrl","Mean_Case")  ,colnames(pure_data))

export_pure_data = export_pure_data[order(abs(as.double(as.character(export_pure_data[,2]))), decreasing = T),]

#write.table("~/Koop_Kulbe/Results/Malignant_Benign/mean_housekeeping_normalized_data_blood.tsv",x=export_pure_data,sep="\t",row.names =F, quote =F)
#write.table("~/Koop_Kulbe/Results/Malignant_Benign/differential_expression_blood.tsv",x=cbind( info,eset$gene.summary.stats.norm ),sep="\t",row.names =F, quote =F)
#write.table("~/Koop_Kulbe/Results/Healthy_Malignant//mean_housekeeping_normalized_data_blood.tsv",x=export_pure_data,sep="\t",row.names =F, quote =F)
#write.table("~/Koop_Kulbe/Results/Healthy_Malignant//differential_expression_blood.tsv",x=cbind( info,eset$gene.summary.stats.norm ),sep="\t",row.names =F, quote =F)
write.table("~/Koop_Kulbe/Results/Healthy_Benign/mean_housekeeping_normalized_data_blood.tsv",x=export_pure_data,sep="\t",row.names =F, quote =F)
write.table("~/Koop_Kulbe/Results/Healthy_Benign//differential_expression_blood.tsv",x=cbind( info,eset$gene.summary.stats.norm ),sep="\t",row.names =F, quote =F)
