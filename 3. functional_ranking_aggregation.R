###Robust rank aggregation of IAV related MAGeCK output

rm(list = ls())
setwd("E:/Research_Project/Drug_screening/genome_wide_screening/Flu")


###For pro-IAV factors screening studies

###1 CRISPR_KO screening
###PMID29642015 (screening round from 1 to 5)
PMID29642015_rd1 = read.table("./PMID29642015/GSE111166/mageck/Rd1.gene_summary.txt", header = T, sep = "\t")
PMID29642015_rd2 = read.table("./PMID29642015/GSE111166/mageck/Rd2.gene_summary.txt", header = T, sep = "\t")
PMID29642015_rd3 = read.table("./PMID29642015/GSE111166/mageck/Rd3.gene_summary.txt", header = T, sep = "\t")
PMID29642015_rd4 = read.table("./PMID29642015/GSE111166/mageck/Rd4.gene_summary.txt", header = T, sep = "\t")
PMID29642015_rd5 = read.table("./PMID29642015/GSE111166/mageck/Rd5.gene_summary.txt", header = T, sep = "\t")

###PMID34871331 (parallel libraries A and B)
PMID34871331_A = read.table("./PMID34871331/PRJNA777037/2.mageck/mageck_A.gene_summary.txt", header = T, sep = "\t")
PMID34871331_B = read.table("./PMID34871331/PRJNA777037/2.mageck/mageck_B.gene_summary.txt", header = T, sep = "\t")

###PMID35354039
PMID35354039 = read.csv("./PMID35354039/MAGECK_full_table.csv", header = T)

###PMID37652009（virus encompass sgRNA strategy）
PMID37652009_p1 = read.table("./PMID37652009/mageck/output/mageck_p1.gene_summary.txt", header = T, sep = "\t")
PMID37652009_p2 = read.table("./PMID37652009/mageck/output/mageck_p2.gene_summary.txt", header = T, sep = "\t")
PMID37652009_p3 = read.table("./PMID37652009/mageck/output/mageck_p3.gene_summary.txt", header = T, sep = "\t")
PMID37652009_p4 = read.table("./PMID37652009/mageck/output/mageck_p4.gene_summary.txt", header = T, sep = "\t")
PMID37652009_p5 = read.table("./PMID37652009/mageck/output/mageck_p5.gene_summary.txt", header = T, sep = "\t")

###PMID32839537
PMID32839537_IAVG = read.csv("./PMID32839537_ko/IAV-G_2_IAV-G.csv", header = T)
PMID32839537_WSN = read.csv("./PMID32839537_ko/IAV-G_2_WSN.csv", header = T)
PMID32839537_cs35m_gfp = read.csv("./PMID32839537_ko/SC35M Flu-GFP.csv", header = T)

###2 RNAi screening
###PMID34556855(two strains with or without type I IFN stimulation)
PMID34556855_H3N2_IFN = read.csv("./PMID34556855_RNAi/H3N2 +IFN.csv", header = T)
PMID34556855_H3N2_noIFN = read.csv("./PMID34556855_RNAi/H3N2 -IFN.csv", header = T)
PMID34556855_H5N1_IFN = read.csv("./PMID34556855_RNAi/H5N1 +IFN.csv", header = T)
PMID34556855_H5N1_noIFN = read.csv("./PMID34556855_RNAi/H5N1 -IFN.csv", header = T)


###3 FACS screening
###PMID31919360
PMID31919360_p = read.csv("./PMID31919360/primary.csv", header = T)
PMID31919360_s = read.csv("./PMID31919360/secondary.csv", header = T)









###For antiIAV factors screening studies

###1 CRISPRa screening
###PMID28813663
PMID28813663_A = read.table("./PMID28813663/GSE89371/2.mageck/1.output/GSE89371_screenA.gene_summary.txt", header = T, sep = "\t")
PMID28813663_B = read.table("./PMID28813663/GSE89371/2.mageck/1.output/GSE89371_screenB.gene_summary.txt", header = T, sep = "\t")
PMID28813663_C = read.table("./PMID28813663/GSE89371/2.mageck/1.output/GSE89371_screenC.gene_summary.txt", header = T, sep = "\t")

###2 FACS screening
###PMID32699298
PMID32699298 = read.csv("./PMID32699298_scirep_anitIAV/mageck.csv", header = T)







###For both pro and anti-IAV screening

###1 FACS screening
###PMID36129973(ranked by logFC)
PMID36129973 = read.csv("./PMID36129973/list.csv", header = T)
PMID36129973_proIAV = PMID36129973[which(PMID36129973$log2FoldChange < 0),]
PMID36129973_antiIAV = PMID36129973[which(PMID36129973$log2FoldChange > 0),]

###2 RNAi screening (META analysis study)
###PMID26651948
zrra = read.csv("./PMID26651948_MetaAnalysis_4Sets/z_score.csv", header = T)
zrra_pro = zrra[which(zrra$Z_RSA < 0),]
zrra_anti = zrra[which(zrra$Z_RSA > 0),]







###Pro-IAV tables aggregation with RRA
library(RobustRankAggreg)

###Build RRA object
aggr = list()

###Add ranking table that ranked by gene functional importance
aggr$PMID29642015_rd1 = PMID29642015_rd1[order(PMID29642015_rd1$pos.score),]$id
aggr$PMID29642015_rd2 = PMID29642015_rd2[order(PMID29642015_rd2$pos.score),]$id
aggr$PMID29642015_rd3 = PMID29642015_rd3[order(PMID29642015_rd3$pos.score),]$id
aggr$PMID29642015_rd4 = PMID29642015_rd4[order(PMID29642015_rd4$pos.score),]$id
aggr$PMID29642015_rd5 = PMID29642015_rd5[order(PMID29642015_rd5$pos.score),]$id

aggr$PMID34871331_A = PMID34871331_A[order(PMID34871331_A$pos.score),]$id
aggr$PMID34871331_B = PMID34871331_B[order(PMID34871331_B$pos.score),]$id

aggr$PMID35354039 = PMID35354039[order(PMID35354039$pos.score),]$id

aggr$PMID37652009_p1 = PMID37652009_p1[order(PMID37652009_p1$pos.score),]$id
aggr$PMID37652009_p2 = PMID37652009_p2[order(PMID37652009_p2$pos.score),]$id
aggr$PMID37652009_p3 = PMID37652009_p3[order(PMID37652009_p3$pos.score),]$id
aggr$PMID37652009_p4 = PMID37652009_p4[order(PMID37652009_p4$pos.score),]$id
aggr$PMID37652009_p5 = PMID37652009_p5[order(PMID37652009_p5$pos.score),]$id

aggr$PMID31919360_p = PMID31919360_p[order(PMID31919360_p$p.value),]$gene
aggr$PMID31919360_s = PMID31919360_s[order(PMID31919360_s$screen.p),]$gene

aggr$PMID36129973_proIAV = PMID36129973_proIAV[order(PMID36129973_proIAV$log2FoldChange),]$GeneSymbol

aggr$PMID32839537_IAVG = PMID32839537_IAVG[order(PMID32839537_IAVG$pos.score),]$id
aggr$PMID32839537_WSN = PMID32839537_WSN[order(PMID32839537_WSN$pos.score),]$id
aggr$PMID32839537_cs35m_gfp = PMID32839537_cs35m_gfp[order(PMID32839537_cs35m_gfp$pos.score),]$id

aggr$PMID34556855_H3N2_IFN = PMID34556855_H3N2_IFN[order(PMID34556855_H3N2_IFN$Average.Z_score.infection),]$Gene_Symbol
aggr$PMID34556855_H3N2_noIFN = PMID34556855_H3N2_noIFN[order(PMID34556855_H3N2_noIFN$Average.Z_score.infection),]$Gene_Symbol
aggr$PMID34556855_H5N1_IFN = PMID34556855_H5N1_IFN[order(PMID34556855_H5N1_IFN$Average.Z_score.infection),]$Gene_Symbol
aggr$PMID34556855_H5N1_noIFN = PMID34556855_H5N1_noIFN[order(PMID34556855_H5N1_noIFN$Average.Z_score.infection),]$Gene_Symbol

aggr$zrra_pro = zrra_pro[order(zrra_pro$Z_RSA, decreasing = F),]$Symbol

###RRA analysis
rra_proIAV = aggregateRanks(glist = aggr)
###Functional score calculation
rra_proIAV$n.log10.RRA = -log10(rra_proIAV$Score)
hist(rra_proIAV$n.log10.RRA, breaks = seq(0,26, by = 0.5))

write.csv(rra_proIAV, "./20240704/1.IAV_ranking_aggregation/RRA_proIAV_score.csv")








###Anti-IAV tables aggregation with RRA

###Build RRA object
aggr_2 = list()

###Add ranking table that ranked by gene functional importance
aggr_2$PMID28813663_A = PMID28813663_A[order(PMID28813663_A$pos.score),]$id
aggr_2$PMID28813663_B = PMID28813663_B[order(PMID28813663_B$pos.score),]$id
aggr_2$PMID28813663_C = PMID28813663_C[order(PMID28813663_C$pos.score),]$id

aggr_2$PMID36129973_antiIAV = PMID36129973_antiIAV[order(PMID36129973_antiIAV$log2FoldChange, decreasing = T),]$GeneSymbol

aggr_2$PMID32699298 = PMID32699298[order(PMID32699298$pos.score),]$id

aggr_2$zrra_anti = zrra_anti[order(zrra_anti$Z_RSA, decreasing = T),]$Symbol


###RRA analysis
rra_antiIAV = aggregateRanks(glist = aggr_2)
###Functional score calculation 
rra_antiIAV$n.log10.RRA = -log10(rra_antiIAV$Score)
hist(rra_antiIAV$n.log10.RRA, breaks = seq(0,26, by = 0.5))

rra_antiIAV$z_RRA = scale(rra_antiIAV$n.log10.RRA)
hist(rra_antiIAV$z_RRA)

write.csv(rra_antiIAV, "./20240704/1.IAV_ranking_aggregation/RRA_antiIAV_score.csv")



###Define the cumulative gene function (functional score) from pro and anti-IAV RRA score
rra_all = merge(rra_proIAV, rra_antiIAV, by = "Name", all = T)
head(rra_all)
colnames(rra_all)[c(2:7)] = c("proIAV_RRA_score", "proIAV_n.log10.RRA", "proIAV_z_RRA", "antiIAV_RRA_score",
                              "antiIAV_n.log10.RRA", "antiIAV_z_RRA")
for (i in 1:nrow(rra_all)) {
  if (is.na(rra_all$antiIAV_n.log10.RRA[i]) == T & is.na(rra_all$proIAV_n.log10.RRA[i]) == F) {
    rra_all$comprehensive_n.log10.RRA[i] = -rra_all$proIAV_n.log10.RRA[i]
  } else if (is.na(rra_all$antiIAV_n.log10.RRA[i]) == F & is.na(rra_all$proIAV_n.log10.RRA[i]) == T) {
    rra_all$comprehensive_n.log10.RRA[i] = rra_all$antiIAV_n.log10.RRA[i]
  } else (rra_all$comprehensive_n.log10.RRA[i] = rra_all$antiIAV_n.log10.RRA[i] - rra_all$proIAV_n.log10.RRA[i])
}
rra_all = rra_all[(order(rra_all$comprehensive_n.log10.RRA, decreasing = T)),]
rownames(rra_all) = 1:nrow(rra_all)

write.csv(rra_all, "./20240704/1.IAV_ranking_aggregation/RRA_score_all.csv")
