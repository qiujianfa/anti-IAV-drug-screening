rm(list = ls())
setwd("E:/Research_Project/Drug_screening/genome_wide_screening/Flu")


###pro-IAV factors screening

###CRISPR_KO
###PMID29642015
PMID29642015_rd1 = read.table("./PMID29642015/GSE111166/mageck/Rd1.gene_summary.txt", header = T, sep = "\t")
PMID29642015_rd2 = read.table("./PMID29642015/GSE111166/mageck/Rd2.gene_summary.txt", header = T, sep = "\t")
PMID29642015_rd3 = read.table("./PMID29642015/GSE111166/mageck/Rd3.gene_summary.txt", header = T, sep = "\t")
PMID29642015_rd4 = read.table("./PMID29642015/GSE111166/mageck/Rd4.gene_summary.txt", header = T, sep = "\t")
PMID29642015_rd5 = read.table("./PMID29642015/GSE111166/mageck/Rd5.gene_summary.txt", header = T, sep = "\t")

###PMID34871331
PMID34871331_A = read.table("./PMID34871331/PRJNA777037/2.mageck/mageck_A.gene_summary.txt", header = T, sep = "\t")
PMID34871331_B = read.table("./PMID34871331/PRJNA777037/2.mageck/mageck_B.gene_summary.txt", header = T, sep = "\t")

###PMID35354039
PMID35354039 = read.csv("./PMID35354039/MAGECK_full_table.csv", header = T)

###PMID37652009
PMID37652009_p1 = read.table("./PMID37652009/mageck/output/mageck_p1.gene_summary.txt", header = T, sep = "\t")
PMID37652009_p2 = read.table("./PMID37652009/mageck/output/mageck_p2.gene_summary.txt", header = T, sep = "\t")
PMID37652009_p3 = read.table("./PMID37652009/mageck/output/mageck_p3.gene_summary.txt", header = T, sep = "\t")
PMID37652009_p4 = read.table("./PMID37652009/mageck/output/mageck_p4.gene_summary.txt", header = T, sep = "\t")
PMID37652009_p5 = read.table("./PMID37652009/mageck/output/mageck_p5.gene_summary.txt", header = T, sep = "\t")

###PMID32839537
PMID32839537_IAVG = read.csv("./PMID32839537_ko/IAV-G_2_IAV-G.csv", header = T)
PMID32839537_WSN = read.csv("./PMID32839537_ko/IAV-G_2_WSN.csv", header = T)
PMID32839537_cs35m_gfp = read.csv("./PMID32839537_ko/SC35M Flu-GFP.csv", header = T)

###RNAi
PMID34556855_H3N2_IFN = read.csv("./PMID34556855_RNAi/H3N2 +IFN.csv", header = T)
PMID34556855_H3N2_noIFN = read.csv("./PMID34556855_RNAi/H3N2 -IFN.csv", header = T)
PMID34556855_H5N1_IFN = read.csv("./PMID34556855_RNAi/H5N1 +IFN.csv", header = T)
PMID34556855_H5N1_noIFN = read.csv("./PMID34556855_RNAi/H5N1 -IFN.csv", header = T)


###FACS screening
###PMID31919360
PMID31919360_p = read.csv("./PMID31919360/primary.csv", header = T)
PMID31919360_s = read.csv("./PMID31919360/secondary.csv", header = T)









###anti-IAV screening

###CRISPRa PMID28813663
PMID28813663_A = read.table("./PMID28813663/GSE89371/2.mageck/1.output/GSE89371_screenA.gene_summary.txt", header = T, sep = "\t")
PMID28813663_B = read.table("./PMID28813663/GSE89371/2.mageck/1.output/GSE89371_screenB.gene_summary.txt", header = T, sep = "\t")
PMID28813663_C = read.table("./PMID28813663/GSE89371/2.mageck/1.output/GSE89371_screenC.gene_summary.txt", header = T, sep = "\t")

###FACSPMID32699298
PMID32699298 = read.csv("./PMID32699298_scirep_anitIAV/mageck.csv", header = T)







###pro and anti screening

###FACS PMID36129973
PMID36129973 = read.csv("./PMID36129973/list.csv", header = T)
PMID36129973_proIAV = PMID36129973[which(PMID36129973$log2FoldChange < 0),]
PMID36129973_antiIAV = PMID36129973[which(PMID36129973$log2FoldChange > 0),]

###RNAi
###PMID26651948
zrra = read.csv("./PMID26651948_MetaAnalysis_4Sets/z_score.csv", header = T)
zrra_pro = zrra[which(zrra$Z_RSA < 0),]
zrra_anti = zrra[which(zrra$Z_RSA > 0),]







###RRA for pro-IAV screening

library(RobustRankAggreg)

aggr = list()

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



###leave-one-out iteration to generate new RRA ranking

rra_proIAV = list()
for (i in 1:length(aggr)) {
  glist = aggr[names(aggr) != names(aggr)[i]]
  rra_proIAV[[i]] = aggregateRanks(glist = glist)
  rra_proIAV[[i]]$proIAV_n.log10.RRA = -log10(rra_proIAV[[i]]$Score)
  print(i)
}
length(rra_proIAV)
rra_proIAV$all = aggregateRanks(glist = aggr)
rra_proIAV$all$proIAV_n.log10.RRA = -log10(rra_proIAV$all$Score)








###RRA for anti-IAV screening

aggr_2 = list()

aggr_2$PMID28813663_A = PMID28813663_A[order(PMID28813663_A$pos.score),]$id
aggr_2$PMID28813663_B = PMID28813663_B[order(PMID28813663_B$pos.score),]$id
aggr_2$PMID28813663_C = PMID28813663_C[order(PMID28813663_C$pos.score),]$id

aggr_2$PMID36129973_antiIAV = PMID36129973_antiIAV[order(PMID36129973_antiIAV$log2FoldChange, decreasing = T),]$GeneSymbol
aggr_2$PMID32699298 = PMID32699298[order(PMID32699298$pos.score),]$id
aggr_2$zrra_anti = zrra_anti[order(zrra_anti$Z_RSA, decreasing = T),]$Symbol


###leave-one-out iteration to generate new RRA ranking

rra_antiIAV = list()
for (i in 1:length(aggr_2)) {
  glist = aggr_2[names(aggr_2) != names(aggr_2)[i]]
  rra_antiIAV[[i]] = aggregateRanks(glist = glist)
  rra_antiIAV[[i]]$antiIAV_n.log10.RRA = -log10(rra_antiIAV[[i]]$Score)
  print(i)
}
rra_antiIAV$all = aggregateRanks(glist = aggr_2)
rra_antiIAV$all$antiIAV_n.log10.RRA = -log10(rra_antiIAV$all$Score)








###generate a merged table with pro and anti RRA score

rra_all = list()
for (i in 1:(length(rra_proIAV)-1)) {
  pro = rra_proIAV[names(rra_proIAV) != names(rra_proIAV)[length(rra_proIAV)]]
  anti = as.data.frame(rra_antiIAV[names(rra_antiIAV) == names(rra_antiIAV)[length(rra_antiIAV)]])
  colnames(anti) = c("Name", "Score", "antiIAV_n.log10.RRA")
  rra_all[[i]] = merge(pro[[i]], anti, by = "Name", all = T)
  colnames(rra_all[[i]]) = c("Name", "pro_Score", "proIAV_n.log10.RRA", "anti_Score", "antiIAV_n.log10.RRA")
}

for (i in 1:(length(rra_antiIAV)-1)) {
  anti = rra_antiIAV[names(rra_antiIAV) != names(rra_antiIAV)[length(rra_antiIAV)]]
  pro = as.data.frame(rra_proIAV[names(rra_proIAV) == names(rra_proIAV)[length(rra_proIAV)]])
  colnames(pro) = c("Name", "Score", "proIAV_n.log10.RRA")
  rra_all[[i+24]] = merge(pro, anti[[i]], by = "Name", all = T)
  colnames(rra_all[[i]]) = c("Name", "pro_Score", "proIAV_n.log10.RRA", "anti_Score", "antiIAV_n.log10.RRA")
}


###calculate the anti-IAV score (comprehensive RRA score) of each gene

for (i in 1:length(rra_all)) {
  for (j in 1:nrow(rra_all[[i]])) {
    if (is.na(rra_all[[i]]$antiIAV_n.log10.RRA[j]) == T & is.na(rra_all[[i]]$proIAV_n.log10.RRA[j]) == F) {
      rra_all[[i]]$comprehensive_n.log10.RRA[j] = -rra_all[[i]]$proIAV_n.log10.RRA[j]
    } else if (is.na(rra_all[[i]]$antiIAV_n.log10.RRA[j]) == F & is.na(rra_all[[i]]$proIAV_n.log10.RRA[j]) == T) {
      rra_all[[i]]$comprehensive_n.log10.RRA[j] = rra_all[[i]]$antiIAV_n.log10.RRA[j]
    } else (rra_all[[i]]$comprehensive_n.log10.RRA[j] = rra_all[[i]]$antiIAV_n.log10.RRA[j] - rra_all[[i]]$proIAV_n.log10.RRA[j])
  }
}


###rank genes by anti-IAV score

for (i in 1:length(rra_all)) {
  rra_all[[i]] = rra_all[[i]][order(rra_all[[i]]$comprehensive_n.log10.RRA, decreasing = T),]
  rra_all[[i]]$Rank = 1:nrow(rra_all[[i]])
}

no_delet = no_delet[order(no_delet$comprehensive_n.log10.RRA, decreasing = T),]
no_delet$Rank = 1:nrow(no_delet)


###read original ranking table

no_delet = read.csv("E:/Research_Project/Drug_screening/genome_wide_screening/Flu/20240704/1.IAV_ranking_aggregation/RRA_score_all.csv", header = T, row.names = 1)


###Calculate the sum of rank differences between leave-one-out and original rank

x = list()
y = c()
for (i in 1:length(rra_all)) {
  x[[i]] = merge(rra_all[[i]], no_delet, by = "Name")
  x[[i]] = abs(x[[i]]$Rank.x - x[[i]]$Rank.y)
  y[i] = sum(x[[i]])
  print(i)
}

names(x) = c(names(aggr[1:24]), names(aggr_2[1:6]))
names(y) = c(names(aggr[1:24]), names(aggr_2[1:6]))
y
plot(y)












###randomly generate 30 ranking table

nrow(no_delet)

random = list()
for (i in 1:30) {
  set.seed(1000+i)
  random[[i]] = sample(no_delet$Name, nrow(no_delet), replace = F)
}

###add one randomly ranking data into a leave-one-out dataset in each iteration for recalculate an perturbed RRA ranking

###for pro-IAV screening datasets

rra_pro_random = list()
for (i in 1:length(aggr)) {
  glist = aggr
  glist$random = random[[i]]
  rra_pro_random[[i]] = aggregateRanks(glist = glist)
  rra_pro_random[[i]]$antiIAV_n.log10.RRA = -log10(rra_pro_random[[i]]$Score)
  print(i)
}

###for anti-IAV screening datasets

rra_anti_random = list()
for (i in 1:length(aggr_2)) {
  glist = aggr_2
  glist$random = random[[i]]
  rra_anti_random[[i]] = aggregateRanks(glist = glist)
  rra_anti_random[[i]]$antiIAV_n.log10.RRA = -log10(rra_anti_random[[i]]$Score)
  print(i)
}

###generate a merged table with pro and anti RRA score

rra_all_random = list()
for (i in 1:length(rra_pro_random)) {
  pro = rra_pro_random[[i]]
  anti = as.data.frame(rra_antiIAV[names(rra_antiIAV) == names(rra_antiIAV)[length(rra_antiIAV)]])
  colnames(anti) = c("Name", "Score", "antiIAV_n.log10.RRA")
  rra_all_random[[i]] = merge(pro, anti, by = "Name", all = T)
  colnames(rra_all_random[[i]]) = c("Name", "pro_Score", "proIAV_n.log10.RRA", "anti_Score", "antiIAV_n.log10.RRA")
}

for (i in 1:length(rra_anti_random)) {
  anti = rra_anti_random[[i]]
  pro = as.data.frame(rra_proIAV[names(rra_proIAV) == names(rra_proIAV)[length(rra_proIAV)]])
  colnames(pro) = c("Name", "Score", "proIAV_n.log10.RRA")
  rra_all_random[[i+24]] = merge(pro, anti, by = "Name", all = T)
  colnames(rra_all_random[[i]]) = c("Name", "pro_Score", "proIAV_n.log10.RRA", "anti_Score", "antiIAV_n.log10.RRA")
}



###calculate the anti-IAV score (comprehensive RRA score) of each gene

for (i in 1:length(rra_all_random)) {
  rra_all_random[[i]]$comprehensive_n.log10.RRA = rra_all_random[[i]]$antiIAV_n.log10.RRA - rra_all_random[[i]]$proIAV_n.log10.RRA
  rra_all_random[[i]] = rra_all_random[[i]][order(rra_all_random[[i]]$comprehensive_n.log10.RRA, decreasing = T),]
  rra_all_random[[i]]$Rank = 1:nrow(rra_all_random[[i]])
}




###Calculate the sum of rank differences between random group and original rank

z = list()
w = c()
for (i in 1:length(rra_all_random)) {
  z[[i]] = merge(rra_all_random[[i]], no_delet, by = "Name")
  z[[i]] = abs(z[[i]]$Rank.x - z[[i]]$Rank.y)
  w[i] = sum(z[[i]])
  print(i)
}

names(z) = c(names(aggr[1:24]), names(aggr_2[1:6]))
names(w) = c(names(aggr[1:24]), names(aggr_2[1:6]))
w
plot(w)






####plot

setwd("E:/Research_Project/Drug_combination/genome_wide_screening/Flu/20240704/5.Bias_test")
library(ggplot2)
###by samples
y = data.frame(rank = y, group = rep("delet", length(y)))
w = data.frame(rank = w, group = rep("random", length(w)))
by_group = rbind(y,w)
write.csv(by_group, "delet_vs_random.csv")


by_group_pro = rbind(y[1:24,], w[1:24,])
by_group_anti = rbind(y[25:30,], w[25:30,])


ggplot(by_group, aes(x = group, y = rank)) +
  geom_boxplot(outliers = F) +
  geom_jitter(width = 0.35, height = 0, size = 2, alpha = 0.6, color = "tomato") +
  theme_classic()
  
ggplot(by_group_pro, aes(x = group, y = rank)) +
  geom_boxplot(outliers = F) +
  geom_jitter(width = 0.35, height = 0, size = 2, alpha = 0.6, color = "tomato") +
  theme_classic()

ggplot(by_group_anti, aes(x = group, y = rank)) +
  geom_boxplot(outliers = F) +
  geom_jitter(width = 0.35, height = 0, size = 2, alpha = 0.6, color = "tomato") +
  theme_classic()


