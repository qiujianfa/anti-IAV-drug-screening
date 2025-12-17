###Screening drugs from cmap database with functional score

rm(list = ls())
setwd("E:/Research_Project/Drug_screening/connectivity_map/level5/20240704/1.ranking/IAV")


###Read functional score table
setwd("E:/Research_Project/Drug_screening/genome_wide_screening/Flu/20240704/1.IAV_ranking_aggregation")
score = read.csv("RRA_score_all.csv", header = T, row.names = 1)


###Transform gene symbol to entrezID 
library(clusterProfiler)

bitr = bitr(score$Name, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Hs.eg.db")
bitr = bitr[!duplicated(bitr$SYMBOL),]
colnames(bitr)[1] = "Name"
score = merge(bitr, score, by = "Name", all = T)
score = score[order(score$comprehensive_n.log10.RRA, decreasing = F),]


###Read level5 z scaled cmap data (including GSE92742 and GSE70138 datasets)

library(cmapR)

###GSE92742 dataset
setwd("E:/Research_Project/Drug_screening/connectivity_map/level5")
level5 = parse_gctx(GSE92742_Broad_LINCS_Level5_COMPZ.MODZ_n473647x12328.gctx)
gc()
level5 = level5@mat
level5 = level5[which(rownames(level5) %in% intersect(score$ENTREZID, rownames(level5))), ]
score_2 = score[which(score$ENTREZID %in% intersect(score$ENTREZID, rownames(level5))), ]
score_2 = score_2[match(score_2$ENTREZID, rownames(level5)), ]

###Calculate the drug score
drug_score = colSums(level5 * score_2$comprehensive_n.log10.RRA)
drug_score = data.frame(drugID = colnames(level5), drug_score = drug_score)
setwd("E:/Research_Project/Drug_screening/connectivity_map/level5/20240704/1.ranking/IAV")
write.csv(drug_score, "drug_score_GSE92742.csv")
gc()




###Match the experimental well ID with drug information
setwd("E:/Research_Project/Drug_screening/connectivity_map/level5")
info = read.table("GSE92742_Broad_LINCS_sig_info.txt", sep = "\t", header = T, quote = "")

colnames(drug_score)[1] = "sig_id"
drug_score = merge(info, drug_score, by = "sig_id")
drug_score = drug_score[order(drug_score$drug_score, decreasing = T),]

summary(as.factor(drug_score$cell_id))
summary(as.factor(drug_score$pert_type))

setwd("E:/Research_Project/Drug_screening/connectivity_map/level5/20240704/1.ranking/IAV")
write.csv(drug_score, "pert_score_list.csv")

###Small moleclue compounds treatment data only
drug_comp = drug_score[which(drug_score$pert_type %in% c("ctl_vehicle", "trt_cp")),]
write.csv(drug_comp, "compound_drug_score_GSE92742.csv")

###A549 cell line data only
drug_comp_A549 = drug_comp[which(drug_comp$cell_id == "A549"),]
write.csv(drug_comp_A549, "A549_compound_drug_score_GSE92742.csv")


###Average drug score of the same compounds treatment
summary(factor(drug_comp_A549$pert_itime))
library(dplyr)
drug_ord = data.frame(score = drug_comp_A549$drug_score,
                      drug_name = drug_comp_A549$pert_iname)
drug_ord = drug_ord %>% group_by(drug_name) %>% summarise(across(.cols = everything(), mean))

count = table(drug_comp_A549$pert_iname)
count = data.frame(count)
colnames(count)[1] = "drug_name"
drug_ord = merge(drug_ord, count, by = "drug_name")

drug_ord = drug_ord[order(drug_ord$score, decreasing = T),]

write.csv(drug_ord, "drug_score_GSE92742.csv", row.names = F)








###GSE70138 dataset
setwd("E:/Research_Project/Drug_screening/connectivity_map/level5")
drug = parse_gctx("GSE70138_Broad_LINCS_Level5_COMPZ_n118050x12328.gctx")
drug = drug@mat
gc()
drug = drug[which(rownames(drug) %in% intersect(score$ENTREZID, rownames(drug))), ]
score_2 = score[which(score$ENTREZID %in% intersect(score$ENTREZID, rownames(drug))), ]
score_2 = score_2[match(score_2$ENTREZID, rownames(drug)), ]

###Calculate the drug score
drug_score = colSums(drug * score_2$comprehensive_n.log10.RRA)
drug_score = data.frame(drugID = colnames(drug), drug_score = drug_score)
drug_score = drug_score[order(drug_score$drug_score, decreasing = T),]
setwd("E:/Research_Project/Drug_screening/connectivity_map/level5/20240704/1.ranking/IAV")
write.csv(drug_score, "S_drug_score.csv")
gc()


###Match the experimental well ID with drug information
setwd("E:/Research_Project/Drug_screening/connectivity_map/level5")
info = read.table("GSE70138_Broad_LINCS_sig_info.txt", sep = "\t", header = T, quote = "")

drug = drug_score
colnames(drug)[1] = "sig_id"
drug = merge(info, drug, by = "sig_id")
drug = drug[order(drug$drug_score, decreasing = T),]

summary(as.factor(drug$cell_id))
summary(as.factor(drug$pert_type))

setwd("E:/Research_Project/Drug_screening/connectivity_map/level5/20240704/1.ranking/IAV")
write.csv(drug, "S_pert_score_list.csv")

###Small moleclue compounds treatment data only
drug_comp = drug_score[which(drug_score$pert_type %in% c("ctl_vehicle", "trt_cp")),]
write.csv(drug_comp, "compound_drug_score_GSE70138.csv")

###A549 cell line data only
drug_comp_A549 = drug_comp[which(drug_comp$cell_id == "A549"),]
write.csv(drug_comp_A549, "A549_compound_drug_score_GSE70138.csv")



###Average drug score of the same compounds treatment
summary(factor(drug_comp_A549$pert_itime))
library(dplyr)
drug_ord = data.frame(score = drug_comp_A549$drug_score,
                      drug_name = drug_comp_A549$pert_iname)
drug_ord = drug_ord %>% group_by(drug_name) %>% summarise(across(.cols = everything(), mean))

count = table(drug_comp_A549$pert_iname)
count = data.frame(count)
colnames(count)[1] = "drug_name"
drug_ord = merge(drug_ord, count, by = "drug_name")

drug_ord = drug_ord[order(drug_ord$score, decreasing = T),]

write.csv(drug_ord, "drug_score_GSE70138.csv", row.names = F)


