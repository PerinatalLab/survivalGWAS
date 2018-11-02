#!/usr/bin/env Rscript

options(stringsAsFactors = F)

library(dplyr)
library(tidyr)

pheno_path= '/mnt/hdd/data/mobaqs/p1724/'
pc_path= '/mnt/hdd/data/geno/harvest-aux/'
kin_path= '/mnt/hdd/common/harvest/relatedness/'
moms_kin= 'moms.kin0'
fets_kin= 'fets.kin0'
outpath= '/mnt/hdd/common/harvest/pheno/'
file_pref= 'harvest_PROM_surv'
final_vars= c('SentrixID_1', 'SVLEN_UL_DG', 'PROM', 'PC1', 'PC2', 'PC3', 'PC4', 'PC5', 'PC6', 'PARITY0','BATCH')

SelectRelated= function(kin_path, sample_list){
  kin= read.table(kin_path, h=T, comment.char = "", sep= '\t')
 kin= kin %>% filter(KINSHIP>0.044)
 kin= kin %>% filter(X.FID1 %in% sample_list & FID2 %in% sample_list)
 kin= kin %>% mutate(ID1= paste(X.FID1,ID1, sep= ":"),
                      ID2= paste(FID2, ID2, sep= ":")) %>% select(ID1, ID2, KINSHIP)
  kin_temp= kin
  colnames(kin_temp)= c("ID2","ID1","KINSHIP")
  kin_temp= rbind(kin_temp, kin)  
  kin_temp= kin_temp %>% add_count(ID1)
  kin_temp= kin_temp %>% add_count(ID2)
  kin_temp= arrange(kin_temp, n, nn)
  to_keep= list()
  
  for (i in 1:nrow(kin_temp)) {
    if (kin_temp[i,"ID1"] %in% unlist(kin_temp[0:i,"ID2"])) {
      kin_temp[i,"ID2"]= "X"
    }
    else
      to_keep[[i]] <- kin_temp[["ID1"]][i]
  }
  to_remove= kin_temp %>% filter(!(ID1 %in% unlist(to_keep))) %>% select(ID1) 
  to_remove= to_remove[!duplicated(to_remove$ID1),]
  to_remove= to_remove %>% separate(ID1, c('FID','ID'), sep=":")

  return(unlist(to_remove[,1]))
}


mfr= read.table(paste0(pheno_path, 'harvest_mfr.csv'), h=T, sep= ';')
pc_moms= read.table(paste0(pc_path, 'plink_covar_mothers'), h=F, sep=' ')
pc_fets= read.table(paste0(pc_path, 'plink_covar_offspring'), h=F, sep= ' ')
flags = read.table(paste0(pc_path, "harvest-flag-list.txt"), h=T)
link = read.table(paste0(pheno_path, 'harvest_linkage.csv'), sep=";", h=T)

names(pc_moms)= c('x','SentrixID_1','PC1','PC2','PC3','PC4','PC5','PC6','PC7','PC8','PC9','PC10')
names(pc_fets)= c('x','SentrixID_1','PC1','PC2','PC3','PC4','PC5','PC6','PC7','PC8','PC9','PC10')

flags= flags %>% filter(genotypesOK=='TRUE', coreOK== 'TRUE') %>% mutate(BATCH = as.numeric(BATCH=="M24")) %>% select(IID, BATCH)

link= inner_join(link, flags[,c('IID','BATCH')], by= c('SentrixID_1'='IID'))
pc_moms= inner_join(pc_moms, link, by= 'SentrixID_1')
pc_fets= inner_join(pc_fets, link, by= 'SentrixID_1')

final= mutate(mfr, PROM = as.numeric(!is.na(VANNAVGANG)), PARITY0= as.numeric(PARITET_5==0))
final= filter(final, FLERFODSEL==0 , DODKAT<6 | DODKAT>10, !is.na(SVLEN_UL_DG), SVLEN_UL_DG<308 & SVLEN_UL_DG>154 & !is.na(PROM) & (VANNAVGANG!=3 | is.na(VANNAVGANG)))

final= final[order(final$PROM, decreasing= T),]

final = filter(final, is.na(IVF),ABRUPTIOP==0,
                           PLACENTA_PREVIA==0,
			FOSTERV_POLYHYDRAMNION==0,
			C00_MALF_ALL==0)

final_sens = filter(final,
                           is.na(PREEKL),
                           is.na(DIABETES_MELLITUS),
                           HYPERTENSJON_ALENE==0 & HYPERTENSJON_KRONISK==0)

final_moms= inner_join(final, pc_moms, by= 'PREG_ID_1724') %>% select(final_vars)
final_fets= inner_join(final, pc_fets, by= 'PREG_ID_1724') %>% select(final_vars)

final_moms= final_moms[!duplicated(final_moms$SentrixID_1),]
final_fets= final_fets[!duplicated(final_fets$SentrixID_1),]

final_sens_moms= inner_join(final_sens, pc_moms, by= 'PREG_ID_1724') %>% select(final_vars)
final_sens_fets= inner_join(final_sens, pc_fets, by= 'PREG_ID_1724') %>% select(final_vars)

final_sens_moms= final_sens_moms[!duplicated(final_sens_moms$SentrixID_1),]
final_sens_fets= final_sens_fets[!duplicated(final_sens_fets$SentrixID_1),]

final_moms= final_moms %>% filter(!(SentrixID_1 %in% SelectRelated(paste0(kin_path, moms_kin), final_moms$SentrixID_1)))

final_sens_moms= final_sens_moms %>% filter(!(SentrixID_1 %in% SelectRelated(paste0(kin_path, moms_kin), final_sens_moms$SentrixID_1)))

final_fets= final_fets %>% filter(!(SentrixID_1 %in% SelectRelated(paste0(kin_path, fets_kin), final_fets$SentrixID_1)))

final_sens_fets= final_sens_fets %>% filter(!(SentrixID_1 %in% SelectRelated(paste0(kin_path, fets_kin), final_sens_fets$SentrixID_1)))


write.table(final_moms, paste0(outpath, file_pref, '_moms'), row.names=F, col.names=T, sep= '\t', quote=F)
write.table(final_sens_moms, paste0(outpath, file_pref, '_moms_sens'), row.names=F, col.names=T, sep= '\t', quote=F)
write.table(final_fets, paste0(outpath, file_pref, '_fets'), row.names=F, col.names=T, sep= '\t', quote=F)
write.table(final_sens_fets, paste0(outpath, file_pref, '_fets_sens'), row.names=F, col.names=T, sep= '\t', quote=F)

