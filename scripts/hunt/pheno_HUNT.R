#!/usr/bin/env Rscript

options(stringsAsFactors = F)

library(dplyr)
library(tidyr)

pheno_path= '/mnt/archive/hunt/phenotypes/mfr/'
pc_path= '/mnt/work/hunt/pca/'
kin_path= '/mnt/work/hunt/relatedness/'
moms_kin= 'mother_samples_related.kin0'
fets_kin= 'fetal_samples_related.kin0'
outpath= '/mnt/work/hunt/pheno/'
file_pref_PROM= 'HUNT_PROM_surv'
file_pref_spont= 'HUNT_spont_surv'
final_vars_PROM_MOR= c('MOR_PID', 'SVLEN_DG', 'PROM', 'PC1', 'PC2', 'PC3', 'PC4', 'PC5', 'PC6', 'FAAR', 'PARITY0')
final_vars_PROM_BARN= c('BARN_PID', 'SVLEN_DG', 'PROM', 'PC1', 'PC2', 'PC3', 'PC4', 'PC5', 'PC6', 'FAAR', 'PARITY0')
final_vars_spont_MOR= c('MOR_PID', 'SVLEN_DG', 'spont', 'PC1', 'PC2', 'PC3', 'PC4', 'PC5', 'PC6', 'FAAR', 'PARITY0')
final_vars_spont_BARN= c('BARN_PID', 'SVLEN_DG', 'spont', 'PC1', 'PC2', 'PC3', 'PC4', 'PC5', 'PC6', 'FAAR', 'PARITY0')
icd_path= '/mnt/archive/hunt/phenotypes/hnt/Gestational age_AP_Avidentifisert.txt'

mfr= read.table(paste0(pheno_path, 'MFR.txt'), h=T, sep= '\t')
pc_moms= read.table(paste0(pc_path, 'mother_pca.sscore'), h=F, sep='\t')
pc_fets= read.table(paste0(pc_path, 'fetal_pca.sscore'), h=F, sep= '\t')

names(pc_moms)= c('MOR_PID', 'IID','X','X1','PC1','PC2','PC3','PC4','PC5','PC6','PC7','PC8','PC9','PC10')
names(pc_fets)= c('BARN_PID', 'IID','X','X1','PC1','PC2','PC3','PC4','PC5','PC6','PC7','PC8','PC9','PC10')

### Spontaneous delivery

SelectRelated= function(kin_path, sample_list, df, var){
  kin= read.table(kin_path, h=T, comment.char = "", sep= '\t')
 kin= kin %>% filter(KINSHIP>0.0884)
 kin= kin %>% filter(X.FID1 %in% sample_list & FID2 %in% sample_list)
 kin= kin %>% mutate(ID1= paste(X.FID1,ID1, sep= ":"),
                      ID2= paste(FID2, ID2, sep= ":")) %>% select(ID1, ID2, KINSHIP)
  kin_temp= kin
  colnames(kin_temp)= c("ID2","ID1","KINSHIP")
  kin_temp= rbind(kin_temp, kin)
  kin_temp= kin_temp %>% add_count(ID1)
  kin_temp= kin_temp %>% add_count(ID2)
df[['ID1']]= paste(df[[var]], df[[var]], sep=':')
kin_temp= inner_join(kin_temp, df[,c('ID1','spont')], by= 'ID1')
  kin_temp= arrange(kin_temp, desc(spont), n, nn)
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

final= mutate(mfr, PARITY0= as.numeric(PARITET_MFR==0), 
		spont= as.numeric(FSTART==1 & (is.na(KSNITT) | KSNITT>1) &
		(is.na(KSNITT_PLANLAGT) | KSNITT_PLANLAGT==1) &
		is.na(INDUKSJON_PROSTAGLANDIN) & is.na(INDUKSJON_ANNET) &
		is.na(INDUKSJON_OXYTOCIN) & is.na(INDUKSJON_AMNIOTOMI)))

final= filter(final, is.na(FLERFODSEL), DODKAT<6 | DODKAT>10, !is.na(SVLEN_DG), SVLEN_DG<308 & SVLEN_DG>154)

final = filter(final, is.na(ART),is.na(ABRUPTIOP),
                           is.na(PLACENTA_PREVIA),
			is.na(FOSTERV_POLYHYDRAMNION),
			is.na(MISD))

final= final[order(final$spont, decreasing= T),]

final_sens= filter(final,
			   is.na(PREEKL),
			   is.na(DIABETES_MELLITUS),
			   is.na(HYPERTENSJON_ALENE) & is.na(HYPERTENSJON_KRONISK))

final_moms= inner_join(final, pc_moms, by= 'MOR_PID') %>% select(final_vars_spont_MOR)
final_fets= inner_join(final, pc_fets, by= 'BARN_PID') %>% select(final_vars_spont_BARN)

final_moms= final_moms[!duplicated(final_moms$MOR_PID),]
final_fets= final_fets[!duplicated(final_fets$BARN_PID),]

final_sens_moms= inner_join(final_sens, pc_moms, by= 'MOR_PID') %>% select(final_vars_spont_MOR)
final_sens_fets= inner_join(final_sens, pc_fets, by= 'BARN_PID') %>% select(final_vars_spont_BARN)

final_sens_moms= final_sens_moms[!duplicated(final_sens_moms$MOR_PID),]
final_sens_fets= final_sens_fets[!duplicated(final_sens_fets$BARN_PID),]

moms= final_moms %>% filter(!(MOR_PID %in% SelectRelated(paste0(kin_path, moms_kin), final_moms$MOR_PID,  final_moms, 'MOR_PID')))
moms_sens= final_sens_moms %>% filter(!(MOR_PID %in% SelectRelated(paste0(kin_path, moms_kin), final_sens_moms$MOR_PID,  final_sens_moms, 'MOR_PID')))

fets= final_fets %>% filter(!(BARN_PID %in% SelectRelated(paste0(kin_path, fets_kin), final_fets$BARN_PID, final_fets,'BARN_PID' )))
fets_sens= final_sens_fets %>% filter(!(BARN_PID %in% SelectRelated(paste0(kin_path, fets_kin), final_sens_fets$BARN_PID, final_sens_fets, 'BARN_PID')))

write.table(moms, paste0(outpath, file_pref_spont, '_moms'), row.names=F, col.names=T, sep= '\t', quote=F)
write.table(moms_sens, paste0(outpath, file_pref_spont, '_moms_sens'), row.names=F, col.names=T, sep= '\t', quote=F)
write.table(fets, paste0(outpath, file_pref_spont, '_fets'), row.names=F, col.names=T, sep= '\t', quote=F)
write.table(fets_sens, paste0(outpath, file_pref_spont, '_fets_sens'), row.names=F, col.names=T, sep= '\t', quote=F)

#### New pheno for PROM. Accorded with Bo and Nils on the 19/01/02
# PROM diagnosis from MFR or ICD codes 9 or 10.
# Following ICD codes:
# ICD9:  658.1 and 658.2
# ICD10: O42

SelectRelated= function(kin_path, sample_list, df, var){
  kin= read.table(kin_path, h=T, comment.char = "", sep= '\t')
 kin= kin %>% filter(KINSHIP>0.0884)
 kin= kin %>% filter(X.FID1 %in% sample_list & FID2 %in% sample_list)
 kin= kin %>% mutate(ID1= paste(X.FID1,ID1, sep= ":"),
                      ID2= paste(FID2, ID2, sep= ":")) %>% select(ID1, ID2, KINSHIP)
  kin_temp= kin
  colnames(kin_temp)= c("ID2","ID1","KINSHIP")
  kin_temp= rbind(kin_temp, kin)  
  kin_temp= kin_temp %>% add_count(ID1)
  kin_temp= kin_temp %>% add_count(ID2)
df[['ID1']]= paste(df[[var]], df[[var]], sep=':')
kin_temp= inner_join(kin_temp, df[,c('ID1','PROM')], by= 'ID1')
  kin_temp= arrange(kin_temp, desc(PROM), n, nn)
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

icd= read.table(icd_path, h=T, sep= '\t', quote='')
icd$FAAR= as.numeric(substring(icd$Diagnosedato, 1, 4))

icd= mutate(icd, 
		PROM_icd9= grepl('658.1|658.2', icd$ICD_9, fixed=F),
		PROM_icd10= grepl('O42', icd$ICD_10, fixed=F))
icd= icd %>% mutate(PROM_icd= ifelse(PROM_icd9== 1 | PROM_icd10== 1, 1, 0)) %>% select(PID.106764, FAAR, PROM_icd)
icd= icd %>% filter(PROM_icd==1)

final= mutate(mfr, PROM = as.numeric(!is.na(VANNAVGANG)), PARITY0= as.numeric(PARITET_MFR==0))

final= left_join(final, icd, by= c('MOR_PID'='PID.106764', 'FAAR'))

final= mutate(final, PROM= ifelse((PROM==0) & (is.na(PROM_icd)), 0, 1))

final= filter(final, is.na(FLERFODSEL), DODKAT<6 | DODKAT>10, !is.na(SVLEN_DG), SVLEN_DG<308 & SVLEN_DG>154 & !is.na(PROM))

final = filter(final, is.na(ART),is.na(ABRUPTIOP),
                           is.na(PLACENTA_PREVIA),
                        is.na(FOSTERV_POLYHYDRAMNION),
                        is.na(MISD))

final= final[order(final$PROM, decreasing= T),]

# PROM sensitivity analysis will be PROM with spontaneous delivery

final_sens= mutate(final, 
		spont= as.numeric(FSTART==1 & (is.na(KSNITT) | KSNITT>1) &
                (is.na(KSNITT_PLANLAGT) | KSNITT_PLANLAGT==1) &
                is.na(INDUKSJON_PROSTAGLANDIN) & is.na(INDUKSJON_ANNET) &
                is.na(INDUKSJON_OXYTOCIN) & is.na(INDUKSJON_AMNIOTOMI)))

final_sens= mutate(final_sens, PROM= ifelse(spont== 1 & PROM==1, 1, 0))
final_sens= final_sens[order(final_sens$PROM, decreasing= T),]

final_moms= inner_join(final, pc_moms, by= 'MOR_PID') %>% select(final_vars_PROM_MOR)
final_fets= inner_join(final, pc_fets, by= 'BARN_PID') %>% select(final_vars_PROM_BARN)

final_moms= final_moms[!duplicated(final_moms$MOR_PID),]
final_fets= final_fets[!duplicated(final_fets$BARN_PID),]

final_sens_moms= inner_join(final_sens, pc_moms, by= 'MOR_PID') %>% select(final_vars_PROM_MOR)
final_sens_fets= inner_join(final_sens, pc_fets, by= 'BARN_PID') %>% select(final_vars_PROM_BARN)

final_sens_moms= final_sens_moms[!duplicated(final_sens_moms$MOR_PID),]
final_sens_fets= final_sens_fets[!duplicated(final_sens_fets$BARN_PID),]

moms= final_moms %>% filter(!(MOR_PID %in% SelectRelated(paste0(kin_path, moms_kin), final_moms$MOR_PID,  final_moms, 'MOR_PID')))
moms_sens= final_sens_moms %>% filter(!(MOR_PID %in% SelectRelated(paste0(kin_path, moms_kin), final_sens_moms$MOR_PID,  final_sens_moms, 'MOR_PID')))

fets= final_fets %>% filter(!(BARN_PID %in% SelectRelated(paste0(kin_path, fets_kin), final_fets$BARN_PID, final_fets,'BARN_PID' )))
fets_sens= final_sens_fets %>% filter(!(BARN_PID %in% SelectRelated(paste0(kin_path, fets_kin), final_sens_fets$BARN_PID, final_sens_fets, 'BARN_PID')))

write.table(moms, paste0(outpath, file_pref_PROM, '_moms'), row.names=F, col.names=T, sep= '\t', quote=F)
write.table(moms_sens, paste0(outpath, file_pref_PROM, '_moms_sens'), row.names=F, col.names=T, sep= '\t', quote=F)
write.table(fets, paste0(outpath, file_pref_PROM, '_fets'), row.names=F, col.names=T, sep= '\t', quote=F)
write.table(fets_sens, paste0(outpath, file_pref_PROM, '_fets_sens'), row.names=F, col.names=T, sep= '\t', quote=F)


