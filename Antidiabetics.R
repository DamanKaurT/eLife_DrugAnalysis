rm(list = ls())
library(data.table)
library(lubridate)
library(dplyr)
library(reshape)
library(naniar)
library(PredictABEL)
library(car)

File <- file.path("C:/Users/Daman Kaur/Desktop/Drug classes DO")
setwd(File)
nacc <- fread("MDallmedscdr.csv") 
df <- nacc 

df <- as.data.frame(df)
df <- df[,c(1,8:10,14,26,49:88,118,126:137,145)]

## Find drug terms based on category
df$SULF <- apply(df, 1, function(x)as.integer(any(grep("CHLORPROPAMIDE|ACETOHEXAMIDE|GLIPIZIDE|
                    GLYBURIDE|TOLAZAMIDE|TOLBUTAMIDE|GLIMEPIRIDE",x))))

df$THI <- apply(df, 1, function(x)as.integer(any(grep("TROGLITAZONE|ROSIGLITAZONE|PIOGLITAZONE",x))))

df$COM <- apply(df, 1, function(x)as.integer(any(grep("metformin-|-metformin|glimepiride-|
                       -sitagliptin|-linagliptin",x, ignore.case = TRUE))))

df$SGLT <- apply(df, 1, function(x)as.integer(any(grep("EMPAGLIFLOZIN|CANAGLIFLOZIN|DAPAGLIFLOZIN",x))))

df$DPP <- apply(df, 1, function(x)as.integer(any(grep("SITAGLIPTIN|SAXAGLIPTIN|LINAGLIPTIN|ALOGLIPTIN",
                                                      x))))

df$INS <- apply(df, 1, function(x)as.integer(any(grep("insulin",x, ignore.case = TRUE))))
df$MET <- apply(df, 1, function(x)as.integer(any(grep("METFORMIN",x))))
df$AGI <- apply(df, 1, function(x)as.integer(any(grep("ACARBOSE|MIGLITOL",x))))
df$IMM <- apply(df, 1, function(x)as.integer(any(grep("EXENATIDE|LIRAGLUTIDE|DULAGLUTIDE",x))))
df$MEG <- apply(df, 1, function(x)as.integer(any(grep("REPAGLINIDE|NATEGLINIDE",x))))

df <- df[,-(7:46)] 
df[df==-4]<-NA

################ Adjusting values ###################################

#### SULF ###
df <- df %>% 
  group_by(NACCID) %>%
  mutate(NSULF = sum(SULF, na.rm=TRUE))

df$SULF[(df$NSULF==0)]<- 0
df$SULF[(df$NSULF==1)]<- 9 
df$SULF[(df$NSULF>=2)]<- 1 
df <- as.data.frame(df)

#### THI ###
df <- df %>% 
  group_by(NACCID) %>%
  mutate(NTHI = sum(THI, na.rm=TRUE))

df$THI[(df$NTHI==0)]<- 0
df$THI[(df$NTHI==1)]<- 9 
df$THI[(df$NTHI>=2)]<- 1 
df <- as.data.frame(df)

#### COM ###
df <- df %>% 
  group_by(NACCID) %>%
  mutate(NCOM = sum(COM, na.rm=TRUE))

df$COM[(df$NCOM==0)]<- 0
df$COM[(df$NCOM==1)]<- 9 
df$COM[(df$NCOM>=2)]<- 1 
df <- as.data.frame(df)

#### SGLT ###
df <- df %>% 
  group_by(NACCID) %>%
  mutate(NSGLT = sum(SGLT, na.rm=TRUE))

df$SGLT[(df$NSGLT==0)]<- 0
df$SGLT[(df$NSGLT==1)]<- 9 
df$SGLT[(df$NSGLT>=2)]<- 1 
df <- as.data.frame(df)

#### DPP ###
df <- df %>% 
  group_by(NACCID) %>%
  mutate(NDPP = sum(DPP, na.rm=TRUE))

df$DPP[(df$NDPP==0)]<- 0
df$DPP[(df$NDPP==1)]<- 9 
df$DPP[(df$NDPP>=2)]<- 1 
df <- as.data.frame(df)

#### INS ###
df <- df %>% 
  group_by(NACCID) %>%
  mutate(NINS = sum(INS, na.rm=TRUE))

df$INS[(df$NINS==0)]<- 0
df$INS[(df$NINS==1)]<- 9 
df$INS[(df$NINS>=2)]<- 1 
df <- as.data.frame(df)

#### MET ###
df <- df %>% 
  group_by(NACCID) %>%
  mutate(NMET = sum(MET, na.rm=TRUE))

df$MET[(df$NMET==0)]<- 0
df$MET[(df$NMET==1)]<- 9 
df$MET[(df$NMET>=2)]<- 1 
df <- as.data.frame(df)

#### AGI ###
df <- df %>% 
  group_by(NACCID) %>%
  mutate(NAGI = sum(AGI, na.rm=TRUE))

df$AGI[(df$NAGI==0)]<- 0
df$AGI[(df$NAGI==1)]<- 9 
df$AGI[(df$NAGI>=2)]<- 1 
df <- as.data.frame(df)

#### IMM ###
df <- df %>% 
  group_by(NACCID) %>%
  mutate(NIMM = sum(IMM, na.rm=TRUE))

df$IMM[(df$NIMM==0)]<- 0
df$IMM[(df$NIMM==1)]<- 9 
df$IMM[(df$NIMM>=2)]<- 1 
df <- as.data.frame(df)

#### MEG ###
df <- df %>% 
  group_by(NACCID) %>%
  mutate(NMEG = sum(MEG, na.rm=TRUE))

df$MEG[(df$NMEG==0)]<- 0
df$MEG[(df$NMEG==1)]<- 9 
df$MEG[(df$NMEG>=2)]<- 1 

file <- as.data.frame(df)
#write.csv(df, file = "MDdiabcdr.csv",row.names=FALSE, na="")

################# Log reg #################################################

df <- file

df$RACE[df$RACE==50] <- NA
df$SEX[df$SEX==2] <- 0  # female 0, male 1
df <- df[df$NEW_VISITNUM == 1, ] 
df <- subset(df, NACCAGE>=40,
             select=c(NACCID:NMEG))
df <- as.data.frame(df)
df <- df %>% filter(SEX == 0)

df$Diagnosis[df$Diagnosis==1 & df$ALZP==1] <- 2
df <- df[df$Diagnosis != 1, ]
df$Diagnosis[df$Diagnosis==2] <- 1

## 9 means medication mentioned only once, decide if need to remove or not
df <- df %>% replace_with_na(replace = list(SULF = 9, THI = 9, COM = 9, SGLT = 9,
                                            DPP = 9, INS = 9, MET = 9, AGI = 9,
                                            IMM = 9, MEG = 9))


df.baseDHH <- df[,c(1:6,9,12:19,21:30)] 
df.baseDHH[ , 16:25 ][ df.baseDHH[ , 16:25 ] == 9 ] <- 1

df.baseDHH.compl <- df.baseDHH[complete.cases(df.baseDHH), ]
lg.baseDHH <- df.baseDHH.compl

lg.baseDHH$RACE <- as.factor(as.character(lg.baseDHH$RACE))
lg.baseDHH$TBI <- as.factor(as.character(lg.baseDHH$TBI))
lg.baseDHH$HEAR <- as.factor(as.character(lg.baseDHH$HEAR))
lg.baseDHH$ALCOHOL <- as.factor(as.character(lg.baseDHH$ALCOHOL))
lg.baseDHH$ABMI <- as.factor(as.character(lg.baseDHH$ABMI))
lg.baseDHH$SEX <- as.factor(as.character(lg.baseDHH$SEX))
lg.baseDHH$DIAB <- as.factor(as.character(lg.baseDHH$DIAB))
lg.baseDHH$HYPERT <- as.factor(as.character(lg.baseDHH$HYPERT))
lg.baseDHH$CVD <- as.factor(as.character(lg.baseDHH$CVD))
lg.baseDHH$DEPRSN <- as.factor(as.character(lg.baseDHH$DEPRSN))
lg.baseDHH$SULF <- as.factor(as.character(lg.baseDHH$SULF))
lg.baseDHH$THI <- as.factor(as.character(lg.baseDHH$THI))
lg.baseDHH$MEG <- as.factor(as.character(lg.baseDHH$MEG))
lg.baseDHH$COM <- as.factor(as.character(lg.baseDHH$COM))
lg.baseDHH$SGLT <- as.factor(as.character(lg.baseDHH$SGLT))
lg.baseDHH$DPP <- as.factor(as.character(lg.baseDHH$DPP))
lg.baseDHH$INS <- as.factor(as.character(lg.baseDHH$INS))
lg.baseDHH$AGI <- as.factor(as.character(lg.baseDHH$AGI))
lg.baseDHH$IMM <- as.factor(as.character(lg.baseDHH$IMM))
lg.baseDHH$MET <- as.factor(as.character(lg.baseDHH$MET))
lg.baseDHH$Diagnosis <- as.factor(as.character(lg.baseDHH$Diagnosis))
lg.baseDHH <- as.data.frame(lg.baseDHH)

riskmodel <- fitLogRegModel(lg.baseDHH, cOutcome=7,
                           cNonGenPreds=c(2:6,8:25), cNonGenPredsCat=c(2,4:6,8:10,13:25),
                          cGenPreds=c(0), cGenPredsCat=c(0))

#gender
riskmodel <- fitLogRegModel(lg.baseDHH, cOutcome=7,
                            cNonGenPreds=c(3:6,8:25), cNonGenPredsCat=c(4:6,8:10,13:25),
                            cGenPreds=c(0), cGenPredsCat=c(0))

vif <- vif(riskmodel)
stat <- cbind(rownames(vif), data.frame(vif, row.names=NULL))
write.csv(stat, file = "vif.csv",row.names=FALSE, na="")


summary(riskmodel)
oddr <- ORmultivariate(riskmodel,filename="medtest.txt")
predRisk <- predRisk(riskmodel)

sex <- cbind(as.data.frame(prop.table(table(data.frame(lg.baseDHH$SEX,
                                                       lg.baseDHH$Diagnosis)), 1)), as.data.frame(with(lg.baseDHH, table(SEX, Diagnosis))))

race <- cbind(as.data.frame(prop.table(table(data.frame(lg.baseDHH$RACE,
                                                        lg.baseDHH$Diagnosis)), 1)), as.data.frame(with(lg.baseDHH, table(RACE, Diagnosis))))

tbi <- cbind(as.data.frame(prop.table(table(data.frame(lg.baseDHH$TBI,
                                                       lg.baseDHH$Diagnosis)), 1)), as.data.frame(with(lg.baseDHH, table(TBI, Diagnosis))))

hear <- cbind(as.data.frame(prop.table(table(data.frame(lg.baseDHH$HEAR,
                                                        lg.baseDHH$Diagnosis)), 1)), as.data.frame(with(lg.baseDHH, table(HEAR, Diagnosis))))

alcohol <- cbind(as.data.frame(prop.table(table(data.frame(lg.baseDHH$ALCOHOL,
                                                           lg.baseDHH$Diagnosis)), 1)), as.data.frame(with(lg.baseDHH, table(ALCOHOL, Diagnosis))))

diab <- cbind(as.data.frame(prop.table(table(data.frame(lg.baseDHH$DIAB,
                                                        lg.baseDHH$Diagnosis)), 1)), as.data.frame(with(lg.baseDHH, table(DIAB, Diagnosis))))

hyper <- cbind(as.data.frame(prop.table(table(data.frame(lg.baseDHH$HYPERT,
                                                         lg.baseDHH$Diagnosis)), 1)), as.data.frame(with(lg.baseDHH, table(HYPERT, Diagnosis))))

sul <- cbind(as.data.frame(prop.table(table(data.frame(lg.baseDHH$SULF,
                                                       lg.baseDHH$Diagnosis)), 1)), as.data.frame(with(lg.baseDHH, table(SULF, Diagnosis))))

thi <- cbind(as.data.frame(prop.table(table(data.frame(lg.baseDHH$THI,
                                                       lg.baseDHH$Diagnosis)), 1)), as.data.frame(with(lg.baseDHH, table(THI, Diagnosis))))

com <- cbind(as.data.frame(prop.table(table(data.frame(lg.baseDHH$COM,
                                                       lg.baseDHH$Diagnosis)), 1)), as.data.frame(with(lg.baseDHH, table(COM, Diagnosis))))

sglt <- cbind(as.data.frame(prop.table(table(data.frame(lg.baseDHH$SGLT,
                                                       lg.baseDHH$Diagnosis)), 1)), as.data.frame(with(lg.baseDHH, table(SGLT, Diagnosis))))

dpp <- cbind(as.data.frame(prop.table(table(data.frame(lg.baseDHH$DPP,
                                                       lg.baseDHH$Diagnosis)), 1)), as.data.frame(with(lg.baseDHH, table(DPP, Diagnosis))))

ins <- cbind(as.data.frame(prop.table(table(data.frame(lg.baseDHH$INS,
                                                       lg.baseDHH$Diagnosis)), 1)), as.data.frame(with(lg.baseDHH, table(INS, Diagnosis))))

met <- cbind(as.data.frame(prop.table(table(data.frame(lg.baseDHH$MET,
                                                       lg.baseDHH$Diagnosis)), 1)), as.data.frame(with(lg.baseDHH, table(MET, Diagnosis))))

agi <- cbind(as.data.frame(prop.table(table(data.frame(lg.baseDHH$AGI,
                                                       lg.baseDHH$Diagnosis)), 1)), as.data.frame(with(lg.baseDHH, table(AGI, Diagnosis))))

imm <- cbind(as.data.frame(prop.table(table(data.frame(lg.baseDHH$IMM,
                                                       lg.baseDHH$Diagnosis)), 1)), as.data.frame(with(lg.baseDHH, table(IMM, Diagnosis))))

meg <- cbind(as.data.frame(prop.table(table(data.frame(lg.baseDHH$MEG,
                                                       lg.baseDHH$Diagnosis)), 1)), as.data.frame(with(lg.baseDHH, table(MEG, Diagnosis))))

bmi <- cbind(as.data.frame(prop.table(table(data.frame(lg.baseDHH$ABMI,
                                                       lg.baseDHH$Diagnosis)), 1)), as.data.frame(with(lg.baseDHH, table(ABMI, Diagnosis))))

cvd <- cbind(as.data.frame(prop.table(table(data.frame(lg.baseDHH$CVD,
                                                       lg.baseDHH$Diagnosis)), 1)), as.data.frame(with(lg.baseDHH, table(CVD, Diagnosis))))

dep <- cbind(as.data.frame(prop.table(table(data.frame(lg.baseDHH$DEPRSN,
                                                       lg.baseDHH$Diagnosis)), 1)), as.data.frame(with(lg.baseDHH, table(DEPRSN, Diagnosis))))

xy <- cbind(sex, alcohol, hear,diab,sul,thi, com,sglt,dpp,ins,met,agi,imm,meg,hyper,dep,tbi,cvd)

write.excel <- function(x,row.names=FALSE,col.names=TRUE,...) {
  write.table(x,"clipboard",sep="\t",row.names=row.names,col.names=col.names,...)}

write.excel(xy)

write.excel(race)
write.excel(bmi)
write.excel(sex)
write.excel(alcohol)
write.excel(hear)
write.excel(diab)
write.excel(sul)
write.excel(thi)
write.excel(com)
write.excel(sglt)
write.excel(dpp)
write.excel(ins)
write.excel(met)
write.excel(agi)
write.excel(imm)
write.excel(meg)
write.excel(hyper)
write.excel(dep)
write.excel(tbi)
write.excel(cvd)
