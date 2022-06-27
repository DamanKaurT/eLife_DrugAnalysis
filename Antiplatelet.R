##### Searching drug terms ####

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
df <- df[,c(1,8:10,14,26,49:88,109,126:137,141)]

#heparin, VitK, Direct Xa, Thrombin IIa, Phosphodiesterase, ADP, Aspirin

## Find drug terms based on category
df$HEPA <- apply(df, 1, function(x)as.integer(any(grep("heparin|enoxaparin|dalteparin|fondaparinux",x))))

df$VITK <- apply(df, 1, function(x)as.integer(any(grep("warfarin|dicumarol",x, ignore.case = TRUE))))

df$XAIN <- apply(df, 1, function(x)as.integer(any(grep("rivaroxaban|apixaban",x, ignore.case = TRUE))))

df$THII <- apply(df, 1, function(x)as.integer(any(grep("dabigatran",x, ignore.case = TRUE))))

df$PHOS <- apply(df, 1, function(x)as.integer(any(grep("dipyridamole|cilostazol",x, ignore.case = TRUE))))
df$ADP <- apply(df, 1, function(x)as.integer(any(grep("ticlopidine|clopidogrel|prasugrel|
                                                       ticagrelor",x, ignore.case = TRUE))))
df$ASPR <- apply(df, 1, function(x)as.integer(any(grep("aspirin|asa/|aspirin-",x, ignore.case = TRUE))))

df <- df[,-(7:46)] 
df[df==-4]<-NA

################ Adjusting values ###################################

#### HEPA ###
df <- df %>% 
  group_by(NACCID) %>%
  mutate(NHEPA = sum(HEPA, na.rm=TRUE))

df$HEPA[(df$NHEPA==0)]<- 0
df$HEPA[(df$NHEPA==1)]<- 9 
df$HEPA[(df$NHEPA>=2)]<- 1 
df <- as.data.frame(df)

#### VITK ###
df <- df %>% 
  group_by(NACCID) %>%
  mutate(NVITK = sum(VITK, na.rm=TRUE))

df$VITK[(df$NVITK==0)]<- 0
df$VITK[(df$NVITK==1)]<- 9 
df$VITK[(df$NVITK>=2)]<- 1 
df <- as.data.frame(df)

#### XAIN ###
df <- df %>% 
  group_by(NACCID) %>%
  mutate(NXAIN = sum(XAIN, na.rm=TRUE))

df$XAIN[(df$NXAIN==0)]<- 0
df$XAIN[(df$NXAIN==1)]<- 9 
df$XAIN[(df$NXAIN>=2)]<- 1 
df <- as.data.frame(df)

#### THII ###
df <- df %>% 
  group_by(NACCID) %>%
  mutate(NTHII = sum(THII, na.rm=TRUE))

df$THII[(df$NTHII==0)]<- 0
df$THII[(df$NTHII==1)]<- 9 
df$THII[(df$NTHII>=2)]<- 1 
df <- as.data.frame(df)

#### PHOS ###
df <- df %>% 
  group_by(NACCID) %>%
  mutate(NPHOS = sum(PHOS, na.rm=TRUE))

df$PHOS[(df$NPHOS==0)]<- 0
df$PHOS[(df$NPHOS==1)]<- 9 
df$PHOS[(df$NPHOS>=2)]<- 1 
df <- as.data.frame(df)

#### ADP ###
df <- df %>% 
  group_by(NACCID) %>%
  mutate(NADP = sum(ADP, na.rm=TRUE))

df$ADP[(df$NADP==0)]<- 0
df$ADP[(df$NADP==1)]<- 9 
df$ADP[(df$NADP>=2)]<- 1 
df <- as.data.frame(df)

#### ASPR ###
df <- df %>% 
  group_by(NACCID) %>%
  mutate(NASPR = sum(ASPR, na.rm=TRUE))

df$ASPR[(df$NASPR==0)]<- 0
df$ASPR[(df$NASPR==1)]<- 9 
df$ASPR[(df$NASPR>=2)]<- 1 
file <- as.data.frame(df)

#write.csv(df, file = "MDaccdr.csv",row.names=FALSE, na="")

################# Log reg #################################################

df <- file
df$RACE[df$RACE==50] <- NA
df$SEX[df$SEX==2] <- 0  # female 0, male 1
df <- df[df$NEW_VISITNUM == 1, ] 
df <- subset(df, NACCAGE>=40,
             select=c(NACCID:NASPR))
df <- as.data.frame(df)
df <- df %>% filter(SEX == 1)

df$Diagnosis[df$Diagnosis==1 & df$ALZP==1] <- 2
df <- df[df$Diagnosis != 1, ]
df$Diagnosis[df$Diagnosis==2] <- 1

## 9 means medication mentioned only once, decide if need to remove or not
df <- df %>% replace_with_na(replace = list(HEPA = 9, VITK = 9, XAIN = 9, PHOS = 9,
                                           THII = 9, ADP = 9, ASPR = 9))


df.baseDHH <- df[,c(1:6,9,12:19,21:27)] 
df.baseDHH[ , 17:22 ][ df.baseDHH[ , 17:22 ] == 9 ] <- 1

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
lg.baseDHH$HEPA <- as.factor(as.character(lg.baseDHH$HEPA))
lg.baseDHH$VITK <- as.factor(as.character(lg.baseDHH$VITK))
lg.baseDHH$XAIN <- as.factor(as.character(lg.baseDHH$XAIN))
lg.baseDHH$PHOS <- as.factor(as.character(lg.baseDHH$PHOS))
lg.baseDHH$THII <- as.factor(as.character(lg.baseDHH$THII))
lg.baseDHH$ADP <- as.factor(as.character(lg.baseDHH$ADP))
lg.baseDHH$ASPR <- as.factor(as.character(lg.baseDHH$ASPR ))
lg.baseDHH$Diagnosis <- as.factor(as.character(lg.baseDHH$Diagnosis))
lg.baseDHH <- as.data.frame(lg.baseDHH)

riskmodel <- fitLogRegModel(lg.baseDHH, cOutcome=7,
                           cNonGenPreds=c(2:6,8:15,17:22), cNonGenPredsCat=c(2,4:6,8:10,13:15,17:22),
                          cGenPreds=c(0), cGenPredsCat=c(0))

#gender
riskmodel <- fitLogRegModel(lg.baseDHH, cOutcome=7,
                            cNonGenPreds=c(3:6,8:15,17:22), cNonGenPredsCat=c(4:6,8:10,13:15,17:22),
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

hepa <- cbind(as.data.frame(prop.table(table(data.frame(lg.baseDHH$HEPA,
                                                          lg.baseDHH$Diagnosis)), 1)), as.data.frame(with(lg.baseDHH, table(HEPA, Diagnosis))))

vitk <- cbind(as.data.frame(prop.table(table(data.frame(lg.baseDHH$VITK,
                                                        lg.baseDHH$Diagnosis)), 1)), as.data.frame(with(lg.baseDHH, table(VITK, Diagnosis))))

xain <- cbind(as.data.frame(prop.table(table(data.frame(lg.baseDHH$XAIN,
                                                        lg.baseDHH$Diagnosis)), 1)), as.data.frame(with(lg.baseDHH, table(XAIN, Diagnosis))))

phos <- cbind(as.data.frame(prop.table(table(data.frame(lg.baseDHH$PHOS,
                                                       lg.baseDHH$Diagnosis)), 1)), as.data.frame(with(lg.baseDHH, table(PHOS, Diagnosis))))

thii <- cbind(as.data.frame(prop.table(table(data.frame(lg.baseDHH$THII,
                                                        lg.baseDHH$Diagnosis)), 1)), as.data.frame(with(lg.baseDHH, table(THII, Diagnosis))))

adp <- cbind(as.data.frame(prop.table(table(data.frame(lg.baseDHH$ADP,
                                                       lg.baseDHH$Diagnosis)), 1)), as.data.frame(with(lg.baseDHH, table(ADP, Diagnosis))))

aspr <- cbind(as.data.frame(prop.table(table(data.frame(lg.baseDHH$ASPR,
                                                        lg.baseDHH$Diagnosis)), 1)), as.data.frame(with(lg.baseDHH, table(ASPR, Diagnosis))))

bmi <- cbind(as.data.frame(prop.table(table(data.frame(lg.baseDHH$ABMI,
                                                       lg.baseDHH$Diagnosis)), 1)), as.data.frame(with(lg.baseDHH, table(ABMI, Diagnosis))))

cvd <- cbind(as.data.frame(prop.table(table(data.frame(lg.baseDHH$CVD,
                                                       lg.baseDHH$Diagnosis)), 1)), as.data.frame(with(lg.baseDHH, table(CVD, Diagnosis))))

dep <- cbind(as.data.frame(prop.table(table(data.frame(lg.baseDHH$DEPRSN,
                                                       lg.baseDHH$Diagnosis)), 1)), as.data.frame(with(lg.baseDHH, table(DEPRSN, Diagnosis))))

xy <- cbind(sex,alcohol,hear,diab,hepa,vitk,xain,thii,phos,adp,aspr,hyper,dep,tbi,cvd)

write.excel <- function(x,row.names=FALSE,col.names=TRUE,...) {
  write.table(x,"clipboard",sep="\t",row.names=row.names,col.names=col.names,...)}

write.excel(xy)

write.excel(race)
write.excel(bmi)

####
write.excel(sex)
write.excel(race)
write.excel(alcohol)
write.excel(hear)
write.excel(hepa)
write.excel(vitk)
write.excel(xain)
write.excel(phos)
write.excel(thii)
write.excel(adp)
write.excel(aspr)
write.excel(diab)
write.excel(hyper) 
write.excel(dep)
write.excel(bmi)
write.excel(tbi)
write.excel(cvd)
