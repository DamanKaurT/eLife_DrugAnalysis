##### Searching drug terms #####

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

df <- df %>%
  group_by(NACCID) %>%
  mutate(check1=is.na(APS1) & is.na(APS2) & is.na(APS3) & is.na(APS4) & is.na(APS5))
df <- df %>% group_by(NACCID) %>% filter(!all(check1=="TRUE"))
df <- df[,-(148)] 
df <- as.data.frame(df)

df <- df %>% 
  group_by(NACCID) %>%
  mutate(DCDC = sum(APS1, na.rm=TRUE))
df$APS1[(df$DCDC==0 & df$NEW_VISITNUM==1)]<- 0 
df$APS1[(df$DCDC>=1 & df$NEW_VISITNUM==1)]<- 1 

df <- df %>% 
  group_by(NACCID) %>%
  mutate(DCDC = sum(APS2, na.rm=TRUE))
df$APS2[(df$DCDC==0 & df$NEW_VISITNUM==1)]<- 0 
df$APS2[(df$DCDC>=1 & df$NEW_VISITNUM==1)]<- 1 

df <- df %>% 
  group_by(NACCID) %>%
  mutate(DCDC = sum(APS3, na.rm=TRUE))
df$APS3[(df$DCDC==0 & df$NEW_VISITNUM==1)]<- 0 
df$APS3[(df$DCDC>=1 & df$NEW_VISITNUM==1)]<- 1 

df <- df %>% 
  group_by(NACCID) %>%
  mutate(DCDC = sum(APS4, na.rm=TRUE))
df$APS4[(df$DCDC==0 & df$NEW_VISITNUM==1)]<- 0 
df$APS4[(df$DCDC>=1 & df$NEW_VISITNUM==1)]<- 1 

df <- df %>% 
  group_by(NACCID) %>%
  mutate(DCDC = sum(APS5, na.rm=TRUE))
df$APS5[(df$DCDC==0 & df$NEW_VISITNUM==1)]<- 0 
df$APS5[(df$DCDC>=1 & df$NEW_VISITNUM==1)]<- 1 

df <- df[,-(148)] 
df <- as.data.frame(df)

df$PSYC[(df$APS1==1 | df$APS2==1 | df$APS3==1 | df$APS4==1 | df$APS5==1)] <- 1
df$PSYC[(df$APS1==0 & df$APS2==0 & df$APS3==0 & df$APS4==0 & df$APS5==0)] <- 0

df <- as.data.frame(df)
df <- df[,c(1,8:10,14,26,49:88,115,126:137,143,148)]

## Find drug terms based on category
df$HALO <- apply(df, 1, function(x)as.integer(any(grep("HALOPERIDOL",x))))

df$LITH <- apply(df, 1, function(x)as.integer(any(grep("LITHIUM",x))))

df$PHEN <- apply(df, 1, function(x)as.integer(any(grep("CHLORPROMAZINE|FLUPHENAZINE|PROCHLORPERAZINE|PROMAZINE|
                                                       TRIFLUPROMAZINE|THIORIDAZINE|PERPHENAZINE|MESORIDAZINE|
                                                       METHOTRIEPRAZINE|TRIFLUOPERAZINE",x))))

df$TRC <- apply(df, 1, function(x)as.integer(any(grep("CLOZAPINE|OLANZAPINE|QUETIAPINE|ASENAPINE|LOXAPINE",x))))

df$BENS <- apply(df, 1, function(x)as.integer(any(grep("RISPERIDONE|ZIPRASIDONE|PALIPERIDONE|
                                                       ILOPERIDONE|LURASIDONE",x))))
df$ARIP <- apply(df, 1, function(x)as.integer(any(grep("ARIPIPRAZOLE",x))))
df$COMBO <- apply(df, 1, function(x)as.integer(any(grep("fluoxetine-|amitriptyline-",x, ignore.case = TRUE))))

df <- df[,-(7:46)] 
df[df==-4]<-NA

################ Adjusting values ###################################

#### HALO ###
df <- df %>% 
  group_by(NACCID) %>%
  mutate(NHALO = sum(HALO, na.rm=TRUE))

df$HALO[(df$NHALO==0)]<- 0
df$HALO[(df$NHALO==1)]<- 9 
df$HALO[(df$NHALO>=2)]<- 1 
df <- as.data.frame(df)

#### LITH ###
df <- df %>% 
  group_by(NACCID) %>%
  mutate(NLITH = sum(LITH, na.rm=TRUE))

df$LITH[(df$NLITH==0)]<- 0
df$LITH[(df$NLITH==1)]<- 9 
df$LITH[(df$NLITH>=2)]<- 1 
df <- as.data.frame(df)

#### PHEN ###
df <- df %>% 
  group_by(NACCID) %>%
  mutate(NPHEN = sum(PHEN, na.rm=TRUE))

df$PHEN[(df$NPHEN==0)]<- 0
df$PHEN[(df$NPHEN==1)]<- 9 
df$PHEN[(df$NPHEN>=2)]<- 1 
df <- as.data.frame(df)

#### TRC ###
df <- df %>% 
  group_by(NACCID) %>%
  mutate(NTRC = sum(TRC, na.rm=TRUE))

df$TRC[(df$NTRC==0)]<- 0
df$TRC[(df$NTRC==1)]<- 9 
df$TRC[(df$NTRC>=2)]<- 1 
df <- as.data.frame(df)

#### BENS ###
df <- df %>% 
  group_by(NACCID) %>%
  mutate(NBENS = sum(BENS, na.rm=TRUE))

df$BENS[(df$NBENS==0)]<- 0
df$BENS[(df$NBENS==1)]<- 9 
df$BENS[(df$NBENS>=2)]<- 1 
df <- as.data.frame(df)

#### ARIP ###
df <- df %>% 
  group_by(NACCID) %>%
  mutate(NARIP = sum(ARIP, na.rm=TRUE))

df$ARIP[(df$NARIP==0)]<- 0
df$ARIP[(df$NARIP==1)]<- 9 
df$ARIP[(df$NARIP>=2)]<- 1 
df <- as.data.frame(df)

#### COMBO ###
df <- df %>% 
  group_by(NACCID) %>%
  mutate(NCOMBO = sum(COMBO, na.rm=TRUE))

df$COMBO[(df$NCOMBO==0)]<- 0
df$COMBO[(df$NCOMBO==1)]<- 9 
df$COMBO[(df$NCOMBO>=2)]<- 1 
file <- as.data.frame(df)

#write.csv(df, file = "NMapsycdr.csv",row.names=FALSE, na="")

################# Log reg #################################################

df <- file
df$RACE[df$RACE==50] <- NA
df$SEX[df$SEX==2] <- 0  # female 0, male 1
df <- df[df$NEW_VISITNUM == 1, ] 
df <- subset(df, NACCAGE>=40,
             select=c(NACCID:NCOMBO))
df <- as.data.frame(df)
df <- df %>% filter(SEX == 1)

df$Diagnosis[df$Diagnosis==1 & df$ALZP==1] <- 2
df <- df[df$Diagnosis != 1, ]
df$Diagnosis[df$Diagnosis==2] <- 1

## 9 means medication mentioned only once, decide if need to remove or not
df <- df %>% replace_with_na(replace = list(HALO = 9, LITH = 9, PHEN = 9, TRC = 9,
                                           BENS = 9, ARIP = 9, COMBO = 9))


df.baseDHH <- df[,c(1:6,9,12:19,21:28)] 
df.baseDHH[ , 17:23][ df.baseDHH[ , 17:23 ] == 9 ] <- 1

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
lg.baseDHH$PSYC <- as.factor(as.character(lg.baseDHH$PSYC))
lg.baseDHH$HALO <- as.factor(as.character(lg.baseDHH$HALO))
lg.baseDHH$LITH <- as.factor(as.character(lg.baseDHH$LITH))
lg.baseDHH$PHEN <- as.factor(as.character(lg.baseDHH$PHEN))
lg.baseDHH$TRC <- as.factor(as.character(lg.baseDHH$TRC))
lg.baseDHH$BENS <- as.factor(as.character(lg.baseDHH$BENS))
lg.baseDHH$ARIP <- as.factor(as.character(lg.baseDHH$ARIP))
lg.baseDHH$COMBO <- as.factor(as.character(lg.baseDHH$COMBO ))
lg.baseDHH$Diagnosis <- as.factor(as.character(lg.baseDHH$Diagnosis))
lg.baseDHH <- as.data.frame(lg.baseDHH)

riskmodel <- fitLogRegModel(lg.baseDHH, cOutcome=7,
              cNonGenPreds=c(2:6,8:16,18:23), cNonGenPredsCat=c(2,4:6,8:10,13:16,18:23),
              cGenPreds=c(0), cGenPredsCat=c(0))

#gender
riskmodel <- fitLogRegModel(lg.baseDHH, cOutcome=7,
                            cNonGenPreds=c(3:6,8:16,18:23), cNonGenPredsCat=c(4:6,8:10,13:16,18:23),
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

psyc <- cbind(as.data.frame(prop.table(table(data.frame(lg.baseDHH$PSYC,
                                                       lg.baseDHH$Diagnosis)), 1)), as.data.frame(with(lg.baseDHH, table(PSYC, Diagnosis))))

hal <- cbind(as.data.frame(prop.table(table(data.frame(lg.baseDHH$HALO,
      lg.baseDHH$Diagnosis)), 1)), as.data.frame(with(lg.baseDHH, table(HALO, Diagnosis))))

lit <- cbind(as.data.frame(prop.table(table(data.frame(lg.baseDHH$LITH,
      lg.baseDHH$Diagnosis)), 1)), as.data.frame(with(lg.baseDHH, table(LITH, Diagnosis))))

phen <- cbind(as.data.frame(prop.table(table(data.frame(lg.baseDHH$PHEN,
    lg.baseDHH$Diagnosis)), 1)), as.data.frame(with(lg.baseDHH, table(PHEN, Diagnosis))))

trc <- cbind(as.data.frame(prop.table(table(data.frame(lg.baseDHH$TRC,
      lg.baseDHH$Diagnosis)), 1)), as.data.frame(with(lg.baseDHH, table(TRC, Diagnosis))))

ben <- cbind(as.data.frame(prop.table(table(data.frame(lg.baseDHH$BENS,
    lg.baseDHH$Diagnosis)), 1)), as.data.frame(with(lg.baseDHH, table(BENS, Diagnosis))))

ari <- cbind(as.data.frame(prop.table(table(data.frame(lg.baseDHH$ARIP,
    lg.baseDHH$Diagnosis)), 1)), as.data.frame(with(lg.baseDHH, table(ARIP, Diagnosis))))

com <- cbind(as.data.frame(prop.table(table(data.frame(lg.baseDHH$COMBO,
                                                       lg.baseDHH$Diagnosis)), 1)), as.data.frame(with(lg.baseDHH, table(COMBO, Diagnosis))))

bmi <- cbind(as.data.frame(prop.table(table(data.frame(lg.baseDHH$ABMI,
                                                       lg.baseDHH$Diagnosis)), 1)), as.data.frame(with(lg.baseDHH, table(ABMI, Diagnosis))))

cvd <- cbind(as.data.frame(prop.table(table(data.frame(lg.baseDHH$CVD,
                                                       lg.baseDHH$Diagnosis)), 1)), as.data.frame(with(lg.baseDHH, table(CVD, Diagnosis))))

dep <- cbind(as.data.frame(prop.table(table(data.frame(lg.baseDHH$DEPRSN,
                                                       lg.baseDHH$Diagnosis)), 1)), as.data.frame(with(lg.baseDHH, table(DEPRSN, Diagnosis))))

xy <- cbind(sex,alcohol,hear,diab,psyc,hal,lit,phen,trc,ben,ari,com,hyper,dep,tbi,cvd)

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
write.excel(hal)
write.excel(lit)
write.excel(phen)
write.excel(trc)
write.excel(ben)
write.excel(ari)
write.excel(com)
write.excel(diab)
write.excel(hyper) 
write.excel(dep)
write.excel(bmi)
write.excel(tbi)
write.excel(cvd) 
