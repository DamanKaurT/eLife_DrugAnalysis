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
  mutate(check1=is.na(PARK1) & is.na(PARK2))
df <- df %>% group_by(NACCID) %>% filter(!all(check1=="TRUE"))
df <- df[,-(148)] 
df <- as.data.frame(df)

df <- df %>% 
  group_by(NACCID) %>%
  mutate(DCDC = sum(PARK1, na.rm=TRUE))
df$PARK1[(df$DCDC==0 & df$NEW_VISITNUM==1)]<- 0 
df$PARK1[(df$DCDC>=1 & df$NEW_VISITNUM==1)]<- 1 

df <- df %>% 
  group_by(NACCID) %>%
  mutate(DCDC = sum(PARK2, na.rm=TRUE))
df$PARK2[(df$DCDC==0 & df$NEW_VISITNUM==1)]<- 0 
df$PARK2[(df$DCDC>=1 & df$NEW_VISITNUM==1)]<- 1 

df <- df[,-(148)] 
df <- as.data.frame(df)

df$PARK[(df$PARK2==1 | df$PARK2==1)] <- 1
df$PARK[(df$PARK2==0 & df$PARK2==0)] <- 0

df <- as.data.frame(df)
df <- df[,c(1,8:10,14,26,49:88,123,126:137,147,148)]

## Find drug terms based on category
df$LEV <- apply(df, 1, function(x)as.integer(any(grep("LEVODOPA",x, ignore.case = TRUE))))

df$DORA <- apply(df, 1, function(x)as.integer(any(grep("bromocriptine|pergolide|cabergoline|
                  apomorphine|pramipexole|ropinirole|rotigotine",x, ignore.case = TRUE))))

df$ACH <- apply(df, 1, function(x)as.integer(any(grep("benzotropine|procyclidine|biperiden|
                                           trihexylphenidyl|diphenhydramine",x, ignore.case = TRUE))))

df$COM <- apply(df, 1, function(x)as.integer(any(grep("TOLCAPONE|ENTACAPONE",x))))
df$MAO <- apply(df, 1, function(x)as.integer(any(grep("selegiline|rasagiline",x, ignore.case = TRUE))))
df$AMA <- apply(df, 1, function(x)as.integer(any(grep("amantidine",x, ignore.case = TRUE))))

df <- df[,-(7:46)] 
df[df==-4]<-NA

################ Adjusting values ###################################

#### LEV ###
df <- df %>% 
  group_by(NACCID) %>%
  mutate(NLEV = sum(LEV, na.rm=TRUE))

df$LEV[(df$NLEV==0)]<- 0
df$LEV[(df$NLEV==1)]<- 9 
df$LEV[(df$NLEV>=2)]<- 1 
df <- as.data.frame(df)

#### DORA ###
df <- df %>% 
  group_by(NACCID) %>%
  mutate(NDORA = sum(DORA, na.rm=TRUE))

df$DORA[(df$NDORA==0)]<- 0
df$DORA[(df$NDORA==1)]<- 9 
df$DORA[(df$NDORA>=2)]<- 1 
df <- as.data.frame(df)

#### ACH ###
df <- df %>% 
  group_by(NACCID) %>%
  mutate(NACH = sum(ACH, na.rm=TRUE))

df$ACH[(df$NACH==0)]<- 0
df$ACH[(df$NACH==1)]<- 9 
df$ACH[(df$NACH>=2)]<- 1 
df <- as.data.frame(df)

#### MAO ###
df <- df %>% 
  group_by(NACCID) %>%
  mutate(NMAO = sum(MAO, na.rm=TRUE))

df$MAO[(df$NMAO==0)]<- 0
df$MAO[(df$NMAO==1)]<- 9 
df$MAO[(df$NMAO>=2)]<- 1 
df <- as.data.frame(df)

#### AMA ###
df <- df %>% 
  group_by(NACCID) %>%
  mutate(NAMA = sum(AMA, na.rm=TRUE))

df$AMA[(df$NAMA==0)]<- 0
df$AMA[(df$NAMA==1)]<- 9 
df$AMA[(df$NAMA>=2)]<- 1 
df <- as.data.frame(df)

#### COM ###
df <- df %>% 
  group_by(NACCID) %>%
  mutate(NCOM = sum(COM, na.rm=TRUE))

df$COM[(df$NCOM==0)]<- 0
df$COM[(df$NCOM==1)]<- 9 
df$COM[(df$NCOM>=2)]<- 1 
file <- as.data.frame(df)

#write.csv(df, file = "NMpdmcdr.csv",row.names=FALSE, na="")

################# Log reg #################################################
df <- file

df$RACE[df$RACE==50] <- NA
df$SEX[df$SEX==2] <- 0  # female 0, male 1
df <- df[df$NEW_VISITNUM == 1, ] 
df <- subset(df, NACCAGE>=40,
             select=c(NACCID:NCOM))
df <- as.data.frame(df)
df <- df %>% filter(SEX == 1)

df$Diagnosis[df$Diagnosis==1 & df$ALZP==1] <- 2
df <- df[df$Diagnosis != 1, ]
df$Diagnosis[df$Diagnosis==2] <- 1

## 9 means medication mentioned only once, decide if need to remove or not
df <- df %>% replace_with_na(replace = list(LEV = 9, DORA = 9, ACH = 9, COM = 9,
                                            MAO = 9, AMA = 9))


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
lg.baseDHH$PARK <- as.factor(as.character(lg.baseDHH$PARK))
lg.baseDHH$LEV <- as.factor(as.character(lg.baseDHH$LEV))
lg.baseDHH$DORA <- as.factor(as.character(lg.baseDHH$DORA))
lg.baseDHH$ACH <- as.factor(as.character(lg.baseDHH$ACH))
lg.baseDHH$COM <- as.factor(as.character(lg.baseDHH$COM))
lg.baseDHH$MAO <- as.factor(as.character(lg.baseDHH$MAO))
lg.baseDHH$AMA <- as.factor(as.character(lg.baseDHH$AMA))
lg.baseDHH$Diagnosis <- as.factor(as.character(lg.baseDHH$Diagnosis))
lg.baseDHH <- as.data.frame(lg.baseDHH)

riskmodel <- fitLogRegModel(lg.baseDHH, cOutcome=7,
                           cNonGenPreds=c(2:6,8:21), cNonGenPredsCat=c(2,4:6,8:10,13:21),
                            cGenPreds=c(0), cGenPredsCat=c(0))

#gender
riskmodel <- fitLogRegModel(lg.baseDHH, cOutcome=7,
                            cNonGenPreds=c(3:6,8:21), cNonGenPredsCat=c(4:6,8:10,13:21),
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

park <- cbind(as.data.frame(prop.table(table(data.frame(lg.baseDHH$PARK,
                                                       lg.baseDHH$Diagnosis)), 1)), as.data.frame(with(lg.baseDHH, table(PARK, Diagnosis))))

lev <- cbind(as.data.frame(prop.table(table(data.frame(lg.baseDHH$LEV,
                                                       lg.baseDHH$Diagnosis)), 1)), as.data.frame(with(lg.baseDHH, table(LEV, Diagnosis))))

dora <- cbind(as.data.frame(prop.table(table(data.frame(lg.baseDHH$DORA,
                                                       lg.baseDHH$Diagnosis)), 1)), as.data.frame(with(lg.baseDHH, table(DORA, Diagnosis))))

ach <- cbind(as.data.frame(prop.table(table(data.frame(lg.baseDHH$ACH,
                                                       lg.baseDHH$Diagnosis)), 1)), as.data.frame(with(lg.baseDHH, table(ACH, Diagnosis))))

com <- cbind(as.data.frame(prop.table(table(data.frame(lg.baseDHH$COM,
                                                       lg.baseDHH$Diagnosis)), 1)), as.data.frame(with(lg.baseDHH, table(COM, Diagnosis))))

mao <- cbind(as.data.frame(prop.table(table(data.frame(lg.baseDHH$MAO,
                                                       lg.baseDHH$Diagnosis)), 1)), as.data.frame(with(lg.baseDHH, table(MAO, Diagnosis))))

ama <- cbind(as.data.frame(prop.table(table(data.frame(lg.baseDHH$AMA,
                                                       lg.baseDHH$Diagnosis)), 1)), as.data.frame(with(lg.baseDHH, table(AMA, Diagnosis))))

bmi <- cbind(as.data.frame(prop.table(table(data.frame(lg.baseDHH$ABMI,
                                                       lg.baseDHH$Diagnosis)), 1)), as.data.frame(with(lg.baseDHH, table(ABMI, Diagnosis))))

cvd <- cbind(as.data.frame(prop.table(table(data.frame(lg.baseDHH$CVD,
                                                       lg.baseDHH$Diagnosis)), 1)), as.data.frame(with(lg.baseDHH, table(CVD, Diagnosis))))

dep <- cbind(as.data.frame(prop.table(table(data.frame(lg.baseDHH$DEPRSN,
                                                       lg.baseDHH$Diagnosis)), 1)), as.data.frame(with(lg.baseDHH, table(DEPRSN, Diagnosis))))

xy <- cbind(sex, alcohol, hear, diab,park,lev,dora,ach,com,mao,ama,hyper,dep,tbi,cvd)

write.excel <- function(x,row.names=FALSE,col.names=TRUE,...) {
  write.table(x,"clipboard",sep="\t",row.names=row.names,col.names=col.names,...)}

write.excel(xy)

write.excel(race)
write.excel(bmi)
write.excel(sex)
write.excel(alcohol)
write.excel(hear)
write.excel(park)
write.excel(lev)
write.excel(dora)
write.excel(ach)
write.excel(com)
write.excel(mao)
write.excel(ama)
write.excel(diab)
write.excel(hyper)
write.excel(dep)
write.excel(tbi)
write.excel(cvd)
