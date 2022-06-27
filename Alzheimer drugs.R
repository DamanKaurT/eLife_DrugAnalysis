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
df <- df[,c(1,8:10,14,26,49:88,112,126:137,146)]

## Find drug terms based on category
df$DON <- apply(df, 1, function(x)as.integer(any(grep("DONEPEZIL",x))))
df$MEM <- apply(df, 1, function(x)as.integer(any(grep("MEMANTINE",x))))
df$GAL <- apply(df, 1, function(x)as.integer(any(grep("GALANTAMINE",x, ignore.case = TRUE))))
df$RIV <- apply(df, 1, function(x)as.integer(any(grep("RIVASTIGMINE",x, ignore.case = TRUE))))
df$TAC <- apply(df, 1, function(x)as.integer(any(grep("TACRINE",x, ignore.case = TRUE))))
df$COM <- apply(df, 1, function(x)as.integer(any(grep("donepezil-memantine|donepezil and memantine",x, 
                                                      ignore.case = TRUE))))

df <- df[,-(7:46)] 
df[df==-4]<-NA

################ Adjusting values ###################################

#### DON ###
df <- df %>% 
  group_by(NACCID) %>%
  mutate(NDON = sum(DON, na.rm=TRUE))

df$DON[(df$NDON==0)]<- 0
df$DON[(df$NDON==1)]<- 9 
df$DON[(df$NDON>=2)]<- 1 
df <- as.data.frame(df)

#### MEM ###
df <- df %>% 
  group_by(NACCID) %>%
  mutate(NMEM = sum(MEM, na.rm=TRUE))

df$MEM[(df$NMEM==0)]<- 0
df$MEM[(df$NMEM==1)]<- 9 
df$MEM[(df$NMEM>=2)]<- 1 
df <- as.data.frame(df)

#### GAL ###
df <- df %>% 
  group_by(NACCID) %>%
  mutate(NGAL = sum(GAL, na.rm=TRUE))

df$GAL[(df$NGAL==0)]<- 0
df$GAL[(df$NGAL==1)]<- 9 
df$GAL[(df$NGAL>=2)]<- 1 
df <- as.data.frame(df)

#### RIV ###
df <- df %>% 
  group_by(NACCID) %>%
  mutate(NRIV = sum(RIV, na.rm=TRUE))

df$RIV[(df$NRIV==0)]<- 0
df$RIV[(df$NRIV==1)]<- 9 
df$RIV[(df$NRIV>=2)]<- 1 
df <- as.data.frame(df)

#### TAC ###
df <- df %>% 
  group_by(NACCID) %>%
  mutate(NTAC = sum(TAC, na.rm=TRUE))

df$TAC[(df$NTAC==0)]<- 0
df$TAC[(df$NTAC==1)]<- 9 
df$TAC[(df$NTAC>=2)]<- 1 
df <- as.data.frame(df)

#### COM ###
df <- df %>% 
  group_by(NACCID) %>%
  mutate(NCOM = sum(COM, na.rm=TRUE))

df$COM[(df$NCOM==0)]<- 0
df$COM[(df$NCOM==1)]<- 9 
df$COM[(df$NCOM>=2)]<- 1 

file <- as.data.frame(df)
#write.csv(df, file = "MDalzcdr.csv",row.names=FALSE, na="")

################# Log reg #################################################
df <- as.data.frame(file) 

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
df <- df %>% replace_with_na(replace = list(DON = 9, MEM = 9, COM = 9, TAC = 9,
                                            GAL = 9, RIV = 9))


df.baseDHH <- df[,c(1:6,9,12:19,21:26)] 
df.baseDHH[ , 16:21 ][ df.baseDHH[ , 16:21 ] == 9 ] <- 1

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
lg.baseDHH$DON <- as.factor(as.character(lg.baseDHH$DON))
lg.baseDHH$MEM <- as.factor(as.character(lg.baseDHH$MEM))
lg.baseDHH$GAL <- as.factor(as.character(lg.baseDHH$GAL))
lg.baseDHH$RIV <- as.factor(as.character(lg.baseDHH$RIV))
lg.baseDHH$TAC <- as.factor(as.character(lg.baseDHH$TAC))
lg.baseDHH$COM <- as.factor(as.character(lg.baseDHH$COM))
lg.baseDHH$Diagnosis <- as.factor(as.character(lg.baseDHH$Diagnosis))
lg.baseDHH <- as.data.frame(lg.baseDHH)

##check column nos.
riskmodel <- fitLogRegModel(lg.baseDHH, cOutcome=7,
                            cNonGenPreds=c(2:6,8:19), cNonGenPredsCat=c(2,4:6,8:10,13:19),
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

don <- cbind(as.data.frame(prop.table(table(data.frame(lg.baseDHH$DON,
                                                       lg.baseDHH$Diagnosis)), 1)), as.data.frame(with(lg.baseDHH, table(DON, Diagnosis))))

mem <- cbind(as.data.frame(prop.table(table(data.frame(lg.baseDHH$MEM,
                                                       lg.baseDHH$Diagnosis)), 1)), as.data.frame(with(lg.baseDHH, table(MEM, Diagnosis))))

gal <- cbind(as.data.frame(prop.table(table(data.frame(lg.baseDHH$GAL,
                                                       lg.baseDHH$Diagnosis)), 1)), as.data.frame(with(lg.baseDHH, table(GAL, Diagnosis))))

riv <- cbind(as.data.frame(prop.table(table(data.frame(lg.baseDHH$RIV,
                                                        lg.baseDHH$Diagnosis)), 1)), as.data.frame(with(lg.baseDHH, table(RIV, Diagnosis))))

tac <- cbind(as.data.frame(prop.table(table(data.frame(lg.baseDHH$TAC,
                                                       lg.baseDHH$Diagnosis)), 1)), as.data.frame(with(lg.baseDHH, table(TAC, Diagnosis))))

com <- cbind(as.data.frame(prop.table(table(data.frame(lg.baseDHH$COM,
                                                       lg.baseDHH$Diagnosis)), 1)), as.data.frame(with(lg.baseDHH, table(COM, Diagnosis))))

bmi <- cbind(as.data.frame(prop.table(table(data.frame(lg.baseDHH$ABMI,
                                                       lg.baseDHH$Diagnosis)), 1)), as.data.frame(with(lg.baseDHH, table(ABMI, Diagnosis))))

cvd <- cbind(as.data.frame(prop.table(table(data.frame(lg.baseDHH$CVD,
                                                       lg.baseDHH$Diagnosis)), 1)), as.data.frame(with(lg.baseDHH, table(CVD, Diagnosis))))

dep <- cbind(as.data.frame(prop.table(table(data.frame(lg.baseDHH$DEPRSN,
                                                       lg.baseDHH$Diagnosis)), 1)), as.data.frame(with(lg.baseDHH, table(DEPRSN, Diagnosis))))

write.excel <- function(x,row.names=FALSE,col.names=TRUE,...) {
  write.table(x,"clipboard",sep="\t",row.names=row.names,col.names=col.names,...)}

write.excel(sex)
write.excel(race)
write.excel(alcohol)
write.excel(hear)
write.excel(don)
write.excel(mem)
write.excel(gal)
write.excel(riv)
write.excel(tac)
write.excel(com)
write.excel(hyper)
write.excel(dep)
write.excel(bmi)
write.excel(tbi)
write.excel(cvd)
