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

df <- df %>%
  group_by(NACCID) %>%
  mutate(check1=is.na(HYPERC1) & is.na(HYPERC2))
df <- df %>% group_by(NACCID) %>% filter(!all(check1=="TRUE"))
df <- df[,-(148)] 
df <- as.data.frame(df)

df <- df %>% 
  group_by(NACCID) %>%
  mutate(DCDC = sum(HYPERC1, na.rm=TRUE))
df$HYPERC1[(df$DCDC==0 & df$NEW_VISITNUM==1)]<- 0 
df$HYPERC1[(df$DCDC>=1 & df$NEW_VISITNUM==1)]<- 1 

df <- df %>% 
  group_by(NACCID) %>%
  mutate(DCDC = sum(HYPERC2, na.rm=TRUE))
df$HYPERC2[(df$DCDC==0 & df$NEW_VISITNUM==1)]<- 0 
df$HYPERC2[(df$DCDC>=1 & df$NEW_VISITNUM==1)]<- 1 

df <- df[,-(148)] 
df <- as.data.frame(df)

df$HYPERC[(df$HYPERC1==1 | df$HYPERC2==1)] <- 1
df$HYPERC[(df$HYPERC1==0 & df$HYPERC2==0)] <- 0

df <- as.data.frame(df)
df <- df[,c(1,8:10,14,26,49:88,121,126:137,139,148)]

## Find drug terms based on category
df$STAT <- apply(df, 1, function(x)as.integer(any(grep("LOVASTATIN|PRAVASTATIN|SIMVASTATIN|
            FLUVASTATIN|ATORVASTATIN|CERIVASTATIN|ROSUVASTATIN|PITAVASTATIN|RED YEAST RICE",x))))

df$NIAC <- apply(df, 1, function(x)as.integer(any(grep("niacin|niacinamide",x, ignore.case = TRUE))))

df$FIBR <- apply(df, 1, function(x)as.integer(any(grep("gemfibrozil|fenofibrate|
                                                       fenofibric acid",x, ignore.case = TRUE))))

df$BILE <- apply(df, 1, function(x)as.integer(any(grep("colestipol|cholestyramine",x, ignore.case = TRUE))))

df$CAI <- apply(df, 1, function(x)as.integer(any(grep("ezetimibe",x, ignore.case = TRUE))))
df$COMB <- apply(df, 1, function(x)as.integer(any(grep("lovastatin-niacin|aspirin-pravastatin|ezetimibe-simvastatin|
                                                       amlodipine-atorvastatin|niacin-simvastatin",x, ignore.case = TRUE))))

df <- df[,-(7:46)] 
df[df==-4]<-NA

################ Adjusting values ###################################

#### STAT ###
df <- df %>% 
  group_by(NACCID) %>%
  mutate(NSTAT = sum(STAT, na.rm=TRUE))

df$STAT[(df$NSTAT==0)]<- 0
df$STAT[(df$NSTAT==1)]<- 9 
df$STAT[(df$NSTAT>=2)]<- 1 
df <- as.data.frame(df)

#### NIAC ###
df <- df %>% 
  group_by(NACCID) %>%
  mutate(NNIAC = sum(NIAC, na.rm=TRUE))

df$NIAC[(df$NNIAC==0)]<- 0
df$NIAC[(df$NNIAC==1)]<- 9 
df$NIAC[(df$NNIAC>=2)]<- 1 
df <- as.data.frame(df)

#### FIBR ###
df <- df %>% 
  group_by(NACCID) %>%
  mutate(NFIBR = sum(FIBR, na.rm=TRUE))

df$FIBR[(df$NFIBR==0)]<- 0
df$FIBR[(df$NFIBR==1)]<- 9 
df$FIBR[(df$NFIBR>=2)]<- 1 
df <- as.data.frame(df)

#### BILE ###
df <- df %>% 
  group_by(NACCID) %>%
  mutate(NBILE = sum(BILE, na.rm=TRUE))

df$BILE[(df$NBILE==0)]<- 0
df$BILE[(df$NBILE==1)]<- 9 
df$BILE[(df$NBILE>=2)]<- 1 
df <- as.data.frame(df)

#### CAI ###
df <- df %>% 
  group_by(NACCID) %>%
  mutate(NCAI = sum(CAI, na.rm=TRUE))

df$CAI[(df$NCAI==0)]<- 0
df$CAI[(df$NCAI==1)]<- 9 
df$CAI[(df$NCAI>=2)]<- 1 
df <- as.data.frame(df)

#### COMB ###
df <- df %>% 
  group_by(NACCID) %>%
  mutate(NCOMB = sum(COMB, na.rm=TRUE))

df$COMB[(df$NCOMB==0)]<- 0
df$COMB[(df$NCOMB==1)]<- 9 
df$COMB[(df$NCOMB>=2)]<- 1 
file <- as.data.frame(df)

#write.csv(df, file = "MDlip.csv",row.names=FALSE, na="")

################# Log reg #################################################

df <- file
df$RACE[df$RACE==50] <- NA
df$SEX[df$SEX==2] <- 0  # female 0, male 1
df <- df[df$NEW_VISITNUM == 1, ] 
df <- subset(df, NACCAGE>=40,
             select=c(NACCID:NCOMB))
df <- as.data.frame(df)
df <- df %>% filter(SEX == 1)

df$Diagnosis[df$Diagnosis==1 & df$ALZP==1] <- 2
df <- df[df$Diagnosis != 1, ]
df$Diagnosis[df$Diagnosis==2] <- 1

## 9 means medication mentioned only once, decide if need to remove or not
df <- df %>% replace_with_na(replace = list(STAT = 9, NIAC = 9, FIBR = 9, BILE = 9,
                                           CAI = 9, COMB = 9))


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
lg.baseDHH$HYPERC <- as.factor(as.character(lg.baseDHH$HYPERC))
lg.baseDHH$STAT <- as.factor(as.character(lg.baseDHH$STAT))
lg.baseDHH$NIAC <- as.factor(as.character(lg.baseDHH$NIAC))
lg.baseDHH$FIBR <- as.factor(as.character(lg.baseDHH$FIBR))
lg.baseDHH$BILE <- as.factor(as.character(lg.baseDHH$BILE))
lg.baseDHH$CAI <- as.factor(as.character(lg.baseDHH$CAI))
lg.baseDHH$COMB <- as.factor(as.character(lg.baseDHH$COMB))
lg.baseDHH$Diagnosis <- as.factor(as.character(lg.baseDHH$Diagnosis))
lg.baseDHH <- as.data.frame(lg.baseDHH)

riskmodel <- fitLogRegModel(lg.baseDHH, cOutcome=7,
                            cNonGenPreds=c(2:6,8:22), cNonGenPredsCat=c(2,4:6,8:10,13:22),
                            cGenPreds=c(0), cGenPredsCat=c(0))

#gender
riskmodel <- fitLogRegModel(lg.baseDHH, cOutcome=7,
                            cNonGenPreds=c(3:6,8:22), cNonGenPredsCat=c(4:6,8:10,13:22),
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

hyperc <- cbind(as.data.frame(prop.table(table(data.frame(lg.baseDHH$HYPERC,
                                                        lg.baseDHH$Diagnosis)), 1)), as.data.frame(with(lg.baseDHH, table(HYPERC, Diagnosis))))

stat <- cbind(as.data.frame(prop.table(table(data.frame(lg.baseDHH$STAT,
                                                       lg.baseDHH$Diagnosis)), 1)), as.data.frame(with(lg.baseDHH, table(STAT, Diagnosis))))

niac <- cbind(as.data.frame(prop.table(table(data.frame(lg.baseDHH$NIAC,
                                                        lg.baseDHH$Diagnosis)), 1)), as.data.frame(with(lg.baseDHH, table(NIAC, Diagnosis))))

fib <- cbind(as.data.frame(prop.table(table(data.frame(lg.baseDHH$FIBR,
                                                       lg.baseDHH$Diagnosis)), 1)), as.data.frame(with(lg.baseDHH, table(FIBR, Diagnosis))))
 
bile <- cbind(as.data.frame(prop.table(table(data.frame(lg.baseDHH$BILE,
                                                       lg.baseDHH$Diagnosis)), 1)), as.data.frame(with(lg.baseDHH, table(BILE, Diagnosis))))

cai <- cbind(as.data.frame(prop.table(table(data.frame(lg.baseDHH$CAI,
                                                       lg.baseDHH$Diagnosis)), 1)), as.data.frame(with(lg.baseDHH, table(CAI, Diagnosis))))

comb <- cbind(as.data.frame(prop.table(table(data.frame(lg.baseDHH$COMB,
                                                       lg.baseDHH$Diagnosis)), 1)), as.data.frame(with(lg.baseDHH, table(COMB, Diagnosis))))

bmi <- cbind(as.data.frame(prop.table(table(data.frame(lg.baseDHH$ABMI,
                                                       lg.baseDHH$Diagnosis)), 1)), as.data.frame(with(lg.baseDHH, table(ABMI, Diagnosis))))

cvd <- cbind(as.data.frame(prop.table(table(data.frame(lg.baseDHH$CVD,
                                                       lg.baseDHH$Diagnosis)), 1)), as.data.frame(with(lg.baseDHH, table(CVD, Diagnosis))))

dep <- cbind(as.data.frame(prop.table(table(data.frame(lg.baseDHH$DEPRSN,
                                                       lg.baseDHH$Diagnosis)), 1)), as.data.frame(with(lg.baseDHH, table(DEPRSN, Diagnosis))))

xy <- cbind(sex,alcohol,hear,diab,hyperc,stat,niac,fib,bile,cai,comb,hyper,dep,tbi,cvd)

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
write.excel(hyperc)
write.excel(stat)
write.excel(niac)
write.excel(fib)
write.excel(bile)
write.excel(cai)
write.excel(comb)
write.excel(diab)
write.excel(hyper)
write.excel(dep)
write.excel(bmi)
write.excel(tbi)
write.excel(cvd)
