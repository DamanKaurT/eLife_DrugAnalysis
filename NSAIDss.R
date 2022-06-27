###### Searching drug terms #####

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
#df <- df %>%
 # group_by(NACCID) %>%
  #mutate(check1=is.na(ARTH1) & is.na(ARTH2))
#df <- df %>% group_by(NACCID) %>% filter(!all(check1=="TRUE"))
#df <- df[,-(148)] 
#df <- as.data.frame(df)

#df <- df %>% 
 # group_by(NACCID) %>%
  #mutate(DCDC = sum(ARTH1, na.rm=TRUE))
#df$ARTH1[(df$DCDC==0 & df$NEW_VISITNUM==1)]<- 0 
#df$ARTH1[(df$DCDC>=1 & df$NEW_VISITNUM==1)]<- 1 

#df <- df %>% 
 # group_by(NACCID) %>%
  #mutate(DCDC = sum(ARTH2, na.rm=TRUE))
#df$ARTH2[(df$DCDC==0 & df$NEW_VISITNUM==1)]<- 0 
#df$ARTH2[(df$DCDC>=1 & df$NEW_VISITNUM==1)]<- 1 

#df <- df[,-(148)] 
#df <- as.data.frame(df)

#df$ARTHRIT[(df$ARTH1==1 | df$ARTH2==1)] <- 1
#df$ARTHRIT[(df$ARTH1==0 & df$ARTH2==0)] <- 0

#df <- as.data.frame(df)
#df <- df[,c(1,8:10,14,26,49:88,122,126:137,140,148)]

#### Find drug terms based on category #####

df <- as.data.frame(df)
df <- df[,c(1,8:10,14,26,49:88,122,126:137,140)]

df$SALC <- apply(df, 1, function(x)as.integer(any(grep("aspirin|asa/|diflunisal|choline sal|salsalate|
                                        sodium sal|magnesium sal|methyl sal|-salicylic|/salicylic|
                                        salicylic acid topical|salicylic acid-",x, ignore.case = TRUE))))

df$PROP <- apply(df, 1, function(x)as.integer(any(grep("IBUPROFEN|NAPROXEN|FENOPROFEN|KETOPROFEN|
                                                       FLURBIPROFEN|OXAPROZIN",x, ignore.case = TRUE))))

df$OXI <- apply(df, 1, function(x)as.integer(any(grep("oxicam",x, ignore.case = TRUE))))

df$COX <- apply(df, 1, function(x)as.integer(any(grep("coxib",x, ignore.case = TRUE))))

df$ACD <- apply(df, 1, function(x)as.integer(any(grep("SULINDAC|TOLMETIN|bromfenac|diclofenac|
                                indomethacin|etodolac|ketorolac|nabumetone",x, ignore.case = TRUE))))

df <- df[,-(7:46)] 
df[df==-4]<-NA

################ Adjusting values ###################################

#### SALC ###
df <- df %>% 
  group_by(NACCID) %>%
  mutate(NSALC = sum(SALC, na.rm=TRUE))

df$SALC[(df$NSALC==0)]<- 0
df$SALC[(df$NSALC==1)]<- 9 
df$SALC[(df$NSALC>=2)]<- 1 
df <- as.data.frame(df)

#### PROP ###
df <- df %>% 
  group_by(NACCID) %>%
  mutate(NPROP = sum(PROP, na.rm=TRUE))

df$PROP[(df$NPROP==0)]<- 0
df$PROP[(df$NPROP==1)]<- 9 
df$PROP[(df$NPROP>=2)]<- 1 
df <- as.data.frame(df)

#### OXI ###
df <- df %>% 
  group_by(NACCID) %>%
  mutate(NOXI = sum(OXI, na.rm=TRUE))

df$OXI[(df$NOXI==0)]<- 0
df$OXI[(df$NOXI==1)]<- 9 
df$OXI[(df$NOXI>=2)]<- 1 
df <- as.data.frame(df)

#### COX ###
df <- df %>% 
  group_by(NACCID) %>%
  mutate(NCOX = sum(COX, na.rm=TRUE))

df$COX[(df$NCOX==0)]<- 0
df$COX[(df$NCOX==1)]<- 9 
df$COX[(df$NCOX>=2)]<- 1 
df <- as.data.frame(df)

#### ACD ###
df <- df %>% 
  group_by(NACCID) %>%
  mutate(NACD = sum(ACD, na.rm=TRUE))

df$ACD[(df$NACD==0)]<- 0
df$ACD[(df$NACD==1)]<- 9 
df$ACD[(df$NACD>=2)]<- 1 
file <- as.data.frame(df)

#write.csv(df, file = "NMnsd.csv",row.names=FALSE, na="")

################# Log reg #################################################

df <- file
df$RACE[df$RACE==50] <- NA
df$SEX[df$SEX==2] <- 0  # female 0, male 1
df <- df[df$NEW_VISITNUM == 1, ] 
df <- subset(df, NACCAGE>=40,
             select=c(NACCID:NACD))
df <- as.data.frame(df)
df <- df %>% filter(SEX == 1)

df$Diagnosis[df$Diagnosis==1 & df$ALZP==1] <- 2
df <- df[df$Diagnosis != 1, ]
df$Diagnosis[df$Diagnosis==2] <- 1

## 9 means medication mentioned only once, decide if need to remove or not
df <- df %>% replace_with_na(replace = list(SALC = 9, PROP = 9, OXI = 9, COX = 9,
                                           ACD = 9))


df.baseDHH <- df[,c(1:6,9,12:19,21:25)] 
df.baseDHH[ , 16:20 ][ df.baseDHH[ , 16:20 ] == 9 ] <- 1

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
#lg.baseDHH$ARTHRIT <- as.factor(as.character(lg.baseDHH$ARTHRIT))
lg.baseDHH$SALC <- as.factor(as.character(lg.baseDHH$SALC))
lg.baseDHH$PROP <- as.factor(as.character(lg.baseDHH$PROP))
lg.baseDHH$OXI <- as.factor(as.character(lg.baseDHH$OXI))
lg.baseDHH$COX <- as.factor(as.character(lg.baseDHH$COX))
lg.baseDHH$ACD <- as.factor(as.character(lg.baseDHH$ACD))
lg.baseDHH$Diagnosis <- as.factor(as.character(lg.baseDHH$Diagnosis))
lg.baseDHH <- as.data.frame(lg.baseDHH)

riskmodel <- fitLogRegModel(lg.baseDHH, cOutcome=7,
                           cNonGenPreds=c(2:6,8:20), cNonGenPredsCat=c(2,4:6,8:10,13:20),
                          cGenPreds=c(0), cGenPredsCat=c(0))

#gender
riskmodel <- fitLogRegModel(lg.baseDHH, cOutcome=7,
                           cNonGenPreds=c(3:6,8:20), cNonGenPredsCat=c(4:6,8:10,13:20),
                          cGenPreds=c(0), cGenPredsCat=c(0))
#riskmodel <- fitLogRegModel(lg.baseDHH, cOutcome=7,
 #                          cNonGenPreds=c(2:6,8:21), cNonGenPredsCat=c(2,4:6,8:10,13:21),
#                          cGenPreds=c(0), cGenPredsCat=c(0))

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

arthrit <- cbind(as.data.frame(prop.table(table(data.frame(lg.baseDHH$ARTHRIT,
  lg.baseDHH$Diagnosis)), 1)), as.data.frame(with(lg.baseDHH, table(ARTHRIT, Diagnosis))))

sal <- cbind(as.data.frame(prop.table(table(data.frame(lg.baseDHH$SALC,
    lg.baseDHH$Diagnosis)), 1)), as.data.frame(with(lg.baseDHH, table(SALC, Diagnosis))))

pro <- cbind(as.data.frame(prop.table(table(data.frame(lg.baseDHH$PROP,
  lg.baseDHH$Diagnosis)), 1)), as.data.frame(with(lg.baseDHH, table(PROP, Diagnosis))))

oxi <- cbind(as.data.frame(prop.table(table(data.frame(lg.baseDHH$OXI,
  lg.baseDHH$Diagnosis)), 1)), as.data.frame(with(lg.baseDHH, table(OXI, Diagnosis))))

cox <- cbind(as.data.frame(prop.table(table(data.frame(lg.baseDHH$COX,
  lg.baseDHH$Diagnosis)), 1)), as.data.frame(with(lg.baseDHH, table(COX, Diagnosis))))

acd <- cbind(as.data.frame(prop.table(table(data.frame(lg.baseDHH$ACD,
lg.baseDHH$Diagnosis)), 1)), as.data.frame(with(lg.baseDHH, table(ACD, Diagnosis))))

bmi <- cbind(as.data.frame(prop.table(table(data.frame(lg.baseDHH$ABMI,
  lg.baseDHH$Diagnosis)), 1)), as.data.frame(with(lg.baseDHH, table(ABMI, Diagnosis))))

cvd <- cbind(as.data.frame(prop.table(table(data.frame(lg.baseDHH$CVD,
lg.baseDHH$Diagnosis)), 1)), as.data.frame(with(lg.baseDHH, table(CVD, Diagnosis))))

dep <- cbind(as.data.frame(prop.table(table(data.frame(lg.baseDHH$DEPRSN,
  lg.baseDHH$Diagnosis)), 1)), as.data.frame(with(lg.baseDHH, table(DEPRSN, Diagnosis))))

xy <- cbind(sex,alcohol,hear,diab,sal,pro,oxi,cox,acd,hyper,dep,tbi,cvd)

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
write.excel(arthrit)
write.excel(sal)
write.excel(pro)
write.excel(oxi)
write.excel(cox)
write.excel(acd)
write.excel(diab)
write.excel(hyper)
write.excel(dep)
write.excel(bmi)
write.excel(tbi)
write.excel(cvd)


###### arthritis #####
df <- df %>%
  group_by(NACCID) %>%
  mutate(check1=is.na(ARTH1) & is.na(ARTH2))
df <- df %>% group_by(NACCID) %>% filter(!all(check1=="TRUE"))
df <- df[,-(148)] 
df <- as.data.frame(df)

df <- df %>% 
  group_by(NACCID) %>%
  mutate(DCDC = sum(ARTH1, na.rm=TRUE))
df$ARTH1[(df$DCDC==0 & df$NEW_VISITNUM==1)]<- 0 
df$ARTH1[(df$DCDC>=1 & df$NEW_VISITNUM==1)]<- 1 

df <- df %>% 
  group_by(NACCID) %>%
  mutate(DCDC = sum(ARTH2, na.rm=TRUE))
df$ARTH2[(df$DCDC==0 & df$NEW_VISITNUM==1)]<- 0 
df$ARTH2[(df$DCDC>=1 & df$NEW_VISITNUM==1)]<- 1 

df <- df[,-(148)] 
df <- as.data.frame(df)

df$ARTHRIT[(df$ARTH1==1 | df$ARTH2==1)] <- 1
df$ARTHRIT[(df$ARTH1==0 & df$ARTH2==0)] <- 0

df <- as.data.frame(df)
df <- df[,c(1,8:10,14,26,49:88,122,126:137,140,148)]
