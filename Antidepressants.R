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
df <- df[,c(1,8:10,14,26,49:88,111,126:137,142)]

## Find drug terms based on category
df$SSRI <- apply(df, 1, function(x)as.integer(any(grep("FLUOXETINE|SERTRALINE|PAROXETINE|FLUVOXAMINE|
                                                       CITALOPRAM|ESCITALOPRAM|5-HYDROXYTRYPTOPHAN",x))))

df$SARI <- apply(df, 1, function(x)as.integer(any(grep("TRAZODONE|NEFAZODONE",x))))

df$SNRI <- apply(df, 1, function(x)as.integer(any(grep("VENLAFAXINE|DULOXETINE|MILNACIPRAN|PRISTIQ|
                                                       DESVENLAFAXINE|LEVOMILNACIPRAN",x))))

df$TCA <- apply(df, 1, function(x)as.integer(any(grep("NORTRIPTYLINE|DESIPRAMINE|AMITRIPTYLINE|DOXEPIN|
                                                      IMIPRAMINE|PROTRIPTYLINE|CLOMIPRAMINE",x))))

df$TECA <- apply(df, 1, function(x)as.integer(any(grep("AMOXAPINE|MAPROTILINE|MIRTAZAPINE",x))))
df$NDRI <- apply(df, 1, function(x)as.integer(any(grep("BUPROPION",x))))
df$SMS <- apply(df, 1, function(x)as.integer(any(grep("VILAZODONE",x))))
df$MAOB <- apply(df, 1, function(x)as.integer(any(grep("SELEGILINE",x))))
df$COMBO <- apply(df, 1, function(x)as.integer(any(grep("fluoxetine-|bupropion-|amitriptyline-",x, ignore.case = TRUE))))

df <- df[,-(7:46)] 
df[df==-4]<-NA

################ Adjusting values ###################################

#### SSRI ###
df <- df %>% 
  group_by(NACCID) %>%
  mutate(NSSRI = sum(SSRI, na.rm=TRUE))

df$SSRI[(df$NSSRI==0)]<- 0
df$SSRI[(df$NSSRI==1)]<- 9 
df$SSRI[(df$NSSRI>=2)]<- 1 
df <- as.data.frame(df)

#### SARI ###
df <- df %>% 
  group_by(NACCID) %>%
  mutate(NSARI = sum(SARI, na.rm=TRUE))

df$SARI[(df$NSARI==0)]<- 0
df$SARI[(df$NSARI==1)]<- 9 
df$SARI[(df$NSARI>=2)]<- 1 
df <- as.data.frame(df)

#### SNRI ###
df <- df %>% 
  group_by(NACCID) %>%
  mutate(NSNRI = sum(SNRI, na.rm=TRUE))

df$SNRI[(df$NSNRI==0)]<- 0
df$SNRI[(df$NSNRI==1)]<- 9 
df$SNRI[(df$NSNRI>=2)]<- 1 
df <- as.data.frame(df)

#### TCA ###
df <- df %>% 
  group_by(NACCID) %>%
  mutate(NTCA = sum(TCA, na.rm=TRUE))

df$TCA[(df$NTCA==0)]<- 0
df$TCA[(df$NTCA==1)]<- 9 
df$TCA[(df$NTCA>=2)]<- 1 
df <- as.data.frame(df)

#### TECA ###
df <- df %>% 
  group_by(NACCID) %>%
  mutate(NTECA = sum(TECA, na.rm=TRUE))

df$TECA[(df$NTECA==0)]<- 0
df$TECA[(df$NTECA==1)]<- 9 
df$TECA[(df$NTECA>=2)]<- 1 
df <- as.data.frame(df)

#### NDRI ###
df <- df %>% 
  group_by(NACCID) %>%
  mutate(NNDRI = sum(NDRI, na.rm=TRUE))

df$NDRI[(df$NNDRI==0)]<- 0
df$NDRI[(df$NNDRI==1)]<- 9 
df$NDRI[(df$NNDRI>=2)]<- 1 
df <- as.data.frame(df)

#### SMS ###
df <- df %>% 
  group_by(NACCID) %>%
  mutate(NSMS = sum(SMS, na.rm=TRUE))

df$SMS[(df$NSMS==0)]<- 0
df$SMS[(df$NSMS==1)]<- 9 
df$SMS[(df$NSMS>=2)]<- 1 
df <- as.data.frame(df)

#### MAOB ###
df <- df %>% 
  group_by(NACCID) %>%
  mutate(NMAOB = sum(MAOB, na.rm=TRUE))

df$MAOB[(df$NMAOB==0)]<- 0
df$MAOB[(df$NMAOB==1)]<- 9 
df$MAOB[(df$NMAOB>=2)]<- 1 
df <- as.data.frame(df)

#### COMBO ###
df <- df %>% 
  group_by(NACCID) %>%
  mutate(NCOMBO = sum(COMBO, na.rm=TRUE))

df$COMBO[(df$NCOMBO==0)]<- 0
df$COMBO[(df$NCOMBO==1)]<- 9 
df$COMBO[(df$NCOMBO>=2)]<- 1 
file <- as.data.frame(df)

#write.csv(df, file = "NMadepcdr.csv",row.names=FALSE, na="")

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
df <- df %>% replace_with_na(replace = list(SSRI = 9, SARI = 9, SNRI = 9, TCA = 9,
                                            TECA = 9, NDRI = 9, SMS = 9, MAOB = 9, COMBO = 9))


df.baseDHH <- df[,c(1:6,9,12:19,21:29)] 
df.baseDHH[ , 16:24][ df.baseDHH[ , 16:24 ] == 9 ] <- 1

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
lg.baseDHH$SSRI <- as.factor(as.character(lg.baseDHH$SSRI))
lg.baseDHH$SARI <- as.factor(as.character(lg.baseDHH$SARI))
lg.baseDHH$SNRI <- as.factor(as.character(lg.baseDHH$SNRI))
lg.baseDHH$TCA <- as.factor(as.character(lg.baseDHH$TCA))
lg.baseDHH$TECA <- as.factor(as.character(lg.baseDHH$TECA))
lg.baseDHH$NDRI <- as.factor(as.character(lg.baseDHH$NDRI))
lg.baseDHH$SMS <- as.factor(as.character(lg.baseDHH$SMS))
lg.baseDHH$MAOB <- as.factor(as.character(lg.baseDHH$MAOB ))
lg.baseDHH$COMBO <- as.factor(as.character(lg.baseDHH$COMBO ))
lg.baseDHH$Diagnosis <- as.factor(as.character(lg.baseDHH$Diagnosis))
lg.baseDHH <- as.data.frame(lg.baseDHH)

###check column levels

riskmodel <- fitLogRegModel(lg.baseDHH, cOutcome=7,
                            cNonGenPreds=c(2:6,8:24), cNonGenPredsCat=c(2,4:6,8:10,13:24),
                            cGenPreds=c(0), cGenPredsCat=c(0))

#gender
riskmodel <- fitLogRegModel(lg.baseDHH, cOutcome=7,
                            cNonGenPreds=c(3:6,8:21,23,24), cNonGenPredsCat=c(4:6,8:10,13:21,23,24),
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

ssr <- cbind(as.data.frame(prop.table(table(data.frame(lg.baseDHH$SSRI,
     lg.baseDHH$Diagnosis)), 1)), as.data.frame(with(lg.baseDHH, table(SSRI, Diagnosis))))

sar <- cbind(as.data.frame(prop.table(table(data.frame(lg.baseDHH$SARI,
   lg.baseDHH$Diagnosis)), 1)), as.data.frame(with(lg.baseDHH, table(SARI, Diagnosis))))

snr <- cbind(as.data.frame(prop.table(table(data.frame(lg.baseDHH$SNRI,
     lg.baseDHH$Diagnosis)), 1)), as.data.frame(with(lg.baseDHH, table(SNRI, Diagnosis))))

tca <- cbind(as.data.frame(prop.table(table(data.frame(lg.baseDHH$TCA,
  lg.baseDHH$Diagnosis)), 1)), as.data.frame(with(lg.baseDHH, table(TCA, Diagnosis))))

teca <- cbind(as.data.frame(prop.table(table(data.frame(lg.baseDHH$TECA,
     lg.baseDHH$Diagnosis)), 1)), as.data.frame(with(lg.baseDHH, table(TECA, Diagnosis))))

ndr <- cbind(as.data.frame(prop.table(table(data.frame(lg.baseDHH$NDRI,
   lg.baseDHH$Diagnosis)), 1)), as.data.frame(with(lg.baseDHH, table(NDRI, Diagnosis))))

sms <- cbind(as.data.frame(prop.table(table(data.frame(lg.baseDHH$SMS,
  lg.baseDHH$Diagnosis)), 1)), as.data.frame(with(lg.baseDHH, table(SMS, Diagnosis))))

mao <- cbind(as.data.frame(prop.table(table(data.frame(lg.baseDHH$MAOB,
      lg.baseDHH$Diagnosis)), 1)), as.data.frame(with(lg.baseDHH, table(MAOB, Diagnosis))))

com <- cbind(as.data.frame(prop.table(table(data.frame(lg.baseDHH$COMBO,
    lg.baseDHH$Diagnosis)), 1)), as.data.frame(with(lg.baseDHH, table(COMBO, Diagnosis))))

bmi <- cbind(as.data.frame(prop.table(table(data.frame(lg.baseDHH$ABMI,
                                                       lg.baseDHH$Diagnosis)), 1)), as.data.frame(with(lg.baseDHH, table(ABMI, Diagnosis))))

cvd <- cbind(as.data.frame(prop.table(table(data.frame(lg.baseDHH$CVD,
                                                       lg.baseDHH$Diagnosis)), 1)), as.data.frame(with(lg.baseDHH, table(CVD, Diagnosis))))

dep <- cbind(as.data.frame(prop.table(table(data.frame(lg.baseDHH$DEPRSN,
                                                       lg.baseDHH$Diagnosis)), 1)), as.data.frame(with(lg.baseDHH, table(DEPRSN, Diagnosis))))

xy <- cbind(sex,alcohol,hear,diab,ssr,sar,snr,tca,teca,ndr,sms,mao,com,hyper,dep,tbi,cvd)

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
write.excel(ssr)
write.excel(sar)
write.excel(snr)
write.excel(tca)
write.excel(teca)
write.excel(ndr)
write.excel(sms)
write.excel(mao)
write.excel(com)
write.excel(diab)
write.excel(hyper) 
write.excel(dep)
write.excel(bmi)
write.excel(tbi)
write.excel(cvd)
