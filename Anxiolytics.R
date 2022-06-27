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
  mutate(check1=is.na(ANX1) & is.na(ANX2) & is.na(ANX3))
df <- df %>% group_by(NACCID) %>% filter(!all(check1=="TRUE"))
df <- df[,-(148)] 
df <- as.data.frame(df)

df <- df %>% 
  group_by(NACCID) %>%
  mutate(DCDC = sum(ANX1, na.rm=TRUE))
df$ANX1[(df$DCDC==0 & df$NEW_VISITNUM==1)]<- 0 
df$ANX1[(df$DCDC>=1 & df$NEW_VISITNUM==1)]<- 1 

df <- df %>% 
  group_by(NACCID) %>%
  mutate(DCDC = sum(ANX2, na.rm=TRUE))
df$ANX2[(df$DCDC==0 & df$NEW_VISITNUM==1)]<- 0 
df$ANX2[(df$DCDC>=1 & df$NEW_VISITNUM==1)]<- 1 

df <- df %>% 
  group_by(NACCID) %>%
  mutate(DCDC = sum(ANX3, na.rm=TRUE))
df$ANX3[(df$DCDC==0 & df$NEW_VISITNUM==1)]<- 0 
df$ANX3[(df$DCDC>=1 & df$NEW_VISITNUM==1)]<- 1 

df <- df[,-(148)] 
df <- as.data.frame(df)

df$ANX[(df$ANX1==1 | df$ANX2==1 | df$ANX3==1)] <- 1
df$ANX[(df$ANX1==0 & df$ANX2==0 & df$ANX3==0)] <- 0

df <- as.data.frame(df)
df <- df[,c(1,8:10,14,26,49:88,108,126:137,144,148)]

## Find drug terms based on category
df$BARB <- apply(df, 1, function(x)as.integer(any(grep("secobarbital|butalbital|amobarbital|
                    pentobarbital|mephobarbital|butabarbital|phenobarbital|thiopental",x, ignore.case = TRUE))))

df$BENZ <- apply(df, 1, function(x)as.integer(any(grep("oxazepam|diazepam|lorazepam|alprazolam|
                          chlordiazepoxide|clonazepam|flurazepam|temazepam|triazolam|
                                            halazepam|estazolam|quazepam|clobazam",x, ignore.case = TRUE))))

df$ATH <- apply(df, 1, function(x)as.integer(any(grep("diphenhydramine|pyrilamine|hydroxyzine|
                                                  doxylamine|promethazine",x, ignore.case = TRUE))))

df$HYP <- apply(df, 1, function(x)as.integer(any(grep("chloral hydrate|zolpidem|melatonin|
                          zaleplon|eszopiclone|ramelteon|sodium oxybate|
                                  acetylcarbromal",x, ignore.case = TRUE))))

df$CARB <- apply(df, 1, function(x)as.integer(any(grep("meprobamate",x, ignore.case = TRUE))))
df$DOX <- apply(df, 1, function(x)as.integer(any(grep("doxepin",x, ignore.case = TRUE))))
df$AZA <- apply(df, 1, function(x)as.integer(any(grep("buspirone",x, ignore.case = TRUE))))

df <- df[,-(7:46)] 
df[df==-4]<-NA

################ Adjusting anx values ###################################

#### BARB ###
df <- df %>% 
  group_by(NACCID) %>%
  mutate(NBARB = sum(BARB, na.rm=TRUE))

df$BARB[(df$NBARB==0)]<- 0
df$BARB[(df$NBARB==1)]<- 9 
df$BARB[(df$NBARB>=2)]<- 1 
df <- as.data.frame(df)

#### BENZ ###
df <- df %>% 
  group_by(NACCID) %>%
  mutate(NBENZ = sum(BENZ, na.rm=TRUE))

df$BENZ[(df$NBENZ==0)]<- 0
df$BENZ[(df$NBENZ==1)]<- 9 
df$BENZ[(df$NBENZ>=2)]<- 1 
df <- as.data.frame(df)

#### ATH ###
df <- df %>% 
  group_by(NACCID) %>%
  mutate(NATH = sum(ATH, na.rm=TRUE))

df$ATH[(df$NATH==0)]<- 0
df$ATH[(df$NATH==1)]<- 9 
df$ATH[(df$NATH>=2)]<- 1 
df <- as.data.frame(df)

#### HYP ###
df <- df %>% 
  group_by(NACCID) %>%
  mutate(NHYP = sum(HYP, na.rm=TRUE))

df$HYP[(df$NHYP==0)]<- 0
df$HYP[(df$NHYP==1)]<- 9 
df$HYP[(df$NHYP>=2)]<- 1 
df <- as.data.frame(df)

#### CARB ###
df <- df %>% 
  group_by(NACCID) %>%
  mutate(NCARB = sum(CARB, na.rm=TRUE))

df$CARB[(df$NCARB==0)]<- 0
df$CARB[(df$NCARB==1)]<- 9 
df$CARB[(df$NCARB>=2)]<- 1 
df <- as.data.frame(df)

#### DOX ###
df <- df %>% 
  group_by(NACCID) %>%
  mutate(NDOX = sum(DOX, na.rm=TRUE))

df$DOX[(df$NDOX==0)]<- 0
df$DOX[(df$NDOX==1)]<- 9 
df$DOX[(df$NDOX>=2)]<- 1 
df <- as.data.frame(df)

#### AZA ###
df <- df %>% 
  group_by(NACCID) %>%
  mutate(NAZA = sum(AZA, na.rm=TRUE))

df$AZA[(df$NAZA==0)]<- 0
df$AZA[(df$NAZA==1)]<- 9 
df$AZA[(df$NAZA>=2)]<- 1 
file <- as.data.frame(df)

#write.csv(df, file = "MDanxcdr.csv",row.names=FALSE, na="")

################# Log reg #################################################

df <- file
df$RACE[df$RACE==50] <- NA
df$SEX[df$SEX==2] <- 0  # female 0, male 1
df <- df[df$NEW_VISITNUM == 1, ] 
df <- subset(df, NACCAGE>=40,
             select=c(NACCID:NAZA))
df <- as.data.frame(df)
df <- df %>% filter(SEX == 1)

df$Diagnosis[df$Diagnosis==1 & df$ALZP==1] <- 2
df <- df[df$Diagnosis != 1, ]
df$Diagnosis[df$Diagnosis==2] <- 1


## 9 means medication mentioned only once, decide if need to remove or not
df <- df %>% replace_with_na(replace = list(BARB = 9, BENZ = 9, ATH = 9, HYP = 9,
                                            CARB = 9, DOX = 9, AZA = 9))


df.baseDHH <- df[,c(1:6,9,12:19,21:28)] 
df.baseDHH[ , 17:23 ][ df.baseDHH[ , 17:23 ] == 9 ] <- 1

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
lg.baseDHH$ANX <- as.factor(as.character(lg.baseDHH$ANX))
lg.baseDHH$BARB <- as.factor(as.character(lg.baseDHH$BARB))
lg.baseDHH$BENZ <- as.factor(as.character(lg.baseDHH$BENZ))
lg.baseDHH$ATH <- as.factor(as.character(lg.baseDHH$ATH))
lg.baseDHH$HYP <- as.factor(as.character(lg.baseDHH$HYP))
lg.baseDHH$CARB <- as.factor(as.character(lg.baseDHH$CARB))
lg.baseDHH$DOX <- as.factor(as.character(lg.baseDHH$DOX))
lg.baseDHH$AZA <- as.factor(as.character(lg.baseDHH$AZA))
lg.baseDHH$Diagnosis <- as.factor(as.character(lg.baseDHH$Diagnosis))
lg.baseDHH <- as.data.frame(lg.baseDHH)

riskmodel <- fitLogRegModel(lg.baseDHH, cOutcome=7,
                            cNonGenPreds=c(2:6,8:20,22,23), cNonGenPredsCat=c(2,4:6,8:10,13:20,22,23),
                            cGenPreds=c(0), cGenPredsCat=c(0))

#gender
riskmodel <- fitLogRegModel(lg.baseDHH, cOutcome=7,
                            cNonGenPreds=c(3:6,8:20,22,23), cNonGenPredsCat=c(4:6,8:10,13:20,22,23),
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

anx <- cbind(as.data.frame(prop.table(table(data.frame(lg.baseDHH$ANX,
                                                           lg.baseDHH$Diagnosis)), 1)), as.data.frame(with(lg.baseDHH, table(ANX, Diagnosis))))

bar <- cbind(as.data.frame(prop.table(table(data.frame(lg.baseDHH$BARB,
                                                       lg.baseDHH$Diagnosis)), 1)), as.data.frame(with(lg.baseDHH, table(BARB, Diagnosis))))

ben <- cbind(as.data.frame(prop.table(table(data.frame(lg.baseDHH$BENZ,
                                                       lg.baseDHH$Diagnosis)), 1)), as.data.frame(with(lg.baseDHH, table(BENZ, Diagnosis))))

ath <- cbind(as.data.frame(prop.table(table(data.frame(lg.baseDHH$ATH,
                                                       lg.baseDHH$Diagnosis)), 1)), as.data.frame(with(lg.baseDHH, table(ATH, Diagnosis))))

hyp <- cbind(as.data.frame(prop.table(table(data.frame(lg.baseDHH$HYP,
                                                       lg.baseDHH$Diagnosis)), 1)), as.data.frame(with(lg.baseDHH, table(HYP, Diagnosis))))

car <- cbind(as.data.frame(prop.table(table(data.frame(lg.baseDHH$CARB,
                                                       lg.baseDHH$Diagnosis)), 1)), as.data.frame(with(lg.baseDHH, table(CARB, Diagnosis))))

dox <- cbind(as.data.frame(prop.table(table(data.frame(lg.baseDHH$DOX,
                                                       lg.baseDHH$Diagnosis)), 1)), as.data.frame(with(lg.baseDHH, table(DOX, Diagnosis))))

aza <- cbind(as.data.frame(prop.table(table(data.frame(lg.baseDHH$AZA,
                                                       lg.baseDHH$Diagnosis)), 1)), as.data.frame(with(lg.baseDHH, table(AZA, Diagnosis))))

bmi <- cbind(as.data.frame(prop.table(table(data.frame(lg.baseDHH$ABMI,
                                                       lg.baseDHH$Diagnosis)), 1)), as.data.frame(with(lg.baseDHH, table(ABMI, Diagnosis))))

cvd <- cbind(as.data.frame(prop.table(table(data.frame(lg.baseDHH$CVD,
                                                       lg.baseDHH$Diagnosis)), 1)), as.data.frame(with(lg.baseDHH, table(CVD, Diagnosis))))

dep <- cbind(as.data.frame(prop.table(table(data.frame(lg.baseDHH$DEPRSN,
                                                       lg.baseDHH$Diagnosis)), 1)), as.data.frame(with(lg.baseDHH, table(DEPRSN, Diagnosis))))

xy <- cbind(sex,alcohol,hear,diab,anx,bar,ben,ath,hyp,car,dox,aza,hyper,dep,tbi,cvd)

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
write.excel(anx)
write.excel(bar)
write.excel(ben)
write.excel(ath)
write.excel(hyp)
write.excel(car)
write.excel(dox)
write.excel(aza)
write.excel(diab)
write.excel(hyper)
write.excel(dep)
write.excel(bmi)
write.excel(tbi)
write.excel(cvd)
