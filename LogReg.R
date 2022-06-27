################## All meds ####################

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

df$SEX[df$SEX==2] <- 0  # female 0, male 1
df$RACE[df$RACE==50] <- NA
df <- df[df$NEW_VISITNUM == 1, ] 
df <- subset(df, NACCAGE>=40,
             select=c(NACCID:NPDM))
df <- as.data.frame(df)
df <- df %>% filter(SEX == 1)

# for analyzing only AD cases
df$Diagnosis[df$Diagnosis==1 & df$ALZP==1] <- 2
df <- df[df$Diagnosis != 1, ]
df$Diagnosis[df$Diagnosis==2] <- 1

## 9 means medication mentioned only once, decide if need to remove or not
df <- df %>% 
  replace_with_na(replace = list(NACCAANX = 9, NACCAC = 9, NACCADEP = 9, NACCADMD = 9,
                                 NACCAHTN = 9, NACCAPSY = 9, NACCDBMD = 9, NACCLIPL = 9, 
                                 NACCNSD = 9, NACCPDMD = 9))


df.baseDHH <- df[,c(1,8:10,14,26,108,109,111:113,115,118,121:123,127,130:137)]
df.baseDHH[ , 7:16 ][ df.baseDHH[ , 7:16 ] == 9 ] <- 1

df.baseDHH.compl <- df.baseDHH[complete.cases(df.baseDHH), ]
lg.baseDHH <- df.baseDHH.compl

lg.baseDHH$DIAB1[(lg.baseDHH$DIAB==0 & lg.baseDHH$NACCDBMD==0)]<- 0 
lg.baseDHH$DIAB1[(lg.baseDHH$DIAB==0 & lg.baseDHH$NACCDBMD==1)]<- 1
lg.baseDHH$DIAB1[(lg.baseDHH$DIAB==1 & lg.baseDHH$NACCDBMD==0)]<- 2 
lg.baseDHH$DIAB1[(lg.baseDHH$DIAB==1 & lg.baseDHH$NACCDBMD==1)]<- 3 

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
lg.baseDHH$NACCAHTN <- as.factor(as.character(lg.baseDHH$NACCAHTN))
lg.baseDHH$NACCLIPL <- as.factor(as.character(lg.baseDHH$NACCLIPL))
lg.baseDHH$NACCNSD <- as.factor(as.character(lg.baseDHH$NACCNSD))
lg.baseDHH$NACCAC <- as.factor(as.character(lg.baseDHH$NACCAC))
lg.baseDHH$NACCADEP <- as.factor(as.character(lg.baseDHH$NACCADEP))
lg.baseDHH$NACCAPSY <- as.factor(as.character(lg.baseDHH$NACCAPSY))
lg.baseDHH$NACCAANX <- as.factor(as.character(lg.baseDHH$NACCAANX))
lg.baseDHH$NACCDBMD <- as.factor(as.character(lg.baseDHH$NACCDBMD))
lg.baseDHH$NACCADMD <- as.factor(as.character(lg.baseDHH$NACCADMD))
lg.baseDHH$NACCPDMD <- as.factor(as.character(lg.baseDHH$NACCPDMD))
lg.baseDHH$Diagnosis <- as.factor(as.character(lg.baseDHH$Diagnosis))
lg.baseDHH <- as.data.frame(lg.baseDHH)

riskmodel <- fitLogRegModel(lg.baseDHH, cOutcome=17,
                            cNonGenPreds=c(2:9,11,12,14:16,19:26), cNonGenPredsCat=c(2,4:9,11,12,14:16,19:20,23:26),
                            cGenPreds=c(0), cGenPredsCat=c(0))
#gender
riskmodel <- fitLogRegModel(lg.baseDHH, cOutcome=17,
                            cNonGenPreds=c(3:9,11:16,18:25), cNonGenPredsCat=c(4:9,11:16,18:20,23:25),
                            cGenPreds=c(0), cGenPredsCat=c(0))

summary(riskmodel)
oddr <- ORmultivariate(riskmodel,filename="medtest.txt")

vif <- vif(riskmodel)
stat <- cbind(rownames(vif), data.frame(vif, row.names=NULL))
write.csv(stat, file = "vif.csv",row.names=FALSE, na="")


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

aanx <- cbind(as.data.frame(prop.table(table(data.frame(lg.baseDHH$NACCAANX,
lg.baseDHH$Diagnosis)), 1)), as.data.frame(with(lg.baseDHH, table(NACCAANX, Diagnosis))))

antc <- cbind(as.data.frame(prop.table(table(data.frame(lg.baseDHH$NACCAC,
lg.baseDHH$Diagnosis)), 1)), as.data.frame(with(lg.baseDHH, table(NACCAC, Diagnosis))))

adep <- cbind(as.data.frame(prop.table(table(data.frame(lg.baseDHH$NACCADEP,
lg.baseDHH$Diagnosis)), 1)), as.data.frame(with(lg.baseDHH, table(NACCADEP, Diagnosis))))

ahtn <- cbind(as.data.frame(prop.table(table(data.frame(lg.baseDHH$NACCAHTN,
lg.baseDHH$Diagnosis)), 1)), as.data.frame(with(lg.baseDHH, table(NACCAHTN, Diagnosis))))

apsy <- cbind(as.data.frame(prop.table(table(data.frame(lg.baseDHH$NACCAPSY,
lg.baseDHH$Diagnosis)), 1)), as.data.frame(with(lg.baseDHH, table(NACCAPSY, Diagnosis))))

dbmd <- cbind(as.data.frame(prop.table(table(data.frame(lg.baseDHH$NACCDBMD,
lg.baseDHH$Diagnosis)), 1)), as.data.frame(with(lg.baseDHH, table(NACCDBMD, Diagnosis))))

lipl <- cbind(as.data.frame(prop.table(table(data.frame(lg.baseDHH$NACCLIPL,
lg.baseDHH$Diagnosis)), 1)), as.data.frame(with(lg.baseDHH, table(NACCLIPL, Diagnosis))))

nsd <- cbind(as.data.frame(prop.table(table(data.frame(lg.baseDHH$NACCNSD,
lg.baseDHH$Diagnosis)), 1)), as.data.frame(with(lg.baseDHH, table(NACCNSD, Diagnosis))))

pdmd <- cbind(as.data.frame(prop.table(table(data.frame(lg.baseDHH$NACCPDMD,
lg.baseDHH$Diagnosis)), 1)), as.data.frame(with(lg.baseDHH, table(NACCPDMD, Diagnosis))))

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
write.excel(aanx)
write.excel(antc)
write.excel(adep)
write.excel(ahtn)
write.excel(apsy)
write.excel(dbmd)
write.excel(lipl)
write.excel(nsd)
write.excel(pdmd)
write.excel(diab)
write.excel(hyper)
write.excel(dep)
write.excel(bmi)
write.excel(tbi)
write.excel(cvd)

#################### Antihypertensives ##################################

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
nacc <- fread("MDhypercdr.csv")
df <- nacc

df$RACE[df$RACE==50] <- NA
df$SEX[df$SEX==2] <- 0  # female 0, male 1
df <- df[df$NEW_VISITNUM == 1, ] 
df <- subset(df, NACCAGE>=40,
             select=c(NACCID:ANGI))
df <- as.data.frame(df)
df <- df %>% filter(SEX == 1)

df$Diagnosis[df$Diagnosis==1 & df$ALZP==1] <- 2
df <- df[df$Diagnosis != 1, ]
df$Diagnosis[df$Diagnosis==2] <- 1

## 9 means medication mentioned only once, decide if need to remove or not
df <- df %>% 
 replace_with_na(replace = list(NACCHTNC = 9, NACCACEI = 9, NACCAAAS = 9,
                               NACCBETA = 9, NACCCCBS = 9, NACCDIUR = 9, 
                              NACCVASD = 9, NACCANGI = 9))


df.baseDHH <- df[,c(1,8:10,14,26,107,110,114,116,117,119,120,124,127,130:137)]
df.baseDHH[ , 7:14 ][ df.baseDHH[ , 7:14 ] == 9 ] <- 1

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
lg.baseDHH$NACCHTNC <- as.factor(as.character(lg.baseDHH$NACCHTNC))
lg.baseDHH$NACCACEI <- as.factor(as.character(lg.baseDHH$NACCACEI))
lg.baseDHH$NACCAAAS <- as.factor(as.character(lg.baseDHH$NACCAAAS))
lg.baseDHH$NACCBETA <- as.factor(as.character(lg.baseDHH$NACCBETA))
lg.baseDHH$NACCCCBS <- as.factor(as.character(lg.baseDHH$NACCCCBS))
lg.baseDHH$NACCDIUR <- as.factor(as.character(lg.baseDHH$NACCDIUR))
lg.baseDHH$NACCVASD <- as.factor(as.character(lg.baseDHH$NACCVASD))
lg.baseDHH$NACCANGI <- as.factor(as.character(lg.baseDHH$NACCANGI))
lg.baseDHH$Diagnosis <- as.factor(as.character(lg.baseDHH$Diagnosis))
lg.baseDHH <- as.data.frame(lg.baseDHH)

riskmodel <- fitLogRegModel(lg.baseDHH, cOutcome=15,
                            cNonGenPreds=c(2:14,16:23), cNonGenPredsCat=c(2,4:14,16:18,21:23),
                            cGenPreds=c(0), cGenPredsCat=c(0))

#gender
riskmodel <- fitLogRegModel(lg.baseDHH, cOutcome=15,
                            cNonGenPreds=c(3:14,16:23), cNonGenPredsCat=c(4:14,16:18,21:23),
                            cGenPreds=c(0), cGenPredsCat=c(0))

summary(riskmodel)
oddr <- ORmultivariate(riskmodel,filename="medtest.txt")

vif <- vif(riskmodel)
stat <- cbind(rownames(vif), data.frame(vif, row.names=NULL))
write.csv(stat, file = "vif.csv",row.names=FALSE, na="")

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

aaas <- cbind(as.data.frame(prop.table(table(data.frame(lg.baseDHH$NACCAAAS,
lg.baseDHH$Diagnosis)), 1)), as.data.frame(with(lg.baseDHH, table(NACCAAAS, Diagnosis))))

acei <- cbind(as.data.frame(prop.table(table(data.frame(lg.baseDHH$NACCACEI,
lg.baseDHH$Diagnosis)), 1)), as.data.frame(with(lg.baseDHH, table(NACCACEI, Diagnosis))))

angi <- cbind(as.data.frame(prop.table(table(data.frame(lg.baseDHH$NACCANGI,
lg.baseDHH$Diagnosis)), 1)), as.data.frame(with(lg.baseDHH, table(NACCANGI, Diagnosis))))

htnc <- cbind(as.data.frame(prop.table(table(data.frame(lg.baseDHH$NACCHTNC,
lg.baseDHH$Diagnosis)), 1)), as.data.frame(with(lg.baseDHH, table(NACCHTNC, Diagnosis))))

beta <- cbind(as.data.frame(prop.table(table(data.frame(lg.baseDHH$NACCBETA,
lg.baseDHH$Diagnosis)), 1)), as.data.frame(with(lg.baseDHH, table(NACCBETA, Diagnosis))))

ccbs <- cbind(as.data.frame(prop.table(table(data.frame(lg.baseDHH$NACCCCBS,
lg.baseDHH$Diagnosis)), 1)), as.data.frame(with(lg.baseDHH, table(NACCCCBS, Diagnosis))))

diur <- cbind(as.data.frame(prop.table(table(data.frame(lg.baseDHH$NACCDIUR,
lg.baseDHH$Diagnosis)), 1)), as.data.frame(with(lg.baseDHH, table(NACCDIUR, Diagnosis))))

vasd <- cbind(as.data.frame(prop.table(table(data.frame(lg.baseDHH$NACCVASD,
lg.baseDHH$Diagnosis)), 1)), as.data.frame(with(lg.baseDHH, table(NACCVASD, Diagnosis))))

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
write.excel(aaas)
write.excel(acei)
write.excel(angi)
write.excel(htnc)
write.excel(beta)
write.excel(ccbs)
write.excel(diur)
write.excel(vasd)
write.excel(diab)
write.excel(hyper)
write.excel(dep)
write.excel(bmi)
write.excel(tbi)
write.excel(cvd)
