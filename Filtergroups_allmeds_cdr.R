rm(list = ls())

library(data.table)
library(lubridate)
library(dplyr)
library(reshape)
library(naniar)
library(tidyverse)

File <- file.path("C:/Users/Daman Kaur/Desktop/Drug classes DO")
setwd(File)
nacc <- fread("kaur10262020.csv")

##Select chosen variables
sub.nacc <- data.frame(NACCID = nacc$NACCID, VISITDAY = nacc$VISITDAY, VISITMO = nacc$VISITMO,
  VISITYR = nacc$VISITYR, MARI = nacc$MARISTAT, VNUM = nacc$NACCVNUM, ALVST = nacc$NACCAVST,
  FAM = nacc$NACCFAM, SEX = nacc$SEX, NACCAGE = nacc$NACCAGE, RACE = nacc$RACE, BMI = nacc$NACCBMI, 
  SMOKE = nacc$SMOKYRS, EDUC = nacc$EDUC, ALCOHOL = nacc$ALCOHOL, HYPERC1 = nacc$HYPERCHO,
  HYPERC2 = nacc$HYPCHOL, DIAB1 = nacc$DIABETES, DIAB2 = nacc$DIABET, HYPERT1 = nacc$HYPERTEN, 
  HYPERT2 = nacc$HYPERT, MEDN = nacc$NACCAMD, TBI1 = nacc$NACCTBI, TBI2 = nacc$TRAUMBRF,
  TBI3 = nacc$BRNINJ, TBI4 = nacc$TBI, HEAR = nacc$HEARING, ARTH1 = nacc$ARTH, ARTH2 = nacc$ARTHRIT,
  CBD1 = nacc$CVDCOG, CBD2 = nacc$CBSTROKE, CBD3 = nacc$CBTIA, CBD4 = nacc$STROKE, 
  CBD5 = nacc$CVD, PARK1 = nacc$PARK, PARK2 = nacc$PD, DEP1 = nacc$DEP, DEP2 = nacc$DEP2YRS,
  DEP3 = nacc$DEPOTHR, ANX1 =nacc$ANXIETY, ANX2 = nacc$ANXIET, ANX3 = nacc$OCD,
  APS1 = nacc$BIPOLAR, APS2 = nacc$BIPOLDX, APS3 = nacc$SCHIZ, APS4 = nacc$SCHIZOP,
  APS5 = nacc$DELIR, UDSD = nacc$NACCUDSD, DIAG = nacc$CDRSUM, ALZP = nacc$NACCALZP, 
  ALZD = nacc$NACCALZD)

sub.nacc <- as.data.frame(sub.nacc)

ds <- nacc[,c(45:84,92,95:103,368:370,372,390:393,401:407,409:414,417:421)]
ds <- as.data.frame(ds)

sub.nacc <- cbind(sub.nacc,ds)
sub.nacc <- as.data.frame(sub.nacc)

sub.nacc <- sub.nacc %>%
  mutate(Date = make_date(VISITYR, VISITMO, VISITDAY))
sub.nacc <- sub.nacc[,-(2)] 

sub.nacc[sub.nacc==-4]<-NA
sub.nacc <- sub.nacc %>% 
  replace_with_na(replace = list(FAM = 9, RACE = 99, BMI = 888.8, EDUC = 99, ALCOHOL = 9,
           HYPERC1 = 9, HYPERC2 = 8,  DIAB1 = 9, DIAB2 = 9,  HYPERT1 = 9, HYPERT2 = 8,                
           TBI1 = 9, TBI2 = 9, TBI4 = 9, HEAR = 9, ARTH1 = 8, ARTH2 = 9, CBD1 = 8,                    
           CBD2 = 9, CBD3 = 9, PARK2 = 9, DEP2 = 9, DEP3 = 9, ANX1 = 9, ANX3 = 9, APS1 = 9,
           APS3 = 9, SMOKE = 99, CVHATT = 9, CVAFIB = 9, CVANGIO = 9, CVBYPASS = 9,
           CVPACDEF = 9, CVPACE = 9, CVCHF = 9, CVANGINA = 9, CVHVALVE = 9, CVOTHR = 9,
           MYOINF = 8, CONGHRT = 8, AFIBRILL = 8, ANGINA = 8, ANGIOCP = 8, ANGIOPCI = 8,
           PACEMAKE = 8, HVALVE = 8))
sub.nacc$SMOKE[sub.nacc$SMOKE == 88] <-  0

#Order according to ID and year of visit

sub.nacc.sort <- sub.nacc[ order(sub.nacc[,1], sub.nacc[,3],sub.nacc[,2]), ]
sub.nacc <- sub.nacc.sort

## Filtering out people with only 1 observation
df <- sub.nacc %>% group_by(NACCID) %>% filter(n()>1)

## new column with new visit number ## 
df <- df %>% group_by(NACCID) %>% mutate(NEW_VISITNUM = row_number())

df$DIAG[(df$DIAG>=0.5 & df$DIAG<=4.0)]<- 3 
df$DIAG[df$DIAG==0] <- 1
df$DIAG[(df$DIAG>=4.5 & df$DIAG<=18.0)] <- 4

## Filter IDs with ascending diagnosis, clear progression trend obtained

is.sorted = Negate(is.unsorted)

dfe <- df %>%
  group_by(NACCID) %>%
  mutate(check1=is.sorted(DIAG)) 

dfe <- dfe %>%
  group_by(NACCID) %>%
  filter(check1=="TRUE")

dfe <- dfe[,-(129)] 

#### FOR healthy to dementia progression, Extract groups with baseline Normal ####

dfe <- dfe %>%
  group_by(NACCID) %>%
  mutate(check1=any(DIAG==1 & NEW_VISITNUM==1)) 

dfe <- dfe %>%
  group_by(NACCID) %>%
  filter(check1=="TRUE")

dfe <- dfe[,-(129)] 

## Now we have stable healthy, healthy to MCI, and healthy to dementia.
## Filter out groups with MCI at last visit.

## MCI
dfe <- dfe %>%
  group_by(NACCID) %>%
  mutate(check1=any(DIAG==3 & NEW_VISITNUM==(n()))) 
dfe <- dfe %>%
  group_by(NACCID) %>%
  filter(check1=="FALSE")
dfe <- dfe[,-(129)] 
nacc <- dfe

## Filter dementia visits, renumber and filter 1st dementia visit 
#(only need 1st dementia visit to identify those who progressed)

dfe <- nacc[nacc$DIAG == 4, ]

dfe <- dfe[,-(128)] 
dfe <- dfe %>% group_by(NACCID) %>% mutate(NEW_VISITNUM = row_number())
dfe <- dfe[dfe$NEW_VISITNUM == 1, ] 

df <- nacc %>%
  group_by(NACCID) %>%
  mutate(check1=any(DIAG==4)) 
xy <- df[df$check1 == "TRUE", ]
xy <- xy[xy$DIAG != 4, ]
xy <- xy[,-(129)] 
mix <- rbind(dfe, xy)

zx <- nacc %>%
  group_by(NACCID) %>%
  mutate(check1=any(DIAG==4)) 
dfe <- zx[zx$check1 == "FALSE", ]
dfe <- dfe[,-(129)] 
nacc <- rbind(dfe, mix) # Now we have stable healthy and healthy to dementia

nacc <- as.data.frame(nacc)
dfe <- nacc[ order(nacc[,1], nacc[,3], nacc[,2]), ]
dfe <- dfe[,-(128)] 

dl <- dfe %>% group_by(NACCID) %>% mutate(NEW_VISITNUM = row_number())
df <- dl %>%
  group_by(NACCID) %>%
  mutate(check1=any(DIAG==4)) 
df$Diagnosis <- df$check1
df$Diagnosis[(df$Diagnosis == 'FALSE' )]<- 0
df$Diagnosis[(df$Diagnosis == 'TRUE' )]<- 1
df <- df[,-(129)] 

write.csv(df, file = "NtoDcdr.csv",row.names=FALSE, na="")

#### FOR healthy to MCI PROGRESSION, Filter groups with baseline Normal ####

dfe <- dfe %>%
  group_by(NACCID) %>%
  mutate(check1=any(DIAG==1 & NEW_VISITNUM==1)) 

dfe <- dfe %>%
  group_by(NACCID) %>%
  filter(check1=="TRUE")

dfe <- dfe[,-(129)] 

## Now we have healthy to healthy,  healthy to MCI,
#and Healthy to dementia. 

dfe <- dfe %>%
  group_by(NACCID) %>%
  mutate(check1=any(DIAG==4)) 

dfe <- dfe %>%
  group_by(NACCID) %>%
  filter(check1=="FALSE")

dfe <- dfe[,-(129)] 

## Now we have healthy to healthy and healthy to MCI

nacc <- dfe
df <- nacc %>%
  group_by(NACCID) %>%
  mutate(check1=any(DIAG==3)) 
dfe <- df[df$check1 == "FALSE", ]
dfe <- dfe[,-(129)] 

xy <- nacc[nacc$DIAG == 3, ]

#renumber and filter 1st MCI visit (Only need first MCI visit to identify those who progressed)
xy <- xy[,-(128)] 
xy <- xy %>% group_by(NACCID) %>% mutate(NEW_VISITNUM = row_number())
xy <- xy[xy$NEW_VISITNUM == 1, ] 

mix <- rbind(dfe, xy)

df <- nacc %>%
  group_by(NACCID) %>%
  mutate(check1=any(DIAG==3)) 
dfe <- df[df$check1 == "TRUE", ]
dfe <- dfe[dfe$DIAG != 3, ]
dfe <- dfe[,-(129)] 
nacc <- rbind(dfe, mix)

nacc <- as.data.frame(nacc)
dfe <- nacc[ order(nacc[,1], nacc[,3], nacc[,2]), ]
dfe <- dfe[,-(128)] 

dl <- dfe %>% group_by(NACCID) %>% mutate(NEW_VISITNUM = row_number())
df <- dl %>%
  group_by(NACCID) %>%
  mutate(check1=any(DIAG==3)) 
df$Diagnosis <- df$check1
df$Diagnosis[(df$Diagnosis == 'FALSE' )]<- 0
df$Diagnosis[(df$Diagnosis == 'TRUE' )]<- 1
df <- df[,-(129)] 
write.csv(df, file = "NtoMcdr.csv",row.names=FALSE, na="")

#### For MCI to Dementia, extract individuals with MCI on 1st visit####
## You get MCI to MCI and MCI to Dementia

df <- dfe %>%
  group_by(NACCID) %>%
  mutate(check1=any(DIAG==3 & NEW_VISITNUM==1)) 
dfe <- df[df$check1 == "TRUE", ]
dfe <- dfe[,-(129)] 

nacc <- dfe

df <- nacc %>%
  group_by(NACCID) %>%
  mutate(check1=any(DIAG==4)) 
dfe <- df[df$check1 == "FALSE", ]
dfe <- dfe[,-(129)] 
xy <- nacc[nacc$DIAG == 4, ]

#renumber and filter 1st dementia visit (only need 1st dementia visit to identify those who progressed)
xy <- xy[,-(128)] 
xy <- xy %>% group_by(NACCID) %>% mutate(NEW_VISITNUM = row_number())
xy <- xy[xy$NEW_VISITNUM == 1, ] 

mix <- rbind(dfe, xy)

df <- nacc %>%
  group_by(NACCID) %>%
  mutate(check1=any(DIAG==4)) 
dfe <- df[df$check1 == "TRUE", ]
dfe <- dfe[dfe$DIAG != 4, ]
dfe <- dfe[,-(129)] 
nacc <- rbind(dfe, mix)

nacc <- as.data.frame(nacc)
dfe <- nacc[ order(nacc[,1], nacc[,3], nacc[,2]), ]
dfe <- dfe[,-(128)] 

dl <- dfe %>% group_by(NACCID) %>% mutate(NEW_VISITNUM = row_number())
df <- dl %>%
  group_by(NACCID) %>%
  mutate(check1=any(DIAG==4)) 
df$Diagnosis <- df$check1
df$Diagnosis[(df$Diagnosis == 'FALSE' )]<- 0
df$Diagnosis[(df$Diagnosis == 'TRUE' )]<- 1
df <- df[,-(129)] 
write.csv(df, file = "MtoDcdr.csv",row.names=FALSE, na="")
