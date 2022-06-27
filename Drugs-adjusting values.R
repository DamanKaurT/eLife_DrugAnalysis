################## All meds ##########################

rm(list = ls())

library(data.table)
library(lubridate)
library(dplyr)
library(reshape)
library(naniar)

File <- file.path("C:/Users/Daman Kaur/Desktop/Drug classes DO")
setwd(File)
nacc <- fread("MtoDcdr.csv")

nacc <- nacc %>% 
  group_by(NACCID) %>%
  mutate(ALZP1 = last(ALZP))
nacc <- nacc %>% 
  group_by(NACCID) %>%
  mutate(ALZD1 = last(ALZD))
nacc <- nacc[,-(49:50)] 

### Only Keep required observations (can remove the 1st MCI (from NtoM) or dementia (from NtoD/MtoD) visit)
df <- subset(nacc, DIAG<4,                      #NtoD/MtoD DIAG<4  #NtoM DIAG<3
             select=c(NACCID:ALZD1))

df <- df %>% group_by(NACCID) %>% filter(!all(is.na(SMOKE)))
df <- df %>% group_by(NACCID) %>% filter(!all(is.na(EDUC)))
df <- df %>% group_by(NACCID) %>% filter(!all(is.na(BMI)))
df <- df %>% group_by(NACCID) %>% filter(!all(is.na(ALCOHOL)))
df <- df %>% group_by(NACCID) %>% filter(!all(is.na(NACCAHTN)))
df <- df %>% group_by(NACCID) %>% filter(!all(is.na(NACCLIPL)))
df <- df %>% group_by(NACCID) %>% filter(!all(is.na(NACCNSD)))
df <- df %>% group_by(NACCID) %>% filter(!all(is.na(NACCAC)))
df <- df %>% group_by(NACCID) %>% filter(!all(is.na(NACCADEP)))
df <- df %>% group_by(NACCID) %>% filter(!all(is.na(NACCAPSY)))
df <- df %>% group_by(NACCID) %>% filter(!all(is.na(NACCAANX)))
df <- df %>% group_by(NACCID) %>% filter(!all(is.na(NACCDBMD)))
df <- df %>% group_by(NACCID) %>% filter(!all(is.na(NACCADMD)))
df <- df %>% group_by(NACCID) %>% filter(!all(is.na(NACCPDMD)))
df <- df %>% group_by(NACCID) %>% filter(!all(is.na(HEAR)))

df <- df %>%
  group_by(NACCID) %>%
  mutate(check1=is.na(HYPERT1) & is.na(HYPERT2))
df <- df %>% group_by(NACCID) %>% filter(!all(check1=="TRUE"))

df <- df %>%
  group_by(NACCID) %>%
  mutate(check1=is.na(DIAB1) & is.na(DIAB2))
df <- df %>% group_by(NACCID) %>% filter(!all(check1=="TRUE"))

df <- df %>%
  group_by(NACCID) %>%
  mutate(check1=is.na(TBI1) & is.na(TBI2) & is.na(TBI3) & is.na(TBI4))
df <- df %>% group_by(NACCID) %>% filter(!all(check1=="TRUE"))

df <- df %>%
  group_by(NACCID) %>%
  mutate(check2=is.na(DEP2) & is.na(DEP3)) 
df <- df %>% group_by(NACCID) %>% filter(!all(check2=="TRUE"))

df <- df %>%
  group_by(NACCID) %>%
  mutate(check2=is.na(CVHATT) & is.na(CVAFIB) & is.na(CVANGIO) & is.na(CVBYPASS) 
         & is.na(CVCHF) & is.na(CVOTHR) & is.na(CVANGINA) & is.na(CVHVALVE) & is.na(CVPACE) 
         & is.na(CVPACDEF) & is.na(MYOINF) & is.na(CONGHRT) & is.na(AFIBRILL) & is.na(ANGINA)
         & is.na(ANGIOCP) & is.na(ANGIOPCI) & is.na(PACEMAKE) & is.na(HVALVE) & is.na(CBD2)
         & is.na(CBD3) & is.na(CBD4) & is.na(CBD5)) 
df <- df %>% group_by(NACCID) %>% filter(!all(check2=="TRUE"))


df <- df[,-(130:131)] 
df <- as.data.frame(df)

df <- df %>% group_by(NACCID) %>% filter(n()>1)
df <- df[ order(df[,1], df[,3], df[,2]), ]
df <- df %>% group_by(NACCID) %>% mutate(NEW_VISITNUM = row_number())

### Cleaning data ###

### DIABETES ###
df <- df %>% 
  group_by(NACCID) %>%
  mutate(DCDC = sum(DIAB1, na.rm=TRUE))

df$DIAB1[(df$DCDC==0 & df$NEW_VISITNUM==1)]<- 0 
df$DIAB1[(df$DCDC>=1 & df$NEW_VISITNUM==1)]<- 1 
df <- df[,-(130)] 
df <- as.data.frame(df)

df <- df %>% 
  group_by(NACCID) %>%
  mutate(DCDC = sum(DIAB2, na.rm=TRUE))

df$DIAB2[(df$DCDC==0 & df$NEW_VISITNUM==1)]<- 0 
df$DIAB2[(df$DCDC>=1 & df$NEW_VISITNUM==1)]<- 1 
df <- df[,-(130)] 
df <- as.data.frame(df)

df$DIAB[(df$DIAB1==1 | df$DIAB2==1)] <- 1
df$DIAB[(df$DIAB1==0 & df$DIAB2==0)] <- 0


### HYPERTENSION ###
df <- df %>% 
  group_by(NACCID) %>%
  mutate(DCDC = sum(HYPERT1, na.rm=TRUE))

df$HYPERT1[(df$DCDC==0 & df$NEW_VISITNUM==1)]<- 0 
df$HYPERT1[(df$DCDC>=1 & df$NEW_VISITNUM==1)]<- 1 
df <- df[,-(131)] 
df <- as.data.frame(df)

df <- df %>% 
  group_by(NACCID) %>%
  mutate(DCDC = sum(HYPERT2, na.rm=TRUE))

df$HYPERT2[(df$DCDC==0 & df$NEW_VISITNUM==1)]<- 0 
df$HYPERT2[(df$DCDC>=1 & df$NEW_VISITNUM==1)]<- 1 
df <- df[,-(131)] 
df <- as.data.frame(df)

df$HYPERT[(df$HYPERT1==1 | df$HYPERT2==1)] <- 1
df$HYPERT[(df$HYPERT1==0 & df$HYPERT2==0)] <- 0

### ALCOHOLISM ###
df <- df %>% 
  group_by(NACCID) %>%
  mutate(DCDC = sum(ALCOHOL, na.rm=TRUE))

df$ALCOHOL[(df$DCDC==0 & df$NEW_VISITNUM==1)]<- 0 
df$ALCOHOL[(df$DCDC>=1 & df$NEW_VISITNUM==1)]<- 1 
df <- df[,-(132)] 
df <- as.data.frame(df)

### DEPRESSION ###
df <- df %>% 
  group_by(NACCID) %>%
  mutate(DCDC = sum(DEP2, na.rm=TRUE))
df$DEP2[(df$DCDC==0 & df$NEW_VISITNUM==1)]<- 0 
df$DEP2[(df$DCDC>=1 & df$NEW_VISITNUM==1)]<- 1 
df <- df[,-(132)] 
df <- as.data.frame(df)

df <- df %>% 
  group_by(NACCID) %>%
  mutate(DCDC = sum(DEP3, na.rm=TRUE))
df$DEP3[(df$DCDC==0 & df$NEW_VISITNUM==1)]<- 0 
df$DEP3[(df$DCDC>=1 & df$NEW_VISITNUM==1)]<- 1 
df <- df[,-(132)] 
df <- as.data.frame(df)

df$DEPRSN[(df$DEP2==1 | df$DEP3==1)] <- 1
df$DEPRSN[(df$DEP2==0 & df$DEP3==0)] <- 0

### HEARING ###
df <- df %>% 
  group_by(NACCID) %>%
  mutate(DCDC = sum(HEAR, na.rm=TRUE))

df$HEAR[(df$DCDC==0 & df$NEW_VISITNUM==1)]<- 0 
df$HEAR[(df$DCDC>=1 & df$NEW_VISITNUM==1)]<- 1 
df <- df[,-(133)] 
df <- as.data.frame(df)

### Number of years smoked ###
df <- df %>% group_by(NACCID) %>% mutate(SMAX = max(SMOKE, na.rm=TRUE))
df <- as.data.frame(df)

### Education ###
df <- df %>% group_by(NACCID) %>% mutate(EDMAX = max(EDUC, na.rm=TRUE))
df <- as.data.frame(df)

### BMI ###
df$DCDC[df$BMI<18.5] <-  1
df$DCDC[df$BMI>=18.5 & df$BMI<=24.99]<-  2
df$DCDC[df$BMI> 24.99] <- 3

df <- df %>%
  group_by(NACCID) %>%
  mutate(check1=mean(DCDC, na.rm = TRUE))

df <- df %>%
  group_by(NACCID) %>% 
  mutate(value = first(DCDC[complete.cases(DCDC)]))

df$ABMI[(df$check1==df$value & df$NEW_VISITNUM==1)]<- 0 ## Stable
df$ABMI[(df$check1<df$value & df$NEW_VISITNUM==1)]<- 1  ## decreasing
df$ABMI[(df$check1>df$value & df$NEW_VISITNUM==1)]<- 2  ## increasing
df <- df[,-(135:137)] 
df <- as.data.frame(df)


### TBI ###
df <- df %>% 
  group_by(NACCID) %>%
  mutate(DCDC = sum(TBI1, na.rm=TRUE))

df$TBI1[(df$DCDC==0 & df$NEW_VISITNUM==1)]<- 0 
df$TBI1[(df$DCDC>=1 & df$NEW_VISITNUM==1)]<- 1 
df <- df[,-(136)] 

df <- df %>% 
  group_by(NACCID) %>%
  mutate(DCDC = sum(TBI2, na.rm=TRUE))

df$TBI2[(df$DCDC==0 & df$NEW_VISITNUM==1)]<- 0 
df$TBI2[(df$DCDC>=1 & df$NEW_VISITNUM==1)]<- 1 
df <- df[,-(136)] 

df <- df %>% 
  group_by(NACCID) %>%
  mutate(DCDC = sum(TBI3, na.rm=TRUE))

df$TBI3[(df$DCDC==0 & df$NEW_VISITNUM==1)]<- 0 
df$TBI3[(df$DCDC>=1 & df$NEW_VISITNUM==1)]<- 1 
df <- df[,-(136)] 

df <- df %>% 
  group_by(NACCID) %>%
  mutate(DCDC = sum(TBI4, na.rm=TRUE))

df$TBI4[(df$DCDC==0 & df$NEW_VISITNUM==1)]<- 0 
df$TBI4[(df$DCDC>=1 & df$NEW_VISITNUM==1)]<- 1 
df <- df[,-(136)] 

df <- as.data.frame(df)

df$TBI[(df$TBI1==1 | df$TBI2==1 | df$TBI3==1 | df$TBI4==1)] <- 1
df$TBI[(df$TBI1==0 & df$TBI2==0 & df$TBI3==0 & df$TBI4==0)] <- 0


### CARDIOVASCULAR DISEASES ###

df <- df %>% 
  group_by(NACCID) %>%
  mutate(DCDC = sum(CVHATT, na.rm=TRUE))
df$CVHATT[(df$DCDC==0 & df$NEW_VISITNUM==1)]<- 0 
df$CVHATT[(df$DCDC>=1 & df$NEW_VISITNUM==1)]<- 1 
df <- df[,-(137)] 

df <- df %>% 
  group_by(NACCID) %>%
  mutate(DCDC = sum(CVAFIB, na.rm=TRUE))
df$CVAFIB[(df$DCDC==0 & df$NEW_VISITNUM==1)]<- 0 
df$CVAFIB[(df$DCDC>=1 & df$NEW_VISITNUM==1)]<- 1 
df <- df[,-(137)] 

df <- df %>% 
  group_by(NACCID) %>%
  mutate(DCDC = sum(CVANGIO, na.rm=TRUE))
df$CVANGIO[(df$DCDC==0 & df$NEW_VISITNUM==1)]<- 0 
df$CVANGIO[(df$DCDC>=1 & df$NEW_VISITNUM==1)]<- 1 
df <- df[,-(137)] 

df <- df %>% 
  group_by(NACCID) %>%
  mutate(DCDC = sum(CVBYPASS, na.rm=TRUE))
df$CVBYPASS[(df$DCDC==0 & df$NEW_VISITNUM==1)]<- 0 
df$CVBYPASS[(df$DCDC>=1 & df$NEW_VISITNUM==1)]<- 1 
df <- df[,-(137)] 

df <- df %>% 
  group_by(NACCID) %>%
  mutate(DCDC = sum(CVCHF, na.rm=TRUE))
df$CVCHF[(df$DCDC==0 & df$NEW_VISITNUM==1)]<- 0 
df$CVCHF[(df$DCDC>=1 & df$NEW_VISITNUM==1)]<- 1 
df <- df[,-(137)] 

df <- df %>% 
  group_by(NACCID) %>%
  mutate(DCDC = sum(CVOTHR, na.rm=TRUE))
df$CVOTHR[(df$DCDC==0 & df$NEW_VISITNUM==1)]<- 0 
df$CVOTHR[(df$DCDC>=1 & df$NEW_VISITNUM==1)]<- 1 
df <- df[,-(137)] 

df <- df %>% 
  group_by(NACCID) %>%
  mutate(DCDC = sum(CVANGINA, na.rm=TRUE))
df$CVANGINA[(df$DCDC==0 & df$NEW_VISITNUM==1)]<- 0 
df$CVANGINA[(df$DCDC>=1 & df$NEW_VISITNUM==1)]<- 1 
df <- df[,-(137)] 

df <- df %>% 
  group_by(NACCID) %>%
  mutate(DCDC = sum(CVHVALVE, na.rm=TRUE))
df$CVHVALVE[(df$DCDC==0 & df$NEW_VISITNUM==1)]<- 0 
df$CVHVALVE[(df$DCDC>=1 & df$NEW_VISITNUM==1)]<- 1 
df <- df[,-(137)] 

df <- df %>% 
  group_by(NACCID) %>%
  mutate(DCDC = sum(CVPACE, na.rm=TRUE))
df$CVPACE[(df$DCDC==0 & df$NEW_VISITNUM==1)]<- 0 
df$CVPACE[(df$DCDC>=1 & df$NEW_VISITNUM==1)]<- 1 
df <- df[,-(137)] 

df <- df %>% 
  group_by(NACCID) %>%
  mutate(DCDC = sum(CVPACDEF, na.rm=TRUE))
df$CVPACDEF[(df$DCDC==0 & df$NEW_VISITNUM==1)]<- 0 
df$CVPACDEF[(df$DCDC>=1 & df$NEW_VISITNUM==1)]<- 1 
df <- df[,-(137)] 

df <- df %>% 
  group_by(NACCID) %>%
  mutate(DCDC = sum(MYOINF, na.rm=TRUE))
df$MYOINF[(df$DCDC==0 & df$NEW_VISITNUM==1)]<- 0 
df$MYOINF[(df$DCDC>=1 & df$NEW_VISITNUM==1)]<- 1 
df <- df[,-(137)] 

df <- df %>% 
  group_by(NACCID) %>%
  mutate(DCDC = sum(CONGHRT, na.rm=TRUE))
df$CONGHRT[(df$DCDC==0 & df$NEW_VISITNUM==1)]<- 0 
df$CONGHRT[(df$DCDC>=1 & df$NEW_VISITNUM==1)]<- 1 
df <- df[,-(137)] 

df <- df %>% 
  group_by(NACCID) %>%
  mutate(DCDC = sum(AFIBRILL, na.rm=TRUE))
df$AFIBRILL[(df$DCDC==0 & df$NEW_VISITNUM==1)]<- 0 
df$AFIBRILL[(df$DCDC>=1 & df$NEW_VISITNUM==1)]<- 1 
df <- df[,-(137)] 

df <- df %>% 
  group_by(NACCID) %>%
  mutate(DCDC = sum(ANGINA, na.rm=TRUE))
df$ANGINA[(df$DCDC==0 & df$NEW_VISITNUM==1)]<- 0 
df$ANGINA[(df$DCDC>=1 & df$NEW_VISITNUM==1)]<- 1 
df <- df[,-(137)] 

df <- df %>% 
  group_by(NACCID) %>%
  mutate(DCDC = sum(ANGIOCP, na.rm=TRUE))
df$ANGIOCP[(df$DCDC==0 & df$NEW_VISITNUM==1)]<- 0 
df$ANGIOCP[(df$DCDC>=1 & df$NEW_VISITNUM==1)]<- 1 
df <- df[,-(137)] 

df <- df %>% 
  group_by(NACCID) %>%
  mutate(DCDC = sum(ANGIOPCI, na.rm=TRUE))
df$ANGIOPCI[(df$DCDC==0 & df$NEW_VISITNUM==1)]<- 0 
df$ANGIOPCI[(df$DCDC>=1 & df$NEW_VISITNUM==1)]<- 1 
df <- df[,-(137)] 

df <- df %>% 
  group_by(NACCID) %>%
  mutate(DCDC = sum(PACEMAKE, na.rm=TRUE))
df$PACEMAKE[(df$DCDC==0 & df$NEW_VISITNUM==1)]<- 0 
df$PACEMAKE[(df$DCDC>=1 & df$NEW_VISITNUM==1)]<- 1 
df <- df[,-(137)] 

df <- df %>% 
  group_by(NACCID) %>%
  mutate(DCDC = sum(HVALVE, na.rm=TRUE))
df$HVALVE[(df$DCDC==0 & df$NEW_VISITNUM==1)]<- 0 
df$HVALVE[(df$DCDC>=1 & df$NEW_VISITNUM==1)]<- 1 
df <- df[,-(137)] 

df <- df %>% 
  group_by(NACCID) %>%
  mutate(DCDC = sum(CBD2, na.rm=TRUE))
df$CBD2[(df$DCDC==0 & df$NEW_VISITNUM==1)]<- 0 
df$CBD2[(df$DCDC>=1 & df$NEW_VISITNUM==1)]<- 1 
df <- df[,-(137)] 

df <- df %>% 
  group_by(NACCID) %>%
  mutate(DCDC = sum(CBD3, na.rm=TRUE))
df$CBD3[(df$DCDC==0 & df$NEW_VISITNUM==1)]<- 0 
df$CBD3[(df$DCDC>=1 & df$NEW_VISITNUM==1)]<- 1 
df <- df[,-(137)] 

df <- df %>% 
  group_by(NACCID) %>%
  mutate(DCDC = sum(CBD4, na.rm=TRUE))
df$CBD4[(df$DCDC==0 & df$NEW_VISITNUM==1)]<- 0 
df$CBD4[(df$DCDC>=1 & df$NEW_VISITNUM==1)]<- 1 
df <- df[,-(137)] 

df <- df %>% 
  group_by(NACCID) %>%
  mutate(DCDC = sum(CBD5, na.rm=TRUE))
df$CBD5[(df$DCDC==0 & df$NEW_VISITNUM==1)]<- 0 
df$CBD5[(df$DCDC>=1 & df$NEW_VISITNUM==1)]<- 1 
df <- df[,-(137)] 

df <- as.data.frame(df)

df$CVD[(df$CVHATT==1 | df$CVAFIB==1 | df$CVANGIO==1 | df$CVBYPASS==1  | df$CVCHF==1 | 
          df$CVOTHR==1 | df$CVANGINA==1  | df$CVHVALVE==1 | df$CVPACE==1 | df$CVPACDEF==1  |
          df$MYOINF==1 | df$CONGHRT==1 | df$AFIBRILL==1  | df$ANGINA==1 | df$ANGIOCP==1 | 
          df$ANGIOPCI==1  | df$PACEMAKE==1 | df$HVALVE==1 | df$CBD2==1  | df$CBD3==1 | 
          df$CBD4==1 | df$CBD5==1)] <- 1
df$CVD[(df$CVHATT==0 & df$CVAFIB==0 & df$CVANGIO==0 & df$CVBYPASS==0  & df$CVCHF==0 & 
          df$CVOTHR==0 & df$CVANGINA==0  & df$CVHVALVE==0 & df$CVPACE==0 & df$CVPACDEF==0  &
          df$MYOINF==0 & df$CONGHRT==0 & df$AFIBRILL==0  & df$ANGINA==0 & df$ANGIOCP==0 & 
          df$ANGIOPCI==0  & df$PACEMAKE==0 & df$HVALVE==0 & df$CBD2==0  & df$CBD3==0 & 
          df$CBD4==0 & df$CBD5==0)] <- 0


#### MEDICATION NUMBERS ##
#### AHTN ###

df <- df %>% 
  group_by(NACCID) %>%
  mutate(NAHTN = sum(NACCAHTN, na.rm=TRUE))

df$NACCAHTN[(df$NAHTN==0)]<- 0
df$NACCAHTN[(df$NAHTN==1)]<- 9 
df$NACCAHTN[(df$NAHTN>=2)]<- 1 
df <- as.data.frame(df)

#### LIPL ###

df <- df %>% 
  group_by(NACCID) %>%
  mutate(NLIPL = sum(NACCLIPL, na.rm=TRUE))

df$NACCLIPL[(df$NLIPL==0)]<- 0
df$NACCLIPL[(df$NLIPL==1)]<- 9 
df$NACCLIPL[(df$NLIPL>=2)]<- 1 
df <- as.data.frame(df)

#### NSD ###

df <- df %>% 
  group_by(NACCID) %>%
  mutate(NNSD = sum(NACCNSD, na.rm=TRUE))

df$NACCNSD[(df$NNSD==0)]<- 0
df$NACCNSD[(df$NNSD==1)]<- 9 
df$NACCNSD[(df$NNSD>=2)]<- 1 
df <- as.data.frame(df)

#### AC ###

df <- df %>% 
  group_by(NACCID) %>%
  mutate(NAC = sum(NACCAC, na.rm=TRUE))

df$NACCAC[(df$NAC==0)]<- 0
df$NACCAC[(df$NAC==1)]<- 9 
df$NACCAC[(df$NAC>=2)]<- 1 
df <- as.data.frame(df)

#### ADEP ###

df <- df %>% 
  group_by(NACCID) %>%
  mutate(NADEP = sum(NACCADEP, na.rm=TRUE))

df$NACCADEP[(df$NADEP==0)]<- 0
df$NACCADEP[(df$NADEP==1)]<- 9 
df$NACCADEP[(df$NADEP>=2)]<- 1 
df <- as.data.frame(df)

#### APSY ###

df <- df %>% 
  group_by(NACCID) %>%
  mutate(NAPSY = sum(NACCAPSY, na.rm=TRUE))

df$NACCAPSY[(df$NAPSY==0)]<- 0
df$NACCAPSY[(df$NAPSY==1)]<- 9 
df$NACCAPSY[(df$NAPSY>=2)]<- 1 
df <- as.data.frame(df)

#### AANX ###

df <- df %>% 
  group_by(NACCID) %>%
  mutate(NANX = sum(NACCAANX, na.rm=TRUE))

df$NACCAANX[(df$NANX==0)]<- 0
df$NACCAANX[(df$NANX==1)]<- 9 
df$NACCAANX[(df$NANX>=2)]<- 1 
df <- as.data.frame(df)

#### DBMD ###

df <- df %>% 
  group_by(NACCID) %>%
  mutate(NDBM = sum(NACCDBMD, na.rm=TRUE))

df$NACCDBMD[(df$NDBM==0)]<- 0
df$NACCDBMD[(df$NDBM==1)]<- 9 
df$NACCDBMD[(df$NDBM>=2)]<- 1 
df <- as.data.frame(df)

#### ADMD ###

df <- df %>% 
  group_by(NACCID) %>%
  mutate(NADM = sum(NACCADMD, na.rm=TRUE))

df$NACCADMD[(df$NADM==0)]<- 0
df$NACCADMD[(df$NADM==1)]<- 9 
df$NACCADMD[(df$NADM>=2)]<- 1 
df <- as.data.frame(df)

#### PDMD ###

df <- df %>% 
  group_by(NACCID) %>%
  mutate(NPDM = sum(NACCPDMD, na.rm=TRUE))

df$NACCPDMD[(df$NPDM==0)]<- 0
df$NACCPDMD[(df$NPDM==1)]<- 9 
df$NACCPDMD[(df$NPDM>=2)]<- 1 
df <- as.data.frame(df)

write.csv(df, file = "MDallmedscdr.csv",row.names=FALSE, na="")


################### Adjust values for antihypertensives ############################################

rm(list = ls())

library(data.table)
library(lubridate)
library(dplyr)
library(reshape)
library(naniar)

File <- file.path("C:/Users/Daman Kaur/Desktop/Drug classes DO")
setwd(File)
nacc <- fread("MDallmedscdr.csv")
df <- nacc

df <- df %>% group_by(NACCID) %>% filter(!all(is.na(NACCHTNC)))
df <- df %>% group_by(NACCID) %>% filter(!all(is.na(NACCACEI)))
df <- df %>% group_by(NACCID) %>% filter(!all(is.na(NACCAAAS)))
df <- df %>% group_by(NACCID) %>% filter(!all(is.na(NACCBETA)))
df <- df %>% group_by(NACCID) %>% filter(!all(is.na(NACCCCBS)))
df <- df %>% group_by(NACCID) %>% filter(!all(is.na(NACCDIUR)))
df <- df %>% group_by(NACCID) %>% filter(!all(is.na(NACCVASD)))
df <- df %>% group_by(NACCID) %>% filter(!all(is.na(NACCANGI)))

#### HTNC ###

df <- df %>% 
  group_by(NACCID) %>%
  mutate(HTNC = sum(NACCHTNC, na.rm=TRUE))

df$NACCHTNC[(df$HTNC==0)]<- 0
df$NACCHTNC[(df$HTNC==1)]<- 9 
df$NACCHTNC[(df$HTNC>=2)]<- 1 
df <- as.data.frame(df)

#### ACEI ###

df <- df %>% 
  group_by(NACCID) %>%
  mutate(ACEI = sum(NACCACEI, na.rm=TRUE))

df$NACCACEI[(df$ACEI==0)]<- 0
df$NACCACEI[(df$ACEI==1)]<- 9 
df$NACCACEI[(df$ACEI>=2)]<- 1 
df <- as.data.frame(df)

#### AAAS ###

df <- df %>% 
  group_by(NACCID) %>%
  mutate(AAAS = sum(NACCAAAS, na.rm=TRUE))

df$NACCAAAS[(df$AAAS==0)]<- 0
df$NACCAAAS[(df$AAAS==1)]<- 9 
df$NACCAAAS[(df$AAAS>=2)]<- 1 
df <- as.data.frame(df)

#### BETA ###

df <- df %>% 
  group_by(NACCID) %>%
  mutate(BETA = sum(NACCBETA, na.rm=TRUE))

df$NACCBETA[(df$BETA==0)]<- 0
df$NACCBETA[(df$BETA==1)]<- 9 
df$NACCBETA[(df$BETA>=2)]<- 1 
df <- as.data.frame(df)

#### CCBS ###

df <- df %>% 
  group_by(NACCID) %>%
  mutate(CCBS = sum(NACCCCBS, na.rm=TRUE))

df$NACCCCBS[(df$CCBS==0)]<- 0
df$NACCCCBS[(df$CCBS==1)]<- 9 
df$NACCCCBS[(df$CCBS>=2)]<- 1 
df <- as.data.frame(df)

#### DIUR ###

df <- df %>% 
  group_by(NACCID) %>%
  mutate(DIUR = sum(NACCDIUR, na.rm=TRUE))

df$NACCDIUR[(df$DIUR==0)]<- 0
df$NACCDIUR[(df$DIUR==1)]<- 9 
df$NACCDIUR[(df$DIUR>=2)]<- 1 
df <- as.data.frame(df)

#### VASD ###

df <- df %>% 
  group_by(NACCID) %>%
  mutate(VASD = sum(NACCVASD, na.rm=TRUE))

df$NACCVASD[(df$VASD==0)]<- 0
df$NACCVASD[(df$VASD==1)]<- 9 
df$NACCVASD[(df$VASD>=2)]<- 1 
df <- as.data.frame(df)

#### ANGI ###

df <- df %>% 
  group_by(NACCID) %>%
  mutate(ANGI = sum(NACCANGI, na.rm=TRUE))

df$NACCANGI[(df$ANGI==0)]<- 0
df$NACCANGI[(df$ANGI==1)]<- 9 
df$NACCANGI[(df$ANGI>=2)]<- 1 
df <- as.data.frame(df)

write.csv(df, file = "MDhypercdr.csv",row.names=FALSE, na="")
