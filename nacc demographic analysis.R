rm(list = ls())
library(data.table)
library(lubridate)
library(dplyr)
library(reshape)
library(naniar)

File <- file.path("C:/Users/Daman Kaur/Desktop/Drug classes DO")
setwd(File)
df <- fread("MDallmedscdr.csv")
#df <- fread("NtoMcdr.csv")

#dh <- which(df$NACCID %in% nacc$NACCID)
#data <- df[ c(dh), ]
#data <- as.data.frame(data)

#write.csv(data, file = "NMc.csv",row.names=FALSE, na="")
data <- nacc
data <- data[data$RACE < 50, ] 

###############################################
stbl <- data %>% filter(Diagnosis == 0)

stbl <- stbl %>%
  mutate(date_visit = as.Date(Date, "%d/%m/%Y"))

ds <- stbl %>% 
  group_by(NACCID) %>%
  summarise(first_visit = min(date_visit), 
            last_visit = max(date_visit), 
            exposure = last_visit - first_visit)

ds$yrs <- (ds$exposure/365.2422)
mean(ds$yrs)
sd(ds$yrs)

############################################

df <- df[df$NEW_VISITNUM == 1, ] 
df <- df[df$RACE < 50, ] 
df <- df[df$NACCAGE >= 40, ] 

#stbl <- df %>% filter(Diagnosis == 0)
prog <- df %>% filter(Diagnosis == 1)

prog <- prog[,c(138:147)]
prog[prog==0]<-NA

xy <- as.data.frame(colMeans(prog, na.rm = TRUE))
write.excel <- function(x,row.names=FALSE,col.names=TRUE,...) {
  write.table(x,"clipboard",sep="\t",row.names=row.names,col.names=col.names,...)}

write.excel(xy)

sd(prog$NAHTN, na.rm = TRUE)
sd(prog$NLIPL, na.rm = TRUE)
sd(prog$NNSD, na.rm = TRUE)
sd(prog$NAC, na.rm = TRUE)
sd(prog$NADEP, na.rm = TRUE)
sd(prog$NAPSY, na.rm = TRUE)
sd(prog$NANX, na.rm = TRUE)
sd(prog$NDBM, na.rm = TRUE)
sd(prog$NADM, na.rm = TRUE)
sd(prog$NPDM, na.rm = TRUE)

#normality
#shapiro.test(stbl$NACCAGE)
#shapiro.test(prog$NACCAGE)


#mean and SD
mean(stbl$NACCAGE)
sd(stbl$NACCAGE)
mean(prog$NACCAGE)
sd(prog$NACCAGE)
t.test(stbl$NACCAGE,prog$NACCAGE)
wilcox.test(stbl$NACCAGE,prog$NACCAGE) 

mean(stbl$ALVST)
sd(stbl$ALVST)
mean(prog$ALVST)
sd(prog$ALVST)

mean(stbl$EDUC)
sd(stbl$EDUC)
mean(prog$EDUC)
sd(prog$EDUC)
t.test(stbl$EDUC,prog$EDUC)
wilcox.test(stbl$EDUC,prog$EDUC) 

mean(stbl$SMOKE, na.rm = TRUE)
sd(stbl$SMOKE, na.rm = TRUE)
mean(prog$SMOKE, na.rm = TRUE)
sd(prog$SMOKE, na.rm = TRUE)
t.test(stbl$SMOKE,prog$SMOKE)
wilcox.test(stbl$SMOKE,prog$SMOKE) 


## marry / fam

stbl1 <- stbl %>% filter(MARI == 1)
stbl1 <- stbl %>% filter(FAM == 1)

stbl1 <- prog %>% filter(MARI == 1)
stbl1 <- prog %>% filter(FAM == 1)
