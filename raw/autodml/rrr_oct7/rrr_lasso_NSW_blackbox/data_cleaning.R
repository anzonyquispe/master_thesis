# prepare data

import.LaLonde.data <- function(file.name=c("nsw_control"), variable.names=names74) {
  data <- read.table(paste(file.name,".txt",sep=""))
  colnames(data) <- variable.names
  return(data)
}
names74 <- c("treat", "age", "education", "black", "hispanic", "married", "nodegree", "re74", "re75", "re78")


import.LaLonde75.data <- function(file.name=c("nsw_control"), variable.names=names75) {
  data <- read.table(paste(file.name,".txt",sep=""))
  colnames(data) <- variable.names
  return(data)
}
names75 <- c("treat", "age", "education", "black", "hispanic", "married", "nodegree", "re75", "re78")

#From Dehejia's webpage:
#treatment indicator (1 if treated, 0 if not treated), age, education, Black (1 if black, 0 otherwise), Hispanic (1 if Hispanic, 0 otherwise), married (1 if married, 0 otherwise), nodegree (1 if no degree, 0 otherwise), RE74 (earnings in 1974), RE75 (earnings in 1975), and RE78 (earnings in 1978)


nswre74_control <- import.LaLonde.data(file.name="nswre74_control",variable.names=names74)
nswre74_treated <- import.LaLonde.data(file.name="nswre74_treated",variable.names=names74)
psid_control1 <- import.LaLonde.data(file.name="psid_controls",variable.names=names74)
cps_control1<- import.LaLonde.data(file.name="cps_controls",variable.names=names74)

nswre75_control <- import.LaLonde75.data(file.name="nsw_control",variable.names=names75)
nswre75_treated <- import.LaLonde75.data(file.name="nsw_treated",variable.names=names75)


##################################
## 			 	##
## 	Prepare data 	##
## 			 	##
##################################


make_dta <- function(treated, control,name){
  df<-rbind(treated,control)
  df_export<-mutate(df,
                    age_2=age^2,
                    education_2=education^2,
                    re74_2=re74^2,
                    re75_2=re75^2,
                    u74=as.numeric(re74==0),
                    u75=as.numeric(re75==0),
                    black_u74=u74*black,
                    black_u75=u75*black,
                    hispanic_u74=u74*hispanic,
                    hispanic_u75=u75*hispanic)
  write.dta(df_export,paste0("/Users/rahul/Documents/research/rrr_lasso_NSW_blackbox/",name,"_export.dta"))
  return(df)
}

nsw<-make_dta(nswre74_treated,nswre74_control,"nsw")
psid<-make_dta(nswre74_treated,psid_control1,"psid")
cps<-make_dta(nswre74_treated,cps_control1,"cps")