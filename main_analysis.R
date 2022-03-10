################################################################################
################################################################################
##
## Code for main analyses in CITS and incidence of infectious diseases
##
################################################################################
################################################################################

##------ LOAD PACKAGES REQUIRED ------
library(dlnm);library(MASS);library(mixmeta);library(scales);library(splines);
library(dplyr)

##------ LOAD Example DATA ------
load("example.RData") ## alldata

##------ GENERATE PROVINCE AND OUTCOME LIST ------
pro_seq <- unique(alldata$prov_id)

##------ LOOP AND STORE FIRST STAGE RESULT ------
coef.stage <- vcov.stage <- vector("list",1)


dt <- alldata

for (i in seq(pro_seq)) {
  tmp <- dt[dt$prov_id==pro_seq[i],]
  tmp <- tmp[order(tmp$year,tmp$month),]
  
  tmp$time <- 1:nrow(tmp) ## LONG-TERM TREND
  
  ## SET NPI INTERVENTION INDICATOR
  tmp$covid <- 0
  tmp$covid[tmp$year==2020 & tmp$month>=2] <- 1
  tmp$covid <- factor(tmp$covid,levels = c("0","1"))
  
  formula <- as.formula(paste0("cases~covid+factor(group)+ns(tmean,df=3)+ns(rh_mean,df=3)+ns(precip,df=3)+factor(month)+time"))
  
  ## model
  fit <- glm(formula,data = tmp,family = quasipoisson,na.action = "na.exclude",offset = log(pop))
  coef.stage[[1]][i] <- coef(fit)["covid1"]
  vcov.stage[[1]][i] <- vcov(fit)["covid1","covid1"]
}


##------ SECOND STAGE: META TO GET POOLED RR AND CALCULATE AF ------
rr <- af <- vector("list",1)

metapost <- mixmeta(coef.stage[[1]]~1, vcov.stage[[1]])
fit <- summary(metapost)

## extract estimation and 95%CI
est_vec <- c(exp(fit$coefficients[1]),exp(fit$coefficients[1]-1.96*fit$coefficients[2]),
             exp(fit$coefficients[1]+1.96*fit$coefficients[2]))

rr[[m]] <- est_vec

## calculate province specific cases
sum_case <- rep(NA, length(pro_seq) + 1)

for (i in seq(pro_seq)) {
  tmp <- alldata[alldata$prov_id==pro_seq[i]&alldata$year==2020&alldata$month>=2,]
  case <- tmp$cases
  sum_case[i] <- sum(case)
}

## calculate province specific avoided cases and fraction
all <- data.frame("Avoided.cases"=(1/est_vec[1]-1)*sum_case[1:(length(sum_case)-1)],
                  "Avoided.cases.upper"=(1/est_vec[2]-1)*sum_case[1:(length(sum_case)-1)],
                  "Avoided.cases.lower"=(1/est_vec[3]-1)*sum_case[1:(length(sum_case)-1)])

all[nrow(all) + 1, ] <- apply(all, 2, sum)
sum_case[length(sum_case)] <- sum(sum_case[1:(length(sum_case)-1)])

all$total_cases <- sum_case
all$percent <- all$Avoided.cases / (all$total_cases+all$Avoided.cases) * 100
all$percent_low <- all$Avoided.cases.lower / (all$total_cases+all$Avoided.cases.lower) * 100
all$percent_high <- all$Avoided.cases.upper / (all$total_cases+all$Avoided.cases.upper) * 100


## re-format the output ####
all$prov_id <- c(pro_seq,NA)
all$prov_id[nrow(all)] <- "99"
all$prov_id <- as.numeric(all$prov_id)
index <- unique(alldata[,c("prov_id","prov_name")])
all <- left_join(all,index)
all$prov_name[nrow(all)] <- "Nation"

all$Case_output <- paste0(formatC(all$Avoided.cases,digits = 0,format = "f")," (",
                        formatC(all$Avoided.cases.lower,digits = 0,format = "f"),", ",
                        formatC(all$Avoided.cases.upper,digits = 0,format = "f"),")")


all$Fraction_output <- paste0(formatC(all$percent,digits = 2,format = "f")," (",
                            formatC(all$percent_low,digits = 2,format = "f"),", ",
                            formatC(all$percent_high,digits = 2,format = "f"),")")

## Final output
Final_table <- all[,c(9,8,10,11)]

## Avoided percentage
Ap <- all[all$prov_name=="Nation",c(5,6,7,11)]






