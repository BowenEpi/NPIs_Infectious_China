################################################################################
################################################################################
##
## Code to generate example data in CITS and incidence of infectious diseases
##
################################################################################
################################################################################

## GENERATE RANDOM EXAMPLE DATA
alldata <- data.frame(prov_name=rep(c("A","B","C","D","E"),each=11*12),
                      prov_id=rep(1:5,each=11*12),
                      year=rep(rep(2010:2020,each=12),5),
                      month=rep(1:12,11*5))

set.seed(1234)
alldata$cases <- round(rpois(nrow(alldata),lambda = 200),0)
alldata$tmean <- rnorm(nrow(alldata),mean = 25,sd = 12)
alldata$precip <- rnorm(nrow(alldata),mean = 25,sd = 12)
alldata$rh_mean <- rnorm(nrow(alldata),mean = 85,sd = 12)
alldata$pop <- runif(nrow(alldata),min = 5000,max = 6000)

alldata$group <- 0
alldata$group[alldata$month==1] <- 1


save(alldata,file = "example.RData")