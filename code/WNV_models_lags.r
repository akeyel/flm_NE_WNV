
# fit and compare models with lagged variables 

library(dplyr)
library(stringr)
library(tidyr)
library(lubridate)

source("code/assembleData_lags_Feb18.r") 

library(gamm4)
library(lme4)
library(mgcv)
library(MuMIn)
library(broom)
library(purrr)

# Adjust data & output file at bottom

#create contrast coefficients for county
allLags$County <- as.factor(allLags$County)
csco <- length(unique(allLags$County))
contrasts(allLags$County) = contr.sum(csco)

# Separate most recent year as out-of-sample data
oosy <- max(allLags$year) 
allLags <- allLags[allLags$year >= 2002,] # lagged variables from earlier years are still there
allLagsT <- allLags[allLags$year < oosy,] # training data
allLagsO <- allLags[allLags$year == oosy,] # out-of-sample data
yrmin <- min(allLagsT$year)
yrmax <- max(allLagsT$year)

# and create contrast coefficient for year
allLagsT$year <- as.factor(allLagsT$year)

csyr <- length(unique(allLagsT$year))
contrasts(allLagsT$year) = contr.sum(csyr)

# no lagged variables

models <- list(cases ~ County + year + offset(log(pop100K)),
               cases ~ CI + County + year + offset(log(pop100K)))

fits <- lapply(models, gam, data=allLagsT, family=nb())

results1 <- map_df(fits, glance)

results1$formula <- map_chr(fits, ~as.character(formula(.x))[3])

write.csv(results1, "temp/results1.csv", row.names = FALSE)

#one lagged variable

tlag = c(12, 18, 24, 30, 36)

models <- c("cases ~ s(lags_tmean%d, by=tmean%d) + County + year + offset(log(pop100K))",
            "cases ~ s(lags_tmean%d, by=tmean%d) + CI + County + year + offset(log(pop100K))",

            "cases ~ s(lags_ppt%d, by=ppt%d) + County + year + offset(log(pop100K))",
            "cases ~ s(lags_ppt%d, by=ppt%d) + CI + County + year + offset(log(pop100K))",

            "cases ~ s(lags_spi%d, by=spi%d) + County + year + offset(log(pop100K))",
            "cases ~ s(lags_spi%d, by=spi%d) + CI + County + year + offset(log(pop100K))",

            "cases ~ s(lags_spei%d, by=spei%d) + County + year + offset(log(pop100K))",
            "cases ~ s(lags_spei%d, by=spei%d) + CI + County + year + offset(log(pop100K))")

allmods_list <- map(tlag,
                    ~sprintf(models, .x, .x))
allmods <- flatten(allmods_list)

allfits <- map(allmods, ~gam(as.formula(.x), data=allLagsT, family=nb()))

results2 <- map_df(allfits,glance)
results2$formula <- map_chr(allfits, ~as.character(formula(.x))[3])

write.csv(results2, "temp/results2.csv", row.names = FALSE)

#temp and ppt, temp & spi, temp & spei


models <- c("cases ~ s(lags_tmean%d, by=tmean%d) + s(lags_ppt%d, by=ppt%d) + County + year + offset(log(pop100K))",
            "cases ~ s(lags_tmean%d, by=tmean%d) + s(lags_ppt%d, by=ppt%d) + CI + County + year + offset(log(pop100K))",

            "cases ~ s(lags_tmean%d, by=tmean%d) + s(lags_spi%d, by=spi%d) + County + year + offset(log(pop100K))",
            "cases ~ s(lags_tmean%d, by=tmean%d) + s(lags_spi%d, by=spi%d) + CI + County + year + offset(log(pop100K))",


            "cases ~ s(lags_tmean%d, by=tmean%d) + s(lags_spei%d, by=spei%d) + County + year + offset(log(pop100K))",
            "cases ~ s(lags_tmean%d, by=tmean%d) + s(lags_spei%d, by=spei%d) + CI + County + year + offset(log(pop100K))")

lag_combinations <- crossing(tlag = c(12, 18, 24, 30, 36), slag = seq(12, 36, 6))

allmods_list <- map2(lag_combinations$tlag, lag_combinations$slag,
                     ~sprintf(models, .x, .x, .y, .y))
allmods <- flatten(allmods_list)
allfits <- map(allmods, ~gam(as.formula(.x), data=allLagsT, family=nb()))
#map_dbl(allfits, AIC)
results3 <- map_df(allfits,glance)
results3$formula <- map_chr(allfits, ~as.character(formula(.x))[3])

write.csv(results3, "temp/results3.csv", row.names = FALSE)

# Same but w/o year

models <- list(cases ~ County + offset(log(pop100K)),
               cases ~ CI + County + offset(log(pop100K)))

fits <- lapply(models, gam, data=allLagsT, family=nb())

results4 <- map_df(fits, glance)

results4$formula <- map_chr(fits, ~as.character(formula(.x))[3])

write.csv(results4, "temp/results4.csv", row.names = FALSE)


#one lagged variable

tlag = c(12, 18, 24, 30, 36)

models <- c("cases ~ s(lags_tmean%d, by=tmean%d) + County + offset(log(pop100K))",
            "cases ~ s(lags_tmean%d, by=tmean%d) + CI + County + offset(log(pop100K))",

            "cases ~ s(lags_ppt%d, by=ppt%d) + County + offset(log(pop100K))",
            "cases ~ s(lags_ppt%d, by=ppt%d) + CI + County + offset(log(pop100K))",

            "cases ~ s(lags_spi%d, by=spi%d) + County + offset(log(pop100K))",
            "cases ~ s(lags_spi%d, by=spi%d) + CI + County + offset(log(pop100K))",

            "cases ~ s(lags_spei%d, by=spei%d) + County + offset(log(pop100K))",
            "cases ~ s(lags_spei%d, by=spei%d) + CI + County + offset(log(pop100K))")

allmods_list <- map(tlag,
                    ~sprintf(models, .x, .x))
allmods <- flatten(allmods_list)

allfits <- map(allmods, ~gam(as.formula(.x), data=allLagsT, family=nb()))

results5 <- map_df(allfits,glance)
results5$formula <- map_chr(allfits, ~as.character(formula(.x))[3])

write.csv(results5, "temp/results5.csv", row.names = FALSE)

# temp and ppt, temp & spi, temp & spei

models <- c("cases ~ s(lags_tmean%d, by=tmean%d) + s(lags_ppt%d, by=ppt%d) + County + offset(log(pop100K))",
            "cases ~ s(lags_tmean%d, by=tmean%d) + s(lags_ppt%d, by=ppt%d) + CI + County + offset(log(pop100K))",

            "cases ~ s(lags_tmean%d, by=tmean%d) + s(lags_spi%d, by=spi%d) + County + offset(log(pop100K))",
            "cases ~ s(lags_tmean%d, by=tmean%d) + s(lags_spi%d, by=spi%d) + CI + County + offset(log(pop100K))",

            "cases ~ s(lags_tmean%d, by=tmean%d) + s(lags_spei%d, by=spei%d) + County + offset(log(pop100K))",
            "cases ~ s(lags_tmean%d, by=tmean%d) + s(lags_spei%d, by=spei%d) + CI + County + offset(log(pop100K))")

lag_combinations <- crossing(tlag = c(12, 18, 24, 30, 36), slag = seq(12, 36, 6))

allmods_list <- map2(lag_combinations$tlag, lag_combinations$slag,
                     ~sprintf(models, .x, .x, .y, .y))
allmods <- flatten(allmods_list)
allfits <- map(allmods, ~gam(as.formula(.x), data=allLagsT, family=nb()))

#map_dbl(allfits, AIC)
results6 <- map_df(allfits,glance)
results6$formula <- map_chr(allfits, ~as.character(formula(.x))[3])

write.csv(results6, "temp/results6.csv", row.names = FALSE)

results <- bind_rows(list(results1, results2, results3, results4, results5, results6))

write.csv(results, "temp/AICresults.csv", row.names = FALSE)



