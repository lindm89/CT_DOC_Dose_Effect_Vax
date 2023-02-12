## Maggie Lind
## Generate the matched dataset for the close contact analysis
  ## matching Cell and Dorm residents based on exposure (known close contact) on facility, room type, vaccination status 

library(dplyr)
library(ggplot2)
library(lubridate)
library(tidyr)
library(hrbrthemes)
library(ggridges)
library(ggsci)
library(survival)
library(forestplot)
library(AF)
library(survminer)
library(multcomp)
library(splines)

## Load data 
survdata <- read.csv(LOAD_DATA.csv)

## subset to the variant of interest (Delta or Omicron in manuscript) (Delta used as exampled here)
survdata <- subset(survdata, variant == "Delta")

## format the survival data 
survdata$exposure <- as.factor(survdata$exposure)
survdata$stfactor <- paste(survdata$Housing.Unit.Code, survdata$LOC_FACILITY_CODE)
survdata$clusterid <- paste(survdata$Housing.Unit.Code, survdata$contact_date)

## Immunity Event Analysis: Infection Risk ########################################################
  ## VACCINATION #### 
    ## run model 
      res.cox <- coxph(Surv(time, outcome) ~ exposure*vax_status + priorinf + ns(Age, df = 4) + 
                         Race + ns(ymd(start_time), df = 4) + roomsize + 
                         strata(stfactor) + cluster(clusterid), data =  survdata, method = "breslow")
      
    ## Save results
      r <- summary(res.cox)
      rc <- r$coefficients
      
      ## vax effect unexposed
      exp0 <- summary(glht(res.cox, linfct = c("vax_statusVaccinated = 0")))
      exp0ci <- confint(glht(res.cox, linfct = c("vax_statusVaccinated = 0")))
      exp0 <- data.frame(coef = exp0$test$coefficients, se = exp0$test$sigma, lower = exp0ci$confint[2], upper = exp0ci$confint[3], pvalue = "-")
      
      ## vax effect cell exposed
      exp1 <- summary(glht(res.cox, linfct = c("vax_statusVaccinated + exposure1:vax_statusVaccinated = 0")))
      exp1ci <- confint(glht(res.cox, linfct = c("vax_statusVaccinated + exposure1:vax_statusVaccinated = 0")))
      exp1 <- data.frame(coef = exp1$test$coefficients, se = exp1$test$sigma, lower = exp1ci$confint[2], upper = exp1ci$confint[3], pvalue = rc[nrow(rc)-1,6])

      ## vax effect cellblock exposed
      exp2 <- summary(glht(res.cox, linfct = c("vax_statusVaccinated + exposure2:vax_statusVaccinated = 0")))
      exp2ci <- confint(glht(res.cox, linfct = c("vax_statusVaccinated + exposure2:vax_statusVaccinated = 0")))
      exp2 <- data.frame(coef = exp2$test$coefficients, se = exp2$test$sigma, lower = exp2ci$confint[2], upper = exp2ci$confint[3], pvalue = rc[nrow(rc),6])
      
      ## Format the saved results
      rt <- data.frame(exposure = seq(0,2), seed = seed, iter = iteration, 
                       analysis = "Sus Immunity Events", group = "vax strat")
      for(ums in 1:3){
        dw <- get(paste("exp", ums-1, "ci", sep = ""))
        dwse <- get(paste("exp", ums-1, sep = ""))
        de <- get(paste("exp", ums-1, sep = ""))
        rt$fit[ums] <- dw$confint[1]
        rt$se[ums] <- dwse$se
        rt$point[ums] <- exp(dw$confint)[1]
        rt$lower[ums] <- exp(dw$confint[2])
        rt$upper[ums] <- exp(dw$confint[3])
        rt$pvalue[ums] <- de$pvalue
        rt$number_vax_0[ums] <- nrow(survdata[which(survdata$exposure == ums-1 & survdata$vax_status == "Unvaccinated"),])
        rt$number_vax_0_cases[ums] <- nrow(survdata[which(survdata$exposure == ums-1 & survdata$vax_status == "Unvaccinated" & survdata$outcome == 1),])
        rt$mean_time_vax_0[ums] <- mean(survdata$time[which(survdata$exposure == ums-1 & survdata$vax_status == "Unvaccinated")])
        rt$sd_time_vax_0[ums] <- sd(survdata$time[which(survdata$exposure == ums-1 & survdata$vax_status == "Unvaccinated")])
        rt$number_vax_1[ums] <- nrow(survdata[which(survdata$exposure == ums-1 & survdata$vax_status == "Vaccinated"),])
        rt$number_vax_1_cases[ums] <- nrow(survdata[which(survdata$exposure == ums-1 & survdata$vax_status == "Vaccinated" & survdata$outcome == 1),])
        rt$mean_time_vax_1[ums] <- mean(survdata$time[which(survdata$exposure == ums-1 & survdata$vax_status == "Vaccinated")])
        rt$sd_time_vax_1[ums] <- sd(survdata$time[which(survdata$exposure == ums-1 & survdata$vax_status == "Vaccinated")])
      }
      
      ## combine the vax results
      vax_results <- rt
      rm(rt)
      
      ## get the scholfelds
      test.ph_sus_vax <- cox.zph(res.cox)$table
      
  ## PRIOR INFECTION #######################
    ## Susceptibility ########################
    res.cox <- coxph(Surv(time, outcome) ~ exposure*priorinf + vax_status + ns(Age, df = 4) + 
                       Race + ns(ymd(start_time), df = 4) +roomsize + 
                       strata(stfactor) + cluster(clusterid), data =  survdata, method = "breslow")
      
    ## Save results
      r <- summary(res.cox)
      rc <- r$coefficients
      
      ## prior infect effect unexposed
      exp0 <- summary(glht(res.cox, linfct = c("priorinf = 0")))
      exp0ci <- confint(glht(res.cox, linfct = c("priorinf = 0")))
      exp0 <- data.frame(coef = exp0$test$coefficients, se = exp0$test$sigma, lower = exp0ci$confint[2], upper = exp0ci$confint[3], pvalue = "-")

      ## prior infect effect cell exposure
      exp1 <- summary(glht(res.cox, linfct = c("priorinf + exposure1:priorinf = 0")))
      exp1ci <- confint(glht(res.cox, linfct = c("priorinf + exposure1:priorinf = 0")))
      exp1 <- data.frame(coef = exp1$test$coefficients, se = exp1$test$sigma, lower = exp1ci$confint[2], upper = exp1ci$confint[3], pvalue = rc[nrow(rc)-1,6])

      ## prior infect effect cellblock exposure
      exp2 <- summary(glht(res.cox, linfct = c("priorinf + exposure2:priorinf = 0")))
      exp2ci <- confint(glht(res.cox, linfct = c("priorinf + exposure2:priorinf = 0")))
      exp2 <- data.frame(coef = exp2$test$coefficients, se = exp2$test$sigma, lower = exp2ci$confint[2], upper = exp2ci$confint[3], pvalue = rc[nrow(rc),6])
      
      ## Format the saved results
      rt <- data.frame(exposure = seq(0,2), seed = seed, iter = iteration, 
                       analysis = "Sus Immunity Events", group = "prior infection strat")
      for(ums in 1:3){
        dw <- get(paste("exp", ums-1, "ci", sep = ""))
        dwse <- get(paste("exp", ums-1, sep = ""))
        de <- get(paste("exp", ums-1, sep = ""))
        rt$fit[ums] <- dw$confint[1]
        rt$se[ums] <- dwse$se
        rt$point[ums] <- exp(dw$confint)[1]
        rt$lower[ums] <- exp(dw$confint[2])
        rt$upper[ums] <- exp(dw$confint[3])
        rt$pvalue[ums] <- de$pvalue
        rt$number_vax_0[ums] <- nrow(survdata[which(survdata$exposure == ums-1 & survdata$priorinf == 0),])
        rt$number_vax_0_cases[ums] <- nrow(survdata[which(survdata$exposure == ums-1 & survdata$priorinf == 0 & survdata$outcome == 1),])
        rt$mean_time_vax_0[ums] <- mean(survdata$time[which(survdata$exposure == ums-1 & survdata$priorinf == 0)])
        rt$sd_time_vax_0[ums] <- sd(survdata$time[which(survdata$exposure == ums-1 & survdata$priorinf == 0)])
        rt$number_vax_1[ums] <- nrow(survdata[which(survdata$exposure == ums-1 & survdata$priorinf == 1),])
        rt$number_vax_1_cases[ums] <- nrow(survdata[which(survdata$exposure == ums-1 & survdata$priorinf == 1 & survdata$outcome == 1),])
        rt$mean_time_vax_1[ums] <- mean(survdata$time[which(survdata$exposure == ums-1 & survdata$priorinf == 1)])
        rt$sd_time_vax_1[ums] <- sd(survdata$time[which(survdata$exposure == ums-1 & survdata$priorinf == 1)])
      }
      
      ## combine the vax results
      pinf_results <- rt
      rm(rt)
      
      ## get the scholfelds
      test.ph_sus_pinf <- cox.zph(res.cox)$table
      
      ## residuals
      resid <- rbind(test.ph, test.ph_sus_vax, test.ph_sus_pinf)

## save results 
  write.csv(vax_results, filename) ## save vax effect results
  write.csv(pinf_results, filename) ## save prior infection
  write.csv(resid, filename) ## save residuals
      