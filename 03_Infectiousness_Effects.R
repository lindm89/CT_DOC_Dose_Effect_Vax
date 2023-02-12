## Maggie Lind
## Feb 10 2023
## Code estimates the effect an index cases vaccination or prior infection history has on transmission within correctional facility residents 

## Load packages
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

## Immunity Event Analysis: INFECTIOUSNESS ########################################################
  ## VACCINATION #### 
    ## run model 
    survdata_vax <- subset(survdata, exposure ==1) %>% mutate(exposure = as.numeric(as.character(exposure)))
    res.cox <- coxph(Surv(time, outcome) ~ vax_status_inf*exposure + priorinf_inf + vax_status + priorinf + ns(Age, df = 4) + 
                       Race + ns(ymd(start_time), df = 4) + roomsize +  
                       strata(stfactor)  + cluster(clusterid), data =  survdata, method = "breslow")

    ## Save results
      r <- summary(res.cox)
      rc <- r$coefficients
      
      ## vax effect cell
      exp1 <- summary(glht(res.cox, linfct = c("vax_status_infVaccinated = 0")))
      exp1ci <- confint(glht(res.cox, linfct = c("vax_status_infVaccinated = 0")))
      exp1 <- data.frame(coef = exp1$test$coefficients, se = exp1$test$sigma, lower = exp1ci$confint[2], upper = exp1ci$confint[3], pvalue = "-")
      
      ## vax effect cellblock
      exp2 <- summary(glht(res.cox, linfct = c("vax_status_infVaccinated + vax_status_infVaccinated:exposure2 = 0")))
      exp2ci <- confint(glht(res.cox, linfct = c("vax_status_infVaccinated + vax_status_infVaccinated:exposure2 = 0")))
      exp2 <- data.frame(coef = exp2$test$coefficients, se = exp2$test$sigma, lower = exp2ci$confint[2], upper = exp2ci$confint[3], pvalue = rc[nrow(rc),6])
      
      ## Format the saved results
      rt <- data.frame(exposure = seq(1,2), seed = seed, iter = iteration, 
                       analysis = "Inf Immunity Events", group = "vax strat")
      for(ums in 1:2){
        dw <- get(paste("exp", ums, "ci", sep = ""))
        dwse <- get(paste("exp", ums, sep = ""))
        de <- get(paste("exp", ums, sep = ""))
        rt$fit[ums] <- dw$confint[1]
        rt$se[ums] <- dwse$se
        rt$point[ums] <- exp(dw$confint)[1]
        rt$lower[ums] <- exp(dw$confint[2])
        rt$upper[ums] <- exp(dw$confint[3])
        rt$pvalue[ums] <- de$pvalue
        rt$number_vax_0[ums] <- nrow(survdata[which(survdata$exposure == ums & survdata$vax_status == "Unvaccinated"),])
        rt$number_vax_0_cases[ums] <- nrow(survdata[which(survdata$exposure == ums & survdata$vax_status == "Unvaccinated" & survdata$outcome == 1),])
        rt$mean_time_vax_0[ums] <- mean(survdata$time[which(survdata$exposure == ums & survdata$vax_status == "Unvaccinated")])
        rt$sd_time_vax_0[ums] <- sd(survdata$time[which(survdata$exposure == ums & survdata$vax_status == "Unvaccinated")])
        rt$number_vax_1[ums] <- nrow(survdata[which(survdata$exposure == ums & survdata$vax_status == "Vaccinated"),])
        rt$number_vax_1_cases[ums] <- nrow(survdata[which(survdata$exposure == ums & survdata$vax_status == "Vaccinated" & survdata$outcome == 1),])
        rt$mean_time_vax_1[ums] <- mean(survdata$time[which(survdata$exposure == ums & survdata$vax_status == "Vaccinated")])
        rt$sd_time_vax_1[ums] <- sd(survdata$time[which(survdata$exposure == ums & survdata$vax_status == "Vaccinated")])
      }
      
      vax_results <- rt
      rm(rt)
      
      ## get the scholfelds
      test.ph_vax <- cox.zph(res.cox)$table
    
  ## PRIOR INFECTION #######################
    ## run model 
    survdata_vax <- subset(survdata, exposure == 1) %>% mutate(exposure = as.numeric(as.character(exposure)), priorinf_inf = as.numeric(priorinf_inf))
    
    res.cox <- coxph(Surv(time, outcome) ~ priorinf_inf*exposure + vax_status_inf + vax_status + priorinf + ns(Age, df = 4) + 
                       Race + ns(ymd(start_time), df = 4) + roomsize + 
                       strata(LOC_FACILITY_CODE)  + cluster(clusterid), data =  survdata, method = "breslow", iter =100)
  
    ## Save results
      r <- summary(res.cox)
      rc <- r$coefficients
      
      ## vax effect cell
      exp1 <- summary(glht(res.cox, linfct = c("priorinf_inf = 0")))
      exp1ci <- confint(glht(res.cox, linfct = c("priorinf_inf = 0")))
      exp1 <- data.frame(coef = exp1$test$coefficients, se = exp1$test$sigma, lower = exp1ci$confint[2], upper = exp1ci$confint[3], pvalue = "-")
      
      ## vax effect cellblock
      exp2 <- summary(glht(res.cox, linfct = c("priorinf_inf + priorinf_inf:exposure2 = 0")))
      exp2ci <- confint(glht(res.cox, linfct = c("priorinf_inf + priorinf_inf:exposure2 = 0")))
      exp2 <- data.frame(coef = exp2$test$coefficients, se = exp2$test$sigma, lower = exp2ci$confint[2], upper = exp2ci$confint[3], pvalue = rc[nrow(rc),6])
      
      ## save regression results HERE 
      rt <- data.frame(exposure = seq(1,2), seed = seed, iter = iteration, 
                       analysis = "Inf Immunity Events", group = "vax strat")
      for(ums in 1:2){
        dw <- get(paste("exp", ums, "ci", sep = ""))
        dwse <- get(paste("exp", ums, sep = ""))
        de <- get(paste("exp", ums, sep = ""))
        rt$fit[ums] <- dw$confint[1]
        rt$se[ums] <- dwse$se
        rt$point[ums] <- exp(dw$confint)[1]
        rt$lower[ums] <- exp(dw$confint[2])
        rt$upper[ums] <- exp(dw$confint[3])
        rt$pvalue[ums] <- de$pvalue
        rt$number_vax_0[ums] <- nrow(survdata[which(survdata$exposure == ums & survdata$vax_status == "Unvaccinated"),])
        rt$number_vax_0_cases[ums] <- nrow(survdata[which(survdata$exposure == ums & survdata$vax_status == "Unvaccinated" & survdata$outcome == 1),])
        rt$mean_time_vax_0[ums] <- mean(survdata$time[which(survdata$exposure == ums & survdata$vax_status == "Unvaccinated")])
        rt$sd_time_vax_0[ums] <- sd(survdata$time[which(survdata$exposure == ums & survdata$vax_status == "Unvaccinated")])
        rt$number_vax_1[ums] <- nrow(survdata[which(survdata$exposure == ums & survdata$vax_status == "Vaccinated"),])
        rt$number_vax_1_cases[ums] <- nrow(survdata[which(survdata$exposure == ums & survdata$vax_status == "Vaccinated" & survdata$outcome == 1),])
        rt$mean_time_vax_1[ums] <- mean(survdata$time[which(survdata$exposure == ums & survdata$vax_status == "Vaccinated")])
        rt$sd_time_vax_1[ums] <- sd(survdata$time[which(survdata$exposure == ums & survdata$vax_status == "Vaccinated")])
      }
  
      pinf_results <- rt
      rm(rt)
      
      ## get the scholfelds
      test.ph_inf <- cox.zph(res.cox)$table
      
      ## residuals
      resid <- rbind(test.ph_vax, test.ph_inf)

## save results 
  write.csv(vax_results, filename) ## save vax effect results
  write.csv(pinf_results, filename) ## save prior infection
  write.csv(resid, filename) ## save residuals


 