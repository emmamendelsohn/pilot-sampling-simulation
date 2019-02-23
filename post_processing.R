library(tidyverse)
library(here)
h <- here::here

sheets <- readxl::excel_sheets(h("C1148 Input files for Match 3 ISM_Prediction Intervals.xlsx"))
dat <- c()

for(sh in sheets){
  tmp <- readxl::read_excel(h("C1148 Input files for Match 3 ISM_Prediction Intervals.xlsx"), sheet = sh)
  tmp$scenario <- sh
  dat <- rbind(dat, tmp)
}

dat$method <- ifelse(grepl("students", dat$scenario,  ignore.case = T), "Student's t", "Chebyshev")
dat$cv <- str_sub(dat$scenario, -3, -1)

pdf(paste0(h("figures", "Figure_percentSampled_ISM_withCI.pdf")), paper="USr", width=10, height=7.25)
par(omi=c(1,0.5, 0, 2.5), xpd=NA)

for (sc in unique(dat$scenario)){
  
  cdat <- dat %>% filter(scenario==sc) %>%
    rename(prob = `P(Total DU>AL)`, pilot = `Pilot Study % of Area`) %>%
    mutate(prob = prob*100, pilot=pilot*100)
  
  mod <- lm(log(prob) ~ log(pilot), data=cdat)
  newdf <- data.frame("pilot"=seq(min(cdat$pilot), max(cdat$pilot), len=1000))
  
  plines <- data.frame(exp(predict(mod, newdf, interval = "prediction", level = 0.95))) %>% rename(p_lwr=lwr, p_upr=upr)
  clines <- data.frame(exp(predict(mod, newdf, interval = "confidence", level = 0.95))) %>% rename(c_lwr=lwr, c_upr=upr)
  outlines <- bind_cols(newdf, full_join(plines, clines))
  
  p <- ggplot() +
    geom_ribbon(data=outlines, aes(x= pilot, ymin = p_lwr, ymax = p_upr), fill = "grey70", alpha=0.7)+
    geom_ribbon(data=outlines, aes(x= pilot, ymin = c_lwr, ymax = c_upr), fill = "grey50", alpha=0.7)+
    geom_point(data=cdat, aes(x=pilot, y=prob), color = "dodgerblue3") +
    geom_line(data=outlines, aes(x=pilot, y=fit), color="red")+
    labs(y = "Percent of Site > AL", 
         x = "Pilot Study Percent of Site",
         title = paste0("ISM Simulation with ", unique(cdat$method), " UCL and CV=", unique(cdat$cv), "\nMatch 3 errors (", nrow(cdat)/10, "%)"))+
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5, size=16),
          axis.text = element_text(size=14),
          axis.title = element_text(size=15))
  
  print(p)
  
}
dev.off()

##########################
dat <- c()

for(cv in c("0.5", "1.5", "3.0")){
  tmp <- read.csv(h(paste0("output/SimulationDecision_CV", cv, ".csv")))
  tmp$cv <- cv
  dat <- rbind(dat, tmp)
}

pdf(paste0(h("figures", "Figure_percentSampled_Discrete_withCI.pdf")), paper="USr", width=10, height=7.25)
par(omi=c(1,0.5, 0, 2.5), xpd=NA)

for (cv_i in unique(dat$cv)){
  
  cdat <- dat %>% filter(cv==cv_i, match==3) %>%
    rename(prob = prop_true_mean_exceeds, pilot = pilot_du_count) %>%
    mutate(prob = prob*100)
  
  mod <- lm(log(prob) ~ log(pilot), data=cdat)
  newdf <- data.frame("pilot"=seq(min(cdat$pilot), max(cdat$pilot), len=1000))
  
  plines <- data.frame(exp(predict(mod, newdf, interval = "prediction", level = 0.95))) %>% rename(p_lwr=lwr, p_upr=upr)
  clines <- data.frame(exp(predict(mod, newdf, interval = "confidence", level = 0.95))) %>% rename(c_lwr=lwr, c_upr=upr)
  outlines <- bind_cols(newdf, full_join(plines, clines))
  
  p <- ggplot() +
    geom_ribbon(data=outlines, aes(x= pilot, ymin = p_lwr, ymax = p_upr), fill = "grey70", alpha=0.7)+
    geom_ribbon(data=outlines, aes(x= pilot, ymin = c_lwr, ymax = c_upr), fill = "grey50", alpha=0.7)+
    geom_point(data=cdat, aes(x=pilot, y=prob), color = "dodgerblue3") +
    geom_line(data=outlines, aes(x=pilot, y=fit), color="red")+
    labs(y = "Percent of Site > AL", 
         x = "Pilot Study Percent of Site",
         title = paste0("Discrete Simulation with CV=", cv_i, "\nMatch 3 errors (", nrow(cdat)/10, "%)"))+
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5, size=16),
          axis.text = element_text(size=14),
          axis.title = element_text(size=15))
  
  print(p)
  
}
dev.off()

