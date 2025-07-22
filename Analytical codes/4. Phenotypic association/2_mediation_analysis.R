#---------------------------------------
# Mediation analysis of height/LHR and CAD, MI, T2D, and HT
#---------------------------------------
pacman::p_load(dplyr,survival,data.table,utils,scales,stats,mediation)

# data
dat0 <- readRDS(file = "/work/data/ukb/phe/height.RDS")
dat <- dat0 %>% filter(ethnic.c=="White")

Ms <- c("bb_BUN", "bb_APOB", "bb_TP", "bc_LYMPH", "bc_MONO", "bb_ALP", "bb_CRE","bc_NEUT", "bc_WBC",
  "bb_ALB", "bb_HBA1C", "bb_CYS", "bb_APOA", "bb_CRP","bb_UA", "bb_LDL", "bb_AST2ALT", "bb_TC",
  "bb_HDL", "bb_BILD", "bb_TBIL","bb_TG", "bb_AST", "bb_GLU", "bb_ALT", "bc_EO", "bb_LPA", "bc_BASO")
Ys <- c('t2dm.2','cad.2','ht.2','mi')
mediation <- data.frame()
for (Y in Ys) {
  for (M in Ms) {
    print(paste("Processing:",M,"->",Y))
    # data process
    dat1 <- dat %>% mutate(
      Y_date = dat[[paste0('icdDate_', Y)]],
      Y_yes = ifelse(is.na(Y_date), 0, 1),
      follow_end_day = fifelse(!is.na(Y_date), Y_date,fifelse(!is.na(date_lost), date_lost, 
                       fifelse(!is.na(date_death), date_death, as.Date("2023-12-31")))),
      follow_years = (as.numeric(follow_end_day) - as.numeric(date_attend)) / 365.25
    ) %>% filter(follow_years > 0) %>% filter(!is.na(leg) & !is.na(height))
    
    dat1$M <- dat1[[M]]
    surv.obj <- Surv(time = dat1$follow_years, event = dat1$Y_yes)
    # mediation
    tryCatch({
      fit.X2Y <- survreg(surv.obj - scale(LHR)+scale(height)+age+sex+whr+college+income+deprivation+
                           smoke_status+alcohol_status+PC1+PC2+PC3+PC4+PC5+days_pa_mod,data = dat1)
      res.X2Y <- summary(fit.X2Y)$table
      fit.X2M <- lm(scale(M) - scale(LHR)+scale(height)+age+sex+whr+college+income+deprivation+
                      smoke_status+alcohol_status+PC1+PC2+PC3+PC4+PC5+days_pa_mod,data = dat1)
      res.X2M <- coef(summary(fit.X2M))
      fit.M2Y <- survreg(surv.obj - scale(LHR)+scale(M)+scale(height)+age+sex+whr+college+income+
                           deprivation+smoke_status+alcohol_status+PC1+PC2+PC3+PC4+PC5+days_pa_mod,data = dat1)
      res.M2Y <- summary(fit.M2Y)$table
      # med model
      med.lhr <- mediation::mediate(fit.X2M, fit.M2Y, treat = "scale(LHR)", mediator = "scale(M)", sims = 100)
      res.lhr <- summary(med.lhr)
      med.height <- mediation::mediate(fit.X2M, fit.M2Y, treat = "scale(height)", mediator = "scale(M)", sims = 100)
      res.height <- summary(med.height)
      
      # extract results
      mediation <- rbind(mediation,
                         data.frame(
                           X = "LHR", M = M, Y = Y,
                           MP = round(res.lhr$n.avg, 4), MP_Lower = round(res.lhr$n.avg.ci[1], 4), MP_Upper = round(res.lhr$n.avg.ci[2], 4), MP_p = formatC(res.lhr$n.avg.p, format = "e", digits = 3),
                           ME = round(res.lhr$d.avg, 4), ME_Lower = round(res.lhr$d.avg.ci[1], 4), ME_Upper = round(res.lhr$d.avg.ci[2], 4), ME_p = formatC(res.lhr$d.avg.p, format = "e", digits = 3),
                           DE = round(res.lhr$z.avg, 4), DE_Lower = round(res.lhr$z.avg.ci[1], 4), DE_Upper = round(res.lhr$z.avg.ci[2], 4), DE_p = formatC(res.lhr$z.avg.p, format = "e", digits = 3),
                           TE = round(res.lhr$tau.coef, 4), TE_Lower = round(res.lhr$tau.ci[1], 4), TE_Upper = round(res.lhr$tau.ci[2], 4), TE_p = formatC(res.lhr$tau.p, format = "e", digits = 3)
                         ),
                         data.frame(
                           X = "height", M = M, Y = Y,
                           MP = round(res.height$n.avg, 4), MP_Lower = round(res.height$n.avg.ci[1], 4), MP_Upper = round(res.height$n.avg.ci[2], 4), MP_p = formatC(res.height$n.avg.p, format = "e", digits = 3),
                           ME = round(res.height$d.avg, 4), ME_Lower = round(res.height$d.avg.ci[1], 4), ME_Upper = round(res.height$d.avg.ci[2], 4), ME_p = formatC(res.height$d.avg.p, format = "e", digits = 3),
                           DE = round(res.height$z.avg, 4), DE_Lower = round(res.height$z.avg.ci[1], 4), DE_Upper = round(res.height$z.avg.ci[2], 4), DE_p = formatC(res.height$z.avg.p, format = "e", digits = 3),
                           TE = round(res.height$tau.coef, 4), TE_Lower = round(res.height$tau.ci[1], 4), TE_Upper = round(res.height$tau.ci[2], 4), TE_p = formatC(res.height$tau.p, format = "e", digits = 3)
                         ))
    }, silent=TRUE) 
  }
}
mediation$MP_FDR <- p.adjust(mediation_df$MP_p, method = "fdr")
write.csv(mediation,'mediation.csv',row.names = F)
