#---------------------------------------
# phenotypic associations of height/LHR and 15 CMD
#---------------------------------------
pacman::p_load(dplyr,survival,data.table,ggplot2,scales,broom)

# data
dat0 <- readRDS(file = "/work/data/ukb/phe/height.RDS")
dat <- dat0 %>% filter(ethnic_cat=="White")

#---------------------------------------
# 1.main phenotypic associations of height/LHR and 15 CMD
#---------------------------------------
Xs <- c('height','LHR')
Ys <- c('af.2','vte.2','cm.2','aa.2','hf.2','ckd','hs.2','is.2','pad.2','as.2','t2dm.2','tia.2','cad.2','ht.2','mi')

X2Y <- data.frame(Model=character(),X=character(),Y=character(),N=numeric(),Event=numeric(),
         HR=numeric(),Lower=numeric(),Upper=numeric(),P=character(),stringsAsFactors=FALSE)
 
for (Y in Ys) {
  for (X in Xs) {
    print(paste("Processing:",X, "->", Y))
    dat1 <- dat %>% mutate(
      Y_date = dat[[paste0('icdDate_', Y)]],
      Y_yes = ifelse(is.na(Y_date), 0, 1),
      follow_end_day = fifelse(!is.na(Y_date), Y_date,fifelse(!is.na(date_lost), date_lost, 
                       fifelse(!is.na(date_death), date_death, as.Date("2023-12-31")))),
      follow_years = (as.numeric(follow_end_day) - as.numeric(date_attend)) / 365.25
    ) %>% filter(follow_years > 0) %>% filter(!is.na(leg) & !is.na(height))
                            
    dat1$X <- dat1[[X]]
    surv.obj <- Surv(time = dat1$follow_years, event = dat1$Y_yes)
    main.X2Y <- coxph(surv.obj - scale(X) +age+sex+whr+college+income+deprivation+
                        smoke_status+alcohol_status+PC1+PC2+PC3+PC4+PC5+days_pa_mod,data=dat1)
    sens.X2Y <- coxph(surv.obj - scale(X) +age+sex+whr+college+income+deprivation+
                        smoke_status+alcohol_status+PC1+PC2+PC3+PC4+PC5+days_pa_mod+bb_HDL+bb_LDL+bb_TC+bb_TG,data=dat1)
    extract_results <- function(model, model_name) {
      coef_summary <- summary(model)$coefficients
      x_row <- coef_summary[grepl("X", rownames(coef_summary)), ]
      sample_size <- model$n;event_count <- model$nevent
      exp_coef <- round(exp(x_row["coef"]), 3)
      lower_95 <- round(exp(x_row["coef"] - 1.96 * x_row["se(coef)"]), 3)
      upper_95 <- round(exp(x_row["coef"] + 1.96 * x_row["se(coef)"]), 3)
      p_value <- formatC(x_row["Pr(>|z|)"], format = "e", digits = 2)
      X2Y <<- rbind(X2Y,data.frame(
          Model=model_name,X=X,Y=Y,N=sample_size,Event=event_count,HR=exp_coef,
          Lower=lower_95,Upper=upper_95,P=p_value,stringsAsFactors=FALSE))
    }
  extract_results(main.X2Y, "main.X2Y"); extract_results(sens.X2Y, "sens.X2Y")
 }
}
write.csv(X2Y, "phenotypic_association.csv", row.names = FALSE)

#---------------------------------------
# 2.grouped phenotypic associations of height/LHR and 15 CMD (3*3)
#---------------------------------------
group_labels <- expand.grid(height = 1:3, lhr = 1:3) %>%
  mutate(label = paste0("H", height, "_L", lhr)) %>%
  pull(label)
results_list <- list()

# loop
for (Y in Ys) {
  dat1 <- dat %>%
    mutate(
      height.cat = cut(height, breaks = quantile(height, probs = c(0, 0.2, 0.8, 1), na.rm = TRUE),
                       labels = c(1, 2, 3), include.lowest = TRUE),
      lhr.cat = cut(LHR, breaks = quantile(LHR, probs = c(0, 0.2, 0.8, 1), na.rm = TRUE),
                       labels = c(1, 2, 3), include.lowest = TRUE),
      Y_date = .data[[paste0("icdDate_", Y)]],
      Y_yes = ifelse(is.na(Y_date), 0, 1),
      follow_end_day = fifelse(!is.na(Y_date), Y_date, fifelse(!is.na(date_lost), date_lost,
                       fifelse(!is.na(date_death), date_death, as.Date("2023-12-31")))),
      follow_years = as.numeric(follow_end_day - date_attend) / 365.25,
      group = paste0("H", height.cat, "_L", lhr.cat)) %>%
      filter(!is.na(height.cat) & !is.na(lhr.cat) & follow_years > 0)
  
  dat1$group <- relevel(factor(dat1$group), ref = "H2_L2")
  surv.obj <- Surv(dat1$follow_years, dat1$Y_yes)
  
  cox.model <- coxph(surv.obj - group +age+sex+whr+college+income+deprivation+
                       smoke_status+alcohol_status+PC1+PC2+PC3+PC4+PC5+days_pa_mod,data=dat1)
  
  # extract results
  res <- tidy(cox.model, exponentiate = TRUE, conf.int = TRUE) %>%
    filter(grepl("^group", term)) %>%
    mutate(group = gsub("^group", "", term)) %>%
    add_row(term = "groupH2_L2", group = "H2_L2", estimate = 1, conf.low = 1, conf.high = 1, p.value = NA) %>%
    mutate(hr_ci = sprintf("%.2f (%.2f–%.2f)", estimate, conf.low, conf.high),
      p = ifelse(is.na(p.value), "", formatC(p.value, format = "e", digits = 2)))
  
  # fomrmat
  res_wide <- res %>% select(group, hr_ci, p) %>%
    pivot_longer(cols = c(hr_ci, p), names_to = "type") %>%
    mutate(colname = paste0(group, ifelse(type == "hr_ci", "_HR_CI", "_P"))) %>%
    select(colname, value) %>% tibble::deframe()
  
  results_list[[Y]] <- res_wide
}

# results
plot_data <- bind_rows(results_list, .id = "Outcome") %>%
  pivot_longer(-Outcome, names_to = "Metric", values_to = "Value") %>%
  filter(str_detect(Metric, "_HR_CI")) %>%
  mutate(group = str_remove(Metric, "_HR_CI"),
    height.cat = as.integer(str_extract(group, "(?<=H)[1-3]")),
    lhr.cat = as.integer(str_extract(group, "(?<=_L)[1-3]")),
    HR = as.numeric(str_extract(Value, "^[0-9.]+")),
    CI_lower = as.numeric(str_extract(Value, "(?<=\\().+?(?=–)")),
    CI_upper = as.numeric(str_extract(Value, "(?<=–).+?(?=\\))"))) %>%
  left_join(bind_rows(results_list, .id = "Outcome") %>%
      pivot_longer(-Outcome, names_to = "Metric", values_to = "Pval") %>%
      filter(str_detect(Metric, "_P")) %>% mutate(group = str_remove(Metric, "_P")),
    by = c("Outcome", "group")) %>%
  mutate(pval = as.numeric(Pval),
    sig = case_when(pval <= 0.001 - "***",pval <= 0.01  - "**",pval <= 0.05  - "*",TRUE - ""))

write.csv(plot_data,'Source Data Fig. 3.csv')

#---------------------------------------
# 3.SD stratified phenotypic associations of height/LHR and 15 CMD 
#---------------------------------------
Y <- 'af.2' # set AF as example
    dat1 <- dat %>% mutate(
      Y_date = dat[[paste0('icdDate_', Y)]],
      Y_yes = ifelse(is.na(Y_date), 0, 1),
      follow_end_day = fifelse(!is.na(Y_date), Y_date,fifelse(!is.na(date_lost), date_lost, 
                       fifelse(!is.na(date_death), date_death, as.Date("2023-12-31")))),
      follow_years = (as.numeric(follow_end_day) - as.numeric(date_attend)) / 365.25
    ) %>% filter(follow_years > 0) %>% filter(!is.na(leg) & !is.na(height))

dat1$X <- dat1[[X]]
surv.obj <- Surv(time = dat1$follow_years, event = dat1$Y_yes)
stratified_hr <- function(var) {
  z <- scale(dat1[[var]])
  g <- cut(z, breaks = c(-Inf, -1.5, -0.5, 0.5, 1.5, Inf), labels = c("-2", "-1", "0", "1", "2"))
  g <- relevel(g, ref = "0")
  dat1$g <- g
  m <- coxph(surv.obj - g +age+sex+whr+college+income+deprivation+smoke_status+alcohol_status+
                        PC1+PC2+PC3+PC4+PC5+days_pa_mod+bb_HDL+bb_LDL+bb_TC+bb_TG, data = dat1)
  broom::tidy(m, exponentiate = TRUE) %>%
    filter(grepl("g", term)) %>%
    left_join(exp(confint(m)) %>% as.data.frame() %>% rownames_to_column("term"),by = "term") %>%
    mutate(group = as.numeric(gsub("g", "", term)),exposure = var) %>%
    select(exposure, group, HR = estimate, Lower = `2.5 %`, Upper = `97.5 %`) %>%
    add_row(exposure = var, group = 0, HR = 1, Lower = 1, Upper = 1)
}

stratified_res <- bind_rows(stratified_hr("height"),stratified_hr("LHR"))
stratified_res$exposure <- factor(stratified_res$exposure,levels = c("height", "LHR"),labels = c("Height", "LHR"))

# Supplementary Fig. 11
ggplot(stratified_res, aes(x = group, y = HR, color = exposure)) +
  geom_point(size = 1) +
  geom_line(aes(group = exposure),size=0.3) +
  geom_errorbar(aes(ymin = Lower, ymax = Upper), width = 0.15, size=0.3) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "gray50",size=0.3) +
  scale_x_continuous(breaks = -2:2, labels = paste0(-2:2)) +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.01)) +
  scale_color_manual(values = c("Height" = "#f74747", "LHR" = "#3b9dff")) +
  theme_minimal(base_size = 10) +
  theme(
    panel.grid = element_blank(),legend.position = "none",
    axis.line = element_line(color = "black", size = 0.3),
    axis.ticks = element_line(color = "black", size = 0.3),
    axis.title = element_blank())