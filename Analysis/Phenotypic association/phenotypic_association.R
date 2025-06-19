#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# phenotypic associations of height/LHR and 15 CMD
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
pacman::p_load(dplyr,survival,data.table,ggplot2,scales,broom)

# data
dat0 <- readRDS(file = "/work/data/ukb/phe/height.RDS")
dat <- dat0 %>% filter(ethnic.c=="White")
# Cox regression loop
Xs <- c('height','lhr.res')
Ys <- c('af.2','vte.2','cm.2','aa.2','hf.2','ckd','hs.2','is.2','pad.2','as.2','t2dm.2','tia.2','cad.2','ht.2','mi')

X2Y <- data.frame(Model=character(),X=character(),Y=character(),N=numeric(),Event=numeric(),
         HR=numeric(),Lower=numeric(),Upper=numeric(),P=character(),stringsAsFactors=FALSE)
 
for (Y in Ys) {
  for (X in Xs) {
    print(paste("Processing:",X, "->", Y))
    dat1 <- dat %>% mutate(Y_date = dat[[paste0('icdDate_', Y)]],Y_yes = ifelse(is.na(Y_date), 0, 1),
      follow_end_day = fifelse(!is.na(Y_date), Y_date,fifelse(!is.na(date_lost),date_lost,fifelse(!is.na(date_death),date_death, as.Date("2023-12-31")))),
      follow_years = (as.numeric(follow_end_day) - as.numeric(date_attend)) / 365.25) %>% filter(follow_years > 0) %>% filter(!is.na(leg) & !is.na(height))
                            
    dat1$X <- dat1[[X]]
    surv.obj <- Surv(time = dat1$follow_years, event = dat1$Y_yes)
    main.X2Y <- coxph(surv.obj ~ scale(X) +age+sex+whr+college+income+deprivation+smoke_status+alcohol_status+
                        PC1+PC2+PC3+PC4+PC5+days_pa_mod,data=dat1)
    sens.X2Y <- coxph(surv.obj ~ scale(X) +age+sex+whr+college+income+deprivation+smoke_status+alcohol_status+
                        PC1+PC2+PC3+PC4+PC5+days_pa_mod+bb_HDL+bb_LDL+bb_TC+bb_TG,data=dat1)
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
# stratified association
Y <- 'af.2' # set AF as example
dat1 <- dat %>% mutate(Y_date = dat[[paste0('icdDate_', Y)]],Y_yes = ifelse(is.na(Y_date), 0, 1),
      follow_end_day = fifelse(!is.na(Y_date), Y_date,fifelse(!is.na(date_lost),date_lost,fifelse(!is.na(date_death),date_death, as.Date("2023-12-31")))),
      follow_years = (as.numeric(follow_end_day) - as.numeric(date_attend)) / 365.25) %>% filter(follow_years > 0) %>% filter(!is.na(leg) & !is.na(height))

dat1$X <- dat1[[X]]
surv.obj <- Surv(time = dat1$follow_years, event = dat1$Y_yes)
stratified_hr <- function(var) {
  z <- scale(dat1[[var]])
  g <- cut(z, breaks = c(-Inf, -1.5, -0.5, 0.5, 1.5, Inf), labels = c("-2", "-1", "0", "1", "2"))
  g <- relevel(g, ref = "0")
  dat1$g <- g
  m <- coxph(surv.obj ~ g +age+sex+whr+college+income+deprivation+smoke_status+alcohol_status+
                        PC1+PC2+PC3+PC4+PC5+days_pa_mod+bb_HDL+bb_LDL+bb_TC+bb_TG, data = dat1)
  broom::tidy(m, exponentiate = TRUE) %>%
    filter(grepl("g", term)) %>%
    left_join(exp(confint(m)) %>% as.data.frame() %>% rownames_to_column("term"),by = "term") %>%
    mutate(group = as.numeric(gsub("g", "", term)),exposure = var) %>%
    select(exposure, group, HR = estimate, Lower = `2.5 %`, Upper = `97.5 %`) %>%
    add_row(exposure = var, group = 0, HR = 1, Lower = 1, Upper = 1)
}

stratified_res <- bind_rows(stratified_hr("height"),stratified_hr("lhr.res"))
stratified_res$exposure <- factor(stratified_res$exposure,levels = c("height", "lhr.res"),labels = c("Height", "LHR"))

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