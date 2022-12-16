library(simstudy)
library(lme4)
library(geepack)
library(gee)
library(ggplot2)
library(rstanarm)
library(brms)
library(sjPlot)
library(JointAI)
library(dplyr)
library(mice)
library(tidyr)
library(jomo)
library(Metrics)
# using sim study package to simulate dataset
# oredictors are not correlated
set.seed(1)
def <- defData(varname = "treatment", dist="binary", formula = 0.5)

def <- defData(def, varname="age", dist="normal", formula=25, variance = 10)

def <- defData(def, varname = "sex", dist="binary", formula = 0.5)

def <- defData(def, varname = "drinker", dist="binary", formula = 0.5)

def <- defData(def, varname = "Y0", dist = "normal", 
               formula = "20 - 0.0125*age + 0.5*drinker + 0.125*sex", variance = 5)

def <- defData(def, varname = "Y1", dist = "normal", 
               formula = "Y0 - 2.5*(treatment+0.125) + 0.75*drinker + 0.075*sex + 0.000025*age", variance = 7)

def <- defData(def, varname = "Y2", dist = "normal", 
               formula = "Y1 - 3*(treatment+0.125) + 1.25*drinker - 0.05*sex + 0.000025*age", variance = 8)

def <- defData(def, varname = "Y3", dist = "normal", 
               formula = "Y2 - 3.25*(treatment+0.125) + 1.75*drinker - 0.0025*sex - 0.00005*age", variance = 10)

# create dataframe
set.seed(18)
wide_df <- genData(30, def)


#wide to long
long_df <- addPeriods(wide_df, nPeriods = 4, idvars = "id", 
                      timevars = c("Y0", "Y1", "Y2", "Y3"), timevarName = "Y")

long_df <- long_df %>%
  dplyr::select(-timeID) %>%
  rename(time = period)

#factor 
long_df <- data.frame(long_df)
needs_factoring <- c("treatment", "sex", "drinker")
long_df[needs_factoring] <- lapply(long_df[needs_factoring] , factor)
levels(long_df$treatment) <- c("control", "intervention")
  

avg_df = data.frame(with(long_df, tapply(Y, list(treatment, time), mean, na.rm=TRUE)) %>% 
                       ftable() %>% 
                       round(2))

names(avg_df)[1] = "treatment"
names(avg_df)[2] = "time"
names(avg_df)[3] = "avg"

# avg plots
p <- ggplot(data = avg_df, aes(x = time, y = avg, group = treatment))+
  geom_point() + 
  geom_line(aes(col = treatment)) +
  geom_point() +
  ggtitle("Avg. cigarettes per week by treatment group") + 
  xlab("Week") + ylab("Cigarettes")
p

# all outcomes
p1 <- ggplot(long_df, aes(x=time, y=Y, color=treatment)) +
  geom_point(alpha=0.5)+
  #guides(color="black") +
  #geom_smooth(method=lm, se=F, color="black") +
  ggtitle("All outcomes by treatment group") + 
  xlab("Week") + ylab("Cigarettes")
p1

p2 <- ggplot(long_df, aes(x=time, y=Y, color=treatment)) +
  geom_smooth(data=long_df, aes(x=time, y=Y), #population trend
              method="lm", inherit.aes=FALSE, se=FALSE, color="black") + 
  geom_smooth(method="lm", se=F) +
  geom_point(alpha=0.5) +
  labs(x="Time", y="Cigarettes") +
  ggtitle("Linear regression") + 
  xlab("Week") + ylab("Cigarettes")
p2

# linear regression for each subject
p3 <- ggplot(long_df, aes(x=time, y=Y, color=id, group=id)) +
  geom_point(alpha=0.5)+
  guides(color=F) +
  geom_smooth(method=lm, se=F) +
  ggtitle("Linear regression for each subject") + 
  xlab("Week") + ylab("Cigarettes")
p3



gee_model <- geeglm(Y ~ time*treatment + sex + drinker + age, id=id, data=long_df, corstr = "ar1")
summary(gee_model)

mixed_model <- lmer(Y ~ time*treatment + sex + drinker + age + (1 + time | id), data=long_df)
summary(mixed_model)

model1 <- brm(Y ~ time*treatment + sex + drinker + age + (1 |id), 
                    data = long_df)

model1 <- brm(Y ~ time*treatment + sex + drinker + age + (1 + time | id), data = long_df, 
          prior = c(
            # for intercept 
            prior(normal(0, 50), class = "Intercept"), 
            # for slope
            prior(normal(0, 10), class = "b"),
            # for interaction
            prior(normal(0, 5), class = "b", 
                  coef = "time:treatmentintervention"),
            # for tau_beta0 and tau_beta1
            prior(gamma(2, 0.2), class = "sd", group = "id"),
            # for sigma
            prior(lkj(1), class = "cor")),
          seed = 18)

plot(model1)
summary(model1)
ground_truth <- summary(model1)$fixed$Estimate
  
# use regular expressions to extract the desired parameters
  # use regular expressions to extract the desired parameters
ran <- as.matrix(model1, pars = "time")
fix <- as.vector(as.matrix(model1, variable = c("sex1", "drinker1"))
comb <- ran + fix
colMeans(comb)  # matches results of coef
fake <- mean(coef(model1)[1]$id[1:30])
fake

comb <- ran + fix
colMeans(comb)  # matches results of coef
<do whatever you want with the samples>

View(coef_post)
# checking the fit
coef_post <- coef(model1, summary = FALSE)

df_lines <- tibble(id = colnames(coef_post$id), 
                   Intercept = colMeans(coef_post$id[ , , "Intercept"]), 
                   time = colMeans(coef_post$id[ , , "time"]))

p4 <- ggplot(long_df, aes(x = time, y = Y)) + 
       geom_point(size = 0.5) + 
       geom_abline(data = df_lines,  aes(intercept = Intercept, slope = time),
              color = "blue", size = 0.8, alpha = 0.5) +
       geom_smooth(se = FALSE, col = "red", size = 0.8, alpha = 0.5) + 
       facet_wrap(~ id, ncol = 6) +
       ggtitle("Posterior mean of parameters for each subject") + 
       xlab("Week") + ylab("Cigarettes")
p4

df_lines=data.frame(df_lines)
df_lines$id <- as.numeric(df_lines$id)

p5 <- ggplot(long_df, aes(x = time, y = Y, col = id)) + 
       geom_jitter(size = 0.5, width = 0.1) + 
      geom_abline(data = df_lines,
              aes(intercept = Intercept, slope = time, col = id), 
              size = 0.8, alpha = 0.5) +
      guides(col = FALSE)


bayesplot::pp_check(model1, type = "ribbon_grouped", group = "id", 
         facet_args = list(ncol = 6))
?pp_check

# poserior predictive check
pp_check(model1, nsamples = 100)
pp_check(model1, type = "stat_2d", stat = c("max", "min"))

# posterior plots with covariates of interest
mcmc_plot(model1, type = "areas", variable = c("time"), 
         prob = 0.90)

mcmc_plot(model1)
?mcmc_plot()

# posterior plots with all covariates of interest
brms::stanplot(model1, type = "areas", pars = c("time:treatment1", "sex1",
                                          "drinker1", "age", "b_time", "treatment"), prob = 0.90)

brms::mcmc_plot(model1, type = "areas", variable = c("time, time:treatment", "sex"), 
          prob = 0.90)
?mcmc_plot


#marginal effects
plot(marginal_effects(model1, effects = "time", 
                   # Request two lines using `conditions`
                   conditions = tibble(treatment = c("0", "1"))), 
  points = TRUE, point_args = list(size = 0.5))

bayes_R2(model1)


#plot showing all models
sjPlot::tab_model(model1, show.icc = FALSE, show.obs = FALSE)







########### should we account for clustering? ####################
ggplot(long_df, aes(x = id, y = Y)) + 
  geom_jitter(width = 0.1, col = "darkgrey") + 
  stat_summary(geom = "point", fun.y = mean, 
               size = 4, shape = 24, fill = "red")


ggplot(long_df, aes(x = time, y = Y)) + 
  geom_point(size = 0.5) + 
  geom_smooth() + 
  # presented by person
  facet_wrap(~ id, ncol = 6)
test_model <- brm(Y ~ (1 | id), data = long_df, 
                  prior = c(# for intercept 
                    prior(normal(0, 20), class = "Intercept"), 
                    # for tau
                    prior(gamma(2, 0.2), class = "sd"), 
                    # for sigma
                    prior(student_t(4, 0, 5), class = "sigma")), 
                  control = list(adapt_delta = .95), 
                  cores = 2L, 
                  seed = 2107)

# Computing ICC
# 1. Obtain posterior draws of tau and sigma
sd_test <- VarCorr(test_model, summary = FALSE)
draws_tau <- sd_test$id$sd[ , "Intercept"]  # tau
draws_sigma <- sd_test$residual__$sd[ , 1]  #sigma
# 2. Compute draws for ICC
draws_icc <- draws_tau^2 / (draws_tau^2 + draws_sigma^2)
# Plot the ICC
qplot(draws_icc, geom = "density", xlab = "ICC", bw = "SJ")

psych::describe(draws_icc)

broom::tidy(test_model, parameters = c("b_Intercept", "sd", "sigma")) %>% 
  knitr::kable()

d_eff <- 1 + (30-1)*(0.51)

############################################################################
# create missing data - monotone only 
#factor 
wide_df <- data.frame(wide_df)
needs_factoring <- c("treatment", "sex", "drinker")
wide_df[needs_factoring] <- lapply(wide_df[needs_factoring] , factor)
levels(wide_df$treatment) <- c("control", "intervention")

# 20%
set.seed(18)
mon_20 <- ampute(wide_df, prop = 0.12, 
                 pattern = data.frame(id = c(1, 1, 1), 
                                      treatment = c(1, 1, 1), 
                                      age = c(1, 0, 1), 
                                      sex = c(1, 1, 0), 
                                      drinker = c(0, 1, 1),
                                      Y0 = c(1, 1, 1),
                                      Y1 = c(1, 1, 0),
                                      Y2 = c(1, 0, 0),
                                      Y3 = c(0, 0, 0)),
                 freq=c(1/3, 1/3, 1/3),
                 mech = "MAR", 
                 bycases = TRUE)

mon_20 <- mon_20$amp
par(mfrow=c(1,1))
md.pattern(mon_20, rotate.names = TRUE)

# wide to long and cleaning
require(dplyr)
mon_20 <- mon_20 %>%
  pivot_longer(!(c(id, treatment, drinker, age, sex)), 
               names_to="week", values_to="Y")

mon_20$time = ave(seq_len(nrow(mon_20)), mon_20$id, FUN= seq_along)
mon_20$time = as.numeric(mon_20$time)
needs_factoring <- c("treatment", "sex", "drinker")
mon_20[needs_factoring] <- lapply(mon_20[needs_factoring] , factor)

mon_20 <- mon_20 %>%
  dplyr::select(-week)
mon_20 <- data.frame(mon_20)
levels(mon_20$treatment) <- c("control", "intervention")


# 40%
set.seed(18)
mon_40 <- ampute(wide_df, prop = 0.45, 
                 pattern = data.frame(id = c(1, 1, 1), 
                                      treatment = c(1, 1, 1), 
                                      age = c(1, 0, 1), 
                                      sex = c(1, 1, 0), 
                                      drinker = c(0, 1, 1),
                                      Y0 = c(1, 1, 1),
                                      Y1 = c(1, 1, 0),
                                      Y2 = c(1, 0, 0),
                                      Y3 = c(0, 0, 0)),
                 freq=c(1/3, 1/3, 1/3),
                 mech = "MAR", 
                 bycases = TRUE)

mon_40 <- mon_40$amp
par(mfrow=c(1,1))
md.pattern(mon_40, rotate.names = TRUE)

# wide to long and cleaning
require(dplyr)
mon_40 <- mon_40 %>%
  pivot_longer(!(c(id, treatment, drinker, age, sex)), 
               names_to="week", values_to="Y")

mon_40$time = ave(seq_len(nrow(mon_40)), mon_40$id, FUN= seq_along)
mon_40$time = as.numeric(mon_40$time)
needs_factoring <- c("treatment", "sex", "drinker")
mon_40[needs_factoring] <- lapply(mon_40[needs_factoring] , factor)

mon_40 <- mon_40 %>%
  dplyr::select(-week)
mon_40 <- data.frame(mon_40)
levels(mon_40$treatment) <- c("control", "intervention")


# 60%
set.seed(18)
mon_60 <- ampute(wide_df, prop = 0.64, 
                 pattern = data.frame(id = c(1, 1, 1), 
                                      treatment = c(1, 1, 1), 
                                      age = c(1, 0, 1), 
                                      sex = c(1, 1, 0), 
                                      drinker = c(0, 1, 1),
                                      Y0 = c(1, 1, 1),
                                      Y1 = c(1, 1, 0),
                                      Y2 = c(1, 0, 0),
                                      Y3 = c(0, 0, 0)),
                 freq=c(1/3, 1/3, 1/3),
                 mech = "MAR", 
                 bycases = TRUE)

mon_60 <- mon_60$amp
par(mfrow=c(1,1))
md.pattern(mon_60, rotate.names = TRUE)

# wide to long and cleaning
require(dplyr)
mon_60 <- mon_60 %>%
  pivot_longer(!(c(id, treatment, drinker, age, sex)), 
               names_to="week", values_to="Y")

mon_60$time = ave(seq_len(nrow(mon_60)), mon_60$id, FUN= seq_along)
mon_60$time = as.numeric(mon_60$time)
needs_factoring <- c("treatment", "sex", "drinker")
mon_60[needs_factoring] <- lapply(mon_60[needs_factoring] , factor)

mon_60 <- mon_60 %>%
  dplyr::select(-week)
mon_60 <- data.frame(mon_60)
levels(mon_60$treatment) <- c("control", "intervention")


# 80%
set.seed(18)
mon_80 <- ampute(wide_df, prop = 0.80, 
                 pattern = data.frame(id = c(1, 1, 1), 
                                      treatment = c(1, 1, 1), 
                                      age = c(1, 0, 1), 
                                      sex = c(1, 1, 0), 
                                      drinker = c(0, 1, 1),
                                      Y0 = c(1, 1, 1),
                                      Y1 = c(1, 1, 0),
                                      Y2 = c(1, 0, 0),
                                      Y3 = c(0, 0, 0)),
                 freq=c(1/3, 1/3, 1/3),
                 mech = "MAR", 
                 bycases = TRUE)

mon_80 <- mon_80$amp
par(mfrow=c(1,1))
md.pattern(mon_80, rotate.names = TRUE)

# wide to long and cleaning
require(dplyr)
mon_80 <- mon_80 %>%
  pivot_longer(!(c(id, treatment, drinker, age, sex)), 
               names_to="week", values_to="Y")

mon_80$time = ave(seq_len(nrow(mon_80)), mon_80$id, FUN= seq_along)
mon_80$time = as.numeric(mon_80$time)
needs_factoring <- c("treatment", "sex", "drinker")
mon_80[needs_factoring] <- lapply(mon_80[needs_factoring] , factor)

mon_80 <- mon_80 %>%
  dplyr::select(-week)
mon_80 <- data.frame(mon_80)
levels(mon_80$treatment) <- c("control", "intervention")



# list-wise deletion
#list for rmse's
mon_rmse_del <- list()
#20%
model_20 <- brm(Y ~ time*treatment + sex + drinker + age + (1 + time | id), 
              data = mon_20, 
               prior = c(# for intercept 
                 prior(normal(0, 50), class = "Intercept"), 
                 # for slope
                 prior(normal(0, 10), class = "b"),
                 # for interaction
                 prior(normal(0, 5), class = "b", 
                       coef = "time:treatmentintervention"),
                 # for tau_beta0 and tau_beta1
                 prior(gamma(2, 0.2), class = "sd", group = "id"),
                 # for sigma
                 prior(lkj(1), class = "cor")),
               seed = 18)
  
mon_rmse_del[1]  <- rmse(ground_truth, 
                        summary(model_20)$fixed$Estimate)

#40%
model_40 <- brm(Y ~ time*treatment + sex + drinker + age + (1 + time | id), 
                data = mon_40, 
                prior = c(# for intercept 
                  prior(normal(0, 50), class = "Intercept"), 
                  # for slope
                  prior(normal(0, 10), class = "b"),
                  # for interaction
                  prior(normal(0, 5), class = "b", 
                        coef = "time:treatmentintervention"),
                  # for tau_beta0 and tau_beta1
                  prior(gamma(2, 0.2), class = "sd", group = "id"),
                  # for sigma
                  prior(lkj(1), class = "cor")),
                seed = 18)

mon_rmse_del[2]  <- rmse(ground_truth, 
                         summary(model_40)$fixed$Estimate)

# 60%
model_60 <- brm(Y ~ time*treatment + sex + drinker + age + (1 + time | id), 
                data = mon_60, 
                prior = c(# for intercept 
                  prior(normal(0, 50), class = "Intercept"), 
                  # for slope
                  prior(normal(0, 10), class = "b"),
                  # for interaction
                  prior(normal(0, 5), class = "b", 
                        coef = "time:treatmentintervention"),
                  # for tau_beta0 and tau_beta1
                  prior(gamma(2, 0.2), class = "sd", group = "id"),
                  # for sigma
                  prior(lkj(1), class = "cor")),
                seed = 18)

mon_rmse_del[3]  <- rmse(ground_truth, 
                         summary(model_60)$fixed$Estimate)



# 80%
model_80 <- brm(Y ~ time*treatment + sex + drinker + age + (1 + time | id), 
                data = mon_80, 
                prior = c(# for intercept 
                  prior(normal(0, 50), class = "Intercept"), 
                  # for slope
                  prior(normal(0, 10), class = "b"),
                  # for interaction
                  prior(normal(0, 5), class = "b", 
                        coef = "time:treatmentintervention"),
                  # for tau_beta0 and tau_beta1
                  prior(gamma(2, 0.2), class = "sd", group = "id"),
                  # for sigma
                  prior(lkj(1), class = "cor")),
                seed = 18)

mon_rmse_del[4]  <- rmse(ground_truth, 
                         summary(model_80)$fixed$Estimate)



# joint AI
# list for rmse's
mon_rmse_jointai <- list()

# 20%
JA_model_20 <- lme_imp(Y ~ time*treatment + sex + drinker + age + (1 + time | id), 
                       data = mon_20, n.iter = 5000)
summary(JA_model_20)
traceplot(JA_model_20)

JA_20 <- c(22.646201, 0.017823,  -4.587823,  1.181247, -0.000663, -0.026534, -2.133610)
mon_rmse_jointai[1] <- rmse(ground_truth, JA_20)


# 40%
JA_model_40 <- lme_imp(Y ~ time*treatment + sex + drinker + age + (1 + time | id), 
                       data = mon_40, n.iter = 5000)
summary(JA_model_40)

JA_40 <- c(21.44184, -0.34708, -3.81575,  1.51500, 0.14289, 0.00652, -1.75974)
mon_rmse_jointai[2] <- rmse(ground_truth, JA_40)


#60%
JA_model_60 <- lme_imp(Y ~ time*treatment + sex + drinker + age + (1 + time | id), 
                       data = mon_60, n.iter = 5000)

summary(JA_model_60)
traceplot(JA_model_20)

JA_60 <- c(22.2692, 0.8394, -4.2306, 1.0625, -0.5017, -0.0584, -2.1007)
mon_rmse_jointai[3] <- rmse(ground_truth, JA_60)


#80%
JA_model_80 <- lme_imp(Y ~ time*treatment + sex + drinker + age + (1 + time | id), 
                       data = mon_80, n.iter = 5000)

summary(JA_model_80)

JA_80 <- c(23.9802, 0.0249, -4.0599, 1.4298, -0.8108, -0.0869, -1.7354)
mon_rmse_jointai[4] <- rmse(ground_truth, JA_80)

traceplot(JA_model_80)
summary(model3)


#jomo
mon_rmse_jomo <- list()

# specify the levels of the variables
lvl_jomolong <- c(id = 1, time = 1, sex = 2, treatment = 2, age = 2, drinker = 2, Y=1)

#20%
t0 <- Sys.time()
jomo_20 <- jomo.lmer(Y ~ time*treatment + sex + drinker + age + (1 + time | id),
                      data = mon_20, level = lvl_jomolong, meth = "common", output = 0,
                      nburn = 5000, nbetween = 500)

library(miceadds)
jomo_20_mids <- datalist2mids(
  split(jomo_20, jomo_20$Imputation)[-1] # split and remove first element
)

jomo_model_20 <- brm_multiple(Y ~ time*treatment + sex + drinker + age + (1 + time | clus), 
                       prior = c(# for intercept 
                         prior(normal(0, 50), class = "Intercept"), 
                         # for slope
                         prior(normal(0, 10), class = "b"),
                         # for interaction
                         prior(normal(0, 5), class = "b", 
                               coef = "time:treatmentintervention"),
                         # for tau_beta0 and tau_beta1
                         prior(gamma(2, 0.2), class = "sd", group = "clus"),
                         # for sigma
                         prior(lkj(1), class = "cor")), 
                       seed = 1955, 
                       data=jomo_20_mids)

mon_rmse_jomo[1] <- rmse(ground_truth, summary(jomo_model_20)$fixed$Estimate)


#40%
t0 <- Sys.time()
jomo_40 <- jomo.lmer(Y ~ time*treatment + sex + drinker + age + (1 + time | id),
                     data = mon_40, level = lvl_jomolong, meth = "common", output = 0,
                     nburn = 5000, nbetween = 500)

library(miceadds)
jomo_40_mids <- datalist2mids(
  split(jomo_40, jomo_40$Imputation)[-1] # split and remove first element
)

jomo_model_40 <- brm_multiple(Y ~ time*treatment + sex + drinker + age + (1 + time | clus), 
                              prior = c(# for intercept 
                                prior(normal(0, 50), class = "Intercept"), 
                                # for slope
                                prior(normal(0, 10), class = "b"),
                                # for interaction
                                prior(normal(0, 5), class = "b", 
                                      coef = "time:treatmentintervention"),
                                # for tau_beta0 and tau_beta1
                                prior(gamma(2, 0.2), class = "sd", group = "clus"),
                                # for sigma
                                prior(lkj(1), class = "cor")), 
                              seed = 1955, 
                              data=jomo_40_mids)

mon_rmse_jomo[2] <- rmse(ground_truth, summary(jomo_model_40)$fixed$Estimate)




#60%
t0 <- Sys.time()
jomo_60 <- jomo.lmer(Y ~ time*treatment + sex + drinker + age + (1 + time | id),
                     data = mon_60, level = lvl_jomolong, meth = "common", output = 0,
                     nburn = 5000, nbetween = 500)
library(miceadds)
jomo_60_mids <- datalist2mids(
  split(jomo_60, jomo_60$Imputation)[-1] # split and remove first element
)

jomo_model_60 <- brm_multiple(Y ~ time*treatment + sex + drinker + age + (1 + time | clus), 
                              prior = c(# for intercept 
                                prior(normal(0, 50), class = "Intercept"), 
                                # for slope
                                prior(normal(0, 10), class = "b"),
                                # for interaction
                                prior(normal(0, 5), class = "b", 
                                      coef = "time:treatmentintervention"),
                                # for tau_beta0 and tau_beta1
                                prior(gamma(2, 0.2), class = "sd", group = "clus"),
                                # for sigma
                                prior(lkj(1), class = "cor")), 
                              seed = 1955, 
                              data=jomo_60_mids)

mon_rmse_jomo[3] <- rmse(ground_truth, summary(jomo_model_60)$fixed$Estimate)


#80%
t0 <- Sys.time()
jomo_80 <- jomo.lmer(Y ~ time*treatment + sex + drinker + age + (1 + time | id),
                     data = mon_80, level = lvl_jomolong, meth = "common", output = 0,
                     nburn = 5000, nbetween = 500)

imp <- jomo.lmer.MCMCchain(Y ~ time*treatment + sex + drinker + age + (1 + time | id),
                    data = mon_80, level = lvl_jomolong, meth = "common", output = 0,
                    nburn = 5000)

library(miceadds)
jomo_80_mids <- datalist2mids(
  split(jomo_80, jomo_80$Imputation)[-1] # split and remove first element
)


jomo_model_80 <- brm_multiple(Y ~ time*treatment + sex + drinker + age + (1 + time | clus), 
                              prior = c(# for intercept 
                                prior(normal(0, 50), class = "Intercept"), 
                                # for slope
                                prior(normal(0, 10), class = "b"),
                                # for interaction
                                prior(normal(0, 5), class = "b", 
                                      coef = "time:treatmentintervention"),
                                # for tau_beta0 and tau_beta1
                                prior(gamma(2, 0.2), class = "sd", group = "clus"),
                                # for sigma
                                prior(lkj(1), class = "cor")), 
                              seed = 1955, 
                              data=jomo_80_mids)

mon_rmse_jomo[4] <- rmse(ground_truth, summary(jomo_model_80)$fixed$Estimate)


# mice
mon_rmse_mice <- list()

#20%
#21only functions (long format)
mon_20$id <- as.numeric(mon_20$id)
mon_20$treatment <- as.numeric(mon_20$treatment)
mon_20$sex <- as.numeric(mon_20$sex)
mon_20$drinker <- as.numeric(mon_20$drinker)


test_mice <- mice(mon_20, maxit = 0)
meth <- test_mice$method
pred <- test_mice$predictorMatrix

meth[c("drinker", "sex", "age")] <- "2lonly.pmm"
meth[c("Y")] <- "2l.norm"
meth

pred[, "id"] <- -2
pred[, "time"] <- 2
pred

mice20 <- mice(mon_20, m=5, predictorMatrix = pred,
                  method=meth, maxit=30, printFlag = FALSE, seed=1874)

plot(mice20)

mice_list <- list(complete(mice20,1), complete(mice20,2), complete(mice20,3), 
                  complete(mice20,4), complete(mice20,5))

new_mice_list <- list()
for(i in 1:5){
  df <- mice_list[[i]]
  df$time = as.numeric(df$time)
  df$treatment <- as.factor(df$treatment)
  df$sex <- as.factor(df$sex)
  df$drinker <- as.factor(df$drinker)
  df$id <- as.integer(df$id)
  levels(df$treatment) <- c("control", "intervention")
  df <- data.frame(df)
  new_mice_list[[i]] <- df
}


mice_model_20 <- brm_multiple(Y ~ time*treatment + sex + drinker + age + (1 + time | id), 
                              prior = c(# for intercept 
                                prior(normal(0, 50), class = "Intercept"), 
                                # for slope
                                prior(normal(0, 10), class = "b"),
                                # for interaction
                                prior(normal(0, 5), class = "b", 
                                      coef = "time:treatmentintervention"),
                                # for tau_beta0 and tau_beta1
                                prior(gamma(2, 0.2), class = "sd", group = "id"),
                                # for sigma
                                prior(lkj(1), class = "cor")), 
                              seed = 1955, 
                              data=new_mice_list)

mon_rmse_mice[1] <- rmse(ground_truth, summary(mice_model_20)$fixed$Estimate)


#40%
#21only functions (long format)
mon_40$id <- as.numeric(mon_40$id)
mon_40$treatment <- as.numeric(mon_40$treatment)
mon_40$sex <- as.numeric(mon_40$sex)
mon_40$drinker <- as.numeric(mon_40$drinker)


mice40 <- mice(mon_40, m=5, predictorMatrix = pred,
               method=meth, maxit=30, printFlag = FALSE, seed=1874)

plot(mice40)

mice_list40 <- list(complete(mice40,1), complete(mice40,2), complete(mice40,3), 
                  complete(mice40,4), complete(mice40,5))

new_mice_list40 <- list()
for(i in 1:5){
  df <- mice_list40[[i]]
  df$time = as.numeric(df$time)
  df$treatment <- as.factor(df$treatment)
  df$sex <- as.factor(df$sex)
  df$drinker <- as.factor(df$drinker)
  df$id <- as.integer(df$id)
  levels(df$treatment) <- c("control", "intervention")
  df <- data.frame(df)
  new_mice_list40[[i]] <- df
}


mice_model_40 <- brm_multiple(Y ~ time*treatment + sex + drinker + age + (1 + time | id), 
                              prior = c(# for intercept 
                                prior(normal(0, 50), class = "Intercept"), 
                                # for slope
                                prior(normal(0, 10), class = "b"),
                                # for interaction
                                prior(normal(0, 5), class = "b", 
                                      coef = "time:treatmentintervention"),
                                # for tau_beta0 and tau_beta1
                                prior(gamma(2, 0.2), class = "sd", group = "id"),
                                # for sigma
                                prior(lkj(1), class = "cor")), 
                              seed = 1955, 
                              data=new_mice_list40)

mon_rmse_mice[2] <- rmse(ground_truth, summary(mice_model_40)$fixed$Estimate)


#60%
#21only functions (long format)
mon_60$id <- as.numeric(mon_60$id)
mon_60$treatment <- as.numeric(mon_60$treatment)
mon_60$sex <- as.numeric(mon_60$sex)
mon_60$drinker <- as.numeric(mon_60$drinker)


mice60 <- mice(mon_60, m=5, predictorMatrix = pred,
               method=meth, maxit=30, printFlag = FALSE, seed=1874)

plot(mice60)

mice_list60 <- list(complete(mice60,1), complete(mice60,2), complete(mice60,3), 
                    complete(mice60,4), complete(mice60,5))

new_mice_list60 <- list()
for(i in 1:5){
  df <- mice_list60[[i]]
  df$time = as.numeric(df$time)
  df$treatment <- as.factor(df$treatment)
  df$sex <- as.factor(df$sex)
  df$drinker <- as.factor(df$drinker)
  df$id <- as.integer(df$id)
  levels(df$treatment) <- c("control", "intervention")
  df <- data.frame(df)
  new_mice_list60[[i]] <- df
}


mice_model_60 <- brm_multiple(Y ~ time*treatment + sex + drinker + age + (1 + time | id), 
                              prior = c(# for intercept 
                                prior(normal(0, 50), class = "Intercept"), 
                                # for slope
                                prior(normal(0, 10), class = "b"),
                                # for interaction
                                prior(normal(0, 5), class = "b", 
                                      coef = "time:treatmentintervention"),
                                # for tau_beta0 and tau_beta1
                                prior(gamma(2, 0.2), class = "sd", group = "id"),
                                # for sigma
                                prior(lkj(1), class = "cor")), 
                              seed = 1955, 
                              data=new_mice_list60)

mon_rmse_mice[3] <- rmse(ground_truth, summary(mice_model_60)$fixed$Estimate)


#80%
#21only functions (long format)
mon_80$id <- as.numeric(mon_80$id)
mon_80$treatment <- as.numeric(mon_80$treatment)
mon_80$sex <- as.numeric(mon_80$sex)
mon_80$drinker <- as.numeric(mon_80$drinker)


mice80 <- mice(mon_80, m=5, predictorMatrix = pred,
               method=meth, maxit=30, printFlag = FALSE, seed=1874)

plot(mice80)

mice_list80 <- list(complete(mice80,1), complete(mice80,2), complete(mice80,3), 
                    complete(mice80,4), complete(mice80,5))

new_mice_list80 <- list()
for(i in 1:5){
  df <- mice_list80[[i]]
  df$time = as.numeric(df$time)
  df$treatment <- as.factor(df$treatment)
  df$sex <- as.factor(df$sex)
  df$drinker <- as.factor(df$drinker)
  df$id <- as.integer(df$id)
  levels(df$treatment) <- c("control", "intervention")
  df <- data.frame(df)
  new_mice_list80[[i]] <- df
}


mice_model_80 <- brm_multiple(Y ~ time*treatment + sex + drinker + age + (1 + time | id), 
                              prior = c(# for intercept 
                                prior(normal(0, 50), class = "Intercept"), 
                                # for slope
                                prior(normal(0, 10), class = "b"),
                                # for interaction
                                prior(normal(0, 5), class = "b", 
                                      coef = "time:treatmentintervention"),
                                # for tau_beta0 and tau_beta1
                                prior(gamma(2, 0.2), class = "sd", group = "id"),
                                # for sigma
                                prior(lkj(1), class = "cor")), 
                              seed = 1955, 
                              data=new_mice_list80)

summary(mice_model_80)
mon_rmse_mice[4] <- rmse(ground_truth, summary(mice_model_80)$fixed$Estimate)

# MAKING PLOT
percentage_missing <- c(20, 40, 60, 80)
rmse_del1 <- c(mon_rmse_del[[1]], mon_rmse_del[[2]], mon_rmse_del[[3]], mon_rmse_del[[4]])
rmse_jointai1 <- c(mon_rmse_jointai[[1]], mon_rmse_jointai[[2]], mon_rmse_jointai[[3]], mon_rmse_jointai[[4]])
rmse_jomo1 <- c(mon_rmse_jomo[[1]], mon_rmse_jomo[[2]], mon_rmse_jomo[[3]], mon_rmse_jomo[[4]])
rmse_mice1 <- c(mon_rmse_mice[[1]], mon_rmse_mice[[2]], mon_rmse_mice[[3]], mon_rmse_mice[[4]])


cov_df <- data.frame(percentage_missing, rmse_del1, rmse_mice1, rmse_jomo1, rmse_jointai1)

require(dplyr)
cov_df <- cov_df %>%
  rename(FCS = rmse_mice1) %>%
  rename(JointAI = rmse_jointai1) %>%
  rename(jomo = rmse_jomo1) %>%
  rename(deletion = rmse_del1)

cov_df <- cov_df %>%
  pivot_longer(!percentage_missing, names_to="Method", values_to="rmse")


rmse_plot <- ggplot(data = cov_df, aes(x = percentage_missing, y = rmse, group = Method))
rmse_plot + geom_point() + 
  geom_line(aes(col = Method)) +
  scale_x_continuous(breaks = seq(min(cov_df$percentage_missing), max(cov_df$percentage_missing),by=20)) +
  ggtitle("RMSE for each method") + 
  xlab("% missing") + ylab("RMSE")


library(ggplot2)
data.estimates = data.frame(
  meth   = c('ground truth', 'deletion'),
  pt = c(7, 7.5),
  lower = c(5,8),
  upper = c(5,8))

p2 <- ggplot(data.estimates, aes(meth, pt, size=10)) + 
  theme_bw(base_size=10)
p2 + geom_point() +
  geom_errorbar(aes(x = meth, ymin = lower, ymax = upper), width = 0.2) + 
  #scale_y_log10(limits=c(0.1, 50), breaks=c(0.1, 0.5, 1, 5, 10, 25, 50)) +
  xlab("Site") + ylab("RR")

# interaction term plots
par(mfrow = c(2, 2))
#20%
d1=data.frame(method=c("grount truth", "jomo", "mice", "jointAI", "deletion"), 
              est=c(-2.47,  -2.59, -2.34, -2.133610, -2.81), 
              lower=c(-3.93,  -3.97, -3.71, -3.157, -4.36), 
              upper=c( -1.02, -1.12, -0.84, -1.128, -1.21))
plot1 <- ggplot() + 
  geom_errorbar(data=d1, mapping=aes(x=method, ymin=upper, ymax=lower), 
                width=0.2, size=1, color="blue") + 
  geom_point(data=d, mapping=aes(x=method, y=est), size=4, shape=21, fill="white") +
  ggtitle("Interaction estimate - 20% missing")


#40%
d2=data.frame(method=c("grount truth", "jomo", "mice", "jointAI", "deletion"), 
              est=c(-2.47,  -2.35, -1.65, -1.75974, -2.27), 
              lower=c(-3.93,  -4.27 , -3.44, -2.992,  -4.23), 
              upper=c( -1.02, -0.66,  0.06, -0.523, -0.25))
plot2 <- ggplot() + 
  geom_errorbar(data=d2, mapping=aes(x=method, ymin=upper, ymax=lower), 
                width=0.2, size=1, color="blue") + 
  geom_point(data=d, mapping=aes(x=method, y=est), size=4, shape=21, fill="white") +
  ggtitle("Interaction estimate - 40% missing")


#60%
d3=data.frame(method=c("grount truth", "jomo", "mice", "jointAI", "deletion"), 
             est=c(-2.47,  -2.73, -2.01, -2.1007, -1.96), 
             lower=c(-3.93,  -4.28, -3.81, -3.1188, -3.91), 
             upper=c( -1.02, -1.32, -0.27, -1.09, -0.05))
plot3 <- ggplot() + 
  geom_errorbar(data=d3, mapping=aes(x=method, ymin=upper, ymax=lower), 
                width=0.2, size=1, color="blue") + 
  geom_point(data=d, mapping=aes(x=method, y=est), size=4, shape=21, fill="white") +
  ggtitle("Interaction estimate - 60% missing")
 


#80%
d4=data.frame(method=c("grount truth", "jomo", "mice", "jointAI", "deletion"), 
              est=c(-2.47,  -2.35, -1.30, -1.7354, -2.50), 
              lower=c(-3.93,  -3.90 , -3.03, -3.094, -7.25), 
              upper=c( -1.02, -0.75, 0.34, -0.367, 2.93))
plot4 <- ggplot() + 
  geom_errorbar(data=d4, mapping=aes(x=method, ymin=upper, ymax=lower), 
                width=0.2, size=1, color="blue") + 
  geom_point(data=d, mapping=aes(x=method, y=est), size=4, shape=21, fill="white") +
  ggtitle("Interaction estimate - 80% missing")


library(patchwork)
plot1+plot2
plot3+plot4




