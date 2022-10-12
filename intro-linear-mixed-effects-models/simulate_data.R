rm(list = ls())

require(dplyr)
require(tidyr)
require(magrittr)
require(ggplot2)
require(lme4)

set.seed(1234)
input_param <- list(m = 10, sigma_a = 3, sigma_e = 1, n=10, b0=10)
naive_var_of_mean <- (input_param$sigma_a^2 + input_param$sigma_e^2)/(input_param$m*input_param$n)
true_var_of_mean <- (input_param$n*input_param$sigma_a^2 + input_param$sigma_e^2)/(input_param$m*input_param$n)

Nsim <- 100

get_variance_estimates <- function(s, input_param) {
  
  m = input_param$m
  n = rpois(m, input_param$n)
  sigma_a = input_param$sigma_a
  sigma_e = input_param$sigma_e
  b0 = input_param$b0
  
  subject_effect <- rnorm(m, mean = 0, sd = sigma_a)
  Y <- vector(mode = "numeric")
  subject_ID <- NULL
  for(i in 1:m) {
    Y %<>% append(., b0 + subject_effect[i] + rnorm(n[i], mean = 0, sd = sigma_e))
    subject_ID %<>% append(., rep(paste0("subjectID_",i), n[i]))
  }
  
  data <- data.frame(subject_ID, Y)
  
  ggplot(data, aes(x=subject_ID, y=Y)) + geom_boxplot()
  
  mean_data <- data %>%
    group_by(subject_ID) %>%
    summarise(mean_i = mean(Y))
  
  sample_n <- data.frame(subject_ID=paste0("subjectID_",1:m), n=n)
  
  mean_data %<>% merge(., sample_n)
  lm_fit <- lm(mean_i ~ 1, data = mean_data)
  slm_fit <- summary(lm_fit)
  mean(slm_fit$residuals^2)
  summary(slm_fit$residuals^2)
  wt_lm_fit <- lm(mean_i ~ 1, data = mean_data, weights = n)
  swlm_fit <- summary(wt_lm_fit)
  mean(swlm_fit$residuals^2)
  summary(swlm_fit$residuals^2)
  simple_obs_var_of_mean <- mean_data %>%
    .$mean_i %>%
    var() %>%
    divide_by(m)
  
  lmer_fit <- lmer(Y ~ (1|subject_ID), data = data)
  slmer_fit <- summary(lmer_fit)
  
  lmer_obs_var_of_mean <- slmer_fit$coefficients[2]^2
  return(c(simple_obs_var_of_mean, lmer_obs_var_of_mean))
}

sim_res <- lapply(1:Nsim, get_variance_estimates, input_param) %>%
  do.call("rbind", .) %>%
  as.matrix() %>%
  as.data.frame() 

colnames(sim_res) <- c("simple_obs_var_of_mean", "lmer_obs_var_of_mean")
col_means <- sim_res %>%
  colMeans()


sim_res %>%
  ggplot(., aes(x=simple_obs_var_of_mean, y=lmer_obs_var_of_mean)) +
  geom_point() +
  geom_abline(slope = 1, intercept = 0) +
  geom_hline(yintercept = true_var_of_mean, lty=2, col="blue") +
  geom_hline(yintercept = naive_var_of_mean, lty=3, col="red") +
  geom_hline(yintercept = col_means[2], lty=2, col="darkblue") 

long_sim_res <- sim_res %>%
  pivot_longer(everything(),names_to = "estimate_type", values_to = "variance")

long_sim_res %>%
  ggplot(., aes(x=estimate_type, y=variance)) +
  geom_boxplot() +
  geom_hline(yintercept = true_var_of_mean, lty=2, col="blue") +
  geom_hline(yintercept = naive_var_of_mean, lty=3, col="red") +
  geom_hline(yintercept = col_means[2], lty=2, col="darkblue") 


##DATASETS
setwd("~/Dropbox (Gladstone)/scripts/Bioinformatics-Workshops/intro-linear-mixed-effects-models")
##simple random effect
str(Dyestuff)
p <- ggplot(Dyestuff, aes(x=Batch, y=Yield)) +
  geom_boxplot() +
  theme_classic() +
  theme(text = element_text(size = 16)) 

p %>% ggsave(plot = .,
             filename = "Dyestuff.pdf",
             width = 7,
             height = 7)
fm1 <- lmer(Yield ~ 1|Batch, Dyestuff)
summary(fm1)

lm1 <- lm(Yield ~ 1, Dyestuff)
summary(lm1)

summ_Dyestuff <- Dyestuff %>% 
  group_by(Batch) %>% 
  summarise(mYield = mean(Yield))

lm2 <- lm(mYield ~ 1, summ_Dyestuff)
summary(lm2)

##very small inter-batch variability
str(Dyestuff2)
##multiple independent random effects
str(Penicillin)

p <- ggplot(Penicillin, aes(x=plate, y=diameter)) +
  geom_boxplot() +
  geom_point(aes(color=sample)) +
  theme(text = element_text(size = 16)) 

p %>% ggsave(plot = .,
             filename = "Penicillin.pdf",
             width = 10,
             height = 7)


fm1 <- lmer(diameter ~ (1|plate) + (1|sample), Penicillin)
summary(fm1)

lm1 <- lm(diameter ~ 1, Penicillin)
summary(lm1)

##hierarchy
str(Pastes)
##random slope and intercept
str(sleepstudy)

require(lattice)
require(lmerTest)
require(ggplot2)
p <- ggplot(sleepstudy, aes(x=Days, y=Reaction)) +
  geom_point() +
  geom_smooth(method = "lm", se=FALSE) +
  ylab("reaction time in ms") +
  theme_classic() +
  theme(text = element_text(size = 16)) 

p %>% ggsave(plot = .,
             filename = "sleepstudy1.pdf",
             width = 7,
             height = 7)

p <- ggplot(sleepstudy, aes(x=Days, y=Reaction, color=Subject)) +
  geom_point() +
  geom_smooth(method = "lm", se=FALSE) +
  ylab("reaction time in ms") +
  theme_classic() +
  theme(text = element_text(size = 16)) 

p %>% ggsave(plot = .,
             filename = "sleepstudy2.pdf",
             width = 10,
             height = 7)

p <- ggplot(sleepstudy, aes(x=Subject, y=Reaction)) +
  geom_boxplot() +
  ylab("reaction time in ms") +
  theme_classic() +
  theme(text = element_text(size = 16)) 

p <- ggplot(sleepstudy, aes(x=Subject, y=Reaction)) +
  geom_boxplot() +
  geom_point(aes(color=Days)) +
  ylab("reaction time in ms") +
  theme_classic() +
  theme(text = element_text(size = 16)) 
p %>% ggsave(plot = .,
             filename = "sleepstudy3.pdf",
             width = 10,
             height = 7)

xyplot(Reaction ~ Days | Subject, sleepstudy, type = c("g","p","r"),
       index = function(x,y) coef(lm(y ~ x))[1],
       xlab = "Days of sleep deprivation",
       ylab = "Average reaction time (ms)", aspect = "xy")

fm1 <- lme4::lmer(Reaction ~ Days + (Days | Subject), sleepstudy)
summary(fm1)

fm2 <- lmerTest::lmer(Reaction ~ Days + (Days || Subject), sleepstudy)
summary(fm2)

lm1 <- lm(Reaction ~ Days, data = sleepstudy)
summary(lm1)
