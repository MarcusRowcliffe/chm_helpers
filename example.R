library(GLMMadaptive)
library(ggplot2)
source("make_chm_data.R")

# Load Ianarilli example data
# observations
obs <- read.csv("data/species_records.csv") %>%
  filter(Species == "BlackBear") %>% 
  droplevels %>% 
  select(-X) %>% 
  mutate(timestamp = ymd_hm(DateTimeOriginal),
         deploymentID = paste(Session, Station, sep="_"))

# deployments
dep <- read.csv("data/CameraTrapProject_CT_data_for_analysis_MASTER.csv", as.is = TRUE) %>% 
  select(-X) %>%
  mutate(Date_setup = mdy(Date_setup),
         Date_retr = mdy(Date_retr),
         Problem1_from = mdy(Problem1_from),
         Problem1_to = mdy(Problem1_to),
         deploymentStart = as.POSIXct(if_else(is.na(Problem1_from), 
                                              Date_setup, 
                                              Problem1_from)) + 
           sample(1:(24*60^2), n(), replace=T),
         deploymentEnd = as.POSIXct(if_else(is.na(Problem1_to), 
                                            Date_retr, 
                                            Problem1_to)) +
           sample(1:(24*60^2), n(), replace=T),
         season = as.factor(substr(Session, 1, 1)),
         deploymentID = paste(Session, Site, sep="_")) %>%
  rename(locationName = Site)


# Generate CHM data, porting season and ghm covariates from deployments
dat <- make_chm_data(dep, obs, covs=c("season", "ghm"))
View(dat)

#Fit basic models
# Unimodal 
unimodal <- mixed_model(fixed = cbind(success, failure) ~ cos(timeRadian) + sin(timeRadian), 
                        random = ~ 1 | locationName,
                        family = binomial(),
                        data = dat)
summary(unimodal)

# Bimodal 
bimodal <- mixed_model(fixed = cbind(success, failure) ~ 
                         cos(timeRadian) + sin(timeRadian) +
                         cos(2*timeRadian) + sin(2*timeRadian), 
                       random = ~ 1 | locationName,
                       family = binomial(),
                       data = dat)
summary(bimodal)

# Cathemeral
null_mod <- mixed_model(fixed = cbind(success, failure) ~ 1, 
                        random = ~ 1 | locationName,
                        family = binomial(),
                        data = dat)
summary(null_mod)

# AIC comparison
AIC(null_mod, unimodal, bimodal)

# Generate predictions 
newdat <- data.frame(timeRadian = seq(0, 24, len=100) * pi / 12)

predict_unimodal <- effectPlotData(unimodal, newdat, marginal = FALSE) %>% 
  mutate(Model = "Unimodal")
predict_bimodal <- effectPlotData(bimodal, newdat, marginal = FALSE) %>% 
  mutate(Model = "Bimodal")
predict_cathemeral <- effectPlotData(null_mod, newdat, marginal = FALSE) %>% 
  mutate(Model = "Cathemeral")

# join and plot results
rbind(predict_unimodal, predict_bimodal, predict_cathemeral) %>% 
  ggplot(., aes(x = timeRadian, y = plogis(pred), group = Model, fill = Model)) +
  geom_line(aes(colour = Model)) +
  geom_ribbon(aes(ymin = plogis(low), ymax = plogis(upp), colour = NULL), alpha = 0.3) +
  labs(x = "Time of Day (Hour)", y = "Predicted Activity Pattern \n (probability)")+
  theme_minimal()+
  theme(legend.position = "top",
        legend.title = element_blank(),
        legend.text = element_text(size=10,face="bold"),
        legend.margin=margin(0,0,0,0),
        legend.box.margin=margin(-5,-10,-10,-10),
        plot.title = element_blank(),
        axis.line = element_line(colour = 'black', linetype = 'solid'),
        axis.ticks = element_line(colour = 'black', linetype = 'solid'),
        axis.title = element_text(size=9,face="bold"),
        panel.grid.minor.y = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.major.x = element_line(colour = 'lightgrey', linetype = 'dashed', linewidth=0.5),
        panel.grid.minor.x = element_blank(),
        strip.text = element_text(size = 9, colour = "black", face = "bold", hjust = 0)
  ) +
  scale_x_continuous(breaks=seq(0, 2*pi, len=5), labels=seq(0, 24, len=5))


# Covariate models 
bi_ssn <- mixed_model(fixed = cbind(success, failure) ~ 
                        cos(timeRadian) * season + 
                        sin(timeRadian) * season +
                        cos(2*timeRadian) * season + 
                        sin(2*timeRadian) * season, 
                      random = ~ 1 | locationName,
                      family = binomial(),
                      data = dat
)

bi_ghm <- mixed_model(fixed = cbind(success, failure) ~ 
                        cos(timeRadian) * ghm + 
                        sin(timeRadian) * ghm +
                        cos(2*timeRadian) * ghm + 
                        sin(2*timeRadian) * ghm, 
                      random = ~ 1 | locationName,
                      family = binomial(),
                      data = dat)

bi_ssn_ghm <- mixed_model(fixed = cbind(success, failure) ~ 
                            cos(timeRadian) * season + cos(timeRadian) * ghm + 
                            sin(timeRadian) * season + sin(timeRadian) * ghm +
                            cos(2*timeRadian) * season + cos(2*timeRadian) * ghm + 
                            sin(2*timeRadian) * season + sin(2*timeRadian) * ghm, 
                          random = ~ 1 | locationName,
                          family = binomial(),
                          data = dat)

summary(bi_ssn_ghm)

bi_ssnXghm <- mixed_model(fixed = cbind(success, failure) ~ 
                            cos(timeRadian) * season * ghm + 
                            sin(timeRadian) * season * ghm +
                            cos(2*timeRadian) * season * ghm + 
                            sin(2*timeRadian) * season * ghm, 
                          random = ~ 1 | locationName,
                          family = binomial(),
                          data = dat)
summary(bi_ssnXghm)
AIC(bimodal, bi_ssn, bi_ghm, bi_ssn_ghm, bi_ssnXghm) %>%
  mutate(dAIC = AIC-min(AIC)) %>%
  arrange(AIC)

# build estimate of activity
newdat <- expand.grid(timeRadian = seq(0, 24, len=100) * pi / 12,
                      season = levels(dat$season))
newdat <- expand.grid(timeRadian = seq(0, 24, len=100) * pi / 12,
                      season = levels(dat$season),
                      ghm = quantile(dat$ghm, c(0.025, 0.975)))
pred_ssn <- effectPlotData(bi_ssn, newdat, marginal = FALSE) %>%
  mutate(ghm = ifelse(ghm==min(ghm), "Low", "High"))
pred_ssn_ghm <- effectPlotData(bi_ssn_ghm, newdat, marginal = FALSE) %>%
  mutate(ghm = ifelse(ghm==min(ghm), "Low", "High"))
pred_ssnXghm <- effectPlotData(bi_ssnXghm, newdat, marginal = FALSE) %>%
  mutate(ghm = ifelse(ghm==min(ghm), "Low", "High"))

# plot
pred_ssnXghm %>%
  ggplot(aes(timeRadian, plogis(pred))) +
  geom_ribbon(aes(ymin = plogis(low), ymax = plogis(upp), 
                  color = interaction(ghm, season), fill = interaction(ghm, season)), 
              alpha = 0.3, linewidth = 0.25) +
  geom_line(aes(color = interaction(ghm, season)), linewidth = 1) +
  coord_cartesian(ylim = c(0, 0.005)) +
  labs(x = "Time of Day (Hour)", y = "Predicted Activity Pattern \n (probability)") +
  theme_minimal() +
  theme(legend.position = "bottom",
        legend.title = element_blank(),
        legend.text = element_text(size=10,face="bold"),
        legend.margin=margin(0,0,0,0),
        legend.box.margin=margin(-5,-10,-10,-10),
        plot.title = element_text(size=10,face="bold"),
        axis.line = element_line(colour = 'black', linetype = 'solid'),
        axis.ticks = element_line(colour = 'black', linetype = 'solid'),
        axis.title = element_text(size=9,face="bold"),
        panel.grid.minor.y = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.major.x = element_line(colour = 'lightgrey', 
                                          linetype = 'dashed', linewidth=0.5),
        panel.grid.minor.x = element_blank(),
        strip.text = element_text(size = 9, colour = "black", 
                                  face = "bold", hjust = 0)
  ) +
  scale_x_continuous(breaks=seq(0, 2*pi,len=5), labels=seq(0, 24, len=5)) +
  facet_grid(season~ghm)
