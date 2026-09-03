install.packages("TMB")
remove.packages("TMB")
install.packages("glmmTMB")
remove.packages("glmmTMB")

library(GLMMadaptive)
library(glmmTMB)
library(ggplot2)
source("https://raw.githubusercontent.com/MarcusRowcliffe/make_chm_data/refs/heads/main/make_chm_data.R")

# Load Ianarilli example data
# observations
obs <- read.csv("data/species_records.csv") %>%
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
         deploymentID = paste(Session, Site, sep="_"),
         cluster = gsub("[^0-9]", "", Site)) %>%
  rename(locationName = Site)


# Generate CHM data, porting season and ghm covariates from deployments
dat <- make_chm_data(dep,
                     subset(obs, Species == "BlackBear"),
                     covs=c("season", "ghm", "locationName", "cluster", "Session"))
View(dat)

# Fit basic models...
# ...using mixed_model
uniform <- mixed_model(fixed = cbind(success, failure) ~ 1,
                       random = ~ 1 | locationName,
                       family = binomial(),
                       data = dat)
unimodal <- mixed_model(fixed = cbind(success, failure) ~ 
                          cos(timeRadian) + sin(timeRadian),
                        random = ~ 1 | locationName,
                        family = binomial(),
                        data = dat)
bimodal <- mixed_model(fixed = cbind(success, failure) ~ 
                         cos(timeRadian) + sin(timeRadian) +
                         cos(2*timeRadian) + sin(2*timeRadian), 
                       random = ~ 1 | locationName,
                       family = binomial(),
                  data = dat)

# AIC comparison
AIC(uniform, unimodal, bimodal)

# Generate predictions 
newdat <- data.frame(timeRadian = seq(0, 24, len=100) * pi / 12)

predict_uniform <- effectPlotData(uniform, newdat, marginal = FALSE) %>% 
  mutate(Model = "Uniform")
predict_unimodal <- effectPlotData(unimodal, newdat, marginal = FALSE) %>% 
  mutate(Model = "Unimodal")
predict_bimodal <- effectPlotData(bimodal, newdat, marginal = FALSE) %>% 
  mutate(Model = "Bimodal")

# join and plot results
rbind(predict_unimodal, predict_bimodal, predict_uniform) %>% 
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


# Fit full fixed effects with various random effects structures
TxSxH_C_SxHxL <- chmTMB(cbind(success, failure) ~ time * season * ghm + (1 | cluster) + (1 + season + ghm | locationName), data=dat)
TxSxH_C_SxL <- chmTMB(cbind(success, failure) ~ time * season * ghm + (1 | cluster) + (1 + season | locationName), data=dat)
TxSxH_C_HxL <- chmTMB(cbind(success, failure) ~ time * season * ghm + (1 | cluster) + (1 + ghm | locationName), data=dat)
TxSxH_C_L <- chmTMB(cbind(success, failure) ~ time * season * ghm + (1 | cluster) + (1 | locationName), data=dat)
TxSxH_SxHxL <- chmTMB(cbind(success, failure) ~ time * season * ghm + (1 + season + ghm | locationName), data=dat)
TxSxH_SxL <- chmTMB(cbind(success, failure) ~ time * season * ghm + (1 + season | locationName), data=dat)
TxSxH_HxL <- chmTMB(cbind(success, failure) ~ time * season * ghm + (1 + ghm | locationName), data=dat)
TxSxH_L <- chmTMB(cbind(success, failure) ~ time * season * ghm + (1 | locationName), data=dat)
# Check model support
BIC(TxSxH_C_SxHxL, TxSxH_C_SxL, TxSxH_C_HxL, TxSxH_C_L, TxSxH_SxHxL, TxSxH_SxL, TxSxH_HxL, TxSxH_L) %>%
  mutate(dBIC = BIC - min(BIC, na.rm=TRUE)) %>%
  arrange(dBIC)

TxSpH_SxL <- chmTMB(cbind(success, failure) ~ time * season + ghm + (1 + season | locationName), data=dat)
TxHpS_SxL <- chmTMB(cbind(success, failure) ~ time * ghm + season + (1 + season | locationName), data=dat)
TpSxH_SxL <- chmTMB(cbind(success, failure) ~ time + season * ghm + (1 + season | locationName), data=dat)
TpSpH_SxL <- chmTMB(cbind(success, failure) ~ time + season + ghm + (1 + season | locationName), data=dat)
TxS_SxL <- chmTMB(cbind(success, failure) ~ time * season + (1 + season | locationName), data=dat)
TxH_L <- chmTMB(cbind(success, failure) ~ time * ghm + (1 | locationName), data=dat)
TpS_SxL <- chmTMB(cbind(success, failure) ~ time + season + (1 + season | locationName), data=dat)
TpH_L <- chmTMB(cbind(success, failure) ~ time + ghm + (1 | locationName), data=dat)
T_L <- chmTMB(cbind(success, failure) ~ time + (1 | locationName), data=dat)
# Check model support
AIC(TxSxH_SxHxL, TxSxH_SxL, TxSpH_SxL, TxHpS_SxL, TpSxH_SxL, TpSpH_SxL, TxS_SxL, TxH_L, TpS_SxL, TpH_L, T_L) %>%
  mutate(dBIC = AIC - min(AIC, na.rm=TRUE)) %>%
  arrange(dBIC)

predict.chm(TxS_SxL, list(season = levels(dat$season))) %>%
  mutate(season = ifelse(season=="F", "Fall", "Spring")) %>%
  ggplot(aes(time, response, col=season)) +
  geom_ribbon(aes(ymin = lcl.response, ymax = ucl.response, fill = season), 
              col = NA, alpha = 0.2) +
  geom_line() +
  scale_x_continuous(breaks=seq(0, 24, len=5)) +
  theme_classic()

predict.chm(TxSpH_SxL, list(season = levels(dat$season),
                        ghm = quantile(dat$ghm, c(0.025, 0.975)))) %>%
  mutate(ghm_level = ifelse(ghm==min(ghm), "Low", "High"),
         season = ifelse(season=="F", "Fall", "Spring")) %>%
  ggplot(aes(time, response, col=season, group=interaction(season, ghm_level))) +
  geom_ribbon(aes(ymin = lcl.response, ymax = ucl.response, fill = season), 
              col = NA, alpha = 0.2) +
  geom_line() +
  scale_x_continuous(breaks=seq(0, 24, len=5)) +
  facet_grid(~ghm_level) +
  theme_classic()

predict.chm(TxSxH_SxHxL, list(season = levels(dat$season),
                              ghm = quantile(dat$ghm, c(0.025, 0.975)))) %>%
  mutate(ghm_level = ifelse(ghm==min(ghm), "Low", "High"),
         season = ifelse(season=="F", "Fall", "Spring")) %>%
  ggplot(aes(time, response, col=season, group=interaction(season, ghm_level))) +
  geom_ribbon(aes(ymin = lcl.response, ymax = ucl.response, fill = season), 
              col = NA, alpha = 0.2) +
  geom_line() +
  scale_x_continuous(breaks=seq(0, 24, len=5)) +
  facet_grid(~ghm_level) +
  theme_classic()

