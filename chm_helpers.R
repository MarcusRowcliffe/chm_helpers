require(dplyr)
require(lubridate)
require(hms)
require(glmmTMB)

#' make_chm_data
#' Function generates data for circular hierarchical modelling from standard 
#' deployment and observation tables. Works by creating daily time bins,
#' defining occasions uniquely identified by deploymentID, date and time bin,
#' and noting whether any observations fall into each occasion.
#' 
#' INPUT
#' deployments: a data frame of deployment data with required fields:
#'  deploymentID: unique deployment identifier, alphanumeric
#'  locationName: deployment name, alphanumeric
#'  deploymentStart, deploymentEnd: POSIXct deployment start and end times
#' observations: a data frame of observation data with required fields:
#'  deploymentID: deployment identifier linking to deployments$deploymentID
#'  timestamp: POSIXct time of observation
#' nBins: integer number of time bins per day (default 24 = hourly)
#' covs: vector of character covariate fields in deployments data
#' collapse: whether to collapse data into location x time x covariate 
#'  categories (see output)
#'
#' OUTPUT
#'  A dataframe of occasions with columns:
#'    locationID: location identifier, taken from deployments
#'    Any covariates defined by the covs argument, taken from deployments
#'    timeHour, timeRadian: respectively hour and radian time of day of 
#'      occasion midpoints
#'    Plus, depending on collapse setting, either:
#'    collapse = TRUE
#'      success, failure: respectively numbers of occasions with and without captures
#'    collapse = FALSE
#'      capture: a binary indicating whether any observations occurred in each occasion
#'      occasionStart, occasionEnd: the start and end date-times of each occasion
#'      
make_chm_data <- function(deployments, observations,
                          nBins = 24,
                          covs = NULL,
                          collapse = TRUE){
  
  required_dep_fields <- c("deploymentID", "locationName", 
                           "deploymentStart", "deploymentEnd",
                           covs)
  required_obs_fields <- c("deploymentID", "timestamp")
  if(!all(required_dep_fields %in% names(deployments)))
    stop(paste("deployments must contain columns:",
               paste(required_dep_fields, collapse=", ")))
  if(!all(required_obs_fields %in% names(observations)))
    stop(paste("observations must contain columns:",
               paste(required_obs_fields, collapse=", ")))
  if(!(inherits(deployments$deploymentStart, "POSIXt") & 
       inherits(deployments$deploymentEnd, "POSIXt") &
       inherits(observations$timestamp, "POSIXt")))
    stop("deploymentStart, deploymentEnd and timestamp must be POSIX values")
  
  timeSeq <- hms::as_hms(seq(0, 60^2*24, len=1+nBins))
  interval <- timeSeq[2]
  
  obs <- observations %>%
    mutate(time = hms::as_hms(timestamp), # time of observation
           bin = findInterval(time, timeSeq), # occasion bin index
           occasionID = paste(deploymentID, 
                              lubridate::as_date(timestamp),
                              timeSeq[bin]))
  
  # time sequence with terminal overshoot
  tseq <- function(from, to, by){
    by <- as.numeric(by)
    from + (0:ceiling(difftime(to, from, units="secs") / by)) * by
  }
  
  occasions <- deployments %>%
    mutate(timeStart = hms::as_hms(deploymentStart), # deployment start time
           binStart = findInterval(timeStart, timeSeq), # sequence starting bin index
           seqStart = lubridate::ymd_hms(paste(date(deploymentStart),  # starting date-time of deployment sequence
                                               timeSeq[binStart]))) %>%
    dplyr::rowwise() %>%
    dplyr::mutate(occasionStart = list(tseq(seqStart, deploymentEnd, by = interval))) %>%
    tidyr::unnest(occasionStart) %>%
    dplyr::mutate(occasionEnd = occasionStart + interval,
                  timeHour = as.numeric(hms::as_hms(occasionStart + interval/2)) / 60^2,
                  timeRadian = timeHour * pi / 12,
                  occasionID = paste(deploymentID, 
                                     occasionStart),
                  capture = occasionID %in% obs$occasionID) %>%
    dplyr::select(locationName, all_of(covs), occasionStart, occasionEnd,
                  timeHour, timeRadian, capture)
  
  if(collapse){
    occasions <- occasions %>%
      dplyr::summarize(success = sum(capture),
                       failure = n() - success,
                       .by = all_of(c("locationName",
                                      covs,
                                      "timeHour", 
                                      "timeRadian")))
  }
  occasions
}

#' chmTMB
#' Fit a circular hierarchical model, specifically a trigonometric binomial 
#' GLMM, fitted using glmmTMB.
#' 
#' INPUT
#'  formula: a formula for the fixed effects, typically:
#'    cbind(success, failure) ~ time + ...
#'    time is a reserved term that is expanded into trigonometric terms internally
#'  type: Activity pattern type (number of activity peaks)
#'  data: a dataframe containing radianTime and the variables named in formula
#'  ...: additional arguments passed to glmmTMB
#'  
#'  OUTPUT
#'    A model object created by glmmTMB::glmmTMB
chmTMB <- function(formula, type = c("bimodal", "unimodal"), data = NULL, ...){
  type = match.arg(type)
  trigTerms <- switch(type,
                      unmodal = "(cos(timeRadian) + sin(timeRadian))",
                      bimodal = "(cos(timeRadian) + sin(timeRadian) + cos(2*timeRadian) + sin(2*timeRadian))")
  formula <- as.formula(gsub("\\btime\\b", 
                             trigTerms, 
                             paste(deparse(formula), collapse = "")))
  mod <- glmmTMB(formula, family=binomial, data=data, ...)
  class(mod) <- c("chmTMB", class(mod))
  mod
}

#' predict.chm
#' Predict activity probabilities from a CHM
#' 
#' INPUT
#' mod: a model fit
#' predictors: a named list of fixed effect values, expanded to all 
#'    combinations internally if multiple elements provided
#'    
#' OUTPUT
#' A dataframe containing, time, radian time, the predictor values, and:
#'  fit: predicted values on the link scale
#'  se.fit: standard errors of fit
#'  response: predicted probabilities
#'  lcl.response, ucl.response: lower and upper probability confidence limits
predict.chmTMB <- function(object, predictors){
  terms <- attr(terms(object), "term.labels")
  terms <- terms[!grepl("timeRadian", terms) & !grepl(":", terms)]
  if(!all(terms %in% names(predictors)))
    stop("Not all required predictors have been provided")
  predictors <- predictors[names(predictors) %in% terms]
  nd <- expand.grid(c(list(timeRadian = seq(0, 2*pi, len=100)), predictors))
  class(object) <- "glmmTMB"
  predict(object, newdata=nd, se.fit=TRUE, re.form=NA) %>%
    as.data.frame() %>%
    dplyr::bind_cols(nd) %>%
    mutate(time = timeRadian * 12 / pi, 
           response = plogis(fit),
           lcl.response = plogis(fit - 1.96 * se.fit),
           ucl.response = plogis(fit + 1.96 * se.fit))
}
