require(dplyr)
require(lubridate)
require(hms)

#' make_chm_data
#' Function generates data for circular hierarchical modelling from standard 
#' deployment and observation tables. Works by creating daily time bins,
#' defining occasions uniquely identified by deploymentID, date and time bin,
#' and noting whether any observations fall into each occasion.
#' 
#' INPUT
#' deployments: a data frame of deployment data with required fields:
#'  deploymentID, locationName, deploymentStart, deploymentEnd
#' observations: a data frame of observation data with required fields:
#'  deploymentID, timestamp
#' nBins: number of time bins (default 24 = hourly)
#' covs: vector of character covariates in deployments data
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
#'      capture: a binary indicating whether any observations occurred in each
#'        occasion
#'      occasionStart, occasionEnd: the start and end date-times of each occasion
#'      
make_chm_data <- function(deployments, observations,
                          nBins = 24,
                          covs = NULL,
                          collapse = TRUE){
  
  required_dep_fields <- c("deploymentID", "locationName", "deploymentStart", "deploymentEnd")
  required_obs_fields <- c("deploymentID", "timestamp")
  stopifnot(all(covs %in% names(deployments)))
  stopifnot(all(required_dep_fields %in% names(deployments)))
  stopifnot(all(required_obs_fields %in% names(observations)))
  
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

