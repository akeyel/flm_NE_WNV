#' Functional Linear Modeling for Vector-borne Disease
#' 
#' Use functional linear modeling statistical approach to predict mosquito-borne disease.
#' Developed in the context of West Nile virus in Nebraska
#' 
#' @docType package
#' @name flm
NULL
#**# NOTE flm has not been checked to see if other packages exist with this name. The package name can easily be changed at this point

# Code for building package (Run during development only)
#devtools::document()
#devtools::check()
#devtools::load_all()
#spelling::spell_check_packages()
#
# Save the .csv data as .rda for use in the package.
# Note the objects need to be read in and are saved with the same object name
#usethis::use_data(NE_county_pops)
#usethis::use_data(sampledat)
#usethis::use_data(NEdat)
#usethis::use_data(spi)
#usethis::use_data(spei)


# Document data sets
#' NE_county_pops
#'
#' Nebraska County populations, 2000-2018, from the U.S. Census Bureau, using annual estimates
#' @docType data
#'
"NE_county_pops"

#' sampledat
#'
#' simulated data on annual numbers of human cases of neuro-invasive and non-neuro-invasive
#' West Nile Virus in Nebraska counties. It is predictions of a model that was trained on
#' actual numbers of cases as recorded in CDC's Arbonet database. It excludes Arthur County,
#' because no cases have been recorded there to date, and we had to exclude it from our
#' modeling to get it to work.
#' @docType data
#'
"sampledat"

#' NEdat
#'
#' temperature and precipitation data for Nebraska counties each month from
#' January 1998 to February 2019 from the National Centers for Environmental
#' Information, National Climatic Data Center 
#' @docType data
#' @source ftp://ftp.ncdc.noaa.gov/pub/data/cirs/climdiv/
#'
"NEdat"

#' spi
#'
#' Extracted monthly values from Westwide Drought Tracker netcdf files through
#' 2019-12-15. See Abatzoglou, J. T., McEvoy, D. J., & Redmond, K. T. (2017).
#' The West Wide Drought Tracker: Drought monitoring at fine spatial scales.
#' Bulletin of the American Meteorological Society, 98(9), 1815–1820.
#' https://doi.org/10.1175/BAMS-D-16-0193.1
#' @docType data
#'
"spi"

#' spei
#'
#' Extracted monthly values from Westwide Drought Tracker netcdf files through
#' 2019-12-15. See Abatzoglou, J. T., McEvoy, D. J., & Redmond, K. T. (2017).
#' The West Wide Drought Tracker: Drought monitoring at fine spatial scales.
#' Bulletin of the American Meteorological Society, 98(9), 1815–1820.
#' https://doi.org/10.1175/BAMS-D-16-0193.1
#' @docType data
#'
"spei"

# Identify dependencies used in dfmip (see also DESCRIPTION file)
# Set up package imports for NAMESPACE
#**# Importing specific functions using @importFrom is a better practice
#**# Watch for trouble with predict - it may be from one of the packages, and not stats
#' @import dplyr
#' @import stringr
#' @import tidyr
#' @importFrom lubridate parse_date_time
#' @import readxl
#' @import broom
#' @import gamm4
#' @import lme4
#' @import mgcv
#' @importFrom MuMIn AICc
#' @import purrr
#' @import ggplot2
#' @importFrom stats contr.sum contrasts<- gaussian na.omit predict rnbinom sd update contrasts formula
#' @importFrom utils write.csv
NULL

#' Main function
#' 
#' @param pop #**# Add documentation
#' @param cases #**# Add documentation
#' @param NEdat #**# Add documentation
#' @param spi #**# Add documentation
#' @param spei #**# Add documentation
#' @param target.date The last date to include for calculation of lags
#' @param start.year The first year to include in the training data
#' @param in.seed The starting number for the random number generator. This makes the results repeatable.
#' @export
call.flm = function(pop, cases, NEdat, spi, spei, target.date = "2018-02-01",
                    start.year = 2002, in.seed = 4872957){

  # Assemble data lags
  message("Assembling Data")
  start.time = Sys.time()
  adl.out = assemble.data.lags(pop, cases, NEdat, spi, spei, target.date, start.year, in.seed)
  allLagsT = adl.out[[1]]
  allLagsO = adl.out[[2]]
  message(sprintf("Elapsed Time: %.2f", Sys.time() - start.time))
  
  
  # Compare models with and without lags
  message("Comparing models with and without lags")

  tlag = c(12, 18, 24, 30, 36)
  models <- c("cases ~ s(lags_tmean%d, by=tmean%d) + County + year + offset(log(pop100K))",
              "cases ~ s(lags_tmean%d, by=tmean%d) + CI + County + year + offset(log(pop100K))",
              
              "cases ~ s(lags_ppt%d, by=ppt%d) + County + year + offset(log(pop100K))",
              "cases ~ s(lags_ppt%d, by=ppt%d) + CI + County + year + offset(log(pop100K))",
              
              "cases ~ s(lags_spi%d, by=spi%d) + County + year + offset(log(pop100K))",
              "cases ~ s(lags_spi%d, by=spi%d) + CI + County + year + offset(log(pop100K))",
              
              "cases ~ s(lags_spei%d, by=spei%d) + County + year + offset(log(pop100K))",
              "cases ~ s(lags_spei%d, by=spei%d) + CI + County + year + offset(log(pop100K))")
  allmods_list <- map(tlag,
                      ~sprintf(models, .x, .x))
  allmods <- flatten(c(list("cases ~ County + year + offset(log(pop100K))",
                            "cases ~ CI + County + year + offset(log(pop100K))"),
                       allmods_list
                       ))
  
  process.start = Sys.time()

  results <- models_lags(allmods[c(1, 3, 5)], allLagsT, allLagsO, results.path) 

  message(sprintf("Elapsed Time: %.2f; Process time: %.2f", (Sys.time() - start.time), (Sys.time() - process.start)))
  
  #**# Update when these are extracted in a format that can be passed to dfmip
  flm.results = results$predictions
  flm.distributions = NA
  flm.other = results$other
  
  # Return model results as a list
  #**# Need to return data in something that can be extracted by the update.df function
  #**# Need to return distribution if available - was that done? Otherwise, can use the point estimate from first output
  #**# Can return anything else as an additional output. This can be returned as a keyed list from dfmip #**# eventually
  return(list(flm.results, flm.distributions, flm.other))  
}

# Functions added to enable compatibility with dfmip package

#' Run the Functional Linear Modeling approach using dfmip inputs
#'
#' @param human.data #**# NEEDS DOCUMENTATION. human.data will be ignored if flm.inputs[[4]] is not NA
#' @param weather.data #**# NEEDS DOCUMENTATION, weather.data will be ignored if flm.inputs[[5]] is not NA.
#' @param flm.inputs Inputs specific to the functional linear modeling approach: \tabular{ll}{
#' pop \tab Human population for each county by year. Required fields are County, cofips, year, pop #**# is cofips needed?\cr
#' spi \tab Drought information corresponding to the spi drought index. #**# Update with link/chunk to other documentation\cr
#' spei \tab Drought information corresponding to the spei drought index #**# Update with link/chunk to other documentation\cr
#' NEdat \tab Weather data in format cofips, year, month, day, temp, ppt; where temp is mean temperature and ppt is mean precipitation.
#' If present, the weather.data argument will be ignored. If absent, this object will be created from the weather.data object.\cr}
#' @param weekinquestion The date for which the forecast is desired
#' @param results.path The path for outputting model results
#' @param in.seed A random seed to make the results repeatable
#'
#' @return flm.out A list of 3 objects: a results data frame, a data frame with random distribution draws, and a list of other outputs.
#'
#' @export run.flm
#' 
run.flm = function(human.data, weather.data, flm.inputs, weekinquestion, results.path, in.seed = 4872957){
  
  # Currently this is handled in dfmip. Can be moved back here if probabilistic forecasts are generated
  # #' @param n.draws The number of draws to include in the returned forecast.distributions object
  #n.draws, 
  
  # For now, have the other inputs come in in the flm.inputs object.
  # Probably want the population inputs to be standardized
  # May be options for the weather data as well
  pop = flm.inputs[[1]]
  spi = flm.inputs[[2]]
  spei = flm.inputs[[3]]
  cases = flm.inputs[[4]]
  NEdat = flm.inputs[[5]]
  start.year = flm.inputs[[6]]
  
  # Convert human.data to cases object
  if (is.na(cases)){
    cases = flm.convert.human.data(human.data)
  }
  
  # Convert weather.data to NEdat object
  if (is.na(NEdat)){
    cofips.lookup = generate.lookup(pop)
    NEdat = flm.convert.weather.data(weather.data, cofips.lookup)
  }
  
  message(weekinquestion)
  
  flm.out = call.flm(pop, cases, NEdat, spi, spei, target.date = weekinquestion,
                     start.year = start.year) #, results.path = results.path
  
  # Adjust forecast.distributions object to have the correct number of draws
  #**# Need to discuss options for probabilistic forecasts
  #**# Currently this step is done in dfmip
  
  # Add zeros for counties that never had a human case, using analysis.districts
  #**# Currently this step is done in dfmip
  
  # For now, assume the outputs are in the correct format
  return(flm.out)
}

#' Convert dfmip standard human case data to flm format
#' 
#' Input data is in format: district date year. Output needs to be: County year cases.
#' @param human.data Human case data by date of symptom onset. Needs a row for each case, and needs three columns: district, date, and year
#' 
#' @return cases Human cases summarized by county and year. Needs one row per county-year. Output columns are County, year, cases.
#' 
#' @noRd
#' 
flm.convert.human.data = function(human.data){
  
  human.data$count = 1 # Add a column for summing in the next step
  cases = stats::aggregate(human.data$count, by = list(human.data$year, human.data$district), sum, na.rm = TRUE)
  cases$County = cases$Group.2
  cases$year = cases$Group.1
  cases$cases = cases$x
  cases$county_year = sprintf("%s_%s", cases$County, cases$year)
  
  # Remove non-required fields
  cases$Group.1 = NULL
  cases$Group.2 = NULL
  cases$x = NULL
  
  # Add zeros for missing years. #**# Assumes no missing internal data!
  districts = unique(cases$County)
  years = seq(min(cases$year, na.rm = TRUE), max(cases$year, na.rm = TRUE))
  
  # Loop through districts
  for (district in districts){
    
    # Loop through years
    for (year in years){
      
      # Check if it is missing
      pattern = sprintf("%s_%s", district, year)
      matches = grep(pattern, cases$county_year)
      # If so, add a zero row
      if (length(matches) == 0){
        zero.row = c(district, year, 0, sprintf("%s_%s", district, year))
        cases = rbind(cases, zero.row)
      }
    }
  }
  
  # Remove county_year field
  cases$county_year = NULL
  
  return(cases)
}


#' Convert weather.data to NEdat object
#' 
#' @param weather.data Input weather data in format district, doy, tmeanc, pr, date, plus other unused columns.
#' @param cofips.lookup A link between cofips and district (using the pop object, to prevent needing another input object)
#' 
#' @return NEdat Output weather data in format: cofips, year, month, temp, ppt
#' 
flm.convert.weather.data = function(weather.data, cofips.lookup){
  
  #**# Need to confirm that pr column is actually precipitation
  
  # Add month column from date #This is slow - should we try to speed it up by converting the date to character?
  month = sapply(weather.data$date, substr, 6, 7)
  
  # Aggregate temperature and precipitation data to the corresponding month #**# Should it be mean or sum? Probably mean, which would be more robust to missing values and differing numbers of days across months
  NEtemp = stats::aggregate(weather.data$tmeanc, by = list(weather.data$year, month, weather.data$district), mean, na.rm = TRUE)
  NEtemp$join.col = sprintf("%s_%s_%s", NEtemp$Group.1, NEtemp$Group.2, NEtemp$Group.3)
  colnames(NEtemp) = c("year", "month", "district", "temp1", 'join.col')
  # Aggregate precipitation data
  NEppt = stats::aggregate(weather.data$pr, by = list(weather.data$year, month, weather.data$district), mean, na.rm = TRUE)
  NEppt$join.col = sprintf("%s_%s_%s", NEppt$Group.1, NEppt$Group.2, NEppt$Group.3)
  colnames(NEppt) = c("year", "month", "district", "ppt1", 'join.col')  
  NEdat = merge(NEtemp, NEppt, by = 'join.col')
  
  # Replace district with the corresponding cofips
  NEdat$County = NEdat$district.x
  NEdat = merge(NEdat, cofips.lookup, by = 'County')
  
  # Add required columns in the appropriate order
  NEdat$year = NEdat$year.x
  NEdat$month = as.numeric(as.character(NEdat$month.x))
  NEdat$temp = NEdat$temp1
  NEdat$ppt = NEdat$ppt1
  
  # Drop all extraneous columns
  NEdat$County = NULL
  NEdat$join.col = NULL
  NEdat$year.x = NULL
  NEdat$month.x = NULL
  NEdat$district.x = NULL
  NEdat$temp1 = NULL
  NEdat$year.y = NULL
  NEdat$month.y = NULL
  NEdat$district.y = NULL
  NEdat$ppt1 = NULL
  
  return(NEdat)
}

#' Create a lookup from the population data
#' 
#' @noRd
generate.lookup = function(pop){
  
  #Get districts
  districts = unique(pop$County)
  
  # Create a data frame to hold the lookup file
  cofips.lookup = data.frame(County = rep(NA, length(districts)), cofips = NA)
  
  # Pull cofips from a corresponding county row
  for (i in 1:length(districts)){
    district = districts[i]
    # Get all instances of county
    instances = grep(district, pop$County)
    if (length(instances) == 0){
      stop("Something went wrong with generating the cofips lookup table")
      #**# Consider if this should be a warning, or if each model should be wrapped in a try statement
    }
    # Cofips should be the same for all instances of the county, so just use the first instance of the county
    # Use the first one to get the corresponding cofips
    cofip = pop$cofips[instances[1]]
    cofips.lookup$County[i] = district
    cofpis.lookup$cofips[i] = cofip
  }
  
  return(cofips.lookup)
}

