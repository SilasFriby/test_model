


#########################################################################
############                HELPFUL FUNCTIONS                ############
#########################################################################



is.wholenumber <- function(x, tol = .Machine$double.eps^0.5)  abs(x - round(x)) < tol


timeFirstDayInHalfYear <- function (charvec, format = "%Y-%m-%d", zone = "", FinCenter = ""){
  if (zone == "") 
    zone = getRmetricsOptions("myFinCenter")
  if (FinCenter == "") 
    FinCenter = getRmetricsOptions("myFinCenter")
  charvec = timeFirstDayInMonth(charvec = charvec, format = format, 
                                FinCenter = FinCenter)
  lt = strptime(charvec, format, tz = "GMT")
  first.quarter = rep(c(1, 7), each = 6) - 1
  lt$mon = first.quarter[1 + lt$mon]
  timeDate(format(lt), format = "%Y-%m-%d", zone = zone, FinCenter = FinCenter)
}

timeFirstDayInYear <- function (charvec, format = "%Y-%m-%d", zone = "", FinCenter = ""){
  if (zone == "") 
    zone = getRmetricsOptions("myFinCenter")
  if (FinCenter == "") 
    FinCenter = getRmetricsOptions("myFinCenter")
  charvec = timeFirstDayInMonth(charvec = charvec, format = format, 
                                FinCenter = FinCenter)
  lt = strptime(charvec, format, tz = "GMT")
  first.quarter = rep(c(1), each = 12) - 1
  lt$mon = first.quarter[1 + lt$mon]
  timeDate(format(lt), format = "%Y-%m-%d", zone = zone, FinCenter = FinCenter)
}





## add years to date
addToYearsDate <- function(startDate, years){
  date     <- as.Date(startDate)
  tmp      <- as.POSIXlt(date)
  tmp$year <- tmp$year + years
  return(as.Date(tmp))
  
}



## get the first date of the follwing period
firstDateOfNextPeriod <- function(date, period){
  
  if (period == 'month'){
    date <- timeFirstDayInMonth(as.Date(date))
    tmp      <- as.POSIXlt(date)
    tmp$mon  <- tmp$mon + 1
  }
  
  if (period == 'quarter'){
    date <- timeFirstDayInQuarter(as.Date(date))
    tmp      <- as.POSIXlt(date)
    tmp$mon  <- tmp$mon + 3
  }
  
  if (period == '6 months'){
    date <- timeFirstDayInHalfYear(as.Date(date))
    tmp      <- as.POSIXlt(date)
    tmp$mon  <- tmp$mon + 6
  }
  
  if (period == 'year'){
    date <- timeFirstDayInYear(as.Date(date))
    tmp      <- as.POSIXlt(date)
    tmp$mon  <- tmp$mon + 12
  }

  return(as.Date(tmp))  
}






#' Compounding rate
#'
#' @param numberInitial - starting number (numeric) 
#' @param dt - step length (numeric)
#' @param rateDt - projection rate (numeric)
#' @param n - number of time points (integer)
#'
#' @return project numbers on the 'n' time points (vector)

addRate <- function(numberInitial, dt, rateDt, n) {

  numberAdjusted    <- rep(NA, n + 1)
  numberAdjusted[1] <- numberInitial
  
  if( n == 0 ) {
    return(numberInitial)
  } else {
    for (i in 2:(n + 1)) {
      numberAdjusted[i] <- numberAdjusted[i - 1] * (1 + rateDt[i - 1]) ^ dt
      
    }
    return(numberAdjusted)
  }
  
}







#' benefit projections with annual rate regulations
#'
#' @param startAnnualBenefit - initial annual benefit (numeric) 
#' @param dates - vector of benfit dates (vector)
#' @param dtDates - step size between benefit dates (numeric)
#' @param rateRegulationAnnual - annual rate regulation (numeric) 
#' @param monthUpdate - month where the rateRegulationAnnual is added (numeric)
#' @param pct - set different than one if benefits are computed as a percentage of a pension base (numeric)
#'
#' @return data frame containing benefit dates and the corresponding regulated benefits

projectBenefitsFct <- function(startAnnualBase, dates, dtDates, rateRegulationAnnual, monthUpdate = 1, pct = 1) {
  
  
  ## project annual benefits
  
  numberOfYears          <- length( unique(format(dates, '%Y')) )  # number of years 
  annualBase             <- addRate(numberInitial = startAnnualBase, 
                                    dt            = 1, 
                                    rateDt        = rateRegulationAnnual, 
                                    n             = numberOfYears )
  annualBenefits         <- annualBase * pct # if benefits are computed as a percentage of a pension base
  
  
  ## add year column
  
  numberOfDates           <- length(dates)
  yearFirstBenefit        <- as.numeric(format(dates[1], '%Y'))
  yearLastBenefit         <- as.numeric(format(dates[numberOfDates], '%Y')) + 1 
  yearColumn              <- seq(yearFirstBenefit, yearLastBenefit)
  annualBenefits          <- data.frame('year'    = yearColumn,
                                        'benefit' = annualBenefits)
  
  ## project benefits per dtBenefit
  
  benefitsPerDt <- rep(NA, numberOfDates)
  for (i in 1:numberOfDates) {
    
    yearDate                  <- as.numeric(format(dates[i], '%Y')) # year of loop date
    monthDate                 <- as.numeric(format(dates[i], '%m')) # month of loop date
    
    # different benefits before and after 'monthUpdate'
    # if (i == 1) {
    #   fooAnnualBenefit          <- annualBenefits$benefit[annualBenefits$year == yearDate]
    # } else 
    
    # different benefits before and after 'monthUpdate' - split cases betwwen the one where the first date is before and after 'monthUpdate'
    if (as.numeric(format(dates[1], '%m')) < monthUpdate) {
      
      if (monthDate < monthUpdate & monthUpdate != 1) {
        fooAnnualBenefit          <- annualBenefits$benefit[annualBenefits$year == yearDate]
      } else if (monthDate >= monthUpdate & monthUpdate != 1) {
        fooAnnualBenefit          <- annualBenefits$benefit[annualBenefits$year == (yearDate + 1)]
      } else {
        fooAnnualBenefit          <- annualBenefits$benefit[annualBenefits$year == yearDate] # if 'monthUpdate' is January
      } 
   
    }
    
    if (as.numeric(format(dates[1], '%m')) >= monthUpdate) { # more tricky than the first case
      
      if (annualBenefits$year[1] == yearDate & monthUpdate != 1) {
        fooAnnualBenefit          <- annualBenefits$benefit[annualBenefits$year == yearDate]
      } else if (monthDate < monthUpdate & annualBenefits$year[1] == (yearDate - 1) & monthUpdate != 1) {
        fooAnnualBenefit          <- annualBenefits$benefit[annualBenefits$year == (yearDate - 1)] 
      } else if (monthDate < monthUpdate & monthUpdate != 1) {
        fooAnnualBenefit          <- annualBenefits$benefit[annualBenefits$year == yearDate] 
      } else if (monthDate >= monthUpdate & monthUpdate != 1) {
        fooAnnualBenefit          <- annualBenefits$benefit[annualBenefits$year == (yearDate + 1)]
      } else {
        fooAnnualBenefit          <- annualBenefits$benefit[annualBenefits$year == yearDate] # if 'monthUpdate' is January
      } 
      
    }
    
    benefitsPerDt[i]          <- fooAnnualBenefit * dtDates
    
  }
  
  benefitsPerDt         <- data.frame('date'    = dates,
                                      'benefit' = benefitsPerDt)
  
  return(benefitsPerDt)
  
} 





#' Convert data.frame object to markdown table
#'
#' @param inFrame - data.frame
#'
#' @return markdown table

tableCat <- function(inFrame) {
  
  outText <- paste(names(inFrame), collapse = " | ")
  outText <- c(outText, paste(rep("---", ncol(inFrame)), collapse = " | "))
  invisible(apply(inFrame, 1, function(inRow) {
    outText <<- c(outText, paste(inRow, collapse = " | "))
  }))
  return(outText)
  
}



#' From integer to month
#'
#' @param monthInteger - month number (numeric) 
#'
#' @return a string containing the corresponding month

fromIntegerToMonth <- function(monthInteger) {
  
  months <- c('January', 'February', 'March', 'April', 'May', 'June', 
              'July', 'August', 'September', 'October', 'November', 'December')
  
  return(months[monthInteger])
  
}





