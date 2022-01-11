####################################################################################################
##
## Script for cleaning possibly overlapping multi-state histories according to a specified hierarchy
## using Norwegian registry data. Individual inclusion dates may be according to some specified 
## criteria occurring between an overall start- and stop-date for the study.
##
## Input: 
## ## dt: a data.table with the columns: id, start, stop and state
##    # id: unique individual cluster identifier (numeric)
##    # start: start date of the state (YYYYMMDD) (numeric)
##    # stop: stop date of the state (YYYYMMDD) (numeric)
##    # state: name of the recorded state (character)
## ## incl.date: inclusion date, overall start of study (YYYYMMDD) (numeric)
## ## end.date: study end date
## ## write_file_name: name of file to write cleaned data to (character)
## ## absorbing: name of absorbing states (character or numeric)
##
## Output: a matrix file with 5 columns: id, start.day, stop.day, state.num, id.incl.date
##
####################################################################################################
library(data.table); library(zoo)

#### Input ####
load("mydata.Rdata")
write_file_name <- "outputfilepath"
incl.date <- 20030101
end.date <- 20101231
absorbing <- "dead"
## Number and merging states
dt$stat <- 0
dt[state == "educ"]$stat <- 1
dt[state == "unem"]$stat <- 2
dt[state == "sohe"]$stat <- 2
dt[state == "work"]$stat <- 3
dt[state == "perm"]$stat <- 3
dt[state == "sick"]$stat <- 4
dt[state == "reha"]$stat <- 4
dt[state == "disa"]$stat <- 4
dt[state == "dead"]$stat <- 5

#### Create empty file to write cleaned state histories to ####
write(NULL, file = paste0(write_file_name), ncolumns =  5,
      append = FALSE, sep = " ")

#### Set missing stop dates to the next start date and remove stays of zero length ####
dt[, start.lead1:= shift(start, type = "lead"), by = id]
dt[is.na(stop), stop:= start.lead1]
dt <- dt[stop > start | state %in% absorbing]

#### Truncate data according to inclusion and end date ####
dt <- dt[!(state %in% absorbing & stop > end.date)]
dt <- dt[stop > incl.date]
dt[start < incl.date, start:= incl.date]
dt <- dt[start < end.date]
dt <- dt[stop >  end.date, stop:= end.date]

#### Convert date variables to real dates ####
dt[, start:= as.Date(as.character(start), format = "%Y%m%d")]
dt[, stop:= as.Date(as.character(stop), format = "%Y%m%d")]
dt <- dt[order(id, start, stop)]

#### Add individual inclusion date and start/stop days since individual inclusion
w <- dt[, .I[1], by = id][[2]]
dt[w, id.incl.date:= start]
dt$id.incl.date <- na.locf(dt$id.incl.date)
dt$start.day <- as.numeric(dt$start - dt$id.incl.date)
dt$stop.day <- as.numeric(dt$stop - dt$id.incl.date)
dt[, id.incl.date:= as.numeric(format(format(id.incl.date, "%Y%m%d")))] # faster write
dt <- dt[, list(id, start.day, stop.day, stat, id.incl.date)]

#### Remove individuals without any stays of length > 0
dt[, max.out:= max(stop.day), by = id]
dt <- dt[max.out > 0]

#### Loop over individuals. Roughly 5 minutes per 1000 ids. Do parallelisation as needed. ####
ids <- unique(dt$id)
options(warn = 2) #stop loop on warnings
for(j in 1:length(ids)) {
  print(j)
  dti <- dt[id==ids[j]]
  #### Vector for storing states. 0: intermediate censoring, 6: right censoring
  msvec <- c(rep(0, dti[1]$max.out + 1), 6)
  #### State Hierarchy. From least to most important (will overwrite entries in msvec)
  ## 1
  Start <- t(dti[stat==1, "start.day"]) + 1
  Stop  <- t(dti[stat==1, "stop.day"]) + 1
  msvec[unlist(mapply(":", Start, Stop))] <- 1
  ## 2:
  Start <- t(dti[stat==2, "start.day"]) + 1
  Stop  <- t(dti[stat==2, "stop.day"]) + 1
  msvec[unlist(mapply(":", Start, Stop))] <- 2
  ## 3:
  Start <- t(dti[stat==3, "start.day"]) + 1
  Stop  <- t(dti[stat==3, "stop.day"]) + 1
  msvec[unlist(mapply(":", Start, Stop))] <- 3
  ## 4:
  Start <- t(dti[stat==4, "start.day"]) + 1
  Stop  <- t(dti[stat==4, "stop.day"]) + 1
  msvec[unlist(mapply(":", Start, Stop))] <- 4
  ## 5 (absorbing):
  Start <- t(dti[stat==5, "start.day"]) + 1
  msvec[Start] <- 5
  ## Limit msvec when censored or absorbed
  msvec <- msvec[1:min(match(c(5,6), msvec), na.rm = rm)]
  ## Create matrix entry for id j
  rlemsvec <- rle(msvec)
  days <- rlemsvec$lengths
  stat <- rlemsvec$values
  start <- c(0,cumsum(days[-length(days)]) + 1)
  stop <- cumsum(days)
  id <- ids[j]
  id.incl.date <- dt[id==ids[j], id.incl.date][1]
  dti <- cbind(id, start, stop, stat, id.incl.date)
  #### Write matrix to file
  write(t(dti),
        file = paste0(write_file_name),
        ncolumns =  5,
        append = TRUE, sep = " ")
}

