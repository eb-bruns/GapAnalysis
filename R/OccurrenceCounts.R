#' @title Generating counts dataframe for taxa
#' @name OccurrenceCounts
#' @description This function creates a summary file counting
#'  the total number of G and H occurrences, including those with coordinates
#'
#' @param species A vector of characters with the species name.
#' @param Occurrence_data A data frame object with the species name, geographical coordinates,
#'  and type of records (G or H) for a given species
#' @return This function returns a data frame object with the following columns:
#'
#' \tabular{lcc}{
#'  species \tab Species name \cr
#'  totalRecords  \tab Total number of records  \cr
#'  hasLat  \tab Number of occurrences with latitude \cr
#'  hasLon  \tab Number of occurrences with longitude \cr
#'  totalUseful  \tab Number of occurrences with coordinates \cr
#'  totalGRecords  \tab Number of G occurrences \cr
#'  totalGUseful  \tab Number of G occurrences with coordinates \cr
#'  totalHRecords  \tab Number of H occurrences \cr
#'  totalHUseful  \tab Number of H occurrences with coordinates \cr
#'  }
#'
#' @examples
#' data(CucurbitaData)
#' ##Obtaining species names from the data
#' Cucurbita_splist <- unique(CucurbitaData$species)
#' sp_counts <- OccurrenceCounts(Cucurbita_splist[[1]],CucurbitaData)
#'
#'@references
#'
#' Ramirez-Villegas et al. (2010) PLOS ONE, 5(10), e13497. doi: 10.1371/journal.pone.0013497
#' Khoury et al. 2019) Ecological Indicators 98: 420-429. doi: 10.1016/j.ecolind.2018.11.016
#'
#' @export
#' @keywords internal


OccurrenceCounts <- function(species,
                             # EBB: need to add raw occurrence data due to data format
                             Occurrence_data_raw,
                             Occurrence_data){

  taxon <- NULL
  latitude <- NULL
  longitude <- NULL
  hasLat <- NULL
  hasLong <- NULL
  hasLatLong <- NULL
  type <- NULL

  #Checking Occurrence_data format
  par_names <- c("species","latitude","longitude","type")

  if(isFALSE(identical(names(Occurrence_data[,1:4]),par_names))){
    stop("Please format the column names in your dataframe as species, latitude, longitude, type")
  }
  if(isFALSE(identical(names(Occurrence_data_raw[,1:4]),par_names))){
    stop("Please format the column names in your dataframe as species, latitude, longitude, type")
  }

  # create an empty dataframe to store counts information
    ##EBB: updating to reflect our workflow with points
  df <- data.frame(matrix(NA, nrow = 1, ncol = 7))
  colNames <- c("Taxon",
                "TotalRecords","TotalCoords",
                "TotalExsituRecords","TotalExsituCoords",
                "TotalRefRecords","TotalRefCoords")
  colnames(df) <- colNames

  #EBB: all species occurrences before filtering for mapping
  speciesOccRaw <- Occurrence_data_raw[which(Occurrence_data_raw$species==species),]

  #speciesOcc$hasLat <- !is.na(speciesOcc$latitude) &
  #  speciesOcc$latitude != "\\N" & speciesOcc$latitude != "" &
  #  !is.null(speciesOcc$latitude) & speciesOcc$latitude != "NULL"

  #speciesOcc$hasLong <- !is.na(speciesOcc$longitude) &
  #  speciesOcc$longitude != "\\N" & speciesOcc$longitude != "" &
  #  !is.null(speciesOcc$longitude) & speciesOcc$longitude != "NULL"

  #speciesOcc$hasLatLong <- speciesOcc$hasLat==speciesOcc$hasLong & speciesOcc$hasLat==TRUE

  tblRaw <- stats::aggregate(speciesOccRaw,list(type = speciesOccRaw$type), length)
  tblRaw <- tblRaw[,c("type","species")]; colnames(tblRaw)[2] <- "total"

  #EBB: species occurrences after filtering for mapping
  speciesOcc <- Occurrence_data[which(Occurrence_data$species==species),]

  tbl <- stats::aggregate(speciesOcc,list(type = speciesOcc$type), length)
  tbl <- tbl[,c("type","species")]; colnames(tbl)[2] <- "total"

   # assign values to the counts dataframe for the species
  df$Taxon <- gsub("_"," ",as.character(species))
  df$TotalRecords <- nrow(speciesOccRaw)
  df$TotalCoords <- nrow(speciesOcc)
  df$TotalExsituRecords <- sum((subset(tblRaw, type == "G"))$total)
  df$TotalExsituCoords <- sum((subset(tbl, type == "G"))$total)
  df$TotalRefRecords <- sum((subset(tblRaw, type == "H"))$total)
  df$TotalRefCoords <- sum((subset(tbl, type == "H"))$total)
#  df$hasLat <- sum(speciesOcc$hasLat)
#  df$hasLong <- sum(speciesOcc$hasLong)

  # returns the counts dataframe
  return(df)
}
