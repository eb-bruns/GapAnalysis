#' @title Final conservation score ex situ
#' @name FCSex
#' @description This function calculates the average of the three ex situ conservation metrics
#'   returning a final conservation score summary table. It also assigns conservation priority categories
#' @param Occurrence_data A data frame object with the species name, geographical coordinates,
#'  and type of records (G or H) for a given species
#' @param Species_list A vector of characters with the species names to calculate the GRSex metrics.
#' @param Raster_list A list of rasters representing the species distribution models for the species list provided
#'  in \var{Species_list}. The order of rasters in this list must match the same order as \var{Species_list}.
#' @param Buffer_distance Geographical distance used to create circular buffers around germplasm.
#'  Default: 50000 (50 km) around germplasm accessions (CA50)
#' @param Ecoregions_shp A shapefile representing Ecoregions information with a field ECO_ID_U representing Ecoregions Ids.
#'  If Ecoregions=NULL the function will use a shapefile provided for use after running GetDatasets()
#' @param Gap_Map logical, if \code{TRUE} the function will calculate gap maps for each species analyzed and will return a list
#'  with three slots: FCSex, GRSex_maps,and ERSex_maps

#' @return This function returns a data frame summarizing the ex-situ gap analysis scores:
#'
#' \tabular{lcc}{
#' species \tab Species name \cr
#' SRSex \tab Sampling representativeness score ex situ \cr
#' GRSex \tab Geographical representativeness score ex situ \cr
#' ERSex \tab Ecological representativeness score ex situ \cr
#' FCSex \tab Final conservation score ex situ  \cr
#' }
#'
#' @examples
#' ##Obtaining occurrences from example
#' data(CucurbitaData)
#' ##Obtaining species names from the data
#' Cucurbita_splist <- unique(CucurbitaData$species)
#' ##Obtaining raster_list
#' data(CucurbitaRasters)
#' CucurbitaRasters <- raster::unstack(CucurbitaRasters)
#' ##Obtaining ecoregions shapefile
#' data(ecoregions)
#' #Running all three Ex-situ gap analysis steps using a unique function
#' FCSex_df <- FCSex(Species_list=Cucurbita_splist,
#'                                       Occurrence_data=CucurbitaData,
#'                                       Raster_list=CucurbitaRasters,
#'                                       Buffer_distance=50000,
#'                                       Ecoregions_shp=ecoregions,
#'                                       Gap_Map=TRUE)
#'
#'@references
#'
#' Khoury et al. (2019) Ecological Indicators 98:420-429. doi: 10.1016/j.ecolind.2018.11.016
#'
#' @export


FCSex <- function(Species_list, Occurrence_data, Raster_list,
                  Buffer_distance=50000, Ecoregions_shp=NULL, Gap_Map=FALSE,
                  #EBB: adding option for separate dataset and filtering
                  #     column for SRS calculation
                  Occurrence_data_raw=Occurrence_data
                  #,Select_database="GBIF"
                  ){

  SRSex_df <- NULL
  GRSex_df <- NULL
  ERSex_df <- NULL
  FCSex_df <- NULL

  #Checking Occurrence_data format
  par_names <- c("species","latitude","longitude","type")

  if(missing(Occurrence_data)){
    stop("Please add a valid data frame with columns: species, latitude, longitude, type")
  }

  if(isFALSE(identical(names(Occurrence_data[,1:4]),par_names))){
    stop("Please format the column names in your dataframe as species, latitude, longitude, type")
  }

  # Load in ecoregions shp
  if(is.null(Ecoregions_shp) | missing(Ecoregions_shp)){
    if(file.exists(system.file("data/preloaded_data/ecoRegion/tnc_terr_ecoregions.shp",
                               package = "GapAnalysis"))){
      Ecoregions_shp <- raster::shapefile(system.file("data/preloaded_data/ecoRegion/tnc_terr_ecoregions.shp",
                                                      package = "GapAnalysis"),encoding = "UTF-8")
    } else {
      stop("Ecoregions file is not available yet. Please run the function GetDatasets() and try again")
    }
  } else{
    Ecoregions_shp <- Ecoregions_shp
  }


  #Checking if Gap_Map option is a boolean or if the parameter is missing left Gap_Map as FALSE
  if(is.null(Gap_Map) | missing(Gap_Map)){ Gap_Map <- FALSE
  } else if(isTRUE(Gap_Map) | isFALSE(Gap_Map)){
    Gap_Map <- Gap_Map
  } else {
    stop("Choose a valid option for GapMap (TRUE or FALSE)")
  }


    #EBB: if raw dataset is provided, filter points to only one database,
    #     to avoid many duplicates
#  par_names_raw <- c("species","latitude","longitude","type","database")
#  if(isFALSE(identical(names(Occurrence_data_raw),par_names_raw))){
#    print("If you'd like to use a different occurrence point dataframe for the SRSex calculation, pass it into the Occurrence_data_raw variable and format as species, latitude, longitude, type, database")
#  } else {
#    print(paste0("For the SRSex calculation, filtering raw points to include those from ",Select_database," only. If you'd like to filter by a different database, pass it into the Select_database variable"))
#    Occurrence_data_raw <- Occurrence_data_raw %>%
#      filter(database == Select_database | database == "Ex_situ") %>%
#      select(-database)
#  }
    #EBB: source edits...
  #source("/Users/emily/Documents/GitHub/GapAnalysis/R/SRSex.R")
  SRSex_df <- SRSex(Species_list = Species_list,
                    Occurrence_data_raw <- Occurrence_data_raw,
                    Occurrence_data = Occurrence_data)
    #EBB: source edits...
  #source("/Users/emily/Documents/GitHub/GapAnalysis/R/GRSex.R")
  GRSex_df <- GRSex(Occurrence_data = Occurrence_data,
                    Species_list = Species_list,
                    Raster_list = Raster_list,
                    Buffer_distance = Buffer_distance,
                    Gap_Map = Gap_Map)
    #EBB: source edits...
  #source("/Users/emily/Documents/GitHub/GapAnalysis/R/ERSex.R")
  ERSex_df <- ERSex(Species_list = Species_list,
                    Occurrence_data = Occurrence_data,
                    Raster_list = Raster_list,
                    Buffer_distance = Buffer_distance,
                    Ecoregions_shp=Ecoregions_shp,
                    Gap_Map = Gap_Map)

  # join the dataframes based on species

  if(is.data.frame(GRSex_df)){
    FCSex_df <- merge(SRSex_df, GRSex_df, by ="species", all.x = TRUE)
  } else {
    FCSex_df <- merge(SRSex_df, GRSex_df$GRSex, by ="species", all.x = TRUE)
  }
  #EBB: check if df for ERS too
  #FCSex_df <- merge(FCSex_df, ERSex_df$ERSex, by = "species", all.x = TRUE)
  if(is.data.frame(ERSex_df)){
    FCSex_df <- merge(FCSex_df, ERSex_df, by ="species", all.x = TRUE)
  } else {
    FCSex_df <- merge(FCSex_df, ERSex_df$ERSex, by ="species", all.x = TRUE)
  }


  # calculate the mean value for each row to determine fcs per species
  FCSex_df$FCSex <- rowMeans(FCSex_df[, c("SRSex", "GRSex", "ERSex")])

  #assign classes (exsitu)
  FCSex_df$FCSex_class <- with(FCSex_df, ifelse(FCSex < 25, "UP",
                                                ifelse(FCSex >= 25 & FCSex < 50, "HP",
                                                       ifelse(FCSex >= 50 & FCSex < 75, "MP",
                                                              "LP"))))
  FCSex_df <- FCSex_df %>% rename(Taxon = species)


  if(isTRUE(Gap_Map)){
    FCSex_df <- list(FCSex=FCSex_df,GRSex_maps=GRSex_df$gap_maps,ERSex_maps=ERSex_df$gap_maps)
  } else{
    FCSex_df <- FCSex_df
  }

  return(FCSex_df)
}
