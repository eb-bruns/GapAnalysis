
# to create species distribution models (SDM) I forked a repo from Dan Carver
# (dcarver1): CWR-of-the-USA-Gap-Analysis

# an alternate or additional way to create an estimated species distribution
# (non-predictive!) in raster form, is to place buffers around the taxon's 
# wild distribution points then rasterize that layer. here we do that!...

################################################################################
# Load libraries
################################################################################

my.packages <- c('raster','tidyverse','terra')
#install.packages(my.packages) #Turn on to install current versions
lapply(my.packages, require, character.only=TRUE)

################################################################################
# Set working directory
################################################################################

# assign main working directory
source("/Users/emily/Documents/GitHub/SDBG_CWR-trees-gap-analysis/0-set_working_directory.R")

# set up additional file paths
path.pts <- file.path(main_dir,occ_dir,"standardized_occurrence_data","taxon_edited_points")
path.rastbuff <- file.path(main_dir,"rasterized_buffers")
  
################################################################################
# Load functions
################################################################################

# source function for filtering occurrence points based on specific flags
source("/Users/emily/Documents/GitHub/SDBG_CWR-trees-gap-analysis/filter_points.R")


# function to create aggregated buffers around occ points, then rasterizing;
# this can be used instead of a species distribution model (SDM)
rasterized.buffer <- function(df,radius,pt_proj,buff_proj,boundary,
                              template_raster){
  # turn occurrence point data into a SpatVector
  spat_pts <- vect(df, geom=c("decimalLongitude","decimalLatitude"), crs=pt_proj)
  # reproject to specified projection
  proj_df <- project(spat_pts,buff_proj)
  # place buffer around each point, then dissolve into one polygon
  buffers <- buffer(proj_df,width=radius)
  buffers <- aggregate(buffers,dissolve = TRUE)
  # clip by boundary so they don't extend into the water
  boundary <- project(boundary,buff_proj)
  buffers_clip <- crop(buffers,boundary)
  # make into simple feature object type
  buffers_clip_sf <- sf::st_as_sf(buffers_clip)
  # rasterize the buffers
  rasterized_buffer <- fasterize::fasterize(buffers_clip_sf, template_raster)
  # trim raster extent to area with values
  rasterized_buffer_trim <- raster::trim(rasterized_buffer,padding=5)
  # return raster
  return(rasterized_buffer_trim)
}

################################################################################
# Read in data for all target taxa
################################################################################


### Protected areas raster 
# **we use this as a template when creating rasterized buffers
# read in PA file downloaded from Dataverse (https://dataverse.harvard.edu/dataverse/GapAnalysis)
protected_areas <- raster(file.path(main_dir,gis_dir,"wdpa_reclass.tif"))


### Ecoregions shapefile 
# **we use this as a boundary to clip the buffers to land
## global ecoregions (TNC terrestrial ecoregions)
eco_global <- vect(file.path(main_dir,gis_dir,
                                    "terr-ecoregions-TNC/tnc_terr_ecoregions.shp"))
# crop to target regions, to make a little smaller
unique(eco_global$WWF_REALM2)
eco_global <- eco_global[eco_global$WWF_REALM2 == "Nearctic" | 
                         eco_global$WWF_REALM2 == "Neotropic",]
eco_global <- project(eco_global, "epsg:4326")
## north america ecoregions (EPA level 1, includes CA and MX)
# global ecoregions don't have separate designation for great lakes, so we
# will get lakes from the north america one
eco_na <- vect(file.path(main_dir,gis_dir,
                         "na_cec_eco_l1/NA_CEC_Eco_Level1.shp"))
# keep just water (great lakes)
eco_water <- eco_na[eco_na$NA_L1NAME == "WATER"]
eco_water <- project(eco_water, "epsg:4326")
## now clip out great lakes from global ecoregions and create aggregated 
# version for clipping buffers
eco_global <- erase(eco_global,eco_water)
boundary.poly <- aggregate(eco_global, dissolve = TRUE)
## save this for use in other scripts
writeVector(boundary.poly,file.path(main_dir,gis_dir,"NorthAm_land_boundary",
                                    "NorthAm_land_boundary.shp"))
#plot(boundary.poly)
rm(eco_global,eco_na,eco_water)


### Occurrence data
# read in flagged occurrence data (ready for spatial analyses)
occ_files <- list.files(path.pts, pattern = ".csv", full.names = T)
occ_dfs <- lapply(occ_files, read.csv)
all_occ <- Reduce(rbind, occ_dfs); rm(occ_files,occ_dfs)
# read in file with manual edits to occurrence data (flagging additional 
# bad points, etc)
manual_pt_edits <- read.csv(
  file.path(main_dir, taxa_dir, "manual_point_edits.csv"),
  na.strings = c("NA",""), colClasses = "character")
# filter points based on flagging columns and manual edits
all_occ <- filter.points(all_occ,manual_pt_edits)

### Select projection to use for points
pt.proj <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"


################################################################################
# Create and save rasterized buffers for each taxon
################################################################################


### Create taxon list from occurrence data
taxa <- sort(unique(all_occ$taxon_name_accepted))


### Create rastertized buffers around occurrence points, to use as 
# species distribution layer in analyses
# can test at multiple buffer sizes (eg 20km, 50km, 100km -- change as desired)
buffer_sizes <- c(20000,50000,100000)

for(i in 1:length(buffer_sizes)){
  buffer_now <- buffer_sizes[i]
  # create output folder for buffer size
  out_fldr <- paste0((buffer_sizes[i]/1000),"km")
  out_dir <- file.path(path.rastbuff,out_fldr)
  if(!dir.exists(out_dir)) dir.create(out_dir, recursive=T)
  # cycle through each target taxon
  for(j in 1:length(taxa)){
    # select taxon occurrence data
    taxon_occ <- all_occ %>% filter(taxon_name_accepted == taxa[j])
    taxon_nm <- gsub(" ","_",taxon_occ$taxon_name_accepted[1])
    # create rasterized buffers
    rast_buff <- rasterized.buffer(taxon_occ,buffer_now,pt.proj,pt.proj,
                                   boundary.poly,protected_areas)
    # save to output directory
    file_nm <- file.path(
      out_dir,paste0(taxon_nm,"-rasterized_buffers_",out_fldr,".tif"))
    writeRaster(rast_buff, file_nm, format="GTiff", overwrite=T)
    print(paste0("Created rasterized buffers for ",taxon_nm," at ",out_fldr))
  }
}


