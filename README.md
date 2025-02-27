## My forked version of the GapAnalysis R package, with small changes to fix issues and bigger changes to reformat the summary HTML based on our project. You can find results at https://NorthAmericanFruitNutTreeCWR.github.io/. Original README below --



# GapAnalysis R package

## Description
The GapAnalysis R package evaluates the ex situ and in situ conservation status of taxa, combines these metrics into an integrated  assessment, and calculates an indicator metric across taxa. GapAnalysis generates quantitative and spatial outputs which demonstrate the state of conservation as well as where gaps in protection exist. The methods are fully described in Carver et al. (2021). Articles by Ramirez-Villegas et al. (2010), Castañeda-Álvarez and Khoury et al. (2016), and Khoury et al. (2019a, b; 2020) describe the main steps toward the current methodology.

The GapAnalysis functions require the user to provide two inputs: a `data.frame` of species occurrences, and a `raster` object of the predicted habitat (species distribution model) for each assessed taxon.

This library consists of 12 functions within 4 families: pre-analysis, ex situ conservation gap analysis, in situ conservation gap analysis, and summary evaluations. In short, the pre-analysis process establishes the file structure and prepares the input data. The ex situ and in situ processes perform the respective conservation strategy gap analyses and produce both quantitative and spatial results. The combined assessment merges the individual assessments, summarizes the results across taxa, calculates the indicator, and generates a summary html document for each taxon, which can be used to evaluate outputs and aid conservation planning.

## Installation
GapAnalysis can be installed as follows
```r
#CRAN
install.packages("GapAnalysis")
#Alternative: GitHub
library(devtools)
remotes::install_github("CIAT-DAPA/GapAnalysis")
```
A full list of libraries needed for the package is included below.

**Dependencies:** `raster`

**Imports:** `base, utils, sp, tmap, data.table, sf, methods, geosphere, data.table, fasterize, rmarkdown`

**Suggests:** `knitr, rgdal, rgeos, kableExtra, DT`


## Usage
We provide the below reproducible example (also available in the package documentation). Please note this example is provided at 10 arc minutes resolution for efficient processing time; the results in the associated published article (Khoury et al. 2019c) and described in the associated R package article (Carver et al. 2021) differ as the analysis was conducted at 2.5 arc minutes resolution. For more details on accessing and utilizing the 2.5 arc minutes dataset see [ecoregions and protected areas](#ecolink).

```r
##Load package
library(raster)
library(GapAnalysis)

##Obtaining occurrences from example
data(CucurbitaData)

##Obtaining species names from the data
speciesList <- unique(CucurbitaData$species)

##Obtaining raster_list
data(CucurbitaRasters)
CucurbitaRasters <- raster::unstack(CucurbitaRasters)

##Obtaining protected areas raster
data(ProtectedAreas)

##Obtaining ecoregions shapefile
data(ecoregions)

#Running all three ex situ gap analysis steps using FCSex function
FCSex_df <- FCSex(Species_list=speciesList,
                  Occurrence_data=CucurbitaData,
                  Raster_list=CucurbitaRasters,
                  Buffer_distance=50000,
                  Ecoregions_shp=ecoregions
)

#Running all three in situ gap analysis steps using FCSin function
FCSin_df <- FCSin(Species_list=speciesList,
                  Occurrence_data=CucurbitaData,
                  Raster_list=CucurbitaRasters,
                  Ecoregions_shp=ecoregions,
                  Pro_areas=ProtectedAreas)

## Combine gap analysis metrics
FCSc_mean_df <- FCSc_mean(FCSex_df = FCSex_df,FCSin_df = FCSin_df)

##Running Conservation indicator across taxa
indicator_df  <- indicator(FCSc_mean_df)

## Generate summary HTML file with all result
GetDatasets()
summaryHTML_file <- SummaryHTML(Species_list=speciesList,
                                Occurrence_data = CucurbitaData,
                                Raster_list=CucurbitaRasters,
                                Buffer_distance=50000,
                                Ecoregions_shp=ecoregions,
                                Pro_areas=ProtectedAreas,
                                Output_Folder=".",
                                writeRasters=FALSE)
```

## Usage with different buffer distances for _ex situ_ gap analysis

```r
#Buffer distances for 5, 10, and 20 km respectively

buffer_distances <- c(5000,10000,20000)

SRSex_df <- SRSex(Species_list = speciesList,
                  Occurrence_data = CucurbitaData)

FCSex_df_list <- list()


#Running all three ex situ gap analysis steps using FCSex function

#Choose if gap maps are calculated for ex situ gap analysis using diferent buffer size
Gap_Map=FALSE

for(i in 1:length(speciesList)){

  FCSex_df_list[[i]] <- FCSex(Species_list=speciesList[i],
                    Occurrence_data=CucurbitaData,
                    Raster_list=CucurbitaRasters[i],
                    Buffer_distance=buffer_distances[i],
                    Ecoregions_shp=ecoregions,
                    Gap_Map=Gap_Map)



};rm(i)

#Returning FCSex object
if(Gap_Map==TRUE){
  FCSex_df <- list(FCSex=do.call(rbind,lapply(FCSex_df_list, `[[`, 1)),
                 GRSex_maps=do.call(c,lapply(FCSex_df_list, `[[`, 2)),
                 ERSex_maps=do.call(c,lapply(FCSex_df_list, `[[`, 3))
                 )
} else {
  FCSex_df <- do.call(rbind,FCSex_df_list)
}


#Running all three in situ gap analysis steps using FCSin function

FCSin_df <- FCSin(Species_list=speciesList,
                  Occurrence_data=CucurbitaData,
                  Raster_list=CucurbitaRasters,
                  Ecoregions_shp=ecoregions,
                  Pro_areas=ProtectedAreas,
                  Gap_Map = NULL)


## Combine gap analysis metrics
FCSc_mean_df <- FCSc_mean(FCSex_df = FCSex_df,FCSin_df = FCSin_df)


##Running Conservation indicator across taxa
indicator_df  <- indicator(FCSc_mean_df)
```

The below sub-sections provide further details on the input data and GapAnalysis steps.

### Data inputs
**_Species occurrences_**

A `data.frame` of species occurrences and record type. This process can handle single or multiple taxa.

species | latitude | longitude | type
------------ | ------------- | -------------| -------------
Cucurbita_cordata | 28.9457 | -113.563 | G
Cucurbita_digitata |  | | H

**species:** this value will be the key for all functions in this library. Ensure it is consistent for all records and is included in the file name of your predicted potential habitat `raster` as well.

**latitude** and **longitude** must be in decimal degrees, preferably with the highest accuracy possible.

**type:** All records must be classified as either a reference observation (typically the main presence data input into the species distribution modeling, labeled H as most records in our previous research source from herbaria), or as a “site of collection” location of an existing ex situ accession from a conservation repository (labeled G, as most records in our previous research source from genebanks). This distinction is significant for multiple evaluations and effort must be taken to ensure the correct assignment of these values.

Digital repositories such as GBIF, EDDmaps, and IDIGBIO contain observed locations of a taxon considered “H” type occurrences in GapAnalysis. From our experience there can be duplication within and between such databases and care should be taken to reduce duplication where possible.

The major sources for G occurrence data that the authors have used in GapAnalysis include [USDA GRIN](https://npgsweb.ars-grin.gov/gringlobal/search.aspx), [GENESYS](https://www.genesys-pgr.org/), [FAO WIEWS](http://www.fao.org/wiews/en/), [PlantSearch](https://tools.bgci.org/plant_search.php), and [GBIF](https://www.gbif.org/) (records designated as 'living specimen'). Generally, data from these sources would be considered a G type if it is still an active and living accession. Duplication between these sources may exist.

More information and examples of how to make the distinction between “H” and “G” points can be found [here](https://doi.org/10.1111/DDI.13008).

<a name="ecolink">
<b><i>Ecoregions and Protected Area </b></i>
</a>
The ecoregion and protected areas datasets are provide through the package via the `GetDatasets()` functions. The files will be downloaded and store at
```r
system.file("data/preloaded_data/ecoRegion/tnc_terr_ecoregions.shp",package = "GapAnalysis")
```
These files can be found accessed directly at the [Dataverse repository](https://dataverse.harvard.edu/dataverse/GapAnalysis) associated with this package.
The original datasets can be found here([ecoregions](http://maps.tnc.org/gis_data.html), [world database of protected areas](https://www.protectedplanet.net/en/thematic-areas/wdpa?tab=WDPA)). The ecoregion dataset is provided in its native vector data type. The package's WDPA layer has been transformed from a vector to a binary raster at 2.5 arc minutes resolution raster.

**_Predicted Habitat_**

The `raster` representing the predicted extent of suitable habitat (species distribution model) is used by multiple functions to represent the maximum potential range of a taxon. This is then compared to what is conserved _ex situ_ or _in situ_. Although a required input, the generation of species distribution models is not included in GapAnalysis because a number of R packages for this process already exist (e.g. packages `sdm`, `wallace`, `dismo` and `maxnet`).


### Workflow
The recommended workflow is as follows:

**Pre-analysis**
 - `GetDatasets` downloads the protected areas and ecoregions datasets from our data repository

**Ex-situ Analysis**
 - `SRSex` calculates the Sampling Representativeness Score for _ex situ_ conservation
 - `GRSex` calculates the Geographic Representativeness Score for _ex situ_ conservation. During this process, an ex situ geographic gap map is also created for each species by subtracting the G buffered areas out of the distribution model of each taxon, leaving only those areas considered not sufficiently sampled for ex situ conservation
 - `ERSex` calculates the Ecological Representativeness Score for _ex situ_ conservation. During this process, an ex situ ecological gap map is also created for each species by mapping only the spatial areas within the distribution model of each taxon which are occupied by ecoregions not represented by G buffers
 - `FCSex` calculates the Final Conservation Score for _ex situ_ conservation as an average of the above 3 scores and assigns a priority category for each taxon based on the final conservation score

**In-situ Analysis**
 - `SRSin` calculates the Sampling Representativeness Score for _in situ_ conservation
 - `GRSin` calculates the Geographic Representativeness Score for _in situ_ conservation. During this process, an in situ geographic gap map is also created for each species by subtracting the protected areas out of the distribution model of each taxon, revealing those areas in the model not currently in protected areas
 - `ERSin` calculates the Ecological Representativeness Score for _in situ_ conservation. During this process, an in situ ecological gap map is also created for each species by by mapping only the spatial areas within the distribution model of each taxon which are occupied by ecoregions not represented at all in protected areas
 - `FCSin` calculates the Final Conservation Score for _in situ_ conservation as an average of the above 3 scores and assigns a priority category for each taxon based on the final conservation score

**Summary evaluations**   
 - `FCSc_mean` computes the mean as well as minimum and maximum of the _ex situ_ and _in situ_ Final Conservation Scores. It also assigns taxa to priority categories based on final conservation scores (high priority (HP) for further conservation action assigned when FCS < 25, medium priority (MP) where 25 ≤ FCS < 50, low priority (LP) where 50 ≤ FCS < 75, and sufficiently conserved (SC) for taxa whose FCS ≥75)
- `indicator` calculates an indicator across assessed taxa, which can be applied at national, regional, global, or any other scale (Khoury et al., 2019). The indicator is calculated separately with regard to ex situ, in situ, min, max, and combined (mean) conservation, by deriving the proportion of taxa categorized as SC or LP out of all taxa.
 - `SummaryHTML` produces a summary HTML output with taxon specific quantitative and spatial results

**Internal functions**
 - `OccurrenceCounts` creates a `data.frame` with counts of G, H, and those record types with coordinates for all taxa, based on input occurrence data
 - `Gbuffer` is an internal function that creates a circular buffer of user-defined size (default is 50 km radius) around each G point for each taxon, which represents the geographic areas already considered to be sufficiently collected for ex situ conservation. The output of this process is a raster. Since this is not an exported function to use it you will need to type `GapAnalysis:::Gbuffer` in R. This is a modified version of [geobuffer_pts.R](https://github.com/valentinitnelav/geobuffer).
 - `ParamTest` checks if occurrence data and distribution models exist for each species and sets the conservation scores to zero if there are no occurrences with coordinates or no model.

Each function can be run as a standalone method and in any order. However, we recommend following this workflow as it will ensure dependencies for individual functions are in place and that the variables are stored correctly to successfully produce the final summary document. For more details on each of these calculations, see the list of references below.

## Authors
Main: Daniel Carver, Chrystian C. Sosa, Colin K. Khoury, and Julian Ramirez-Villegas

Other contributors: Harold A. Achicanoy, Maria Victoria Diaz, Steven Sotelo, Nora P. Castaneda-Alvarez

## References

Carver D, Sosa CC, Khoury CK, Achicanoy HA, Diaz MV, Sotelo S, Castañeda-Álvarez NP, and Ramírez-Villegas JR (2021) GapAnalysis: an R package to calculate conservation indicators using spatial information. Ecography. doi: 10.1111/ecog.05430. (https://doi.org/10.1111/ecog.05430)

Castañeda-Álvarez NP, Khoury CK, Achicanoy HA, Bernau V, Dempewolf H, Eastwood RJ, Guarino L, Harker RH, Jarvis A, Maxted N, Mueller JV, Ramirez-Villegas J, Sosa CC, Struik PC, Vincent H, and Toll J (2016) Global conservation priorities for crop wild relatives. Nature Plants 2(4): 16022. doi: [10.1038/nplants.2016.22](http://www.nature.com/articles/nplants201622)

Khoury CK, Amariles D, Soto JS, Diaz MV, Sotelo S, Sosa CC, Ramirez-Villegas J, Achicanoy HA, Velásquez-Tibata J, Guarino L, Leon B, Navarro-Racines C, Castañeda-Álvarez NP, Dempewolf H, Wiersema JH, and Jarvis A (2019a) Comprehensiveness of conservation of useful wild plants: an operational indicator for biodiversity and sustainable development targets. Ecological Indicators 98: 420-429. doi: [10.1016/j.ecolind.2018.11.016](https://doi.org/10.1016/j.ecolind.2018.11.016)

Khoury CK, Carver D, Barchenger DW, Barboza G, van Zonneweld M, Jarret R, Bohs L, Kantar MB, Uchanski M, Mercer K, Nabhan GP, Bosland PW, and Greene SL (2019b) Modeled distributions and conservation status of the wild relatives of chile peppers (Capsicum L). Diversity and Distributions 26(2): 209-225. doi: 10.1111/DDI.13008. https://doi.org/10.1111/DDI.13008

Khoury CK, Carver D, Greene SL, Williams KA, Achicanoy HA, Schori M, León B, Wiersema JH, and Frances A (2020) Crop wild relatives of the United States require urgent conservation action. Proc Natl Acad Sci USA 117(52): 33351-33357. doi: 10.1073/pnas.2007029117. https://doi.org/10.1073/pnas.2007029117

Khoury CK, Carver D, Kates HR, Achicanoy HA, van Zonneweld M, Thomas E, Heinitz C, Jarret R, Labate JA, Reitsma K, Nabhan GP, and Greene SL (2019c) Distributions, conservation status, and abiotic stress tolerance potential of wild cucurbits (Cucurbita L.). Plants, People, Planet 2(3): 269-283. doi: 10.1002/ppp3.10085. https://doi.org/10.1002/ppp3.10085

Ramirez-Villegas J, Khoury CK, Jarvis A, Debouck DG, Guarino L (2010) A gap analysis methodology for collecting crop genepools: a case study with Phaseolus beans. PLoS One 5, e13497. [doi:10.1371/journal.pone.0013497](http://dx.doi.org/10.1371%2Fjournal.pone.0013497)

## License
GNU GENERAL PUBLIC LICENSE Version 3
