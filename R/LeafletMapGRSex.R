

create.buffers <- function(df,radius,pt_proj,buff_proj,boundary,
                           lat.col,long.col){
  # turn occurrence point data into a SpatVector
  spat_pts <- vect(df, geom=c(long.col, lat.col), crs=pt_proj)
  # reproject to specified projection
  proj_df <- project(spat_pts,buff_proj)
  # place buffer around each point, then dissolve into one polygon
  buffers <- buffer(proj_df,width=radius)
  buffers <- aggregate(buffers,dissolve = TRUE)
  # clip by boundary so they don't extend into the water
  boundary <- project(boundary,buff_proj)
  buffers_clip <- crop(buffers,boundary)
  # make into object that can be mapped in leaflet
  buffers_clip_sf <- sf::st_as_sf(buffers_clip)
  # return buffer polygons
  return(buffers_clip_sf)
}

# map when all elements are available
create.full.map <- function(spp.sdm,spp.pts.buffer,spp.pts,
                            ex.pts.buffer,ex.pts){
    # color of SDM
    sdm.pal <- colorNumeric("#0e6187",values(spp.sdm),na.color="transparent")
    ## create leaflet map
    map <- leaflet(width = "100%") %>%
    # Base layer groups
      addProviderTiles(providers$Esri.WorldGrayCanvas,group = "Esri.WorldGrayCanvas",
                       options = providerTileOptions(maxZoom = 10)) %>%
      addProviderTiles(providers$OpenStreetMap,group = "OpenStreetMap",
                       options = providerTileOptions(maxZoom = 10)) %>%
      addProviderTiles(providers$Esri.WorldTopoMap,group = "Esri.WorldTopoMap",
                       options = providerTileOptions(maxZoom = 10)) %>%
      addProviderTiles(providers$CartoDB.Positron,group = "CartoDB.Positron",
                       options = providerTileOptions(maxZoom = 10)) %>%
      ### Overlay groups (can toggle)
      addRasterImage(spp.sdm, colors=sdm.pal, opacity = 0.8,
                 group = "Potential distribution") %>%
      # In situ buffers
      addPolygons(data = spp.pts.buffer,
                fillColor = "#965da7", fillOpacity = 0.45,
                weight = 0, opacity = 0, color = "white",
                smoothFactor = 0,
                group = "Estimated distribution") %>%
      # Ex situ buffers
      addPolygons(data = ex.pts.buffer,
                  fillColor = "#45B320", fillOpacity = 0.7, #80bf30
                  weight = 0, color = "white", opacity = 0,
                  smoothFactor = 0,
                  group = "Potential distribution conserved ex situ") %>%
      # In situ points
      addCircleMarkers(data = spp.pts, ~longitude, ~latitude,
                     popup = ~paste0(
                      "<b>Accepted taxon name:</b> ",species,"<br/>",
                      "<b>Verbatim taxon name:</b> ",taxon_name,"<br/>",
                      "<b>Source database:</b> ",database,"<br/>",
                      "<b>All databases with record:</b> ",all_source_databases,"<br/>",
                      "<b>Year:</b> ",year,"<br/>",
                      "<b>Basis of record:</b> ",basisOfRecord,"<br/>",
                      "<b>Dataset name:</b> ",datasetName,"<br/>",
                      "<b>Establishment means:</b> ",establishmentMeans,"<br/>",
                      "<b>Latitude:</b> ",latitude,"<br/>",
                      "<b>Longitude:</b> ",longitude,"<br/>",
                      "<b>Coordinate uncertainty:</b> ",coordinateUncertaintyInMeters,"<br/>",
                      "<b>ID:</b> ",UID),
                     color = "#56116b", radius = 1, fillOpacity = 0.9, #stroke = T,
                     group = "Occurrence points") %>%
      # Ex situ points
      addCircleMarkers(data = ex.pts, ~longitude, ~latitude,
                     popup = ~paste0(
                      "<b>Accepted taxon name:</b> ",species,"<br/>",
                      "<b>Verbatim taxon name:</b> ",taxon_name,"<br/>",
                      "<b>Source database:</b> ",database,"<br/>",
                      "<b>Institution name:</b> ",datasetName,"<br/>",
                      "<b>Accession number:</b> ",references,"<br/>",
                      "<b>Collection date:</b> ",year,"<br/>",
                      "<b>Number of individuals:</b> ",individualCount,"<br/>",
                      "<b>Coordinate uncertainty:</b> ",coordinateUncertaintyInMeters,"<br/>",
                      "<b>ID:</b> ",UID),
                     color = "#77d111", radius = 4, fillOpacity = 0.9, #stroke = T,
                     group = "Wild source localities for ex situ material") %>%
      ### Legend and notes
      addLegend(labels = c(
        "Potential distribution",
        "Occurrence points",
        "Estimated distribution",
        "Wild source localities for ex situ material",
        "Potential distribution conserved ex situ"),
        colors = c("#0e6187","#56116b","#965da7","#77d111","#45B320"),
        #title = "Map of GRS ex situ",
        position = "topright", opacity = 0.8) %>%
      #addControl(
      #  "Use the layers iconon the left to change<br/>the basemap and toggle layers on/off",
      #  position = "topright") %>%
      addControl(
        "Click each point for more information",
        position = "topright") %>%
     ### Layers control
      addLayersControl(
        baseGroups = c("Esri.WorldGrayCanvas","OpenStreetMap",
                       "Esri.WorldTopoMap","CartoDB.Positron"),
        overlayGroups = c(
          "Potential distribution",
          "Occurrence points",
          "Estimated distribution",
          "Wild source localities for ex situ material",
          "Potential distribution conserved ex situ"),
        options = layersControlOptions(collapsed = TRUE),
        position = "topleft") %>%
      #addControl(
      #  "Please do not share this map outside the North American Fruit and Nut Tree Crop Wild Relatives Working Group",
      #  position = "bottomleft") %>%
      setView(-98, 40, zoom = 4)
  return(map)
}

# may when no ex situ points
create.noex.map <- function(spp.sdm,spp.pts.buffer,spp.pts){

  # color of SDM
  sdm.pal <- colorNumeric("#0e6187",values(spp.sdm),na.color="transparent")
  ## create leaflet map
  map <- leaflet(width = "100%") %>%
    # Base layer groups
    addProviderTiles(providers$Esri.WorldGrayCanvas,group = "Esri.WorldGrayCanvas",
                     options = providerTileOptions(maxZoom = 10)) %>%
    addProviderTiles(providers$OpenStreetMap,group = "OpenStreetMap",
                     options = providerTileOptions(maxZoom = 10)) %>%
    addProviderTiles(providers$Esri.WorldTopoMap,group = "Esri.WorldTopoMap",
                     options = providerTileOptions(maxZoom = 10)) %>%
    addProviderTiles(providers$CartoDB.Positron,group = "CartoDB.Positron",
                     options = providerTileOptions(maxZoom = 10)) %>%
    ### Overlay groups (can toggle)
    addRasterImage(spp.sdm, colors=sdm.pal, opacity = 0.8,
                   group = "Potential distribution") %>%
    # In situ buffers
    addPolygons(data = spp.pts.buffer,
                fillColor = "#965da7", fillOpacity = 0.45,
                weight = 0, opacity = 0, color = "white",
                smoothFactor = 0,
                group = "Estimated distribution") %>%
    # In situ points
    addCircleMarkers(data = spp.pts, ~longitude, ~latitude,
                     popup = ~paste0(
                       "<b>Accepted taxon name:</b> ",species,"<br/>",
                       "<b>Verbatim taxon name:</b> ",taxon_name,"<br/>",
                       "<b>Source database:</b> ",database,"<br/>",
                       "<b>All databases with record:</b> ",all_source_databases,"<br/>",
                       "<b>Year:</b> ",year,"<br/>",
                       "<b>Basis of record:</b> ",basisOfRecord,"<br/>",
                       "<b>Dataset name:</b> ",datasetName,"<br/>",
                       "<b>Establishment means:</b> ",establishmentMeans,"<br/>",
                       "<b>Latitude:</b> ",latitude,"<br/>",
                       "<b>Longitude:</b> ",longitude,"<br/>",
                       "<b>Coordinate uncertainty:</b> ",coordinateUncertaintyInMeters,"<br/>",
                       "<b>ID:</b> ",UID),
                     color = "#56116b", radius = 1, fillOpacity = 0.9, #stroke = T,
                     group = "Occurrence points") %>%
    ### Legend and notes
    addLegend(labels = c(
      "Potential distribution",
      "Occurrence points",
      "Estimated distribution",
      "Wild source localities for ex situ material",
      "Potential distribution conserved ex situ"),
      colors = c("#0e6187","#56116b","#965da7","#77d111","#45B320"),
      #title = "Map of GRS ex situ",
      position = "topright", opacity = 0.8) %>%
    #addControl(
    #  "Use the layers iconon the left to change<br/>the basemap and toggle layers on/off",
    #  position = "topright") %>%
    addControl(
      "Click each point for more information",
      position = "topright") %>%
    ### Layers control
    addLayersControl(
      baseGroups = c("Esri.WorldGrayCanvas","OpenStreetMap",
                     "Esri.WorldTopoMap","CartoDB.Positron"),
      overlayGroups = c(
        "Potential distribution",
        "Occurrence points",
        "Estimated distribution"),
      options = layersControlOptions(collapsed = TRUE),
      position = "topleft") %>%
    #addControl(
    #  "Please do not share this map outside the North American Fruit and Nut Tree Crop Wild Relatives Working Group",
    #  position = "bottomleft") %>%
    setView(-98, 40, zoom = 4)
  return(map)
}

# map when no ex situ points and no SDM
create.noex.nosdm.map <- function(eco_clip,spp.pts.buffer,spp.pts,pro_areas){
  ## create leaflet map
  map <- leaflet(width = "100%") %>%
    # Base layer groups
    addProviderTiles(providers$Esri.WorldGrayCanvas,group = "Esri.WorldGrayCanvas",
                     options = providerTileOptions(maxZoom = 10)) %>%
    addProviderTiles(providers$OpenStreetMap,group = "OpenStreetMap",
                     options = providerTileOptions(maxZoom = 10)) %>%
    addProviderTiles(providers$Esri.WorldTopoMap,group = "Esri.WorldTopoMap",
                     options = providerTileOptions(maxZoom = 10)) %>%
    addProviderTiles(providers$CartoDB.Positron,group = "CartoDB.Positron",
                     options = providerTileOptions(maxZoom = 10)) %>%
    ### Overlay groups (can toggle)
    # In situ buffers
    addPolygons(data = spp.pts.buffer,
                fillColor = "#965da7", fillOpacity = 0.45,
                weight = 0, opacity = 0, color = "white",
                smoothFactor = 0,
                group = "Estimated distribution") %>%
    # In situ points
    addCircleMarkers(data = spp.pts, ~longitude, ~latitude,
                     popup = ~paste0(
                       "<b>Accepted taxon name:</b> ",species,"<br/>",
                       "<b>Verbatim taxon name:</b> ",taxon_name,"<br/>",
                       "<b>Source database:</b> ",database,"<br/>",
                       "<b>All databases with record:</b> ",all_source_databases,"<br/>",
                       "<b>Year:</b> ",year,"<br/>",
                       "<b>Basis of record:</b> ",basisOfRecord,"<br/>",
                       "<b>Dataset name:</b> ",datasetName,"<br/>",
                       "<b>Establishment means:</b> ",establishmentMeans,"<br/>",
                       "<b>Latitude:</b> ",latitude,"<br/>",
                       "<b>Longitude:</b> ",longitude,"<br/>",
                       "<b>Coordinate uncertainty:</b> ",coordinateUncertaintyInMeters,"<br/>",
                       "<b>ID:</b> ",UID),
                     color = "#56116b", radius = 1, fillOpacity = 0.9, #stroke = T,
                     group = "Occurrence points") %>%
    ### Legend and notes
    addControl(
      "There were not enough points to model a potential distribution<br>All analyses in this report use the estimated distribution instead",
      position = "topright") %>%
    addLegend(labels = c(
      "Occurrence points",
      "Estimated distribution",
      "Wild source localities for ex situ material",
      "Potential distribution conserved ex situ"),
      colors = c("#56116b","#965da7","#77d111","#45B320"),
      #title = "Map of GRS ex situ",
      position = "topright", opacity = 0.8) %>%
    #addControl(
    #  "Use the layers iconon the left to change<br/>the basemap and toggle layers on/off",
    #  position = "topright") %>%
    addControl(
      "Click each point for more information",
      position = "topright") %>%
    ### Layers control
    addLayersControl(
      baseGroups = c("Esri.WorldGrayCanvas","OpenStreetMap",
                     "Esri.WorldTopoMap","CartoDB.Positron"),
      overlayGroups = c(
        "Occurrence points",
        "Estimated distribution"),
      options = layersControlOptions(collapsed = TRUE),
      position = "topleft") %>%
    #addControl(
    #  "Please do not share this map outside the North American Fruit and Nut Tree Crop Wild Relatives Working Group",
    #  position = "bottomleft") %>%
    setView(-98, 40, zoom = 4)
  return(map)
}

# map when no SDM
create.nosdm.map <- function(eco_clip,spp.pts.buffer,spp.pts,ex.pts.buffer,
                             ex.pts,pro_areas){
  map <- leaflet(width = "100%") %>%
      # Base layer groups
    addProviderTiles(providers$Esri.WorldGrayCanvas,group = "Esri.WorldGrayCanvas",
                     options = providerTileOptions(maxZoom = 10)) %>%
    addProviderTiles(providers$OpenStreetMap,group = "OpenStreetMap",
                     options = providerTileOptions(maxZoom = 10)) %>%
    addProviderTiles(providers$Esri.WorldTopoMap,group = "Esri.WorldTopoMap",
                     options = providerTileOptions(maxZoom = 10)) %>%
    addProviderTiles(providers$CartoDB.Positron,group = "CartoDB.Positron",
                     options = providerTileOptions(maxZoom = 10)) %>%
    ### Overlay groups (can toggle)
    # In situ buffers
    addPolygons(data = spp.pts.buffer,
                fillColor = "#965da7", fillOpacity = 0.45,
                weight = 0, opacity = 0, color = "white",
                smoothFactor = 0,
                group = "Estimated distribution") %>%
    # Ex situ buffers
    addPolygons(data = ex.pts.buffer,
                fillColor = "#45B320", fillOpacity = 0.7, #80bf30
                weight = 0, color = "white", opacity = 0,
                smoothFactor = 0,
                group = "Potential distribution conserved ex situ") %>%
    # In situ points
    addCircleMarkers(data = spp.pts, ~longitude, ~latitude,
                     popup = ~paste0(
                       "<b>Accepted taxon name:</b> ",species,"<br/>",
                       "<b>Verbatim taxon name:</b> ",taxon_name,"<br/>",
                       "<b>Source database:</b> ",database,"<br/>",
                       "<b>All databases with record:</b> ",all_source_databases,"<br/>",
                       "<b>Year:</b> ",year,"<br/>",
                       "<b>Basis of record:</b> ",basisOfRecord,"<br/>",
                       "<b>Dataset name:</b> ",datasetName,"<br/>",
                       "<b>Establishment means:</b> ",establishmentMeans,"<br/>",
                       "<b>Latitude:</b> ",latitude,"<br/>",
                       "<b>Longitude:</b> ",longitude,"<br/>",
                       "<b>Coordinate uncertainty:</b> ",coordinateUncertaintyInMeters,"<br/>",
                       "<b>ID:</b> ",UID),
                     color = "#56116b", radius = 1, fillOpacity = 0.9, #stroke = T,
                     group = "Occurrence points") %>%
    # Ex situ points
    addCircleMarkers(data = ex.pts, ~longitude, ~latitude,
                     popup = ~paste0(
                       "<b>Accepted taxon name:</b> ",species,"<br/>",
                       "<b>Verbatim taxon name:</b> ",taxon_name,"<br/>",
                       "<b>Source database:</b> ",database,"<br/>",
                       "<b>Institution name:</b> ",datasetName,"<br/>",
                       "<b>Accession number:</b> ",references,"<br/>",
                       "<b>Collection date:</b> ",year,"<br/>",
                       "<b>Number of individuals:</b> ",individualCount,"<br/>",
                       "<b>Coordinate uncertainty:</b> ",coordinateUncertaintyInMeters,"<br/>",
                       "<b>ID:</b> ",UID),
                     color = "#77d111", radius = 4, fillOpacity = 0.9, #stroke = T,
                     group = "Wild source localities for ex situ material") %>%
    ### Legend and notes
    addControl(
      "There were not enough points to model a potential distribution<br>All analyses in this report use the estimated distribution instead",
      position = "topright") %>%
    addLegend(labels = c(
      "Occurrence points",
      "Estimated distribution",
      "Wild source localities for ex situ material",
      "Potential distribution conserved ex situ"),
      colors = c("#56116b","#965da7","#77d111","#45B320"),
      #title = "Map of GRS ex situ",
      position = "topright", opacity = 0.8) %>%
    #addControl(
    #  "Use the layers iconon the left to change<br/>the basemap and toggle layers on/off",
    #  position = "topright") %>%
    addControl(
      "Click each point for more information",
      position = "topright") %>%
    ### Layers control
    addLayersControl(
      baseGroups = c("Esri.WorldGrayCanvas","OpenStreetMap",
                     "Esri.WorldTopoMap","CartoDB.Positron"),
      overlayGroups = c(
        "Occurrence points",
        "Estimated distribution",
        "Wild source localities for ex situ material",
        "Potential distribution conserved ex situ"),
      options = layersControlOptions(collapsed = TRUE),
      position = "topleft") %>%
    #addControl(
    #  "Please do not share this map outside the North American Fruit and Nut Tree Crop Wild Relatives Working Group",
    #  position = "bottomleft") %>%
    setView(-98, 40, zoom = 4)
  return(map)
}
