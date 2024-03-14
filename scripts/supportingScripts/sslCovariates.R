getCovariates <- function(extent=c(-2499.76, -999.6641, -0.09183017, 999.9718)){
  
  ssl_site_info <- readr::read_csv("data/ssl_site_info.csv") %>% 
    dplyr::select(PARENTSITE,SUBSITE,REGION,RCA,E_W_144,ROOK,LAT,LON) %>% 
    `colnames<-`(tolower(colnames(.))) %>% 
    st_as_sf(coords=c("lon","lat"), crs=4326) %>% st_transform(3338)
  
  bathy <- getNOAA.bathy(120,-100,30,89,resolution=1.2179,
                         antimeridian = TRUE) %>% marmap::as.raster() %>% terra::rast()
  
  slope <- terra::terrain(bathy)
  bathy <- terra::project(bathy, "epsg:3338") %>% `names<-`("depth")
  slope <- terra::project(slope, "epsg:3338") %>% `names<-`("slope")
  
  # clip to extent
  bathy <- crop(bathy,extent(extent*1000))
  slope <- crop(slope,extent(extent*1000))
  
  Grid_r_c <- terra::xyFromCell(bathy,1:ncell(bathy)) %>% as.data.frame() %>% 
    st_as_sf(coords=c('x','y'), crs=3338)
  dist_site <- bathy[[1]] %>% 
    `values<-`(nngeo::st_nn(Grid_r_c, ssl_site_info, k=1, returnDist=TRUE)$dist %>% unlist()) %>% 
    `names<-`("dist_site")
  
  hab_cov_expand <- stack(
    raster(bathy),
    raster(dist_site),
    raster(slope)
  )
  
  # convert to km ######
  hab_cov_expand$depth <- hab_cov_expand$depth/1000
  hab_cov_expand$dist_site = hab_cov_expand$dist_site/1000
  
  covlist <- list(depth = hab_cov_expand$depth,
                  slope = hab_cov_expand$slope,
                  d2site = hab_cov_expand$dist_site,
                  d2site2 = hab_cov_expand$dist_site^2)
  
  for(i in 1:length(covlist)) {
    extent(covlist[[i]]) <- extent(c(xmin(covlist[[i]]), xmax(covlist[[i]]), 
                                     ymin(covlist[[i]]), ymax(covlist[[i]]))/1000)
    projection(covlist[[i]]) <- gsub("units=m", "units=km", projection(covlist[[i]]))
  }
  
  # rescale depth
  absmin <- abs(min(getValues(covlist$depth),na.rm=TRUE))
  covlist$depth <- covlist$depth+absmin
  covlist[["depthslope"]] <- covlist$depth * covlist$slope
  attr(covlist$depth,"absmin") <- absmin
  
  return(covlist[c("depth","slope","depthslope","d2site","d2site2")])
  
}