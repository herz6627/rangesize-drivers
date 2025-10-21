
"mcp" <- function(xy, percent=95, unin=c("m", "km"), unout=c("ha", "km2", "m2"))
{
  
  library(tidyverse)
  library(sf)
  library(rnaturalearth)
  library(maps)
  
  ## Verifications
  if (!inherits(xy, "sf"))
    stop("xy should be of class sf")
  if (ncol(st_coordinates(xy))>2)
    stop("xy should be defined in two dimensions")
  pfs <- st_crs(xy)
  
  if (length(percent)>1)
    stop("only one value is required for percent")
  if (percent>100) {
    warning("The MCP is estimated using all relocations (percent>100)")
    percent<-100
  }
  #unin <- match.arg(unin)
  #unout <- match.arg(unout)
  
  
  if (inherits(xy, "sf")) {
    if (!"species" %in% colnames(xy)) {
      warning("xy should have 'species' column, species ignored")
      id <- factor(rep("a", nrow(as.data.frame(xy))))
    } else {
      id <- xy[["species"]]
    }
  } else {
    id <- factor(rep("a", nrow(as.data.frame(xy))))
  }
  
  if (percent>100) {
    warning("The MCP is estimated using all relocations (percent>100)")
    percent<-100
  }
  
  
  
  if (min(table(id))<5)
    stop("At least 5 relocations are required to fit an home range")
  id<-factor(id)
  
  xy <- as.data.frame(st_coordinates(xy))
  
  
  ## Computes the centroid of points for each taxon
  r<-split(xy, id)
  est.cdg<-function(xy) apply(xy, 2, mean)
  cdg<-lapply(r,est.cdg)
  levid<-levels(id)
  
  res <- lapply(1:length(r), function(i) {
    k<-levid[i]
    df.t<-r[[levid[i]]]
    cdg.t<-cdg[[levid[i]]]
    
    ## Distances from points to the centroid: we keep
    ## the "percent" closest
    dist.cdg <-function(xyt) {
      d<-sqrt( ( (xyt[1]-cdg.t[1])^2 ) + ( (xyt[2]-cdg.t[2])^2 ) )
      return(d)
    }
    
    di<-apply(df.t, 1, dist.cdg)
    key<-c(1:length(di))
    
    acons<-key[di<=quantile(di,percent/100)]
    xy.t<-df.t[acons,]
    
    ## Coordinates of the MCP
    coords.t<-chull(xy.t[,1], xy.t[,2])
    xy.bord<-xy.t[coords.t,]
    xy.bord <- rbind(xy.bord[nrow(xy.bord),], xy.bord)
    

    polygon <-  cbind(xy.bord$X, xy.bord$Y) %>% # https://stackoverflow.com/questions/49266736/clip-spatial-polygon-by-world-map-in-r # need to reattach them together for some reason
      st_linestring() %>%
      st_cast("POLYGON") %>%
      st_sfc(crs = "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84", check_ring_dir = TRUE) %>% # crs = 4328 = WGS84
      st_sf()
    
    land <- rnaturalearth::ne_countries(returnclass = "sf", continent = c("North America", "South America")) %>%
      #filter(name %in% c("United States of America", "Canada", "Mexico")) %>% 
      st_union() %>% 
      st_transform(crs = "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84")

    clipped_range <- st_intersection(polygon, land)
    clipped_range$sp <- k
    
    #print(paste("done with i=", i))

    return(clipped_range)
  })
  
  are <- unlist(lapply(1:length(res), function(i) {
    st_area(res[[i]])
  }))

  if (unin == "m") {
    if (unout == "ha")
      are <- are/10000
    if (unout == "km2")
      are <- are/1e+06
  }
  if (unin == "km") {
    if (unout == "ha")
      are <- are * 100
    if (unout == "m2")
      are <- are * 1e+06
  }

  df <- data.frame(area = are, sp = levid)
  
  res <- list(res, df)
  names(res) <- c("polygons", "areas")

  print('clipping all points not over terrestrial America, incliuding islands')
  return(res)
}

