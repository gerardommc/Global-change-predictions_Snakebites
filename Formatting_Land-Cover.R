#Formatting data for projecting snake models
library(raster); library(rgdal)
library(doParallel)

#Land cover
#Original categories:
#Forest           0
#Degraded         1
#Agriculture      2
#Urban            3
#Tea              4
#Paddy to ldv     5

#Land cover data

proj <- raster("pop_2010.asc")

ssps <- paste0("SSP", c(1, 2, 5))
gcms <- c("CNRM-CM5", "GFDL-CM3", "MPI-ESM-LR")

lc.files <- list.files("Land cover", ".asc", full.names = T, recursive = T)
lc.r <- lapply(lc.files, function(x){
      X <- raster(x)
      proj4string(X) <- CRS(proj4string(proj))
      X <- projectRaster(X, crs = CRS(SRS_string="EPSG:5235"), method = "ngb")
      return(X)
})
lc <- stack(lc.r) + 1
lc.s <- reclassify(lc, rcl = matrix(c(6, 2), ncol = 2))

folders <- substr(lc.files, 1, nchar(lc.files) - 15)
un.fold <- unique(folders)
sapply(un.fold, function(x){dir.create(paste0(x, "tif"))})

for(i in 1:nlayers(lc.s)){
   writeRaster(lc.s[[i]], paste0(folders[i], "tif", 
                                 substr(lc.files[i], nchar(lc.files[i]) - 15,
                                                     nchar(lc.files[i]) - 4)),
               "GTiff", overwrite= T)
   }

#Computing the variables we need to project the snake models
#Variables needed: Distance to forests,proportion of agricultural land

forest <- lc == 1
forest <- reclassify(forest, rcl = matrix(c(0, NA), ncol = 2))

library(doParallel); registerDoParallel(cores = 4)
dist.forest <- foreach(i = 1:nlayers(forest))%dopar%{distance(forest[[i]])}
dist.forest <- stack(dist.forest)
dist.forest <- mask(dist.forest, lc)

#Writing up distance to forests data
dir.create("Distance-forests")
sapply(ssps, function(x){dir.create(paste0("Distance-forests", "/", x))})

fold.dist <- paste0("Distance-forests", substr(un.fold, 11, nchar(un.fold)))
sapply(fold.dist, function(x){dir.create(x)})

fold.dist.files <- rep(fold.dist, each = 41)

years <- rep(2010:2050, 41)

for(i in seq_along(fold.dist.files)){
   writeRaster(dist.forest[[i]], 
               paste0(fold.dist.files[i], "Dist-for-", years[i]),
               "GTiff", overwrite = T)
}

##Proportion of agriculture
agr <- lc == 3
agr.5 <- aggregate(agr, 5, "mean")
agr.1 <- resample(agr.5, agr, method = "bilinear")
prop.agr <- mask(agr.1, agr)

#Writing up distance to proportion of agricultural land
dir.create("Proportion-agriculture")
sapply(ssps, function(x){dir.create(paste0("Proportion-agriculture", "/", x))})

fold.agr <- paste0("Proportion-agriculture", substr(un.fold, 11, nchar(un.fold)))
sapply(fold.agr, function(x){dir.create(x)})

fold.agr.files <- rep(fold.agr, each = 41)

for(i in seq_along(fold.dist.files)){
   writeRaster(prop.agr[[i]], 
               paste0(fold.agr.files[i], "Prop-agr-", years[i]),
               "GTiff", overwrite = T)
}

########################
####Tree cover data#####
tc.files <- list.files("Tree cover", "asc", full.names = T, recursive = T)
tc.paths <- substr(tc.files, 1, nchar(tc.files) - 15)

tc <- stack(tc.files)
proj4string(tc) <- CRS(proj4string(proj))
tc <- (tc - 5) * (-1 / 5)
#tc <- tc * 90

#Computing a correction factor
tc.2010 <- readRDS("Data-objects/Snakes-env-data.rds")
tc.2010.r <- rasterFromXYZ(tc.2010[,c("x", "y", "tree")])
proj4string(tc.2010.r) <- CRS("+init=epsg:5235")
tc.2010.r <- tc.2010.r/cellStats(tc.2010.r, max)

tc <- projectRaster(tc, crs = CRS("+init=epsg:5235"), method = "bilinear")
tc <- resample(tc, tc.2010.r, method = "bilinear")

#Computng correction factor
year <- rep(2010:2050, 7)

tc.drop <- which(year != 2010)
tc.2010.r1 <- dropLayer(tc, tc.drop)
tc.2010.r1 <- reclassify(tc.2010.r1, rcl = matrix(c(0.96, 1, 0.99), ncol = 3, byrow = T))

tc.2010.logit <- log(tc.2010.r1/(1-tc.2010.r1)) #Logit transform of model projections for data-years

tc.2010.l <- log(tc.2010.r/(1 - tc.2010.r)) #Logit transform of data-years

tc.corr.factor <- mean(tc.2010.l - tc.2010.logit) #Average difference between projections and data in logit scale
   
## Logit transformation of all projections
tc.res <- reclassify(tc, rcl = matrix(c(0.96, 1, 0.99), ncol = 3, byrow = T))
tc.logit <- log(tc.res/(1 - tc.res))

corr <- tc.logit - tc.corr.factor
tc.corr <- exp(corr)/(1 + exp(corr))
zeroes <- tc.2010.r > 100
tc.corr <- merge(zeroes, tc.corr) * 100

tc.path.uniq <- unique(tc.paths)

sapply(tc.path.uniq, function(x){dir.create(paste0(x, "tif"))})

registerDoParallel(cores = 8)
foreach(i = seq_along(tc.paths)) %dopar% {
   writeRaster(tc.corr[[i]],
               paste0(tc.paths[i], "/tif/",
                      "Tree-cover-", year[i]),
               "GTiff", overwrite = T)
}



#Resampling land cover to 5km

ref.5k <- readRDS("Data-objects/Incidence models data-Apr2020.rds")
ref.id <- 1:nrow(ref.5k)
ref.id.r <- rasterFromXYZ(data.frame(ref.5k[, c("x", "y")], ref.id))

library(dplyr)

registerDoParallel(cores = 4)
lc.5k <- foreach(i = 1:nlayers(lc.s)) %dopar% {
      lr <- stack(lapply(1:5, function(x){return(lc.s[[i]] == x)}))
      l.df <- data.frame(rasterToPoints(lr))
      names(l.df) <- c("x", "y", "a", "b", "c", "d", "e")
      
      zones <- extract(ref.id.r, l.df[, c("x", "y")])
      l.df$zones <- zones
      
      l.df.5k <- l.df %>% group_by(zones) %>% summarise(a = mean(a, na.rm = T),
                                                        b = mean(b, na.rm = T),
                                                        c = mean(c, na.rm = T),
                                                        d = mean(d, na.rm = T),
                                                        e = mean(e, na.rm = T))
      l.5k <- data.frame(l.df.5k[, c("a", "b", "c", "d", "e")])
      l.5k.s <- apply(l.5k, 1, which.max)
      
      l.5k.r <- rasterFromXYZ(data.frame(ref.5k[, c("x", "y")], lc = l.5k.s))
      return(l.5k.r)
}
lc.5k <- stack(lc.5k)

un.fold.5k <- paste0("Land-cover-5km",substr(un.fold, 11, nchar(un.fold)))
dir.create("Land-cover-5km")

for(i in 1:nlayers(lc.5k)){
   writeRaster(lc.5k[[i]], paste0("Land-cover-5km", substr(folders[i], 11, nchar(folders[i])), 
                                 substr(lc.files[i], nchar(lc.files[i]) - 15,
                                        nchar(lc.files[i]) - 4)),
               "GTiff", overwrite = T)
}
