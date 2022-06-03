#Projecting snake ppms
library(raster); library(spatstat); library(doParallel)

snakes <- lapply(paste0("Snake-presences/", 
                        c("Bungarus caeruleus", "Bungarus ceylonicus",
                          "Daboia russellii", "Echis carinatus",
                          "Hypnale spp", "Naja naja",
                          "Trimeresurus trigonocephalus"),
                          "-5500.csv"),read.csv)
names(snakes) <- c("Bungarus caeruleus", "Bungarus ceylonicus",
                   "Daboia russelii", "Echis carinatus",
                   "Hypnale spp", "Naja naja",
                   "Trimeresurus trigonocephalus")
arb.df <- readRDS("Data-objects/Snakes-env-data.rds")

topo <- rasterFromXYZ(arb.df[, c("x", "y", "topo")])
roads <- rasterFromXYZ(arb.df[, c("x", "y", "roads")])
roads <- roads/roads; roads <- roads != 1

coords <- data.frame(rasterToPoints(topo))

#Spatstat formatting
ux = sort(unique(coords$x))
uy = sort(unique(coords$y))
nx = length(ux)
ny = length(uy)
col.ref = match(coords$x, ux)
row.ref = match(coords$y, uy)
all.vec = rep(NA, max(row.ref)*max(col.ref))
vec.ref = (col.ref - 1)*max(row.ref) + row.ref
all.vec[vec.ref] = 1
Lka.mask = matrix(all.vec, max(row.ref), max(col.ref), dimnames = list(uy, ux))

#Setting up the study window

Lka.win = as.owin(im(Lka.mask, xcol = ux, yrow = uy))

#Preparing presences for the quadrature scheme

X <- lapply(snakes, function(x){coordinates(x)[,1]})
Y <- lapply(snakes, function(x){coordinates(x)[,2]})

ppp.dat <- foreach(i = seq_along(snakes)) %do%{
      ppp(X[[i]], Y[[i]], window = Lka.win, check = FALSE)
} 

quads <- ppp(arb.df$x, arb.df$y, window = Lka.win)

Q <- foreach(i = seq_along(snakes)) %do% {
      quadscheme(data = ppp.dat[[i]], dummy = quads, method = "grid",
                 ntile = c(nx, ny), npix = c(nx, ny))
}

stackToImList <- function(x){
      X <- with(x, 
                    cbind(roads,
                          topo,
                          tc, pa, df,
                          Bungarus_caeruleus,
                          Bungarus_ceylonicus,
                          Daboia_russelii, 
                          Echis_carinatus,
                          Hypnale_spp,
                          Naja_naja,
                          Trimeresurus_trigonocephalus
                    )
      ) #There's an issue with the number of pixels to fill

      im.list <- list()
      for (i in 1:dim(X)[2]){
            all.vec = rep(NA, max(row.ref)*max(col.ref))
            vec.ref = (col.ref - 1)*max(row.ref) + row.ref
            all.vec[vec.ref] = X[,i]
            im.list[[i]] = im(matrix(all.vec, max(row.ref), max(col.ref),
                                          dimnames = list(uy, ux)),
                              xcol = ux, yrow = uy)
      }
      names(im.list) <- c("dist.cit.roads",
                          "topo",
                          "tree",
                          "prop.agric",
                          "dist.forest",
                          "Bungarus_caeruleus.SLD99",
                          "Bungarus_ceylonicus.SLD99",
                          "Daboia_russelii.SLD99",
                          "Echis_carinatus.SLD99",
                          "Hypnale_hypnale.SLD99",
                          "Naja_naja.SLD99",
                          "Trimeresurus_trigonocephalus.SLD99")
      return(im.list)
}

## SSP and RCP setup
ssp <- c("SSP1", "SSP2", "SSP5")
rcp <- c("NA", "RCP 4.5", "RCP 8.5")
gcm <- c("CNRM-CM5", "GFDL-CM3", "MPI-ESM-LR")
years <- 2010:2050

#Fixed variables


ref <- topo
proj4string(ref) <- CRS("+init=epsg:5235")

proj4string(roads) <- CRS("+init=epsg:5235")
roads <- resample(roads, ref)

coords.ref <- data.frame(rasterToPoints(ref))[, c("x", "y")]

dir.create("Snakes")

snake.models <- paste0("../SNAKEBITE MODELLING/Snakes Fundamental niches/PPMs/AIMs-all-species.rds")
ppm <- readRDS(snake.models)
echis <- readRDS("../SNAKEBITE MODELLING/Snakes Fundamental niches/PPMs/Echis-best-ppms.rds")
ppm[[4]] <- echis[[2]]
spp <- c("Bungarus caeruleus",
         "Bungarus ceylonicus",
         "Daboia russelii",
         "Echis carinatus",
         "Hypnale spp",
         "Naja naja", 
         "Trimeresurus trigonocephalus")

# SSP2 & 5
for(i in 2:3){
   dir.create(paste0("Snakes/", ssp[i]))
      
   for (j in 1:3) {
         dir.create(paste0("Snakes/", ssp[i], "/", gcm[j]))

         for(k in seq_along(years)){
              dir.create(paste0("Snakes/", ssp[i], "/", gcm[j], "/", years[k]))
              
              lc <- list(tc = raster(paste0("Tree cover/", ssp[i], "/", gcm[j],  "/tif/Tree-cover-", years[k], ".tif")),
                         df = raster(paste0("Distance-forests/", ssp[i], "/", gcm[j],  "/Dist-for-", years[k], ".tif"))/1000,
                         pa = raster(paste0("Proportion-agriculture/", ssp[i], "/", gcm[j],  "/Prop-agr-", years[k], ".tif")))
              
              lc <- lapply(lc, function(x){
                    proj4string(x) <- CRS("+init=epsg:5235")
                    x <- resample(x, ref)
              })
              
              snakes.1 <- stack(list.files(paste0("../SNAKEBITE MODELLING/Snakes Fundamental niches/Projections/", rcp[i], "/", gcm[j], "/", years[k] - 2), "tif", full.names = T))
              snakes.2 <- stack(list.files(paste0("../SNAKEBITE MODELLING/Snakes Fundamental niches/Projections/", rcp[i], "/", gcm[j], "/", years[k] - 1), "tif", full.names = T))
              snakes.3 <- stack(list.files(paste0("../SNAKEBITE MODELLING/Snakes Fundamental niches/Projections/", rcp[i], "/", gcm[j], "/", years[k] ), "tif", full.names = T))
              
              sn <- (snakes.1 + snakes.2 + snakes.3)/3
              proj4string(sn) <- CRS("+init=epsg:5235")
              sn <- resample(sn, ref)
            
              lc$snakes <- sn
              lc$roads <- roads
              lc$topo <- ref
              
              lc.r <- stack(lc)
              
              lc.df <- data.frame(extract(lc.r, coords.ref))
              
              lc.im <- stackToImList(lc.df)
              
              for(l in seq_along(ppm)){
                           r <- raster(predict(ppm[[l]], covariates = lc.im, type = "trend", ngrid = c(469, 267)))
                           writeRaster(r, paste0("Snakes/", ssp[i], "/", gcm[j], "/", years[k],"/", spp[l]),
                                       "GTiff", overwrite = T)
              }
              
             }
      }
}

## SSP1

sn <- stack(list.files("../SNAKEBITE MODELLING/Snakes Fundamental niches/Niches/Final-dists", "tif", full.names = T, recursive = F))
proj4string(sn) <- CRS("+init=epsg:5235")
sn <- resample(sn, ref)
names(sn) <- c("Bungarus_caeruleus",
               "Bungarus_ceylonicus",
               "Daboia_russelii", 
               "Echis_carinatus",
               "Hypnale_spp",
               "Naja_naja",
               "Trimeresurus_trigonocephalus")

dir.create(paste0("Snakes/", ssp[1]))

for(k in seq_along(years)){
      dir.create(paste0("Snakes/", ssp[1], "/", "/", years[k]))
      
      lc <- list(tc = raster(paste0("Tree cover/", ssp[1], "/tif/Tree-cover-", years[k], ".tif")),
                 df = raster(paste0("Distance-forests/", ssp[1], "/Dist-for-", years[k], ".tif"))/1000,
                 pa = raster(paste0("Proportion-agriculture/", ssp[1],  "/Prop-agr-", years[k], ".tif")))
      
      lc <- lapply(lc, function(x){
         proj4string(x) <- CRS("+init=epsg:5235")
         x <- resample(x, ref)
      })
      
      lc$snakes <- sn
      lc$roads <- roads
      lc$topo <- ref
      
      lc.r <- stack(lc)
      
      lc.df <- data.frame(extract(lc.r, coords.ref))
      
      lc.im <- stackToImList(lc.df)
      
      for(l in seq_along(ppm)){
         r <- raster(predict(ppm[[l]], covariates = lc.im, type = "trend", ngrid = c(469, 267)))
         writeRaster(r, paste0("Snakes/", ssp[1], "/", years[k],"/", spp[l]),
                     "GTiff", overwrite = T)
      }
      
}
