rcp <- c("RCP 4.5", "RCP 8.5")
gcm <- c("CNRM-CM5", "GFDL-CM3", "MPI-ESM-LR")
spp <- c("Bungarus caeruleus", "Bungarus ceylonicus",
         "Daboia russelii", "Echis carinatus",
         "Hypnale spp", "Naja naja",
         "Trimeresurus trigonocephalus")

data <- expand.grid(rcp = rcp, 
                    gcm = gcm, 
                    spp = spp)

library(doParallel)
library(terra)

registerDoParallel(cores = 24)

DNCs <- foreach(i = 1:nrow(data), .combine = rbind) %dopar% {
    f <- paste0("Snake-DNCs/",
                data$rcp[i], "/",
                data$gcm[i], "/",
                2006:2050, "/",
                data$spp[i], ".tif")
  r <- rast(f)
  m <- global(r, function(x){median(x, na.rm = T)})
  iqr.1 <- global(r, function(x){quantile(x, 0.25, na.rm = T)})
  iqr.2 <- global(r, function(x){quantile(x, 0.75, na.rm = T)})
  
  df <- data.frame(spp = data$spp[i],
                   rcp = data$rcp[i],
                   gcm = data$gcm[i],
                   Median = m,
                   IQ.25 = iqr.1,
                   IQ.75 = iqr.2,
                   year = 2006:2050)
  return(df)
}

DNCs <- data.frame(DNCs)

names(DNCs)[4] <- c("Median") 

library(ggplot2)

ggplot(DNCs) + geom_point(aes(x = year, y = log10(Median), colour = gcm, 
                               shape = rcp), size = 0.7) +
  geom_smooth(aes(x = year, y = log10(global), colour = gcm, linetype = rcp, fill = gcm), alpha = 0.2) +
  facet_wrap(~spp)

pdf("Plots/DNC-snakes.pdf", width = 8, height = 8)
ggplot(DNCs) + geom_point(aes(x = year, y = log10(Median), colour = gcm, 
                              shape = rcp), size = 0.7) +
  geom_smooth(aes(x = year, y = log10(Median), colour = gcm, linetype = rcp, fill = gcm), alpha = 0.2) +
  facet_wrap(~spp)
dev.off()
