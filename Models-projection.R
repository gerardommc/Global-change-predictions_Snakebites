
# Projection functions

format.proj.data <- function(file = "", max.s){
      require(readODS)
      l <- readRDS(file)
      idx <- read.csv("../Snakebite-zoonotic-transmission/Data/Agressivenes-indices.csv")
      rel.abund <- read_ods("../Questionnaires/Parameters-questions.ods",
                             sheet = 2)$Density_5k
      
      df <- lapply(l, function(x){
            S <- as.matrix(x[, 3:9])
            S <- sweep(S, 2, max.s, "/")
            S <- sweep(S, 2, rel.abund, "*")
            
            log.hum <- log(x$hum.pop)
            x$log.hum <- log.hum
            
            return(data.frame(hum.pop = x[, "hum.pop"],
                              log.hum = log.hum,
                              land.cover = x[, "land.cover"],
                              S))  
      })
      return(df)
}

sb.func <- function(S, log.hum, hum.pop, land.cover, beta0, beta1, ind, aggr, rho){
   S.1 <- sweep(as.matrix(S), 2, aggr * ind, "*")
   S.sum <- rowSums(S.1)
   Beta <- exp(beta0[land.cover] + beta1[land.cover] * log.hum^2)
   P <- (1 - exp(-Beta * S.sum))*exp(rho)
   return(P)
}

env.func <- function(S, H.bit, hum.pop, log.hum, land.cover, s.effects, int, sever, beta0, beta1, rho){
   S.1 <- sweep(S, 2, sever, "*")
   Beta.i <- exp(beta0[land.cover] + beta1[land.cover] * log.hum^2)
   S.2 <- sweep(S.1, 1, Beta.i, "*")
   df <- data.frame(forest = as.numeric(land.cover==1),
                    degraded = as.numeric(land.cover==2),
                    agriculture = as.numeric(land.cover==3),
                    urban = as.numeric(land.cover==4),
                    tea = as.numeric(land.cover==5), S.2)
   pars <- c(int, s.effects)
   P <- binomial()$linkinv(pars  %*% t(as.matrix(df))) 
   H.env <- c(P) * H.bit
   P.env <- (H.env/hum.pop)*exp(rho)
   return(P.env)
}

snake.env.contrib <- function(S, H.bit, hum.pop, log.hum, land.cover, 
                              s.effects, int, sever, beta0, beta1){
   S.1 <- sweep(S, 2, sever, "*")
   Beta.i <- exp(beta0[land.cover] + beta1[land.cover] * log.hum^2)
   S.2 <- sweep(S.1, 1, Beta.i, "*")
   max.sp <- apply(S.2, 1, which.max)
   return(max.sp)
}

#Incidence data
incid.data <- readRDS("../SNAKEBITE MODELLING/Data objects/Incid-models-data-Aug-2020/Incidence models data-Apr2020.rds")
S <- incid.data[, 3:9]
max.s <- apply(S, 2, max)
#Loading up models

library(coda)
sb.model <- readRDS("../Snakebite-zoonotic-transmission/Model-results/Snakebites-mass-action-NIMBLE-model.rds")
env.model <- readRDS("../Snakebite-zoonotic-transmission/Model-results/Envenoming-NIMBLE-model.rds")

sb.mcmc <- rbind(sb.model[[1]], sb.model[[2]], sb.model[[3]])
env.mcmc <- rbind(env.model[[1]], env.model[[2]], env.model[[3]])

sb.summary <- list(
 beta0 = data.frame(
    median = apply(sb.mcmc[, paste0("beta0[", 1:5, "]")], 2, median),
    ints = t(apply(sb.mcmc[, paste0("beta0[", 1:5, "]")], 2, function(x){
       HPDinterval(as.mcmc(x), 0.95)
    }))
 ),
 beta1 = data.frame(
    median = apply(sb.mcmc[, paste0("beta1[", 1:5, "]")], 2, median),
    ints = t(apply(sb.mcmc[, paste0("beta1[", 1:5, "]")], 2, function(x){
       HPDinterval(as.mcmc(x), 0.95)
    }))
 ),
 ind = data.frame(
    median = apply(sb.mcmc[, paste0("indices[", 1:7, "]")], 2, median),
    ints = t(apply(sb.mcmc[, paste0("indices[", 1:7, "]")], 2, function(x){
       HPDinterval(as.mcmc(x), 0.95)
    }))
 ),
 rho = data.frame(
    median = apply(sb.mcmc[, paste0("rho[", 1:3057, "]")], 2, median),
    ints = t(apply(sb.mcmc[, paste0("rho[", 1:3057, "]")], 2, function(x){
       HPDinterval(as.mcmc(x), 0.95)
    }))
 )
)

env.summary <- list(
   int = data.frame(
      median = apply(env.mcmc[, paste0("inter[", 1:5, "]")], 2, median),
      ints = t(apply(env.mcmc[, paste0("inter[", 1:5, "]")], 2, function(x){
         HPDinterval(as.mcmc(x), 0.95)
      }))
   ),
   s.effects = data.frame(
      median = apply(env.mcmc[, paste0("s.effects[", 1:7, "]")], 2, median),
      ints = t(apply(env.mcmc[, paste0("s.effects[", 1:7, "]")], 2, function(x){
         HPDinterval(as.mcmc(x), 0.95)
      }))
   ),
   rho = data.frame(
      median = apply(env.mcmc[, paste0("rho[", 1:3057, "]")], 2, median),
      ints = t(apply(env.mcmc[, paste0("rho[", 1:3057, "]")], 2, function(x){
         HPDinterval(as.mcmc(x), 0.95)
      }))
   )
)

#Indices
agr.ind <- read.csv("../Snakebite-zoonotic-transmission/Data/Agressivenes-indices.csv")

proj.files <- list.files("Incid-models", "rds", full.names = T)

proj.names <- list.files("Incid-models", "rds", full.names = F)
proj.names <- paste0("Projection", substr(proj.names, 10, nchar(proj.names)))

dir.create("../Population projections/Incidence-projections")

library(doParallel); registerDoParallel(cores = 8)
for(i in seq_along(proj.files)){
   dat <- format.proj.data(proj.files[i], max.s = max.s)
   projection <- foreach(j = seq_along(dat)) %dopar% {
      sb.incid <- foreach(k = 1:3, .combine = cbind) %do% {
         sb.func(S = dat[[j]][, 4:10],
                          log.hum = dat[[j]]$log.hum,
                          hum.pop = dat[[j]]$hum.pop,
                          land.cover = dat[[j]]$land.cover,
                          beta0 = sb.summary$beta0[, k],
                          beta1 = sb.summary$beta1[, k],
                          ind = sb.summary$ind[, k],
                          aggr = agr.ind$Agressiveness/10,
                 rho = sb.summary$rho[, k])
         } #Working fine
      
         env.incid <- foreach(k = 1:3, .combine = cbind) %do% { 
            env.func(S = dat[[j]][, 4:10],
                            H.bit = sb.incid[, 1] * dat[[j]]$hum.pop,
                            hum.pop = dat[[j]]$hum.pop,
                            log.hum = dat[[j]]$log.hum,
                            land.cover = dat[[j]]$land.cover,
                            s.effects = env.summary$s.effects[, k],
                            int = env.summary$int[, k],
                            sever = agr.ind$Severity/10,
                            beta0 = sb.summary$beta0[, k],
                            beta1 = sb.summary$beta1[, k],
                            rho = env.summary$rho[, k]
                     )
         }
         
         k = 1
         snakes.contr <- snake.env.contrib(S = dat[[j]][, 4:10],
                                                   H.bit = sb.incid[, 1] * dat[[j]]$hum.pop,
                                                   hum.pop = dat[[j]]$hum.pop,
                                                   log.hum = dat[[j]]$log.hum,
                                                   land.cover = dat[[j]]$land.cover,
                                                   s.effects = env.summary$s.effects[, k],
                                                   int = env.summary$int[, k],
                                                   sever = agr.ind$Severity/10,
                                                   beta0 = sb.summary$beta0[, k],
                                                   beta1 = sb.summary$beta1[, k])
         
      df <- data.frame(cbind(sb.incid,env.incid, snakes.contr))
      names(df) <- c("sb.med", "sb.025", "sb.97",
                     "env.med", "env.025", "env.97", "snakes")
      return(df)
   }
   saveRDS(projection, paste0("Incidence-projections/",
                              proj.names[i]))
}


## Envenoming snakes

