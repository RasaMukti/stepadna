list.of.packages <- c("deSolve","ReacTran","RColorBrewer","geosphere","fields","stats4","GEOmap","geomapdata","plsgenomics","paramtest","bbmle","parallel","maps","mapdata","mapplots","rworldmap","scales","tidyverse","FME","optparse","readxl")
new.packages <-
  list.of.packages[!(list.of.packages %in% installed.packages()[, "Package"])]
if (length(new.packages)) {
  install.packages(new.packages,repos = "http://cran.us.r-project.org")
}

library(optparse)
option_list = list(
  make_option(c("-f", "--file"), type="character", default=NULL,
              help="dataset file name"),
  make_option(c("-o", "--out"), type="character", default="out.csv",
              help="output file name [default= %default]"),
  make_option(c("-a", "--age"), type="numeric", default=NULL,
              help="age of the allele"),
  make_option(c("-c", "--cores"), type="integer", default=1,
              help="number of cores to us in parallel"),
  make_option(c("-i", "--init"), type="integer", default=2,
              help="number of initial points for simulated annealing"),
  make_option(c("-l", "--loglik"), type="character", default="ph",
              help="whether to use pseudohaploid genotypes (ph) or genotype likelihoods (gl)"),
  make_option(c("--lon"), type="numeric", default=NULL,
              help="longitude of the origin of the allele. Has to be in the range from -10 to 80"),
  make_option(c("--lat"), type="numeric", default=NULL,
              help="latitude of the origin of the allele. Has to be in the range from 30 to 75")
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

library(deSolve)
library(ReacTran) # set up grid
library(RColorBrewer)
library(geosphere)
library(fields)
library(stats4)
library(GEOmap) # getETOPO
library(geomapdata) #for ETOPO data
library(plsgenomics)
library(paramtest)
library(bbmle)
library(parallel)
library(maps)
library(mapdata)
library(mapplots)
library(rworldmap)
library(scales)
library(tidyverse)
library(FME)
library("readxl")
source("./stepfunctions.R")
options(scipen = 999)

# Processing present LCT data
present_lct <- read.table("LCT_data/Lactase_all_GP.txt", header = TRUE)

present_lct <- present_lct[(present_lct["Long"] <= 80 &
                              present_lct["Long"] >= -10 &
                              present_lct["Lat"] <= 75 &
                              present_lct["Lat"] >= 30),]

present_lct <- present_lct[,c(6,5,4,3)]
colnames(present_lct) <- c("longitude", "latitude", "chrom_count", "frequency")

# Processing ancient LCT data
ancient_lct <- read_excel("LCT_data/ancient_data.xlsx")

ancient_lct <- data.frame(ancient_lct[,c(9,8,11,5)])
colnames(ancient_lct) <- c("longitude", "latitude", "pseudo_haplotype", "TimeBP")
ancient_lct <- na.omit(ancient_lct)

ancient_lct <- ancient_lct[(ancient_lct["longitude"] <= 80 &
                              ancient_lct["longitude"] >= -10 &
                              ancient_lct["latitude"] <= 75 &
                              ancient_lct["latitude"] >= 30),]
# Convert to generations
ancient_lct$TimeBP <- ancient_lct$TimeBP/29

# Round and scale times
ancient_lct$TimeBP= round(ancient_lct$TimeBP/10)

# Define map boundaries
la1 <- 30
la2 <- 75
lo1 <- -10
lo2 <- 80

# Create topographical data grid
topo  <- createTopo(la1, la2, lo1, lo2)
topoLon <- topo$topoLon
topoLat <- topo$topoLat
topo <- topo$topo

# Points with potential initial allele
lonseq <- seq(-5, 100, 10)
latseq <- seq(37, 70, 5)
Lonpoints <- rep(lonseq, length(latseq))
Latpoints <- rep(latseq, each=length(lonseq))
ptsremoved <- c(2,3,4,6,8,10,11,13,14,16,18,20,21,22,23,26,28,30,32,33,34,36,38,40,42,43,44,45,46,48,50,52,54,55,56,57,58,60,62,64,65,66,67,68,70,72,74,76,77)
Lonpoints <- Lonpoints[-ptsremoved]
Latpoints <- Latpoints[-ptsremoved]
Lonpoints[25] <- Lonpoints[25]+1
Lonpoints[14] <- Lonpoints[14]+1

# Initialise parameters
Nx <- nrow(topo)
Ny <- ncol(topo)
D = 2.5 # Initial population density
timeP = round(opt$age/290) 
timeSplit = round(5000/290)

x0_fixed <- Nx+1-which(topoLat == closest(topoLat, opt$lat))
y0_fixed <- which(topoLon == closest(topoLon, opt$lon))

if (is.null(opt$lat) == FALSE & is.null(opt$lon) == FALSE){
  print(c(opt$lat, opt$lon))
} else{print("is null")}

######################################## MLE ######################################## 

# Function for calculating likelihoods
LL_LCT <- function(iter, Nx=Nx, Ny=Ny, s=s, sigmax=sigmax, sigmay=sigmay, D=D, x0=x0, y0=y0, xtran=xtran, ytran=ytran, timeP=timeP, timeSplit=timeSplit, oldnew=NULL){
  
  # Set boundaries for SANN algorithm
  s.low <- 0
  s.high <- 0.1
  sigma.low <- 1
  sigma.high <- 100
  tran.low <- -2.5
  tran.high <- 2.5
  if (s < s.low | s > s.high | sigmax < sigma.low | sigmay < sigma.low |
      sigmax > sigma.high | sigmay > sigma.high |
      xtran < tran.low | ytran < tran.low | xtran > tran.high | ytran > tran.high
      | x0 < 1 | y0 < 1 | x0 > Nx | y0 > Ny | (topo[x0,y0] < 0) == TRUE){
    return(c(99999))
  } else{

    x0 <- round(x0)
    y0 <- round(y0)
    
    if (oldnew == 0){
      out <- FischerSolver(Nx, Ny, s, sigmax, sigmay, D, x0, y0, xtran, ytran, (timeP-timeSplit), oldnew = "old")
      out <- out[,2:ncol(out)]
      
      time_periods <- seq(timeP,timeSplit,-1)
      timevect <- timeP+1 - time_periods
    }

    
    if (oldnew == 1){
      out <- FischerSolver(Nx, Ny, s, sigmax, sigmay, D, x0, y0, xtran, ytran, timeSplit, oldnew="new", initMat=initMat)
      out <- out[,2:ncol(out)]
      
      time_periods <- seq(timeSplit,0,-1)
      timevect <- timeSplit+1 - time_periods
    }

    
    l <- 0

    # Optimizer crashes with high s values
    if (nrow(out) != length(time_periods)){
      print("ERRORRRRR")
      return(c(100000))
    }
    else{

      for (tp in 1:length(timevect)){

        mat_out <- matrix(unlist(out[timevect[tp],]),nrow=Nx,ncol=Ny)
        p <- t(mat_out[nrow(mat_out):1,])

        if (time_periods[tp] != 0){

          current_time <- ancient_lct[which(ancient_lct$TimeBP == time_periods[tp]),]
          # Pseudo haploid calculation
          if (nrow(current_time) != 0){
            
            for (i in 1:nrow(current_time)){

              # Extract spatial points matching observed point
              y_topo <- which(topoLon == closest(topoLon, current_time$longitude[i]))[1]
              x_topo <- which(topoLat == closest(topoLat, current_time$latitude[i]))[1]
  
              k <- (2-current_time$pseudo_haplotype[i])/2
              small_num <- 1e-5
              if (p[y_topo,x_topo] == 0 | round(p[y_topo,x_topo], 5) == 0){l <- l + log(small_num^k*(1-small_num)^(1-k))}
              else if (p[y_topo,x_topo] == 1 | round(p[y_topo,x_topo], 5) == 1){l <- l + log((1-small_num)^k*(small_num)^(1-k))}
              else{l <- l + log(p[y_topo,x_topo]^k*(1-p[y_topo,x_topo])^(1-k))}
            }
          }
        }
        # Likelihood calculation for present day
        if (time_periods[tp] == 0){
          for (i in 1:nrow(present_lct)){
            y_topo <- which(topoLon == closest(topoLon, present_lct$longitude[i]))[1]
            x_topo <- which(topoLat == closest(topoLat, present_lct$latitude[i]))[1]
            n <- present_lct$chrom_count[i]
            k <- round(present_lct$chrom_count[i]*present_lct$frequency[i])
            small_num <- 1e-5
            if (p[y_topo,x_topo] == 0 | round(p[y_topo,x_topo], 5) == 0){l <- l+bignchoosek(n,k)+k*log(small_num)+(n-k)*log(1-small_num)}
            else if (p[y_topo,x_topo] == 1 | round(p[y_topo,x_topo], 5) == 1){l <- l+bignchoosek(n,k)+k*log(1-small_num)+(n-k)*log(small_num)}
            else {l <- l+bignchoosek(n,k)+k*log(p[y_topo,x_topo])+(n-k)*log(1-p[y_topo,x_topo])}
            
          }
        }
      }

      print(round(c(s, sigmax, sigmay, xtran, ytran, -l), 3))

      if (is.nan(-l)){
        return(c(99999))
      } else{
        return(c(-l))}
    }
  }
}


#########################################################################################
if(opt$loglik == "ph"){
  func <- LL_LCT
} else if (opt$loglik == "gl"){
  func <- LL_GL
}

start_time <- Sys.time()

if (is.null(opt$lat) == FALSE & is.null(opt$lon) == FALSE){
  ResOld <- param_optim_fixedInit(init_points = opt$init, seed = 5, num_cores = opt$cores, LL=func, oldnew="old")
} else{ResOld <- param_optim_old(init_points = opt$init, seed = 5, num_cores = opt$cores, LL=func)}

CIres <- getCI(ResOld, oldnew="old", LL=func)
bestResOld <- CIres[1][[1]]
oldCIvect <- CIres[2:length(CIres)]

outOld <- FischerSolver(Nx, Ny, coef(bestResOld)[["s"]], coef(bestResOld)[["sigmax"]], coef(bestResOld)[["sigmay"]], D, coef(bestResOld)[["x0"]], coef(bestResOld)[["y0"]], coef(bestResOld)[["xtran"]], coef(bestResOld)[["ytran"]], (timeP-timeSplit), oldnew="old")
outOld <- outOld[,2:ncol(outOld)]
initMat <- matrix(unlist(outOld[nrow(outOld),]),nrow=Nx,ncol=Ny)
initMat[which(topo < 0, arr.ind = TRUE)] <- 0

if (is.null(opt$lat) == FALSE & is.null(opt$lon) == FALSE){
  ResNew <- param_optim_fixedInit(init_points = opt$init, seed = 5, num_cores = opt$cores, LL=func, oldnew="new")
} else{ResNew <- param_optim_new(init_points = opt$init, seed = 5, num_cores = opt$cores, LL=func)}

CIresNew <- getCI(ResNew, oldnew="new", LL=func)
bestResNEW <- CIresNew[1][[1]]
newCIvect <- CIresNew[2:length(CIresNew)]
end_time <- Sys.time()
print(c("Run time: ", (end_time-start_time)))

outNew <- FischerSolver(Nx, Ny, coef(bestResNEW)[["s"]], coef(bestResNEW)[["sigmax"]], coef(bestResNEW)[["sigmay"]], D, coef(bestResNEW)[["x0"]], coef(bestResNEW)[["y0"]], coef(bestResNEW)[["xtran"]], coef(bestResNEW)[["ytran"]], timeSplit, oldnew="new", initMat=initMat)
outNew <- outNew[,2:ncol(outNew)]

outDF <- data.frame(. = c(paste("before 5000"), "CI low", "CI high"), s = c(round(coef(bestResOld)[["s"]], 4), oldCIvect[[1]], oldCIvect[[2]]),
                    sigmax=c(round(coef(bestResOld)[["sigmax"]], 3), max(oldCIvect[[3]], 0), oldCIvect[[4]]),
                    sigmay=c(round(coef(bestResOld)[["sigmay"]], 3), max(oldCIvect[[5]], 0), oldCIvect[[6]]),
                    xtran=c(round(coef(bestResOld)[["xtran"]], 3), oldCIvect[[7]], oldCIvect[[8]]),
                    ytran=c(round(coef(bestResOld)[["ytran"]], 3), oldCIvect[[9]], oldCIvect[[10]]),
                    longitude=round(topoLon[coef(bestResOld)[["y0"]]]), latitude=round(topoLat[Nx+1-coef(bestResOld)[["x0"]]]),
                    alleleAge.years=timeP*290)
 
outDF <- rbind(outDF, data.frame(. = c(paste("after 5000"), "CI low", "CI high"), s = c(round(coef(bestResNEW)[["s"]], 4), newCIvect[[1]], newCIvect[[2]]),
                    sigmax=c(round(coef(bestResNEW)[["sigmax"]], 3), max(newCIvect[[3]], 0), newCIvect[[4]]),
                    sigmay=c(round(coef(bestResNEW)[["sigmay"]], 3), max(newCIvect[[5]], 0), newCIvect[[6]]),
                    xtran=c(round(coef(bestResNEW)[["xtran"]], 3), newCIvect[[7]], newCIvect[[8]]),
                    ytran=c(round(coef(bestResNEW)[["ytran"]], 3), newCIvect[[9]], newCIvect[[10]]),
                    longitude=round(topoLon[coef(bestResNEW)[["y0"]]]), latitude=round(topoLat[Nx+1-coef(bestResNEW)[["x0"]]]), alleleAge.years=timeP*290))

write.csv(outDF, opt$out, quote=FALSE)
save.image(paste0(opt$out,".RData"))


