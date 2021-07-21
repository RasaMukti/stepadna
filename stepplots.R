############################### Plot the results #################################
library(optparse)
option_list = list(
  make_option(c("-f", "--file"), type="character", default=NULL,
              help="dataset file name"),
  make_option(c("-o", "--out"), type="character", default="outplot.gif",
              help="output file name [default= %default]"),
  make_option(c("-n", "--nimage"), type="integer", default=20,
              help="number of images in GIF")
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

library(GEOmap)
library(geomapdata)
library(ReacTran)
library(geosphere)
library(deSolve)

library(magrittr)
library(maptools)
library(ggplot2)
library(cowplot)
library(animation)
source("./stepfunctions.R")

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

# Initialize some variables
readyFile <- read.csv(opt$file)

Nx <- nrow(topo)
Ny <- ncol(topo)
x0 <- Nx+1-which(round(topoLat) == readyFile[["latitude"]][1])[1]
y0 <- which(round(topoLon) == readyFile[["longitude"]][1])[1]
D = 2.5
timeP = round(readyFile[["alleleAge.years"]][1]/290) 
timeSplit = round(5000/290)

# Obtain the spatial distribution for time periods before 5000 years
outOld <- FischerSolver(Nx, Ny, readyFile[["s"]][1], readyFile[["sigmax"]][1], readyFile[["sigmay"]][1], D, x0, y0, readyFile[["xtran"]][1], readyFile[["ytran"]][1], (timeP-timeSplit), oldnew="old")
outOld <- outOld[,2:ncol(outOld)]
initMat <- matrix(unlist(outOld[nrow(outOld),]),nrow=Nx,ncol=Ny)
initMat[which(topo < 0, arr.ind = TRUE)] <- 0

# Obtain the spatial distribution for time periods after 5000 years
outNew <- FischerSolver(Nx, Ny, readyFile[["s"]][4], readyFile[["sigmax"]][4], readyFile[["sigmay"]][4], D, x0, y0, readyFile[["xtran"]][4], readyFile[["ytran"]][4], timeSplit, oldnew="new", initMat=initMat)
outNew <- outNew[,2:ncol(outNew)]

# Create the GIF and save it
saveGif(col_num = 20, image_num = opt$nimage, outName = opt$out, diffusion_mat = rbind(outOld, outNew), timeP = timeP, Nx = Nx, Ny = Ny, lo1, lo2, la1, la2, readyFile[["longitude"]][1], readyFile[["latitude"]][1])



