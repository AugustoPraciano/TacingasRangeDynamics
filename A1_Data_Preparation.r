#######################################################################################
###  Written by Sampaio, A. C. P. (1) & Cavalcante, A. de M. B. (2), 2022
###  1, 2 Instituto Nacional de Pesquisas Espaciais (INPE), Eusébio, Brazil
###  Published in: 
#######################################################################################
##############################  DATA PREPARATION  #####################################
#######################################################################################
### This script is used to prepare and select of species' presence data and 
### environmental variables to be used in the modeling step.
### Applied to R software version 4.2.1
#######################################################################################
###################################  Index  ###########################################
### 1 - Load packages
### 2 - Presence records preparation (check for duplicates & reduce autocorrelation)
### 3 - Analysis of multicollinearity & selection of predictors variables
#######################################################################################

#######################################################################################
################################## 1 - Load packages ##################################
#######################################################################################

rm(list = ls(all = TRUE))
# Set a directory
setwd ("C:/Tinamoena")
getwd()

library(raster) # stack(), scale(), crop(), crs(), writeRaster() & raster() functions
library(sp) # required for raster packages, coordinates() function
library(dismo) # for evaluate() gridsample() & maxent() functions
library(rJava) # required for dismo
library(maptools) # readShapeSpatial() function
library(rgeos) # required for maptools
library(Hmisc) # rcorr() function
library(ade4) # dudi.pca() function
library(factoextra) # get_eigenvalue() function
library(effsize) #cohen.d() function
data(wrld_simpl) # countries boundaries

#######################################################################################
### 2 - Presence records preparation (check for duplicates & reduce autocorrelation)###
#######################################################################################

# this is the species' sample file we will use:
file <- paste("C:/Tinamoena/tacinga_inamoena.csv", sep="")
# read it
inamoena <- read.table(file, header=TRUE, sep=",")
# inspect the values of the file (first rows)
head(inamoena)
# we only need columns 2 and 3:
inamoena <- inamoena[,2:3]
head(inamoena)

# plot species records map
plot(wrld_simpl, xlim=c(-80,-30), ylim=c(-40,10), axes=TRUE, col="light grey")
box() # restore the box around the map
points(inamoena$dd_long, inamoena$dd_lat, col='red', pch=20, cex=0.75)

# check for duplicate records
dups2 <- duplicated(inamoena[, c('dd_long', 'dd_lat')])
sum(dups2) # number of duplicates
# keep the records that are not duplicated
ina <- inamoena[!dups2, ]

# use the coordinates function to create a SpatialPointsDataFrame
coordinates(ina) <- ~dd_long+dd_lat
crs(ina) <- crs(wrld_simpl)
class(ina)
# create a SpatialPolygonsDataFrame
class(wrld_simpl)

# spatially rarefy occurence data (reduce spatial autocorrelation)
# create a RasterLayer with the extent of T.inamoena records
r <- raster(ina)
#set the resolution of the cells to ~ 1 minute (value in degrees)
res(r) <- 0.016
# expand (extend) the extent of the RasterLayer a little
r <- extend(r, extent(r)+1)
# resample:
inaSelected <- gridSample(ina, r, n=1)
# to illustrate the method and show the result
plot(ina)
points(inaSelected, cex=1, col='red', pch='x') # selected points in red
# save the rarefied data set
write.csv(x = inaSelected, file = "C:/Tinamoena/T_inamoena.csv")

#######################################################################################
######## 3 - Analysis of multicollinearity & selection of predictors variables ########
#######################################################################################

# Define the extent of the study area: lon 45W 34W; lat 2S 17S
ext <- extent(-45, -34, -17, -2)

# load all environmental variables (predictors)

# finds all the files with extension "asc" in the directory
files <- list.files(path=paste('C:/Tinamoena/Layers', 
sep=''), pattern='asc', full.names=TRUE )
filesStack <- stack(files) # create a raster stack
names(filesStack) # get the predictor names
# plot a single layer in a RasterStack, and plot some additional data on top of it
plot(filesStack, 1) # first layer of the RasterStack
points(inaSelected, col='blue', pch=20)

filesVector <- na.omit(values(filesStack)) # As vector (NA's are omitted)
filesDf <- data.frame(filesVector) # put in a dataframe
# convert to matrix, required by fucntion rcorr()
filesMatrix <- as.matrix(filesDf)
colnames(filesMatrix)

# calculate the Spearman rank correlation between all predictors variables
setwd ("C:/Tinamoena/Layers")
predRcorr <- rcorr(filesMatrix, type = "spearman")
# export the Spearman r values and significance levels
write.table(predRcorr$r,
         file = "predSpearmanR.txt",
         sep = ",",
         quote = FALSE,
         append = FALSE,
         na = "NA",
         qmethod = "escape")
write.table(predRcorr$P,
         file = "predSpearmanSignificance.txt",
         sep = ",",
         quote = FALSE,
         append = FALSE,
         na = "NA",
         qmethod = "escape")

# outside R, identify correlated predictors based on Spearman rank > 0.7

# calculate PCA for predictors variables
# standardize data, required for PCA analysis
filesPredScale <- scale(filesMatrix, center = TRUE, scale = TRUE)
# calculate PCA
pcaPred <- dudi.pca(filesPredScale, center = TRUE, scale = TRUE, scannf = FALSE, nf = 4)
pcaPred$co # the column coordinates, which is the same as vector load for each PC
write.csv(pcaPred$co, 
file = "C:/Tinamoena/Layers/pcaPredVariables.csv")

eigVal <- get_eigenvalue(pcaPred)# variance of each PCA axis
head(eigVal)
summary(pcaPred) #summary of results

# plot a graph of variance
barplot(eigVal[, 2], names.arg=1:nrow(eigVal), 
       main = "Variances",
       xlab = "Principal Components",
       ylab = "Percentage of variances",
       col ="steelblue")
# Add connected line segments to the plot
lines(x = 1:nrow(eigVal), eigVal[, 2], 
      type="b", pch=19, col = "red")

# outside R, keep only non correlated variables (Spearman rank > 0.7) and
# with larger variance contribution for PCA axis in Layers folder

#######################################################################################
################################  END OF CODE  ########################################
#######################################################################################