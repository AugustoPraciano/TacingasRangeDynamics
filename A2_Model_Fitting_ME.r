#######################################################################################
###  Written by Sampaio, A. C. P. (1) & Cavalcante, A. de M. B. (2), 2022
###  1, 2 Instituto Nacional de Pesquisas Espaciais (INPE), Eusébio, Brazil
###  Published in:
#######################################################################################
##############################   MODEL FITTING ########################################
#######################################################################################
### This is the main script used to fit a Maxent distribution model and transfer it to 
### future climate change scenarios. AUC is used to evaluate the models. For this
### the species data are partitioned into training and test data
#######################################################################################
###################################  Index  ###########################################
### 1 - Load packages
### 2 - Load data (predictors layers, presences & background)
### 3 - Model fitting, evaluate & tranfers
### 4 - Combining models predictions
### 5 - Export results
#######################################################################################

#######################################################################################
################################## 1 - Load packages ##################################
#######################################################################################

rm(list = ls(all = TRUE))
# Set a directory
setwd ("C:/Tinamoena")
getwd()

#Load Packages
library(raster) # stack(), scale(), crop(), writeRaster() & raster() functions
library(sp) # required for raster packages, coordinates() function
library(dismo) # for evaluate() & maxent() functions
library(rJava) # required for dismo
library(maptools) # readShapeSpatial() function
library(rgeos) # required for maptools
library(Hmisc) # rcorr() function
library(ade4) # dudi.pca() function
library(factoextra) # get_eigenvalue() function
library(effsize) #cohen.d() function
data(wrld_simpl) # countries boundaries

#######################################################################################
############## 2 - Load data(predictors layers, presences & background) ###############
#######################################################################################
##### Data that was prepared in "A1 Data Preparation" must be contained in their ###### 
##################### directories, so first define the directories ####################

# Define the extent of a rectangular study area: lon 45W 34W; lat 2S 17S
ext = extent(-45, -34, -17, -2)

############################## Predictors Layers ######################################

### load predictors layers used in model trainning

# finds all the files with extension "asc" in the directory 
files <- list.files(path=paste('C:/Tinamoena/Layers', sep=''), 
pattern='asc', full.names=TRUE)
predictors <- stack(files) # create a raster stack
projection(predictors) <- CRS('+proj=longlat +datum=WGS84') # Project stack
names(predictors)
plot(predictors)

### load predictors layers used in model tranfers 

### SSP 245 for 2050 scenario
# finds all the files with extension "asc" in the directory 
files1 <- list.files(path=paste('C:/Tinamoena/Layers_1', sep=''), 
pattern='asc', full.names=TRUE)
transfer1 <- stack(files1) # create a raster stack
projection(transfer1) <- CRS('+proj=longlat +datum=WGS84') # Project stack

### SSP 585 for 2050 scenario
# finds all the files with extension "asc" in the directory
files2 <- list.files(path=paste('C:/Tinamoena/Layers_2', sep=''), 
pattern='asc', full.names=TRUE)
transfer2 <- stack(files2) # create a raster stack
projection(transfer2) <- CRS('+proj=longlat +datum=WGS84') # Project stack

### SSP 245 for 2070 scenario
# finds all the files with extension "asc" in the directory
files3 <- list.files(path=paste('C:/Tinamoena/Layers_3', sep=''), 
pattern='asc', full.names=TRUE)
transfer3 <- stack(files3) # create a raster stack
projection(transfer3) <- CRS('+proj=longlat +datum=WGS84') # Project stack

### SSP 585 for 2070 scenario
# finds all the files with extension "asc" in the directory
files4 <- list.files(path=paste('C:/Tinamoena/Layers_4', sep=''), 
pattern='asc', full.names=TRUE)
transfer4 <- stack(files4) # create a raster stack
projection(transfer4) <- CRS('+proj=longlat +datum=WGS84') # Project stack

#################################### Presence Data ####################################

# this is the file wich presence records we will use:
file <- paste("C:/Tinamoena/T_inamoena.csv", sep="")
inamoena <- read.table(file, header=TRUE, sep=',')
# we do not need the first column
inamoena <- inamoena[,-1]
# extract values of the predictors at the presence points
presValues <- extract(predictors, inamoena)
# plot first layer of the RasterStack
plot(predictors, 1)
points(inamoena, col='red', pch='+')

################################# Background Data #####################################

set.seed(0) # setting random seed to always create the same random set
backgr <- randomPoints(predictors, 10000) # select 10,000 random points
# extract values from predictors to background points
backValues <- extract(predictors, backgr)

#######################################################################################
####################### 3 - Model fitting, evaluate & tranfers ########################
#######################################################################################

repetition <- c(1:20) # the number of repetitions is defined
# now we can fitt and test our model twenty times

# create an empty vector to store everything

e <- list() # models evaluate are stored in a list called 'e'

maxResults <- list() # Maxent results are stored in a list called 'maxResults'

# models tranfers are stored in a list called 'pmx'
pmx <- list() # current scenario
pmx1 <- list() # SSP 245 for 2050 scenario
pmx2 <- list() # SSP 585 for 2050 scenario
pmx3 <- list() # SSP 245 for 2070 scenario
pmx4 <- list() # SSP 585 for 2070 scenario


for(w in repetition) {
# a set of 75% of randon select presence records are used to fit the model
samp <- sample(nrow(presValues), round(0.75 * nrow(presValues)))
trainPres <- presValues[samp,]
# the others 25% is only used to evaluate the model
testPres <- presValues[-samp,]

trainData <- data.frame(rbind(trainPres, backValues)) # create a dataframe for train
# differentiate presence and background values
pb <- c(rep(1, nrow(trainPres)), rep(0, nrow(backValues)))

# define 'soil' as categorical variable (called a 'factor' in R)
trainData[,'soil'] = as.factor(trainData[,'soil'])

# differentiate presence and abscense (background) values
pa <- c(rep(1, nrow(testPres)), rep(0, nrow(backValues)))
# create a dataframe for test
testData <- data.frame(cbind(pa, rbind(testPres, backValues)))

# define 'soil' as categorical variable (called a 'factor' in R)
testData[,'soil'] = as.factor(testData[,'soil'])

### Fitting a model
# function maxent() is used(Hijmans & al., 2013)
mx <- maxent(trainData, 
             pb,
			 path = "C:/Tinamoena/Maxent/output",
			 args = c("redoifexists", "notooltips", "noautofeature", "linear", 
			 "quadratic", "hinge", "product", "threshold"))

# models evaluate
e[[w]] <- evaluate(testData[testData==1,], testData[testData==0,], mx)
# Maxent results are exported as .asc file
maxResults[[w]] <- read.csv("C:/Tinamoena/Maxent/output/maxentResults.csv") 
# models tranfers
pmx[[w]] <- predict(predictors, mx, ext=ext, progress='') # current scenario
pmx1[[w]] <- predict(transfer1, mx, ext=ext, progress='') # SSP 245 for 2050 scenario
pmx2[[w]] <- predict(transfer2, mx, ext=ext, progress='') # SSP 585 for 2050 scenario
pmx3[[w]] <- predict(transfer3, mx, ext=ext, progress='') # SSP 245 for 2070 scenario 
pmx4[[w]] <- predict(transfer4, mx, ext=ext, progress='') # SSP 585 for 2070 scenario
 }

# extract AUC values
auc <- sapply( e, function(x){slot(x, 'auc')} )

# extract "Maximum of the sum of the sensitivity and specificity" threshold
mst <- sapply( e, function(x){ x@t[which.max(x@TPR + x@TNR)] } )
threshold <- mean(mst)

# create a dataframe with Maxent results
maxentResults <- data.frame(rbind(maxResults[[1]], maxResults[[2]], maxResults[[3]], 
maxResults[[4]], maxResults[[5]], maxResults[[6]], maxResults[[7]], maxResults[[8]], 
maxResults[[9]], maxResults[[10]], maxResults[[11]], maxResults[[12]], 
maxResults[[13]], maxResults[[14]], maxResults[[15]], maxResults[[16]], 
maxResults[[17]], maxResults[[18]], maxResults[[19]], maxResults[[20]]))

# AUC values are exported as .asc file
write.csv(x = auc, file = "C:/Tinamoena/Maxent/testAucValues.csv")

# threshold values are exported as .asc file
write.csv(x = mst, file = "C:/Tinamoena/Maxent/msthresholdValues.csv")

# Maxent results are exported as .asc file
write.csv(x = maxentResults, file = "C:/Tinamoena/Maxent/maxentResults.csv")

#######################################################################################
########################## 4 - Combining models predictions ###########################
#######################################################################################

# create a raster stack (predictions for current scenario)
models <- stack(pmx[[1]], pmx[[2]],pmx[[3]], pmx[[4]], pmx[[5]], pmx[[6]], 
pmx[[7]],pmx[[8]], pmx[[9]], pmx[[10]], pmx[[11]], pmx[[12]],pmx[[13]], 
pmx[[14]], pmx[[15]], pmx[[16]], pmx[[17]],pmx[[18]], pmx[[19]], pmx[[20]])

# create a raster stack (predictions for SSP 245 for 2050 scenario)
models1 <- stack(pmx1[[1]], pmx1[[2]],pmx1[[3]], pmx1[[4]], pmx1[[5]], pmx1[[6]], 
pmx1[[7]],pmx1[[8]], pmx1[[9]], pmx1[[10]], pmx1[[11]], pmx1[[12]],pmx1[[13]], 
pmx1[[14]], pmx1[[15]], pmx1[[16]], pmx1[[17]],pmx1[[18]], pmx1[[19]], pmx1[[20]])

# create a raster stack (predictions for SSP 585 for 2050 scenario) 
models2 <- stack(pmx2[[1]], pmx2[[2]],pmx2[[3]], pmx2[[4]], pmx2[[5]], pmx2[[6]], 
pmx2[[7]],pmx2[[8]], pmx2[[9]], pmx2[[10]], pmx2[[11]], pmx2[[12]],pmx2[[13]], 
pmx2[[14]], pmx2[[15]], pmx2[[16]], pmx2[[17]],pmx2[[18]], pmx2[[19]], pmx2[[20]])

# create a raster stack (predictions for SSP 245 for 2070 scenario)
models3 <- stack(pmx3[[1]], pmx3[[2]],pmx3[[3]], pmx3[[4]], pmx3[[5]], pmx3[[6]], 
pmx3[[7]],pmx3[[8]], pmx3[[9]], pmx3[[10]], pmx3[[11]], pmx3[[12]],pmx3[[13]], 
pmx3[[14]], pmx3[[15]], pmx3[[16]], pmx3[[17]],pmx3[[18]], pmx3[[19]], pmx3[[20]])

# create a raster stack (predictions for SSP 585 for 2070 scenario)
models4 <- stack(pmx4[[1]], pmx4[[2]],pmx4[[3]], pmx4[[4]], pmx4[[5]], pmx4[[6]], 
pmx4[[7]],pmx4[[8]], pmx4[[9]], pmx4[[10]], pmx4[[11]], pmx4[[12]],pmx4[[13]], 
pmx4[[14]], pmx4[[15]], pmx4[[16]], pmx4[[17]],pmx4[[18]], pmx4[[19]], pmx4[[20]])

# compute the simple average for predictions:
m <- mean(models) # current scenario
m1 <- mean(models1) # SSP 245 for 2050 scenario
m2 <- mean(models2) # SSP 585 for 2050 scenario
m3 <- mean(models3) # SSP 245 for 2070 scenario
m4 <- mean(models4) # SSP 585 for 2070 scenario

# plot models results 
par(mfrow=c(2,3))
plot(m, main='Current')
plot(m1, main='2050 SSP245')
plot(m2, main='2050 SSP585')
plot(m3, main='2070 SSP245')
plot(m4, main='2070 SSP585')

# create a binary map of presence and abscence
presAbsc <- (m > threshold) # current scenario
presAbsc1 <- (m1 > threshold)  # SSP 245 for 2050 scenario
presAbsc2 <- (m2 > threshold) # SSP 585 for 2050 scenario
presAbsc3 <- (m3 > threshold) # SSP 245 for 2070 scenario
presAbsc4 <- (m4 > threshold) # SSP 585 for 2070 scenario

#######################################################################################
################################ 5 - Export results ###################################
#######################################################################################

# suitability map (current scenario)
writeRaster(m,
            filename  = "C:/Tinamoena/Maxent/current_ME.asc",
            format    = 'ascii',
            NAflag    = -9999,
            overwrite = TRUE)
			
writeRaster(presAbsc,
            filename  = "C:/Tinamoena/Maxent/current_ME_pa.asc",
            format    = 'ascii',
            NAflag    = -9999,
            overwrite = TRUE)			

# suitability map (SSP 245 for 2050 scenario)
writeRaster(m1,
            filename  = "C:/Tinamoena/Maxent/2050_SSP245_ME.asc",
            format    = 'ascii',
            NAflag    = -9999,
            overwrite = TRUE)
			
writeRaster(presAbsc1,
            filename  = "C:/Tinamoena/Maxent/2050_SSP245_ME_pa.asc",
            format    = 'ascii',
            NAflag    = -9999,
            overwrite = TRUE)	

# suitability map (SSP 585 for 2050 scenario)
writeRaster(m2,
            filename  = "C:/Tinamoena/Maxent/2050_SSP585_ME.asc",
            format    = 'ascii',
            NAflag    = -9999,
            overwrite = TRUE)		
			
writeRaster(presAbsc2,
            filename  = "C:/Tinamoena/Maxent/2050_SSP585_ME_pa.asc",
            format    = 'ascii',
            NAflag    = -9999,
            overwrite = TRUE)	

# suitability map (SSP 245 for 2070 scenario) 	
writeRaster(m3,
            filename  = "C:/Tinamoena/Maxent/2070_SSP245_ME.asc",
            format    = 'ascii',
            NAflag    = -9999,
            overwrite = TRUE)	
			
writeRaster(presAbsc3,
            filename  = "C:/Tinamoena/Maxent/2070_SSP245_ME_pa.asc",
            format    = 'ascii',
            NAflag    = -9999,
            overwrite = TRUE)	

# suitability map (SSP 585 for 2070 scenario)
writeRaster(m4,
            filename  = "C:/Tinamoena/Maxent/2070_SSP585_ME.asc",
            format    = 'ascii',
            NAflag    = -9999,
            overwrite = TRUE)	
			
writeRaster(presAbsc4,
            filename  = "C:/Tinamoena/Maxent/2070_SSP585_ME_pa.asc",
            format    = 'ascii',
            NAflag    = -9999,
            overwrite = TRUE)			

#######################################################################################
################################  END OF CODE  ########################################
#######################################################################################