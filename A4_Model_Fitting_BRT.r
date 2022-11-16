#######################################################################################
###  Written by Sampaio, A. C. P. (1) & Cavalcante, A. de M. B. (2), 2022
###  1, 2 Instituto Nacional de Pesquisas Espaciais (INPE), Eusébio, Brazil
###  Published in:
#######################################################################################
##############################   MODEL FITTING ########################################
#######################################################################################
### This is the main script, used to fit a specie distribution model with Boosted
### Regression Trees (BRT) method and transfer it to future climate change scenarios.
### AUC is used to evaluate the models. For this the species data are partitioned into 
### training and test data
#######################################################################################
###################################  Index  ###########################################
### 1 - Load packages
### 2 - Load data (predictors layers, presences & abscenses)
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
library(gbm)
library(sp) # required for raster packages, coordinates() function
library(raster) # stack(), scale(), crop(), writeRaster() & raster() functions
library(dismo) # for evaluate() & gbm.step() functions

#######################################################################################
############## 2 - Load data(predictors layers, presences & abscenses) ################
#######################################################################################
##### os dados preparados no apêndice 1 devem estar contidos nos seus respectivos ##### 
################# diretórios, portanto, primeiro defina os diretórios #################

# Define the extent of a rectangular study area: lon 45W 34W; lat 2S 17S
ext = extent(-45, -34, -17, -2)

############################## Predictors Layers ######################################

### load predictors layers used in model trainning

# finds all the files with extension "asc" in the directory 
files <- list.files(path=paste('C:/Tinamoena/Layers', sep=''), 
pattern='asc', full.names=TRUE)
predictors <- stack(files) # create a raster stack
# Project stack if layers is not project
projection(predictors) <- CRS('+proj=longlat +datum=WGS84') # Project stack
predictors
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

################################## Presence Data ######################################

# this is the file wich presence records we will use:
file <- paste("C:/Tinamoena/T_inamoena.csv", sep="")
inamoena <- read.table(file, header=TRUE, sep=',')
# we do not need the first column
inamoena <- inamoena[,-1]
# extract values of the predictors at the presence points
presValues <- extract(predictors, inamoena)
# first layer of the RasterStack
plot(predictors, 1)
points(inamoena, col='red', pch='+')

################################# Absences Data #######################################

# Define a mask layer for absence sampling
abscMask <- raster("C:/Tinamoena/absence_mask_inamoena.asc")
projection(abscMask) <- CRS('+proj=longlat +datum=WGS84') # Project

# Create absence sample
set.seed(0) # setting random seed to always create the same random set
absence <- randomPoints(abscMask, 273) # select 273 random points
# extract values from predictors to absences points
abscValues <- extract(predictors, absence)
points(absence, col='blue', pch='+')

#######################################################################################
####################### 3 - Model fitting, evaluate & tranfers ########################
#######################################################################################

repetition <- c(1:20) # the number of repetitions is defined
# now we can fitt and test our model twenty times

# create an empty vector to store everything

# models evaluate are stored in a list called 'e'
e <- list() 

# Contributions of predictors are stored in a list called 'varContribution'
varCon <- list() 

# models tranfers are stored in a list called 'brt'
brt <- list() # current scenario
brt1 <- list() # SSP 245 for 2050 scenario
brt2 <- list() # SSP 585 for 2050 scenario
brt3 <- list() # SSP 245 for 2070 scenario
brt4 <- list() # SSP 585 for 2070 scenario

for(w in repetition) {

# a set of 75% of randon select presence records are used to fit the model
sampPres <- sample(nrow(presValues), round(0.75 * nrow(presValues)))
trainPres <- presValues[sampPres,]
# the others 25% is only used to evaluate the model
testPres <- presValues[-sampPres,]

# a set of 75% of randon select absence records are used to fit the model
sampAbsc <- sample(nrow(abscValues), round(0.75 * nrow(abscValues)))
trainAbsc <- abscValues[sampAbsc,]
# the others 25% is only used to evaluate the model
testAbsc <- abscValues[-sampAbsc,]

# differentiate presence and abscense values
train <- c(rep(1, nrow(trainPres)), rep(0, nrow(trainAbsc)))
# create a dataframe for train
train.data <- data.frame(cbind(train, rbind(trainPres, trainAbsc)))

# define 'soil' as categorical variable (called a 'factor' in R )
train.data[,'soil'] = as.factor(train.data[,'soil'])

# differentiate presence and abscense values
test <- c(rep(1, nrow(testPres)), rep(0, nrow(testAbsc)))
# create a dataframe for test
test.data <- data.frame(cbind(test, rbind(testPres, testAbsc)))

# define 'soil' as categorical variable (called a 'factor' in R )
test.data[,'soil'] = as.factor(test.data[,'soil'])

# function gbm.step()is used (Ridgeway, 2006 and Elith et al. 2008)
inamoena.tc2.lr002 <- gbm.step(data=train.data,
gbm.x = 2:9,
gbm.y = 1,
family = "bernoulli",
tree.complexity = 2,
learning.rate = 0.002,
bag.fraction = 0.5)

### Test data
### Predict values for test data
preds <- predict.gbm(inamoena.tc2.lr002, test.data, 
n.trees=inamoena.tc2.lr002$gbm.call$best.trees, type="response")

# Calculate deviance
calc.deviance(obs=test.data$test, pred=preds, calc.mean=TRUE)

# models evaluate
d <- cbind(test.data$test, preds)
pres <- d[d[,1]==1, 2]
absc <- d[d[,1]==0, 2]
e[[w]] <- evaluate(p=pres, a=absc)

varCon[[w]] <- inamoena.tc2.lr002$contributions

# models tranfers
brt[[w]] <- predict(predictors, inamoena.tc2.lr002, 
n.trees=inamoena.tc2.lr002$gbm.call$best.trees, type="response")# current 
brt1[[w]] <- predict(transfer1, inamoena.tc2.lr002, 
n.trees=inamoena.tc2.lr002$gbm.call$best.trees, type="response")# SSP 245 for 2050 
brt2[[w]] <- predict(transfer2, inamoena.tc2.lr002, 
n.trees=inamoena.tc2.lr002$gbm.call$best.trees, type="response")# SSP 585 for 2050 
brt3[[w]] <- predict(transfer3, inamoena.tc2.lr002, 
n.trees=inamoena.tc2.lr002$gbm.call$best.trees, type="response")# SSP 245 for 2070 
brt4[[w]] <- predict(transfer4, inamoena.tc2.lr002, 
n.trees=inamoena.tc2.lr002$gbm.call$best.trees, type="response")# SSP 585 for 2070 
}

# extract AUC values
auc <- sapply( e, function(x){slot(x, 'auc')} )

# extract "Maximum of the sum of the sensitivity and specificity" threshold
mst <- sapply( e, function(x){ x@t[which.max(x@TPR + x@TNR)] } )
threshold <- mean(mst)

# create a dataframe with variables contributions
varContributions <- data.frame(varCon[[1]], varCon[[2]], varCon[[3]], varCon[[4]], 
varCon[[5]], varCon[[6]], varCon[[7]], varCon[[8]], varCon[[9]], varCon[[10]], 
varCon[[11]], varCon[[12]], varCon[[13]], varCon[[14]], varCon[[15]], varCon[[16]], 
varCon[[17]], varCon[[18]], varCon[[19]], varCon[[20]])

# AUC values are exported as .asc file
write.csv(x = auc, file = "C:/Tinamoena/BRT/testAucValues.csv")

# threshold values are exported as .asc file
write.csv(x = mst, file = "C:/Tinamoena/BRT/msthresholdValues.csv")

# Variables contributions values are exported as .asc file
write.csv(x = varContributions, file = "C:/Tinamoena/BRT/variables_contributions.csv")

#######################################################################################
########################## 4 - Combining models predictions ###########################
#######################################################################################

# create a raster stack (predictions for current scenario)
models <- stack(brt[[1]], brt[[2]],brt[[3]], brt[[4]], brt[[5]], brt[[6]], 
brt[[7]],brt[[8]], brt[[9]], brt[[10]], brt[[11]], brt[[12]],brt[[13]], 
brt[[14]], brt[[15]], brt[[16]], brt[[17]],brt[[18]], brt[[19]], brt[[20]])

# create a raster stack (predictions for SSP 245 for 2050 scenario)
models1 <- stack(brt1[[1]], brt1[[2]],brt1[[3]], brt1[[4]], brt1[[5]], brt1[[6]], 
brt1[[7]],brt1[[8]], brt1[[9]], brt1[[10]], brt1[[11]], brt1[[12]],brt1[[13]], 
brt1[[14]], brt1[[15]], brt1[[16]], brt1[[17]],brt1[[18]], brt1[[19]], brt1[[20]])

# create a raster stack (predictions for SSP 585 for 2050 scenario) 
models2 <- stack(brt2[[1]], brt2[[2]],brt2[[3]], brt2[[4]], brt2[[5]], brt2[[6]], 
brt2[[7]],brt2[[8]], brt2[[9]], brt2[[10]], brt2[[11]], brt2[[12]],brt2[[13]], 
brt2[[14]], brt2[[15]], brt2[[16]], brt2[[17]],brt2[[18]], brt2[[19]], brt2[[20]])

# create a raster stack (predictions for SSP 245 for 2070 scenario)
models3 <- stack(brt3[[1]], brt3[[2]],brt3[[3]], brt3[[4]], brt3[[5]], brt3[[6]], 
brt3[[7]],brt3[[8]], brt3[[9]], brt3[[10]], brt3[[11]], brt3[[12]],brt3[[13]], 
brt3[[14]], brt3[[15]], brt3[[16]], brt3[[17]],brt3[[18]], brt3[[19]], brt3[[20]])

# create a raster stack (predictions for SSP 585 for 2070 scenario)
models4 <- stack(brt4[[1]], brt4[[2]],brt4[[3]], brt4[[4]], brt4[[5]], brt4[[6]], 
brt4[[7]],brt4[[8]], brt4[[9]], brt4[[10]], brt4[[11]], brt4[[12]],brt4[[13]], 
brt4[[14]], brt4[[15]], brt4[[16]], brt4[[17]],brt4[[18]], brt4[[19]], brt4[[20]])

# compute the simple average for predictions:
m <- mean(models) # current scenario
m1 <- mean(models1) # SSP 245 for 2050 scenario
m2 <- mean(models2) # SSP 585 for 2050 scenario
m3 <- mean(models3) # SSP 245 for 2070 scenario
m4 <- mean(models4) # SSP 285 for 2070 scenario

# plot models results 
par(mfrow=c(2,3))
plot(m, main='Current')
plot(m1, main='2050 SSP245')
plot(m2, main='2050 SSP585')
plot(m3, main='2070 SSP245')
plot(m4, main='2070 SSP585')

# create a binary map of presence and abscence
presAbsc <- (m > threshold) # current scenario
presAbsc1 <- (m1 > threshold) # SSP 245 for 2050 scenario
presAbsc2 <- (m2 > threshold) # SSP 585 for 2050 scenario
presAbsc3 <- (m3 > threshold) # SSP 245 for 2070 scenario
presAbsc4 <- (m4 > threshold) # SSP 285 for 2070 scenario

#######################################################################################
################################ 5 - Export results ###################################
#######################################################################################

# suitability map (current scenario)			
writeRaster(m,
            filename  = "C:/Tinamoena/BRT/current_BRT.asc",
            format    = 'ascii',
            NAflag    = -9999,
            overwrite = TRUE)	
			
writeRaster(presAbsc,
            filename  = "C:/Tinamoena/BRT/current_BRT_pa.asc",
            format    = 'ascii',
            NAflag    = -9999,
            overwrite = TRUE)	
			
# suitability map (SSP 245 for 2050 scenario)
writeRaster(m1,
            filename  = "C:/Tinamoena/BRT/2050_SSP245_BRT.asc",
            format    = 'ascii',
            NAflag    = -9999,
            overwrite = TRUE)	
			
writeRaster(presAbsc1,
            filename  = "C:/Tinamoena/BRT/2050_SSP245_BRT_pa.asc",
            format    = 'ascii',
            NAflag    = -9999,
            overwrite = TRUE)	

# suitability map (SSP 585 for 2050 scenario)
writeRaster(m2,
            filename  = "C:/Tinamoena/BRT/2050_SSP585_BRT.asc",
            format    = 'ascii',
            NAflag    = -9999,
            overwrite = TRUE)
			
writeRaster(presAbsc2,
            filename  = "C:/Tinamoena/BRT/2050_SSP585_BRT_pa.asc",
            format    = 'ascii',
            NAflag    = -9999,
            overwrite = TRUE)	


# suitability map (SSP 245 for 2070 scenario) 			
writeRaster(m3,
            filename  = "C:/Tinamoena/BRT/2070_SSP245_BRT.asc",
            format    = 'ascii',
            NAflag    = -9999,
            overwrite = TRUE)	
			
writeRaster(presAbsc3,
            filename  = "C:/Tinamoena/BRT/2070_SSP245_BRT_pa.asc",
            format    = 'ascii',
            NAflag    = -9999,
            overwrite = TRUE)	


# suitability map (SSP 585 for 2070 scenario)	
writeRaster(m4,
            filename  = "C:/Tinamoena/BRT/2070_SSP585_BRT.asc",
            format    = 'ascii',
            NAflag    = -9999,
            overwrite = TRUE)	
			
writeRaster(presAbsc4,
            filename  = "C:/Tinamoena/BRT/2070_SSP585_BRT_pa.asc",
            format    = 'ascii',
            NAflag    = -9999,
            overwrite = TRUE)	
			
