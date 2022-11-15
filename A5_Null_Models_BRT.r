#######################################################################################
################################### 1 - Null Models ###################################
#######################################################################################
### Esse script é usado para criar e avaliar 1.000 modelos de distribuição nulos, e 
### deve ser execultado após o script principal "A2_Model_Fitting". Os modelos nulos 
### são usados para testar a significância dos modelos criados com dados de espécies 
### reais

# setting random seed to always create the same random set
set.seed(1993)

# fitt and test one thousand null models
nullRep <- c(1:1000) # the number of repetitions is defined

# null models evaluate are stored in a list called 'nullEval'
nullEval <- list()

for(i in nullRep) {
# a set of 25% of randon select presence records are used to evaluate the model
presSamp <- sample(nrow(presValues), round(0.25 * nrow(presValues)))
testPres <- presValues[presSamp,]
# a set of 25% of randon select absence records are used to evaluate the model
abscSamp <- sample(nrow(abscValues), round(0.25 * nrow(abscValues)))
testAbsc <- abscValues[abscSamp,]

# a set of randon points are used like presence records to fit the model
# same sample size of specie presence records (75% of total)
presNull <- randomPoints(predictors, 21) 
trainPresNull <- extract(predictors, presNull)
# a set of randon points are used like absence records to fit the model
# same sample size of absence records (75% of total)
abscNull <- randomPoints(predictors, 150)
trainAbscNull <- extract(predictors, abscNull)

# differentiate null presences and null absence values 
trainNull <- c(rep(1, nrow(trainPresNull)), rep(0, nrow(trainAbscNull)))
# create a dataframe for train
train.dataNull <- data.frame(cbind(trainNull, rbind(trainPresNull, trainAbscNull)))

# differentiate presence and abscense values
test <- c(rep(1, nrow(testPres)), rep(0, nrow(testAbsc)))
# create a dataframe for test
test.data <- data.frame(cbind(test, rbind(testPres, testAbsc)))

#######################################################################################
# function gbm.fixed()is used (Ridgeway, 2006 and Elith et al. 2008)
null.tc2.lr002 <- gbm.fixed(data=train.dataNull,
gbm.x = 2:9,
gbm.y = 1,
tree.complexity = 2,
learning.rate = 0.002,
n.trees=1000)

#Predict values for test data
predsNull <- predict.gbm(null.tc2.lr002, test.data, n.trees=1000, "response")

# null models evaluate
d <- cbind(test.data$test, predsNull)
pres <- d[d[,1]==1, 2]
absc <- d[d[,1]==0, 2]
nullEval[[i]] <- evaluate(p=pres, a=absc)
}

# extract null AUC values
nullAuc <- sapply( nullEval, function(x){slot(x, 'auc')} )
# null AUC values are exported as .asc file
write.csv(x = nullAuc, file = "C:/Tfunalis/BRT/testNullAucValues.csv")

#######################################################################################
################################  END OF CODE  ########################################
#######################################################################################