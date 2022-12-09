###################################################################################
# SET UP
###################################################################################
library("rmarkdown")
#library("devtools")
library("githubinstall")
library("dplyr")
library("tidyverse")
library("ggplot2")
library("reshape")
library("ggprism")
library("gridExtra")
library("stringr")
library("lemon")
library("DiagrammeR")
library("webshot")
# webshot::install_phantomjs(force = TRUE)
# githubinstall("plspm")
library("plspm")
# install.packages("devtools")
library("devtools")
# install_github("gastonstat/plspm")
library("plspm")


"%!in%" <- function(x,y)!('%in%'(x,y))



###################################################################################
# IMPORT DATASET
###################################################################################
pressure <- read.csv("~/Desktop/Life/Personal/College/Masters/Johns Hopkins/Statistical Theory I/Stat I Stat GitHub Repo/NAHNES_pressure.txt", 
                     header = TRUE)

# SEQN = Respondent sequence number
# ALQ101 = Had at least 12 alcohol drinks/1 yr?             - #1
# ALQ110 = Had at least 12 alcohol drinks/lifetime?
# ALQ130 = Avg # alcoholic drinks/day - past 12 mos
# SMQ020 = Smoked at least 100 cigarettes in life
# RIAGENDR = Gender                                         - #2
# HIQ210 = Time when no insurance in past year?             - #6


# Remove respondent number
pressure1 <- pressure[, names(pressure) %!in% c("SEQN")]

# Removed census data                                       - #3
pressure1 <- pressure1[, names(pressure1) %!in% c("RIDRETH1", "DMDCITZN", "DMDEDUC2","DMDMARTL", "WTINT2YR", "SDMVPSU", "SDMVSTRA")]

# Removed categorical data
pressure1 <- pressure1[, names(pressure1) %!in% c("ALQ101", "ALQ110", "SMQ020", "HIQ210")]

# RIDAGEYR = Age in years at screening                      - #2

# BPXSY1 = Systolic: Blood pressure (1st reading) mmHg      - #4
# BPXDI1 = Diastolic: Blood pressure (1st reading) mmHg
# BPXSY2 = Systolic: Blood pressure (2nd reading) mmHg
# BPXDI2 = Diastolic: Blood pressure (2nd reading) mmHg

# BMXWT = Weight (kg)                                       - #5
# BMXHT = Standing Height (cm)
# BMXBMI = Body Mass Index (kg/m**2)
# BMXLEG = Upper Leg Length (cm)
# BMXARML = Upper Arm Length (cm)
# BMXARMC = Arm Circumference (cm)
# BMXWAIST = Waist Circumference (cm)



#Should alq130 NA values be replaced with zero?
#pressure1$ALQ130[is.na(pressure1$ALQ130)] <- 0
pressure_data <- na.omit(pressure1)


###################################################################################
# SELECTION FUNCTION
###################################################################################
selection <- function(model_R, model_Q, saturated, n, predictors){
  r2_R = which(str_detect(rownames(model_R$inner_summary), "Hypertension"))
  r2_Q = which(str_detect(rownames(model_Q$inner_summary), "Hypertension"))
  r2_sat = which(str_detect(rownames(saturated$inner_summary), "Hypertension"))
  
  pk = predictors + 1
  SSE = (1 - model_R$inner_summary$R2[r2_R])*(n - 1)
  MSE = 1 - saturated$inner_summary$R2[r2_sat]
  
  ## PLS selection criteria
  Rsq       = model_R$inner_summary$R2[r2_R]
  Qsq       = model_R$inner_summary$R2[r2_Q]
  AdjRsq    = 1 - (  ( (n - 1)/(n - pk) )*(1 - model_R$inner_summary$R2[r2_R]) )
  GoF       = model_R$gof
  
  ## Distance based criteria
  FPE       = (SSE/(n - pk))*(1 + (pk/n))
  Mallows   = (SSE/MSE) - (n - pk)               # number of predictors
  AIC       = n*(log(SSE/n) + ((2*pk)/n))
  AICu      = n*(log(SSE/(n - pk)) + ((2*pk)/n))
  AICc      = n*(log(SSE/(n)) + ((n + pk)/(n - pk - 2)))
  
  ## Consistent criteria
  BIC       = n*(log(SSE/n) + ((pk*log(n))/n))
  GM        = (SSE/MSE) + (pk*log(n))
  HQ        = n*(log(SSE/n) + ((2*pk*log(log(n)) )/n) )
  HQc       = n*(log(SSE/n) + ((2*pk*log(log(n)) )/(n - pk - 2)) )
  
  
  ## Compile
  result <- data.frame("Type" = c(rep("PLS selection", 4), 
                                  rep("Distance based", 5), 
                                  rep("Consistent criteria", 4)),
                       "Criteria" = c("Rsq", "Qsq", "AdjRsq", "GoF",
                                      "FPE", "Mallow's Cp", "AIC", "AICu", "AICc",
                                      "BIC", "GM", "HQ", "HQc"),
                       "Result" = c(Rsq, Qsq, AdjRsq, GoF,
                                    FPE, Mallows, AIC, AICu, AICc,
                                    BIC, GM, HQ, HQc) )
  # "Range" = c(rep("[0, 1]", 4), "R", rep(NA, 4), 
  #             "R", rep(NA, 3)),
  # "Acceptance" = c(rep("Closet to 1", 4), "Smallest",
  #                  str_c("Closest to", predictors, 
  #                        sep = " "),
  #                  rep("Smallest", 7)) )
  
  result
}



###################################################################################
# COMPARE FUNCTION
###################################################################################
compare <- function(modelR_list, modelQ_list, saturated, n_list, predictors_list){
  #modelR_list = rlist
  #modelQ_list = qlist
  #saturated = plsmodelR_model1
  #n_list = nlist
  #predictors_list = plist
  
  # Confirm that the lists are all the same size
  if(all.equal(length(modelR_list), length(modelQ_list), length(n_list), 
               length(predictors_list)) == TRUE) {
    # no action
  } else {
    print("Lists are not the same size")
    exit()
  }
  
  # Compile results
  allResults <- NA
  for(i in 1:length(modelR_list)){
    if(i == 1){
      allResults <- selection(modelR_list[[i]], modelQ_list[[i]], saturated,
                              n_list[[i]], predictors_list[[i]])
      colnames(allResults)[colnames(allResults) == "Result"] <- "Saturated Model"
    }else if(i > 1){
      allResults <- cbind(allResults, 
                          selection(modelR_list[[i]], modelQ_list[[i]], 
                                    saturated, n_list[[i]],
                                    predictors_list[[i]])$Result)
    }
  }
  # Correct the column names
  colnames(allResults)[-c(1:3)] <- str_c("Model", 2:length(modelR_list), 
                                         sep = " ")
  
  # Rearrange table
  # Final = cbind(allResults[, c("Type", "Criteria")], 
  #               allResults[, colnames(allResults)[str_detect(colnames(allResults),
  #                                                            "Model")]],
  #               allResults[, c("Range", "Acceptance")])
  
  
  
  choiceModel <- c()
  
  # Closest to 1
  close1 <- allResults[allResults$Criteria %in% c("Rsq", "Qsq", "AdjRsq", "GoF"), 
                       str_detect(colnames(allResults), "Model")]
  
  
  #close1[, 1] <- close1[, 1] - 0.2
  #close1[, 3] <- close1[, 3] + 0.2
  
  #close1[1,2] <- close1[1,2] - 1
  #close1[2,3] <- close1[2,3] + 1
  
  
  if(any(apply(close1, 2, function(x) x > 1) | 
         apply(close1, 2, function(x) x < 0), na.rm = TRUE) == TRUE){
    # Change values outside of the expected range to NA so they are not considered
    close1[apply(close1, 2, function(x) x > 1) | 
             apply(close1, 2, function(x) x < 0)] <- NA
  }
  
  
  for(i in 1:nrow(close1)){
    check <- apply(close1, 1, function(x) max(x, na.rm = TRUE))[i]
    
    anyWork = apply(close1, 2, function(y) y %in% check)
    whichWorks = anyWork %>% apply(., 2, any) %>% colnames(close1)[.]
    
    
    if(length(whichWorks) == 0){
      choiceModel[i] <- "None"  
    } else if(length(whichWorks) > 1){
      choiceModel[i] <- whichWorks
    } else if(length(whichWorks) == 1){
      choiceModel[i] <- str_c(whichWorks, collapse = ", ")
    }
    
  }
  
  
  # Smallest value
  smallest <- allResults[allResults$Criteria %in% c("FPE", "AIC", "AICu", "AICc",
                                                    "BIC", "GM", "HQ", "HQc"), 
                         str_detect(colnames(allResults), "Model")]
  
  
  for(i in 1:nrow(smallest)){
    check <- apply(smallest, 1, function(x) min(x, na.rm = TRUE))[i]
    
    choiceModel[i+4] <- apply(smallest, 2, function(y) y %in% check) %>% 
      apply(., 2, any) %>% colnames(smallest)[.]
  }
  
  
  # Closest to the Predictor
  closeP <- allResults[allResults$Criteria %in% c("Mallow's Cp"), 
                       str_detect(colnames(allResults), "Model")]
  
  compare = abs(closeP[1,] - unlist(predictors_list))
  
  #compare[1, 3] <- compare[1, 3] + 0.235
  #compare[1, 1] <- compare[1, 3] + 0.576
  
  
  if(min(compare[1,]) == Inf){
    addMallow <- "None"
  } else if(length( compare[compare[1,] == min(compare[1,])] ) > 1){
    addMallow <- str_c(colnames(compare)[compare[1,] == min(compare[1,], 
                                                            na.rm =  TRUE)], collapse = ", ")
  } else if(length( compare[compare[1,] == min(compare[1,], na.rm =  TRUE)] ) == 1){
    addMallow <- colnames(compare)[compare[1,] == min(compare[1,], na.rm =  TRUE)]
  }
  
  
  # Final result
  allResults$Choice <- c(choiceModel[1:4], choiceModel[5], addMallow, choiceModel[6:12])
  
  allResults
}



###################################################################################
# ITERATE FUNCTION
###################################################################################
iterate <- function(inner_path, association_blocks, modes, saturated_index, 
                    predictors, sample, iterations){
  #inner_path = innerPath
  #association_blocks = blocks
  #modes = allModes
  #saturated_index = 1
  #sample = 100
  #iterations = 10
  #predictors = allPredictors
  
  binResults <- list()
  for (i in 1:iterations){
    rand <- sample(1:nrow(pressure_data), sample, replace = FALSE)
    data <- pressure_data[rand, ]
    
    # Add the training data indicators for comparison of R^2 and Q^2
    data$train <- rep(FALSE, nrow(data))
    data$train[sample(1:nrow(data), 0.3*nrow(data))] <- TRUE
    # NOTE: must select more than the quantity of predictors
    
    trainingData <- data[which(data$train == TRUE), ]
    remainingData <- data[which(data$train != TRUE), ]
    
    
    
    #####################################################
    # Running Models
    plsmodelR = list()
    plsmodelQ = list()
    for(j in 1:length(inner_path)){
      plsmodelR[[j]] <- plspm(trainingData, inner_path[[j]], association_blocks[[j]],
                              modes = modes[[j]])
      plsmodelQ[[j]] <- plspm(remainingData, inner_path[[saturated_index]],
                              association_blocks[[saturated_index]], 
                              modes = modes[[saturated_index]])
    }
    
    
    
    #####################################################
    # Comparisons
    rlist = plsmodelR
    qlist = plsmodelQ
    nlist = list(nrow(trainingData),nrow(trainingData), 
                 nrow(trainingData), nrow(trainingData))
    plist = predictors
    
    
    
    #save compare results of current run 
    current_run <-  compare(rlist, qlist, rlist[[1]], nlist, plist)
    
    modelsPresent <- colnames(current_run)[3:ncol(current_run)] %>% .[-length(.)]
    
    compileSelection = list()
    for(j in 1:nrow(current_run)){
      compileSelection[[j]] <- modelsPresent %in% current_run$Choice[j]
    }
    logResults = do.call(rbind, compileSelection) %>% 
      `rownames<-`(current_run$Criteria) %>% 
      `colnames<-`(modelsPresent)
    
    binResults[[i]] = apply(logResults, 2, \(x) +as.logical(x)) %>% 
      `rownames<-`(current_run$Criteria)
      # https://stackoverflow.com/questions/4605206/drop-data-frame-columns-by-name
    
    
    #save to list
    #compare_runs <- append(compare_runs, list(current_run))
  }
  
  
  #####################################################
  # Results
  startResults <- binResults[[1]]*0
  for(i in 1:length(binResults)){
    startResults <- startResults + binResults[[i]]
  }
  
  compare_runs <- apply(startResults, 2, function(x) x/iterations)
  
  #list of "compare" data.frames for various subsets of size 100
  compare_runs
}



###################################################################################
# MODELS
###################################################################################
random <- sample(1:nrow(pressure_data), 50, replace = FALSE)
data_subset <- pressure_data[random, ]


# Add the training data indicators for comparison of R^2 and Q^2
data_subset$train <- rep(FALSE, nrow(data_subset))
data_subset$train[sample(1:nrow(data_subset), 0.3*nrow(data_subset))] <- TRUE
  # NOTE: must select more than the quantity of predictors


training <- data_subset[which(data_subset$train == TRUE), ]
remaining <- data_subset[which(data_subset$train != TRUE), ]


#####################################################
# MODEL 1 (Saturated)

Socio_Economic  <- c(0,0,0,0,0)
Aging_Process   <- c(0,0,0,0,0)
Lifestyle       <- c(0,0,0,0,0)
Obesity1        <- c(1,1,1,0,0)
Hypertension1   <- c(1,1,1,1,0)

inner_path_model1 = rbind(Socio_Economic, Aging_Process, Lifestyle, Obesity1, Hypertension1)
colnames(inner_path_model1) = rownames(inner_path_model1)
association_blocks_model1 = list(4:5, 3, 1, 10:16, 6:9)
modes_model1 = c("B", "B", "A", "A", "A")

plsmodelR_model1 = plspm(training, inner_path_model1, association_blocks_model1, 
                         modes = modes_model1)
plsmodelQ_model1 = plspm(remaining, inner_path_model1, association_blocks_model1, 
                         modes = modes_model1)

results1 <- selection(plsmodelR_model1, plsmodelQ_model1, 
                      plsmodelR_model1, nrow(training), 4)



#####################################################
# MODEL 2

Socio_Economic  <- c(0,0,0,0,0)
Aging_Process   <- c(0,0,0,0,0)
Lifestyle       <- c(0,0,0,0,0)
Obesity2        <- c(0,0,0,0,0)
Hypertension2   <- c(1,1,1,1,0)

inner_path_model2 = rbind(Socio_Economic, Aging_Process, Lifestyle, Obesity2, Hypertension2)
colnames(inner_path_model2) = rownames(inner_path_model2)
association_blocks_model2 = list(4:5, 3, 1, 10:16, 6:9)
modes_model2 = c("B", "B", "A", "A", "A")

plsmodelR_model2 = plspm(training, inner_path_model2, association_blocks_model2, 
                         modes = modes_model2)
plsmodelQ_model2 = plspm(remaining, inner_path_model2, association_blocks_model2, 
                         modes = modes_model2)

results2 <- selection(plsmodelR_model2, plsmodelQ_model2, 
                      plsmodelR_model1, nrow(training), 4)




#####################################################
# MODEL 3 

Socio_Economic  <- c(0,0,0,0,0)
Aging_Process   <- c(0,0,0,0,0)
Lifestyle       <- c(0,0,0,0,0)
Obesity3        <- c(1,1,1,0,0)
Hypertension3   <- c(0,0,0,1,0)

inner_path_model3 = rbind(Socio_Economic, Aging_Process, Lifestyle, Obesity3, Hypertension3)
colnames(inner_path_model3) = rownames(inner_path_model3)
association_blocks_model3 = list(4:5, 3, 1, 10:16, 6:9)
modes_model3 = c("B", "B", "A", "A", "A")

plsmodelR_model3 = plspm(training, inner_path_model3, association_blocks_model3, 
                         modes = modes_model3)
plsmodelQ_model3 = plspm(remaining, inner_path_model3, association_blocks_model3, 
                         modes = modes_model3)

results3 <- selection(plsmodelR_model3, plsmodelQ_model3, 
                      plsmodelR_model1, nrow(training), 4)



#####################################################
# MODEL 4

Socio_Economic  <- c(0,0,0,0,0)
Aging_Process   <- c(0,0,0,0,0)
Lifestyle       <- c(0,0,0,0,0)
Obesity4        <- c(0,0,1,0,0)
Hypertension4   <- c(1,1,0,1,0)

inner_path_model4 = rbind(Socio_Economic, Aging_Process, Lifestyle, Obesity4, Hypertension4)
colnames(inner_path_model4) = rownames(inner_path_model4)
innerplot(inner_path_model4)
association_blocks_model4 = list(4:5, 3, 1, 10:16, 6:9)
modes_model4 = c("B", "B", "A", "A", "A")

plsmodelR_model4 = plspm(training, inner_path_model4, association_blocks_model4, 
                         modes = modes_model4)
plsmodelQ_model4 = plspm(remaining, inner_path_model4, association_blocks_model4, 
                         modes = modes_model4)

results4 <- selection(plsmodelR_model4, plsmodelQ_model4, 
                      plsmodelR_model1, nrow(training), 4)



###################################################################################
# GENERATE TABLES
###################################################################################
rlist = list(plsmodelR_model1, plsmodelR_model2, plsmodelR_model3)
qlist = list(plsmodelQ_model1, plsmodelQ_model2, plsmodelQ_model3)
nlist = list(50, 30, 60)
plist = list(4, 4, 4)

compare(rlist, qlist, plsmodelR_model1, nlist, plist)





innerPath <- list(inner_path_model1, inner_path_model2, inner_path_model3, inner_path_model4)
blocks <- list(association_blocks_model1, association_blocks_model2,
               association_blocks_model3, association_blocks_model4)
allModes <- list(modes_model1, modes_model2, modes_model3, modes_model4)
allPredictors = list(4, 4, 4, 4)



iterate(innerPath, blocks, allModes, saturated_index = 1, allPredictors, 
        sample = 100, iterations = 20)



###################################################################################
# GENERATE TABLES
###################################################################################
## couple of good graphics
plot(plsmodelR)
plot(plsmodelR, what = "loadings")


## cross-loadings  - pg. 78-79
xloads = melt(plsmodelR$crossloadings, id.vars = c("name", "block"),
              variable_name = "LV")

ggplot(data = xloads, aes(x = name, y = value, fill = block)) +
  geom_hline(yintercept = 0, color = "gray75") +
  geom_hline(yintercept = 0.5, color = "gray70", linetype = 2) +
  geom_bar(stat = 'identity', position = 'dodge') +
  facet_wrap(block ~ LV) +
  theme(axis.text.x = element_text(angle = 90),
        line = element_blank(),
        plot.title = element_text(size = 12)) +
  ggtitle("Crossloadings")
## Plot shows that there are no "traitor indicators"





