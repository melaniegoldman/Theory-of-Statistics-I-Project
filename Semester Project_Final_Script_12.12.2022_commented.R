###################################################################################
# SET UP
###################################################################################
library("rmarkdown")        # version 2.16
library("reshape")          # version 0.8.9
library("reshap2")          # version 1.4.4
library("githubinstall")    # version 0.2.2
library("dplyr")            # version 1.0.10
library("tidyverse")        # version 1.3.2
library("ggplot2")          # version 3.3.6
library("stringr")          # version 1.5.0
library("lemon")            # version 0.4.5
library("DiagrammeR")       # version 1.0.9
library("webshot")          # version 0.5.4
# webshot::install_phantomjs(force = TRUE)
# install.packages("devtools")
# library("devtools")       # version 2.4.5
# githubinstall("plspm")
library("plspm")            # version 0.4.9
set.seed(22)


"%!in%" <- function(x,y)!('%in%'(x,y))



###################################################################################
# IMPORT DATASET
###################################################################################
pressure <- read.csv(file.choose())

# Remove respondent number
pressure1 <- pressure[, names(pressure) %!in% c("SEQN")]

# Removed census data                                       - #3
pressure1 <- pressure1[, names(pressure1) %!in% c("RIDRETH1", "DMDCITZN", "DMDEDUC2","DMDMARTL", "WTINT2YR", "SDMVPSU", "SDMVSTRA")]

# Removed categorical data
pressure1 <- pressure1[, names(pressure1) %!in% c("ALQ101", "ALQ110", "SMQ020", "HIQ210")]

# Remove all rows with an NA
pressure_data <- na.omit(pressure1)


###################################################################################
# SELECTION FUNCTION
###################################################################################
selection <- function(model_R, model_Q, saturated, n, predictors){
  # Finds the row associated with the hypertension latent variable
  r2_R = which(str_detect(rownames(model_R$inner_summary), "Hypertension"))
  r2_Q = which(str_detect(rownames(model_Q$inner_summary), "Hypertension"))
  r2_sat = which(str_detect(rownames(saturated$inner_summary), "Hypertension"))
  
  # Calculated values used in the criteria equations
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
  Mallows   = (SSE/MSE) - (n - pk)
  AIC       = n*(log(SSE/n) + ((2*pk)/n))
  AICu      = n*(log(SSE/(n - pk)) + ((2*pk)/n))
  AICc      = n*(log(SSE/(n)) + ((n + pk)/(n - pk - 2)))
  
  ## Consistent criteria
  BIC       = n*(log(SSE/n) + ((pk*log(n))/n))
  GM        = (SSE/MSE) + (pk*log(n))
  HQ        = n*(log(SSE/n) + ((2*pk*log(log(n)) )/n) )
  HQc       = n*(log(SSE/n) + ((2*pk*log(log(n)) )/(n - pk - 2)) )
  
  
  ## Compiles the results into a table
  result <- data.frame("Type" = c(rep("PLS selection", 4), 
                                  rep("Distance based", 5), 
                                  rep("Consistent criteria", 4)),
                       "Criteria" = c("Rsq", "Qsq", "AdjRsq", "GoF",
                                      "FPE", "Mallow's Cp", "AIC", "AICu", "AICc",
                                      "BIC", "GM", "HQ", "HQc"),
                       "Result" = c(Rsq, Qsq, AdjRsq, GoF,
                                    FPE, Mallows, AIC, AICu, AICc,
                                    BIC, GM, HQ, HQc) )
  
  result
}



###################################################################################
# COMPARE FUNCTION
###################################################################################
compare <- function(modelR_list, modelQ_list, saturated, n_list, predictors_list){

  # Quick check to confirm that the lists are all the same size
  if(all.equal(length(modelR_list), length(modelQ_list), length(n_list), 
               length(predictors_list)) == TRUE) {
    # no action
  } else {
    print("Lists are not the same size")
    exit()
  }
  
  # Compile the selection() results for all models into one table
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
  
  
  
  #####################################################
  # Ranking models by criteria conditions
  
  choiceModel <- c()
  
  # Closest to 1
    close1 <- allResults[allResults$Criteria %in% c("Rsq", "Qsq", "AdjRsq", "GoF"), 
                         str_detect(colnames(allResults), "Model")]
    
    # Change values outside of the expected range to NA so they are not considered
    if(any(apply(close1, 2, function(x) x > 1) | 
           apply(close1, 2, function(x) x < 0), na.rm = TRUE) == TRUE){
      
      close1[apply(close1, 2, function(x) x > 1) | 
               apply(close1, 2, function(x) x < 0)] <- NA
      
    }
    
    # 
    for(i in 1:nrow(close1)){
      # Find the value that meets the condition best
      check <- apply(close1, 1, function(x) max(x, na.rm = TRUE))[i]
      
      # Locates this value within the matrix
      anyWork = apply(close1, 2, function(y) y %in% check)
      whichWorks = anyWork %>% apply(., 2, any) %>% colnames(close1)[.]
      
      # Prepares the final notation based on how many models match the criteria best
      if(length(whichWorks) == 0){
        choiceModel[i] <- "None"                      # No models were closest to 1
      } else if(length(whichWorks) == 1){
        choiceModel[i] <- whichWorks                  # Only one model was closes to 1
      } else if(length(whichWorks) > 1){
        choiceModel[i] <- str_c(whichWorks, 
                                collapse = ", ")      # Multiple models tied
      }
      
    }
  
  
  # Smallest value
    smallest <- allResults[allResults$Criteria %in% c("FPE", "AIC", "AICu", "AICc",
                                                      "BIC", "GM", "HQ", "HQc"), 
                           str_detect(colnames(allResults), "Model")]
    
    # Finds the value that is smallest by row (aka by criteria)
    for(i in 1:nrow(smallest)){
      check <- apply(smallest, 1, function(x) min(x, na.rm = TRUE))[i]
      
      choiceModel[i+4] <- apply(smallest, 2, function(y) y %in% check) %>% 
        apply(., 2, any) %>% colnames(smallest)[.]
    }
  
  
  # Closest to the Predictor
    closeP <- allResults[allResults$Criteria %in% c("Mallow's Cp"), 
                         str_detect(colnames(allResults), "Model")]
    
    # Subtracts the predictor value from each Mallow's Cp to compare the min absolute value
    compare = abs(closeP[1,] - unlist(predictors_list))
    
    # Prepares the final notation based on how many models match the criteria best
    # Similar set-up as in the "Closest to 1" section
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
  
  binResults <- list()
  for (i in 1:iterations){
    # Randomly select a subset from the dataset of the size "sample"
    rand <- sample(1:nrow(pressure_data), sample, replace = FALSE)
    data <- pressure_data[rand, ]
    
    # Add the training data indicators for comparison of R^2 and Q^2
    data$train <- rep(FALSE, nrow(data))
    data$train[sample(1:nrow(data), 0.7*nrow(data))] <- TRUE
      # NOTE: must select more than the quantity of predictors
    
    trainingData <- data[which(data$train == TRUE), ]
    remainingData <- data[which(data$train != TRUE), ]
    
    
    
    #####################################################
    # Running Models
    plsmodelR = list()
    plsmodelQ = list()
    
    # Generate the model for each set-up indicated by the user
    for(j in 1:length(inner_path)){
      plsmodelR[[j]] <- plspm(trainingData, inner_path[[j]], association_blocks[[j]],
                              modes = modes[[j]])
      plsmodelQ[[j]] <- plspm(remainingData, inner_path[[j]],
                              association_blocks[[j]], 
                              modes = modes[[j]])
    }
    
    
    
    #####################################################
    # Comparisons
    rlist = plsmodelR
    qlist = plsmodelQ
    nlist = as.list(rep(nrow(trainingData),length(rlist)))
    plist = predictors
    
    
    # Save compare() results of current iteration 
    current_run <-  compare(rlist, qlist, rlist[[1]], nlist, plist)
    
    # Save the model names for the ones being evaluated. Used for boolean, match testing
    modelsPresent <- colnames(current_run)[3:ncol(current_run)] %>% .[-length(.)]
    
    # Boolean, name matching to count selection rate
    compileSelection = list()
    for(j in 1:nrow(current_run)){
      # Accounts for multiple models chosen
      compileSelection[[j]] <- modelsPresent %in% str_split_1(current_run$Choice[j], ", ")
    }
    logResults = do.call(rbind, compileSelection) %>% 
      `rownames<-`(current_run$Criteria) %>% 
      `colnames<-`(modelsPresent)
    
    # Changes the logical TRUE/FALSE to a 1/0 binary for rate calculation
    binResults[[i]] = apply(logResults, 2, \(x) +as.logical(x)) %>% 
      `rownames<-`(current_run$Criteria)
      # Source: https://stackoverflow.com/questions/4605206/drop-data-frame-columns-by-name
    
    
    # Save to list
    save_compares <<- binResults
  }
  
  
  #####################################################
  # Results
  
  # Adds each bin result matrices to show the number of times a model was
  # selected by criteria
  startResults <- binResults[[1]]*0
  for(i in 1:length(binResults)){
    startResults <- startResults + binResults[[i]]
  }
  
  # Divides the selection rate by total of iterations for the final percent hits
  compare_runs <- apply(startResults, 2, function(x) x/iterations)
  
  # List of "compare" data.frames for various subsets of size 100
  compare_runs
}



###################################################################################
# MODELS
###################################################################################
random <- sample(1:nrow(pressure_data), 50, replace = FALSE)
data_subset <- pressure_data[random, ]


# Add the training data indicators for comparison of R^2 and Q^2
data_subset$train <- rep(FALSE, nrow(data_subset))
data_subset$train[sample(1:nrow(data_subset), 0.7*nrow(data_subset))] <- TRUE
  # NOTE: must select more than the quantity of predictors


training <- data_subset[which(data_subset$train == TRUE), ]
remaining <- data_subset[which(data_subset$train != TRUE), ]


#####################################################
# MODEL 1 (Saturated)

# Matrix that dictates the interactions of each latent
# variable on the others where column j affects row i.
# i.e. Obesity affects Hypertension
Socio_Economic1  <- c(0,0,0,0,0)
Aging_Process1   <- c(0,0,0,0,0)
Lifestyle1       <- c(0,0,0,0,0)
Obesity1         <- c(1,1,1,0,0)
Hypertension1    <- c(1,1,1,1,0)

inner_path_model1 = rbind(Socio_Economic1, Aging_Process1, Lifestyle1, Obesity1, Hypertension1)
colnames(inner_path_model1) = rownames(inner_path_model1)

# Defines the manifest variables that are associated with a 
# given latent variable.
# i.e. Hypertension manifest variables are columns 6-9
association_blocks_model1 = list(4:5, 3, 1, 10:16, 6:9)

# Defines each mode that the latent its respective manifest variables
# interact: A = reflective or B = formative
modes_model1 = c("B", "B", "A", "A", "A")

# Generates the model using the set-up defined above
# NOTE: one is for R^2 calculation and the other for Q^2
plsmodelR_model1 = plspm(training, inner_path_model1, association_blocks_model1, 
                         modes = modes_model1)
plsmodelQ_model1 = plspm(remaining, inner_path_model1, association_blocks_model1, 
                         modes = modes_model1)

# Model results
results1 <- selection(plsmodelR_model1, plsmodelQ_model1, 
                      plsmodelR_model1, nrow(training), 4)

# Plots the latent variable associations for the inner model
innerplot(inner_path_model1)



#####################################################
# MODEL 2

Socio_Economic2  <- c(0,0,0,0,0)
Aging_Process2   <- c(0,0,0,0,0)
Lifestyle2       <- c(0,0,0,0,0)
Obesity2         <- c(0,0,0,0,0)
Hypertension2    <- c(1,1,1,1,0)

inner_path_model2 = rbind(Socio_Economic2, Aging_Process2, Lifestyle2, Obesity2, Hypertension2)
colnames(inner_path_model2) = rownames(inner_path_model2)
association_blocks_model2 = list(4:5, 3, 1, 10:16, 6:9)
modes_model2 = c("B", "B", "A", "A", "A")


plsmodelR_model2 = plspm(training, inner_path_model2, association_blocks_model2, 
                         modes = modes_model2)
plsmodelQ_model2 = plspm(remaining, inner_path_model2, association_blocks_model2, 
                         modes = modes_model2)


results2 <- selection(plsmodelR_model2, plsmodelQ_model2, 
                      plsmodelR_model1, nrow(training), 4)


innerplot(inner_path_model2)



#####################################################
# MODEL 3 

Socio_Economic3  <- c(0,0,0,0,0)
Aging_Process3   <- c(0,0,0,0,0)
Lifestyle3       <- c(0,0,0,0,0)
Obesity3         <- c(1,1,1,0,0)
Hypertension3    <- c(0,0,0,1,0)

inner_path_model3 = rbind(Socio_Economic3, Aging_Process3, Lifestyle3, Obesity3, Hypertension3)
colnames(inner_path_model3) = rownames(inner_path_model3)
association_blocks_model3 = list(4:5, 3, 1, 10:16, 6:9)
modes_model3 = c("B", "B", "A", "A", "A")

plsmodelR_model3 = plspm(training, inner_path_model3, association_blocks_model3, 
                         modes = modes_model3)
plsmodelQ_model3 = plspm(remaining, inner_path_model3, association_blocks_model3, 
                         modes = modes_model3)

results3 <- selection(plsmodelR_model3, plsmodelQ_model3, 
                      plsmodelR_model1, nrow(training), 4)


innerplot(inner_path_model3)



#####################################################
# MODEL 4

Socio_Economic4  <- c(0,0,0,0,0)
Aging_Process4   <- c(0,0,0,0,0)
Lifestyle4       <- c(0,0,0,0,0)
Obesity4         <- c(0,0,1,0,0)
Hypertension4    <- c(1,1,0,1,0)

inner_path_model4 = rbind(Socio_Economic4, Aging_Process4, Lifestyle4, Obesity4, Hypertension4)
colnames(inner_path_model4) = rownames(inner_path_model4)
association_blocks_model4 = list(4:5, 3, 1, 10:16, 6:9)
modes_model4 = c("B", "B", "A", "A", "A")

plsmodelR_model4 = plspm(training, inner_path_model4, association_blocks_model4, 
                         modes = modes_model4)
plsmodelQ_model4 = plspm(remaining, inner_path_model4, association_blocks_model4, 
                         modes = modes_model4)

results4 <- selection(plsmodelR_model4, plsmodelQ_model4, 
                      plsmodelR_model1, nrow(training), 4)


innerplot(inner_path_model4)




#####################################################
# MODEL 5 (manifest variable change of Model 2)

Aging_Process5   <- c(0,0,0,0)
Lifestyle5       <- c(0,0,0,0)
Obesity5         <- c(0,0,0,0)
Hypertension5    <- c(1,1,1,0)

inner_path_model5 = rbind(Aging_Process5, Lifestyle5, Obesity5, Hypertension5)
colnames(inner_path_model5) = rownames(inner_path_model5)
association_blocks_model5 = list(3, 1, 10:16, 6:9)
modes_model5 = c("B", "A", "A", "A")


plsmodelR_model5 = plspm(training, inner_path_model5, association_blocks_model5, 
                         modes = modes_model5)
plsmodelQ_model5 = plspm(remaining, inner_path_model5, association_blocks_model5, 
                         modes = modes_model5)

results5 <- selection(plsmodelR_model5, plsmodelQ_model5, 
                      plsmodelR_model2, nrow(training), 4)


innerplot(inner_path_model5)



#####################################################
# MODEL 6 (manifest variable change of Model 2)

Socio_Economic6  <- c(0,0,0,0)
Lifestyle6       <- c(0,0,0,0)
Obesity6         <- c(0,0,0,0)
Hypertension6    <- c(1,1,1,0)

inner_path_model6 = rbind(Socio_Economic6, Lifestyle6, Obesity6, Hypertension6)
colnames(inner_path_model6) = rownames(inner_path_model6)
association_blocks_model6 = list(4:5, 1, 10:16, 6:9)
modes_model6 = c("B", "A", "A", "A")


plsmodelR_model6 = plspm(training, inner_path_model6, association_blocks_model6, 
                         modes = modes_model6)
plsmodelQ_model6 = plspm(remaining, inner_path_model6, association_blocks_model6, 
                         modes = modes_model6)

results6 <- selection(plsmodelR_model6, plsmodelQ_model6, 
                      plsmodelR_model2, nrow(training), 4)


innerplot(inner_path_model6)



#####################################################
# MODEL 7 (manifest variable change of Model 2)

Socio_Economic7  <- c(0,0,0,0)
Aging_Process7   <- c(0,0,0,0)
Obesity7         <- c(0,0,0,0)
Hypertension7    <- c(1,1,1,0)

inner_path_model7 = rbind(Socio_Economic7, Aging_Process7, Obesity7, Hypertension7)
colnames(inner_path_model7) = rownames(inner_path_model7)
association_blocks_model7 = list(4:5, 3, 10:16, 6:9)
modes_model7 = c("B", "B", "A", "A")


plsmodelR_model7 = plspm(training, inner_path_model7, association_blocks_model7, 
                         modes = modes_model7)
plsmodelQ_model7 = plspm(remaining, inner_path_model7, association_blocks_model7, 
                         modes = modes_model7)

results7 <- selection(plsmodelR_model7, plsmodelQ_model7, 
                      plsmodelR_model2, nrow(training), 4)


innerplot(inner_path_model7)




###################################################################################
# GENERATE TABLES
###################################################################################
rlist = list(plsmodelR_model1, plsmodelR_model2, plsmodelR_model3, 
             plsmodelR_model4, plsmodelR_model5)
qlist = list(plsmodelQ_model1, plsmodelQ_model2, plsmodelQ_model3,
             plsmodelQ_model4, plsmodelQ_model5)
nlist = list(nrow(training), nrow(training), nrow(training), 
             nrow(training), nrow(training))
plist = list(4, 4, 4, 4, 4)


compare(rlist, qlist, plsmodelR_model1, nlist, plist)



###################################################################################
# First round of model comparisons - latent variables
innerPath <- list(inner_path_model1, inner_path_model2, inner_path_model3, 
                  inner_path_model4)
blocks <- list(association_blocks_model1, association_blocks_model2,
               association_blocks_model3, association_blocks_model4)
allModes <- list(modes_model1, modes_model2, modes_model3, modes_model4)
allPredictors = list(4, 4, 4, 4)



it <- iterate(innerPath, blocks, allModes, saturated_index = 1, allPredictors, 
        sample = 100, iterations = 100)




###################################################################################
# Second round of model comparisons - manifest variables
innerPath <- list(inner_path_model2, inner_path_model5, inner_path_model6, inner_path_model7)
blocks <- list(association_blocks_model2, association_blocks_model5,
               association_blocks_model6, association_blocks_model7)
allModes <- list(modes_model2, modes_model5, modes_model6, modes_model7)
allPredictors = list(4, 3, 3, 3)



it <- iterate(innerPath, blocks, allModes, saturated_index = 1, allPredictors, 
              sample = 100, iterations = 100)




###################################################################################
# OTHER
###################################################################################
## couple of good graphics
plot(plsmodelR_model1)
plot(plsmodelR_model1, what = "loadings", arr.width = 0.1, lcol = "black")
plot(plsmodelR_model2, what = "loadings", arr.width = 0.1)
plot(plsmodelR_model3, what = "loadings", arr.width = 0.1)


## cross-loadings  - pg. 78-79
xloads = melt(plsmodelR_model1$crossloadings, id.vars = c("name", "block"),
              variable_name = "LV")

ggplot(data = xloads, aes(x = name, y = value, fill = block)) +
  geom_hline(yintercept = 0, color = "gray75") +
  geom_hline(yintercept = 0.5, color = "gray70", linetype = 2) +
  geom_bar(stat = 'identity', position = 'dodge') +
  facet_wrap(block ~ variable) +
  theme(axis.text.x = element_text(angle = 90),
        line = element_blank(),
        plot.title = element_text(size = 12)) +
  ggtitle("Crossloadings")
## Plot shows that there are no "traitor indicators"


