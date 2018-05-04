##############################################
##############################################
#
# R Code supplementing the paper
# The mechanics of predator-prey interactions: first principles of physics predict predator-prey size ratios
# by Portalier, Fussmann, Loreau, Cherif

# April 2018
#
##############################################
##############################################
#

#### READ ME #################################
#
# The following code is used to test accuracy of the model by computing the True Skill Statistics (TSS)
# over a subset of the database (data points that come from food web studies)
# 
# The code calls the FoodWebs.RData file that contains:    
# 1. body sizes and names for all species in each food web used in the analysis
# 2. the observed community matrices for each food web 
# 3. the predicted prey size ranges from the model (see main text)
# 
# For a given predator, the predicted size range for a prey includes all prey sizes 
# that lead to a positive net gain (see main text)
# 
# Then, the code builds the predicted community matrices, 
# compares observed and predicted community matrices, and computes the TSS
#
# The code shows the procedure for the first food web
# It can be repeated for the others
# 
##############################################


#### Function to build the predicted community matrix
BuildPredicted = function(n,minsize,maxsize) {
  
  NumberofSpecies = length(n)   	
  PredictedMatrix = matrix(0,nrow=NumberofSpecies,ncol=NumberofSpecies)
  
  for(i in 1:NumberofSpecies){
    for(j in 1:NumberofSpecies){
      if(n[j] > minsize[i] && n[j] < maxsize[i]){
        PredictedMatrix[i,j] = 1
      }
    }
  }
  return(PredictedMatrix)	
}

load('FoodWebs.RData')

#### Retrieve data for a given food web (food web #1)
ObservedFoodWeb=ObservedMatrices$ObservedFoodWeb1
Species=SpeciesSizes$SizeFoodWeb1
PredictedRange=PredictedRangesOfPrey$PredictedRangeOfPrey1

NumberofSpecies=nrow(Species)
BodySize=Species$BodyMass_kg
MinPreySize=PredictedRange$MinPreySize_kg
MaxPreySize=PredictedRange$MaxPreySize_kg

PredictedFoodWeb=BuildPredicted(BodySize,MinPreySize,MaxPreySize)

#### TSS 
a=0 # link observed and predicted
b=0 # link not observed but predicted
c=0 # link observed but not predicted
d=0 # link not observed and not predicted

for (i in 1:NumberofSpecies){
  for (j in 1:NumberofSpecies){
    if (ObservedFoodWeb[i,j]==1){
      if (PredictedFoodWeb[i,j]==1){ # link observed and predicted
        a=a+1
      }else{ # link observed but not predicted
        c=c+1 
      }
    }else{
      if (PredictedFoodWeb[i,j]==1){ # link not observed but predicted
        b=b+1
      }else{ # link not observed and not predicted
        d=d+1
      }
    }
  }
}

tss=(a*d-b*c)/((a+c)*(b+d))


