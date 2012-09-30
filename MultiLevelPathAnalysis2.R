# initialize the pth object
makePathObj <- function(paths, DF, covs=NULL, RFX=NULL, intercepts = TRUE, slopes = TRUE, nBootReps = 0, glmerText=NULL){
  pth = list()
  pth$paths = paths
  pth$DF = DF
  pth$varNames = NULL
  pth$nVars = 0
  pth$connectionMatrix = matrix()
  pth$coefMatrix = matrix()
  pth$IP = list()
  pth$covs = covs
  pth$RFX = RFX
  pth$intercepts = intercepts
  pth$slopes = slopes
  pth$glmerText = glmerText
  pth$modelList = list()
  pth$nBootReps = nBootReps
  return (pth)
}

# establish a list of variables that have connections among them.
getVarNames <- function(pth){
  varNames<-NULL
  for( i in 1:length(pth$paths) ){
    thisPath<-pth$paths[i]
    theseVars<-strsplit(thisPath, "->")[[1]]
    V1<-theseVars[1]
    V2<-theseVars[2]
    varNames = union(varNames, V1)
    varNames = union(varNames, V2)
  }
  pth$varNames = varNames
  pth$nVars = length(varNames)
  return(pth)
}

#create matrix of direct path connections
makeConnectionMatrix <- function(pth){
  
  connectionMatrix = matrix(0, pth$nVars, pth$nVars)
  rownames(connectionMatrix) = pth$varNames
  colnames(connectionMatrix) = pth$varNames
  
  for( i in 1:length(pth$paths) ){
    thisPath<-pth$paths[i]
    theseVars<-strsplit(thisPath, "->")[[1]]
    V1<-which(theseVars[1] == pth$varNames)
    V2<-which(theseVars[2] == pth$varNames)
    connectionMatrix[V1,V2] = 1
  }
  pth$connectionMatrix = connectionMatrix
  return(pth)
}

# create matrix of direct path coefficients	
directPathCoeffs <- function(pth){
  coefMatrix = pth$connectionMatrix
  coefMatrix[,] = 0
  
  modelList<-list()
  eq = 0
  for( i in 1:ncol(pth$connectionMatrix) ){    
    thisCol<-pth$connectionMatrix[,i]
    
    if(sum(thisCol)>0){
      IVs<-pth$varNames[thisCol==1]
      DV = pth$varNames[i]      
      IVsCombined	<- paste( c(IVs) , collapse = "+" )
      formulaText <- paste( DV, IVsCombined , sep = "~" )
      eq = eq+1
      
      for(j in 1:length(pth$covs)){
        formulaText <- paste(formulaText, pth$covs[j], sep = "+")
      }
      
      if(!is.null(pth$glmerText)){
        formulaText <- paste(formulaText, pth$glmerText, sep = "+")
      }
           
      for(j in 1:length(pth$RFX)){
       
        if(pth$intercepts){
          thisRandEffect<-paste('(1|',pth$RFX[j],')', sep = "")
          formulaText <- paste(formulaText, thisRandEffect, sep = "+")
        }
        
        if(pth$slopes){
          for(k in 1:length(IVs)){
            thisSlope<-paste('(', IVs[k], '+0|', pth$RFX[k], ')', sep = "")
            formulaText <- paste(formulaText, thisSlope, sep = "+")
          }
        }
      }
                  
      thisFormula <- as.formula(formulaText)
      
      if(!is.null(pth$glmerText) | pth$intercepts | pth$slopes){    
        thisMod<-glmer(thisFormula, pth$DF) # run a mixed glm
        theseCoef<-attr(thisMod, "fixef")
      } else {
        thisMod<-lm(thisFormula, pth$DF)
        theseCoef <- coefficients(summary(thisMod))[,"Estimate"] #otherwise run a normal linear model
      }
      modelList[[eq]] = thisMod
      colsToReplace = which(thisCol==1)
      coefMatrix[colsToReplace,i] = theseCoef[2:(1+length(colsToReplace))]
    }
  }
  pth$coefMatrix = coefMatrix
  pth$modelList = modelList
  return(pth)
}

# a wrapper to call the recursive function "findIndirectPaths"
findIndirectPathsWrapper <- function(pth) {
  pth$ix<-0
  pth = findIndirectPaths(pth, pth$varNames)
  pth$thisIP<-NULL
  pth$ix<-NULL  
  return(pth)
}

# recursively find all the indirect paths
findIndirectPaths <- function(pth, varsToSearch){  
  
  for(v in varsToSearch){	
    
    i = which(pth$varNames==v)
    thisRow<-pth$coefMatrix[i,]
    thisVar<-pth$varNames[i]
    pth$thisIP = c(pth$thisIP, thisVar)
    
    if (length(unique(pth$thisIP)) != length(pth$thisIP)) {
      badPath = paste(pth$thisIP, collapse = "->")
      errText = paste("Cyclic path found: ", badPath, ". Path connections must be cyclic.", sep="")
      stop(errText)
    }
    
    if (length(pth$thisIP) > 2){      
      pth$ix <- pth$ix + 1 # update path index
      pth$IP[[pth$ix]] = list()
      pth$IP[[pth$ix]]$path = pth$thisIP
    }
    
    if(sum(abs(thisRow))>0){
      newVarsToSearch = pth$varNames[which(thisRow!=0)]
      newCoefMatrix = pth$coefMatrix
      newCoefMatrix[i,] = 0
      pth = findIndirectPaths(pth, newVarsToSearch)
    }
    pth$thisIP = pth$thisIP[-length(pth$thisIP)]
  }
  return(pth)
}


# compute all indirect path coefficients
indirectPathCoefs <- function(pth){

  allIndirectPathCoef = list()
  for( i in 1:length(pth$IP) ){
    thisPath = pth$IP[[i]]$path
    thisIndCoef = 1
    for( j in 1:(length(thisPath)-1) ){
      
      V1<-which(thisPath[j] == pth$varNames)
      V2<-which(thisPath[j+1] == pth$varNames)
      
      thisDirCoef = pth$coefMatrix[V1,V2];
      thisIndCoef = thisIndCoef * thisDirCoef      
    }
    allIndirectPathCoef[[i]] = thisIndCoef
    pth$IP[[i]]$coef = thisIndCoef
  }
  return(pth)
}	

# Find bootstrap confidence intervals and p-values
bootStrapCoefs <- function(pth){
  dir_paths = array(0, c(dim(pth$coefMatrix), pth$nBootReps))
  ind_paths = NULL
  for(i in 1:pth$nBootReps){
    boot_DF 	<- pth$DF[sample( 1:nrow(pth$DF), replace=TRUE),]
    boot_pth = pth
    boot_pth$DF = boot_DF
    boot_pth = directPathCoeffs(boot_pth)
    boot_pth = indirectPathCoefs(boot_pth)
    theseIndCoefs = sapply(boot_pth$IP, function(x) getElement(x,"coef"))
    ind_paths 	<- as.data.frame(rbind(ind_paths, theseIndCoefs))
    dir_paths[,,i]<-boot_pth$coefMatrix
  }
  idxCI = c(round(.025*nrow(ind_paths)+.5), round(.975*nrow(ind_paths)+.5))
  
  #pth$bootDirect$dist = dir_paths
  pth$bootDirect$CI95Pct = array(list(NULL), dim(pth$coefMatrix))
  rownames(pth$bootDirect$CI95Pct) = rownames(pth$coefMatrix)
  colnames(pth$bootDirect$CI95Pct) = colnames(pth$coefMatrix)
  #pth$bootDirect$p = array(0, dim(pth$coefMatrix))
  pth$bootDirect$p = pth$coefMatrix
  pth$bootDirect$p[,] = 0

  for(i in 1:pth$nVars){
    for(j in 1:pth$nVars){
      if (pth$connectionMatrix[i,j]==1){
      thisSortedVec = sort(dir_paths[i,j,])
      pth$bootDirect$CI95Pct[[i,j]] = c(thisSortedVec[idxCI[1]],thisSortedVec[idxCI[2]])
      dist_prop<- as.numeric(pth$coefMatrix[i,j]>thisSortedVec)
      p_h = 2*abs(.5 - sum(dist_prop)/length(dist_prop))
      p = sapply(p_h, function(x) min(c(x, 1-x)))
      if (p==0) {
        p = 1/pth$nBootReps
      }
      pth$bootDirect$p[i,j] = p
      }
    }
  }

  for(i in 1:length(ind_paths)){
    sortedDists = sort(ind_paths[[i]])
    distProp = as.numeric(pth$IP[[i]]$coef > sortedDists)
    p_h = 2*abs(.5 - sum(distProp)/length(distProp))
    
    if (p_h>0) {
      pth$IP[[i]]$p = p
    }else {
      pth$IP[[i]]$p = 1/pth$nBootReps
    }
    #pth$IP[[i]]$dist = ind_paths[[i]]
    pth$IP[[i]]$CI95pct = sortedDists[idxCI]
  }

  pth$bootDirect$dist = NULL

return(pth)
}


# The main function.  Run a path analysis!
pathAnalysis <- function(paths, DF, covs=NULL, RFX=NULL, intercepts = TRUE, slopes = TRUE, nBootReps = 0, glmerText=NULL){
  
  pth = makePathObj(paths, DF, covs, RFX, intercepts, slopes, nBootReps, glmerText)
  pth = getVarNames(pth)
  pth = makeConnectionMatrix(pth)
  pth = directPathCoeffs(pth)
  pth = findIndirectPathsWrapper(pth)  
  pth = indirectPathCoefs(pth)
  pth = bootStrapCoefs(pth)
  
  return(pth)
}

#Usage Case
n	<- 100
DV	<- rnorm( n )
MV2	<- DV + rnorm( n , 2 )
MV1	<- MV2 + rnorm( n , 2 )
c1  <- MV2 + rnorm( n , 2 )
IV	<- MV1 + rnorm( n , 2 )
subs <- round(2*runif(100))
DF	<- data.frame( DV , MV2, MV1, IV, c1, subs)


#paths = c('IV->MV1', 'IV->MV2', 'MV1->DV', 'MV2->DV')
paths = c('IV->MV1', 'MV1->MV2', 'MV2->DV')

covs = ("c1")
RFX = "subs"
pathRes<-pathAnalysis(paths, DF, covs, RFX, nBootReps = 100)

##TO DO: 

# make sure the p-values for the bootstrap are correct.
# organize direct path stuff
# check for crazy and non-semantic input
# allow for lists of covariates and glmerText, in case each direct path requires something special. 
# thoroughly check other test cases
# is there a way to suppress the matrix of bootstrapped values from appearing when the pathRes variable is typed?  (e.g. slots?) If not, consider not storing the bootstrapped values...