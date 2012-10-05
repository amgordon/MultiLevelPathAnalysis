# MultiLevel Path Analysis
# Alan Gordon, September 2012

require("lme4")

# initialize the pth object
makePathObj <- function(paths, DF, covs=NULL, RFX=NULL, intercepts = TRUE, slopes = TRUE, nBootReps = 0, glmerText=NULL){
  pth = list()
  pth$init = list()
  pth$init$paths = paths
  pth$init$DF = DF
  pth$init$varNames = NULL
  pth$init$nVars = 0
  pth$DP = list()
  pth$DP$connectionMatrix = matrix()
  pth$DP$coefMatrix = matrix()
  pth$IP = list()
  pth$init$covs = covs
  pth$init$RFX = RFX
  pth$init$intercepts = intercepts
  pth$init$slopes = slopes
  pth$init$glmerText = glmerText
  pth$DP$modelList = list()
  pth$init$nBootReps = nBootReps
  return (pth)
}

# establish a list of variables that have connections among them.
getVarNames <- function(pth){
  varNames<-NULL
  for( i in 1:length(pth$init$paths) ){
    thisPath<-pth$init$paths[i]
    theseVars<-strsplit(thisPath, "->")[[1]]
    V1<-theseVars[1]
    V2<-theseVars[2]
    varNames = union(varNames, V1)
    varNames = union(varNames, V2)
  }
  pth$init$varNames = varNames
  pth$init$nVars = length(varNames)
  return(pth)
}

#create matrix of direct path connections
makeConnectionMatrix <- function(pth){
  
  connectionMatrix = matrix(FALSE, pth$init$nVars, pth$init$nVars)
  rownames(connectionMatrix) = pth$init$varNames
  colnames(connectionMatrix) = pth$init$varNames
  
  for( i in 1:length(pth$init$paths) ){
    thisPath<-pth$init$paths[i]
    theseVars<-strsplit(thisPath, "->")[[1]]
    V1<-which(theseVars[1] == pth$init$varNames)
    V2<-which(theseVars[2] == pth$init$varNames)
    connectionMatrix[V1,V2] = TRUE
  }
  pth$DP$connectionMatrix = connectionMatrix
  return(pth)
}

# run a regression
runRegression <- function(pth, DV, IVs){
  
  IVsCombined  <- paste( c(IVs) , collapse = "+" )
  formulaText <- paste( DV, IVsCombined , sep = "~" )  
  
  for(j in 1:length(pth$init$covs)){
    formulaText <- paste(formulaText, pth$init$covs[j], sep = "+")
  }
  
  for(j in 1:length(pth$init$RFX)){
    
    if(pth$init$intercepts){
      thisRandEffect<-paste('(1|',pth$init$RFX[j],')', sep = "")
      formulaText <- paste(formulaText, thisRandEffect, sep = "+")
    }
    
    if(pth$init$slopes){
      for(k in 1:length(IVs)){

        thisSlope<-paste('(', IVs[k], '+0|', pth$init$RFX[j], ')', sep = "")
        formulaText <- paste(formulaText, thisSlope, sep = "+")
        
      }
    }
  }
  
  if(!is.null(pth$init$glmerText)){
    formulaText <- paste(formulaText, pth$init$glmerText, sep = "+")
  }

  thisFormula <- as.formula(formulaText)

  if(!is.null(pth$init$glmerText) | pth$init$intercepts | pth$init$slopes){    
    thisMod<-do.call("glmer", args=list(thisFormula, data = as.name("DF"), REML=FALSE))
    theseCoef<-attr(thisMod, "fixef")
  } else {
    thisMod<-lm(thisFormula, pth$init$DF)
    theseCoef <- coefficients(summary(thisMod))[,"Estimate"] #otherwise run a normal linear model
  }

  output = list()
  output$mod = thisMod
  output$coef = theseCoef
  return(output)
}


# create matrix of direct path coefficients	
directPathCoeffs <- function(pth){
  coefMatrix = pth$DP$connectionMatrix
  coefMatrix[,] = 0
  modelList<-list()
  eq = 0
  for( i in 1:ncol(pth$DP$connectionMatrix) ){    
    thisCol<-pth$DP$connectionMatrix[,i]
    if(any(thisCol)){
      IVs<-pth$init$varNames[thisCol]
      DV = pth$init$varNames[i]      
      
      thisReg = runRegression(pth, DV, IVs)
      eq = eq+1
      modelList[[eq]] = thisReg$mod
      colsToReplace = which(thisCol)
      coefMatrix[colsToReplace,i] = thisReg$coef[2:(1+length(colsToReplace))]
    }
  }
  pth$DP$coefMatrix = coefMatrix
  pth$DP$modelList = modelList
  return(pth)
}

# use a log-likelihood ratio test to determine the significance of a coefficient in a mer object
LRT <- function(theModel, idx=1){
 
  modFull = theModel
  modelTerms_h<-terms(theModel)  
  modelTerms = attr(modelTerms_h, "term.labels")
  
  updateStr<-paste(".~.-",as.character(modelTerms[idx]))
  modNested = update(theModel, updateStr)  
  llFull = logLik(modFull)
  llNested = logLik(modNested)
  
  chsqval<-2*(llFull[1] - llNested[1])
  dfDiff<- attr(llFull,"df") - attr(llNested,"df") 
  
  pVal<-1-pchisq(chsqval,dfDiff)
 
  return(pVal)
}


# perform a d-sep test of model fit
dSepTest <- function(pth){
  idxBf = 0
  pVals = 0
  pth$MF = list()
  for( i in 2:ncol(pth$DP$connectionMatrix) ){
    for( j in 1:(i-1) ){
      if((!pth$DP$connectionMatrix[i,j]) & (!pth$DP$connectionMatrix[j,i])){
        idxBf = idxBf+1
        V1 = pth$init$varNames[i]
        V2 = pth$init$varNames[j]
        theseVars = c(V1, V2)
        precedingVars1 = pth$init$varNames[pth$DP$connectionMatrix[,i]]
        precedingVars2 = pth$init$varNames[pth$DP$connectionMatrix[,j]]
        precedingVars = union(precedingVars1, precedingVars2)

        indPaths = mapply(function(x) x$path, pth$IP)
        indPathsContainsVars = mapply(function(x) is.element(V1,x) & is.element(V2,x), indPaths)
        
        # do any indirect paths contain both variables?
        if (any(indPathsContainsVars)){
          
          # if so, the variable located earlier in the path is the IV
          # and the variable located later is the DV
          thisPath = unlist(indPaths[which(indPathsContainsVars)[1]])
          match1 = match(V1,thisPath)
          match2 = match(V2,thisPath)
          if(match1 > match2){
            thisReg = runRegression(pth, V1, c(V2, precedingVars))
          }else{
            thisReg = runRegression(pth, V2, c(V1, precedingVars))
          }
          pVals[idxBf] = LRT(thisReg$mod)
        }  else {
          # if no indirect path contains both variables, do the regression both ways
          # and take the mean p-value across the regressions.
          thisReg1 = runRegression(V1, c(V2, precedingVars))
          thisReg2 = runRegression(V2, c(V1, precedingVars))
          pVals[idxBf] = mean(c(LRT(thisReg1), LRT(thisReg2)))
        }
      }
    }
  }
  C = -2*sum(log(pVals))
  pth$MF$basisPVals = pVals
  pth$MF$C = C
  pth$MF$df = 2*length(pVals)
  pth$MF$p = 1-pchisq(C,2*length(pVals))
  return(pth) 
}

# a wrapper to call the recursive function "findIndirectPaths"
findIndirectPathsWrapper <- function(pth) {
  pth$ix<-0
  pth = findIndirectPaths(pth, pth$init$varNames)
  pth$thisIP<-NULL
  pth$ix<-NULL  
  return(pth)
}

# recursively find all the indirect paths
findIndirectPaths <- function(pth, varsToSearch){  
  
  for(v in varsToSearch){	
    
    i = which(pth$init$varNames==v)
    thisRow<-pth$DP$coefMatrix[i,]
    thisVar<-pth$init$varNames[i]
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
      newVarsToSearch = pth$init$varNames[which(thisRow!=0)]
      newCoefMatrix = pth$DP$coefMatrix
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
      
      V1<-which(thisPath[j] == pth$init$varNames)
      V2<-which(thisPath[j+1] == pth$init$varNames)
      
      thisDirCoef = pth$DP$coefMatrix[V1,V2];
      thisIndCoef = thisIndCoef * thisDirCoef      
    }
    allIndirectPathCoef[[i]] = thisIndCoef
    pth$IP[[i]]$coef = thisIndCoef
  }
  return(pth)
}	

# Find bootstrap confidence intervals and p-values
bootStrapCoefs <- function(pth){
  dir_paths = array(0, c(dim(pth$DP$coefMatrix), pth$init$nBootReps))
  ind_paths = NULL
  for(i in 1:pth$init$nBootReps){
    boot_DF 	<- pth$init$DF[sample( 1:nrow(pth$init$DF), replace=TRUE),]
    boot_pth = pth
    boot_pth$init$DF = boot_DF
    boot_pth = directPathCoeffs(boot_pth)
    boot_pth = indirectPathCoefs(boot_pth)
    theseIndCoefs = sapply(boot_pth$IP, function(x) getElement(x,"coef"))
    ind_paths 	<- as.data.frame(rbind(ind_paths, theseIndCoefs))
    dir_paths[,,i]<-boot_pth$DP$coefMatrix
  }
  idxCI = c(floor(.025*nrow(ind_paths)), ceiling(.975*nrow(ind_paths)))
  
  #pth$DP$bootDirect$dist = dir_paths
  pth$DP$bootDirect$CI95Pct = array(list(NULL), dim(pth$DP$coefMatrix))
  rownames(pth$DP$bootDirect$CI95Pct) = rownames(pth$DP$coefMatrix)
  colnames(pth$DP$bootDirect$CI95Pct) = colnames(pth$DP$coefMatrix)
  #pth$DP$bootDirect$p = array(0, dim(pth$DP$coefMatrix))
  pth$DP$bootDirect$p = pth$DP$coefMatrix
  pth$DP$bootDirect$p[,] = 0

  for(i in 1:pth$init$nVars){
    for(j in 1:pth$init$nVars){
      if (pth$DP$connectionMatrix[i,j]){
      thisSortedVec = sort(dir_paths[i,j,])
      pth$DP$bootDirect$CI95Pct[[i,j]] = c(thisSortedVec[idxCI[1]],thisSortedVec[idxCI[2]])
      dist_prop<- as.numeric(pth$DP$coefMatrix[i,j]>thisSortedVec)
      p_h = 2*abs(.5 - sum(dist_prop)/length(dist_prop))
      p = sapply(p_h, function(x) min(c(x, 1-x)))
      if (p==0) {
        p = 1/pth$init$nBootReps
      }
      pth$DP$bootDirect$p[i,j] = p
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
      pth$IP[[i]]$p = 1/pth$init$nBootReps
    }
    #pth$IP[[i]]$dist = ind_paths[[i]]
    pth$IP[[i]]$CI95pct = sortedDists[idxCI]
  }

  pth$DP$bootDirect$dist = NULL

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
  pth = dSepTest(pth)
  
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

#pathRes<-pathAnalysis(paths, DF, covs, RFX, nBootReps = 100)

##TO DO: 

# check for crazy and non-semantic input
# allow for lists of covariates and glmerText, in case each direct path requires something special. 
# thoroughly check other test cases, including glmerText
# is there a way to suppress the matrix of bootstrapped values from appearing when the pathRes variable is typed?  (e.g. slots?) If not, consider not storing the bootstrapped values...