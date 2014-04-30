# MultiLevel Path Analysis
# Alan Gordon, Stanford University Dept. of Psychology, September 2012
# The d-sep test of goodness-of-fit for multilevel models is taken from Shipley, B. "Confirmatory path analysis in a generalized multilevel context," Ecology 90:363-368, 2009.

# the main function is:
# pathAnalysis(paths, DF,  covs=NULL, RFX=NULL, intercepts = TRUE, slopes = TRUE, nBootReps = 0, dichotVars = NULL, stdCoeffs = FALSE)

# paths: A list of unidirectional connections, each in the form 'var1->var2,' indicating that var1 acts on var2.
# DF: A data frame containing data from all variables, covariates, and random effects.
# covs: A list of covariates variable names.  These variables are controlled for in the path equations, but are not included in the paths.
# RFX: A list of random effect variable names.
# intercepts: TRUE = model all random intercepts.  FALSE = do not model random intercepts.
# slopes: TRUE = model all random slopes.  FALSE = do not model random slopes.
# nBootReps: The number of bootstrapping iterations used to determine significance for indirect paths.
# dichotVars: A list of all variables that are dichotomous.  When these varialbes are dependent variables, logistic regressions will be used.
# stdCoeffs: do we want to standardize the coefficients of logistic models, as in the fully standardized regression coefficients mentioned in Scott Menard (2004) Six Approaches to Calculating Standardized Logistic Regression Coefficients, The American Statistician, 58:3, 218-223
#
#----------------Usage Cases----------------
#

# 1) Using a bad model to predicting tree death with Shipley.dat 
# DF = read.delim("Shipley.dat", sep=" ")
# paths1 = c('lat->DD', 'DD->Date', 'DD->Growth', 'Growth->Live')
# pathRes1 = pathAnalysis(paths1, DF, dichotVars = 'Live', RFX=c('site', 'tree'), intercepts = TRUE, slopes = FALSE)
# pathRes1$MF #check model fit. p = .0012, therefore the data significantly diverges from the model
# 
# 2) Using a better model to predicting tree death with Shipley.dat 
# DF = read.delim("Shipley.dat", sep=" ")
# paths2 = c('lat->DD', 'DD->Date', 'Date->Growth', 'Growth->Live')
# pathRes2 = pathAnalysis(paths2, DF, dichotVars = 'Live', RFX=c('site', 'tree'), intercepts = TRUE, slopes = FALSE)
# pathRes2$MF #check model fit. p = .588, therefore we can retain the model
# 
# 3) Use randomly generated data to test significance of indirect paths
# n  <- 100
# DV  <- rnorm( n )
# MV2	<- DV + rnorm( n , 2 )
# MV1	<- MV2 + rnorm( n , 2 )
# c1  <- MV2 + rnorm( n , 2 )
# IV	<- MV1 + rnorm( n , 2 )
# subs <- round(2*runif(100))
# DF	<- data.frame( DV , MV2, MV1, IV, c1, subs)
# paths = c('IV->MV1', 'MV1->MV2', 'MV2->DV')
# covs = ("c1")
# RFX = "subs"
# 
# pathRes3<-pathAnalysis(paths, DF, covs, RFX, nBootReps = 2000) # Run the path Analysis
# pathRes3$MF  # Check the model fit, derived with a d-sep test
# pathRes3$DP # List the coefficients, p-values, and models for direct connections between variables
# pathRes3$IP # List all indirect paths, their coefficients, bootstrapped p-vals, and 95% confidence intervals


require("lme4")
require("abind")

# initialize the pth object
makePathObj <- function(paths, DF, covs=NULL, RFX=NULL, intercepts = TRUE, slopes = TRUE, nBootReps = 0, dichotVars = NULL, stdCoeffs = FALSE){
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
  pth$DP$modelList = list()
  pth$init$nBootReps = nBootReps
  pth$init$dichotVars= dichotVars
  pth$init$stdCoeffs = stdCoeffs
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
  
  #subset the data frame to only include relevant variables
  allVarNames = c(varNames, pth$init$RFX, pth$init$covs)
  pth$init$DF = pth$init$DF[,allVarNames]
  
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
    connectionMatrix[V1,V2] = TRUE #turn this matrix TRUE if there is a connection between V1 and V2
  }
  pth$DP$connectionMatrix = connectionMatrix
  return(pth)
}

# run a regression
runRegression <- function(pth, DV, IVs, covs){
  
  IVsCombined  <- paste( c(IVs) , collapse = "+" )
  formulaText <- paste( DV, IVsCombined , sep = "~" )  # basic regression formula
  
  # add covariates
  if (!is.null(covs)){
    for(j in 1:length(covs)){
      formulaText <- paste(formulaText, covs[j], sep = "+")
    }
  }
  
  for(j in 1:length(pth$init$RFX)){    
    # add random intercept terms
    if(pth$init$intercepts){
      thisRandEffect<-paste('(1|',pth$init$RFX[j],')', sep = "")
      formulaText <- paste(formulaText, thisRandEffect, sep = "+")
    }
    #add random slope terms
    if(pth$init$slopes){
      for(k in 1:length(IVs)){
        thisSlope<-paste('(', IVs[k], '+0|', pth$init$RFX[j], ')', sep = "")
        formulaText <- paste(formulaText, thisSlope, sep = "+")       
      }
    }
  }
  
  thisFormula <- as.formula(formulaText)
  thisDF = pth$init$DF
  
  # if the DV is dichotomous, specify that a binomial family will be used
  if(is.element(DV, pth$init$dichotVars)){
    familyText <-  "binomial"
    modText <-"glmer"
    lme4Args=list(thisFormula, data = thisDF, family=familyText)
  }else{
    modText <- "lmer"
    familyText <- "gaussian"
    lme4Args=list(formula=thisFormula, data = thisDF, REML=FALSE)
    lme4Txt = paste(formulaText, "data=thisDF", "REML=FALSE", sep = ', ')
  }
  
  # if random slope or intercepts are included, use glmer.  Otherwise, use glm.
  if(pth$init$intercepts | pth$init$slopes){    
        
    thisMod<-do.call(modText, args=lme4Args)
    theseCoef<-attr(thisMod, "fixef")
  } else {
    thisMod<-glm(thisFormula, pth$init$DF, family=familyText)
    theseCoef <- coefficients(summary(thisMod))[,"Estimate"] #otherwise run a normal linear model
  }

  if(is.element(DV, pth$init$dichotVars) & pth$init$stdCoeffs){
    for(k in 2:length(theseCoef)){
      sdIV = sd(getElement(pathRes$init$DF, names(theseCoef[k])))
      sd_logit_Yhat = sd(logit(fitted(thisMod)))
      Rsquared = summary(lm(thisMod@y ~ fitted(thisMod)))$r.squared
      theseCoef[k] = theseCoef[k] * sdIV * sqrt(Rsquared) / sd_logit_Yhat
    }
  }
  
  output = list()
  output$mod = thisMod
  output$coef = theseCoef
  return(output)
}


# create matrix of direct path (DP) and simple effect (SE) coefficients
directPathCoeffs <- function(pth){
  coefMatrix = pth$DP$connectionMatrix
  coefMatrix[,] = 0
  coefMatrixSimple = coefMatrix
  p = pth$DP$connectionMatrix
  p[,] = -1
  pSimple = p
  modelList<-list()
  modelListSimple = modelList
  eq = 0
  eqS = 0
  
  for( i in 1:ncol(pth$DP$connectionMatrix) ){    
    thisCol<-pth$DP$connectionMatrix[,i]
    if(any(thisCol)){
      IVs<-pth$init$varNames[thisCol]
      DV = pth$init$varNames[i]      
      thisReg = runRegression(pth, DV, IVs, pth$init$covs)
      eq = eq+1

      colsToReplace = which(thisCol)
      coefMatrix[colsToReplace,i] = fixef(thisReg$mod)[2:(1+length(colsToReplace))]
      p[colsToReplace,i] = LRT(thisReg$mod, names(colsToReplace))
      
      # remove the data frame from being reported in the call slot
      thisMod = thisReg$mod
      thisMod@call = thisMod@call[1:2]
      modelList[[eq]] = thisMod
      
      for( j in which(thisCol==1) ){
        eqS = eqS+1
        IV<-pth$init$varNames[j]
        DV = pth$init$varNames[i]      
        thisReg = runRegression(pth, DV, IV, pth$init$covs)

        colsToReplace = j

        coefMatrixSimple[colsToReplace,i] = fixef(thisReg$mod)[2]

        pSimple[colsToReplace,i] = LRT(thisReg$mod, 1)
        
        # remove the data frame from being reported in the call slot
        thisMod = thisReg$mod
        thisMod@call = thisMod@call[1:2]
        modelListSimple[[eqS]] = thisMod
      }
    }
  }
  pth$DP$coefMatrix = coefMatrix
  pth$DP$pMatrix = p
  pth$DP$modelList = modelList
  
  pth$SE$coefMatrix = coefMatrixSimple
  pth$SE$pMatrix = pSimple
  pth$SE$modelList = modelListSimple
  return(pth)
}

# use a log-likelihood ratio test to determine the significance of a coefficient in a mer object
# idx can be specified as a list of variable names, or a numerical list of variable indices in the model.
LRT <- function(theModel, idx=1){
  
  modFull = theModel
  pVal = -1;

  modelTerms_h<-terms(theModel)  
  modelTerms = attr(modelTerms_h, "term.labels")
  
  for (i in 1:length(idx)){
    if (is.character(idx)) {
      updateStr<-paste(".~.-",idx[i])
    } else {
      thisIdx = idx[i]
      updateStr<-paste(".~.-",as.character(modelTerms[thisIdx]))
    }
    modNested = update(theModel, updateStr)  
    llFull = logLik(modFull)
    llNested = logLik(modNested)
    
    chsqval<-2*(llFull[1] - llNested[1])
    dfDiff<- attr(llFull,"df") - attr(llNested,"df") 
    
    pVal[i]<-1-pchisq(chsqval,dfDiff)
  }
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
        
        # do any indirect paths contain both V1 and V2?
        indPathsContainsVars = mapply(function(x) is.element(V1,x) & is.element(V2,x), indPaths)
        
        # if a path exists that contains both V1 and V2

        if (any(indPathsContainsVars)){
          
          # if a path exists that contains both V1 and V2, the variable 
          # located earlier in the path is the IV and the variable 
          # located later is the DV.
          thisPath = unlist(indPaths[which(indPathsContainsVars)[1]])
          match1 = match(V1,thisPath)
          match2 = match(V2,thisPath)
          if(match1 > match2){
            thisReg = runRegression(pth, V1, V2, precedingVars)
          }else{
            thisReg = runRegression(pth, V2, V1, precedingVars)
          }
          pVals[idxBf] = LRT(thisReg$mod)
        }  else {
          # if no indirect path contains both variables, do the regression 
          # both with V1 as the DV and V2 as the DV, and take the mean 
          # p-value across the regressions.
          thisReg1 = runRegression(pth, V1, V2, precedingVars)
          thisReg2 = runRegression(pth, V2, V1, precedingVars)
          pVals[idxBf] = .5*(LRT(thisReg1$mod) + LRT(thisReg2$mod))
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
    
    # if a variable appears twice in an indirect path, the path is cyclical,
    # and cannot be solved with this algorithm
    if (length(unique(pth$thisIP)) != length(pth$thisIP)) {
      badPath = paste(pth$thisIP, collapse = "->")
      errText = paste("Cyclic path found: ", badPath, ". Path connections must be cyclic.", sep="")
      stop(errText)
    }
    
    if (length(pth$thisIP) > 2){      
      pth$ix <- pth$ix + 1 # update path index
      pth$IP[[pth$ix]] = list()
      pth$IP[[pth$ix]]$path = pth$thisIP # update the path list
    }
    
    if(sum(abs(thisRow))>0){
      newVarsToSearch = pth$init$varNames[which(thisRow!=0)]
      newCoefMatrix = pth$DP$coefMatrix
      newCoefMatrix[i,] = 0
      pth = findIndirectPaths(pth, newVarsToSearch) # recursively call findIndirectPaths
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

bootIterFun <-function(pth){
  
  boot_DF   <- pth$init$DF[sample( 1:nrow(pth$init$DF), replace=TRUE),]
  boot_pth = pth
  boot_pth$init$DF = boot_DF
  
  boot_pth = directPathCoeffs(boot_pth)
  boot_pth = indirectPathCoefs(boot_pth)
  theseIndCoefs = sapply(boot_pth$IP, function(x) getElement(x,"coef"))
  
  out = list()
  out$ind_paths <- theseIndCoefs
  out$dir_paths<-boot_pth$DP$coefMatrix
  out$simple_paths<-boot_pth$SE$coefMatrix
  
  return(out)
}

# Find bootstrap confidence intervals and p-values
bootStrapCoefs <- function(pth){
  dir_paths = array(0, c(dim(pth$DP$coefMatrix), pth$init$nBootReps))
  simple_paths = array(0, c(dim(pth$SE$coefMatrix), pth$init$nBootReps))
  ind_paths = NULL
  
  cat("bootstrapping", pth$init$nBootReps, "iterations: ");    flush.console()
  for(i in 1:pth$init$nBootReps){    
    #if (i %% round(pth$init$nBootReps/10)==0){
    #  cat(paste(round(100*i/pth$init$nBootRep), "% ", sep = ""));    flush.console()
    #}
    
    boot_DF 	<- pth$init$DF[sample( 1:nrow(pth$init$DF), replace=TRUE),]
    boot_pth = pth
    boot_pth$init$DF = boot_DF
    
    boot_pth = directPathCoeffs(boot_pth)
    boot_pth = indirectPathCoefs(boot_pth)
    theseIndCoefs = sapply(boot_pth$IP, function(x) getElement(x,"coef"))
    ind_paths 	<- as.data.frame(rbind(ind_paths, theseIndCoefs))
    dir_paths[,,i]<-boot_pth$DP$coefMatrix
    simple_paths[,,i]<-boot_pth$SE$coefMatrix
  }
  
  #bootRes = replicate(pth$init$nBootReps,bootIterFun(pth))
  
  #ind_paths = (bootRes["ind_paths",])
  #dir_paths = abind(bootRes["dir_paths",], along = 3)
  #simple_paths = abind(bootRes["simple_paths",], along = 3)
  
  idxCI = c(floor(.025*nrow(ind_paths)), ceiling(.975*nrow(ind_paths)))
  
  if (idxCI[1]==0){
    idxCI[1]=1
  }
  
  #The code below performs bootstrapping for testing the significance of the direct paths
  pth$DP$bootDirect$dist = dir_paths
  #pth$DP$bootDirect$CI95Pct = array(list(NULL), dim(pth$DP$coefMatrix))
  #rownames(pth$DP$bootDirect$CI95Pct) = rownames(pth$DP$coefMatrix)
  #colnames(pth$DP$bootDirect$CI95Pct) = colnames(pth$DP$coefMatrix)
  #pth$DP$bootDirect$p = pth$DP$coefMatrix
  #pth$DP$bootDirect$p[,] = NaN
  pth$DP$bootDirect$CMinusCPrime$p = pth$DP$coefMatrix
  pth$DP$bootDirect$CMinusCPrime$p[,] = NaN
  pth$DP$bootDirect$CMinusCPrime$DiffScore = pth$DP$bootDirect$CMinusCPrime$p
  
  for(i in 1:pth$init$nVars){
    for(j in 1:pth$init$nVars){
      if (pth$DP$connectionMatrix[i,j]){
        #thisSortedVec = sort(dir_paths[i,j,])
        #pth$DP$bootDirect$CI95Pct[[i,j]] = c(thisSortedVec[idxCI[1]],thisSortedVec[idxCI[2]])
        #dist_prop<- as.numeric(thisSortedVec>0)
        #p_h = 2*min(c(mean(dist_prop), 1-mean(dist_prop)))
        #
        #if (p_h==0) {
        #  p = 1/length(dist_prop)
        #}else{
        #  p = p_h
        #}
        #pth$DP$bootDirect$p[i,j] = p
        
        # C minus Cprime bootstrapped
        if (any(dir_paths[i,j,]!=simple_paths[i,j,])){
          thisVec = simple_paths[i,j,] - dir_paths[i,j,]
          dist_prop<- as.numeric(thisVec>0)
          p_h = 2*min(c(mean(dist_prop), 1-mean(dist_prop)))
          
          if (p_h==0) {
            p = 1/length(dist_prop)
          }else{
            p = p_h
          }
          pth$DP$bootDirect$CMinusCPrime$p[i,j] = p
          pth$DP$bootDirect$CMinusCPrime$DiffScore[i,j] = mean(thisVec)
        }
      }
    }
  }
  
  # bootstrap testing for indirect paths.
  for(i in 1:length(ind_paths)){
    sortedDists = sort(ind_paths[[i]])
    distProp = mean(as.numeric(sortedDists > 0))
    p_h = 2*min(c(mean(distProp), 1-mean(distProp)))
    
    if (p_h==0) {
      p = 1/length(sortedDists)
    }else{
      p = p_h
    }
    
    pth$IP[[i]]$p = p
    
    #pth$IP[[i]]$dist = ind_paths[[i]]
    pth$IP[[i]]$CI95pct = sortedDists[idxCI]
  }
  
  pth$DP$bootDirect$dist = NULL
  
  return(pth)
}


# The main function.  Run a path analysis!
pathAnalysis <- function(paths, DF,  covs=NULL, RFX=NULL, intercepts = TRUE, slopes = TRUE, nBootReps = 0, dichotVars = NULL, stdCoeffs = FALSE){
  
  pth = makePathObj(paths, DF, covs, RFX, intercepts, slopes, nBootReps, dichotVars, stdCoeffs)
  pth = getVarNames(pth)
  pth = makeConnectionMatrix(pth)
  pth = directPathCoeffs(pth)
  pth = findIndirectPathsWrapper(pth)  
  pth = indirectPathCoefs(pth)
  if(nBootReps>0){
    pth = bootStrapCoefs(pth)
  }
  pth = dSepTest(pth)
  
  return(pth)
}





##TO DO: 

# allow for path connections to be entered as a matrix
# check more thoroughly for nonsense input, use other test cases.
# allow for lists of covariates and glmerText, in case each direct path requires different model params.