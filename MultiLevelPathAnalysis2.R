# initialize the pth object
makePathObj <- function(paths, DF, covs=NULL, RFX=NULL, intercepts = TRUE, slopes = TRUE, nBootReps = 0, glmerText=NULL){
  pth = list()
  pth$paths = paths
  pth$DF = DF
  pth$varNames = attributes(DF)$names
  pth$nVars = length(pth$varNames)
  pth$connectionMatrix = matrix()
  pth$coefMatrix = matrix()
  pth$indirectPaths = list()
  pth$indirectPathCoefs = list()
  pth$covs = covs
  pth$RFX = RFX
  pth$intercepts = intercepts
  pth$slopes = slopes
  pth$glmerText = glmerText
  pth$nBootReps = nBootReps
  return (pth)
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
  for( i in 1:ncol(pth$connectionMatrix) ){    
    thisCol<-pth$connectionMatrix[,i]
    
    if(sum(thisCol)>0){
      IVs<-pth$varNames[thisCol==1]
      DV = pth$varNames[i]      
      IVsCombined	<- paste( c(IVs) , collapse = "+" )
      formulaText <- paste( DV, IVsCombined , sep = "~" )
      
      for(j in 1:length(pth$covs)){
        formulaText <- paste(formulaText, pth$covs[j], sep = "+")
      }
      
      if(!is.null(pth$glmerText)){
        formulaText <- paste(formulaText, pth$glmerText, sep = "+")
      }
           
      for(j in 1:length(pth$RFX)){
        print(j)
        if(pth$intercepts){
          print(pth$RFX)
          thisRandEffect<-paste('(1|',pth$RFX[j],')', sep = "")
          formulaText <- paste(formulaText, thisRandEffect, sep = "+")
        }
        
        if(pth$slopes){
          print("slopes")
          for(k in 1:length(IVs)){
            thisSlope<-paste('(', IVs[k], '+0|', pth$RFX[k], ')', sep = "")
            print(thisSlope)
            formulaText <- paste(formulaText, thisSlope, sep = "+")
          }
        }
      }
              
      print("text is")
      print(formulaText)
      
      thisFormula <- as.formula(formulaText)
      
      if(!is.null(pth$glmerText)){    
        thisMod<-glmer(thisFormula, pth$DF) # run a mixed glm
        theseCoef<-attr(thisMod, "fixef")
      } else {
        theseCoef <- coefficients(summary(lm(thisFormula, pth$DF)))[,"Estimate"] #otherwise run a normal linear model
      }
      colsToReplace = which(thisCol==1)
      coefMatrix[colsToReplace,i] = theseCoef[2:(1+length(colsToReplace))]
    }
  }
  pth$coefMatrix = coefMatrix
  return(pth)
}

# a wrapper to call the function "findIndirectPaths"
findIndirectPathsWrapper <- function(pth) {
  pth$ix<-0
  pth = findIndirectPaths(pth, pth$varNames)
  pth$thisIP<-NULL
  pth$ix<-NULL
  
  #trim off repeat paths (unclear if this is still needed)
  pth$indirectPaths<-unique(pth$indirectPaths)
  
  return(pth)
}

# recursively find all the indirect paths
findIndirectPaths <- function(pth, varsToSearch){  
  
  for(v in varsToSearch){	
    
    i = which(pth$varNames==v)
    thisRow<-pth$coefMatrix[i,]
    thisVar<-pth$varNames[i]
    pth$thisIP = c(pth$thisIP, thisVar) 	
    
    if (length(pth$thisIP) > 2){      
      pth$ix <- pth$ix + 1 # update path index
      pth$indirectPaths[[pth$ix]]<-pth$thisIP # include this path in the list 
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
  for( i in 1:length(pth$indirectPaths) ){
    thisPath = pth$indirectPaths[[i]]
    thisIndCoef = 1
    for( j in 1:(length(thisPath)-1) ){
      
      V1<-which(thisPath[j] == pth$varNames)
      V2<-which(thisPath[j+1] == pth$varNames)
      
      thisDirCoef = pth$coefMatrix[V1,V2];
      thisIndCoef = thisIndCoef * thisDirCoef      
    }
    allIndirectPathCoef[[i]] = thisIndCoef
  }
  pth$indirectPathCoefs = allIndirectPathCoef
  return(pth)
}	

bootStrapCoefs <- function(pth){
  dir_paths = array(0, c(dim(pth$coefMatrix), pth$nBootReps))
  ind_paths = NULL
  for(i in 1:pth$nBootReps){
    boot_DF 	<- pth$DF[sample( 1:nrow(pth$DF), replace=TRUE),]
    boot_pth = pth
    boot_pth$DF = boot_DF
    boot_pth = directPathCoeffs(boot_pth)
    boot_pth = indirectPathCoefs(boot_pth)
    ind_paths 	<- as.data.frame(rbind(ind_paths, unlist(boot_pth$indirectPathCoefs)))
    dir_paths[,,i]<-boot_pth$coefMatrix
  }
  idxCI = c(round(.025*nrow(ind_paths)+.5), round(.975*nrow(ind_paths)+.5))
  
  pth$bootDirect$dist = dir_paths
  pth$bootDirect$CI95Pct = array(list(NULL), dim(pth$coefMatrix))
  pth$bootDirect$p = array(-1, dim(pth$coefMatrix))
  for(i in 1:pth$nVars){
    for(j in 1:pth$nVars){
      thisSortedVec = sort(dir_paths[i,j,])
      pth$bootDirect$CI95Pct[[i,j]] = list(thisSortedVec[idxCI[1]],thisSortedVec[idxCI[2]])
      dist_prop<- as.numeric(pth$coefMatrix[i,j]>thisSortedVec)
      p_h = 2*abs(.5 - sum(dist_prop)/length(dist_prop))
      p = sapply(p_h, function(x) min(c(x, 1-x)))
      if (p==0) {
        p = 1/pth$nBootReps
      }
      pth$bootDirect$p[i,j] = p
    }
  }
  
  ind = list()
  ind$sorted_dists = sapply(ind_paths, sort)
  ind$dist_prop<-sapply(as.matrix(ind_paths), function(x) as.numeric(x>unlist(pth$indirectPathCoefs)))
  ind$p_h = 2*abs(.5 - rowSums(ind$dist_prop)/ncol(ind$dist_prop))
  
  p = sapply(ind$p_h, function(x) min(c(x, 1-x)))
  if (p>0) {
    ind$p = p
  }else {
    ind$p = 1/pth$nBootReps
  }

pth$bootIndirect$dist = ind_paths
pth$bootIndirect$CI95pct = ind$sorted_dists[idxCI,]
pth$bootIndirect$p = ind$p

return(pth)
}


# the main function.  given some paths and a data frame, run a path analysis
pathAnalysis <- function(paths, DF, covs=NULL, RFX=NULL, intercepts = TRUE, slopes = TRUE, nBootReps = 0, glmerText=NULL){
  
  pth = makePathObj(paths, DF, covs, RFX, intercepts, slopes, nBootReps, glmerText)
  
  pth = makeConnectionMatrix(pth)
  pth = directPathCoeffs(pth)
  pth = findIndirectPathsWrapper(pth)  
  pth = indirectPathCoefs(pth)
  
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
# better connect the ind coefficients and their descriptions
# report full equations
# limit coef matrix etc. only to variables in the path
# check for crazy input, non-semantic input, and cyclic connections (which would invalidate this procedure)
# implement bootstrapping for significance testing of direct coefficients
# allow for lists of covariates and glmerText, in case each direct path requires something special. 
# thoroughly check other test cases
# allow entry for slopes = TRUE and intercepts = TRUE
# label variables for bootIndirect and bootDirect