makePathObj <- function(paths, DF){
  pth = list()
  pth$paths = paths
  pth$DF = DF
  pth$varNames = attributes(DF)$names
  pth$nVars = length(pth$varNames)
  pth$connectionMatrix = matrix()
  pth$coefMatrix = matrix()
  pth$indirectPaths = list()
  pth$indirectPathCoefs = list()
  return (pth)
}

#create matrix of direct connections
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
estDirectCoefficients <- function(pth){
  coefMatrix = pth$connectionMatrix
  coefMatrix[,] = 0
  for( i in 1:ncol(pth$connectionMatrix) ){    
    thisCol<-pth$connectionMatrix[,i]
    
    if(sum(thisCol)>0){
      IVs<-pth$varNames[thisCol==1]
      DV = pth$varNames[i]      
      IVsCombined	<- paste( c(IVs) , collapse = "+" )
      thisFormula <- as.formula( paste( DV, IVsCombined , sep = "~" ) )
      theseCoef <- coefficients( summary( lm(thisFormula, pth$DF ) ) )[,"Estimate"]
      
      coefMatrix[which(thisCol==1),i] = theseCoef[2:length(theseCoef)]
    }
  }
  pth$coefMatrix = coefMatrix
  return(pth)
}

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
      print(pth$ix)
    }
    
    if(sum(abs(thisRow))>0){
      newVarsToSearch = pth$varNames[which(thisRow!=0)]
      newCoefMatrix = coefMatrix
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

# The main function.  Given some paths and a data frame, run a path analysis
pathAnalysis <- function(paths, DF)	{
  
  pth = makePathObj(paths, DF)
  
  pth = makeConnectionMatrix(pth)
  pth = estDirectCoefficients(pth)
  pth = findIndirectPathsWrapper(pth)  
  pth = indirectPathCoefs(pth)
  
  return(pth)
}

#Usage Case
n	<- 100
DV	<- rnorm( n )
MV2	<- DV + rnorm( n , 2 )
MV1	<- MV2 + rnorm( n , 2 )
IV	<- MV1 + rnorm( n , 2 )
DF	<- data.frame( DV , MV2, MV1, IV )


paths = c('IV->MV1', 'IV->MV2', 'MV1->DV', 'MV2->DV')
#paths = c('IV->MV1', 'MV1->MV2', 'MV2->DV')

pathRes<-pathAnalysis(paths, DF)

##TO DO: 
# better connect the ind coefficients and their descriptions
# check for crazy input, non-semantic input, and cyclic connections (which would invalidate this procedure)
#implement bootstrapping for significance testing of ind coefficients
#make it work for lmer