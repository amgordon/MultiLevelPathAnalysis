#create matrix of direct connections
makeConnectionMatrix <- function(paths, DF){
  
  varNames = attributes(DF)$names
  nVars = length(varNames)
  connectionMatrix = matrix(0, nVars, nVars)
  rownames(connectionMatrix) = varNames
  colnames(connectionMatrix) = varNames
  
  for( i in 1:length(paths) ){
    varNames = attributes(DF)$names
    thisPath<-paths[i]
    theseVars<-strsplit(thisPath, "->")[[1]]
    V1<-which(theseVars[1] == varNames)
    V2<-which(theseVars[2] == varNames)
    connectionMatrix[V1,V2] = 1
  }
  return(connectionMatrix)
}

# create matrix of direct path coefficients	
estDirectCoefficients <- function(connectionMatrix, DF){
  coefMatrix = connectionMatrix
  coefMatrix[,] = 0
  varNames = colnames(connectionMatrix)
  for( i in 1:ncol(connectionMatrix) ){
    
    thisCol<-connectionMatrix[,i]
    
    if(sum(thisCol)>0){
      IVs<-varNames[thisCol==1]
      DV = varNames[i]
      
      IVsCombined	<- paste( c(IVs) , collapse = "+" )
      thisFormula <- as.formula( paste( DV, IVsCombined , sep = "~" ) )
      theseCoef <- coefficients( summary( lm(thisFormula, DF ) ) )[,"Estimate"]
      
      coefMatrix[which(thisCol==1),i] = theseCoef[2:length(theseCoef)]
    }
  }
  return(coefMatrix)
}


# recursively find all the indirect paths
findIndirectPaths <- function(coefMatrix, IPs, varsToSearch){
  
  allVars<-rownames(coefMatrix)
  
  for(v in varsToSearch){	
    
    i = which(allVars==v)
    thisRow<-coefMatrix[i,]
    thisVar<-rownames(coefMatrix)[i]
    IPs$thisIP = c(IPs$thisIP, thisVar) 	
    
    if (length(IPs$thisIP) > 2){
      
      IPs$ix <- IPs$ix + 1 # update path index
      IPs$IPlist[[IPs$ix]]<-IPs$thisIP # include this path in the list 
    }
    
    if(sum(abs(thisRow))>0){
      newVarsToSearch = allVars[which(thisRow!=0)]
      newCoefMatrix = coefMatrix
      newCoefMatrix[i,] = 0
      IPs = findIndirectPaths(coefMatrix, IPs, newVarsToSearch)
    }
    IPs$thisIP = IPs$thisIP[-length(IPs$thisIP)]
  }
  return(IPs)
}


# compute all indirect path coefficients
indirectPathCoefs <- function(coefMatrix, allIndirectPathDes){
  varNames = colnames(coefMatrix)
  allIndirectPathCoef = list()
  for( i in 1:length(allIndirectPathDes) ){
    thisPath = allIndirectPathDes[[i]]
    thisIndCoef = 1
    for( j in 1:(length(thisPath)-1) ){
      
      V1<-which(thisPath[j] == varNames)
      V2<-which(thisPath[j+1] == varNames)
      
      thisDirCoef = coefMatrix[V1,V2];
      thisIndCoef = thisIndCoef * thisDirCoef      
    }
    allIndirectPathCoef[[i]] = thisIndCoef
  }
  return(allIndirectPathCoef)
}	

# The main function.  Given some paths and a data frame, run a path analysis
pathAnalysis <- function(paths, DF)	{
  
  conMatrix = makeConnectionMatrix(paths,DF)
  coefMatrix = estDirectCoefficients(conMatrix,DF)
  
  IPs = list()
  IPs$ix = 0
  IPs$IPlist = list()
  IPs$thisIP = NULL
  
  IPs<-findIndirectPaths(coefMatrix, IPs, attributes(DF)$names)
  
  #trim off repeat paths (though I don't think these should exist)
  allIndirectPathDes<-unique(IPs$IPlist)
  
  indCoefs<-indirectPathCoefs(coefMatrix, allIndirectPathDes)
  
  res = list()
  res$coefMatrix = coefMatrix
  res$allIndirectPathDes = allIndirectPathDes
  res$indCoefs = indCoefs
  return(res)
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

##TO DO: pass around an object with things like varnams and nVars already included in it
# better connect the ind coefficients and their descriptions
# check for crazy input, non-semantic input, and cyclic connections (which would invalidate this procedure)
#implement bootstrapping for significance testing of ind coefficients
#make it work for lmer