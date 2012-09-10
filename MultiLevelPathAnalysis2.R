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
	

 
findIndirectPaths <- function(coefMatrix, indirectPaths, varsToSearch){
ix_h = 0 
allVars<-rownames(coefMatrix)
 	
for(v in varsToSearch){	
	i = which(allVars==v)
	thisRow<-coefMatrix[i,]
	thisVar<-rownames(coefMatrix)[i]
	indirectPaths = c(indirectPaths, thisVar) 	
	
	if (length(indirectPaths) > 2){
		#print(indirectPaths)
		ix_h  = ix + 1
		ix <<- ix_h # path index
		
		#print(ix)
		indirectPathsSet[[ix]]<<-indirectPaths # include this path 
		#print(indirectPathsSet)
	}
 	
	if(sum(abs(thisRow))>0){
		newVarsToSearch = allVars[which(thisRow!=0)]
		newCoefMatrix = coefMatrix
		newCoefMatrix[i,] = 0
		findIndirectPaths(coefMatrix, indirectPaths, newVarsToSearch)
	}
	indirectPaths = indirectPaths[-length(indirectPaths)]
}
 	return(indirectPaths)
}

 
#run the function
#findIndirectPaths(coefMatrix, indirectPaths, varNames)



# compute all Indirect Path Coefficients
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


pathAnalysis <- function(paths, DF)	{

	indirectPathsSet <<- list()
	ix <<- 0	
	
	conMatrix = makeConnectionMatrix(paths,DF)
	coefMatrix = estDirectCoefficients(conMatrix,DF)
	findIndirectPaths(coefMatrix, NULL, attributes(DF)$names)
	
	#trim off repeat paths (though I don't think these should exist)
	allIndirectPathDes<-unique(indirectPathsSet)

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

#indirectPaths = NULL
#varNames = attributes(DF)$names
#varsToSearch = varNames

paths = c('IV->MV1', 'IV->MV2', 'MV1->DV', 'MV2->DV')
#paths = c('IV->MV1', 'MV1->MV2', 'MV2->DV')

pathRes<-pathAnalysis(paths, DF)

##TO DO: pass around an object with things like varnams and nVars already included in it
# figure out a possibly less ghetto way to use global variables
# better connect the ind coefficients and their descriptions
# check for crazy input, non-semantic input, and cyclic connections (which would invalidate this procedure)
#implement bootstrapping for significance testing of ind coefficients
#make it work for lmer