MultiLevel Path Analysis

Alan Gordon, Stanford University Dept. of Psychology, September 2012

This package enables statistical path analysis in a multilevel framework.  Given a path structure and a data frame, it computes overall model fit and direct effects.  Additionally, it automatically computes all indirect effects.  Model fit is tested using a d-sep goodness of fit test.  Significance of indirect effects is tested with bootstrapping methods.  

The d-sep test of goodness-of-fit for multilevel models is taken from Shipley, B. "Confirmatory path analysis in a generalized multilevel context," Ecology 90:363-368, 2009.

Mixed linear models are solved using R package lme4: http://cran.r-project.org/web/packages/lme4/index.html

An example of how this method can be used in cognitive neuroscience can be found in Gordon, A. M. et al. (in press). "Cortical reinstatement mediates the relationship between content-specific encoding activity and subsequent recollection decisions." Cerebral Cortex 

the main function is:
pathAnalysis(paths, DF,  covs=NULL, RFX=NULL, intercepts = TRUE, slopes = TRUE, nBootReps = 0, dichotVars = NULL, stdCoeffs = FALSE)

**inputs:**
* `paths` : A list of unidirectional connections, each in the form 'var1->var2,' indicating that var1 acts on var2.
* `DF` : A data frame containing data from all variables, covariates, and random effects.
* `covs` : A list of covariates variable names.  These variables are controlled for in the path equations, but are not included in the paths.
* `RFX` : A list of random effect variable names.
* `intercepts` : TRUE = model all random intercepts.  FALSE = do not model random intercepts.
* `slopes` : TRUE = model all random slopes.  FALSE = do not model random slopes.
* `nBootReps` : The number of bootstrapping iterations used to determine significance for indirect paths.
* `dichotVars` : A list of all variables that are dichotomous.  When these varialbes are dependent variables, logistic regressions will be used.
* `stdCoeffs` : Do we want to standardize the coefficients of logistic models, as in the fully standardized regression coefficients mentioned in Scott Menard (2004) Six Approaches to Calculating 
Standardized Logistic Regression Coefficients, The American Statistician, 58:3, 218-223

**outputs:**
* `MF` : Display the overall model fit.
  * MF$basisPVals : Component p values of non-connections among variables, used to calculate model fit.
  * MF$C : Chi-squared statistic of model fit.
  * MF$df : Degrees of freedom for test of model fit.
  * MF$p : Model fit p-value.
* `DP` : List the coefficients, p-values, and models for direct connections between variables.
  * DP$connectionMatrix : Matrix of modeled connections among variables.  IVs are rows, DVs are columns.
  * DP$coefMatrix: Matrix of direct connections among variables.
  * DP$modelList: List of solved models used to generate direct connection coefficients and p-values.
  * DP$pMatrix: Matrix of pValues of relationships among variables.
* `IP` : List all indirect paths, their coefficients, bootstrapped p-vals, and 95% confidence intervals.

Usage Cases:


1) Using a bad model to predict tree death with Shipley.dat 
```
{
DF = read.delim("Shipley.dat", sep=" ")
paths1 = c('lat->DD', 'DD->Date', 'DD->Growth', 'Growth->Live')
pathRes1 = pathAnalysis(paths1, DF, dichotVars = 'Live', RFX=c('site', 'tree'), intercepts = TRUE, slopes = FALSE)
pathRes1$MF #check model fit. p = .0012, therefore the data significantly diverges from the model
}
```


2) Using a better model to predict tree death with Shipley.dat 
```
{
DF = read.delim("Shipley.dat", sep=" ")
paths2 = c('lat->DD', 'DD->Date', 'Date->Growth', 'Growth->Live')
pathRes2 = pathAnalysis(paths2, DF, dichotVars = 'Live', RFX=c('site', 'tree'), intercepts = TRUE, slopes = FALSE)
pathRes2$MF #check model fit. p = .597, therefore we can retain the model
}
```


3) Use randomly generated data to test significance of indirect paths
```
{
n  <- 100
DV  <- rnorm( n )
MV2	<- DV + rnorm( n , 2 )
MV1	<- MV2 + rnorm( n , 2 )
c1  <- MV2 + rnorm( n , 2 )
IV	<- MV1 + rnorm( n , 2 )
subs <- round(2*runif(100))
DF	<- data.frame( DV , MV2, MV1, IV, c1, subs)
paths = c('IV->MV1', 'MV1->MV2', 'MV2->DV')
covs = ("c1")
RFX = "subs"
pathRes3<-pathAnalysis(paths, DF, covs, RFX, nBootReps = 2000)
}
```
