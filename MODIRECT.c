#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <limits.h>

//#include "global.h"
#include "paretofiltering.h"
#include "myarray.h"
#include "wrapbbob.h"
#include "myrandom.h"


// define

#define MIN(a,b) (((a)<(b))?(a):(b))
#define MAX(a,b) (((a)>(b))?(a):(b))

// to sort the distances
int cmp(const void *v1, const void *v2) {
  const double d1 = **(const double**)v1;
  const double d2 = **(const double**)v2;
  
  return (d1 < d2)? -1:(d1 > d2);
}
// returns hmax
size_t max_depth(size_t minLevel, size_t nVars) {
	return minLevel + nVars * nVars + 20;
}



// MODIRECT algorithm
void MODIRECT(size_t nVars, size_t nObjs, size_t maxFuncEvals, size_t *func_id, size_t *inst_id, double loBound, double upBound, double pDepth) {
			   
	
	// MOSOO parameters
	size_t totalNumPoints = maxFuncEvals;
	
	char isVerbose = 0;
	//char isMaxDepth = 0;
	size_t  curDim;
	size_t splitDim;
	size_t  curLevel;
	//size_t  startIdx = 0; // for efficient loop over the variables
	double curScale = 0.;
	double offset = 0.;
	double deltaBound = upBound - loBound;
	double tempDist=0.; // to store distance between obje vectors temp
	// MODIRECT data-structures (pointers of pointers)
	struct array1Dp pPoints;
	struct array1Dp pObjs;
	double *objVectors;
	double *pointVectors;
	double *distVector;
	double **distPtr;// pointer to distVector to compute the offset
	int *pointLevels;
	int *dimSeq;
	// a set of flags for the sample points
	// TO DO : for even K we can replace its item
	bool *isExpandable ; // is it expandable currently
	//bool *isV; // is it part of the 

	
	// allocate
	objVectors = malloc(sizeof(double) * nObjs * totalNumPoints);
	pointVectors = malloc(sizeof(double) * nVars * totalNumPoints);
	distVector = malloc(sizeof(double) * nVars);
	distPtr = malloc(sizeof(double*) * nVars);
	isExpandable = malloc(sizeof(bool) * totalNumPoints);
	dimSeq = malloc(sizeof(int) * nVars * totalNumPoints);
	pointLevels  = malloc(sizeof(int) * totalNumPoints);	
	array1Dp_construct(&pPoints, totalNumPoints);
	array1Dp_construct(&pObjs, totalNumPoints);
	
	// initialization 
	size_t numPoints = 1;
	size_t pointIdx = 0;
	//size_t curPoint;
	// to keep track of the present levels (depths) in the tree
	size_t minLevel = 1;
	size_t nExpMinLvl = 0;
	size_t maxLevel = 0;
	size_t iter = 0;
	// point-wise initialization
	for(size_t p=0; p < totalNumPoints; p++) {
	  *array1Dp_element(&pPoints, p) = &(pointVectors[p * nVars]);
	  *array1Dp_element(&pObjs, p) = &(objVectors[p * nObjs]);
	  isExpandable[p] = false;
	}
	
	isExpandable[0] = true;
	pointLevels[0] = 0;
	// sample the center
	for(size_t j =0; j<nVars; j++) {
	  pointVectors[pointIdx * nVars + j] = 0.5 * deltaBound + loBound;
	  dimSeq[pointIdx * nVars + j] = j;
	}
	mobbob_eval_testLogging_mosoo(pObjs.elements, pPoints.elements, func_id, inst_id, nVars, numPoints, nObjs, pointIdx);
	
	//curLevel=1;	   
	// Core of the algorithm
	while (true) {
	  iter = iter + 1;
    // 1. Expand Q nodes sequentially and sample their representative sets
    // a) Expand
    pointIdx = numPoints;

    if (isVerbose) printf("=== A NEW ROUND===\n");
    for(size_t p=0; p < numPoints; p++) if (isExpandable[p]) {
      //curPoint = p;
		  curLevel = pointLevels[p];
		  curDim = (curLevel) % nVars;
		  curScale = deltaBound / pow(3, (curLevel) / nVars + 1);
		  // online update minLevel
		  if (pointLevels[p] == minLevel) {
		    nExpMinLvl = nExpMinLvl + 1;
		    if (nExpMinLvl == (size_t) pow(3, (minLevel))) minLevel = minLevel + 1;
		  }
		  // update the parameters of p then settle its kids
		  isExpandable[p] = false;
		  pointLevels[p] = pointLevels[p] + (nVars - curDim);
		  // generate the new nodes
		  // expand depends on the number of dimension the sequence
		  for (size_t k= curDim; k < nVars; k++){
			  splitDim = dimSeq[p * nVars + k];
			  offset = curScale ;
			  // one +
		    for(size_t j =0; j<nVars; j++)  {
		      if (j == splitDim) pointVectors[pointIdx * nVars + j] = pointVectors[p * nVars + j] + offset;
		      else pointVectors[pointIdx * nVars + j] = pointVectors[p * nVars + j];
		    }
		    // update point-wise info
		    pointLevels[pointIdx] = curLevel + (k - curDim + 1);
			  isExpandable[pointIdx] = false;
			  // force exit if you run out of budget (but sample before exit)
			  pointIdx = pointIdx + 1;
		    if (pointIdx == totalNumPoints) goto eval_children;	
			  // one -
		    for(size_t j =0; j<nVars; j++)  {
		      if (j == splitDim) pointVectors[pointIdx * nVars + j] = pointVectors[p * nVars + j] - offset;
		      else pointVectors[pointIdx * nVars + j] = pointVectors[p * nVars + j];
		    }
		    // update point-wise info
		    pointLevels[pointIdx] = curLevel + (k - curDim + 1);
			  isExpandable[pointIdx] = false;
			  // force exit if you run out of budget (but sample before exit)
			  pointIdx = pointIdx + 1;
		    if (pointIdx == totalNumPoints) goto eval_children;	
		   }
		
		  // b) Sample and update their levels and dim sequence
	    //printf("quick test: start:%zu, end:%zu\n", pointIdx, numPoints + pointIdx);
	    eval_children:
	      mobbob_eval_testLogging_mosoo(pObjs.elements, pPoints.elements, func_id, inst_id, nVars, pointIdx - numPoints, nObjs, numPoints);	
	    // Display some info
      if (isVerbose){
        printf("Exp node. %zu: (", p);
        for (size_t j=0; j < nVars- 1; j ++) printf("%f,", pointVectors[p * nVars +j]);
        printf("%f), depth=%zu, splitDim=%zu, Obj. vector:(", pointVectors[p * nVars + nVars  - 1], curLevel, curDim);
        for (size_t j=0; j < nObjs- 1; j ++) printf("%f,", objVectors[p * nObjs +j]);
        printf("%f).\n", objVectors[p * nObjs + nObjs - 1]);
       }
	    // check for termination  
		  if (pointIdx == totalNumPoints) break;	  
	    // computes the distance
	    for (size_t k= curDim; k < nVars; k++){
	      tempDist = 0.;
	      distVector[k] = 0.;
	      for(size_t j=0; j < nObjs; j++) {
	        tempDist = tempDist + abs(objVectors[p * nObjs +j]-objVectors[(numPoints + 2 * (k - curDim)) * nObjs +j]);
	        distVector[k]= distVector[k] + abs(objVectors[p * nObjs +j]-objVectors[(numPoints + 2 * (k - curDim) + 1) * nObjs +j]);
	      }
	      distVector[k] = 1.0/ MIN(tempDist, distVector[k]);
	      // re-intialize the distPtr
	      distPtr[k] = &distVector[k];
	    }
	    // sort the distances
	    qsort(distPtr+curDim, nVars - curDim, sizeof(double*), cmp);
	    // update levels
	    for (size_t k= curDim; k < nVars; k++){
	    
	      pointLevels[numPoints + 2 * (k - curDim)]= curLevel + 1 + (*(distPtr+k)- &distVector[curDim]);
	      pointLevels[numPoints + 2 * (k - curDim) + 1]= curLevel + 1 + (*(distPtr+k)- &distVector[curDim]);
	      
	    }
	    // update dimension sequence
	    for (size_t q= numPoints; q < pointIdx; q++){
	      // update the maxlevel
		    if (maxLevel < pointLevels[q]) maxLevel = pointLevels[q];
	      for (size_t k= curDim; k < nVars; k++){
		      dimSeq[q  * nVars + k]= dimSeq[p* nVars +(*(distPtr+k)- &distVector[curDim])];
	      }   
	    }
	    // Display some info about the child of the expanded nodes
      if (isVerbose){
        for (size_t q= numPoints; q < pointIdx; q++){
          printf("\t Child node. %zu: (", q- numPoints);
          for (size_t j=0; j < nVars- 1; j ++) printf("%f,", pointVectors[q * nVars +j]);
          printf("%f), depth=%zu, Obj. vector:(", pointVectors[q * nVars + nVars  - 1], pointLevels[q]);
          for (size_t j=0; j < nObjs- 1; j ++) printf("%f,", objVectors[q * nObjs +j]);
          printf("%f). DimSeq:[", objVectors[q * nObjs + nObjs - 1]);
          for (size_t j=0; j < nVars- 1; j ++) printf("%zu,", dimSeq[q * nVars +j]);
          printf("%zu).\n", dimSeq[q * nVars + nVars - 1]);
        }
       }
	    

      numPoints = pointIdx;

    }
	
    
	
	if (pointIdx == totalNumPoints) break;
  
    
    
    
    // 2. Identify potentially optimal nodes at this iteration
	// a. update isV to only those nodes of depth <= h_max
	for (size_t p = 0; p < numPoints; p++) if (pointLevels[p] <= max_depth(minLevel, nVars)) isExpandable[p] = true;
    // b. Filter the nodes of depth less than h_max
    nd_direct(isExpandable, objVectors, pointLevels, numPoints, nObjs);
    //for(size_t p = startIdx; p < numPoints; p++) if (pointLevels[p] == curLevel) isExpandable[p] = isV[p];
    
    // monitor the sampling
    if (isVerbose) sleep(3);
				
	}
			 
	// free memory:
	free(isExpandable);
	free(pointLevels);
	free(objVectors);
	free(pointVectors);
	free(dimSeq);
	array1Dp_destruct(&pPoints);
	array1Dp_destruct(&pObjs);  		   
			   
}
