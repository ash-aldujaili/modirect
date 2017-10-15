// paretofront returns the logical Pareto membership of a set of points
// synopsis:  frontFlag = paretofront(objMat)
// Created by Yi Cao: y.cao@cranfield.ac.uk
// for compiling type:
//   mex paretofront.c

#include <stdio.h>
#include <stdlib.h>  // memory, e.g. malloc
#include <stdbool.h> // to use the bool datatype, required C99
#include <math.h>
#include "paretofiltering.h"

#ifdef __cplusplus
extern "C" {
#endif

void paretoFiltering(struct SolutionsArchive *archive) {
    // Create the objective vectors and frontFlag of appropriate format for paretofront()
    size_t len = archive->size;
    size_t nObjs = archive->numobj;
    bool *frontFlag = (bool*) malloc(len * sizeof(bool));
    double *obj = (double*) malloc(len * nObjs * sizeof(double));
    for (size_t i=0; i < len; i++) {
        for (size_t k=0; k < nObjs; k++) {
            obj[i + k*len] = archive->active[i]->obj[k];
        }
        frontFlag[i] = false;
    }
    
    // Call the non-dominated sorting engine
    paretofront(frontFlag, obj, len, nObjs);
    
    // Mark non-dominated solutions and filter out dominated ones
    size_t s = 0;
    for (size_t i=0; i < len; i++) {
        if (frontFlag[i] == true) {
            archive->active[i]->status = 1;
            if (i != s)
                archive->active[s] = archive->active[i];
            s++;
        } else {
            archive->active[i]->status = 0; // filter out dominated solutions
        }
    }
    archive->size = s;
    
    free(obj);
    free(frontFlag);
}


void paretofront(bool *frontFlag, double *obj, unsigned nrow, unsigned ncol) {
    unsigned t, s, i, j, j1, j2;
    bool *checklist, colDominatedFlag;
    
    checklist = (bool*)malloc(nrow*sizeof(bool));
    
    for(t=0; t<nrow; t++)
        checklist[t] = true;
    for(s=0; s<nrow; s++) {
        t = s;
        if (!checklist[t])
            continue;
        checklist[t] = false;
        colDominatedFlag = true;
        for(i=t+1; i<nrow; i++) {
            if (!checklist[i])
                continue;
            checklist[i] = false;
            for (j=0,j1=i,j2=t; j<ncol; j++,j1+=nrow,j2+=nrow) {
                if (obj[j1] < obj[j2]) {
                    checklist[i] = true;
                    break;
                }
            }
            if (!checklist[i])
                continue;
            colDominatedFlag = false;
            for (j=0,j1=i,j2=t; j<ncol; j++,j1+=nrow,j2+=nrow) {
                if (obj[j1] > obj[j2]) {
                    colDominatedFlag = true;
                    break;
                }
            }
            if (!colDominatedFlag) { //swap active index continue checking
                frontFlag[t] = false;
                checklist[i] = false;
                colDominatedFlag = true;
                t = i;
            }
        }
        frontFlag[t] = colDominatedFlag;
        if (t>s) {
            for (i=s+1; i<t; i++) {
                if (!checklist[i])
                    continue;
                checklist[i] = false;
                for (j=0,j1=i,j2=t; j<ncol; j++,j1+=nrow,j2+=nrow) {
                    if (obj[j1] < obj[j2]) {
                        checklist[i] = true;
                        break;
                    }
                }
            }
        }
    }
    free(checklist); 
}


void selectiveparetofront(bool *frontFlag, double *obj, size_t nPoints, size_t nObjs) {
    
    //bool *checklist, colDominatedFlag;
    size_t nPosDiff;
    size_t nNegDiff;
    //checklist = (bool*)malloc(nrow*sizeof(bool));
    
    for (size_t p = 0; p < nPoints; p++) {
      if (!frontFlag[p]) continue; // skip this point if not of interest
      // loop over other points
      for (size_t q = p + 1; q < nPoints; q++) {
        // skip this point if it is dominated or of no interest
        if (!frontFlag[q]) continue;
        // otherwise compare their objective-differences
        nPosDiff = 0;
        nNegDiff = 0;
        for (size_t j = 0; j < nObjs; j++) {
          //if (obj[p * nObjs + j] < obj[q * nObjs + j]) nPosDiff = nPosDiff + 1; to filter out identical items
          if (obj[p * nObjs + j] > obj[q * nObjs + j]) nNegDiff = nNegDiff + 1;
          else nPosDiff = nPosDiff + 1;
          // speed up step
          if (nNegDiff > 0 && nPosDiff > 0) break; // incomparable
        }
        if (nNegDiff > 0 && nPosDiff > 0) continue; // incomparable
        else if (nNegDiff > 0) { // p is dominated
          frontFlag[p] = false;
          break; 
        }
        else if (nPosDiff > 0) frontFlag[q] = false; // q is dominated
        
      }
    }
    
    //free(checklist); 
}


// returns the potentially optimal hyper-rectangle
void nd_direct(bool *frontFlag, double *obj, int *level, size_t nPoints, size_t nObjs) {
    
    //bool *checklist, colDominatedFlag;
    size_t nPosDiff;
    size_t nNegDiff;
    //checklist = (bool*)malloc(nrow*sizeof(bool));
    
    for (size_t p = 0; p < nPoints; p++) {
      if (!frontFlag[p]) continue; // skip this point if not of interest
      // loop over other points
      for (size_t q = p + 1; q < nPoints; q++) {
        // skip this point if it is dominated or of no interest
        if (!frontFlag[q]) continue;
        // otherwise compare their objective-differences
        nPosDiff = 0;
        nNegDiff = 0;
		// compare their sizes
        if (level[p] > level[q]) nNegDiff = nNegDiff + 1;	
        else nPosDiff = nPosDiff + 1;		  
		// compare the values
        for (size_t j = 0; j < nObjs; j++) {
		  if (obj[p * nObjs + j] > obj[q * nObjs + j]) nNegDiff = nNegDiff + 1;
          else nPosDiff = nPosDiff + 1;

          // speed up step
          if (nNegDiff > 0 && nPosDiff > 0) break; // incomparable
        }
        if (nNegDiff > 0 && nPosDiff > 0) continue; // incomparable
        else if (nNegDiff > 0) { // p is dominated
          frontFlag[p] = false;
          break; 
        }
        else if (nPosDiff > 0) frontFlag[q] = false; // q is dominated
        
      }
    }
    
    //free(checklist); 
}


#ifdef __cplusplus
}
#endif
