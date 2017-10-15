#ifndef PARETOFILTERING_H
#define PARETOFILTERING_H

#include <stdbool.h> // to use the bool datatype, required C99
#ifndef RECORDER_H
#include "recorder.h"
#endif

#ifdef __cplusplus
extern "C" {
#endif

void paretoFiltering(struct SolutionsArchive *archive);
void paretofront(bool *frontFlag, double *obj, unsigned nrow, unsigned ncol);
// perform pareto front among selected elements of *obj as specified by frontFlag
void selectiveparetofront(bool *frontFlag, double *obj, size_t nPoints, size_t nObjs);
void nd_direct(bool *frontFlag, double *obj, int *level, size_t nPoints, size_t nObjs);

#ifdef __cplusplus
}
#endif

#endif /* PARETOFILTERING_H */
