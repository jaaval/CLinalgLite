#ifndef CHOLESKY_H
#define CHOLESKY_H

#include "matrixOps.h"

#ifdef __cplusplus
extern "C" {
#endif


typedef struct CholDecomp {
	int size;
	Matrix* D; //ldlt diagonal vector
	Matrix* L; //triangular
	Matrix* X; //solver solution
} CholDecomp;

CholDecomp* newCholDecomp(int size);
void initCholDecomp(CholDecomp* dec, int size);
void resetCholDecomp(CholDecomp* dec, int size);
void deleteCholDecomp(CholDecomp* dec);

void cholDecompose(CholDecomp* chol, Matrix* A);
void cholSolve(CholDecomp* chol, Matrix* b);
void ldlDecompose(CholDecomp* ldlt, Matrix* A);
void ldlSolve(CholDecomp* ldlt, Matrix* b);
//void ldlInverse(DecompData* ldlt, Matrix* A);

#ifdef __cplusplus
}
#endif

#endif
