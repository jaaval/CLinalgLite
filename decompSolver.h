#ifndef DECOMPSOLVER_H
#define DECOMPSOLVER_H

#include "matrixOps.h"

#ifdef __cplusplus
extern "C" {
#endif


typedef struct DecompData {
	int size;
	int cols;
	Matrix* D; //ldlt diagonal vector
	Matrix* L; //triangular
	Matrix* X; //solver solution
} DecompData;

void initDecomp(DecompData* dec, int size); //TODO
void resetDecomp(DecompData* dec);
//void deleteDecomp(DecompData* dec);

void cholDecompose(DecompData* chol, Matrix* A);
void cholSolve(DecompData* chol, Matrix* b);
void ldlDecompose(DecompData* ldlt, Matrix* A);
void ldlSolve(DecompData* ldlt, Matrix* b);
//void ldlInverse(DecompData* ldlt, Matrix* A);

#ifdef __cplusplus
}
#endif

#endif