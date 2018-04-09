#ifndef EIGEN_H
#define EIGEN_H

#include "matrixOps.h"

#ifdef __cplusplus
extern "C" {
#endif

//TODO
typedef struct EigenDecomp {
	int size;
  Matrix* vectors;
  Matrix* lambda;
} EigenDecomp;


EigenDecomp* newEigenDecomp(int size);
void initEigen(EigenDecomp* dec, int size);
void resetEigen(EigenDecomp* dec, int size);
void deleteEigen(EigenDecomp* dec);

void eigenDecompose(EigenDecomp* dec, Matrix* A);

void __eigenDecomp2(EigenDecomp* dec, Matrix* A);
void __eigenDecomp3(EigenDecomp* dec, Matrix* A);


#endif
