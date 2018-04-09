
#include "eigen.h"
#include <stdlib.h>
#include <math.h> //sqrt

void initEigen(EigenDecomp* dec, int size) {
	dec->size = size;
	dec->vectors = (Matrix*)malloc(sizeof(Matrix));
	dec->lambda = (Matrix*)malloc(sizeof(Matrix));
	createMatrix(dec->vectors,size,size);
	createMatrix(dec->lambda,size,1);
}

void resetEigen(EigenDecomp* dec, int size) {
  if (size != dec->size) {
    dec->size = size;
    resizeMatrix(dec->vectors, size, size);
    resizeMatrix(dec->lambda, size, 1);
  }
  setZero(dec->vectors);
  setZero(dec->lambda);
}

void deleteEigen(EigenDecomp* dec) {
	deleteMatrix(dec->vectors);
	deleteMatrix(dec->lambda);
	free(dec);
}


void eigenDecompose(EigenDecomp* dec, Matrix* A){
  // special cases for size 2 and 3 decompositions
  if (dec->size == 2) {
    __eigenDecomp2(dec, A);
  }
  else if (dec->size == 3) {
    __eigenDecomp3(dec, A);
  }
  else {

  }
}

void __eigenDecomp2(EigenDecomp* dec, Matrix* A) {
  double temp;
  double trA = A->array[0] + A->array[3];
  double detA = A->array[0] * A->array[3] - A->array[1] * A->array[2];
  double sqdet = sqrt(trA*trA - 4*detA);
  dec->lambda->array[0] = (trA + sqdet) / 2.0;
  dec->lambda->array[1] = (trA - sqdet) / 2.0;
  dec->vectors->array[0] = A->array[1];
  dec->vectors->array[1] = A->array[0] - dec->lambda->array[0];
  dec->vectors->array[2] = A->array[3] - dec->lambda->array[1];
  dec->vectors->array[3] = A->array[2];
  // normalization to unity length
  temp = sqrt(pow(dec->vectors->array[0],2) + pow(dec->vectors->array[2], 2));
  dec->vectors->array[0] /= temp;
  dec->vectors->array[2] /= temp;
  temp = sqrt(pow(dec->vectors->array[1],2) + pow(dec->vectors->array[3], 2));
  dec->vectors->array[1] /= temp;
  dec->vectors->array[3] /= temp;
}

void __eigenDecomp3(EigenDecomp* dec, Matrix* A) {
  return;
}
