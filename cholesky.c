

#include "cholesky.h"
#include <stdlib.h>
#include <math.h>


void initDecomp(CholDecomp* dec, int size) {
	dec->size = size;
	dec->L = (Matrix*)malloc(sizeof(Matrix));
	dec->D = (Matrix*)malloc(sizeof(Matrix));
	dec->X = (Matrix*)malloc(sizeof(Matrix));
	createMatrix(dec->L,size,size);
	createMatrix(dec->D,size,1);
	createMatrix(dec->X,size,size);
}

void resetDecomp(CholDecomp* dec, int size) {
	if (size != dec->size) {
		dec->size = size;
		resizeMatrix(dec->L, size, size);
		resizeMatrix(dec->D, size, 1);
		resizeMatrix(dec->X, size, size);
	}
	setZero(dec->L);
	setZero(dec->D);
	setZero(dec->X);
}

void deleteDecomp(CholDecomp* dec) {
	deleteMatrix(dec->L);
	deleteMatrix(dec->D);
	deleteMatrix(dec->X);
	free(dec);
}

void cholDecompose(CholDecomp* choldec, Matrix* A) {
	int n = choldec->size;
	double sum;
	for (int i = 0; i < n; i++) {
		for (int j = i; j < n; j++) {
			sum = get(A, i, j);
			for (int k = 0; k < i; k++) {
				sum -= (get(choldec->L, k, i)) * (get(choldec->L, k, j));
			}
			if (i == j) {
				set(choldec->L, i, j, sqrt(sum));
			} else {
				set(choldec->L, i, j, sum / get(choldec->L, i,i));
				set(choldec->L, j, i, 0);
			}
		}
	}
}

void cholSolve(CholDecomp* choldec, Matrix* b) {
	int temp_index;
	int n = choldec->size;
	//forward solution
	for (int c = 0; c < b->cols; c++) {//a bit of a hack for now to solve for matrices
		for (int i = 0; i < n; i++) {
			set(choldec->X, i, c, get(b, i, c));
			temp_index = i * choldec->X->cols + c;
			for (int k = 0; k < i; k++) {
				choldec->X->array[temp_index] -= get(choldec->X, k,c) * get(choldec->L, k,i); // check the indeX order k/i
			}
			choldec->X->array[temp_index] /= get(choldec->L, i,i);
		}
		//backward solution
		for (int i = n-1; i >= 0; i--) {
			temp_index = i * choldec->X->cols + c;
			for (int k = i+1; k < n; k++) {
				choldec->X->array[temp_index] -= get(choldec->X, k,c) * get(choldec->L, i,k); // check the indeX order k/i
			}
			choldec->X->array[temp_index] /= get(choldec->L, i,i);
		}
	}
}

void ldlDecompose(CholDecomp* ldlt, Matrix* A) {
	int n = ldlt->size;
	int ij, ji;
	double sum;
	for (int i = 0; i < n; i++) {
		for (int j = 0; j <= i; j++) {
			ij = i * ldlt->L->cols + j;
			ji = j * ldlt->L->cols + i;
			sum = A->array[ij];
			for (int k = 0; k < j; k++) {
				sum -= (get(ldlt->L, i,k)) * (get(ldlt->L, j,k)) * (ldlt->D->array[k]);
			}
			if (i == j) {
				ldlt->D->array[i] = sum;
				ldlt->L->array[ij] = 1;
			} else {
				ldlt->L->array[ij] = sum / ldlt->D->array[j];
				ldlt->L->array[ji] = 0;
			}
		}
	}
}

void ldlSolve(CholDecomp* ldlt, Matrix* b) {
	int temp_index;
	int n = ldlt->size;
	//forward solution

	for (int c = 0; c < b->cols; c++) {//a bit of a hack for now to solve for matrices
		for (int i = 0; i < n; i++) {
			temp_index = i * ldlt->X->cols + c;
			ldlt->X->array[temp_index] = b->array[temp_index];
			for (int k = 0; k < i; k++) {
				ldlt->X->array[temp_index] -= get(ldlt->X, k,c) * get(ldlt->L, i,k);
			}
		}
		//diagonal stuff
		for (int i = 0; i < n; i++) {
			temp_index = i * ldlt->X->cols + c;
			ldlt->X->array[temp_index] /= ldlt->D->array[i];
		}
		//backward solution
		for (int i = n-1; i >= 0; i--) {
			temp_index = i * ldlt->X->cols + c;
			for (int k = i+1; k < n; k++) {
				ldlt->X->array[temp_index] -= get(ldlt->X, k,c) * get(ldlt->L, k,i);
			}
		}
	}
}


//void ldlInverse(CholDecomp* ldlt, Matrix* A) {}
