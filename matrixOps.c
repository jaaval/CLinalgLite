/**
	Some basic linear algebra functionality defined
	Only the core functionality of the system implemented - basically no error checking
	i.e. if you use it wrong it will crash - you are supposed to know what you put in to the functions
*/

#include "matrixOps.h"
#include <math.h> //sqrt
#include <stdlib.h>

//for now lets just assume memory allocation succeeds always
void createMatrix(Matrix* mat, int r, int c) {
	mat->cols = c;
	mat->rows = r;
	mat->transposed = FALSE;
	mat->array = (double *)malloc(c*r*sizeof(double));
}

Matrix* newMatrix(int r, int c) {
	Matrix* mat = (Matrix *)malloc(sizeof(Matrix));
	createMatrix(mat, r, c);
	return mat;
}

// "destructor" releases all the memory and invalidates the matrix pointer
void deleteMatrix(Matrix* mat) {
	free(mat->array);
	free(mat);
}

void resizeMatrix(Matrix* mat, int r, int c) {
	if (mat->cols != c || mat->rows != r) {
		free(mat->array);
		createMatrix(mat, r, c);
	}
}


// define some basic matrix functionality
// all matrices are in row major format

// scalar mult in place
void scalarMultInPlace(Matrix* matrix, double scalar) {
	for (int i = 0; i < matrix->rows * matrix->cols; i++) {
		matrix->array[i] *= scalar;
	}
}
// scalar mult
void scalarMult(Matrix* matA, Matrix* matB, double scalar) {
	for (int i = 0; i < matA->rows * matA->cols; i++) {
			matB->array[i] = matA->array[i]*scalar;
	}
}

// matrix multiplication
void matrixMult(Matrix* matA, Matrix* matB, Matrix* matC) {
	double sum;
	for (int i = 0; i < matA->rows; i++) {
		for (int j = 0; j < matB->cols; j++) {
			sum = 0;
			for (int k = 0; k < matA->cols; k++) {
				sum += get(matA, i, k) * get(matB, k, j);
			}
			set(matC, i, j, sum);
		}
	}
}

// matrix multiplication with a transpose on matB
void matrixMultTranspose(Matrix* matA, Matrix* matB, Matrix* matC) {
	double sum;
	for (int i = 0; i < matA->rows; i++) {
		for (int j = 0; j < matB->rows; j++) {
			sum = 0;
			for (int k = 0; k < matA->cols; k++) {
				sum += get(matA, i, k) * get(matB, j, k);
			}
			set(matC, i, j, sum);
		}
	}
}

//when it is not an array
double normXYZ(double x, double y, double z) {
	return sqrt(x*x+y*y+z*z);
}

void cross(Matrix* vecA, Matrix* vecB, Matrix* vecC) {
	vecC->array[0] = vecA->array[2]*vecB->array[3] - vecA->array[3]*vecB->array[2];
	vecC->array[1] = vecA->array[3]*vecB->array[1] - vecA->array[1]*vecB->array[3];
	vecC->array[2] = vecA->array[1]*vecB->array[2] - vecA->array[2]*vecB->array[1];
}

// matrix addition
void plusMatrixInPlace(Matrix* matA, Matrix* matB) {
	for (int i = 0; i < matA->rows*matA->cols; i++) {
		matA->array[i] += matB->array[i];
	}
}

void minusMatrixInPlace(Matrix* matA, Matrix* matB) {
	for (int i = 0; i < matA->rows*matA->cols; i++) {
		matA->array[i] -= matB->array[i];
	}
}

void plusMatrix(Matrix* matA, Matrix* matB, Matrix* matC) {
	for (int i = 0; i < matA->rows*matA->cols; i++) {
		matC->array[i] = matA->array[i] + matB->array[i];
	}
}

// scalar addition to matrix
void plusScalarInPlace(Matrix* matA, double scalar) {
	for (int i = 0; i < matA->rows*matA->cols; i++) {
		matA->array[i] += scalar;
	}
}

// set to
void setZero(Matrix* matrix) {
	setTo(matrix, 0);
}
void setTo(Matrix* matrix, double value) {
	for (int i = 0; i < matrix->rows*matrix->cols; i++) {
		matrix->array[i] = value;
	}
}

void setIdentity(Matrix* matrix) {
	setDiag(matrix, 1.0);
}

void setDiag(Matrix* matrix, double value) {
	for (int i = 0; i < matrix->rows; i++) {
		for (int j = 0; j < matrix->cols; j++) {
			if (i == j) {
				set(matrix, i, j, value);
			}
			else {
				set(matrix, i, j, 0);
			}
		}
	}
}

void copyMat(Matrix* matA, Matrix* matB) {
	for (int i = 0; i < matA->rows*matA->cols; i++) {
		matB->array[i] = matA->array[i];
	}
}

void transpose(Matrix* matA, Matrix* matB) {
	int indexA, indexB;
	for (int i = 0; i < matA->rows; i++) {
		for (int j = 0; j < matA->cols; j++) {
			indexA = i*matA->cols+j;
			indexB = j*matB->cols+i;
			matB->array[indexB] = matA->array[indexA];
		}
	}
}


void extractVector(Matrix* matA, Matrix* matB, int ind, int col) {
	if (col) {

	} else {

	}
}

void copyBlock(Matrix* matA, Matrix* matB, int y0, int x0, int sizey, int sizex) {
	int indexA;
	int indexB = 0;
	for (int i = y0; i < sizey+y0; i++) {
		for (int j = x0; j < sizex+x0; j++, indexB++) {
			indexA = i * matA->cols + j;
			matB->array[indexB] = matA->array[indexA];
		}
	}
}

void setBlock(Matrix* matA, Matrix* matB, int y0, int x0, int sizey, int sizex) {
	int indexB;
	int indexA = 0;
	for (int i = y0; i < sizey+y0; i++) {
		for (int j = x0; j < sizex+x0; j++, indexA++) {
			indexB = i * matA->cols + j;
			matB->array[indexB] = matA->array[indexA];
		}
	}
}

int sizeEquals(Matrix* matA, Matrix* matB) {
	if (matA->cols == matB->cols && matA->rows == matB->rows) {
		return 1;
	}
	return 0;
}

// setter and getter
double get(Matrix* mat, int i, int j) {
	int index;
	if (FALSE) { //if transposed. not implemented. do not use
		index = j * mat->rows + i;
	} else {
		index = i * mat->cols + j;
	}
	return mat->array[index];
}
void set(Matrix* mat, int i, int j, double val) {
	int index;
	if (FALSE) { //if transposed. not implemented. do not use
		index = j * mat->rows + i;
	} else {
		index = i * mat->cols + j;
	}
	mat->array[index] = val;
}
