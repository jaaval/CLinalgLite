/**
	Some basic linear algebra functionality defined
	Only the core functionality of the system implemented - basically no error checking
	i.e. if you use it wrong it will crash - you are supposed to know what you put in to the functions
*/

#ifndef __MATRIXOPS__
#define __MATRIXOPS__


#ifdef __cplusplus
extern "C" {
#endif

#define TRUE 1
#define FALSE 0

typedef struct Matrix {
	int cols;
	int rows;
	int transposed; //TODO this can be used for transpose in the future
	double* array;
} Matrix;

//constructors
//for now lets just assume memory allocation succeeds always
void createMatrix(Matrix* mat, int r, int c);
Matrix* newMatrix(int r, int c);

// "destructor" releases all the memory and invalidates the matrix pointer
void deleteMatrix(Matrix* mat);

//this can be used to resize matrices in other functions
//should not be needed if things are done correctly
void resizeMatrix(Matrix* mat, int rows, int cols);


// define some basic matrix functionality
// all matrices are in row major format

// scalar mult in place
void scalarMultInPlace(Matrix* matrix, double scalar);
// scalar mult
void scalarMult(Matrix* matA, Matrix* matB, double scalar);

// matrix multiplication
void matrixMult(Matrix* matA, Matrix* matB, Matrix* matC);

// matrix multiplication with a transpose on matB
void matrixMultTranspose(Matrix* matA, Matrix* matB, Matrix* matC);

//when it is not an array
double normXYZ(double x, double y, double z);

//cross product of two length 3 vectors
void cross(Matrix* vecA, Matrix* vecB, Matrix* vecC);

// matrix addition
void plusMatrixInPlace(Matrix* matA, Matrix* matB);

void minusMatrixInPlace(Matrix* matA, Matrix* matB);

void plusMatrix(Matrix* matA, Matrix* matB, Matrix* matC);

// scalar addition to matrix
void plusScalarInPlace(Matrix* matA, double scalar);

// set to
void setIdentity(Matrix* matrix);

void setTo(Matrix* matrix, double value);
void setZero(Matrix* matrix);

void setDiag(Matrix* matrix, double value);

void copyMat(Matrix* matA, Matrix* matB);

void transpose(Matrix* matA, Matrix* matB);


//copies a vector from index ind. if col -> copy a column else copy a row
void extractVector(Matrix* matA, Matrix* matB, int ind, int col);

//copies a block from matrix A to matriB that is of correct size
void copyBlock(Matrix* matA, Matrix* matB, int y0, int x0, int sizey, int sizex);

//copies matrix A to a block of matrix B
void setBlock(Matrix* matA, Matrix* matB, int y0, int x0, int sizey, int sizex);

//helper for error checking
int sizeEquals(Matrix* matA, Matrix* matB);

// setter and getter
double get(Matrix* mat, int i, int j);
void set(Matrix* mat, int i, int j, double val);

// extra
double det(Matrix* mat);
double diagSum(Matrix* mat);


// inversion for matrices. this is easier than writing an actual solver
// input has to be square matrix
//void stupidInv(Matrix* matA, Matrix* matB);

#ifdef __cplusplus
}
#endif

#endif
