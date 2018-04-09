#define TESTING // undefine to not mess with arduino build
#ifdef TESTING
#include "matrixOps.h"
#include "cholesky.h"
#include "eigen.h"
#include <stdlib.h>
#include <stdio.h>

int main() {
	printf("Hello!\n");
	Matrix* A = newMatrix(3,3);
	Matrix* b = newMatrix(3,3);
	CholDecomp* ldlt = newCholDecomp(3);
	printf("Cholesky (and LDLT) tests:\n");
	setZero(A);
	printf("--------\n");
	set(A,0,0, 2);
	set(A,0,1, -1);
	set(A,0,2, 0);
	set(A,1,0, -1);
	set(A,1,1, 2);
	set(A,1,2, -1);
	set(A,2,0, 0);
	set(A,2,1, -1);
	set(A,2,2, 1);
	printf("--------\n");
	set(b,0,0, 3);
	set(b,0,1, -1);
	set(b,0,2, 1);
	set(b,1,0, -2);
	set(b,1,1, 3);
	set(b,1,2, 0);
	set(b,2,0, 1);
	set(b,2,1, -2);
	set(b,2,2, 0);
	printf("--------\n");
	resetCholDecomp(ldlt, 3);
	printf("--------\n");
	ldlDecompose(ldlt, A);
	printf("%f  %f  %f\n%f  %f  %f\n%f  %f  %f\n",
		get(ldlt->L, 0,0),get(ldlt->L, 0,1),get(ldlt->L, 0,2),
		get(ldlt->L, 1,0),get(ldlt->L, 1,1),get(ldlt->L, 1,2),
		get(ldlt->L, 2,0),get(ldlt->L, 2,1),get(ldlt->L, 2,2));
	printf("--------\n");
	printf("%f  %f  %f\n",
		ldlt->D->array[0],ldlt->D->array[1],ldlt->D->array[2]);
	printf("--------\n");
	ldlSolve(ldlt, b);
	printf("%f  %f  %f\n%f  %f  %f\n%f  %f  %f\n",
		get(ldlt->X, 0,0),get(ldlt->X, 0,1),get(ldlt->X, 0,2),
		get(ldlt->X, 1,0),get(ldlt->X, 1,1),get(ldlt->X, 1,2),
		get(ldlt->X, 2,0),get(ldlt->X, 2,1),get(ldlt->X, 2,2));
	printf("--------\n");
	resetCholDecomp(ldlt, 3);
	printf("--------\n");
	cholDecompose(ldlt, A);
	printf("%f  %f  %f\n%f  %f  %f\n%f  %f  %f\n",
		get(ldlt->L, 0,0),get(ldlt->L, 0,1),get(ldlt->L, 0,2),
		get(ldlt->L, 1,0),get(ldlt->L, 1,1),get(ldlt->L, 1,2),
		get(ldlt->L, 2,0),get(ldlt->L, 2,1),get(ldlt->L, 2,2));
	printf("--------\n");
	printf("%f  %f  %f\n",
		ldlt->D->array[0],ldlt->D->array[1],ldlt->D->array[2]);
	printf("--------\n");
	cholSolve(ldlt, b);
	printf("%f  %f  %f\n%f  %f  %f\n%f  %f  %f\n",
		get(ldlt->X, 0,0),get(ldlt->X, 0,1),get(ldlt->X, 0,2),
		get(ldlt->X, 1,0),get(ldlt->X, 1,1),get(ldlt->X, 1,2),
		get(ldlt->X, 2,0),get(ldlt->X, 2,1),get(ldlt->X, 2,2));
	printf("--------\n");

	deleteMatrix(A);
	deleteMatrix(b);
	deleteCholDecomp(ldlt);

	printf("Eigen tests:\n");
	printf("2x2:\n");
	A = newMatrix(2,2);
	EigenDecomp* eig = newEigenDecomp(2);

	set(A,0,0, 2);
	set(A,0,1, 1);
	set(A,1,0, 3);
	set(A,1,1, 6);
	eigenDecompose(eig, A);
	printf("--------\n");
	printf("%f  %f\n",
		eig->lambda->array[0],eig->lambda->array[1]);
	printf("--------\n");
	printf("%f  %f\n%f  %f\n",
		get(eig->vectors, 0,0),get(eig->vectors, 0,1),
		get(eig->vectors, 1,0),get(eig->vectors, 1,1));
	printf("--------\n");


	deleteMatrix(A);
	deleteEigen(eig);
	return 0;
}

#endif
