#define TESTING // undefine to not mess with arduino build
#ifdef TESTING
#include "matrixOps.h"
#include "decompSolver.h"
#include <stdlib.h>
#include <stdio.h>

int main() {
	printf("Hello!\n");
	Matrix* A = (Matrix *)malloc(sizeof(Matrix));
	Matrix* b = (Matrix *)malloc(sizeof(Matrix));;
	DecompData* ldlt = (DecompData *)malloc(sizeof(DecompData));;
	printf("--------\n");
	createMatrix(A, 3,3);
	printf("--------\n");
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
	createMatrix(b, 3,3);
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
	initDecomp(ldlt,3);
	resetDecomp(ldlt);
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
	resetDecomp(ldlt);
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
	return 0;
}

#endif
