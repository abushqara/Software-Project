/*
 * a.h
 *
 *  Created on: Aug 14, 2020
 *  	Author: Muhammad Watad && Shadi Abu Shqara
 */

#ifndef A_H_
#define A_H_

typedef struct A_t {

	/* Matrix size (n*n) */
	int		n;
	int *colind;
	int *rowptr;
	unsigned int first_empty;
	unsigned int nnz;
	int m;
	int *k;

} A;



/* Allocates a new arrays sparse matrix of size n with nnz non-zero elements */
A* A_create(int n, int nnz);

void A_destroy(A *a);

void A_add_row(A *a, const int *row, int size, int i);

void A_multiply(A *a, const double *v, double *result, int *g, int n_g);


#endif /* A_H_ */
#pragma once
