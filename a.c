/*
 * a.c
 *
 *  Created on: Aug 14, 2020
 *  	Author: Muhammad Watad && Shadi Abu Shqara
 */



#include <stdlib.h>
#include "a.h"
#include "error.h"


A* A_create(int n, int nnz) {

	A* a = (A*)malloc(sizeof(A));
	exit_if_no_memory(a);

	a->n = n;

	a->colind = (int*)malloc(nnz * sizeof(int));
	exit_if_no_memory(a->colind);

	a->rowptr = (int*)malloc((n + 1) * sizeof(int));
	exit_if_no_memory(a->rowptr);
	(a->rowptr)[n] = nnz;

	a->first_empty = 0;
	a->nnz = nnz;

	return a;

}


void A_add_row(A* a, const int *row, int size, int i) {

	unsigned int j;
	register int col;
	register const int *ptr = row;
	register int *ptr_1;

	j = a->first_empty;
	(a->rowptr)[i] = j;
	ptr_1 = a->colind + j;

	for (col = 0; col < size; col++) {

		*ptr_1++ = *ptr;
		a->first_empty++;
		ptr++;

	}

}

void A_destroy(A *a) {

	free(a->colind);
	free(a->rowptr);
	free(a->k);
	free(a);

}

void A_multiply(A* a, const double *v, double *res, int *g, int n_g) {

	register int i;
	int g_i;
	register int *row_ptr, *i_ptr;
	int *i_ptr_last, *g_ptr_last = g + n_g, *colind = a->colind;
	register int *g_ptr = g, *ptr;

	row_ptr = (a->rowptr);
	ptr = g;

	for (i = 0; i < n_g; i++) {

		g_i = *ptr;

		i_ptr = colind + row_ptr[g_i];
		i_ptr_last = colind + row_ptr[g_i + 1];
		g_ptr = g;

		while (i_ptr < i_ptr_last && g_ptr < g_ptr_last) {

			if (*i_ptr < *g_ptr) {
				i_ptr++;
			}
			else if (*g_ptr < *i_ptr) {
				g_ptr++;
			}
			else {
				res[g_i] += v[*g_ptr];
				g_ptr++;
				i_ptr++;
			}
		}

		ptr++;
	}

}
