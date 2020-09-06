/*
 * b.h
 *
 *  Created on: Aug 14, 2020
 *  	Author: Muhammad Watad && Shadi Abu Shqara
 */

#ifndef B_H_
#define B_H_

#include "a.h"

#define EPSILON 0.00001
#define IS_POSITIVE(X) ((X) > EPSILON)


typedef struct B_t {

	A *a;
	int *g;
	int n_g;

	/* constant vector which is calculated only once */
	double *f;
	double b_norm;

} B;

typedef enum { false, true } boolean;
typedef enum { moved, unmoved } status;


B* B_create(A *a, int n_g);

void B_destroy(B *b);

void B_initialize_g(B *b);

void B_calculate_f(B *b);

double B_calculate_1_norm(B *b, double *counters, double *aux_array);

void B_multiply(B *b, double *v, double *res);

void B_shifted_multiply(B *b, double *v, double *res);

boolean B_calculate_next(B *b, double *v_prev, double *v_next);

double B_find_leading_eigenpair(B *b, double *leading_v, double *counters, double *aux_array);

boolean B_divide_into_two_groups(B *b, B **b1, B **b2, double *leading_v, double *counters, double *aux_array);

void B_maximize_modularity(B *b, B *b1, B *b2, double *s, status *unmoved,
	double *aux_array, double *improve, int *indices, boolean b_was_used);

double B_calculate_modularity(B *b, double *s, double *res);

#endif /* B_H_ */
