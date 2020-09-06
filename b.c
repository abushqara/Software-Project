/*
 * b.c
 *
 *  Created on: Aug 14, 2020
 *  	Author: Muhammad Watad && Shadi Abu Shqara
 */

#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "error.h"
#include "division.h"



int num_of_shared_values(int *arr1, int *arr2, int m, int n) {

	int cnt = 0;

	register int *last_arr1 = arr1 + m;
	register int *last_arr2 = arr2 + n;

	while (arr1 < last_arr1 && arr2 < last_arr2) {

		if (*arr1 < *arr2)
			arr1++;

		else if (*arr2 < *arr1)
			arr2++;

		else {
			arr1++;
			arr2++;
			cnt++;
		}
	}

	return cnt;

}


B* B_create(A* a, int n_g) {

	B* b = (B*)malloc(sizeof(B));
	exit_if_no_memory(b);

	if (n_g == 0) {
		b->g = NULL;
		b->f = NULL;
	}
	else {
		b->g = (int*)malloc(sizeof(int) * n_g);
		exit_if_no_memory(b->g);

		b->f = (double*)calloc(n_g, sizeof(double));
		exit_if_no_memory(b->f);
	}

	b->a = a;
	b->n_g = n_g;
	b->b_norm = 0;

	return b;

}

void B_destroy(B* b) {

	free(b->g);
	free(b->f);
	free(b);

}

void B_initialize_g(B* b) {

	int n_g = b->n_g, *g = b->g;
	register int i;

	for (i = 0; i < n_g; i++) {
		*g = i;
		g++;
	}

	B_calculate_f(b);

}

void B_calculate_f(B *b) {

	int n_g, *k_array, *g, *colind, *rowptr, g_i, k_i;
	double m;
	register int *g_ptr_i, *g_ptr_j;
	register int i, j;
	register double *f_ptr;

	k_array = b->a->k;
	f_ptr = b->f;
	m = b->a->m;
	n_g = b->n_g;
	g = b->g;
	g_ptr_i = g;
	g_ptr_j = g;

	colind = b->a->colind;
	rowptr = b->a->rowptr;

	exit_if_dividing_by_zero(m);

	for (i = 0; i < n_g; i++) {

		g_ptr_j = g;
		g_i = *g_ptr_i;
		k_i = k_array[g_i];

		for (j = 0; j < n_g; j++) {

			*f_ptr -= (k_i * k_array[*g_ptr_j]) / m;
			g_ptr_j++;

		}

		*f_ptr += num_of_shared_values(g, colind + rowptr[g_i], n_g, k_i);
		g_ptr_i++;
		f_ptr++;

	}

}


/*
 * Assumes size > 0 and element are non negative
 */
double Max(double *array, int size) {

	register int i;
	register double *ptr = array;
	double max = 0;
	double val;

	for (i = 0; i < size; i++) {
		val = *ptr;
		max = val > max ? val : max;
		ptr++;
	}

	return max;

}

/* counters is a pre-allocated array of size b->a->n of zeros */
/* aux_array is a pre-allocated array of size b->a->n of zeros */
double B_calculate_1_norm(B *b, double *counters, double *aux_array) {

	int n_g, *g, *colind, *rowptr, *k_array, len, g_i, *g_ptr_last, *i_ptr_last;
	register int *i_ptr, *g_ptr1, *g_ptr2;
	register int i;
	double m;
	register double *cnt_ptr, *aux_ptr, *f_ptr;
	register double val = 0;

	n_g = b->n_g;
	g_ptr1 = b->g;
	g_ptr2 = b->g;
	g = b->g;
	colind = b->a->colind;
	rowptr = b->a->rowptr;
	k_array = b->a->k;
	m = b->a->m;
	cnt_ptr = counters;
	aux_ptr = aux_array;

	g_ptr_last = g + n_g;

	for (i = 0; i < n_g; i++) {

		g_i = *g_ptr2;

		len = rowptr[g_i + 1] - rowptr[g_i];
		i_ptr = colind + rowptr[g_i];
		i_ptr_last = i_ptr + len;
		g_ptr1 = b->g;

		while (i_ptr < i_ptr_last && g_ptr1 < g_ptr_last) {

			val = 0;

			if (*i_ptr < *g_ptr1)
				i_ptr++;

			else if (*g_ptr1 < *i_ptr) {
				val = -(k_array[g_i] * k_array[*g_ptr1]) / m;

				if (g_i == *g_ptr1)
					*aux_ptr = val;
				else
					*cnt_ptr += fabs(val);

				g_ptr1++;
			}

			else {
				val = -(k_array[g_i] * k_array[*g_ptr1]) / m;

				if (g_i == *g_ptr1)
					*aux_ptr = val;
				else
					*cnt_ptr += fabs(1 + val);

				g_ptr1++;
				i_ptr++;
			}

		}

		while (g_ptr1 < g_ptr_last) {
			val = -(k_array[g_i] * k_array[*g_ptr1]) / m;

			if (g_i == *g_ptr1)
				*aux_ptr = val;
			else
				*cnt_ptr += fabs(val);

			g_ptr1++;
		}


		g_ptr2++;
		cnt_ptr++;
		aux_ptr++;

	}

	cnt_ptr = counters;
	aux_ptr = aux_array;
	f_ptr = b->f;

	for (i = 0; i < n_g; i++) {
		*cnt_ptr += fabs(*aux_ptr - *f_ptr);
		cnt_ptr++;
		aux_ptr++;
		f_ptr++;
	}

	return Max(counters, n_g);

}

void B_multiply(B *b, double *v, double *res) {

	int n_g = b->n_g, *k_array = b->a->k, *g = b->g, m = b->a->m, g_i;
	register double k_v = 0;
	register int i;
	register int *g_ptr;
	register double *f_ptr;

	g_ptr = b->g;

	exit_if_dividing_by_zero(m);

	memset(res, 0, b->a->n * sizeof(double));
	A_multiply(b->a, v, res, g, n_g);

	for (i = 0; i < n_g; i++) {
		g_i = *g_ptr;
		k_v += k_array[g_i] * v[g_i];
		g_ptr++;
	}

	g_ptr = b->g;
	for (i = 0; i < n_g; i++) {
		g_i = *g_ptr;
		res[g_i] -= (k_array[g_i] * k_v) / m;
		g_ptr++;
	}

	g_ptr = b->g;
	f_ptr = b->f;

	for (i = 0; i < n_g; i++) {
		g_i = *g_ptr;
		res[g_i] -= (*f_ptr) * v[g_i];
		g_ptr++;
		f_ptr++;
	}

}

void B_shifted_multiply(B *b, double *v, double *res) {

	int n_g = b->n_g, *g_ptr = b->g, g_i;
	double b_norm = b->b_norm;
	register int i;

	B_multiply(b, v, res);

	for (i = 0; i < n_g; i++) {
		g_i = *g_ptr;
		res[g_i] += v[g_i] * b_norm;
		g_ptr++;
	}

}


boolean B_calculate_next(B *b, double *v_prev, double *v_next) {

	int n_g = b->n_g, g_i;
	register double norm = 0, val = 0, norm_inverse;
	register int i, *g_ptr = b->g;
	boolean res = true;

	B_shifted_multiply(b, v_prev, v_next);

	for (i = 0; i < n_g; i++) {
		val = v_next[*g_ptr];
		norm += val * val;
		g_ptr++;
	}

	norm = sqrt(norm);
	norm_inverse = 1 / norm;

	g_ptr = b->g;

	for (i = 0; i < n_g; i++) {

		g_i = *g_ptr;
		v_next[g_i] *= norm_inverse;

		if (res && fabs(v_next[g_i] - v_prev[g_i]) >= EPSILON)
			res = false;

		g_ptr++;

	}

	return res;

}

void create_random_vector(double *init_vector, int n_g, int *g) {

	register int i;
	register int *g_ptr = g;

	for (i = 0; i < n_g; i++) {
		init_vector[*g_ptr] = rand();
		g_ptr++;
	}

}

double vector_dot_product(double *v_1, double *v_2, int *g, int n_g) {

	int g_i;
	register int i, *g_ptr = g;
	register double dot_product = 0;

	for (i = 0; i < n_g; i++) {
		g_i = *g_ptr;
		dot_product += v_1[g_i] * v_2[g_i];
		g_ptr++;
	}

	return dot_product;
}

/*
 *
 * Assumes init_vector is a ZERO vector
 *
 */
double B_find_leading_eigenpair(B *b, double *leading_v, double *counters, double *aux_array) {

	double *v_prev, *v_next, *temp;
	double numerator, denominator;

	int *g = b->g;
	int n_g = b->n_g;

	b->b_norm = B_calculate_1_norm(b, counters, aux_array);

	create_random_vector(aux_array, n_g, g);

	v_prev = aux_array;
	v_next = leading_v;

	while (1) {

		if (!B_calculate_next(b, v_prev, v_next)) {
			temp = v_prev;
			v_prev = v_next;
			v_next = temp;
		}

		else
			break;
	}

	B_shifted_multiply(b, leading_v, aux_array);

	numerator = vector_dot_product(leading_v, aux_array, g, n_g);
	denominator = vector_dot_product(leading_v, leading_v, g, n_g);
	exit_if_dividing_by_zero(denominator);

	return numerator / denominator - b->b_norm;
}

int num_of_positive_entries(double *leading_v, int *g, int n_g) {

	register int i, cnt = 0;
	register int *g_ptr = g;

	for (i = 0; i < n_g; i++) {
		if (IS_POSITIVE(leading_v[*g_ptr]))
			cnt++;
		g_ptr++;
	}

	return cnt;

}

boolean B_divide_into_two_groups(B *b, B **b1, B **b2, double *leading_v, double *counters, double *aux_array) {

	int n_g = b->n_g, *g, g_i;
	register int i;
	register int *g_ptr;
	int n_g1, *ptr_1, *ptr_2;
	double modularity;
	double eigen_val = B_find_leading_eigenpair(b, leading_v, counters, aux_array);

	g = b->g;
	g_ptr = g;

	for (i = 0; i < n_g; i++) {

		g_i = *g_ptr;

		if (IS_POSITIVE(leading_v[g_i]))
			leading_v[g_i] = 1;
		else
			leading_v[g_i] = -1;

		g_ptr++;

	}

	modularity = B_calculate_modularity(b, leading_v, aux_array);

	if (eigen_val <= 0 || !IS_POSITIVE(modularity)) {

		g_ptr = g;
		for (i = 0; i < n_g; i++) {
			leading_v[*g_ptr] = 1;
			g_ptr++;
		}

		*b1 = b;
		*b2 = B_create(b->a, 0);

		return true;

	}

	else {

		n_g1 = num_of_positive_entries(leading_v, g, n_g);
		*b1 = B_create(b->a, n_g1);
		*b2 = B_create(b->a, n_g - n_g1);

		ptr_1 = (*b1)->g;
		ptr_2 = (*b2)->g;
		g_ptr = g;

		for (i = 0; i < n_g; i++) {

			g_i = *g_ptr;

			if (leading_v[g_i] == 1) {
				*ptr_1 = g_i;
				ptr_1++;
			}
			else {
				*ptr_2 = g_i;
				ptr_2++;
			}

			g_ptr++;
		}

		return false;
	}

}

double B_calculate_modularity(B *b, double *s, double *res) {

	int n_g, g_i;
	register int i, *g_ptr;
	register double sum = 0;

	B_multiply(b, s, res);
	n_g = b->n_g;
	g_ptr = b->g;

	for (i = 0; i < n_g; i++) {
		g_i = *g_ptr;
		sum += s[g_i] * res[g_i];
		g_ptr++;
	}

	return sum;

}


void B_update_groups(B *b, B *b1, B *b2, double *s, boolean b_eq_b1) {

	int n_g1, n_g, *ptr_1, *ptr_2, g_i, *g = b->g;
	register int i, *g_ptr = b->g;

	n_g = b->n_g;
	n_g1 = num_of_positive_entries(s, b->g, n_g);

	if (!b_eq_b1) {
		free(b1->g);
	}

	free(b1->f);
	free(b2->g);
	free(b2->f);

	b1->g = (int*)malloc(n_g1 * sizeof(int));
	exit_if_no_memory(b1->g);
	b1->f = (double*)calloc(n_g1, sizeof(double));
	exit_if_no_memory(b1->f);
	b2->g = (int*)malloc((n_g - n_g1) * sizeof(int));
	exit_if_no_memory(b2->g);
	b2->f = (double*)calloc((n_g - n_g1), sizeof(double));
	exit_if_no_memory(b2->f);

	b1->n_g = n_g1;
	b2->n_g = n_g - n_g1;

	ptr_1 = b1->g;
	ptr_2 = b2->g;

	for (i = 0; i < n_g; i++) {

		g_i = *g_ptr;

		if (s[g_i] == 1) {
			*ptr_1 = g_i;
			ptr_1++;
		}
		else {
			*ptr_2 = g_i;
			ptr_2++;
		}

		g_ptr++;
	}

	if (b_eq_b1) {
		free(g);
	}
	B_calculate_f(b1);
	B_calculate_f(b2);

}

double delta_modularity(B *b, double *d, int i) {

	int *rowptr, *i_last, *g_last, *k_array, *g, n_g, g_j, g_i, k_g_i;
	double d_g_i;
	register int *i_ptr, *g_ptr;
	register int j;
	double m, res;

	rowptr = b->a->rowptr;
	k_array = b->a->k;
	m = b->a->m;
	n_g = b->n_g;
	g = b->g;

	g_ptr = g;
	g_i = g[i];
	d_g_i = d[g_i];
	k_g_i = k_array[g_i];
	i_ptr = b->a->colind + rowptr[g_i];
	i_last = i_ptr + k_g_i;
	g_last = g + b->n_g;

	res = (4 * k_g_i * k_g_i) / m;

	for (j = 0; j < n_g; j++) {
		g_j = *g_ptr;
		res -= 4 * d_g_i * d[g_j] * ((k_g_i * k_array[g_j]) / m);
		g_ptr++;
	}

	g_ptr = g;

	while (g_ptr < g_last && i_ptr < i_last) {

		if (*i_ptr < *g_ptr) {
			i_ptr++;
		}
		else if (*g_ptr < *i_ptr) {
			g_ptr++;
		}
		else {
			res += 4 * d_g_i * d[*g_ptr];
			g_ptr++;
			i_ptr++;
		}
	}

	return res;

}

void B_maximize_modularity(B *b, B *b1, B *b2, double *s, status *unmoved_arr,
	double *aux_array, double *improve, int *indices, boolean b_was_used) {

	int n_g, j_tag = -1, i_tag, *g, g_k;
	register int i, j, k;
	register status* um_ptr;
	register int *g_ptr, *ind_ptr, iterations = 0;
	register double *scr_ptr, *imp_ptr;
	double *score, max, delta_Q, val;

	delta_Q = 0;
	score = aux_array;

	n_g = b->n_g;
	g = b->g;

	do {

		um_ptr = unmoved_arr;
		for (i = 0; i < n_g; i++) {
			*um_ptr = unmoved;
			um_ptr++;
		}

		ind_ptr = indices;
		imp_ptr = improve;

		for (i = 0; i < n_g; i++) {

			g_ptr = g;
			um_ptr = unmoved_arr;
			scr_ptr = score;

			for (k = 0; k < n_g; k++) {
				if (*um_ptr == unmoved) {
					g_k = *g_ptr;
					s[g_k] *= -1;
					*scr_ptr = delta_modularity(b, s, k);
					max = *scr_ptr;
					j_tag = k;
					s[g_k] *= -1;
				}
				g_ptr++;
				um_ptr++;
				scr_ptr++;
			}

			scr_ptr = score;
			um_ptr = unmoved_arr;
			for (j = 0; j < n_g; j++) {

				if (*scr_ptr > max && *um_ptr == unmoved) {
					max = *scr_ptr;
					j_tag = j;
				}
				scr_ptr++;
				um_ptr++;
			}

			s[g[j_tag]] *= -1;
			*ind_ptr = j_tag;

			if (i == 0)
				*imp_ptr = score[j_tag];
			else
				*imp_ptr = *(imp_ptr - 1) + score[j_tag];

			unmoved_arr[j_tag] = moved;

			ind_ptr++;
			imp_ptr++;

		}

		max = improve[0];
		i_tag = 0;

		imp_ptr = improve;

		for (i = 0; i < n_g; i++) {

			val = *imp_ptr;

			if (val > max) {
				max = val;
				i_tag = i;
			}

			imp_ptr++;
		}

		ind_ptr = indices + n_g - 1;

		for (i = n_g - 1; i > i_tag; i--) {
			j = *ind_ptr;
			s[g[j]] *= -1;
			ind_ptr--;
		}

		if (i_tag == n_g - 1)
			delta_Q = 0;
		else
			delta_Q = improve[i_tag];

		iterations++;

	} while (IS_POSITIVE(delta_Q));

	if (iterations == 1 && !IS_POSITIVE(delta_Q))
		return;

	B_update_groups(b, b1, b2, s, b_was_used);

}

