/*
 * main.c
 *
 *  Created on: Aug 14, 2020
 *  	Author: Muhammad Watad && Shadi Abu Shqara
 */


#include <stdlib.h>
#include <string.h>
#include <time.h>
#include "error.h"
#include "parser.h"
#include "SPBufferset.h"



int main(int argc, char **argv) {

	B *b, *b1, *b2;
	A *a;
	int n, *indices;
	double *leading_v, *counters, *aux_array, *improve;
	status *unmoved_array;
	Division *P, *O;
	clock_t start, end;
	boolean b_was_used;

	start = clock();

	if (argc != 3)
		exit_insufficient_arguments();

	SP_BUFF_SET();
	srand(time(NULL));

	b = parser_read_input(argv[1]);
	a = b->a;
	n = b->a->n;

	P = division_create();
	O = division_create();

	leading_v = (double*)malloc(sizeof(double) * b->a->n);
	exit_if_no_memory(leading_v);

	counters = (double*)calloc(n, sizeof(double));
	exit_if_no_memory(counters);

	aux_array = (double*)calloc(n, sizeof(double));
	exit_if_no_memory(aux_array);

	unmoved_array = (status*)calloc(n, sizeof(status));
	exit_if_no_memory(unmoved_array);

	improve = (double*)calloc(n, sizeof(double));
	exit_if_no_memory(improve);

	indices = (int*)calloc(n, sizeof(int));
	exit_if_no_memory(indices);

	division_add_group(P, b);

	while (P->size != 0) {

		b = division_extract_group(P);
		b_was_used = B_divide_into_two_groups(b, &b1, &b2, leading_v, counters, aux_array);
		B_maximize_modularity(b, b1, b2, leading_v, unmoved_array, aux_array, improve, indices, b_was_used);

		if (!b_was_used)
			B_destroy(b);

		if (b1->n_g == 0 || b2->n_g == 0) {

			if (b1->n_g == 0) {
				B_destroy(b1);
				division_add_group(O, b2);
			}
			else {
				B_destroy(b2);
				division_add_group(O, b1);
			}
		}

		else {

			if (b1->n_g == 1)
				division_add_group(O, b1);
			else
				division_add_group(P, b1);

			if (b2->n_g == 1)
				division_add_group(O, b2);
			else
				division_add_group(P, b2);

		}

		memset(leading_v, 0, n * sizeof(double));
		memset(counters, 0, n * sizeof(double));
		memset(aux_array, 0, n * sizeof(double));
		memset(unmoved_array, 0, n * sizeof(status));
		memset(improve, 0, n * sizeof(double));
		memset(indices, 0, n * sizeof(int));


	}

	parser_write_output(argv[2], O);

	/*division_print(O);*/

	free(leading_v);
	free(counters);
	free(aux_array);
	free(unmoved_array);
	free(improve);
	free(indices);

	division_destroy(P);
	division_destroy(O);
	A_destroy(a);

	end = clock();
	printf("Program took: %f seconds\n", (double)(end - start) / CLOCKS_PER_SEC);

	return 0;

}


