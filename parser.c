/*
 * parser.c
 *
 *  Created on: Aug 13, 2020
 *  	Author: Muhammad Watad && Shadi Abu Shqara
 */


#include "parser.h"
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include "error.h"



int num_of_non_zero_entries(char *fileName, int n) {

	int num_of_ints;
	struct stat st;

	stat(fileName, &st);

	num_of_ints = st.st_size / sizeof(int);

	return num_of_ints - n - 1;

}


/* Should we handle n = 0, n = 1... ? */
B* parser_read_input(char *fileName) {

	int n, res, rank_of_i, *neighbors_i, m = 0, *k_array, nnz;
	register int i;
	FILE *input_file;
	A *a;
	B *b;

	input_file = fopen(fileName, "rb");
	exit_if_file_failure(input_file);

	res = fread(&n, sizeof(int), 1, input_file);
	exit_if_read_write_failed(res, 1);

	nnz = num_of_non_zero_entries(fileName, n);

	k_array = (int*)malloc(sizeof(int) * n);
	exit_if_no_memory(k_array);

	a = A_create(n, nnz);

	neighbors_i = (int*)malloc(sizeof(int) * n);
	exit_if_no_memory(neighbors_i);

	for (i = 0; i < n; i++) {

		res = fread(&rank_of_i, sizeof(int), 1, input_file);
		exit_if_read_write_failed(res, 1);

		m += rank_of_i;
		k_array[i] = rank_of_i;

		res = fread(neighbors_i, sizeof(int), rank_of_i, input_file);
		exit_if_read_write_failed(res, rank_of_i);

		A_add_row(a, neighbors_i, rank_of_i, i);

		memset(neighbors_i, 0, n * sizeof(int));

	}

	a->m = m;
	a->k = k_array;

	b = B_create(a, n);
	B_initialize_g(b);

	free(neighbors_i);
	fclose(input_file);

	return b;
}


void parser_write_output(char *fileName, Division *division) {

	Node *ptr = division->first;
	int size = division->size, res, n_g, *g;
	FILE *output_file;

	output_file = fopen(fileName, "wb");
	exit_if_file_failure(output_file);

	res = fwrite(&size, sizeof(int), 1, output_file);
	exit_if_read_write_failed(res, 1);

	while (ptr != NULL) {

		n_g = ptr->b->n_g;

		res = fwrite(&n_g, sizeof(int), 1, output_file);
		exit_if_read_write_failed(res, 1);

		g = ptr->b->g;

		res = fwrite(g, sizeof(int), n_g, output_file);
		exit_if_read_write_failed(res, n_g);

		ptr = ptr->next;

	}

	fclose(output_file);

}

