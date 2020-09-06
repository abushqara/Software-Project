/*
 * division.c
 *
 *  Created on: Aug 13, 2020
 *  	Author: Muhammad Watad && Shadi Abu Shqara
 */


#include <stdlib.h>
#include "error.h"
#include "division.h"



Division* division_create() {

	Division* division = (Division*)malloc(sizeof(Division));
	exit_if_no_memory(division);

	division->first = NULL;
	division->last = NULL;
	division->size = 0;

	return division;

}

void division_destroy(Division *division) {

	Node *ptr, *toBeDeleted;

	if (division == NULL)
		return;

	ptr = division->first;

	while (ptr != NULL) {
		toBeDeleted = ptr;
		ptr = ptr->next;
		B_destroy(toBeDeleted->b);
		free(toBeDeleted);
	}

	free(division);

}

/* add element to the end */
void division_add_group(Division *division, B *b) {

	Node *node = (Node*)malloc(sizeof(Node));
	exit_if_no_memory(node);

	node->b = b;
	node->next = NULL;

	if (division->size == 0) {

		division->first = node;
		division->last = node;

	}

	else {

		division->last->next = node;
		division->last = node;

	}

	division->size++;
}


B* division_extract_group(Division *division) {

	B *b;
	Node *node;

	if (division->size == 0)
		return NULL;

	node = division->first;
	division->first = division->first->next;
	division->size--;

	if (division->size == 0)
		division->last = NULL;

	b = node->b;

	free(node);

	return b;

}

void division_print(Division *division) {

	Node *ptr = division->first;
	int *g, n_g, n;
	register int i, j;

	n = division->size;

	printf("n = %d\n", n);

	printf("====================\n");

	for (i = 0; i < n; i++) {

		printf("group[%d]: Size = %d\n---\n", i, ptr->b->n_g);

		g = ptr->b->g;
		n_g = ptr->b->n_g;

		for (j = 0; j < n_g; j++) {
			printf("%d\n", *g);
			g++;
		}

		printf("====================\n");

		ptr = ptr->next;

	}

}




