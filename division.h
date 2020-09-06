/*
 * division.h
 *
 *  Created on: Aug 13, 2020
 *  	Author: Muhammad Watad && Shadi Abu Shqara
 */

#ifndef DIVISION_H_
#define DIVISION_H_

#include "b.h"

struct node_t {
	B *b;
	struct node_t *next;
};

typedef struct node_t Node;

typedef struct division_t {

	Node *first;
	Node *last;

	int size;

} Division;


Division* division_create();

void division_destroy(Division *division);

void division_add_group(Division *division, B *b);

B* division_extract_group(Division *division);

void division_print(Division *division);


#endif /* DIVISION_H_ */
#pragma once
