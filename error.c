/*
 * error.c
 *
 *  Created on: Aug 13, 2020
 *  	Author: Muhammad Watad && Shadi Abu Shqara
 */

#include <stdio.h>
#include <stdlib.h>

void exit_if_no_memory(void *ptr) {

	if (ptr == NULL) {
		printf("Error: memory allocation failed!\n");
		exit(-1);
	}

}


void exit_if_file_failure(FILE *file) {

	if (file == NULL) {
		printf("Error: can't open a file!\n");
		exit(-1);
	}

}

void exit_if_read_write_failed(int a, int b) {

	if (a != b) {
		printf("Error: I/O operation failed!\n");
		exit(-1);
	}

}

void exit_if_dividing_by_zero(double m) {

	if (m == 0) {
		printf("Error: can't divide by zero!\n");
		exit(-1);
	}

}

void exit_insufficient_arguments() {

	printf("Error: insufficient number of arguments were supplied!\n");
	exit(-1);

}
