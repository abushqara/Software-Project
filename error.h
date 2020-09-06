/*
 * error.h
 *
 *  Created on: Aug 13, 2020
 *  	Author: Muhammad Watad && Shadi Abu Shqara
 */

#ifndef ERROR_H_
#define ERROR_H_

#include <stdio.h>

void exit_if_no_memory(void *ptr);
void exit_if_file_failure(FILE *file);
void exit_if_read_write_failed(int a, int b);
void exit_if_dividing_by_zero(double m);
void exit_insufficient_arguments();

#endif /* ERROR_H_ */
