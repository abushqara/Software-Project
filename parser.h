/*
 * parser.h
 *
 *  Created on: Aug 13, 2020
 *  	Author: Muhammad Watad && Shadi Abu Shqara
 */



#ifndef PARSER_H_
#define PARSER_H_

#define _CRT_SECURE_NO_WARNINGS
#include "division.h"

B* parser_read_input(char *fileName);

void parser_write_output(char *fileName, Division *division);


#endif /* PARSER_H_ */
