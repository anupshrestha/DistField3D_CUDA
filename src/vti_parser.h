/*
 * vti_parser.h
 *
 *  Created on: Jun 04, 2015
 *      Author: Anup Shrestha
 *
 *  Audit Trail:
 *     Date: 
 *  Changes: 
 */

#ifndef VTI_Parser_H_
#define VTI_Parser_H_

#include <stdio.h>
#include <string.h>

/*
 * Structure for storing the dimension
 * and spacing of the grid
 */
typedef struct {
	int    x, y, z;
} Grid3D;

// public method declarations
Grid3D make_grid3D(int x, int y, int z);
void vti_get_dimensions(FILE *vti, double *d);
void vti_get_data(FILE *vti, int *l, int b_l, double *d, double b_d, Grid3D g);

#endif /* VTI_Parser_H_ */
