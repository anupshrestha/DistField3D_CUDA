/*
 * vti_parser.c
 * This file parses information from the VTI file
 * and stores them in the appropriate output
 * parameters. Also responsible for setting the 
 * border values for the arrays.
 *
 *  Created on: Jun 04, 2015
 *      Author: Anup Shrestha
 *
 *  Audit Trail:
 *     Date: 
 *  Changes: 
 */
#include "vti_parser.h"

// private method declarations
static void move_file_pointer(FILE *file_ptr, int lineNumber, int r);
static void get_location(FILE *vti, int *l, int b_l, Grid3D g);
static void get_distance(FILE *vti, double *d, double b_d, Grid3D g);

/*
 * Creates a Grid3D struct.
 * Arguments:
 *   int x, y, z [in] - dimensions
 * Returns:
 *   Grid3D - struct representing 3D grid
 */
Grid3D make_grid3D(int x, int y, int z){
	Grid3D g;
	g.x = x; g.y = y; g.z = z;

	return g;
}

/* 
 * Moves the pointer of a FILE to specified line.
 * Argument(s):
 *   FILE* [in/out] - file pointer that is to be moved
 * 	   int [in]     - line number to move the pointer to
 * 	   int [in]     - 0 - rewind; 1 - no rewind
 * Returns:
 */
static void move_file_pointer(FILE *file_ptr, int lineNumber, int r) {
	char tmpStr[512];
	if(r) rewind(file_ptr);
	while (lineNumber > 0){
		fgets (tmpStr, 511, file_ptr);
		lineNumber--;
	}
}

/*
 * Searches for x, y, z, dx, dy, dz information
 * in the VTI file and puts it into the array
 * taken as the second argument.
 * Argument(s):
 * 	  FILE*  [in]  - vti file to be parsed
 * 	 double* [out] - array to store parsed data
 * Returns:
 */
void vti_get_dimensions(FILE *vti, double *d) {
	char tmpStr[512];
	rewind(vti);
	while (1) {
		fgets (tmpStr, 511, vti);
		if ( strstr(tmpStr, "ImageData WholeExtent") ) {
			sscanf(tmpStr, "    <ImageData WholeExtent=\"0 %lf 0 %lf 0 %lf\" Spacing=\"%lf %lf %lf\">",
					&d[0], &d[1], &d[2], &d[3], &d[4], &d[5]);
			break;
		}
	}
}

/*
 * Calls the appropriate methods to parse the
 * file for location and distance values and
 * stores them in an int and double array
 * respectively.
 * Arguments:
 * 	   FILE* [in]  - vti file to be parsed
 *      int* [out] - array to store location values
 *       int [in]  - border value for location
 * 	 double* [out] - array to store distance values
 * 	  double [in]  - border value for distance
 *   Grid3D* [in]  - dimensions of the grid
 * Returns:
 */
void vti_get_data(FILE *vti, int *l, int b_l, double *d, double b_d, Grid3D g) {
	
	// move the file pointer to
	// line 6 from beginning
	move_file_pointer(vti, 6, 1);

	get_location(vti, l, b_l, g);

	// move the file pointer 2 lines
	// forward from its last position
	move_file_pointer(vti, 2, 0);

	get_distance(vti, d, b_d, g);
	
}

/*
 * Parses the file for location values and
 * stores them in the int array. Also adds
 * the default border values for location
 * in the array.
 * Arguments:
 * 	   FILE* [in]  - vti file to be parsed
 *      int* [out] - array to store location values
 *       int [in]  - border value for location
 *   Grid3D* [in]  - dimensions of the grid
 * Returns:
 */
static void get_location(FILE *vti, int *l, int b_l, Grid3D g) {
	int i, j, k, *t = &l[0];
	for (i = 0; i < g.z; i++){
		for (j = 0; j < g.y; j++) {
			for (k = 0; k < g.x; k++) {
				// Border
				if (k == 0 || k == g.x-1 || j == 0 || j == g.y-1 || i == 0 || i == g.z-1 ) {
					*(t++) = b_l;
				}
				else{ // Interior
					fscanf(vti, "%d ", t++);
				}
			}
		}
	}
}

/*
 * Parses the file for distance values and
 * stores them in the double array. Also adds
 * the default border values for distance
 * in the array.
 * Arguments:
 * 	   FILE* [in]  - vti file to be parsed
 * 	 double* [out] - array to store distance values
 * 	  double [in]  - border value for distance
 *   Grid3D* [in]  - dimensions of the grid
 * Returns:
 */
static void get_distance(FILE *vti, double *d, double b_d, Grid3D g) {
	int i, j, k;
	double *t = &d[0];
	for (i = 0; i < g.z; i++){
		for (j = 0; j < g.y; j++) {
			for (k = 0; k < g.x; k++) {
				// Border distance
				if (k == 0 || k == g.x-1 || j == 0 || j == g.y-1 || i == 0 || i == g.z-1 ) {
					*(t++) = b_d;
				}
				else{ // Interior distance
					fscanf(vti, "%lf ", t++);
				}
			}
		}
	}
}