#include <stdio.h>
#include <math.h>
#include "Utilities.h"

// private methods only used within the 
// this file
static void get_location(FILE *vti, int *l, Phi3D_Dim sz);
static void get_distance(FILE *vti, int *l, double *d, Phi3D_Dim sz);
static void move_file_pointer(FILE *file_ptr, int lineNumber, int r);

/*
 * Searches for x, y, z, dx, dy, dz information
 * in the VTI file and puts it into the array
 * taken as the second argument.
 * Argument(s):
 * 	 FILE *vti [in]: pointer to a file
 * 	double *d [out]: pointer to an array of double
 * Returns:
 */
void VTI_Get_dimensions(FILE *vti, double *d) {
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
 * Prints for x, y, z, dx, dy, dz information
 * represented in the double array taken in as
 * a parameter.
 * Argument(s):
 * 	const double *d [in]: pointer to an array of double
 * Returns:
 */
void VTI_Print_dimensions(const double *d) {
	printf("\nDimensions: \n");
	printf("x: %f\tdx: %f\n", d[0], d[3]);
	printf("y: %f\tdy: %f\n", d[1], d[4]);
	printf("z: %f\tdz: %f\n", d[2], d[5]);
	printf("-----------------------------\n\n");
}

void VTI_Get_locAndDist(FILE *vti, int *l, double *d, Phi3D_Dim sz) {
	
	// move the file pointer to
	// line 6 from beginning
	move_file_pointer(vti, 6, 1);

	get_location(vti, l, sz);

	// move the file pointer to 2 lines
	// forward from its last position
	move_file_pointer(vti, 2, 0);

	get_distance(vti, l, d, sz);
	
}

static void get_location(FILE *vti, int *l, Phi3D_Dim sz) {
	int i, j, k, index;
	for (i = 0; i < sz.z; i++){
		for (j = 0; j < sz.y; j++) {
			for (k = 0; k < sz.x; k++) {
				index = linear3dIndex(k, j, i, sz.x, sz.y);
				if (k == 0 || k == sz.x-1 || j == 0 || j == sz.y-1 || i == 0 || i == sz.z-1 ) {
					l[index] = -1;
				}
				fscanf(vti, "%d ", &l[index]);
			}
		}
	}
}

static void get_distance(FILE *vti, int *l, double *d, Phi3D_Dim sz) {
	int i, j, k, index;
	for (i = 0; i < sz.z; i++){
		for (j = 0; j < sz.y; j++) {
			for (k = 0; k < sz.x; k++) {
				index = linear3dIndex(k, j, i, sz.x, sz.y);
				// Border distance
				if (k == 0 || k == sz.x-1 || j == 0 || j == sz.y-1 || i == 0 || i == sz.z-1 ) {
					d[index] = INFINITY;
				}
				else{ // Interior distance
					double temp;
					fscanf(vti, "%lf ", &temp);
					if(l[index] == 1 && d[index] == INFINITY) {
						d[index] = -1;
					}
					else if ( (temp > 0.0) || ( temp < 0.0) ) {
						d[index] = temp;
					}
					else{
						d[index] = 90000;
					}
				}
			}
		}
	}
}

/* 
 * Moves the pointer of a FILE to specified line.
 * Argument(s):
 * 	FILE *file_ptr [in/out]: pointer to a file
 * 	    int lineNumber [in]: line number to move the pointer to
 * 	             int r [in]: 0 - rewind; 1 - no rewind
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
