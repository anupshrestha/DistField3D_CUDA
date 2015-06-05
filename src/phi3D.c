/*
 * phi3D.c
 * This is a wrapper file that calls appropriate
 * methods from other different files that:
 *   1) Creates and initializes the phi function
 *      variables
 *   2) Calculates the distance field
 *   3) Writes the data to file
 *   4) Performs cleanup
 *
 *  Created on: Jun 04, 2015
 *      Author: Anup Shrestha
 *
 *  Audit Trail:
 *     Date: 
 *  Changes: 
 */
#include "phi3D.h"

// private method declarations
static void adjust_boundary(Phi *p, Grid3D gridEx);
static int linear3dIndex(int x, int y, int z, int max_x, int max_y);


// external method declarations
extern void init_phiFncn(Phi *p, double *d);
extern void free_phiFncn(Phi *p);
extern void calc_distfield(Phi *p, Grid3D g3d);


/*
 * Creates phi function, calls the file parser
 * methods to get the dimensions, location
 * and distance data from the file. Also,
 * allocates memory and initializes variables
 * required by phi function.
 * Arguments:
 *     FILE* [in] - pointer to a file
 * Returns:
 *   Phi* - pointer to phi function
 */
Phi* phi_create(FILE *vti) {
	
	Phi *p = (Phi *) malloc(sizeof(Phi));

	// parsing the information about the dimensions
	// and the spacing of the nodes from the file and
	// storing them in an array
	double *dims = (double *) malloc(sizeof(double) * 6);
	vti_get_dimensions(vti, dims);

	// now using the information about the dimensions
	// and the spacing create a phi function by
	// initializing all the attributes and allocating
	// space for the location and distance data
	init_phiFncn(p, dims);

	// the phi function stores the internal dimensions
	// so store the external dimensions in a global
	// Grid3D struct variable for future use.
	// Grid3D and make_grid3D are defined in vti_parser.h
	Grid3D gridEx = make_grid3D(p->x + 2, p->y + 2, p->z + 2);

	// now parse the location and distance data from 
	// the file and store it in appropriate arrays
	// inside the phi function.
	vti_get_data(vti, p->location, DEFAULT_BORDER_LOCATION,
		                p->distance, DEFAULT_BORDER_DISTANCE, gridEx);

	return p;
}

/*
 * Calls the method that deallocates
 * memory allocated for the phi function.
 * Arguments:
 *    Phi* [in] - pointer to phi function
 * Returns:
 */
void phi_destroy(Phi *p) {

	free_phiFncn(p);

}

/*
 * Calls the method that calculates the distance
 * field.
 * Arguments:
 *    Phi* [in/out] - pointer to phi function
 * Returns:
 */
void phi_calc_distance_field(Phi *p) {

	Grid3D gridEx = make_grid3D(p->x + 2, p->y + 2, p->z + 2);
	calc_distfield(p, gridEx);

}


/*
 * Updates the boundary values of the distance
 * array, then calls the method that writes
 * distance field to different files based on
 * the FileOut parameter
 * Arguments:
 *      Phi* [in] - pointer to phi function
 *   FileOut [in] - Output file type enumerator
 *     char* [in] - name of the output file
 * Returns:
 */
void phi_gen_file(Phi *p, FileOut out, char *fileName) {

  Grid3D gridEx = make_grid3D(p->x + 2, p->y + 2, p->z + 2);
  adjust_boundary(p, gridEx);
  FileGrid fg = make_fileGrid(gridEx.x, gridEx.y, gridEx.z, p->dx, p->dy, p->dz);
  FileType ft = make_fileType(fileName, p->distance, out, fg);
  
  // generate file(s)
  file_generate(&ft);

}

/*
 * Adjusts the boundary values
 * Arguments:
 *     Phi* [in/out] - pointer to phi function
 *   Grid3D [in]     - dimensions of grid
 * Returns:
 */
static void adjust_boundary(Phi *p, Grid3D gridEx) {
	// before writing to file we need to update the 
  // boundary values
  int x, y, z, i, j, k;
  x = gridEx.x;
  y = gridEx.y;
  z = gridEx.z;

  for(i = 0; i < z; i++){
    for(j = 0; j < y; j++){
      for(k = 0; k < x; k++){
      	int I = i, J = j, K = k;
      	I = (i == z-1) ? I-1 : (!i) ? I+1 : I;
      	J = (j == y-1) ? J-1 : (!j) ? J+1 : J;
      	K = (k == x-1) ? K-1 : (!k) ? K+1 : K;
      	if( i != I || j != J || k != K) {
      		p->distance[linear3dIndex(k, j, i, x, y)] = p->distance[linear3dIndex(K, J, I, x, y)];
        }
      }
    }
  }
}

/*
 * Convert 3D indexing to 1D indexing 
 * Arguments:
 *   int x, y, z [in] - 3D coordinate
 *   int max_x   [in] - size of x-dimension
 *   int max_y   [in] - size of y-dimension
 * Returns:
 */
static int linear3dIndex(int x, int y, int z, int max_x, int max_y) {
  return z * max_y * max_x + y * max_x + x;
}