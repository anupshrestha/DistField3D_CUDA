/*
 * Utilities.c
 * Contains methods that are used throughout
 * the distance field calculation.
 *
 *  Created on: May 14, 2015
 *      Author: Anup Shrestha
 *
 *  Audit Trail:
 *     Date: May 19, 2015
 *  Changes: Added a new function writeToNetCDF, that
 *           write a file in NetCDF format - A.S.
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <netcdf.h>
#include "Utilities.h"

/* Calculates the required size of the grid for
 * launching CUDA kernels. This method is used
 * during the creation of blocks per grid before
 * the launch of any CUDA kernel.
 * Argument(s):
 * 	int a: total number of threads needed
 * 	       in a single dimension
 * 	int b: total number of threads per block
 * 	       in a single dimension
 * Returns:
 * 	int: Number of grids required to launch the
 * 	     kernel with specified number of total
 * 	     threads in a single dimension
 */
int iDivUp(int a, int b) {
	return ( (a % b) != 0 ) ? (a / b + 1) : (a / b);
}

/* Converts 3D indexing to 1D indexing.
 * Example : 0,0,0 - 0
 *           0,0,1 - 1
 * Argument(s):
 * 	int x, y, z: 3D-coordinate
 * 	  int max_x: maximum x-coordinate
 * 	  int max_y: maximum y-coordinate
 * Returns:
 * 	int: Linear 1D index of the 3D index
 */
int linear3dIndex(int x, int y, int z, int max_x, int max_y) {
  return z * max_y * max_x + y * max_x + x;
}


/* Prints a message to the console
 * and exits the program.
 * Argument(s):
 * 	char *msg: message to be printed
 * Returns:
 */
void printMsgAndExit(char *msg) {
	printf(msg); exit(EXIT_FAILURE);
}

/* Creates a file with the specified filename
 * and writes the data from the distance array
 * to the file.
 * Argument(s):
 * 	  Phi *phiFncn: pointer to phi
 * 	char *fileName: name of the file
 * Returns:
 */
void writeToFile(Phi *phiFncn, char *fileName) {

	int SIZE_X = phiFncn->size_x; double dx = phiFncn->dx;
	int SIZE_Y = phiFncn->size_y; double dy = phiFncn->dy;
	int SIZE_Z = phiFncn->size_z; double dz = phiFncn->dz;

	FILE *fp = fopen(fileName, "w");

	fprintf(fp, "<?xml version=\"1.0\"?>\n");
	fprintf(fp, "<VTKFile type=\"ImageData\" version=\"0.1\" byte_order=\"LittleEndian\">\n");
	fprintf(fp, "\t<ImageData WholeExtent=\"0 %d 0 %d 0 %d\" Origin=\"0 0 0\" Spacing=\" %f %f %f  \">\n", SIZE_X+1, SIZE_Y+1, SIZE_Z+1, dx, dy, dz);
	fprintf(fp, "\t\t<Piece Extent=\"0 %d 0 %d 0 %d\">\n", SIZE_X+1, SIZE_Y+1, SIZE_Z+1);
	fprintf(fp, "\t\t\t<PointData Scalars=\"Distance Field\">\n");
	fprintf(fp, "\t\t\t\t<DataArray type=\"Float32\" Name=\"Distance Field\" format=\"ascii\">\n");

	int iend = SIZE_Z + 2;
	int jend = SIZE_Y + 2;
	int kend = SIZE_X + 2;
	int i,j,k;
	for(i=0; i<iend;i++){
		for(j=0; j<jend;j++){
			for(k=0; k<kend; k++) {
			        int I = i, J = j, K = k;
				I = (i == iend-1) ? I-1 : (!i) ? I+1 : I;
				J = (j == jend-1) ? J-1 : (!j) ? J+1 : J;
				K = (k == kend-1) ? K-1 : (!k) ? K+1 : K;
			        fprintf(fp, "%f  ", phiFncn->distance[linear3dIndex(K, J, I, kend, jend)]);
			}
			fprintf(fp, "\n");
		}
		fprintf(fp, "\n");
	}

	fprintf(fp, "\t\t\t\t</DataArray>\n");
	fprintf(fp, "\t\t\t</PointData>\n");
	fprintf(fp, "\t\t\t<CellData>\n");
	fprintf(fp, "\t\t\t</CellData>\n");
	fprintf(fp, "\t\t</Piece>\n");
	fprintf(fp, "\t</ImageData>\n");
	fprintf(fp, "</VTKFile>\n");

	fclose(fp);
	printf("VTI file created: %s\n", fileName);
}

/* Creates a NetCDF file with the specified filename
 * and writes the data from the distance array
 * to the file.
 * Argument(s):
 * 	  Phi *phiFncn: pointer to phi
 * 	char *fileName: name of the file
 * Returns:
 */
void writeToNetCDF(Phi *phiFncn, char *fileName) {
  
  int ndims = 3;
  
  // dimensions
  int nx = phiFncn->size_x + 2;
  int ny = phiFncn->size_y + 2;
  int nz = phiFncn->size_z + 2;
  
  // spacing 
  double dx = phiFncn->dx;
  double dy = phiFncn->dy;
  double dz = phiFncn->dz;
  
  // Making sure the boundaries have 
  // the correct value
  int i,j,k;
  for(i=0; i<nz;i++){
    for(j=0; j<ny;j++){
      for(k=0; k<nx; k++){
	int I = i, J = j, K = k;
	I = (i == nz-1) ? I-1 : (!i) ? I+1 : I;
	J = (j == ny-1) ? J-1 : (!j) ? J+1 : J;
	K = (k == nx-1) ? K-1 : (!k) ? K+1 : K;
	if( i != I || j != J || k != K)
	  phiFncn->distance[linear3dIndex(k, j, i, nx, ny)] = phiFncn->distance[linear3dIndex(K, J, I, nx, ny)];
      }
    }
  }
  
  // IDs for the netCDF file, dimensions, and variables
  int ncid, x_dimid, y_dimid, z_dimid, retval;
  int dist_varid, x_varid, y_varid, z_varid;
  int dimids[ndims];
  
  // Data for x, y, z variables
  float lats[ny], lons[nx], deps[nz];
  
  for(k = 0; k < nx; k++)
    lons[k] = dx * k;
  for(j = 0; j < ny; j++)
    lats[j] = dy * j;
  for(i = 0; i < nz; i++)
    deps[i] = dz * i;
  
  // Create the file. The NC_NETCDF4 flag tells netCDF to
  // create a netCDF-4/HDF5 file.
  if((retval = nc_create(fileName, NC_NETCDF4  | NC_CLOBBER, &ncid)))
    ERR(retval);

  // Define the dimensions in the root group. Dimensions are visible
  // in all subgroups
  if ((retval = nc_def_dim(ncid, "x", nx, &x_dimid)))
    ERR(retval);
  
  if ((retval = nc_def_dim(ncid, "y", ny, &y_dimid)))
    ERR(retval);
  
  if ((retval = nc_def_dim(ncid, "z", nz, &z_dimid)))
    ERR(retval);
   
  
  // The dimids passes the IDs of
  // the dimensions of the variable
  dimids[0] = z_dimid;
  dimids[1] = y_dimid;
  dimids[2] = x_dimid;
  
  // Defining the coordinate variables
  if ((retval = nc_def_var(ncid, "x", NC_FLOAT, 1, &x_dimid, &x_varid)))
      ERR(retval);
  
  if ((retval = nc_def_var(ncid, "y", NC_FLOAT, 1, &y_dimid, &y_varid)))
      ERR(retval);
  
  if ((retval = nc_def_var(ncid, "z", NC_FLOAT, 1, &z_dimid, &z_varid)))
      ERR(retval);
  
  // Assign delta attribute to coordinate variables
  if ((retval = nc_put_att_double(ncid, x_varid, "delta", NC_DOUBLE, 1, &dx)))
    ERR(retval);
  
  if ((retval = nc_put_att_double(ncid, y_varid, "delta", NC_DOUBLE, 1, &dy)))
    ERR(retval);
  
  if ((retval = nc_put_att_double(ncid, z_varid, "delta", NC_DOUBLE, 1, &dz)))
    ERR(retval);
  
  // Define a float variable in root, using dimensions
  // in the root group.
  if ((retval = nc_def_var(ncid, "Distance Field", NC_FLOAT, ndims, dimids, &dist_varid)))
      ERR(retval);
   
  // End define mode
  if ((retval = nc_enddef(ncid)))
      ERR(retval);
  
  // Put the data values for the dimensions and distance
  if ((retval = nc_put_var_float(ncid, x_varid, &lons[0])))
      ERR(retval);
  
  if ((retval = nc_put_var_float(ncid, y_varid, &lats[0])))
      ERR(retval);
  
  if ((retval = nc_put_var_float(ncid, z_varid, &deps[0])))
      ERR(retval);
  
  if ((retval = nc_put_var_double(ncid, dist_varid, phiFncn->distance)))
      ERR(retval);
  
  // Close the file.
  if ((retval = nc_close(ncid)))
      ERR(retval);
  
  printf("NetCDF file created: %s\n", fileName);
}
