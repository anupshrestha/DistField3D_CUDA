/*
 * file_writer.c
 * Generates file of different types based
 * on the information provided through the
 * FileType structure. The types of file that
 * can be generated are:
 *  1) VTI - VTK Image Data
 *  2) DAT - Data file
 *  3) NCF - NetCDF file
 *
 *  Created on: Jun 04, 2015
 *      Author: Anup Shrestha
 *
 *  Audit Trail:
 *     Date: 
 *  Changes:
 */
#include "file_writer.h"

// private method declarations
static void file_gen_dat(double *data);
static void file_gen_vti(double *data);
static void file_gen_ncf(double *data);

// global variables
static char fileName[255]; // stores the full filename
static int x, y, z;        // grid dimensions
static double dx, dy, dz;  // node spacing

/*
 * Creates a FileGrid struct.
 * Arguments:
 *         int x, y, z [in] - 3D dimension of grid
 *   double dx, dy, dz [in] - spacing of the nodes on
 *                            each dimension of grid
 * Returns:
 *   FileGrid - struct with dimension and spacing 
 *              information
 */
FileGrid make_fileGrid(int x, int y, int z, double dx, double dy, double dz) {
	FileGrid g;
	g.x = x; g.dx = dx;
	g.y = y; g.dy = dy;
	g.z = z; g.dz = dz;

	return g;
}

/*
 * Creates a FileType struct.
 * Arguments:
 *      char* [in] - name of the file
 *    double* [in] - data to be written
 *    FileOut [in] - type of file to generate
 *   FileGrid [in] - dimensions/spacing of the grid
 * Returns:
 *   FileType - struct containing all the information
 *              to generate a file
 */
FileType make_fileType( char *n, double *d, FileOut o, FileGrid g) {
	FileType f;
	f.name = n;
	f.data = d;
	f.out  = o;
	f.grid = g;

	return f;
}

/*
 * Determines what types of file to output
 * generates filename for each type and 
 * calls the appropriate methods to generate
 * the file.
 * Arguments:
 *   FileType* [in] - pointer to FileType
 * Returns:
 */
void file_generate(FileType *f) {
	
	// initialize global variables
	x = f->grid.x; dx = f->grid.dx;
	y = f->grid.y; dy = f->grid.dy;
	z = f->grid.z; dz = f->grid.dz;

	// check if DAT bit is set
	if(f->out & DAT) {
		sprintf(fileName, "%s.dat", f->name);
		file_gen_dat(f->data);
		printf("DAT file created: %s\n", fileName);
	}

	// check if VTI bit is set
	if(f->out & VTI) {
		sprintf(fileName, "%s.vti", f->name);
		file_gen_vti(f->data);
		printf("VTI file created: %s\n", fileName);
	}

	// check if NCF bit is set
	if(f->out & NCF) {
		sprintf(fileName, "%s.nc", f->name);
		file_gen_ncf(f->data);
		printf("NCF file created: %s\n", fileName);
	}
}

/*
 * Generates .dat output file
 * Arguments:
 *   double* - data to be written
 * Returns:
 */
static void file_gen_dat(double *data) {
	file_gen_vti(data);
}

/*
 * Generates .vti output file
 * Arguments:
 *   double* - data to be written
 * Returns:
 */
static void file_gen_vti(double *data) {
	
	FILE *fp = fopen(fileName, "w");
	
	// Header Information + Data for VTI file
	fprintf(fp, "<?xml version=\"1.0\"?>\n");
	fprintf(fp, "<VTKFile type=\"ImageData\" version=\"0.1\" byte_order=\"LittleEndian\">\n");
	fprintf(fp, "\t<ImageData WholeExtent=\"0 %d 0 %d 0 %d\" Origin=\"0 0 0\" Spacing=\" %f %f %f  \">\n", x-1, y-1, z-1, dx, dy, dz);
	fprintf(fp, "\t\t<Piece Extent=\"0 %d 0 %d 0 %d\">\n", x-1, y-1, z-1);
	fprintf(fp, "\t\t\t<PointData Scalars=\"Distance Field\">\n");
	fprintf(fp, "\t\t\t\t<DataArray type=\"Float32\" Name=\"Distance Field\" format=\"ascii\">\n");

	int i,j,k;
	double *t = &data[0];
	for(i = 0; i < z; i++){
		for(j = 0; j < y; j++){
			for(k = 0; k < x; k++) {
			  fprintf(fp, "%f  ", *(t++));
			}
			fprintf(fp, "\n");
		}
		fprintf(fp, "\n");
	}

	// Closing matching entities
	fprintf(fp, "\t\t\t\t</DataArray>\n");
	fprintf(fp, "\t\t\t</PointData>\n");
	fprintf(fp, "\t\t\t<CellData>\n");
	fprintf(fp, "\t\t\t</CellData>\n");
	fprintf(fp, "\t\t</Piece>\n");
	fprintf(fp, "\t</ImageData>\n");
	fprintf(fp, "</VTKFile>\n");

	fclose(fp);
}

/*
 * Generates .nc output file
 * Arguments:
 *   double* - data to be written
 * Returns:
 */
static void file_gen_ncf(double *data){
	
	int ndims = 3;
  
  // IDs for the netCDF file, dimensions, and variables
  int ncid, x_dimid, y_dimid, z_dimid, retval;
  int dist_varid, x_varid, y_varid, z_varid;
  int dimids[ndims];
  
  // Data for x, y, z variables
  float lats[y], lons[x], deps[z];
  int i,j,k;
  for(k = 0; k < x; k++)
    lons[k] = dx * k;
  for(j = 0; j < y; j++)
    lats[j] = dy * j;
  for(i = 0; i < z; i++)
    deps[i] = dz * i;
  
  // Create the file. The NC_NETCDF4 flag tells netCDF to
  // create a netCDF-4/HDF5 file.
  if((retval = nc_create(fileName, NC_NETCDF4  | NC_CLOBBER, &ncid)))
    ERR(retval);

  // Define the dimensions in the root group. Dimensions are visible
  // in all subgroups
  if ((retval = nc_def_dim(ncid, "x", x, &x_dimid)))
    ERR(retval);
  
  if ((retval = nc_def_dim(ncid, "y", y, &y_dimid)))
    ERR(retval);
  
  if ((retval = nc_def_dim(ncid, "z", z, &z_dimid)))
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
  
  if ((retval = nc_put_var_double(ncid, dist_varid, data)))
    ERR(retval);
  
  // Close the file.
  if ((retval = nc_close(ncid)))
    ERR(retval);
}