#include <stdio.h>
#include <netcdf.h>
#include "Phi3D_Writer.h"

static write_to_dat(Phi3D *phiFncn, Grid3D *g3d, char *fileName);
static write_to_vti(Phi3D *phiFncn, Grid3D *g3d, char *fileName);
static write_to_ncf(Phi3D *phiFncn, Grid3D *g3d, char *fileName);

static int linear3dIndex(int x, int y, int z, int max_x, int max_y) {
  return z * max_y * max_x + y * max_x + x;
}

void write_to_file(Phi3D *phiFncn, Grid3D *g3d, FileType type, char *prefix) {

	char fileName[255];

	// check if DAT bit is set
	if(type & DAT) {
		sprintf(fileName, "%s_distfield.dat", prefix);
		write_to_dat(phiFncn, g3d, fileName);
		printf("VTI file created: %s\n", fileName);
	}

	// check if VTI bit is set
	if(type & VTI) {
		sprintf(fileName, "%s_distfield.vti", prefix);
		write_to_vti(phiFncn, g3d, fileName);
		printf("DAT file created: %s\n", fileName);
	}

	// check if NCF bit is set
	if(type & NCF) {
		sprintf(fileName, "%s_distfield.nc", prefix);
		write_to_ncf(phiFncn, g3d, fileName);
		printf("NCF file created: %s\n", fileName);
	}

} 

static write_to_dat(Phi3D *phiFncn, Grid3D *g3d, char *fileName) {
	write_to_vti(phiFncn, g3d, fileName);
}


static write_to_vti(Phi3D *phiFncn, Grid3D *g3d, char *fileName) {
	
	int SIZE_X = phiFncn->x; double dx = phiFncn->dx;
	int SIZE_Y = phiFncn->y; double dy = phiFncn->dy;
	int SIZE_Z = phiFncn->z; double dz = phiFncn->dz;

	FILE *fp = fopen(fileName, "w");

	fprintf(fp, "<?xml version=\"1.0\"?>\n");
	fprintf(fp, "<VTKFile type=\"ImageData\" version=\"0.1\" byte_order=\"LittleEndian\">\n");
	fprintf(fp, "\t<ImageData WholeExtent=\"0 %d 0 %d 0 %d\" Origin=\"0 0 0\" Spacing=\" %f %f %f  \">\n", SIZE_X+1, SIZE_Y+1, SIZE_Z+1, dx, dy, dz);
	fprintf(fp, "\t\t<Piece Extent=\"0 %d 0 %d 0 %d\">\n", SIZE_X+1, SIZE_Y+1, SIZE_Z+1);
	fprintf(fp, "\t\t\t<PointData Scalars=\"Distance Field\">\n");
	fprintf(fp, "\t\t\t\t<DataArray type=\"Float32\" Name=\"Distance Field\" format=\"ascii\">\n");

	int i,j,k;
	for(i = 0; i < g3d->z; i++){
		for(j = 0; j < g3d->y; j++){
			for(k = 0; k < g3d->x; k++) {
			  fprintf(fp, "%f  ", phiFncn->distance[linear3dIndex(k, j, i, g3d->x, g3d->y)]);
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
}


static write_to_ncf(Phi3D *phiFncn, Grid3D *g3d, char *fileName) {
	
	int ndims = 3;
  
  // dimensions    // spacing 
  int nx = g3d->x; double dx = phiFncn->dx;
  int ny = g3d->y; double dy = phiFncn->dy;
  int nz = g3d->z; double dz = phiFncn->dz;
  
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
}