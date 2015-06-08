/*
 * file_writer.h
 *
 *  Created on: Jun 04, 2015
 *      Author: Anup Shrestha
 *
 *  Audit Trail:
 *     Date: 
 *  Changes:
 */

#ifndef File_Writer_H_
#define File_Writer_H_

#include <stdio.h>
#include <stdlib.h>
#include <netcdf.h>

/*
 * Enumerator for types of output files
 *  1) VTI - VTK Image Data
 *  2) DAT - Data file
 *  3) NCF - NetCDF file
 */
typedef enum { 
	VTI = 0x01,
	DAT = 0x02,
	NCF = 0x04
} FileOut;

/*
 * Structure for storing the dimension
 * and spacing of the grid
 */
typedef struct {
	int    x, y, z;
	double dx, dy, dz;
} FileGrid;

/*
 * Structure for information required
 * to generate an output file
 */
typedef struct {
	char      *name; // name of the output file
	double    *data; // data to be put into file
	FileOut   out;   // type of file to output
	FileGrid  grid;  // grid dimensions/spacing
} FileType;
	
// NetCDF Macros
// Handle errors by printing an error message and exiting with a
// non-zero status.
#define ERRCODE 2
#define ERR(e) {printf("Error: %s\n", nc_strerror(e)); exit(ERRCODE);}
 
// public method declarations
FileGrid make_fileGrid(int x, int y, int z, double dx, double dy, double dz);
FileType make_fileType(char *n, double *d, FileOut o, FileGrid g);
void file_generate(FileType *f);

#endif /* File_Writer_H_ */
