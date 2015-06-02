/*
 * DistField3D.c
 * Main module that calculates in parallel using CUDA,
 * the distance field for an interface in a given VTI
 * file using the
 *     - Fast Sweeping Method
 * The algorithm implemented for parallel fast sweeping
 * method is from a paper in the Journal of Computational
 * Physics titled "A parallel fast sweeping method for
 * the Eikonal equation" by Miles Detrixhe, Federic Gibou,
 * and Chohong Min.
 * DOI: http://www.sciencedirect.com/science/article/pii/S002199911200722X
 *
 *  Created on: May 14, 2015
 *      Author: Anup Shrestha
 *
 *  Audit Trail:
 *     Date: May 19, 2015
 *  Changes: Creating NetCDF file format - A.S.
 *
 */

#include <time.h>
#include "Phi3D.h"
#include "Utilities.h"
#include "VTI_Parser.h"

#define SPEED 1

char *fileToRead;
char fileToWrite1[255]; // distance field in NetCDF format
char fileToWrite2[255]; // distance field in VTI format
char fileToWrite3[255]; // distance field in regular file format

static void parse_cmdLineArgs(int argc, char *argv[]) {
	if ( argc != 3 )
		printMsgAndExit("\nUsage: distance_field <file_to_process> <output_prefix>\n\n");

	fileToRead = argv[1];
	char *prefix = argv[2];
	sprintf(fileToWrite1, "%s_distfield.nc", prefix);
	sprintf(fileToWrite2, "%s_distfield.vti", prefix);
	sprintf(fileToWrite3, "%s_distfield.dat", prefix);
}

int main (int argc, char *argv[]) {

	time_t start, stop;
	time(&start);

	parse_cmdLineArgs(argc, argv);

	// opening the file as read-only
	FILE *inputFile = fopen( fileToRead, "r" );
	if ( inputFile == NULL ) printMsgAndExit("Error opening file!\n");

	// fetching and printing the dimension
	// and spacing values from the input
	// file
	double *dims = (double *) malloc(sizeof(double) * 6);
	VTI_Get_dimensions(inputFile, dims);
	VTI_Print_dimensions(dims);

	Phi *phiFncn = Phi3D_Init(dims);
	Phi3D_Calc_distField(phiFncn, inputFile);
	
	printf("-----------------------------\n\n");
	
	writeToNetCDF(phiFncn, fileToWrite1);
	writeToFile(phiFncn, fileToWrite2);
	writeToFile(phiFncn, fileToWrite3);
	
	Phi3D_Finalize(phiFncn);

	free(dims);
	fclose(inputFile);

	time(&stop);
	printf("\n-----------------------------\n");
	printf("Finished in about %f s \n", difftime(stop, start));
	printf("-----------------------------\n\n");
}


