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

#define SPEED 1

extern void setupAndRun_phiFSM(Phi* phiFncn, double *dims, FILE *vti);

char *fileToRead;
char fileToWrite1[255]; // distance field in NetCDF format
char fileToWrite2[255]; // distance field in VTI format
char fileToWrite3[255]; // distance field in regular file format

int main (int argc, char *argv[]) {

	time_t start, stop;
	time(&start);

	parse_cmdLineArgs(argc, argv);

	FILE * vti = fopen( argv[1], "r" );
	if ( vti == NULL ) printMsgAndExit("Error opening file!\n");

	// parses the dimension data from the file
	// and stores it in the dims pointer
	double* dims = (double *) malloc(sizeof(double) * 6);
	getDimensions(vti, dims);

	// print dimensions along with the spacing
	printf("\nDimensions: \n");
	printf("x: %f\tdx: %f\n", dims[0], dims[3]);
	printf("y: %f\tdy: %f\n", dims[1], dims[4]);
	printf("z: %f\tdz: %f\n", dims[2], dims[5]);
	printf("-----------------------------\n\n");


	Phi *phiFncn = create_phiFncn(dims);
	setupAndRun_phiFSM(phiFncn, dims, vti);
	
	printf("-----------------------------\n\n");
	
	writeToNetCDF(phiFncn, fileToWrite1);
	writeToFile(phiFncn, fileToWrite2);
	writeToFile(phiFncn, fileToWrite3);
	
	destroy_phiFncn(phiFncn);
	free(dims);
	fclose(vti);

	time(&stop);
	printf("\n-----------------------------\n");
	printf("Finished in about %f s \n", difftime(stop, start));
	printf("-----------------------------\n\n");
}

void parse_cmdLineArgs(int argc, char *argv[]) {
	if ( argc != 3 )
		printMsgAndExit("\nUsage: distance_field <file_to_process> <output_prefix>\n\n");

	fileToRead = argv[1];
	char *prefix = argv[2];
	sprintf(fileToWrite1, "%s_distfield.nc", prefix);
	sprintf(fileToWrite2, "%s_distfield.vti", prefix);
	sprintf(fileToWrite3, "%s_distfield.dat", prefix);
}

/* Creates the phi function for fast sweeping method
 * Initializes all variables for the phi function
 * and allocates memory required by the location and
 * distance arrays.
 * Argument(s):
 * 	double *dims: Pointer to a double with the
 * 	              dimension data
 * Returns:
 * 	Phi* phiFncn: pointer to the Phi (phi function)
 */
Phi* create_phiFncn(double *dims) {

	Phi *phiFncn = ( Phi * ) malloc ( sizeof( Phi ) );

	// initialize the variables
	phiFncn->F = SPEED; // speed of propagation [Equation (1) in the report]
	phiFncn->size_x = (int) dims[0] + 1; phiFncn->dx = dims[3];
	phiFncn->size_y = (int) dims[1] + 1; phiFncn->dy = dims[4];
	phiFncn->size_z = (int) dims[2] + 1; phiFncn->dz = dims[5];

	int totalGridNodes = (phiFncn->size_x + 2) * (phiFncn->size_y + 2) * (phiFncn->size_z + 2);

	// allocating host memory to store location
	// and distance data
	phiFncn->location = (int *) malloc( sizeof(int) * totalGridNodes );
	phiFncn->distance = (double *) malloc( sizeof(double) * totalGridNodes );
	return phiFncn;
}

/* Deallocates all memory allocated by
 * the create_phiFncn()
 * Argument(s):
 * 	Phi* phiFncn: pointer to the Phi (phi function)
 * 	Returns:
 */
void destroy_phiFncn(Phi *phiFncn) {
	free(phiFncn->location);
	free(phiFncn->distance);
	free(phiFncn);
}

/*
 * Searches for x, y, z, dx, dy, dz information
 * in the VTI file and puts it into the array
 * taken as the second argument.
 * Argument(s):
 * 	   FILE *vti: pointer to the file
 * 	double *dims: pointer to an array of double
 * Returns:
 */
void getDimensions(FILE *vti, double *dims) {
	char tmpStr[512];
	rewind(vti);
	while (1) {
		fgets (tmpStr, 511, vti);
		if ( strstr(tmpStr, "ImageData WholeExtent") ) {
			sscanf(tmpStr, "    <ImageData WholeExtent=\"0 %lf 0 %lf 0 %lf\" Spacing=\"%lf %lf %lf\">",
					&dims[0], &dims[1], &dims[2], &dims[3], &dims[4], &dims[5]);
			break;
		}
	}
}
