/*
 * main.c
 * Main driver file that executes methods for calculating
 * the distance field for an interface in a given VTI
 * file using the - Fast Sweeping Method
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

void run_dfp_fsm(FILE *vti);
void parse_cmdLineArgs(int argc, char *argv[]);

char *fileToRead;
char *prefix;

int main (int argc, char *argv[]) {

	time_t start, stop;
	time(&start);

	parse_cmdLineArgs(argc, argv);

	// open the file as read-only
	FILE *inputFile = fopen( fileToRead, "r" );
	if ( inputFile == NULL ) exit(1);

	// run the distance field pre-processor
	// using the fast sweeping method
	run_dfp_fsm(inputFile);

	fclose(inputFile);

	time(&stop);
	printf("\n-----------------------------\n");
	printf("Finished in about %f s \n", difftime(stop, start));
	printf("-----------------------------\n\n");

	return EXIT_SUCCESS;
}

void run_dfp_fsm(FILE *vti) {

	Grid3D *g3d = (Grid3D *) malloc(sizeof(Grid3D));
	
	Phi3D *phiFncn = Phi3D_Init(vti, g3d);
	
	Phi3D_Calc_distfield(phiFncn, g3d);
	
	Phi3D_Write(phiFncn, g3d, VTI | DAT, prefix);

	Phi3D_Finalize(phiFncn);	
	
	free(g3d);

}

void parse_cmdLineArgs(int argc, char *argv[]) {
	if ( argc != 3 ){
		printf("\nUsage: distance_field <file_to_process> <output_prefix>\n\n");
		exit(EXIT_FAILURE);
	}

	fileToRead = argv[1];
	prefix = argv[2];
}