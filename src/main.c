#include <time.h>
#include "phiADT.h"

// Local Method Declarations
void parse_cmdLineArgs(int argc, char *argv[]);

// Global variables
char *fileToRead; // name of the file to read
char *prefix;     // prefix for the output file

/*
 * Main method requires two command line arguments
 * Usage: program <file_to_process> <output_prefix>
 */
int main (int argc, char *argv[]) {

	time_t start, stop;
	time(&start);

	parse_cmdLineArgs(argc, argv);

	// open the file as read-only
	FILE *inputFile = fopen( fileToRead, "r" );
	if ( inputFile == NULL ) exit(1);

	char fileName[255];
	sprintf(fileName, "%s_distfield", prefix);

	Phi * phiFncn;
	phiFncn = phi_create(inputFile);
	phi_calc_distance_field(phiFncn);
	phi_gen_file(phiFncn, VTI|DAT, fileName);
	phi_destroy(phiFncn);

	fclose(inputFile);
	
	time(&stop);
	printf("\n-----------------------------\n");
	printf("Finished in about %f s \n", difftime(stop, start));
	printf("-----------------------------\n\n");

	return EXIT_SUCCESS;
}

/*
 * Method to parse and check the command line
 * arguments passed when running the code.
 * Arguments:
 *        int [in] - total number of arguments passed
 *   char* [] [in] - pointer to char array
 * Returns:
 */
void parse_cmdLineArgs(int argc, char *argv[]) {
	if ( argc != 3 ){
		printf("\nUsage: distance_field <file_to_process> <output_prefix>\n\n");
		exit(EXIT_FAILURE);
	}

	fileToRead = argv[1];
	prefix = argv[2];
}