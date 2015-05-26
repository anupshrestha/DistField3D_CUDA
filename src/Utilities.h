/*
 * Utilities.h
 * Contains methods that are used throughout
 * the distance field calculation.
 *
 *  Created on: May 14, 2015
 *      Author: Anup Shrestha
 */

#ifndef UTILITIES_H_
#define UTILITIES_H_
	int iDivUp(int a, int b);
	void printMsgAndExit(char* msg);
	void move_file_pointer(FILE *file_ptr, int lineNumber, int r);
	int linear3dIndex(int x, int y, int z, int max_x, int max_y) ;
	void writeToFile(Phi *phiFncn, char* fileName);
	void writeToNetCDF(Phi *phiFncn, char *fileName);
	
	// Handle errors by printing an error message and exiting with a
	// non-zero status.
	#define ERRCODE 2
	#define ERR(e) {printf("Error: %s\n", nc_strerror(e)); exit(ERRCODE);}

#endif /* UTILITIES_H_ */
