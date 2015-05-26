/*
 * Phi3D.h
 * Main header file.
 *
 *  Created on: May 14, 2015
 *      Author: Anup Shrestha
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifndef PHI3D_H_
#define PHI3D_H_

	// dx, dy, dz are used for propagation.
	typedef struct  {
		int *  location;
		double * distance;
		int size_x, size_y, size_z;
		double dx, dy, dz, F;
	} Phi;

	// Host Method Definitions
	// Methods implemented in DistField3D.c
	Phi* create_phiFncn(double *dims);
	void getDimensions(FILE *vti, double *dims);
	void destroy_phiFncn(Phi *phiFncn);
	void parse_cmdLineArgs(int argc, char *argv[]);

#endif /* PHI3D_H_ */
