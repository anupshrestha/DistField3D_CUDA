/*
 * Phi3D.h
 * Main header file.
 *
 *  Created on: May 14, 2015
 *      Author: Anup Shrestha
 */
#include <stdio.h>
#include <math.h>

#ifndef PHI3D_H_
#define PHI3D_H_

 	#define SPEED 1
  #define DEFAULT_BORDER_LOCATION -1
  #define DEFAULT_BORDER_DISTANCE INFINITY

  // holds the total size of the
  // 3D grid in each dimension
 	typedef struct {
		int x, y, z;
	} Phi3D_Dim;

	// structure of the Phi Function
	// dx, dy, dz are used for propagation.
	typedef struct  {
		int    * location;
		double * distance;
		int      size_x, size_y, size_z;
		double   dx, dy, dz, F;
	} Phi;

	// Method Declarations
	Phi* Phi3D_Init(double *d);
	void Phi3D_Finalize(Phi *phiFncn);
	void Phi3D_Calc_distField(Phi *phiFncn, FILE *vti);

#endif /* PHI3D_H_ */
