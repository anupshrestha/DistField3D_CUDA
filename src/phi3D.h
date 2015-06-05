/*
 * phi3D.h
 * Main header file for 3D Phi Function.
 *
 *  Created on: Jun 04, 2015
 *      Author: Anup Shrestha
 *
 *  Audit Trail:
 *     Date: 
 *  Changes: 
 */

#ifndef PHI3D_H_
#define PHI3D_H_

#include "phiADT.h"
#include "vti_parser.h"
#include <math.h>

/* speed of propagation [Equation (1) in the report] */
#define SPEED                     1
#define DEFAULT_BORDER_LOCATION   -1
#define DEFAULT_BORDER_DISTANCE   INFINITY
#define DEFAULT_INTERIOR_DISTANCE 90000

/*
 * Completing the struct declared in the 
 * phiADT.h This definition is for a 3D phi.
 *  x, y, z    - interior dimensions
 *  dx, dy, dz - node spacing
 *  F          - speed
 *  location   - location values
 *  distance   - distance values
 */
struct phi_type {
  int      x, y, z;
  double   dx, dy, dz, F;
  int    * location;
  double * distance;
};

#endif /* PHI3D_H_ */
