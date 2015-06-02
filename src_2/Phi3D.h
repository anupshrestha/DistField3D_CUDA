/*
* Phi3D.h
* Main header file.
*
*  Created on: May 14, 2015
*      Author: Anup Shrestha
*/

#include <stdlib.h>
#include <stdio.h>
#include "Phi3D_Writer.h"

#ifndef PHI3D_H_
#define PHI3D_H_
  
  /* speed of propagation [Equation (1) in the report] */
  #define SPEED 1
  #define DEFAULT_BORDER_LOCATION -1
  #define DEFAULT_BORDER_DISTANCE INFINITY
  #define DEFAULT_INTERIOR_DISTANCE 90000

  /*
   * Phi3D structure stores the interior dimensions
   * of the mesh(grid) as x, y, z and the spacing
   * of the nodes as dx, dy, dz which are used for
   * propagation. It also stores a pointer to an
   * int array for location values and a pointer to
   * a double array for distance values.
   */
  typedef struct  {
    int    * location;
    double * distance;
    int      x, y, z;
    double   dx, dy, dz, F;
  } Phi3D;

  /*
   * Grid3D structure stores the dimensions of 
   * the mesh(grid) and the total number of nodes
   */
  typedef struct {
    int x, y, z;
    int totalNodes;
  } Grid3D;
  
  extern Phi3D* Phi3D_Init(FILE *vti, Grid3D *g3d);
  extern void Phi3D_Calc_distfield(Phi3D *phiFncn, Grid3D *g3d);
  //extern void Phi3D_Write(Phi3D *phiFncn, Grid3D *g3d, FileType type, char *prefix);
  extern void Phi3D_Finalize(Phi3D *phiFncn);

  
#endif /* PHI3D_H_ */
