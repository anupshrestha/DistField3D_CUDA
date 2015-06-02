/*
 * VTI_Parser.h
 *
 *  Created on: Jun 01, 2015
 *      Author: Anup Shrestha
 */

#ifndef VTI_Parser_H_
#define VTI_Parser_H_
  
  #include "Phi3D.h"
  
  void VTI_Get_dimensions(FILE *vti, double *d);
  void VTI_Get_locAndDist(FILE *vti, int *l, double *d, Grid3D *sz);

#endif /* VTI_Parser_H_ */
