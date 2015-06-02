/*
 * Phi3D_Writer.h
 *
 *  Created on: Jun 02, 2015
 *      Author: Anup Shrestha
 */

#ifndef Phi3D_Writer_H_
#define Phi3D_Writer_H_
  
 #include "Phi3D.h"

 typedef enum { VTI = 0x01, DAT = 0x02, NCF = 0x04 } FileType;

 void write_to_file(Phi3D *phiFncn, Grid3D *g3d, FileType type, char *prefix);

#endif /* Phi3D_Writer_H_ */
