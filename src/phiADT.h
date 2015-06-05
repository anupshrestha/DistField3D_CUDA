/*
 * phiADT.h
 * Abstract Data Type for the Distance Field
 * Preprocessor.
 *
 *  Created on: Jun 04, 2015
 *      Author: Anup Shrestha
 *
 *  Audit Trail:
 *     Date: 
 *  Changes: 
 */
#ifndef PHIADT_H_
#define PHIADT_H_

#include <stdio.h>
#include <stdlib.h>
#include "file_writer.h"

// Type Definitions
typedef struct phi_type Phi;

// Method Declarations 
Phi *phi_create(FILE *);
void phi_destroy(Phi *);
void phi_calc_distance_field(Phi *);
void phi_gen_file(Phi *, FileOut, char *);

#endif /* PHIADT_H_ */
