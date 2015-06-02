#include "Phi3D.h"
#include "VTI_Parser.h"
#include "Phi3D_Writer.h"

static int linear3dIndex(int x, int y, int z, int max_x, int max_y) {
  return z * max_y * max_x + y * max_x + x;
}

Phi3D* Phi3D_Init(FILE *vti, Grid3D *g3d) {
	
	Phi3D *phiFncn = (Phi3D *) malloc(sizeof(Phi3D));

	// get the dimensions from the vti file
	double *dims = (double *) malloc(sizeof(double) * 6);
	VTI_Get_dimensions(vti, dims);

	// initialize variables and allocate memory
	// for the phi function and the 3D grid
	init_phiFncn(phiFncn, g3d, dims);

	// parse and fetch the location and distance
	// data from the vti file and store in arrays
	// in phi function
	VTI_Get_locAndDist(vti, phiFncn->location, phiFncn->distance, g3d);

	return phiFncn;
}

void Phi3D_Calc_distfield(Phi3D *phiFncn, Grid3D *g3d) {
	calculate_distfield(phiFncn, g3d);
}

void Phi3D_Write(Phi3D *phiFncn, Grid3D *g3d, FileType type, char *prefix) {
  
  // before writing to file we need to update the 
  // boundary values
  // dimensions                
  int x = g3d->x;
  int y = g3d->y;
  int z = g3d->z;
  int i,j,k;
  for(i = 0; i < z; i++){
    for(j = 0; j < y; j++){
      for(k = 0; k < x; k++){
      	int I = i, J = j, K = k;
      	I = (i == z-1) ? I-1 : (!i) ? I+1 : I;
      	J = (j == y-1) ? J-1 : (!j) ? J+1 : J;
      	K = (k == x-1) ? K-1 : (!k) ? K+1 : K;
      	if( i != I || j != J || k != K)
      		phiFncn->distance[linear3dIndex(k, j, i, x, y)] = phiFncn->distance[linear3dIndex(K, J, I, x, y)];
        }
      }
    }
  }
  
  // write to file
  write_to_file(phiFncn, g3d, type, prefix);
}

void Phi3D_Finalize(Phi3D *phiFncn) {
  free_phiFncn(phiFncn);
}