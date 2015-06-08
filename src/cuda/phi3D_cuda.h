/*
 * phi3D_cuda.h
 * Phi CUDA header file.
 *
 *  Created on: Jun 04, 2015
 *      Author: Anup Shrestha
 *
 *  Audit Trail:
 *     Date: 
 *  Changes: 
 */

#ifndef PHI3D_CUDA_H_
#define PHI3D_CUDA_H_
  
  #include "../phi3D.h"
  
  /*
   * Structure for storing the cudaExtent and the
   * cudaPitchedPtr for the 3D array memory in 
   * device (GPU) for location and distance.
   */
  typedef struct  {
    cudaExtent loc_ext, dst_ext;
    cudaPitchedPtr loc_dPitchPtr, dst_dPitchPtr;
  } Phi3D_d;

  /* 
   * Structure for storing information used during
   * sweeping to manage internal grid dimensions,
   * sweep directions, position of the node on the
   * array and its offset in each kernel block
   */
  typedef struct {
    int level;
    int xDim, yDim, zDim;
    int xOffSet, yOffset;
    int xSweepOff, ySweepOff, zSweepOff;
    double dx, dy, dz;
  } SweepInfo;

  // Macro for checking CUDA errors following a CUDA launch or API call
  #define cudaCheckError() {\
	  cudaError_t e = cudaGetLastError();\
	  if( e != cudaSuccess ) {\
		  printf("\nCuda failure %s:%d: '%s'\n",__FILE__,__LINE__,cudaGetErrorString(e));\
		  exit(EXIT_FAILURE);\
	  }\
  }

  // External Linkage Method Declarations
  extern "C" {
    void init_phiFncn(Phi *p, double *d);
    void free_phiFncn(Phi *p);
    void calc_distfield(Phi *p, Grid3D g3d);
  }
  
#endif /* PHI3D_CUDA_H_ */
