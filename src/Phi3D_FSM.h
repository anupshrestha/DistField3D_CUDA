/*
 * Phi3D_FSM.h
 *
 *  Created on: May 14, 2015
 *      Author: Anup Shrestha
 */

#include "Phi3D.h"

#ifndef PHI3D_FSM_H_
#define PHI3D_FSM_H_

	typedef struct  {
		cudaExtent loc_ext, dist_ext;
		cudaPitchedPtr loc_dPitchPtr, dist_dPitch_ptr;
	} d_Phi;

	/* Structure for storing information used during
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

	extern "C" void setupAndRun_phiFSM(Phi* phiFncn, double *dims, FILE *vti);
	void initializeGrid(Phi* phiFncn, d_Phi* d_phiFncn);
	void update_phiFSM(FILE* vti, Phi* phiFncn);
	void fastSweep(Phi* phiFncn, int noOfTimesToSweep, cudaPitchedPtr dist_dPitch_ptr);
	void cudaMemcpy3D_D2H(Phi *phiFncn, d_Phi *d_phiFncn);
	void cudaMemcpy3D_H2D(Phi *phiFncn, d_Phi *d_phiFncn);

	// Device Kernel/Method Definitions
	__global__ void setInteriorAndBorderDistance(cudaPitchedPtr dist_dPitch_ptr, int w, int h, int d);
	__global__ void setNegativeDistanceInside(cudaPitchedPtr dist_dPitch_ptr, cudaPitchedPtr loc_dPitchPtr, int w, int h, int d);
	__global__ void parallel_fastSweep(cudaPitchedPtr dist_dPitch_ptr, SweepInfo s);
	__device__ double solveEikonal(double cur_dist, double minX, double minY, double minZ, double dx, double dy, double dz);

#endif /* PHI3D_FSM_H_ */
