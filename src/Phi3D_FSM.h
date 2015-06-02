/*
 * Phi3D_FSM.h
 *
 *  Created on: May 14, 2015
 *      Author: Anup Shrestha
 */

#include "Phi3D.h"

#ifndef PHI3D_FSM_H_
#define PHI3D_FSM_H_

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

	// Device Kernel/Method Declarations
	__global__ void setNegativeDistanceInside(cudaPitchedPtr dist_dPitch_ptr, cudaPitchedPtr loc_dPitchPtr, int w, int h, int d);
	__global__ void parallel_fastSweep(cudaPitchedPtr dist_dPitch_ptr, SweepInfo s);
	__device__ double solveEikonal(double cur_dist, double minX, double minY, double minZ, double dx, double dy, double dz);

#endif /* PHI3D_FSM_H_ */
