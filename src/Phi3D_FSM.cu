/*
 * Phi3D_FSM.cu
 * Parallel 3D implementation of Fast Sweeping Method
 * using CUDA C/C++.
 * Takes the information about locations and IB distances
 * from a VTI file and propagates it.
 * The algorithm implemented for parallel fast sweeping
 * method is from a paper in the Journal of Computational
 * Physics titled "A parallel fast sweeping method for
 * the Eikonal equation" by Miles Detrixhe, Federic Gibou,
 * and Chohong Min.
 * DOI: http://www.sciencedirect.com/science/article/pii/S002199911200722X
 *
 *  Created on: May 14, 2015
 *      Author: Anup Shrestha
 *
 *  Audit Trail:
 *          Date:
 *  Changes:
 *
 */
#include "Phi3D.h"
#include "Phi3D_FSM.h"
extern "C"{
  #include "Utilities.h"
}

// CUDA  Kernels

__global__ void setNegativeDistanceInside(cudaPitchedPtr dist_dPitch_ptr, cudaPitchedPtr loc_dPitchPtr, int w, int h, int d) {
	int x = blockIdx.x*blockDim.x+threadIdx.x + 1;
	int y = blockIdx.y*blockDim.y+threadIdx.y + 1;
	int z = blockIdx.z*blockDim.z+threadIdx.z + 1;

	if (x < w-1 &&  y < h-1 && z < d-1 )  {

		char *dist_devPtr = (char *) dist_dPitch_ptr.ptr;
		size_t dist_pitch = dist_dPitch_ptr.pitch;
		size_t dist_slicePitch = dist_pitch * h;
		char* dist_slice = dist_devPtr + z * dist_slicePitch;
		double* dist_row = (double *) (dist_slice + y * dist_pitch);

		char *loc_devPtr = (char *) loc_dPitchPtr.ptr;
		size_t loc_pitch = loc_dPitchPtr.pitch;
		size_t loc_slicePitch = loc_pitch * h;
		char* loc_slice = loc_devPtr + z * loc_slicePitch;
		int* loc_row = (int *) (loc_slice + y * loc_pitch);

		if(loc_row[x] == 1) dist_row[x] = -1;
	}
}

__global__ void parallel_fastSweep(cudaPitchedPtr dist_dPitch_ptr, SweepInfo s) {
	int x = (blockIdx.x * blockDim.x + threadIdx.x) + s.xOffSet;
	int y = (blockIdx.y * blockDim.y + threadIdx.y) + s.yOffset;
	if(x <= s.xDim && y <= s.yDim) {
		int z = s.level - (x+y);
		if(z > 0 && z <= s.zDim){
			int i = abs(z-s.zSweepOff);
			int j = abs(y-s.ySweepOff);
			int k = abs(x-s.xSweepOff);

			char *devPtr = (char *) dist_dPitch_ptr.ptr;
			size_t pitch = dist_dPitch_ptr.pitch;
			size_t slicePitch = pitch * (s.yDim + 2);

			double *c_row = (double *) ( (devPtr + i * slicePitch) + j * pitch);          // center row
			double center = c_row[k];                                                     // center distance
			double left = c_row[k-1];                                                     // left distance
			double right = c_row[k+1];                                                    // right distance
			double up = ((double *) ( (devPtr + i * slicePitch) + (j-1) * pitch) )[k];    // upper distance
			double down = ((double *) ( (devPtr + i * slicePitch) + (j+1) * pitch))[k];   // lower distance
			double front = ((double *) ( (devPtr + (i-1) * slicePitch) + j * pitch) )[k]; // front distance
			double back = ((double *) ( (devPtr + (i+1) * slicePitch) + j * pitch) )[k];  // back distance

			double minX = min(left, right);
			double minY = min(up, down);
			double minZ = min(front, back);
			c_row[k] = solveEikonal(center, minX, minY, minZ, s.dx, s.dy, s.dz);
		}
	}
}

__device__ double solveEikonal(double cur_dist, double minX, double minY, double minZ, double dx, double dy, double dz) {
	double dist_new = 0;
	double m[] = { minX, minY, minZ} ;
	double d[] = { dx, dy, dz};
	// sort the mins
	for(int i = 1; i < 3; i++){
		for(int j = 0; j < 3-i; j++) {
			if(m[j] > m[j+1]) {
				double tmp_m = m[j];
				double tmp_d = d[j];
				m[j] = m[j+1]; d[j] = d[j+1];
				m[j+1] = tmp_m; d[j+1] = tmp_d;
			}
		}
	}

	double m_0 = m[0], m_1 = m[1], m_2 = m[2];
	double d_0 = d[0], d_1 = d[1], d_2 = d[2];
	double m2_0 = m_0 * m_0, m2_1 = m_1 * m_1, m2_2 = m_2 * m_2;
	double d2_0 = d_0 * d_0, d2_1 = d_1 * d_1, d2_2 = d_2 * d_2;

	dist_new = m_0 + d_0;
	if(dist_new > m_1) {
	  dist_new = ( m_1 * d2_0 + m_0 * d2_1 + d_0 * d_1 * sqrt(- m2_0 + 2 * m_0 * m_1 - m2_1 + d2_0 + d2_1) ) / (d2_0 + d2_1);
	  if(dist_new > m_2) {
	    double a = sqrt(- m2_0 * d2_1 - m2_0 * d2_2 + 2 * m_0 * m_1 * d2_2
	                    - m2_1 * d2_0 - m2_1 * d2_2 + 2 * m_0 * m_2 * d2_1
	                    - m2_2 * d2_0 - m2_2 * d2_1 + 2 * m_1 * m_2 * d2_0
	                    + d2_0 * d2_1 + d2_0 * d2_2 + d2_1 * d2_2);
	    dist_new = (m_2 * d2_0 * d2_1 + m_1 * d2_0 * d2_2 + m_0 * d2_1 * d2_2 + d_0 * d_1 * d_2 * a) /
	               (d2_0 * d2_1 + d2_0 * d2_2 + d2_1 * d2_2);
	  }
	}
	return min(cur_dist, dist_new);
}
