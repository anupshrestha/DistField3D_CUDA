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

/* Sets up the size of the grid and device memory
 * allocations for location and distance arrays
 * using cudaMalloc3D. Then calls the appropriate
 * methods to run the fast sweeping method.
 * Argument(s):
 * 	  Phi* phiFncn: pointer to the Phi (phi function)
 * 	  double *dims: Pointer to a double with the
 * 	                dimension data
 * 	FILE *file_ptr: pointer to a file
 * 	Returns:
 */
void setupAndRun_phiFSM(Phi* phiFncn, double *dims, FILE *vti) {
  
	d_Phi* d_phiFncn = ( d_Phi * ) malloc ( sizeof( d_Phi ) );

	// setup the total grid size including boundaries
	int grid_x = phiFncn->size_x + 2;
	int grid_y = phiFncn->size_y + 2;
	int grid_z = phiFncn->size_z + 2;
	
	// setup for device memory allocation of
	// location and distance arrays
	d_phiFncn->loc_ext = make_cudaExtent(grid_x * sizeof(int), grid_y, grid_z);
	d_phiFncn->dist_ext = make_cudaExtent(grid_x * sizeof(double), grid_y, grid_z);

	// allocating device memory using cudaMalloc3D
	cudaMalloc3D(&d_phiFncn->loc_dPitchPtr, d_phiFncn->loc_ext); cudaCheckError();
	cudaMalloc3D(&d_phiFncn->dist_dPitch_ptr, d_phiFncn->dist_ext); cudaCheckError();

	initializeGrid(phiFncn, d_phiFncn);

	cudaMemcpy3D_D2H(phiFncn, d_phiFncn);
	update_phiFSM(vti, phiFncn);
	cudaMemcpy3D_H2D(phiFncn, d_phiFncn);

	// CUDA Timer
	cudaEvent_t start, stop;
	float elapsedTime;
	cudaEventCreate(&start);
	
	cudaEventCreate(&stop);
	cudaEventRecord(start, 0);

	fastSweep(phiFncn, 3, d_phiFncn->dist_dPitch_ptr);

	cudaEventRecord(stop, 0);
	cudaEventSynchronize(stop);
	cudaEventElapsedTime(&elapsedTime,start,stop);
	printf("Elapsed Time for FSM: %f s\n", elapsedTime/1000.0);

	// Setup 3D-Grid and 3D-Block for Kernel launch
        // Running 256 threads per block	
	dim3 bs(8, 8, 8);
	dim3 gs(iDivUp(grid_x, bs.x), iDivUp(grid_y , bs.y), iDivUp(grid_z, bs.z));
	setNegativeDistanceInside<<<gs, bs>>>(d_phiFncn->dist_dPitch_ptr,  d_phiFncn->loc_dPitchPtr, grid_x, grid_y, grid_z);
	cudaThreadSynchronize(); cudaCheckError();

	cudaMemcpy3D_D2H(phiFncn, d_phiFncn);

	cudaFree(d_phiFncn->loc_dPitchPtr.ptr);
	cudaFree(d_phiFncn->dist_dPitch_ptr.ptr);
	free(d_phiFncn);
}

/* Initializes the points in the grid by setting the
 * location to -1. Then calls the kernel to set the
 * distance of the interior and border points on the
 * grid.
 * Argument(s):
 * 	    Phi* phiFncn: pointer to the Phi (Host)
 * 	d_Phi* d_phiFncn: pointer to the d_Phi (Device)
 * Returns:
 */
void initializeGrid(Phi* phiFncn, d_Phi* d_phiFncn) {

	int grid_x = phiFncn->size_x + 2;
	int grid_y = phiFncn->size_y + 2;
	int grid_z = phiFncn->size_z + 2;

	// set every value in the location array in device to -1
	cudaMemset3D(d_phiFncn->loc_dPitchPtr, -1, d_phiFncn->loc_ext); cudaCheckError();

	// Setup 3D-Grid and 3D-Block for Kernel launch
	// Running 256 threads per block
	dim3 bs(8, 8, 8);
	dim3 gs(iDivUp(grid_x, bs.x), iDivUp(grid_y , bs.y), iDivUp(grid_z, bs.z));
	setInteriorAndBorderDistance<<<gs, bs>>>(d_phiFncn->dist_dPitch_ptr, grid_x, grid_y, grid_z);
	cudaThreadSynchronize(); cudaCheckError();
}

/* Set up for the parallel FSM implementation
 * Determines the number of levels and the nodes
 * on each level that are to be run in parallel.
 * Determines the direction of the sweep and
 * chooses and offset to translate the coordinate
 * for sweeping from different directions. Then
 * calls the CUDA kernel that calculates the
 * distance field for a 3D grid using Gauss-Seidal
 * iterations.
 * (1) i = 1:I, j = 1:J, k = 1:K
 * (2) i = I:1, j = 1:J, k = K:1
 * (3) i = I:1, j = 1:J, k = 1:K
 * (4) i = 1:I, j = 1:J, k = K:1
 * (5) i = I:1, j = J:1, k = K:1
 * (6) i = 1:I, j = J:1, k = 1:K
 * (7) i = 1:I, j = J:1, k = K:1
 * (8) i = I:1, j = J:1, k = 1:K
 *
 * This is Sweeping step discussed in section 2.3 of the report.
 * Argument(s):
 * 	                  Phi* phiFncn: pointer to the Phi
 * 	          int noOfTimesToSweep: number of sweep iterations
 * 	cudaPitchedPtr dist_dPitch_ptr: pointer to distance array in
 *                                  device memory
 * 	Returns:
 */
void fastSweep(Phi* phiFncn, int noOfTimesToSweep, cudaPitchedPtr dist_dPitch_ptr) {
	// Information regarding sweeping and linear indexing
	int meshDim = 3;
	SweepInfo sw;
	sw.xDim = phiFncn->size_x;
	sw.yDim = phiFncn->size_y;
	sw.zDim = phiFncn->size_z;
	sw.dx = phiFncn->dx;
	sw.dy = phiFncn->dy;
	sw.dz = phiFncn->dz;

	int totalLevels = sw.xDim + sw.yDim + sw.zDim;

	// loop till the number of times to sweep
	int fastSweepLoopCount = 1;
	while( fastSweepLoopCount <= noOfTimesToSweep){
		printf("Please wait. Sweeping...[%d/%d]\n", fastSweepLoopCount, noOfTimesToSweep);
		for(int swCount = 1; swCount <=8; ++swCount){
			int start = (swCount == 2 || swCount == 5 || swCount == 7 || swCount == 8 ) ? totalLevels : meshDim;
			int end = ( start == meshDim ) ? totalLevels + 1 : meshDim-1;
			int incr = ( start == meshDim ) ? true : false;

			// sweep offset is used for translating the 3D coordinates
			// to perform sweeps from different directions
			sw.xSweepOff = (swCount == 4 || swCount == 8 ) ? sw.xDim + 1 : 0;
			sw.ySweepOff = (swCount == 2 || swCount == 6 ) ? sw.yDim + 1 : 0;
			sw.zSweepOff = (swCount == 3 || swCount == 7 ) ? sw.zDim + 1 : 0;

			for(int level = start; level != end; level = (incr) ? level+1 : level-1){
				int xs = max(1, level-(sw.yDim + sw.zDim)), ys = max(1,level-(sw.xDim + sw.zDim));
				int xe = min(sw.xDim, level-(meshDim-1)),   ye = min(sw.yDim, level-(meshDim-1));

				int xr = xe-xs + 1, yr = ye-ys + 1;
				int tth = xr * yr; // Total number of threads needed

				dim3 bs(16, 16, 1);
				if(tth < 256){
					bs.x = xr;
					bs.y = yr;
				}
				dim3 gs(iDivUp(xr, bs.x), iDivUp(yr , bs.y), 1);

				sw.level = level;
				sw.xOffSet = xs;
				sw.yOffset = ys;

				parallel_fastSweep<<<gs, bs>>>(dist_dPitch_ptr, sw);
				cudaThreadSynchronize();
				cudaCheckError();
			}
		}
		printf("Sweeping finished!......[%d/%d]\n", fastSweepLoopCount, noOfTimesToSweep);
		++fastSweepLoopCount;
	}
}

/* Parses the location and distance data in the
 * input file and updates the data array.
 * Called before starting the sweep.
 * Argument(s):
 * 	  Phi* phiFncn: pointer to the Phi (phi function)
 * 	FILE *file_ptr: pointer to a file
 * 	Returns:
 */
void update_phiFSM(FILE* vti, Phi* phiFncn) {
	int intGrid_x = phiFncn->size_x + 1;
	int intGrid_y = phiFncn->size_y + 1;
	int intGrid_z = phiFncn->size_z + 1;

	move_file_pointer(vti, 6, 1);
	for (int i =1; i < intGrid_z; i++){
		for (int j =1; j < intGrid_y; j++) {
			for (int k =1; k < intGrid_x; k++) {
				int index = linear3dIndex(k, j, i, intGrid_x + 1, intGrid_y + 1);
				fscanf(vti, "%d ", &phiFncn->location[index]);
			}
		}
	}

	move_file_pointer(vti, 2, 0);

	for (int i =1; i < intGrid_z; i++){
		for (int j =1; j < intGrid_y; j++) {
			for (int k =1; k < intGrid_x; k++) {
				double temp;
				int index = linear3dIndex(k, j, i, intGrid_x + 1, intGrid_y + 1);
				fscanf(vti, "%lf ", &temp);
				if(phiFncn->location[index] == 1 && phiFncn->distance[index] == INFINITY) {
					phiFncn->distance[index] = -1;
				}
				else if ( (temp > 0.0) || ( temp < 0.0) ) {
					phiFncn->distance[index] = temp;
				}
			}
		}
	}
}

// CUDA  Kernels
__global__ void setInteriorAndBorderDistance(cudaPitchedPtr dist_dPitch_ptr, int w, int h, int d){
	int x = blockIdx.x*blockDim.x+threadIdx.x;
	int y = blockIdx.y*blockDim.y+threadIdx.y;
	int z = blockIdx.z*blockDim.z+threadIdx.z;

	if ( x < w && y < h && z < d) {
		char *devPtr = (char *) dist_dPitch_ptr.ptr;
		size_t pitch = dist_dPitch_ptr.pitch;
		size_t slicePitch = pitch * h;
		char* slice = devPtr + z * slicePitch;
		double* row = (double *) (slice + y * pitch);

		// Border distance
		if (x == 0 || x == w-1 || y == 0 || y == h-1 || z == 0 || z == d-1 )  {
			row[x] = INFINITY;
		}
		else{ // Interior distance
			row[x] = 90000;
		}
	}
}

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

/* Copies location and distance array memory
 * from device to host
 * Argument(s):
 * 	    Phi* phiFncn: pointer to the Phi (Host - Destination)
 * 	d_Phi* d_phiFncn: pointer to the d_Phi (Device - Source)
 * Returns:
 */
void cudaMemcpy3D_D2H(Phi *phiFncn, d_Phi *d_phiFncn) {
	cudaMemcpy3DParms mcpy3D_p = { 0 };
	mcpy3D_p.kind = cudaMemcpyDeviceToHost;

	// copy parameters for location
	mcpy3D_p.srcPtr = d_phiFncn->loc_dPitchPtr;
	mcpy3D_p.dstPtr.ptr = phiFncn->location;
	mcpy3D_p.dstPtr.pitch = d_phiFncn->loc_ext.width;
	mcpy3D_p.dstPtr.xsize = (size_t) (d_phiFncn->loc_ext.width/sizeof(int));
	mcpy3D_p.dstPtr.ysize = d_phiFncn->loc_ext.height;
	mcpy3D_p.extent = d_phiFncn->loc_ext;
	cudaMemcpy3D(&mcpy3D_p); cudaCheckError();

	// copy parameters for distance
	mcpy3D_p.srcPtr = d_phiFncn->dist_dPitch_ptr;
	mcpy3D_p.dstPtr.ptr = phiFncn->distance;
	mcpy3D_p.dstPtr.pitch = d_phiFncn->dist_ext.width;
	mcpy3D_p.dstPtr.xsize = (size_t) (d_phiFncn->dist_ext.width/sizeof(double));
	mcpy3D_p.dstPtr.ysize = d_phiFncn->dist_ext.height;
	mcpy3D_p.extent = d_phiFncn->dist_ext;
	cudaMemcpy3D(&mcpy3D_p); cudaCheckError();
}

/* Copies location and distance array memory
 * from host to device
 * Argument(s):
 * 	    Phi* phiFncn: pointer to the Phi (Host - Source)
 * 	d_Phi* d_phiFncn: pointer to the d_Phi (Device - Destination)
 * Returns:
 */
void cudaMemcpy3D_H2D(Phi *phiFncn, d_Phi *d_phiFncn) {
	cudaMemcpy3DParms mcpy3D_p = { 0 };
	mcpy3D_p.kind = cudaMemcpyHostToDevice;

	// copy parameters for location
	mcpy3D_p.dstPtr = d_phiFncn->loc_dPitchPtr;
	mcpy3D_p.srcPtr.ptr = phiFncn->location;
	mcpy3D_p.srcPtr.pitch = d_phiFncn->loc_ext.width;
	mcpy3D_p.srcPtr.xsize = (size_t) (d_phiFncn->loc_ext.width/sizeof(int));
	mcpy3D_p.srcPtr.ysize = d_phiFncn->loc_ext.height;
	mcpy3D_p.extent = d_phiFncn->loc_ext;
	cudaMemcpy3D(&mcpy3D_p); cudaCheckError();

	// copy parameters for distance
	mcpy3D_p.dstPtr = d_phiFncn->dist_dPitch_ptr;
	mcpy3D_p.srcPtr.ptr = phiFncn->distance;
	mcpy3D_p.srcPtr.pitch = d_phiFncn->dist_ext.width;
	mcpy3D_p.srcPtr.xsize = (size_t) (d_phiFncn->dist_ext.width/sizeof(double));
	mcpy3D_p.srcPtr.ysize = d_phiFncn->dist_ext.height;
	mcpy3D_p.extent = d_phiFncn->dist_ext;
	cudaMemcpy3D(&mcpy3D_p); cudaCheckError();
}
