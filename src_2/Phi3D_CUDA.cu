#include "Phi3D_CUDA.h"

static void update_distance(int *loc_d, double *dst_d, Grid3D *g3d);
static void fast_sweep(Phi3D* phiFncn, int noOfTimesToSweep, cudaPitchedPtr dst_dPitchPtr);
static void set_distance_negative_inside(Phi3D_d *d_phiFncn, Grid3D *g3d);
static void cudaMemcpy3D_P2D(Phi3D_d *d_phiFncn, double *dst, int *loc);
static void cudaMemcpy3D_D2P(Phi3D_d *d_phiFncn, double *dst, int *loc);

__global__ void update_distance_kernel(int *loc, double *dst, Grid3D g3d);
__global__ void set_distance_negative_inside_kernel(cudaPitchedPtr dst_dPitchPtr, cudaPitchedPtr loc_dPitchPtr, Grid3D g3d);
__global__ void fast_sweep_kernel(cudaPitchedPtr dst_dPitchPtr, SweepInfo s);
__device__ double solve_eikonal(double cur_dist, double minX, double minY, double minZ, double dx, double dy, double dz);


int iDivUp(int a, int b) {
	return ( (a % b) != 0 ) ? (a / b + 1) : (a / b);
}

void init_phiFncn(Phi3D *phiFncn, Grid3D *g3d, double *d) {
	
	g3d->x = (int) d[0] + 3;
	g3d->y = (int) d[1] + 3;
	g3d->z = (int) d[2] + 3;
	g3d->totalNodes = g3d->x * g3d->y * g3d->z;
	phiFncn->F = SPEED;
	phiFncn->x = (int) d[0] + 1; phiFncn->dx = d[3];
	phiFncn->y = (int) d[1] + 1; phiFncn->dy = d[4];
	phiFncn->z = (int) d[2] + 1; phiFncn->dz = d[5];

	// allocating pinned memory for the 
	// location and distance arrays
	cudaHostAlloc((void **)&phiFncn->location, sizeof(int) * g3d->totalNodes, cudaHostAllocMapped);
	cudaCheckError();
	cudaHostAlloc((void **)&phiFncn->distance, sizeof(double) * g3d->totalNodes, cudaHostAllocMapped);
	cudaCheckError();
}

void free_phiFncn(Phi3D *phiFncn) {
	cudaFreeHost(phiFncn->location); cudaCheckError();
	cudaFreeHost(phiFncn->distance); cudaCheckError();
	free(phiFncn);
}

void calculate_distfield(Phi3D *phiFncn, Grid3D *g3d) {
	
	// get the device pointers to the pinned memory
	int *loc_d; double *dst_d;
	cudaHostGetDevicePointer(&loc_d, phiFncn->location, 0); cudaCheckError();
	cudaHostGetDevicePointer(&dst_d, phiFncn->distance, 0); cudaCheckError();

	update_distance(loc_d, dst_d, g3d);
	
	// Setup for running the calculation on device (GPU)
	Phi3D_d* d_phiFncn = (Phi3D_d *) malloc (sizeof(Phi3D_d));
	
	// device memory allocation of location and
	// distance arrays using cudaMalloc3D
	d_phiFncn->loc_ext = make_cudaExtent(g3d->x * sizeof(int), g3d->y, g3d->z);
	d_phiFncn->dst_ext = make_cudaExtent(g3d->x * sizeof(double), g3d->y, g3d->z);

	cudaMalloc3D(&d_phiFncn->loc_dPitchPtr, d_phiFncn->loc_ext); cudaCheckError();
	cudaMalloc3D(&d_phiFncn->dst_dPitchPtr, d_phiFncn->dst_ext); cudaCheckError();

	// copying the distance and location array in pinned memory
	// to device memory allocated using cudaMalloc3D
	cudaMemcpy3D_P2D(d_phiFncn, dst_d, loc_d);

	// run parallel fast sweeping method to solve the
	// Eikonal equation and update distance
	fast_sweep(phiFncn, 3, d_phiFncn->dst_dPitchPtr);

	// set the inside distance to negative
	set_distance_negative_inside(d_phiFncn, g3d);

	// copying the distance and location array in device memory
	// to pinned memory
	cudaMemcpy3D_D2P(d_phiFncn, dst_d, loc_d);
	
	// freeing all device allocated memory
	cudaFree(d_phiFncn->loc_dPitchPtr.ptr);
	cudaFree(d_phiFncn->dst_dPitchPtr.ptr);
	free(d_phiFncn);
}

static void update_distance(int *loc_d, double *dst_d, Grid3D *g3d) {

	// Setup 3D-Grid and 3D-Block for Kernel launch
	// Running 256 threads per block
	dim3 bs(8, 8, 8);
	dim3 gs(iDivUp(g3d->x, bs.x), iDivUp(g3d->y , bs.y), iDivUp(g3d->z, bs.z));

	update_distance_kernel<<<gs, bs>>>(loc_d, dst_d, *g3d);
	cudaThreadSynchronize(); cudaCheckError();

}

/* 
 * Set up for the parallel FSM implementation
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
 * 	             Phi* phiFncn [in]: pointer to the Phi
 * 	     int noOfTimesToSweep [in]: number of sweep iterations
 * 	cudaPitchedPtr dstPtr [in/out]: pointer to distance array in
 *                                  device memory
 * 	Returns:
 */
static void fast_sweep(Phi3D* phiFncn, int noOfTimesToSweep, cudaPitchedPtr dst_dPitchPtr) {
	// Information regarding sweeping and linear indexing
	int meshDim = 3;
	SweepInfo sw;
	sw.xDim = phiFncn->x; sw.dx = phiFncn->dx;
	sw.yDim = phiFncn->y; sw.dy = phiFncn->dy;
	sw.zDim = phiFncn->z; sw.dz = phiFncn->dz;
	
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

				fast_sweep_kernel<<<gs, bs>>>(dst_dPitchPtr, sw);
				cudaThreadSynchronize();
				cudaCheckError();
			}
		}
		printf("Sweeping finished!......[%d/%d]\n", fastSweepLoopCount, noOfTimesToSweep);
		++fastSweepLoopCount;
	}
}

static void set_distance_negative_inside(Phi3D_d *d_phiFncn, Grid3D *g3d) {

	// Setup 3D-Grid and 3D-Block for Kernel launch
	// Running 256 threads per block
	dim3 bs(8, 8, 8);
	dim3 gs(iDivUp(g3d->x, bs.x), iDivUp(g3d->y , bs.y), iDivUp(g3d->z, bs.z));
	set_distance_negative_inside_kernel<<<gs, bs>>>(d_phiFncn->dst_dPitchPtr, d_phiFncn->loc_dPitchPtr, *g3d);
	cudaThreadSynchronize(); cudaCheckError();

}

__global__ void update_distance_kernel(int *loc, double *dst, Grid3D g3d) {
	int blockId = blockIdx.x + blockIdx.y * gridDim.x + gridDim.x * gridDim.y * blockIdx.z;
	int tid = blockId * (blockDim.x * blockDim.y * blockDim.z) + 
	                    (threadIdx.z * (blockDim.x * blockDim.y))+ 
	                    (threadIdx.y * blockDim.x) + threadIdx.x;

	//int tid = blockIdx.x * blockDim.x + threadIdx.x;

	if (tid < g3d.totalNodes)  {
		int l = loc[tid];
		double d = dst[tid];
		if(l != DEFAULT_BORDER_LOCATION &&
		   d != DEFAULT_BORDER_DISTANCE ) {
			dst[tid] = (l == 1 && d == INFINITY) ? -1 : (d > 0.0 || d < 0.0) ? d : DEFAULT_INTERIOR_DISTANCE;
		}
	}
}

__global__ void fast_sweep_kernel(cudaPitchedPtr dst_dPitchPtr, SweepInfo s) {
	int x = (blockIdx.x * blockDim.x + threadIdx.x) + s.xOffSet;
	int y = (blockIdx.y * blockDim.y + threadIdx.y) + s.yOffset;
	if(x <= s.xDim && y <= s.yDim) {
		int z = s.level - (x+y);
		if(z > 0 && z <= s.zDim){
			int i = abs(z-s.zSweepOff);
			int j = abs(y-s.ySweepOff);
			int k = abs(x-s.xSweepOff);

			char *devPtr = (char *) dst_dPitchPtr.ptr;
			size_t pitch = dst_dPitchPtr.pitch;
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
			c_row[k] = solve_eikonal(center, minX, minY, minZ, s.dx, s.dy, s.dz);
		}
	}
}

__device__ double solve_eikonal(double cur_dist, double minX, double minY, double minZ, double dx, double dy, double dz) {
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


__global__ void set_distance_negative_inside_kernel(cudaPitchedPtr dst_dPitchPtr, cudaPitchedPtr loc_dPitchPtr, Grid3D g3d) {
	int x = blockIdx.x*blockDim.x+threadIdx.x + 1;
	int y = blockIdx.y*blockDim.y+threadIdx.y + 1;
	int z = blockIdx.z*blockDim.z+threadIdx.z + 1;

	if (x < g3d.x-1 &&  y < g3d.y-1 && z < g3d.z-1 )  {

		char *dist_devPtr = (char *) dst_dPitchPtr.ptr;
		size_t dist_pitch = dst_dPitchPtr.pitch;
		size_t dist_slicePitch = dist_pitch * g3d.y;
		char* dist_slice = dist_devPtr + z * dist_slicePitch;
		double* dist_row = (double *) (dist_slice + y * dist_pitch);

		char *loc_devPtr = (char *) loc_dPitchPtr.ptr;
		size_t loc_pitch = loc_dPitchPtr.pitch;
		size_t loc_slicePitch = loc_pitch * g3d.y;
		char* loc_slice = loc_devPtr + z * loc_slicePitch;
		int* loc_row = (int *) (loc_slice + y * loc_pitch);

		if(loc_row[x] == 1) dist_row[x] = -1;
	}
}


/* 
 * Copies location and distance array memory
 * from device pinned memory to device memory
 * Argument(s):
 * 	           int * loc [in]: pointer to location (Pinned - Source)
 *          double * dst [in]: pointer to distance (Pinned - Source)
 * 	Phi_d* d_phiFncn [in/out]: pointer to d_Phi (Device - Destination)
 * Returns:
 */
static void cudaMemcpy3D_P2D(Phi3D_d *d_phiFncn, double *dst, int *loc) {
	cudaMemcpy3DParms mcpy3D_p = { 0 };
	mcpy3D_p.kind = cudaMemcpyDeviceToDevice;

	// copy parameters for distance
	mcpy3D_p.dstPtr = d_phiFncn->dst_dPitchPtr;
	mcpy3D_p.srcPtr.ptr = dst;
	mcpy3D_p.srcPtr.pitch = d_phiFncn->dst_ext.width;
	mcpy3D_p.srcPtr.xsize = (size_t) (d_phiFncn->dst_ext.width/sizeof(double));
	mcpy3D_p.srcPtr.ysize = d_phiFncn->dst_ext.height;
	mcpy3D_p.extent = d_phiFncn->dst_ext;
	cudaMemcpy3D(&mcpy3D_p); cudaCheckError();

	// copy parameters for location
	mcpy3D_p.dstPtr = d_phiFncn->loc_dPitchPtr;
	mcpy3D_p.srcPtr.ptr = loc;
	mcpy3D_p.srcPtr.pitch = d_phiFncn->loc_ext.width;
	mcpy3D_p.srcPtr.xsize = (size_t) (d_phiFncn->loc_ext.width/sizeof(int));
	mcpy3D_p.srcPtr.ysize = d_phiFncn->loc_ext.height;
	mcpy3D_p.extent = d_phiFncn->loc_ext;
	cudaMemcpy3D(&mcpy3D_p); cudaCheckError();
}

static void cudaMemcpy3D_D2P(Phi3D_d *d_phiFncn, double *dst, int *loc) {
	cudaMemcpy3DParms mcpy3D_p = { 0 };
	mcpy3D_p.kind = cudaMemcpyDeviceToDevice;

	// copy parameters for distance
	mcpy3D_p.srcPtr = d_phiFncn->dst_dPitchPtr;
	mcpy3D_p.dstPtr.ptr = dst;
	mcpy3D_p.dstPtr.pitch = d_phiFncn->dst_ext.width;
	mcpy3D_p.dstPtr.xsize = (size_t) (d_phiFncn->dst_ext.width/sizeof(double));
	mcpy3D_p.dstPtr.ysize = d_phiFncn->dst_ext.height;
	mcpy3D_p.extent = d_phiFncn->dst_ext;
	cudaMemcpy3D(&mcpy3D_p); cudaCheckError();

	// copy parameters for location
	mcpy3D_p.srcPtr = d_phiFncn->loc_dPitchPtr;
	mcpy3D_p.dstPtr.ptr = loc;
	mcpy3D_p.dstPtr.pitch = d_phiFncn->loc_ext.width;
	mcpy3D_p.dstPtr.xsize = (size_t) (d_phiFncn->loc_ext.width/sizeof(int));
	mcpy3D_p.dstPtr.ysize = d_phiFncn->loc_ext.height;
	mcpy3D_p.extent = d_phiFncn->loc_ext;
	cudaMemcpy3D(&mcpy3D_p); cudaCheckError();
}