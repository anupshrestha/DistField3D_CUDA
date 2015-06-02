#include "Utilities.h"
#include "VTI_Parser.h"
#include "Phi3D_FSM.h"
#include "Phi3D.h"

typedef struct  {
	cudaExtent loc_ext, dist_ext;
	cudaPitchedPtr loc_dPitchPtr, dist_dPitch_ptr;
} Phi_d;

static Phi3D_Dim p3d_sz;
static void fastSweep(Phi* phiFncn, int noOfTimesToSweep, cudaPitchedPtr dist_dPitch_ptr);
static void cudaMemcpy3D_D2H(Phi *phiFncn, Phi_d *d_phiFncn);
static void cudaMemcpy3D_D2D(Phi_d *d_phiFncn, int *loc, double *dst);

/* 
 * Creates the phi function for fast sweeping method
 * Initializes all variables for the phi function
 * and allocates memory required by the location and
 * distance arrays.
 * Argument(s):
 * 	double *d [in]: Pointer to a double with the
 * 	                dimension data
 * Returns:
 * 	Phi* phiFncn: pointer to Phi (phi function)
 */
Phi* Phi3D_Init(double *d) {

	// p3d_sz can be accessed from 
	// anywhere within this file
	p3d_sz.x = d[0]+3;
	p3d_sz.y = d[1]+3;
	p3d_sz.z = d[2]+3;

	Phi *phiFncn = ( Phi * ) malloc ( sizeof( Phi ) );

	// initialize the variables
	phiFncn->F = 1; // speed of propagation [Equation (1) in the report]
	phiFncn->size_x = (int) d[0] + 1; phiFncn->dx = d[3];
	phiFncn->size_y = (int) d[1] + 1; phiFncn->dy = d[4];
	phiFncn->size_z = (int) d[2] + 1; phiFncn->dz = d[5];

	int totalGridNodes = p3d_sz.x * p3d_sz.y * p3d_sz.z;

	// allocating host memory to store location
	// and distance data
	phiFncn->location = (int *) malloc( sizeof(int) * totalGridNodes );
	phiFncn->distance = (double *) malloc( sizeof(double) * totalGridNodes );

	return phiFncn;
}

/* Deallocates all memory allocated by
 * the create_phiFncn()
 * Argument(s):
 * 	Phi* phiFncn [in/out]: pointer to Phi (phi function)
 * 	Returns:
 */
void Phi3D_Finalize(Phi *phiFncn) {
	free(phiFncn->location);
	free(phiFncn->distance);
	free(phiFncn);
}

void Phi3D_Calc_distField(Phi *phiFncn, FILE *vti) {
	int totalGridNodes = p3d_sz.x * p3d_sz.y * p3d_sz.z;
	cudaSetDeviceFlags(cudaDeviceMapHost); cudaCheckError();
	
	// pointers to pinned memory
	int *loc_h, *loc_d;
	double *dst_h, *dst_d;
	
	// allocating pinned memory
	cudaHostAlloc((void **)&loc_h, sizeof(int) * totalGridNodes, cudaHostAllocMapped);
	cudaCheckError();
	cudaHostAlloc((void **)&dst_h, sizeof(double) * totalGridNodes, cudaHostAllocMapped);
	cudaCheckError();

	// parsing the location and distance data
	// from the file and storing them in the 
	// pinned memory array
	VTI_Get_locAndDist(vti, loc_h, dst_h, p3d_sz);

	// once we have everything for our distance
	// field calculation, we copy it to device (GPU)
	// where the calculations are done
	// for increased performance of the calculation
	// we will allocate memory in device using the 
	// cudaMalloc3D API, which aligns our 3D grid in
	// device memory using a pitch, for faster
	// access to memory
	Phi_d* d_phiFncn = ( Phi_d * ) malloc ( sizeof( Phi_d ) );
	d_phiFncn->loc_ext = make_cudaExtent(p3d_sz.x * sizeof(int), p3d_sz.y, p3d_sz.z);
	d_phiFncn->dist_ext = make_cudaExtent(p3d_sz.x * sizeof(double), p3d_sz.y, p3d_sz.z);

	// allocating device memory using cudaMalloc3D
	cudaMalloc3D(&d_phiFncn->loc_dPitchPtr, d_phiFncn->loc_ext); cudaCheckError();
	cudaMalloc3D(&d_phiFncn->dist_dPitch_ptr, d_phiFncn->dist_ext); cudaCheckError();

	// the device pointers to the pinned memory
	cudaHostGetDevicePointer(&loc_d, loc_h, 0); cudaCheckError();
	cudaHostGetDevicePointer(&dst_d, dst_h, 0); cudaCheckError();
	
	// copying the pinned memory to device memory
	// allocated using cudaMalloc3D
	cudaMemcpy3D_D2D(d_phiFncn, loc_d, dst_d);

	// free the pinned memory
	cudaFreeHost(loc_h); cudaCheckError();
	cudaFreeHost(dst_h); cudaCheckError();

	// CUDA Timer for FSM
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

	// free events
	cudaEventDestroy(start); cudaCheckError();
	cudaEventDestroy(stop);  cudaCheckError();


	// Setup 3D-Grid and 3D-Block for Kernel launch
  // Running 256 threads per block	
	dim3 bs(8, 8, 8);
	dim3 gs(iDivUp(p3d_sz.x, bs.x), iDivUp(p3d_sz.y , bs.y), iDivUp(p3d_sz.z, bs.z));
	setNegativeDistanceInside<<<gs, bs>>>(d_phiFncn->dist_dPitch_ptr,  d_phiFncn->loc_dPitchPtr, p3d_sz.x, p3d_sz.y, p3d_sz.z);
	cudaThreadSynchronize(); cudaCheckError();

	cudaMemcpy3D_D2H(phiFncn, d_phiFncn);

	cudaFree(d_phiFncn->loc_dPitchPtr.ptr);
	cudaFree(d_phiFncn->dist_dPitch_ptr.ptr);
	free(d_phiFncn);

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
static void fastSweep(Phi* phiFncn, int noOfTimesToSweep, cudaPitchedPtr dist_dPitch_ptr) {
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


/* 
 * Copies location and distance array memory
 * from device to host
 * Argument(s):
 * 	    Phi* phiFncn: pointer to Phi (Host - Destination)
 * 	Phi_d* d_phiFncn: pointer to d_Phi (Device - Source)
 * Returns:
 */
static void cudaMemcpy3D_D2H(Phi *phiFncn, Phi_d *d_phiFncn) {
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

/* 
 * Copies location and distance array memory
 * from device pinned memory to device memory
 * Argument(s):
 * 	           int * loc [in]: pointer to location (Pinned - Source)
 *          double * dst [in]: pointer to distance (Pinned - Source)
 * 	Phi_d* d_phiFncn [in/out]: pointer to d_Phi (Device - Destination)
 * Returns:
 */
static void cudaMemcpy3D_D2D(Phi_d *d_phiFncn, int *loc, double *dst) {
	cudaMemcpy3DParms mcpy3D_p = { 0 };
	mcpy3D_p.kind = cudaMemcpyDeviceToDevice;

	// copy parameters for location
	mcpy3D_p.dstPtr = d_phiFncn->loc_dPitchPtr;
	mcpy3D_p.srcPtr.ptr = loc;
	mcpy3D_p.srcPtr.pitch = d_phiFncn->loc_ext.width;
	mcpy3D_p.srcPtr.xsize = (size_t) (d_phiFncn->loc_ext.width/sizeof(int));
	mcpy3D_p.srcPtr.ysize = d_phiFncn->loc_ext.height;
	mcpy3D_p.extent = d_phiFncn->loc_ext;
	cudaMemcpy3D(&mcpy3D_p); cudaCheckError();

	// copy parameters for distance
	mcpy3D_p.dstPtr = d_phiFncn->dist_dPitch_ptr;
	mcpy3D_p.srcPtr.ptr = dst;
	mcpy3D_p.srcPtr.pitch = d_phiFncn->dist_ext.width;
	mcpy3D_p.srcPtr.xsize = (size_t) (d_phiFncn->dist_ext.width/sizeof(double));
	mcpy3D_p.srcPtr.ysize = d_phiFncn->dist_ext.height;
	mcpy3D_p.extent = d_phiFncn->dist_ext;
	cudaMemcpy3D(&mcpy3D_p); cudaCheckError();

}