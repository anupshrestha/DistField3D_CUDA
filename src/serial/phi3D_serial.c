#include "phi3D_serial.h"
#define min(a,b)(a>b)?b:a

static void fast_sweep(Phi *p, int noOfTimesToSweep);
static double solveEikonal(Phi *p, int i, int j, int k);
static void update_distance(int *loc, double *dst, Grid3D g3d);
static void set_distance_negative_inside(Phi *p, Grid3D g3d);
static int linear3dIndex(int x, int y, int z, int max_x, int max_y);

void init_phiFncn(Phi *p, double *d) {
	p->x = (int) d[0] + 1; p->dx = d[3];
	p->y = (int) d[1] + 1; p->dy = d[4];
	p->z = (int) d[2] + 1; p->dz = d[5];
	p->F = SPEED;

	// allocating memory for the 
	// location and distance arrays
	int totalNodes = (p->x + 2) * (p->y + 2) * (p->z + 2);
	size_t l_arr = sizeof(int) * totalNodes;
	size_t d_arr = sizeof(double) * totalNodes;
	p->location = (int *) malloc(l_arr);
	p->distance = (double *) malloc(d_arr);

}

void free_phiFncn(Phi *p) {
	free(p->location);
	free(p->distance);
	free(p);
}

void calc_distfield(Phi *p, Grid3D g3d) {
	update_distance(p->location, p->distance, g3d);
	fast_sweep(p, 3);
	set_distance_negative_inside(p, g3d);
}

static void update_distance(int *loc, double *dst, Grid3D g3d) {
	int i, j, k;
	int *lptr = &loc[0]; double *dptr = &dst[0];

	for (i = 0; i < g3d.z; i++){
		for (j = 0; j < g3d.y; j++) {
			for (k = 0; k < g3d.x; k++) {
				int l = *lptr; double d = *dptr;
				if(l != DEFAULT_BORDER_LOCATION &&
		       d != DEFAULT_BORDER_DISTANCE ) {
					*dptr = (l == 1 && d == INFINITY) ? -1 : (d > 0.0 || d < 0.0) ? d : DEFAULT_INTERIOR_DISTANCE;
				}
				lptr++; dptr++;
			}
		}
	}
}

static void fast_sweep(Phi *p, int noOfTimesToSweep) {	
  int s,i,j,k;

  int sweeps[8][3] = {{1,1,1}, {0,1,0}, {0,1,1}, {1,1,0}, {0,0,0}, {1,0,1}, {1,0,0}, {0,0,1}};

  // loop till the number of times to sweep
	int fastSweepLoopCount = 1;
	while( fastSweepLoopCount <= noOfTimesToSweep) {
	 printf("Please wait. Sweeping...[%d/%d]\n", fastSweepLoopCount, noOfTimesToSweep);
	  for(s = 0; s < 8; ++s) {
		  int iStart = (sweeps[s][0]) ? 1 : p->z;
		  int iEnd   = (sweeps[s][0]) ? p->z + 1 : 0;

		  int jStart = (sweeps[s][1]) ? 1 : p->y;
		  int jEnd   = (sweeps[s][1]) ? p->y + 1 : 0;

		  int kStart = (sweeps[s][2]) ? 1 : p->x;
		  int kEnd   = (sweeps[s][2]) ? p->x + 1 : 0;

		  for (i = iStart; i != iEnd; i = (sweeps[s][0]) ? i+1 : i-1) {
			  for (j = jStart; j != jEnd; j = (sweeps[s][1]) ? j+1 : j-1) {
				  for (k = kStart; k != kEnd; k = (sweeps[s][2]) ? k+1 : k-1) {
				  	int index = linear3dIndex(k, j, i, p->x+2, p->y+2);
					  p->distance[index] = solveEikonal(p, i, j, k);
				  }
			  }
		  }
	  }
	printf("Sweeping finished!......[%d/%d]\n", fastSweepLoopCount++, noOfTimesToSweep);
	}
}

static double solveEikonal(Phi *p, int i, int j, int k) {
	int max_x = p->x+2, max_y = p->y+2;

	double dist_new = 0;
	double dist_old = p->distance[linear3dIndex(k, j, i, max_x, max_y)];

	double dx = p->dx, dy = p->dy, dz = p->dz;
	double minX = min(p->distance[linear3dIndex(k-1,j,i,max_x,max_y)], p->distance[linear3dIndex(k+1,j,i,max_x,max_y)]);
	double minY = min(p->distance[linear3dIndex(k,j-1,i,max_x,max_y)], p->distance[linear3dIndex(k,j+1,i,max_x,max_y)]);
	double minZ = min(p->distance[linear3dIndex(k,j,i-1,max_x,max_y)], p->distance[linear3dIndex(k,j,i+1,max_x,max_y)]);

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

	// simplifying the variables
	double m_0 = m[0], m_1 = m[1], m_2 = m[2];
	double d_0 = d[0], d_1 = d[1], d_2 = d[2];
	double m2_0 = m_0 * m_0, m2_1 = m_1 * m_1, m2_2 = m_2 * m_2;
	double d2_0 = d_0 * d_0, d2_1 = d_1 * d_1, d2_2 = d_2 * d_2;

	dist_new = m_0 + d_0;
	if(dist_new > m_1) {
	  
	  double s = sqrt(- m2_0 + 2 * m_0 * m_1 - m2_1 + d2_0 + d2_1);
	  dist_new = ( m_1 * d2_0 + m_0 * d2_1 + d_0 * d_1 * s) / (d2_0 + d2_1);
	  
	  if(dist_new > m_2) {
	    
	    double a = sqrt(- m2_0 * d2_1 - m2_0 * d2_2 + 2 * m_0 * m_1 * d2_2
	                    - m2_1 * d2_0 - m2_1 * d2_2 + 2 * m_0 * m_2 * d2_1
	                    - m2_2 * d2_0 - m2_2 * d2_1 + 2 * m_1 * m_2 * d2_0
	                    + d2_0 * d2_1 + d2_0 * d2_2 + d2_1 * d2_2);

	    dist_new = (m_2 * d2_0 * d2_1 + m_1 * d2_0 * d2_2 + m_0 * d2_1 * d2_2 + d_0 * d_1 * d_2 * a) /
	               (d2_0 * d2_1 + d2_0 * d2_2 + d2_1 * d2_2);
	  }
	}

	return min(dist_old, dist_new);
}

static void set_distance_negative_inside(Phi *p, Grid3D g3d) {
	int i, j, k, index;
	for (i = 1; i < g3d.z-1; i++){
		for (j = 1; j < g3d.y-1; j++) {
			for (k = 1; k < g3d.x-1; k++) {	
				index = linear3dIndex(k, j, i, g3d.x, g3d.y);
				if(p->location[index] == 1) p->distance[index] = -1;
			}
		}
	}
}

/*
 * Convert 3D indexing to 1D indexing 
 * Arguments:
 *   int x, y, z [in] - 3D coordinate
 *   int max_x   [in] - size of x-dimension
 *   int max_y   [in] - size of y-dimension
 * Returns:
 */
static int linear3dIndex(int x, int y, int z, int max_x, int max_y) {
  return z * max_y * max_x + y * max_x + x;
}