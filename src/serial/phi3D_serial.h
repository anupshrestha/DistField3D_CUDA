#ifndef PHI3D_SERIAL_H_
#define PHI3D_SERIAL_H_
  
  #include "../phi3D.h"
  
  // External Linkage Method Declarations
  extern void init_phiFncn(Phi *p, double *d);
  extern void free_phiFncn(Phi *p);
  extern void calc_distfield(Phi *p, Grid3D g3d);
  
#endif /* PHI3D_SERIAL_H_ */
