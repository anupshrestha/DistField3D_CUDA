# DistField3D_CUDA - GIN3D-DFP

The distance field propagation calculation for an interface
in a given VTI file using the Fast Sweeping Method. The Fast
Sweeping Method can be calculated using Serial or Parallel (CUDA)
implementation based on specifying the target while running make.
The algorithm implemented for parallel fast sweeping
method is from a paper in the Journal of Computational Physics
titled "A parallel fast sweeping method for the Eikonal equation"
by Miles Detrixhe, Federic Gibou, and Chohong Min.
DOI: http://www.sciencedirect.com/science/article/pii/S002199911200722X

# Build GIN3D-DFP
Clean the project

    make clobber // delete directories and files
    make clean   // delete .o, .vti, .nc, .dat files
Build the project

    // serial (CPU) implementation
    make serial
    // parallel (CUDA) implementation
    make parallel
Execute the program

    cd bin
    ./GIN3D-DFP <filename> <prefix>
