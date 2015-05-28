# DistField3d_CUDA - GIN3D-DFP

Calculates in parallel using CUDA, the distance field for an interface
in a given VTI file using the Fast Sweeping Method.
The algorithm implemented for parallel fast sweeping
method is from a paper in the Journal of Computational Physics
titled "A parallel fast sweeping method for the Eikonal equation"
by Miles Detrixhe, Federic Gibou, and Chohong Min.
DOI: http://dx.doi.org/10.1016j.jcp.2012.112042

# Build GIN3D-DFP
Clean the project

    make clobber // delete directories and files
    make clean   // delete .o, .vti, .nc, .dat files
Build the project

    make
Execute the program

    cd bin
    ./GIN3D-DFP <filename> <prefix>
