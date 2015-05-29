# DistField3D_CUDA - GIN3D-DFP

Calculates in parallel using CUDA, the distance field for an interface
in a given VTI file using the Fast Sweeping Method.
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

    make
Execute the program

    cd bin
    ./GIN3D-DFP <filename> <prefix>
