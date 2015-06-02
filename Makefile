EXEC := GIN3D-DFP
SRCDIR := src
BINDIR := bin
OBJDIR := obj

OBJ := $(patsubst $(SRCDIR)/%.c,$(OBJDIR)/%.o,$(wildcard $(SRCDIR)/*.c))
CUOBJ := $(patsubst $(SRCDIR)/%.cu,$(OBJDIR)/%.o,$(wildcard $(SRCDIR)/*.cu))
TARGET :=  $(BINDIR)/$(EXEC)

# Compiler options
CC = gcc -Wall -g -std=c99
CCFLAGS = -I/usr/local/wrfUtils/netcdf-c/include -I/usr/local/wrfUtils/szip/include -I/usr/local/wrfUtils/hdf5/include
LDFLAGS = -L/usr/local/wrfUtils/netcdf-c/lib -L/usr/local/wrfUtils/hdf5/lib -L/usr/local/wrfUtils/szip/lib -lm -lhdf5 -lsz -lnetcdf

NVCC = nvcc 
NVCCFLAGS = -arch=sm_20

$(OBJDIR)/%.o: $(SRCDIR)/%.c
	$(CC) $(CCFLAGS) -o $@ -c $< $(LDFLAGS)
	
$(OBJDIR)/%.o: $(SRCDIR)/%.cu
	$(NVCC) $(NVCCFLAGS) $(CCFLAGS) $(LDFLAGS) -o $@ -c $< 

$(TARGET): makedirectories $(OBJ) $(CUOBJ)
	$(NVCC) $(NVCCFLAGS) $(CCFLAGS) $(LDFLAGS) -o $@ $(OBJ) $(CUOBJ)

makedirectories:
	mkdir -p $(OBJDIR) $(BINDIR)

clean:
	rm -f $(OBJDIR)/*.o $(TARGET) $(BINDIR)/*_distfield.*

clobber: clean
	rm -rf $(OBJDIR)
