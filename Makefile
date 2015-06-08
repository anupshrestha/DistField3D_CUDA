EXEC := GIN3D-DFP
SRCDIR := src
BINDIR := bin
OBJDIR := obj
SERDIR := src/serial
PARDIR := src/cuda

OBJ := $(patsubst $(SRCDIR)/%.c,$(OBJDIR)/%.o,$(wildcard $(SRCDIR)/*.c))
SEROBJ := $(patsubst $(SERDIR)/%.c,$(OBJDIR)/%.o,$(wildcard $(SERDIR)/*.c))
PAROBJ := $(patsubst $(PARDIR)/%.cu,$(OBJDIR)/%.o,$(wildcard $(PARDIR)/*.cu))
TARGET :=  $(BINDIR)/$(EXEC)

# Compiler options
CC = gcc -Wall -g -std=c99
CCFLAGS = -I/usr/local/wrfUtils/netcdf-c/include -I/usr/local/wrfUtils/szip/include -I/usr/local/wrfUtils/hdf5/include
LDFLAGS = -L/usr/local/wrfUtils/netcdf-c/lib -L/usr/local/wrfUtils/hdf5/lib -L/usr/local/wrfUtils/szip/lib -lm -lhdf5 -lsz -lnetcdf

NVCC = nvcc 
NVCCFLAGS = -arch=sm_20

parallel: makedirectories $(OBJ) $(PAROBJ)
	$(NVCC) $(NVCCFLAGS) $(CCFLAGS) $(LDFLAGS) -o $(TARGET) $(OBJ) $(PAROBJ)
	
serial: makedirectories $(OBJ) $(SEROBJ)
	$(CC) $(CCFLAGS) $(LDFLAGS) -o $(TARGET) $(OBJ) $(SEROBJ)

$(OBJDIR)/%.o: $(SRCDIR)/%.c
	$(CC) $(CCFLAGS) -o $@ -c $< $(LDFLAGS)

$(OBJDIR)/%.o: $(SERDIR)/%.c
	$(CC) $(CCFLAGS) -o $@ -c $< $(LDFLAGS)

$(OBJDIR)/%.o: $(PARDIR)/%.cu
	$(NVCC) $(NVCCFLAGS) $(CCFLAGS) $(LDFLAGS) -o $@ -c $<

makedirectories:
	mkdir -p $(OBJDIR) $(BINDIR)

clean:
	rm -f $(OBJDIR)/*.o $(TARGET) $(BINDIR)/*_distfield.*

clobber: clean
	rm -rf $(OBJDIR)
