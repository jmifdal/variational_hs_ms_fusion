

hdrdir=   -I/usr/local/include/
libdir=    -L/usr/local/lib/     

COPT	= -O3 -funroll-loops -fomit-frame-pointer  -ffast-math -ftree-vectorize
CXXFLAGS  +=  -g -no-pie -std=c++11 $(COPT)     -Wall -Wextra  -Wno-write-strings  -Wno-deprecated -fopenmp  $(hdrdir)
LDFLAGS +=  $(CXXFLAGS) $(libdir) -lpng -ltiff -ljpeg  -lfftw3 -lfftw3f -fopenmp -lpthread -lX11

LIBMX= ../library/*.o

OBJ = src_filter_convol.o  src_filter_convol_gaussian.o src_create_gaussian_kernel.o src_filter_nlmeans.o src_filter_patch_max.o src_filter_patch_median.o src_filter_fft_convol_gaussian.o src_filter_variational_TVNL_L2w.o

BIN = src_filter_convol src_filter_convol_gaussian  src_create_gaussian_kernel src_filter_nlmeans src_filter_patch_max src_filter_patch_median src_filter_fft_convol_gaussian src_filter_variational_TVNL_L2w

LIBBIN=../bin

all: $(OBJ)	$(BIN) 


$(OBJ) : %.o : %.cpp
	$(CXX) -c    $< -o $@ $(CXXFLAGS)


$(BIN) : % : %.o  $(LIBMX)
	$(CXX)   -o $(LIBBIN)/$@  $^ $(LDFLAGS)



.PHONY : clean
clean:
	$(RM) $(OBJ) ; cd $(LIBBIN); rm $(BIN)





