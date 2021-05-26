hdrdir=   -I/usr/local/include/
libdir=    -L/usr/local/lib/     

COPT	= -O3 -funroll-loops -fomit-frame-pointer  -ffast-math -ftree-vectorize
CXXFLAGS  +=  -g -no-pie -std=c++11 $(COPT)     -Wall -Wextra  -Wno-write-strings  -Wno-deprecated -fopenmp  $(hdrdir)
LDFLAGS +=  $(CXXFLAGS) $(libdir) -lpng -ltiff -ljpeg  -lfftw3 -lfftw3f -fopenmp -lpthread -lX11

LIBMX= ../library/*.o

OBJ = src_spline_zoom.o  src_fft.o	src_fft_sample.o src_fft_translate.o src_sample.o src_fft_rot.o	src_fft_shear.o	src_fft_zoom.o src_transform_homography.o src_spline_translate.o

BIN = src_spline_zoom  src_fft	src_fft_sample src_fft_translate src_sample src_fft_rot src_fft_shear src_fft_zoom src_transform_homography src_spline_translate

LIBBIN=../bin

all: $(OBJ)	$(BIN) 


$(OBJ) : %.o : %.cpp
	$(CXX) -c    $< -o $@ $(CXXFLAGS)


$(BIN) : % : %.o  $(LIBMX)
	$(CXX)   -o $(LIBBIN)/$@  $^ $(LDFLAGS)



.PHONY : clean
clean:
	$(RM) $(OBJ) ; cd $(LIBBIN); rm $(BIN)





