
hdrdir=   -I/usr/local/include/
libdir=    -L/usr/local/lib/     

COPT	= -O3 -funroll-loops -fomit-frame-pointer  -ffast-math -ftree-vectorize
CXXFLAGS  +=  -g -no-pie -std=c++11 $(COPT) -std=c99 -Wall -Wextra  -Wno-write-strings  -Wno-deprecated -fopenmp  $(hdrdir)
LDFLAGS +=  $(CXXFLAGS) $(libdir) -lpng -ltiff -ljpeg  -lfftw3 -lfftw3f -fopenmp -lpthread -lX11

LIBMX= ../library/*.o

OBJ = src_hyperspectral_fusion_L1_weights_double.o src_hyperspectral_generate_hyperpan.o src_hyperspectral_generate_multispectral.o src_hyperspectral_fusion_L2.o src_hyperspectral_fusion_L2_weights.o src_hyperspectral_fusion_L1_weights.o src_hyperspectral_fusion_L1.o src_hyperspectral_fusion_L1_double.o
BIN = src_hyperspectral_fusion_L1_weights_double src_hyperspectral_generate_hyperpan src_hyperspectral_generate_multispectral src_hyperspectral_fusion_L2 src_hyperspectral_fusion_L2_weights src_hyperspectral_fusion_L1_weights src_hyperspectral_fusion_L1 src_hyperspectral_fusion_L1_double

LIBBIN=../bin

all: $(OBJ)	$(BIN) 


$(OBJ) : %.o : %.cpp
	$(CXX) -c    $< -o $@ $(CXXFLAGS)


$(BIN) : % : %.o  $(LIBMX)
	$(CXX)   -o $(LIBBIN)/$@  $^ $(LDFLAGS)



.PHONY : clean
clean:
	$(RM) $(OBJ) ; cd $(LIBBIN); rm $(BIN)

    
