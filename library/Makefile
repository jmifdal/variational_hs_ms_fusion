
hdrdir=	-I/usr/local/include/
libdir=	-L/usr/local/lib/	

COPT	= -lm -O3 -funroll-loops -fomit-frame-pointer -ffast-math -ftree-vectorize
CXXFLAGS 	+= 	-g -std=c++11 $(COPT)   -pedantic -Wall -Wextra    -Wno-write-strings  -Wno-deprecated  $(hdrdir) -fopenmp
LDFLAGS	+=	-g $(CXXFLAGS) $(libdir) -lpng -ltiff -ljpeg -lfftw3 -lfftw3f -fopenmp


CFLAGS	=	-O3 -std=c99 -I/usr/local/include/   -L/usr/local/lib/

OBJC	=	libImageFormats.o
OBJ	=	libBasic.o libImageFormatPM.o libImage.o libHyperspectral.o libHyperspectralDOUBLE.o libHyperspectralWeights.o libHyperspectralWeightsDouble.o libFilteringNLmeans.o libFilteringVariational.o 

all:	$(OBJC)	$(OBJ)	$(BIN)

$(OBJ) : %.o : %.cpp
	$(CXX) -c $(CXXFLAGS)   $< -o $@


$(BIN) : % : %.o
	$(CXX) $(LDFLAGS)  -o $(LIBBIN)/$@  $^

.PHONY : clean
clean:
	$(RM) $(OBJ) $(OBJC) ;
