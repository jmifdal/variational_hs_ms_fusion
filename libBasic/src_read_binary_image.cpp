#include <iostream>
#include <fstream>
#include "../library/libImage.h"
#include "../library/libBasic.h"


int main(int argc, char **argv) {
    
	std::vector <OptStruct *> options;
	
	std::vector<ParStruct *> parameters;
	ParStruct pinput = {"input", NULL, "input file"}; parameters.push_back(&pinput);
	ParStruct poutput = {"output", NULL, "output file"}; parameters.push_back(&poutput);
    ParStruct plines = {"lines", NULL, "number of lines of input file"}; parameters.push_back(&plines);
    ParStruct pcolumns = {"columns", NULL, "number of columns of input file"}; parameters.push_back(&pcolumns);
    ParStruct pbytes = {"bytes", NULL, "number of bytes per pixel"}; parameters.push_back(&pbytes);
    
	
    if (!parsecmdline("src_read_binary_image", "Reads binary image", argc, argv, options, parameters))
        return 0;
    
    //! Parameters
    int lines = atoi(plines.value);
    int columns = atoi(pcolumns.value);
    int dim = lines * columns;
    int bytes = atoi(pbytes.value);
    
    //! Output
    libUSTG::flimage output(columns, lines);
    
    
    std::ifstream file(pinput.value, std::ios::binary | std::ios::in);
    
    
    
    
}

/*
    
    int *vector = new int[dim];
    
    
    FILE *out = fopen("outfile_long.txt", "w+");
    
    
    //std::streampos begin, end;
 
    long x;
    
    while(file.read((char*)&x,sizeof(x)).gcount()==sizeof(x))
    {
        fprintf(out, "%i\n",x);
    }
    
    
    file.close();
    
    return EXIT_SUCCESS;
*/