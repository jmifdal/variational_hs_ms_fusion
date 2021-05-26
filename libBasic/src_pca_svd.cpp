#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "../library/libImage.h"
#include "../library/libBasic.h"


using namespace std;

int main(int argc, char **argv)
{
	
    vector <OptStruct *> options;
	 
	vector<ParStruct *> parameters;
	ParStruct pinput = {"input", NULL, "input file"}; parameters.push_back(&pinput);
	
	if(!parsecmdline("src_pca_svd", "Apply PCA and provide the right singular vectors V^T in SUV^T decomposition", argc, argv, options, parameters))
        return EXIT_FAILURE;
	
	
	// Input
    libUSTG::cflimage input;
    input.load(pinput.value);
    
    
    // Process
    libUSTG::laVector S(input.c());
    libUSTG::laMatrix US(input.wh(), input.c()), V(input.c(), input.c()), X(input.wh(), input.c());

    for(int ii = 0; ii < input.wh(); ii++)
        for(int iC = 0; iC < input.c(); iC++)
            X[ii][iC] = input[iC * input.wh() + ii];
    
    libUSTG::compute_pca_svd(X, S, V, US);
    
    
    // Output
    for(int i = 0; i < input.c(); i++, printf("\n"))
        for(int j = 0; j < input.c(); j++)
            printf("%f \t", V[i][j]);
    
    
	return EXIT_SUCCESS;
}



