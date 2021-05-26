#include "../library/libBasic.h"
#include "../library/libImage.h"



int main(int argc, char **argv)
{
    std::vector <OptStruct *> options;
	
    std::vector<ParStruct *> parameters;
	ParStruct pinput = {"input", NULL, "input txt file"}; parameters.push_back(&pinput);
	ParStruct pout = {"output", NULL, "output image"}; parameters.push_back(&pout);
    ParStruct pwidth = {"width", NULL, "width"}; parameters.push_back(&pwidth);
    ParStruct pheight = {"height", NULL, "height"}; parameters.push_back(&pheight);

	if(!parsecmdline("ustg_convert_txt2image", "Convert txt to image", argc, argv, options, parameters))
		return EXIT_FAILURE;

    // Image size
    int width = atoi(pwidth.value);
    int height = atoi(pheight.value);

    
	// Read input info
 	std::ifstream f(pinput.value);
	if(!f.good())
    {
        printf("error :: file not found or unreadable");
        return EXIT_FAILURE;
    }
	
	// Output memory
	libUSTG::flimage out(width, height);
	
	
	for (int i = 0; i < out.whc(); i++)
		f >> out[i];
	
	f.close();
    out.save(pout.value);
	
    
    return EXIT_SUCCESS;
}
