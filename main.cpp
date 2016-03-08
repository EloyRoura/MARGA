/**
* MARGA: Multispectral Adaptive Region Growing Algorithm for brain extraction on axial MRI
*
*
* Implemented by Eloy Roura Perez
* mailto: eloyroura@eia.udg.edu
*/

#include <iostream>
#include <cstdlib>

// ITK IO includes
//#include "itkOrientedImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"

// ITK Registration includes 
#include "itkCastImageFilter.h" 
#include "itkWarpImageFilter.h" 
#include "itkLinearInterpolateImageFunction.h"

#include "itkCheckerBoardImageFilter.h"

//Own includes
#include <brainio.h>
#include <brainExtraction.h>

typedef struct{
        char *inputImageName;
        char *outputResultName;
        float t_lower;
        float t_upper;
        int minSlice;
        int maxSlice;
        int minSm;
        int maxSm;
        int dilateRadius;
        int erodeRadius;
}PARAM;

typedef struct{
        bool outputResultFlag;
        bool t_lowerFlag;
        bool t_upperFlag;
        bool minSliceFlag;
        bool maxSliceFlag;
        bool minSmFlag;
        bool maxSmFlag;
        bool dilateRadiusFlag;
        bool erodeRadiusFlag;
        bool inputImageFlag;
}FLAG;

void Usage(char *exec){
    std::cout << "********************************************************************************" << std::endl;
    std::cerr << "Usage: " << std::endl;
    std::cerr << "\tSkullStrip -input <filename> [OPTIONS] " << std::endl;
    std::cerr << "" <<std::endl;
    std::cerr << "\tOPTIONS:" <<std::endl;
    std::cerr << "\t\t -output <filename> -> mask image filename (default: MARGA_mask.nii)" <<std::endl;
    std::cerr << "\t\t -t_lower <float> -> Lower threshold (default: 0.35)" <<std::endl;
    std::cerr << "\t\t -t_upper <float> -> Upper threshold (default: 0.80)" <<std::endl;
    std::cerr << "\t\t -minSlice <int> ->  Minimum slice where to cut the process (default: 0)" <<std::endl;
    std::cerr << "\t\t -maxSlice <int> ->  Maximum slice where to cut the process <numberOfSlices - maxSlice> (default: 0)" <<std::endl;
    std::cerr << "\t\t -minSm <int> ->  Minimum slice where to reduce the threshold in 0.03 (default: 0)" <<std::endl;
    std::cerr << "\t\t -maxSm <int> ->  Maximum slice where to reduce the threshold in 0.03 <numberOfSlices - maxSm> (default: 0)" <<std::endl;
    std::cerr << "\t\t -dilateRadius <int> -> Radius in mm used in the dilation morphological operation (default: 5mm)" <<std::endl;
    std::cerr << "\t\t -erodeRadius <int> ->  Radius in mm used in the erosion morphological operation (default: 5mm)" <<std::endl;
    std::cout << std::endl << "********************************************************************************" << std::endl;
    std::cerr << "\tDefault command: " << std::endl;
    std::cerr << "\t\tSkullStrip -input T1.nii -output MARGA_mask.nii -t_lower 0.35 -t_upper 0.8 -minSlice 0 -maxSlice 0 -minSm 0 -maxSm 0 -dilateRadius 5 -erodeRadius 5" << std::endl;
    std::cout << std::endl << "********************************************************************************" << std::endl;
}


int main( int argc, char * argv [] )
{
    PARAM *param = (PARAM *)calloc(1,sizeof(PARAM));
    FLAG *flag = (FLAG *)calloc(1,sizeof(FLAG));

    std::cout << std::endl << "********************************************************************************" << std::endl << std::endl;

    // Verify the number of parameters in the command line
    if( argc < 2)
    {
        Usage(argv[0]);
        return EXIT_FAILURE;
    }


    /* read the input parameter */
    for(int i=1;i<argc;i++){
            if(strcmp(argv[i], "-help")==0 || strcmp(argv[i], "-Help")==0 ||
                            strcmp(argv[i], "-h")==0 || strcmp(argv[i], "-HELP")==0 ||
                            strcmp(argv[i], "--h")==0 || strcmp(argv[i], "--help")==0){
                Usage(argv[0]);
                return 0;
            }
            else if(strcmp(argv[i], "-input") == 0){
                param->inputImageName=argv[++i];
                flag->inputImageFlag=1;
            }
            else if(strcmp(argv[i], "-output") == 0){
                param->outputResultName=argv[++i];
                flag->outputResultFlag=1;
            }
            else if(strcmp(argv[i], "-t_lower") == 0){
                param->t_lower = atof(argv[++i]);
                flag->t_lowerFlag=1;
            }
            else if(strcmp(argv[i], "-t_upper") == 0){
                param->t_upper = atof(argv[++i]);
                flag->t_upperFlag=1;
            }
            else if(strcmp(argv[i], "-minSlice") == 0){
                param->minSlice=atoi(argv[++i]);
                flag->minSliceFlag=1;
            }
            else if(strcmp(argv[i], "-maxSlice") == 0){
                param->maxSlice=atoi(argv[++i]);
                flag->maxSliceFlag=1;
            }
            else if(strcmp(argv[i], "-maxSm") == 0){
                param->maxSm=atoi(argv[++i]);
                flag->maxSmFlag=1;
            }
            else if(strcmp(argv[i], "-minSm") == 0){
                param->minSm=atoi(argv[++i]);
                flag->minSmFlag=1;
            }
            else if(strcmp(argv[i], "-dilateRadius") == 0){
                param->dilateRadius=atoi(argv[++i]);
                flag->dilateRadiusFlag=1;
            }
            else if(strcmp(argv[i], "-erodeRadius") == 0){
                param->erodeRadius=atoi(argv[++i]);
                flag->erodeRadiusFlag=1;
            }
    }

    if(!flag->inputImageFlag){
            fprintf(stderr,"ERROR:\tThe input images have to be defined.\n");
            return EXIT_FAILURE;
    }

    /*Initialization of the default parameters*/
    if(!flag->t_lowerFlag)
        param->t_lower = 0.35;

    if(!flag->t_upperFlag)
        param->t_upper = 0.8;

    if(!flag->minSliceFlag)
        param->minSlice = 0;

    if(!flag->maxSliceFlag)
        param->maxSlice = 0;

    if(!flag->minSmFlag)
        param->minSm = 0;

    if(!flag->maxSmFlag)
        param->maxSm = 0;

    if(!flag->dilateRadiusFlag)
        param->dilateRadius = 5;

    if(!flag->erodeRadiusFlag)
        param->erodeRadius = 5;
    /*********************************************************************************************************************************************************************/


    std::cout << "\t\t - Skull Striping (based on Region Growing) - " << std::endl << std::endl;

    std::cout << "Extracting the brain from: " << param->inputImageName << std::endl;

    BrainIO* brainio = new BrainIO();
    Image inputImage = brainio->ReadImage(param->inputImageName);

    BrainExtraction* brainExtract = new BrainExtraction();
    MaskImage rgMask = brainExtract->BE_RR(inputImage, param->t_lower, param->t_upper, param->minSlice, param->maxSlice, param->minSm, param->maxSm, param->dilateRadius, param->erodeRadius);

    brainio->WriteMaskImage(param->outputResultName,rgMask);


    std::cout<<"The mask has been saved in: " << param->outputResultName << std::endl;
    std::cout << std::endl << "********************************************************************************" << std::endl;

    return EXIT_SUCCESS;
}



