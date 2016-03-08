/**
* Statistics class
*/

#include <iostream>
#include <cstdlib>

// ITK IO includes
//#include "itkOrientedImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"

//My classes includes
#include "brainio.h"
#include "imageDefinitions.h"

int main( int argc, char * argv [] )
{

  // Verify the number of parameters in the command line
  if( argc < 3 )
    {
    std::cout<<std::endl<< "************************************************************************************************"<<std::endl<<std::endl;
    std::cerr << "Usage: " << std::endl;
    std::cerr << "\t" << argv[0] << " SourceImage Output" << std::endl;
    std::cerr << std::endl << "************************************************************************************************"<<std::endl<<std::endl;
    return EXIT_FAILURE;
    }



  BrainIO* brainio = new BrainIO();
  MaskImage sourceImage = brainio->ReadMaskImage(argv[1]);


  //Inverse the original mask
  NotImageFilterType::Pointer filterInvers1 = NotImageFilterType::New();
  filterInvers1->SetInput(sourceImage);
  filterInvers1->Update();



  //Getting the largest connected component (e.g. the background)
  ConnectedComponentImageFilterType::Pointer labelFilter = ConnectedComponentImageFilterType::New ();
  labelFilter->SetInput(filterInvers1->GetOutput());
  labelFilter->Update();

  LabelShapeKeepNObjectsImageFilterType::Pointer labelShapeKeepNObjectsImageFilter = LabelShapeKeepNObjectsImageFilterType::New();
  labelShapeKeepNObjectsImageFilter->SetInput( labelFilter->GetOutput() );
  labelShapeKeepNObjectsImageFilter->SetBackgroundValue( 0 );
  labelShapeKeepNObjectsImageFilter->SetNumberOfObjects( 1 );
  labelShapeKeepNObjectsImageFilter->SetAttribute( LabelShapeKeepNObjectsImageFilterType::LabelObjectType::NUMBER_OF_PIXELS);
  labelShapeKeepNObjectsImageFilter->Update();

  RescaleFilterType::Pointer rescaleFilter = RescaleFilterType::New();
  rescaleFilter->SetOutputMinimum(0);
  //rescaleFilter->SetOutputMaximum(itk::NumericTraits<PixelType>::max());
  rescaleFilter->SetOutputMaximum(255);
  rescaleFilter->SetInput(labelShapeKeepNObjectsImageFilter->GetOutput());
  rescaleFilter->Update();



  //Inverse the background
  NotImageFilterType::Pointer filterInvers2 = NotImageFilterType::New();
  filterInvers2->SetInput(rescaleFilter->GetOutput());
  filterInvers2->Update();



  brainio->WriteMaskImage(argv[2],filterInvers2->GetOutput());

  return EXIT_SUCCESS;
}

