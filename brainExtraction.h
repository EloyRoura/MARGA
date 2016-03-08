#ifndef BRAINEXTRACTION_H
#define BRAINEXTRACTION_H

#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>
#include <itkMinimumMaximumImageCalculator.h>
#include <itkBinaryDilateImageFilter.h>
#include <itkBinaryBallStructuringElement.h>
#include <itkCastImageFilter.h>
#include <itkBinaryThresholdImageFilter.h>
#include "itkOtsuThresholdImageFilter.h"
#include <itkVotingBinaryHoleFillingImageFilter.h>
#include <itkVotingBinaryIterativeHoleFillingImageFilter.h>
#include "itkBinaryFunctorImageFilter.h"
#include "itkAndImageFilter.h"
#include "itkOrImageFilter.h"
#include <itkNotImageFilter.h>

#include "itkBinaryBallStructuringElement.h"
#include "itkGrayscaleErodeImageFilter.h"
#include "itkGrayscaleDilateImageFilter.h"
#include "itkRescaleIntensityImageFilter.h"
#include "itkRelabelComponentImageFilter.h"
#include "itkScalarToRGBColormapImageFilter.h"
#include "itkConnectedComponentImageFilter.h"
#include "itkLabelShapeKeepNObjectsImageFilter.h"
#include "itkSubtractImageFilter.h"

#include "itkOrientImageFilter.h"


#include "itkScalarImageToHistogramGenerator.h"

#include <brainio.h>
#include <imageDefinitions.h>
#include <iostream>
#include <ctime>
#include <string>

class BrainExtraction
{

    typedef itk::MinimumMaximumImageCalculator< ImageType > MinimumMaximumImageCalculatorType;
    typedef itk::MinimumMaximumImageCalculator< SliceType > MinimumMaximumSliceCalculatorType;

    typedef itk::Statistics::ScalarImageToHistogramGenerator< ImageType > HistogramGeneratorType;
    typedef itk::Statistics::ScalarImageToHistogramGenerator< SliceType > SliceHistogramGeneratorType;
    typedef HistogramGeneratorType::HistogramType HistogramType;
    typedef SliceHistogramGeneratorType::HistogramType SliceHistogramType;

    typedef itk::MaskImageFilter< ImageType, MaskImageType, ImageType > MaskImageFilterType;
    typedef itk::MaskImageFilter< SliceType, MaskSliceType, SliceType > MaskSliceFilterType;
    typedef itk::MaskImageFilter< MaskSliceType, MaskSliceType, MaskSliceType > MaskSliceMaskFilterType;

    typedef itk::BinaryThresholdImageFilter< ImageType, MaskImageType > ThresholdingFilterType;
    typedef itk::BinaryThresholdImageFilter< SliceType, MaskSliceType > ThresholdingSliceFilterType;
    typedef itk::BinaryThresholdImageFilter< MaskSliceType, MaskSliceType > ThresholdingMaskSliceFilterType;
    typedef itk::OtsuThresholdImageFilter < SliceType, MaskSliceType > OtsuThresholdImageFilterType;
    typedef itk::OrImageFilter <MaskSliceType> OrImageFilterType;
    typedef itk::NotImageFilter <MaskSliceType, MaskSliceType> NotSliceFilterType;
    typedef itk::NotImageFilter <MaskImageType, MaskImageType> NotImageFilterType;

    typedef itk::VotingBinaryHoleFillingImageFilter< MaskImageType, MaskImageType > HoleFillingImageFilterType;
    typedef itk::VotingBinaryHoleFillingImageFilter< MaskSliceType, MaskSliceType > HoleFillingSliceFilterType;

    typedef itk::BinaryBallStructuringElement < MaskSliceType::PixelType, 2> StructuringElementType;
    typedef itk::GrayscaleErodeImageFilter <MaskSliceType, MaskSliceType, StructuringElementType> GrayscaleErodeSliceFilterType;
    typedef itk::GrayscaleDilateImageFilter <MaskSliceType, MaskSliceType, StructuringElementType> GrayscaleDilateSliceFilterType;


    typedef itk::BinaryBallStructuringElement < MaskImageType::PixelType, 3> Structuring3DElementType;
    typedef itk::GrayscaleErodeImageFilter <MaskImageType, MaskImageType, Structuring3DElementType> GrayscaleErodeImageFilterType;
    typedef itk::GrayscaleDilateImageFilter <MaskImageType, MaskImageType, Structuring3DElementType> GrayscaleDilateImageFilterType;


    typedef itk::ConnectedComponentImageFilter < MaskSliceType, MaskSliceType > ConnectedComponentImageFilterType;
    typedef itk::LabelShapeKeepNObjectsImageFilter < MaskSliceType > LabelShapeKeepNObjectsImageFilterType;
    typedef itk::RescaleIntensityImageFilter< MaskSliceType, MaskSliceType > RescaleFilterType;
    typedef itk::RelabelComponentImageFilter<MaskSliceType, MaskSliceType> RelabelFilterType;

    typedef itk::SubtractImageFilter < MaskSliceType, MaskSliceType > SubtractImageFilterType;
    typedef itk::OrientImageFilter < MaskImageType,MaskImageType > OrientImageFilterType;


public:
    BrainExtraction();

    MaskImage BE_RR(Image inputImage, float t_lower, float t_upper, int minSlice, int maxSlice, int minSmSlice, int maxSmSlice, int dilateRadius, int erodeRadius);
    const HistogramGeneratorType::Pointer getHistogram(Image inputImage);
    const double * getCumulativeHistogram(Image imROI, int max, int min);
    const double * getCumulativeHistogramSlice(Slice sliceROI, int max, int min);
    unsigned int lowerHist(const HistogramType * histogram);
    unsigned int upperHist(const HistogramType * histogram);
    unsigned int getSizeM(MaskSlice M);
    MaskImage getRoi(Image inputImage, float t2, float t98);
    MaskSlice otsuThreshold(Slice slice);
    MaskSlice imThreshold(Slice slice, double threshold);
    Image maskImage(Image inputImage, MaskImage roi);
    Slice maskSlice(Slice inputSlice, MaskSlice roi);
    MaskSlice maskMaskSlice(MaskSlice inputSlice, MaskSlice M);
    double meanImage(Image imROI);
    double getMeanSlice(Slice slice);
    Slice getSlice(Image inputImage, int slice);
    MaskSlice getSlice(MaskImage inputImage, int slice);
    MaskImage setSlice(MaskImage image, MaskSlice slice, int nSlice);
    MaskSlice getM(Slice slice);
    double getTlower(float t_lower, Slice slice, float meanROI);
    double getTupper(float t_upper, Slice img_inter);
    MaskSlice erodeSlice(MaskSlice inputSl, unsigned int radius);
    MaskImage erodeIm(MaskImage inputIm, unsigned int radius);
    MaskSlice fillSlice(MaskSlice slice);
    MaskSlice getLergestConnectedComponent(MaskSlice maskSl);
    MaskSlice dilateSlice(MaskSlice inputSl, unsigned int radius);
    MaskImage dilateIm(MaskImage inputIm, unsigned int radius);
    MaskSlice seedRegion(MaskSlice thrISlice, MaskSlice M);
    MaskSlice voxelAssignment(MaskSlice assignedVoxels, MaskSlice newVoxels, MaskSlice M);
    MaskSlice componentAssignment(MaskSlice assignedVoxels, MaskSlice M);
    MaskSlice imageSubtract(MaskSlice image1, MaskSlice image2);
    MaskSlice extractPixelsWithValue(MaskSlice slice, int pValue);
    MaskImage removeBrainStem(MaskImage image);
};

#endif // BRAINIO_H
