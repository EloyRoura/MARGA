#ifndef IMAGEDEFINITIONS_H
#define IMAGEDEFINITIONS_H
#include <itkImage.h>
#include "itkImageDuplicator.h"

#include <itkImageRegionIterator.h>
#include <itkNeighborhoodIterator.h>

#include <itkVector.h>
#include <itkMatrix.h>
#include <itkVariableSizeMatrix.h>
//#include <itkCompose2DVectorImageFilter.h>
//#include <itkCompose3DVectorImageFilter.h>
//#include <itkImageToVectorImageFilter.h>

#include <itkNotImageFilter.h>
#include "itkRescaleIntensityImageFilter.h"
#include "itkRelabelComponentImageFilter.h"
#include "itkScalarToRGBColormapImageFilter.h"
#include "itkConnectedComponentImageFilter.h"
#include "itkLabelShapeKeepNObjectsImageFilter.h"
#include <itkBinaryThresholdImageFilter.h>

#include <itkExtractImageFilter.h>

#include "itkCropImageFilter.h"

#include <itkSTAPLEImageFilter.h>


typedef std::string FileNameType;

typedef double PixelType;
typedef int MaskPixelType;
typedef unsigned short shortPixelType;

typedef itk::Image< PixelType, 3 > ImageType;
typedef ImageType::Pointer Image;

typedef itk::Image< PixelType, 2 > SliceType;
typedef SliceType::Pointer Slice;

typedef itk::Image< MaskPixelType, 3 > MaskImageType;
typedef MaskImageType::Pointer MaskImage;

typedef itk::Image< MaskPixelType, 2 > MaskSliceType;
typedef MaskSliceType::Pointer MaskSlice;

typedef itk::Image< shortPixelType, 2 > LabelsSliceType;
typedef LabelsSliceType::Pointer LabelsSlice;

typedef itk::BinaryThresholdImageFilter< ImageType, ImageType > ThresholdingFilterType;


typedef ImageType::IndexType IndexType;
typedef ImageType::RegionType RegionType;
typedef ImageType::SizeType SizeType;
typedef ImageType::DirectionType DirectionType;

typedef MaskSliceType::RegionType SliceRegionType;

typedef itk::ImageRegionIterator<ImageType> Iterator;
typedef itk::ImageRegionIterator<MaskImageType> MaskImageIterator;
typedef itk::NeighborhoodIterator<ImageType> NeighborhoodIterator;
typedef itk::ImageRegionIterator<SliceType> SliceIterator;
typedef itk::ImageRegionIterator<MaskSliceType> MaskSliceIterator;
typedef itk::NeighborhoodIterator<SliceType> SliceNeighborhoodIterator;
typedef itk::NeighborhoodIterator<MaskSliceType> MaskSliceNeighborhoodIterator;

//slicers
typedef itk::ExtractImageFilter< ImageType, SliceType > SlicerType;
typedef itk::ExtractImageFilter< MaskImageType, MaskSliceType > MaskSlicerType;

typedef itk::CropImageFilter <ImageType, ImageType> CropImageFilterType;

typedef MaskImageType::RegionType MaskRegionType;
typedef MaskImageType::SizeType MaskSizeType;
typedef MaskImageType::DirectionType MaskDirectionType;
typedef MaskImageType::SizeType::SizeValueType MaskSizeValueType;
typedef MaskImageType::IndexType MaskIndexType;
typedef itk::ImageRegionIterator<MaskImageType> MaskIterator;
typedef itk::NeighborhoodIterator<MaskImageType> MaskNeighborhoodIterator;

typedef itk::ImageDuplicator< MaskSliceType > DuplicatorType;

typedef itk::NotImageFilter <MaskImageType, MaskImageType> NotImageFilterType;
typedef itk::ConnectedComponentImageFilter < MaskImageType, MaskImageType > ConnectedComponentImageFilterType;
typedef itk::LabelShapeKeepNObjectsImageFilter < MaskImageType > LabelShapeKeepNObjectsImageFilterType;
typedef itk::RescaleIntensityImageFilter< MaskImageType, MaskImageType > RescaleFilterType;
typedef itk::RelabelComponentImageFilter<MaskImageType, MaskImageType> RelabelFilterType;

typedef itk::STAPLEImageFilter < ImageType, ImageType > STAPLEImageType;

enum { MI_MATTES, MI, MI_HIST, NMI, CORR, NORM_CORR, GRADIENT_DIFF, KULLBACK_LEIBLER, KAPPA, MEAN_SQUARES, MEAN_SQUARES_HIST, MEAN_SQUARES_RECIPROCAL };

template <class T> const T& min ( const T& a, const T& b ) {
  return (a<b)?a:b;     // or: return comp(a,b)?a:b; for the comp version
}
template <class T, class Tt> const T& min ( const T& a, const Tt& b ) {
  return (a<b)?a:b;     // or: return comp(a,b)?a:b; for the comp version
}
template <class T> const T& max ( const T& a, const T& b ) {
  return (a>b)?a:b;     // or: return comp(a,b)?a:b; for the comp version
}
template <class T, class Tt> const T& max ( const T& a, const Tt& b ) {
  return (a>b)?a:b;     // or: return comp(a,b)?a:b; for the comp version
}


#endif // IMAGEDEFINITIONS_H
