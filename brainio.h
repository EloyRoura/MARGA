#ifndef BRAINIO_H
#define BRAINIO_H

#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>
#include <itkConnectedComponentImageFilter.h>
#include <itkMinimumMaximumImageCalculator.h>
#include <itkBinaryDilateImageFilter.h>
#include <itkBinaryBallStructuringElement.h>
#include <itkCastImageFilter.h>

#include <imageDefinitions.h>
#include <iostream>
#include <ctime>
#include <string>

class BrainIO
{

    // Image
    typedef itk::ImageFileReader< ImageType > ReaderType;
    typedef itk::ImageFileWriter< ImageType > WriterType;

    //Slice
    typedef itk::ImageFileReader< SliceType > ReaderSliceType;
    typedef itk::ImageFileWriter< SliceType > WriterSliceType;

    // Masks
    typedef itk::ImageFileReader< MaskImageType > MaskReaderType;
    typedef itk::ImageFileWriter< MaskImageType > MaskWriterType;


    typedef itk::ImageFileReader< MaskSliceType > MaskSliceReaderType;
    typedef itk::ImageFileWriter< MaskSliceType > MaskSliceWriterType;

public:
    BrainIO();

    // Read
    Image ReadImage(FileNameType fileName);
    MaskImage ReadMaskImage(FileNameType fileName);

    // Write
    void WriteImage(FileNameType fileName, Image image);
    void WriteMaskImage(FileNameType fileName, MaskImage bimage);
    void WriteSlice(FileNameType fileName, Slice slice);
    void WriteMaskSlice(FileNameType fileName, MaskSlice slice);
};

#endif // BRAINIO_H
