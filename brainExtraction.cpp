/******************************************************


  Methods specification:

Get M:
    1. ThresholdingFilterType
    2. Erode
    3. largest component
    4. dilate
    5. Fill holes

Assignment Rules (B->3 U->2 N->1):
1. Voxel

2. Components

Contact: Eloy Roura  mailto:eloyroura@eia.udg.edu
******************************************************/



#include <brainExtraction.h>
//#include "QuickView.h"

//using namespace BrainExtraction;

BrainExtraction::BrainExtraction()
{
}

MaskImage BrainExtraction::BE_RR(Image inputImage, float t_lower, float t_upper, int minSlice, int maxSlice, int minSmSlice, int maxSmSlice, int dilateRadius, int erodeRadius) {
    /***************** Init *******************/
    int nSlices, slice, t2, t98, anteriorSize=0;
    double iTupper, iTlower, increment;
    MaskSlice M, thrISlice, seedsISlice, assignedVoxels, assignedComponents, voxelsProcessed, newVoxels;
    Slice middleSlice, sliceBackward, sliceForward;
    BrainIO* brainio = new BrainIO();

    //Create the output image
    MaskImageType::Pointer rgMask = MaskImageType::New();
    rgMask->SetRegions(inputImage->GetLargestPossibleRegion());
    rgMask->SetSpacing(inputImage->GetSpacing());
    rgMask->SetOrigin(inputImage->GetOrigin());
    rgMask->SetDirection(inputImage->GetDirection());
    rgMask->Allocate();
    rgMask->FillBuffer(0);
    rgMask->Update();

    nSlices = inputImage->GetLargestPossibleRegion().GetSize()[2];
    //std::cout << "The volume contains " << nSlices << " slices" << std::endl;

    HistogramGeneratorType::Pointer histogramGenerator = getHistogram(inputImage);
    const HistogramType * histogram = histogramGenerator->GetOutput();

    /* Determination of ROI */
    // 98% and 2% estimation from the cumulative histogram
    t2 = lowerHist(histogram);
    t98 = upperHist(histogram);

    MaskImageType::Pointer roi = getRoi(inputImage,t2,t98);
    brainio->WriteMaskImage("roiITK.nii.gz",roi);
    ImageType::Pointer imROI = maskImage(inputImage, roi);


    brainio->WriteImage("imRoi.nii.gz", imROI);

    double meanROI = meanImage(imROI);
    //std::cout << "Mean ROI = " << meanROI << std::endl;

    increment = (t_upper-t_lower) / (maxSmSlice - maxSlice);
    std::cout << "increment: " << increment << std::endl;


    /***************** Initial Slice Processing *******************/
    slice = nSlices/2;
    //slice=165;

    std::cout << "Middle slice = " << slice << std::endl;
    std::cout << "\tThreshold t2-t98: [" << t2 << "-" << t98 << "]" << std::endl;
    middleSlice = getSlice(imROI, slice);

    M = getM(middleSlice);

    Slice img_inter = maskSlice(middleSlice, M);
    MinimumMaximumSliceCalculatorType::Pointer maxIntensityCalculator = MinimumMaximumSliceCalculatorType::New();
    maxIntensityCalculator->SetImage( img_inter );
    maxIntensityCalculator->Compute();


    if (maxIntensityCalculator->GetMaximum()>0){       

        iTlower = getTlower(t_lower, middleSlice, meanROI);
        iTupper = getTupper(t_upper, img_inter);

        std::cout << "\tThreshold Up2Low: [" << iTupper << "-" << iTlower << "]" << std::endl;

        thrISlice = imThreshold(middleSlice,iTupper);
        seedsISlice = seedRegion(thrISlice, M);
        assignedComponents = seedsISlice;

        for(int dinamicTupper=iTupper-1; dinamicTupper>=iTlower; dinamicTupper--){
            voxelsProcessed = thrISlice;
            thrISlice = imThreshold(middleSlice,dinamicTupper);
            newVoxels = imageSubtract(thrISlice, voxelsProcessed);

            assignedVoxels = voxelAssignment(assignedComponents, newVoxels, M);
            assignedComponents = componentAssignment(assignedVoxels,M);
        }
        M = extractPixelsWithValue(assignedComponents, 3);
        M = dilateSlice(M, dilateRadius);
        M = erodeSlice(M, erodeRadius);
        M = fillSlice(M);
    }

    DuplicatorType::Pointer duplicator = DuplicatorType::New();
    duplicator->SetInputImage(M);
    duplicator->Update();
    MaskSliceType::Pointer initialM = duplicator->GetOutput();

    rgMask = setSlice(rgMask, M, slice);
    //std::cout << "Initial slice processed" << std::endl;
    /******************************************/


    /***************** Backward *******************/
    float originalT_Lower = t_lower;
    for (int nSlice = slice-1; nSlice>=minSlice; nSlice--){
        sliceBackward = getSlice(imROI, nSlice);


        //anteriorSize = size de M
        anteriorSize = getSizeM(M);

        MinimumMaximumSliceCalculatorType::Pointer maxIntensityCalculator = MinimumMaximumSliceCalculatorType::New();
        maxIntensityCalculator->SetImage( sliceBackward );
        maxIntensityCalculator->Compute();

        Slice img_inter = maskSlice(sliceBackward, M);

        MinimumMaximumSliceCalculatorType::Pointer maxIntensityCalculator2 = MinimumMaximumSliceCalculatorType::New();
        maxIntensityCalculator2->SetImage( img_inter );
        maxIntensityCalculator2->Compute();

        if (maxIntensityCalculator->GetMaximum()>0 && maxIntensityCalculator2->GetMaximum()>0){

            iTlower = getTlower(t_lower, sliceBackward, meanROI);

            std::cout << "Slice: " << nSlice << std::endl;
            std::cout << "\tt_lower: " << t_lower << std::endl;

            if(iTlower<0){
                iTlower=0;
            }

                iTupper = getTupper(t_upper, img_inter);

                std::cout << "\tThreshold Up2Low: [" << iTupper << "-" << iTlower << "]" << std::endl;

                thrISlice = imThreshold(sliceBackward,iTupper);
                assignedComponents = seedRegion(thrISlice, M);

                for(int dinamicTupper=iTupper-1;dinamicTupper>=iTlower;dinamicTupper--){
                    voxelsProcessed = thrISlice;
                    thrISlice = imThreshold(sliceBackward,dinamicTupper);
                    newVoxels = imageSubtract(thrISlice, voxelsProcessed);

                    assignedVoxels = voxelAssignment(assignedComponents, newVoxels, M);
                    assignedComponents = componentAssignment(assignedVoxels, M);
                }
                M = extractPixelsWithValue(assignedComponents, 3);
                M = dilateSlice(M, dilateRadius);
                M = erodeSlice(M, erodeRadius);
                M = getLergestConnectedComponent(M);
                //if(nSlice != 100)
                M = fillSlice(M);

                //if(nSlice <= 79 && nSlice >= 77)
                    //std::cout << t_lower << std::endl;

                /*if(nSlice <= minSmSlice && t_lower < (t_upper - 0.03)){
                    t_lower = t_lower + 0.03;
                }*/

                if(nSlice <= minSmSlice && t_lower < (t_upper - increment)){
                    t_lower = t_lower + increment;
                }

                //if (anteriorSize > sizeM*0.80)
                /*if (anteriorSize > getSizeM(M)*0.75)
                    rgMask = setSlice(rgMask, M, nSlice);
                else
                    break;*/

                //rgMask = setSlice(rgMask, M, nSlice);
                rgMask = setSlice(rgMask, maskMaskSlice(M, getSlice(rgMask, nSlice+1)), nSlice);
            //}
        }
    }

    t_lower = originalT_Lower;
    //std::cout << "Backward processed" << std::endl;
    /******************************************/




    /***************** Forward *******************/
    M = initialM;

    for (int nSlice = slice+1; nSlice<nSlices-maxSlice; nSlice++){
        sliceForward = getSlice(imROI, nSlice);

        //anteriorSize = size de M
        anteriorSize = getSizeM(M);

        MinimumMaximumSliceCalculatorType::Pointer maxIntensityCalculator = MinimumMaximumSliceCalculatorType::New();
        maxIntensityCalculator->SetImage( sliceForward );
        maxIntensityCalculator->Compute();

        Slice img_inter = maskSlice(sliceForward, M);
        MinimumMaximumSliceCalculatorType::Pointer maxIntensityCalculator2 = MinimumMaximumSliceCalculatorType::New();
        maxIntensityCalculator2->SetImage( img_inter );
        maxIntensityCalculator2->Compute();

        if (maxIntensityCalculator->GetMaximum()>0 && maxIntensityCalculator2->GetMaximum()>0){
            iTlower = getTlower(t_lower, sliceForward, meanROI);


            if(iTlower<0){
                iTlower=0;
            }

                iTupper = getTupper(t_upper, img_inter);

                //std::cout << "\tThreshold Up2Low "<< nSlice << ": [" << iTupper << "-" << iTlower << "(" << t_lower << ")]" << std::endl;

                thrISlice = imThreshold(sliceForward,iTupper);
                assignedComponents = seedRegion(thrISlice, M);

                for(int dinamicTupper=iTupper-1;dinamicTupper>=iTlower;dinamicTupper--){
                   voxelsProcessed = thrISlice;
                   thrISlice = imThreshold(sliceForward,dinamicTupper);
                   newVoxels = imageSubtract(thrISlice, voxelsProcessed);

                   assignedVoxels = voxelAssignment(assignedComponents, newVoxels, M);
                   assignedComponents = componentAssignment(assignedVoxels, M);
                }
                M = extractPixelsWithValue(assignedComponents, 3);
                M = dilateSlice(M, dilateRadius);
                M = erodeSlice(M, erodeRadius);
                M = fillSlice(M);

                /*if(nSlice >= nSlices-maxSmSlice && t_lower < (t_upper - 0.03)){
                    t_lower = t_lower + 0.03;
                }*/

                if(nSlice >= nSlices-maxSmSlice && t_lower < (t_upper - increment)){
                    t_lower = t_lower + increment;
                }

                //rgMask = setSlice(rgMask, M, nSlice);

                //if (anteriorSize > sizeM*0.80)
                /*if (anteriorSize > getSizeM(M)*0.75)
                    rgMask = setSlice(rgMask, M, nSlice);
                else
                    break;*/

                rgMask = setSlice(rgMask, maskMaskSlice(M, getSlice(rgMask, nSlice-1)), nSlice);
                //rgMask = setSlice(rgMask, M, nSlice);
            //}
        }
    }

    //std::cout << "Forward processed" << std::endl;
    /******************************************/

    rgMask = removeBrainStem(rgMask);
    return rgMask;
}


unsigned int BrainExtraction::getSizeM(MaskSlice M){
    int nVoxels=0;
    MaskSliceIterator itMask = MaskSliceIterator( M, M->GetRequestedRegion() );
    for (itMask.GoToBegin(); !itMask.IsAtEnd(); ++itMask ) {
        if(itMask.Get() > 0)
            nVoxels++;
    }
    return nVoxels;
}

const BrainExtraction::HistogramGeneratorType::Pointer BrainExtraction::getHistogram(Image inputImage){
    MinimumMaximumImageCalculatorType::Pointer maxIntensityCalculator = MinimumMaximumImageCalculatorType::New();
    maxIntensityCalculator->SetImage(inputImage);
    maxIntensityCalculator->Compute();

    HistogramGeneratorType::Pointer histogramGenerator = HistogramGeneratorType::New();
    histogramGenerator->SetInput(inputImage);
    histogramGenerator->SetNumberOfBins( maxIntensityCalculator->GetMaximum() );
    histogramGenerator->SetHistogramMin( maxIntensityCalculator->GetMinimum() );
    histogramGenerator->SetHistogramMax( maxIntensityCalculator->GetMaximum() );
    histogramGenerator->Compute();

    maxIntensityCalculator = NULL;

    return histogramGenerator;
}


unsigned int BrainExtraction::lowerHist(const HistogramType * histogram){
    const float lowerTh = 2;

    const unsigned int histogramSize = histogram->Size();

    double totalSum = histogram->GetTotalFrequency();
    unsigned int currentValue = 0;
    unsigned int sum = 0;

    while ( currentValue < histogramSize && double(sum)/double(totalSum) < (lowerTh/100))//0.02 )
    {
        sum += histogram->GetFrequency( currentValue, 0 );
        currentValue++;
    }

    return currentValue;
}

unsigned int BrainExtraction::upperHist(const HistogramType * histogram){
    const float upperTh = 98;

    const unsigned int histogramSize = histogram->Size();
    double totalSum = histogram->GetTotalFrequency();
    unsigned int currentValue = 0;
    unsigned int sum = 0;

    while ( currentValue < histogramSize && double(sum)/double(totalSum) < (upperTh/100))//0.98 )
    {
        sum += histogram->GetFrequency( currentValue, 0 );
        currentValue++;
    }

    return currentValue;
}

MaskImage BrainExtraction::getRoi(Image inputImage, float t2, float t98){
    float t = 0.1*(t98-t2) + t2;
    //std::cout << t << std::endl;

    ThresholdingFilterType::Pointer thresholder = ThresholdingFilterType::New();
    thresholder->SetLowerThreshold(t);
    thresholder->SetInput(inputImage);
    thresholder->SetInsideValue(1);
    thresholder->SetOutsideValue(0);
    thresholder->Update();


    //Filling image holes
    /*NotImageFilterType::Pointer filterInvers1 = NotImageFilterType::New();
    filterInvers1->SetInput(thresholder->GetOutput());
    filterInvers1->Update();

    MaskImageType::Pointer background = getImageLergestConnectedComponent(filterInvers1->GetOutput());

    NotImageFilterType::Pointer filterInvers2 = NotImageFilterType::New();
    filterInvers2->SetInput(background);
    filterInvers2->Update();

    return filterInvers2->GetOutput();*/

    //Filling image holes by dilation and erosion
    //BrainIO* brainio = new BrainIO();
    //brainio->WriteMaskImage("maskROI_1.nii", thresholder->GetOutput());
    MaskImage outputImage = dilateIm(thresholder->GetOutput(), 5);
    outputImage = erodeIm(outputImage, 5);

    //brainio->WriteMaskImage("maskROI_2.nii", outputImage);

    //MaskImage outputImage = thresholder->GetOutput();
    //return thresholder->GetOutput();
    return outputImage;


    /*HoleFillingImageFilterType::Pointer filler = HoleFillingImageFilterType::New();
    ImageType::SizeType radius;
     radius[0] = holeSize;
     radius[1] = holeSize;
     radius[2] = 1;

    filler->SetInput(thresholder->GetOutput());
    filler->SetForegroundValue(1);
    filler->SetBackgroundValue(0);
    filler->SetRadius(radius);
    filler->Update();


    MaskImageType::Pointer outputImage = filler->GetOutput();
    //MaskImageType::Pointer outputImage = thresholder->GetOutput();
    thresholder = NULL;
    //radius = NULL;
    //filler = NULL;

    return outputImage;*/
}

MaskSlice BrainExtraction::otsuThreshold(Slice slice){
    OtsuThresholdImageFilterType::Pointer otsuFilter = OtsuThresholdImageFilterType::New();
    otsuFilter->SetInput(slice);

    otsuFilter->SetInsideValue(0);
    otsuFilter->SetOutsideValue(1);
    otsuFilter->Update();

    MaskSliceType::Pointer outputSlice = otsuFilter->GetOutput();
    otsuFilter = NULL;

    return outputSlice;
}

MaskSlice BrainExtraction::imThreshold(Slice slice, double threshold){

    ThresholdingSliceFilterType::Pointer thresholder = ThresholdingSliceFilterType::New();
    thresholder->SetLowerThreshold(threshold);
    thresholder->SetInput(slice);
    thresholder->SetInsideValue(1);
    thresholder->SetOutsideValue(0);
    thresholder->Update();

    MaskSliceType::Pointer outputSlice = thresholder->GetOutput();
    thresholder = NULL;

    return outputSlice;
}

MaskSlice BrainExtraction::fillSlice(MaskSlice slice){

    NotSliceFilterType::Pointer filterInvers1 = NotSliceFilterType::New();
    filterInvers1->SetInput(slice);
    filterInvers1->Update();

    MaskSliceType::Pointer background = getLergestConnectedComponent(filterInvers1->GetOutput());

    NotSliceFilterType::Pointer filterInvers2 = NotSliceFilterType::New();
    filterInvers2->SetInput(background);
    filterInvers2->Update();

    return filterInvers2->GetOutput();
}



Image BrainExtraction::maskImage(Image inputImage, MaskImage roi){
    MaskImageFilterType::Pointer maskImageFilter = MaskImageFilterType::New();
    maskImageFilter->SetInput( inputImage );
    maskImageFilter->SetOutsideValue( 0 );
    maskImageFilter->SetInput2( roi );
    maskImageFilter->Update();

    ImageType::Pointer outputImage = maskImageFilter->GetOutput();
    maskImageFilter = NULL;

    return outputImage;
}

Slice BrainExtraction::maskSlice(Slice inputSlice, MaskSlice M){
    MaskSliceFilterType::Pointer maskSliceFilter = MaskSliceFilterType::New();
    maskSliceFilter->SetInput( inputSlice );
    maskSliceFilter->SetOutsideValue( -1 );
    maskSliceFilter->SetInput2( M );
    maskSliceFilter->Update();

    SliceType::Pointer outputSlice = maskSliceFilter->GetOutput();
    maskSliceFilter = NULL;

    return outputSlice;
}

MaskSlice BrainExtraction::maskMaskSlice(MaskSlice inputSlice, MaskSlice M){
    MaskSliceMaskFilterType::Pointer maskSliceFilter = MaskSliceMaskFilterType::New();
    maskSliceFilter->SetInput( inputSlice );
    maskSliceFilter->SetOutsideValue( 0 );
    maskSliceFilter->SetInput2( M );
    maskSliceFilter->Update();

    MaskSliceType::Pointer outputSlice = maskSliceFilter->GetOutput();
    maskSliceFilter = NULL;

    return outputSlice;
}

double BrainExtraction::meanImage(Image imROI){
    double meanU = 0;
    int nMaskedVoxels = 0;

    Iterator itROI = Iterator( imROI, imROI->GetRequestedRegion() );

    for (itROI.GoToBegin(); !itROI.IsAtEnd(); ++itROI ) {
        if (itROI.Get()>0) {
            meanU += itROI.Get();
            nMaskedVoxels++;
        }
    }

    return meanU/nMaskedVoxels;
}

double BrainExtraction::getMeanSlice(Slice slice){
    double meanS = 0;
    int nMaskedVoxels = 0;

    SliceIterator itSlice = SliceIterator(slice, slice->GetRequestedRegion());

    for (itSlice.GoToBegin(); !itSlice.IsAtEnd(); ++itSlice ) {
        if (itSlice.Get()>=0) {
            meanS += itSlice.Get();
            nMaskedVoxels++;
        }
    }

    //To control the div 0

    return meanS/nMaskedVoxels;
}

Slice BrainExtraction::getSlice(Image inputImage, int slice){

    SlicerType::Pointer slicer = SlicerType::New();
    RegionType segRegion = inputImage->GetLargestPossibleRegion();
    RegionType desiredRegion;
    SizeType size = segRegion.GetSize();
    size[2] = 0;
    IndexType start = segRegion.GetIndex();
    start[2] = slice;

    desiredRegion.SetSize(  size  );
    desiredRegion.SetIndex( start );

    slicer->SetExtractionRegion( desiredRegion );
    slicer->SetInput( inputImage );
    slicer->SetDirectionCollapseToIdentity();
    slicer->Update();

    SliceType::Pointer img = slicer->GetOutput();
    //img->Update();

    slicer = NULL;

    return img;
}

MaskSlice BrainExtraction::getM(Slice slice){

    MaskSliceType::Pointer M = otsuThreshold(slice);
    M = erodeSlice(M, 2);
    M = getLergestConnectedComponent(M);

    M = dilateSlice(M, 2);

    return M;
}

double BrainExtraction::getTlower(float t_lower, Slice slice, float meanROI){

    MinimumMaximumSliceCalculatorType::Pointer maxIntensityCalculator = MinimumMaximumSliceCalculatorType::New();
    maxIntensityCalculator->SetImage( slice );
    maxIntensityCalculator->Compute();

    SliceHistogramGeneratorType::Pointer histogramGenerator = SliceHistogramGeneratorType::New();
    histogramGenerator->SetInput(slice);
    histogramGenerator->SetNumberOfBins( maxIntensityCalculator->GetMaximum());
    histogramGenerator->SetHistogramMin( 0 );
    histogramGenerator->SetHistogramMax( maxIntensityCalculator->GetMaximum() );
    histogramGenerator->Compute();

    const SliceHistogramType * histogram = histogramGenerator->GetOutput();
    const unsigned int histogramSize = histogram->Size();

    double totalSum = histogram->GetTotalFrequency()-histogram->GetFrequency( 0, 0 );

    //std::cout << histogram->GetTotalFrequency() << "-" << totalSum << " - " << histogramSize << std::endl;

    unsigned int currentValue = 1;
    unsigned int sum = 0;
    while ( currentValue < histogramSize && sum < t_lower*totalSum)
    {
        sum += histogram->GetFrequency( currentValue, 0 );
        currentValue++;
    }

    double meanSlice = getMeanSlice(slice);

    maxIntensityCalculator = NULL;
    histogramGenerator = NULL;

    return currentValue + (meanSlice - meanROI);
}

double BrainExtraction::getTupper(float t_upper, Slice img_inter){

    MinimumMaximumSliceCalculatorType::Pointer maxIntensityCalculator = MinimumMaximumSliceCalculatorType::New();
    maxIntensityCalculator->SetImage( img_inter );
    maxIntensityCalculator->Compute();

    SliceHistogramGeneratorType::Pointer histogramGenerator = SliceHistogramGeneratorType::New();
    histogramGenerator->SetInput(img_inter);
    histogramGenerator->SetNumberOfBins( maxIntensityCalculator->GetMaximum());
    histogramGenerator->SetHistogramMin( 0 );
    histogramGenerator->SetHistogramMax( maxIntensityCalculator->GetMaximum() );
    histogramGenerator->Compute();
    const SliceHistogramType * histogram = histogramGenerator->GetOutput();
    const unsigned int histogramSize = histogram->Size();

    double totalSum = histogram->GetTotalFrequency()-histogram->GetFrequency( 0, 0 );
    unsigned int currentValue = 1;
    unsigned int sum = 0;
    while ( currentValue < histogramSize && sum < t_upper*totalSum)
    {
        sum += histogram->GetFrequency( currentValue, 0 );
        currentValue++;
    }

    maxIntensityCalculator = NULL;
    histogramGenerator = NULL;

    return currentValue;
}

MaskSlice BrainExtraction::erodeSlice(MaskSlice inputSl, unsigned int radius){
    StructuringElementType structuringElement;
    structuringElement.SetRadius(radius);
    structuringElement.CreateStructuringElement();

    GrayscaleErodeSliceFilterType::Pointer erodeFilter = GrayscaleErodeSliceFilterType::New();
    erodeFilter->SetInput(inputSl);
    erodeFilter->SetKernel(structuringElement);
    erodeFilter->Update();

    MaskSliceType::Pointer slice = erodeFilter->GetOutput();

    erodeFilter = NULL;

    return slice;
}

MaskImage BrainExtraction::erodeIm(MaskImage inputIm, unsigned int radius){
    Structuring3DElementType structuringElement;
    structuringElement.SetRadius(radius);
    structuringElement.CreateStructuringElement();

    GrayscaleErodeImageFilterType::Pointer erodeFilter = GrayscaleErodeImageFilterType::New();
    erodeFilter->SetInput(inputIm);
    erodeFilter->SetKernel(structuringElement);
    erodeFilter->Update();

    MaskImageType::Pointer image = erodeFilter->GetOutput();

    erodeFilter = NULL;

    return image;
}

MaskSlice BrainExtraction::getLergestConnectedComponent(MaskSlice maskSl){
    ConnectedComponentImageFilterType::Pointer labelFilter = ConnectedComponentImageFilterType::New ();
    labelFilter->SetInput(maskSl);
    labelFilter->Update();

    LabelShapeKeepNObjectsImageFilterType::Pointer labelShapeKeepNObjectsImageFilter = LabelShapeKeepNObjectsImageFilterType::New();
    labelShapeKeepNObjectsImageFilter->SetInput( labelFilter->GetOutput() );
    labelShapeKeepNObjectsImageFilter->SetBackgroundValue( 0 );
    labelShapeKeepNObjectsImageFilter->SetNumberOfObjects( 1 );
    labelShapeKeepNObjectsImageFilter->SetAttribute( LabelShapeKeepNObjectsImageFilterType::LabelObjectType::NUMBER_OF_PIXELS);
    labelShapeKeepNObjectsImageFilter->Update();

    RescaleFilterType::Pointer rescaleFilter = RescaleFilterType::New();
    rescaleFilter->SetOutputMinimum(0);
    rescaleFilter->SetOutputMaximum(255);
    rescaleFilter->SetInput(labelShapeKeepNObjectsImageFilter->GetOutput());
    rescaleFilter->Update();

    MaskSliceType::Pointer labelSlice = rescaleFilter->GetOutput();

    labelFilter = NULL;
    labelShapeKeepNObjectsImageFilter = NULL;
    rescaleFilter = NULL;

    return labelSlice;
}

MaskSlice BrainExtraction::dilateSlice(MaskSlice inputSl, unsigned int radius){
    StructuringElementType structuringElement;
    structuringElement.SetRadius(radius);
    structuringElement.CreateStructuringElement();

    GrayscaleDilateSliceFilterType::Pointer dilateFilter = GrayscaleDilateSliceFilterType::New();
    dilateFilter->SetInput(inputSl);
    dilateFilter->SetKernel(structuringElement);
    dilateFilter->Update();

    MaskSliceType::Pointer slice = dilateFilter->GetOutput();

    dilateFilter = NULL;

    return slice;
}

MaskImage BrainExtraction::dilateIm(MaskImage inputIm, unsigned int radius){
    Structuring3DElementType structuringElement;
    structuringElement.SetRadius(radius);
    structuringElement.CreateStructuringElement();

    GrayscaleDilateImageFilterType::Pointer dilateFilter = GrayscaleDilateImageFilterType::New();
    dilateFilter->SetInput(inputIm);
    dilateFilter->SetKernel(structuringElement);
    dilateFilter->Update();

    MaskImageType::Pointer image = dilateFilter->GetOutput();

    dilateFilter = NULL;

    return image;
}

MaskSlice BrainExtraction::seedRegion(MaskSlice thrISlice, MaskSlice M){
    //Get connected components
    ConnectedComponentImageFilterType::Pointer labelFilter = ConnectedComponentImageFilterType::New ();
    labelFilter->SetInput(thrISlice);
    labelFilter->Update();

    //Processing
    int * labelCounts = new int[labelFilter->GetObjectCount()];
    int * labelAssignments = new int[labelFilter->GetObjectCount()];
    std::fill_n(labelCounts, labelFilter->GetObjectCount(), 0);
    std::fill_n(labelAssignments, labelFilter->GetObjectCount(), 0);

    MaskSliceType::Pointer connected = labelFilter->GetOutput();

    MaskSliceIterator itComponents = MaskSliceIterator( connected, connected->GetRequestedRegion() );
    MaskSliceIterator itMask = MaskSliceIterator( M, M->GetRequestedRegion() );
    for (itComponents.GoToBegin(), itMask.GoToBegin(); !itComponents.IsAtEnd(); ++itComponents, ++itMask ) {
        if(itComponents.Get()>0){
            int label = itComponents.Get();
            labelCounts[label-1]=labelCounts[label-1]+1;
            if(itMask.Get() > 0)
                labelAssignments[label-1]=1;
        }
    }

   //Assignment
   MaskSliceIterator itComponentAssignment = MaskSliceIterator( connected, connected->GetRequestedRegion() );
   for (itComponentAssignment.GoToBegin(); !itComponentAssignment.IsAtEnd(); ++itComponentAssignment ) {
        if(itComponentAssignment.Get()>0){
            int label = itComponentAssignment.Get();
            if(labelCounts[label-1]<=5){
                itComponentAssignment.Set(2);
            }else
                if(labelAssignments[label-1]==1){
                    itComponentAssignment.Set(3);
                }else{
                    itComponentAssignment.Set(1);                    
                }
        }
    }

    labelFilter = NULL;

    return connected;
}

MaskSlice BrainExtraction::voxelAssignment(MaskSlice assignedVoxels, MaskSlice newVoxels, MaskSlice M){
    /*
     *
     *
     *
     */
    MaskSliceType::SizeType radius;
    radius[0] = 1;
    radius[1] = 1;

    MaskSliceNeighborhoodIterator nbIt = MaskSliceNeighborhoodIterator(radius, assignedVoxels, assignedVoxels->GetRequestedRegion());

    MaskSliceIterator itNewVoxels = MaskSliceIterator( newVoxels, newVoxels->GetRequestedRegion() );
    MaskSliceIterator itAssigned = MaskSliceIterator( assignedVoxels, assignedVoxels->GetRequestedRegion() );
    MaskSliceIterator itM = MaskSliceIterator( M, M->GetRequestedRegion() );
    for (itNewVoxels.Begin(), nbIt.Begin(), itAssigned.Begin(), itM.Begin(); !itNewVoxels.IsAtEnd(); ++itNewVoxels, ++nbIt, ++itAssigned, ++itM) {
        if(itNewVoxels.Get()>0){
            if(nbIt.GetPixel(1)==1 || nbIt.GetPixel(3)==1 || nbIt.GetPixel(5)==1 || nbIt.GetPixel(7)==1)
                itAssigned.Set(1);
            else if((nbIt.GetPixel(1)==3 || nbIt.GetPixel(3)==3 || nbIt.GetPixel(5)==3 || nbIt.GetPixel(7)==3))// && itM.Get()>0)
                itAssigned.Set(3);
            else
                itAssigned.Set(2);
        }
    }

    return assignedVoxels;
}

MaskSlice BrainExtraction::componentAssignment(MaskSlice assignedVoxels, MaskSlice M){

    //Extract undetermined components
    MaskSlice undComp = extractPixelsWithValue(assignedVoxels, 2);

    //Connected components once it is dilated to check and know wich is the component we are checking each time
    ConnectedComponentImageFilterType::Pointer labelFilter = ConnectedComponentImageFilterType::New ();
    labelFilter->SetInput(undComp);
    labelFilter->Update();
    MaskSlice indexComponent = labelFilter->GetOutput();
    int * labelAssignments = new int[labelFilter->GetObjectCount()];
    std::fill_n(labelAssignments, labelFilter->GetObjectCount(), 2);

    //Get boundaries of undetermined components via a rombus dilation and update labelAssignments
    MaskSliceType::SizeType radius;
    radius[0] = 1;
    radius[1] = 1;

    MaskSliceNeighborhoodIterator nbIt = MaskSliceNeighborhoodIterator(radius, assignedVoxels, assignedVoxels->GetRequestedRegion());
    MaskSliceIterator itUndComp = MaskSliceIterator( undComp, undComp->GetRequestedRegion() );
    MaskSliceIterator itIndComp = MaskSliceIterator( indexComponent, indexComponent->GetRequestedRegion() );
    for (itUndComp.Begin(), itIndComp.Begin(), nbIt.Begin(); !itUndComp.IsAtEnd(); ++itUndComp, ++itIndComp, ++nbIt) {
        if(itUndComp.Get()>0){
            if(nbIt.GetPixel(1)==1 || nbIt.GetPixel(3)==1 || nbIt.GetPixel(5)==1 || nbIt.GetPixel(7)==1)
                labelAssignments[itIndComp.Get()-1]=1;
            else if(nbIt.GetPixel(1)==3 || nbIt.GetPixel(3)==3 || nbIt.GetPixel(5)==3 || nbIt.GetPixel(7)==3){
                if(labelAssignments[itIndComp.Get()-1]!=1){
                    labelAssignments[itIndComp.Get()-1]=3;
                }
            }

        }
    }

    //Update assignedVoxels by iterating labelAssignments
    MaskSliceIterator itAssigned = MaskSliceIterator( assignedVoxels, assignedVoxels->GetRequestedRegion() );
    MaskSliceIterator itM = MaskSliceIterator( M, M->GetRequestedRegion() );
    for (itAssigned.Begin(), itIndComp.GoToBegin(), itM.GoToBegin(); !itAssigned.IsAtEnd(); ++itAssigned, ++itIndComp, ++itM){
        if (itIndComp.Get() > 0 )//&& itM.Get()>0){
            itAssigned.Set(labelAssignments[itIndComp.Get()-1]);        
    }


    labelFilter = NULL;

    return assignedVoxels;
}

MaskSlice BrainExtraction::imageSubtract(MaskSlice image1, MaskSlice image2){

    SubtractImageFilterType::Pointer subtractFilter = SubtractImageFilterType::New ();
    subtractFilter->SetInput1(image1);
    subtractFilter->SetInput2(image2);
    subtractFilter->Update();

    MaskSlice imSubstracted = subtractFilter->GetOutput();

    subtractFilter = NULL;

    return imSubstracted;
}

MaskSlice BrainExtraction::extractPixelsWithValue(MaskSlice slice, int pValue){
    ThresholdingMaskSliceFilterType::Pointer thresholder = ThresholdingMaskSliceFilterType::New();
    thresholder->SetUpperThreshold(pValue);
    thresholder->SetLowerThreshold(pValue);
    thresholder->SetInput(slice);
    thresholder->SetInsideValue(1);
    thresholder->SetOutsideValue(0);
    thresholder->Update();

    MaskSliceType::Pointer outputSlice = thresholder->GetOutput();
    thresholder = NULL;

    return outputSlice;
}

MaskImage BrainExtraction::setSlice(MaskImage image, MaskSlice slice, int nSlice){
    RegionType region = image->GetLargestPossibleRegion();
    IndexType start = region.GetIndex();
    start[2] = nSlice;
    SizeType size = region.GetSize();
    size[2] = 0;
    region.SetIndex(start);
    region.SetSize(size);

    MaskSliceIterator itSlice = MaskSliceIterator( slice, slice->GetRequestedRegion() );
    MaskImageIterator itImage = MaskImageIterator( image, region );
    for (itSlice.Begin(), itImage.GoToBegin(); !itSlice.IsAtEnd(); ++itSlice, ++itImage) {
        itImage.Set(itSlice.Get());
    }

    return image;
}

MaskSlice BrainExtraction::getSlice(MaskImage inputImage, int slice){

    MaskSlicerType::Pointer slicer = MaskSlicerType::New();
    RegionType segRegion = inputImage->GetLargestPossibleRegion();
    RegionType desiredRegion;
    SizeType size = segRegion.GetSize();
    size[2] = 0;
    IndexType start = segRegion.GetIndex();
    start[2] = slice;

    desiredRegion.SetSize(  size  );
    desiredRegion.SetIndex( start );

    slicer->SetExtractionRegion( desiredRegion );
    slicer->SetInput( inputImage );
    slicer->SetDirectionCollapseToIdentity();
    slicer->Update();

    MaskSliceType::Pointer img = slicer->GetOutput();
    //img->Update();

    slicer = NULL;

    return img;
}

MaskImage BrainExtraction::removeBrainStem(MaskImage image){
    int * labelAssignments = new int[image->GetLargestPossibleRegion().GetSize()[2]];
    std::fill_n(labelAssignments, image->GetLargestPossibleRegion().GetSize()[2]/2, 0);

    RegionType region = image->GetLargestPossibleRegion();

    MaskImageIterator itImage = MaskImageIterator( image, region );
    for (itImage.GoToBegin(); itImage.GetIndex()[2]<=image->GetLargestPossibleRegion().GetSize()[2]/2; ++itImage){
        labelAssignments[itImage.GetIndex()[2]]+=itImage.Get();
    }

    int min = 0;
    int i;
    for (i=1; i<image->GetLargestPossibleRegion().GetSize()[2]/2; i++){
        if(labelAssignments[i] < labelAssignments[min])
            min = i;
    }

    for (itImage.GoToBegin(); !itImage.IsAtEnd(); ++itImage) {
        if(itImage.GetIndex()[2]<=min)
            itImage.Set(0);
    }

    return image;
}
