#include <brainio.h>

BrainIO::BrainIO()
{
}

Image BrainIO::ReadImage(FileNameType fileName) {
    // Init
    ReaderType::Pointer intreader = ReaderType::New();;
    ImageType::Pointer intimage = NULL;

    if (fileName.compare("")!=0) {
        intreader->SetFileName( fileName );
        // Reading the image.
        try {
            intreader->Update();
            intimage = intreader->GetOutput();
        }
        catch (itk::ExceptionObject & e) {
            std::cerr << "exception in file reader " << std::endl;
            std::cerr << e.GetDescription() << std::endl;
            std::cerr << e.GetLocation() << std::endl;
        }
    }
    return(intimage);
}

MaskImage BrainIO::ReadMaskImage(FileNameType fileName) {
    // Init
    MaskReaderType::Pointer breader = MaskReaderType::New();
    MaskImageType::Pointer bimage = NULL;

    if (fileName.compare("")!=0) {
        breader->SetFileName( fileName );
        // Reading the image.
        try {
            breader->Update();
            bimage = breader->GetOutput();
        }
        catch (itk::ExceptionObject & e) {
            std::cerr << "exception in file reader " << std::endl;
            std::cerr << e.GetDescription() << std::endl;
            std::cerr << e.GetLocation() << std::endl;
        }
    }
    return(bimage);
}

void BrainIO::WriteImage(FileNameType fileName, Image image) {

    WriterType::Pointer writer = WriterType::New();

    if (fileName.compare("")!=0) {
        writer->SetFileName( fileName );
        writer->SetInput( image );
        // Writting the image.
        try {
            writer->Update();
        }
        catch (itk::ExceptionObject & e) {
            std::cerr << "exception in file writer " << std::endl;
            std::cerr << e.GetDescription() << std::endl;
            std::cerr << e.GetLocation() << std::endl;
        }
    }
}

void BrainIO::WriteSlice(FileNameType fileName, Slice slice) {

    WriterSliceType::Pointer writer = WriterSliceType::New();

    if (fileName.compare("")!=0) {
        writer->SetFileName( fileName );
        writer->SetInput( slice );
        // Writting the slice.
        try {
            writer->Update();
        }
        catch (itk::ExceptionObject & e) {
            std::cerr << "exception in file writer " << std::endl;
            std::cerr << e.GetDescription() << std::endl;
            std::cerr << e.GetLocation() << std::endl;
        }
    }
}

void BrainIO::WriteMaskImage(FileNameType fileName, MaskImage bimage) {

    MaskWriterType::Pointer writer = MaskWriterType::New();

    if (fileName.compare("")!=0) {
        writer->SetFileName( fileName );
        writer->SetInput( bimage );
        // Writting the image.
        try {
            writer->Update();
        }
        catch (itk::ExceptionObject & e) {
            std::cerr << "exception in file writer " << std::endl;
            std::cerr << e.GetDescription() << std::endl;
            std::cerr << e.GetLocation() << std::endl;
        }
    }
}

void BrainIO::WriteMaskSlice(FileNameType fileName, MaskSlice slice) {

    MaskSliceWriterType::Pointer writer = MaskSliceWriterType::New();

    if (fileName.compare("")!=0) {
        writer->SetFileName( fileName );
        writer->SetInput( slice );
        // Writting the slice.
        try {
            writer->Update();
        }
        catch (itk::ExceptionObject & e) {
            std::cerr << "exception in file writer " << std::endl;
            std::cerr << e.GetDescription() << std::endl;
            std::cerr << e.GetLocation() << std::endl;
        }
    }
}




