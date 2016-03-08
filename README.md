# MARGA (Multispectral Adaptive Region Growing Algorithm)

## About: 

**MARGA** (Multispectral Adaptive Region Growing Algorithm) is a skull stripping algorithm that uses the complementary information provided by conventional MRI such as T1-weighted and T2-weighted to perform the brain segmentation in axial views and low quality images. The seed regions are spread using a 2D RG algorithm which behaves differently in specific zones of the brain, allowing to deal with the fact that middle MRI slices have better image contrast between the brain and non-brain regions than superior and inferior brain slices where the contrast is smaller.

Please browse the SALEM webpage for more information (http://eia.udg.edu/salem/margaToolbox).


## Installation:

MARGA has been developed using the Insight ToolKit. In order to compile some of the code included, you should have a working itk installation (detailed installation instructions for windows and linux can be found at http://www.itk.org/CourseWare/Training/GettingStarted-I.pdf).

It can be compiled using the included cmakelist.txt file with the cmake cross platform compilation tool.

Known issues:

 *  The code provided has been compiled with itk 4.x .
	
 *  For some reason, the cxx compiler g++ is not installed by default in some ubuntu distributions. Fix this with : sudo apt-get install g++
	
 *  Also itk needs the uuid tool but does not say so, you can install it from the ubuntu software center , then you will also need the development files so you also have to run: sudo apt-get install uuid-dev
 
 *  When configuring itk with cmake (after having followed the guidelines), you need to set two flags (press "t" to see the flags and set ITK_USE_OPTIMIZED_REGISTRATION and ITK_USE_REVIEW to ON.


To re-compile this code run the following commands:

```
$ cd MARGA/Build
$ ccmake ../
$ make -j 4 //depending on your number of cores
$ make install -j 4 //depending on your number of cores
```

## Usage:

Data formats: All tests have been run using the .nii (NIFTI) format. If you need to perform conversions from nifti, you might want to use the "dcm2nii" tool (can be found in the ubuntu software centre).

Running the shell script Process.sh:

```
input source path
input T1 image
input T2 image
output image

./Process ./SourcePath ./sourceT1_image ./sourceT2_image ./output_image"

*Attention: All the images must be entered without the extension!

```


However, MARGA can be applied directly to a single image with different parameters to play around:

```
Marga -input <filename> [OPTIONS] 

	OPTIONS:
		-output <filename> -> mask image filename (default: MARGA_mask.nii)
		-t_lower <float> -> Lower threshold (default: 0.35)
		-t_upper <float> -> Upper threshold (default: 0.80)
		-minSlice <int> ->  Minimum slice where to cut the process (default: 0)
		-maxSlice <int> ->  Maximum slice where to cut the process <numberOfSlices - maxSlice> (default: 0)
		-minSm <int> ->  Minimum slice where to reduce the threshold in 0.03 (default: 0)
		-maxSm <int> ->  Maximum slice where to reduce the threshold in 0.03 <numberOfSlices - maxSm> (default: 0)
		-dilateRadius <int> -> Radius in mm used in the dilation morphological operation (default: 5mm)
		-erodeRadius <int> ->  Radius in mm used in the erosion morphological operation (default: 5mm)

********************************************************************************
	Default command: 
		Marga -input T1.nii -output MARGA_mask.nii -t_lower 0.35 -t_upper 0.8 -minSlice 0 -maxSlice 0 -minSm 0 -maxSm 0 -dilateRadius 5 -erodeRadius 5

```



With this last release, althoug originally it was thought to be applied on T1w and T2w image of 1.5T, it has been tested over 3T FLAIR images obtaining a good performance. 
