#!/bin/sh 

# USAGE:
if [ "$#" != "4" ];
then
    echo "*********************************************************************************"
    echo "Usage:"
    echo "		./Process dataPath"
    echo
    echo "Example: ./Process ./SourcePath(/marga ./sourceT1_image ./sourceT2_image ./output_image"
    echo " "
    echo "**Atention**"
    echo "      The only format accepted is .nii!"
    echo "      All the images must be entered without the extension!"
    echo " "
    echo "Current entry: $0 $*"
    echo "*********************************************************************************"
    echo " "

    exit
fi

#Get parameters
sourcePath=$1
subjectT1=$2
subjectT2=$3
subjectout=$4

echo
echo "	**************************************************************"
echo "	* Current entry:"
echo "	*	$0 $*"
echo " "
echo "  **Atention**"
echo "      The only format accepted is .nii!"
echo "      All the images must be entered without the extension!"
echo "	**************************************************************"
echo

timeIn=$(date +%s)

datapath=$1
	
	#########################################################################################
	################  - SKULL STRIPPING -  ###################################################
	#########################################################################################
	echo 'Processing -> ' $subjectT1	
	$sourcePath -input $subjectT1.nii -output $subjectout.tmp.t1.nii -t_lower 0.35 -t_upper 0.8 -minSlice 0 -maxSlice 0 -minSm 0 -maxSm 0 -dilateRadius 5 -erodeRadius 5

	echo 'Processing -> ' $subjectT2
	$sourcePath -input $subjectT2.nii -output $subjectout.tmp.t2.nii -t_lower 0.35 -t_upper 0.8 -minSlice 0 -maxSlice 0 -minSm 0 -maxSm 0 -dilateRadius 5 -erodeRadius 5

	echo 'Combining -> ' $subjectT1 ' and ' $subjectT2
	fslmaths $subjectout.tmp.t1.nii -add $subjectout.tmp.t2.nii $subjectout.tmp.nii
	fslmaths $subjectout.tmp.nii -bin $subjectout.tmp.nii
	fslmaths $subjectout.tmp.nii -kernel sphere 5 -dilM $subjectout.tmp.nii
	fslmaths $subjectout.tmp.nii -kernel sphere 5 -ero $subjectout.tmp.nii
	$sourcePath/fillInHoles $subjectout.tmp.nii $subjectout.nii
	rm $subjectout.tmp.nii.gz


timeOut=$(date +%s)
diffTime=$(($timeOut - $timeIn ))

echo "All data set processed in $diffTime seconds"




