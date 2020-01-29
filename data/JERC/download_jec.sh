#!/bin/bash

# This script downloads JEC tarballs from the official repository, unpacks them, deletes unneeded
# files, and copies remaining ones into the current directory.

set -e


year="$1"

if [ "$year" = "2016" ]; then
  dataVersion="11"
  simVersion="11"
  date="Summer16_07Aug2017"
  periods="BCD EF GH"
elif [ "$year" = "2017" ]; then
  dataVersion="32"
  simVersion="32"
  date="Fall17_17Nov2017"
  periods="B C DE F"
else
  echo "unrecognized year \"${year}\""
  exit 1
fi

for period in $periods; do
  blocks+=("${date}${period}_V${dataVersion}_DATA")
done

blocks+=("${date}_V${simVersion}_MC")

dirSource=`pwd`

# Temporary directory to store intermediate files
dirTemp=`mktemp -d`
cd $dirTemp
echo $dirTemp


# Download and unpack all data
dirUnpacked=unpacked
mkdir $dirUnpacked

for block in ${blocks[*]}
do
    wget https://github.com/cms-jet/JECDatabase/raw/master/tarballs/${block}.tar.gz
    tar -xf ${block}.tar.gz
    
    if [ -d textFiles/${block} ]; then
        mv textFiles/${block}/* $dirUnpacked/
    elif [ -d ${block} ]; then
        mv ${block}/* $dirUnpacked/
    else
        mv ${block}* $dirUnpacked/
    fi
done


# Remove all files that are not needed
cd $dirUnpacked
find . -not -name "*AK4PFchs*" -delete
find . -name "*DataMcSF*" -delete
find . \( -name "*_L1RC_*" -or -name "*_L2Residual_*" \) -delete
find . -name "*_MC_L2L3Residual*" -delete
find . -name "*_DATA_Uncertainty*" -delete
find . -name "*_UncertaintySources_*" -delete

chmod 644 *.txt


# Move remaining files to the original current directory and remove temporary
# directory
mv *.txt $dirSource
cd $dirSource
rm -r $dirTemp

