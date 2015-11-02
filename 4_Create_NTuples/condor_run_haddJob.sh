#!/bin/bash

#untar the tarball fed to condor via the submit script
#echo "untarring tarball"
tar -zxvf run_ppNTuples.tar > /dev/null

#cmsenv
echo "cmsenv"
#export SCRAM_ARCH=slc6_amd64_gcc472
source /osg/app/cmssoft/cms/cmsset_default.sh

#something i'll leave here for now
echo "Job started on `date` at WN: `hostname` "
echo "Job is running on `uname -a`"

#what compiler is being used
gcc --version

#debugging
echo "listing contents.."
ls

#arguments used by macro fed to runscript by submit script
echo "Processing..."

flavor=$1

fileName=""
dir=""

if [ $flavor -eq 0 ]; then
    echo "hadd'ing Data in Data Folder..."
    fileName="data_NTuple_11.1.root"
#    dir="/mnt/hadoop/cms/store/user/ilaflott/Leos_Analysis_NTuples/Data_reForest_11.1"
    dir="/mnt/hadoop/cms/store/user/ilaflott/Leos_Analysis_NTuples/Data_MaybeBadAndOldFiles_11.1"
fi
if [ $flavor -eq 1 ]; then
    echo "hadd'ing QCDJets files in QCDJets Folder..."
    fileName="QCDJets_NTuple_noVsJets_11.1.root"
    dir="/mnt/hadoop/cms/store/user/ilaflott/Leos_Analysis_NTuples/QCDJets_Official_11.1"
fi
if [ $flavor -eq 2 ]; then
    echo "hadd'ing BJets files in BJets Folder..."
    fileName="BJets_NTuple_addStat_11.1.root"
    dir="/mnt/hadoop/cms/store/user/ilaflott/Leos_Analysis_NTuples/BJets_addStat_11.1"
fi
if [ $flavor -eq 3 ]; then
    echo "hadd CJets file in CJets Folder..."
    fileName="CJets_NTuple_addStat_11.1.root"
    dir="/mnt/hadoop/cms/store/user/ilaflott/Leos_Analysis_NTuples/CJets_addStat_11.1"
fi

hadd -f "${fileName}" "${dir}/"*.root
echo "moving hadd'd output to folder"
mkdir "${dir}/haddOutput"
mv "${fileName}" "${dir}/haddOutput/."

