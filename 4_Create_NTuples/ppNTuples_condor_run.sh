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
job=$1
flavor=$2

#mkdir /mnt/hadoop/cms/store/user/ilaflott/Leos_Analysis_NTuples
#chmod a+rxw /mnt/hadoop/cms/store/user/ilaflott/Leos_Analysis_NTuples
#mkdir /mnt/hadoop/cms/store/user/ilaflott/Leos_Analysis_weight_info
#chmod a+rxw /mnt/hadoop/cms/store/user/ilaflott/Leos_Analysis_weight_info
#rm /mnt/hadoop/cms/store/user/ilaflott/Leos_Analysis_NTuples/*.root

#execute the script like this...
root -b -l <<EOF
.x bTagNTuple.C+(${job}, ${flavor})
.q
EOF

mv *.root /mnt/hadoop/cms/store/user/ilaflott/Leos_Analysis_NTuples
mv *.txt  /mnt/hadoop/cms/store/user/ilaflott/Leos_Analysis_weight_info
