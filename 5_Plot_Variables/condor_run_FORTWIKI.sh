B#export SCRAM_ARCH=slc6_amd64_gcc472
export SCRAM_ARCH=slc5_amd64_gcc462
#scram arch for slc6 is defined above
s0;95;cource /osg/app/cmssoft/cms/cmsset_default.sh
#source /apps/02/cmssoft/cms/cmsset_default.sh
#source /osg/app/glite/etc/profile.d/grid_env.sh

cd /net/hisrv0001/home/ilaflott/Leos_Analysis/CMSSW_5_3_20_FOREST_PLOTS/

#eval `scram list CMSSW`
eval `scramv1 runtime -sh`

cd /net/hisrv0001/home/ilaflott/Leos_Analysis/CMSSW_5_3_20_FOREST_PLOTS/src/For_Ian/5_Plot_Variables/

echo "root directory: $ROOTSYS"

gcc --version

echo "Processing..."

root -b -q "bTagPlots_IanEdit.C+(\"condor_test\")"

echo "Done!"

echo "Copying output files to " $destination