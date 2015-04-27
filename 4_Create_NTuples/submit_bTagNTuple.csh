#tell CERN batch where your proxy is if you need to (this one needs to)
export X509_USER_PROXY=/afs/cern.ch/user/u/ilaflott/x509up_u49069
#tell CERN batch to enter project directory and run cmsenv equivalent
cd /afs/cern.ch/user/i/ilaflott/Leos_Analysis/CMSSW_5_3_20_FOREST_PLOTS/src
eval `scramv1 runtime -sh`
#go to relvant directory
cd /afs/cern.ch/user/i/ilaflott/Leos_Analysis/CMSSW_5_3_20_FOREST_PLOTS/src/For_Ian/4_Create_NTuples
#run scripts
root -b -q "bTagNTuple_Original.C+(0)"
echo "Done!"
