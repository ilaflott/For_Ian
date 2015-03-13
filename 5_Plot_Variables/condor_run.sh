#export SCRAM_ARCH=slc6_amd64_gcc472
export SCRAM_ARCH=slc5_amd64_gcc462
#scram arch for slc6 is defined above
source /osg/app/cmssoft/cms/cmsset_default.sh
#source /apps/02/cmssoft/cms/cmsset_default.sh
#source /osg/app/glite/etc/profile.d/grid_env.sh

cd /net/hisrv0001/home/ilaflott/Leos_Analysis/CMSSW_5_3_20_FOREST_PLOTS/

#eval `scram list CMSSW`
eval `scramv1 runtime -sh`

cd /net/hisrv0001/home/ilaflott/Leos_Analysis/CMSSW_5_3_20_FOREST_PLOTS/src/For_Ian/5_Plot_Variables/

echo "root directory: $ROOTSYS"

gcc --version

echo "Processing..."

root -b -q "bTagPlots_IanEdit.C+(\"vz15_jteta2_jtpt40_HLTPAMu12v1_muCut_WCut_3.12.15\",0,\"vz<15&&vz>-15&&jteta<2&&jteta>-2&&jtpt>40&&HLT_PAMu12_v1&&mupt!=0&&mupt/rawpt<0.95\",1)"
root -b -q "bTagPlots_IanEdit.C+(\"vz15_jteta2_jtpt30_HLTPAMu3v1_muCut_WCut_3.12.15\",0,\"vz<15&&vz>-15&&jteta<2&&jteta>-2&&jtpt>30&&HLT_PAMu3_v1&&mupt!=0&&mupt/rawpt<0.95\",1)"
root -b -q "bTagPlots_IanEdit.C+(\"vz15_jteta2_jtpt30_HLTPAMu7v1_muCut_WCut_3.12.15\",0,\"vz<15&&vz>-15&&jteta<2&&jteta>-2&&jtpt>30&&HLT_PAMu7_v1&&mupt!=0&&mupt/rawpt<0.95\",1)"
root -b -q "bTagPlots_IanEdit.C+(\"vz15_jteta2_jtpt30_HLTPAMu12v1_muCut_WCut_3.12.15\",0,\"vz<15&&vz>-15&&jteta<2&&jteta>-2&&jtpt>30&&HLT_PAMu12_v1&&mupt!=0&&mupt/rawpt<0.95\",1)"
root -b -q "bTagPlots_IanEdit.C+(\"vz15_jteta2_jtpt40_HLTPAMu3PFJet40v1_muCut_WCut_3.12.15\",0,\"vz<15&&vz>-15&&jteta<2&&jteta>-2&&jtpt>40&&HLT_PAMu3PFJet40_v1&&mupt!=0&&mupt/rawpt<0.95\",1)"
root -b -q "bTagPlots_IanEdit.C+(\"vz15_jteta2_jtpt30_HLTPAMu3PFJet40v1_muCut_WCut_3.12.15\",0,\"vz<15&&vz>-15&&jteta<2&&jteta>-2&&jtpt>30&&HLT_PAMu3PFJet40_v1&&mupt!=0&&mupt/rawpt<0.95\",1)"
root -b -q "bTagPlots_IanEdit.C+(\"vz15_jteta2_jtpt30_HLTPAMu3v1_muCut_WCut_svtxCut_3.12.15_oneCutAtATime\",0,\"vz<15&&vz>-15&&jteta<2&&jteta>-2&&jtpt>30&&HLT_PAMu3_v1&&mupt!=0&&mupt/rawpt<0.95&&svtxdl>0.01&&svtxdl<2.5&&svtxdls>3.0&&svtxm<6.5\",1)"
root -b -q "bTagPlots_IanEdit.C+(\"vz15_jteta2_jtpt30_HLTPAMu3v1_muCut_svtxCut_3.12.15_oneCutAtATime\",0,\"vz<15&&vz>-15&&jteta<2&&jteta>-2&&jtpt>30&&HLT_PAMu3_v1&&mupt!=0&&svtxdl>0.01&&svtxdl<2.5&&svtxdls>3.0&&svtxm<6.5\",1)"
root -b -q "bTagPlots_IanEdit.C+(\"vz15_jteta2_jtpt30_HLTPAMu3v1_svtxCut_3.12.15_oneCutAtATime\",0,\"vz<15&&vz>-15&&jteta<2&&jteta>-2&&jtpt>30&&HLT_PAMu3_v1&&svtxdl>0.01&&svtxdl<2.5&&svtxdls>3.0&&svtxm<6.5\",1)"
root -b -q "bTagPlots_IanEdit.C+(\"vz15_jteta2_jtpt30_HLTPAMu3v1_3.12.15_oneCutAtATime\",0,\"vz<15&&vz>-15&&jteta<2&&jteta>-2&&jtpt>30&&HLT_PAMu3_v1\",1)"
root -b -q "bTagPlots_IanEdit.C+(\"vz15_jteta2_jtpt30_3.12.15_oneCutAtATime\",0,\"vz<15&&vz>-15&&jteta<2&&jteta>-2&&jtpt>30\",1)"
root -b -q "bTagPlots_IanEdit.C+(\"jtpt30_3.12.15_oneCutAtATime\",0,\"jtpt>30&&vz<15&&vz>-15\",1)"
root -b -q "bTagPlots_IanEdit.C+(\"noCuts_3.12.15_oneCutAtATime\",0,\"vz<15&&vz>-15\",1)"

echo "Done!"

echo "Copying output files to " $destination