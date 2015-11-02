rm BJets*.txt
#OfficialLowPt
ls /mnt/hadoop/cms/store/user/ilaflott/Leos_Analysis/runForest/BJets_pthat15_7.13.15/Pythia_BJet_Pt15_TuneZ2_2760GeV/crab_pthat15/150720_205250/0000/*.root >> BJets_OfficialLowPt_Forests.txt
ls /mnt/hadoop/cms/store/user/ilaflott/Leos_Analysis/runForest/BJets_pthat30_7.13.15/Pythia_BJet_Pt30_TuneZ2_2760GeV/crab_pthat30/150713_214353/0000/*.root >> BJets_OfficialLowPt_Forests.txt
ls /mnt/hadoop/cms/store/user/ilaflott/Leos_Analysis/runForest/BJets_pthat50_7.13.15/Pythia_BJet_Pt50_TuneZ2_2760GeV/crab_pthat50/150713_214536/0000/*.root >> BJets_OfficialLowPt_Forests.txt
ls /mnt/hadoop/cms/store/user/ilaflott/Leos_Analysis/runForest/BJets_pthat80_7.13.15/Pythia_BJet_Pt80_TuneZ2_2760GeV/crab_pthat80/150713_214559/0000/*.root >> BJets_OfficialLowPt_Forests.txt
ls /mnt/hadoop/cms/store/user/ilaflott/Leos_Analysis/runForest/BJets_pthat120_7.13.15/Pythia_BJet_Pt120_TuneZ2_2760GeV/crab_pthat120/150713_214620/0000/*.root >> BJets_OfficialLowPt_Forests.txt
ls /mnt/hadoop/cms/store/user/ilaflott/Leos_Analysis/runForest/BJets_pthat170_7.13.15/Pythia_BJet_Pt170_TuneZ2_2760GeV/crab_pthat170_pt1/150713_214645/0000/*.root >> BJets_OfficialLowPt_Forests.txt
#unOfficialHighPt
ls /mnt/hadoop/cms/store/user/ilaflott/Leos_Analysis/runForest/BJets_pthat170_7.13.15/ppBJetMC_pthat280_MultiCrab_2760GeV_STEP1_GEN/crab_pthat170_pt2/150715_204631/0000/*.root >> BJets_unOfficialHighPt_Forests.txt
ls /mnt/hadoop/cms/store/user/ilaflott/Leos_Analysis/runForest/BJets_pthat220_7.13.15/ppBJetMC_pthat220_MultiCrab_2760GeV_STEP1_GEN/crab_pthat220/150715_165457/0000/*.root >> BJets_unOfficialHighPt_Forests.txt
ls /mnt/hadoop/cms/store/user/ilaflott/Leos_Analysis/runForest/BJets_pthat370_7.13.15/ppBJetMC_pthat370_MultiCrab_2760GeV_STEP1_GEN/crab_pthat370/150714_173857/0000/*.root >> BJets_unOfficialHighPt_Forests.txt
ls /mnt/hadoop/cms/store/user/ilaflott/Leos_Analysis/runForest/BJets_pthat460_7.13.15/ppBJetMC_pthat460_MultiCrab_2760GeV_STEP1_GEN/crab_pthat460/150714_173911/0000/*.root >> BJets_unOfficialHighPt_Forests.txt
ls /mnt/hadoop/cms/store/user/ilaflott/Leos_Analysis/runForest/BJets_pthat540_7.13.15/ppBJetMC_pthat540_MultiCrab_2760GeV_STEP1_GEN/crab_pthat540/150714_173924/0000/*.root >> BJets_unOfficialHighPt_Forests.txt
#unOfficialAddStat
ls /mnt/hadoop/cms/store/user/ilaflott/Leos_Analysis/runForest_AddStat/pthat80/BJetMC_pthat80_STEP1_GEN/crab_pthat80/151010_213401/000?/*.root >> BJets_unOfficialHighPt_addStat_Forests.txt
ls /mnt/hadoop/cms/store/user/ilaflott/Leos_Analysis/runForest_AddStat/pthat120/BJetMC_pthat120_STEP1_GEN/crab_pthat120/151010_213427/000?/*.root >> BJets_unOfficialHighPt_addStat_Forests.txt
ls /mnt/hadoop/cms/store/user/ilaflott/Leos_Analysis/runForest_AddStat/pthat170/BJetMC_pthat170_STEP1_GEN/crab_pthat170/151010_213452/000?/*.root >> BJets_unOfficialHighPt_addStat_Forests.txt
ls /mnt/hadoop/cms/store/user/ilaflott/Leos_Analysis/runForest_AddStat/pthat220/BJetMC_pthat220_STEP1_GEN/crab_pthat220/151010_213510/000?/*.root >> BJets_unOfficialHighPt_addStat_Forests.txt
ls /mnt/hadoop/cms/store/user/ilaflott/Leos_Analysis/runForest_AddStat/pthat280/BJetMC_pthat280_STEP1_GEN/crab_pthat280/151010_213526/000?/*.root >> BJets_unOfficialHighPt_addStat_Forests.txt
ls /mnt/hadoop/cms/store/user/ilaflott/Leos_Analysis/runForest_AddStat/pthat370/BJetMC_pthat370_STEP1_GEN/crab_pthat370/151010_213542/000?/*.root >> BJets_unOfficialHighPt_addStat_Forests.txt
ls /mnt/hadoop/cms/store/user/ilaflott/Leos_Analysis/runForest_AddStat/pthat460/BJetMC_pthat460_STEP1_GEN/crab_pthat460/151010_213558/000?/*.root >> BJets_unOfficialHighPt_addStat_Forests.txt
ls /mnt/hadoop/cms/store/user/ilaflott/Leos_Analysis/runForest_AddStat/pthat540/BJetMC_pthat540_STEP1_GEN/crab_pthat540/151010_213615/000?/*.root >> BJets_unOfficialHighPt_addStat_Forests.txt
cat BJets*.txt >> BJets_allAvailable_Forests.txt