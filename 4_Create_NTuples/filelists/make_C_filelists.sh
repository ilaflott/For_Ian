rm CJets*.txt
#OfficialLowPt
ls /mnt/hadoop/cms/store/user/ilaflott/Leos_Analysis/runForest/CJets_pthat15_7.13.15/Pythia_CJet_Pt15_TuneZ2_2760GeV/crab_pthat15/150713_214332/0000/*.root >> CJets_OfficialLowPt_Forests.txt
ls /mnt/hadoop/cms/store/user/ilaflott/Leos_Analysis/runForest/CJets_pthat30_7.13.15/Pythia_CJet_Pt30_TuneZ2_2760GeV/crab_pthat30/150713_214349/0000/*.root >> CJets_OfficialLowPt_Forests.txt
ls /mnt/hadoop/cms/store/user/ilaflott/Leos_Analysis/runForest/CJets_pthat50_7.13.15/Pythia_CJet_Pt50_TuneZ2_2760GeV/crab_pthat50/150713_214410/0000/*.root >> CJets_OfficialLowPt_Forests.txt
ls /mnt/hadoop/cms/store/user/ilaflott/Leos_Analysis/runForest/CJets_pthat80_7.13.15/Pythia_CJet_Pt80_TuneZ2_2760GeV/crab_pthat80/150713_214426/0000/*.root >> CJets_OfficialLowPt_Forests.txt
ls /mnt/hadoop/cms/store/user/ilaflott/Leos_Analysis/runForest/CJets_pthat120_7.13.15/Pythia_CJet_Pt120_TuneZ2_2760GeV/crab_pthat120/150713_214443/0000/*.root >> CJets_OfficialLowPt_Forests.txt
ls /mnt/hadoop/cms/store/user/ilaflott/Leos_Analysis/runForest/CJets_pthat170_7.13.15/Pythia_CJet_Pt170_TuneZ2_2760GeV/crab_pthat170/150713_214459/0000/*.root >> CJets_OfficialLowPt_Forests.txt
#unOfficialHighPt
ls /mnt/hadoop/cms/store/user/ilaflott/Leos_Analysis/runForest/CJets_pthat220_7.13.15/ppCJetMC_pthat220_MultiCrab_2760GeV_STEP1_GEN/crab_pthat220/150714_174106/0000/*.root >> CJets_unOfficialHighPt_Forests.txt
ls /mnt/hadoop/cms/store/user/ilaflott/Leos_Analysis/runForest/CJets_pthat280_7.13.15/ppCJetMC_pthat280_MultiCrab_2760GeV_STEP1_GEN/crab_pthat280/150714_174129/0000/*.root >> CJets_unOfficialHighPt_Forests.txt
ls /mnt/hadoop/cms/store/user/ilaflott/Leos_Analysis/runForest/CJets_pthat370_7.13.15/ppCJetMC_pthat370_MultiCrab_2760GeV_STEP1_GEN/crab_pthat370/150714_174144/0000/*.root >> CJets_unOfficialHighPt_Forests.txt
ls /mnt/hadoop/cms/store/user/ilaflott/Leos_Analysis/runForest/CJets_pthat460_7.13.15/ppCJetMC_pthat460_MultiCrab_2760GeV_STEP1_GEN/crab_pthat460/150714_174158/0000/*.root >> CJets_unOfficialHighPt_Forests.txt
ls /mnt/hadoop/cms/store/user/ilaflott/Leos_Analysis/runForest/CJets_pthat540_7.13.15/ppCJetMC_pthat540_MultiCrab_2760GeV_STEP1_GEN/crab_pthat540/150714_174215/0000/*.root >> CJets_unOfficialHighPt_Forests.txt
#unOfficialAddStat
ls    /mnt/hadoop/cms/store/user/ilaflott/Leos_Analysis/runForest_AddStat/pthat80/CJetMC_pthat80_STEP1_GEN/crab_pthat80/151011_170935/000?/*.root >> CJets_unOfficialHighPt_addStat_Forests.txt
ls /mnt/hadoop/cms/store/user/ilaflott/Leos_Analysis/runForest_AddStat/pthat120/CJetMC_pthat120_STEP1_GEN/crab_pthat120/151011_170952/000?/*.root >> CJets_unOfficialHighPt_addStat_Forests.txt
ls /mnt/hadoop/cms/store/user/ilaflott/Leos_Analysis/runForest_AddStat/pthat170/CJetMC_pthat170_STEP1_GEN/crab_pthat170/151011_171008/000?/*.root >> CJets_unOfficialHighPt_addStat_Forests.txt
ls /mnt/hadoop/cms/store/user/ilaflott/Leos_Analysis/runForest_AddStat/pthat220/CJetMC_pthat220_STEP1_GEN/crab_pthat220/151011_171025/000?/*.root >> CJets_unOfficialHighPt_addStat_Forests.txt
ls /mnt/hadoop/cms/store/user/ilaflott/Leos_Analysis/runForest_AddStat/pthat280/CJetMC_pthat280_STEP1_GEN/crab_pthat280/151011_171042/000?/*.root >> CJets_unOfficialHighPt_addStat_Forests.txt
ls /mnt/hadoop/cms/store/user/ilaflott/Leos_Analysis/runForest_AddStat/pthat370/CJetMC_pthat370_STEP1_GEN/crab_pthat370/151011_171058/000?/*.root >> CJets_unOfficialHighPt_addStat_Forests.txt
ls /mnt/hadoop/cms/store/user/ilaflott/Leos_Analysis/runForest_AddStat/pthat460/CJetMC_pthat460_STEP1_GEN/crab_pthat460/151011_171114/000?/*.root >> CJets_unOfficialHighPt_addStat_Forests.txt
ls /mnt/hadoop/cms/store/user/ilaflott/Leos_Analysis/runForest_AddStat/pthat540/CJetMC_pthat540_STEP1_GEN/crab_pthat540/151011_171130/000?/*.root >> CJets_unOfficialHighPt_addStat_Forests.txt
cat CJets*.txt >> CJets_allAvailable_Forests.txt