root -l -q "bTagNTuple.C+(2,1)"
root -l -q "bTagNTuple.C+(2,2)"
root -l -q "bTagNTuple.C+(2,3)"
hadd -rf good_NTuples/TotalMCNTuple_WithWeights.root good_NTuples/CJets_NTuple_8.22.15_withWeights.root good_NTuples/BJets_NTuple_8.22.15_withWeights.root good_NTuples/QCDJets_NTuple_8.22.15_withWeights.root 