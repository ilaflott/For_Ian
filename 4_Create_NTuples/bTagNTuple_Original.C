/*
  LAST EDITED: 02.11.2014

  DESCRIPTION:
  Produces an ntuple from data crab output or MC pthat files with
  some selections. Run once for each: data, QCD MC, BJet MC, CJet MC

  Input: crab output/MC pthat bin files
  Output: a ntuple/tree

  EVENT SELECTION:
  |hiEvt.vz| <= 15

  DATA:
  skim.pPAcollisionEventSelectionPA
  skim.pHBHENoiseFilter
  trackMax/jtpt > 0.01            << Actual collision event

  JET SELECTION:
  jtpt >= 20
  |jteta| <= 2

  MUON SELECTION:
  (mupt!=0 || mueta!=0 || muphi!=0)
  mupt/rawpt < 0.95                   << W event filtering
*/

// C/C++ library includes
#include <iostream>
#include <fstream>
#include <string>
#include <cmath>

// Root includes
#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TH1D.h"

using namespace std;

// Function declarations
int makeNTuple(int type);
static double MCWeights(double pthat);
static void heavyJetWeights(double *pthatEntries);
static inline void newBranches(TTree *newTree);
static inline void branchAddresses(TTree *akPu3);
//int mergeMCSamples();
//int calculateWeights(int type);
//int forestStatistics(int type);

///////////////////////////////////////////////////////////////////////////////////////////////////
//
//
//
//
//              GLOBAL VARIABLES
//
//
//
//
///////////////////////////////////////////////////////////////////////////////////////////////////

// Switches
//int dataType = 0;
string weights_file;
const bool doTracks=true;

// Macro settings/constants
//FileLists
const string fileListPath = "/net/hisrv0001/home/ilaflott/Leos_Analysis/CMSSW_5_3_20_FOREST_PLOTS/src/For_Ian/4_Create_NTuples/filelists/";
const string dataFileList = "/afs/cern.ch/user/i/ilaflott/Leos_Analysis/CMSSW_5_3_20_FOREST_PLOTS/src/For_Ian/4_Create_NTuples/filelists/ppMuon_data_filelist.txt";
const string QCDFileList  = "/net/hisrv0001/home/ilaflott/Leos_Analysis/CMSSW_5_3_20_FOREST_PLOTS/src/For_Ian/4_Create_NTuples/filelists/pp_MC_QCD_filelist.txt";
const string BJetFileList = "/net/hisrv0001/home/ilaflott/Leos_Analysis/CMSSW_5_3_20_FOREST_PLOTS/src/For_Ian/4_Create_NTuples/filelists/pp_MC_BJet_2760GeV_halfOfficial_filelist.txt";
const string CJetFileList = "/net/hisrv0001/home/ilaflott/Leos_Analysis/CMSSW_5_3_20_FOREST_PLOTS/src/For_Ian/4_Create_NTuples/filelists/pp_MC_CJet_2760GeV_halfOfficial_filelist.txt";

//Weight Files
const string weightFilePath = "/net/hisrv0001/home/ilaflott/pp_MC_2760GeV_bTag_forests_ntuples/NTuples/";
const string dataWeightsFile = "/afs/cern.ch/work/i/ilaflott/bTagNTuples_ppMC_2760GeV/weights_data_updated_4.13.15.txt";			 
const string QCDWeightsFile  = "/net/hisrv0001/home/ilaflott/pp_MC_2760GeV_bTag_forests_ntuples/NTuples/kurts_QCDMC/weights_QCD_updated_4.13.15.txt";		 
const string BJetWeightsFile = "/net/hisrv0001/home/ilaflott/pp_MC_2760GeV_bTag_forests_ntuples/NTuples/halfOfficial_HFMC/weights_BJet_halfOfficial_updated_4.13.15_schemeA.txt";
const string CJetWeightsFile = "/net/hisrv0001/home/ilaflott/pp_MC_2760GeV_bTag_forests_ntuples/NTuples/halfOfficial_HFMC/weights_CJet_halfOfficial_updated_4.13.15_schemeA.txt";

//Output Files
const string outputFilePath = "/net/hisrv0001/home/ilaflott/pp_MC_2760GeV_bTag_forests_ntuples/NTuples";
const char* dataOutFile   = "/afs/cern.ch/work/i/ilaflott/bTagNTuples_ppMC_2760GeV/data_updated_4.13.15.root";
const char* QCDOutFile    = "/net/hisrv0001/home/ilaflott/pp_MC_2760GeV_bTag_forests_ntuples/NTuples/kurts_QCDMC/QCD_updated_4.13.15.root";
const char* BJetOutFile   = "/net/hisrv0001/home/ilaflott/pp_MC_2760GeV_bTag_forests_ntuples/NTuples/halfOfficial_HFMC/BJet_halfOfficial_updated_4.13.15_schemeA.root";
const char* CJetOutFile   = "/net/hisrv0001/home/ilaflott/pp_MC_2760GeV_bTag_forests_ntuples/NTuples/halfOfficial_HFMC/CJet_halfOfficial_updated_4.13.15_schemeA.root";

const int weightsMode = 1; //1 for weight scheme A, anything else for scheme B
//const int weightsMode = -1;

const int QCDBins = 10;
//const int HFBins = 5;

const int pthatBin[] = {30,50,80,120,170,220,280,370,460,540};

const string pthatCut[] = 
  {
    "pthat<30",
    "pthat>=30 && pthat<50",
    "pthat>=50 && pthat<80",
    "pthat>=80 && pthat<120",
    "pthat>=120 && pthat<170",
    "pthat>=170 && pthat<220",
    "pthat>=220 && pthat<280",
    "pthat>=280 && pthat<370",
    "pthat>=370 && pthat<460",
    "pthat>=460 && pthat<540",
    "pthat>=540"
  };

// Cross sections for weighting in mb
const double xsection[] = 
  {
    0.2034,    // 15
    1.075e-2,  // 30
    1.025e-3,  // 50
    9.865e-5,  // 80
    1.129e-5,  // 120
    1.465e-6,  // 170
    2.837e-7,  // 220
    5.323e-8,  // 280
    5.934e-9,  // 370
    8.125e-10, // 460
    1.467e-10  // 540
  }; 

// Data storage variables
// akPu3PFJetAnalyzer/t
int   nref;
float rawpt[1000];
float jtpt[1000];
float jteta[1000];
float jtphi[1000];
float trackMax[1000];
float mupt[1000];
int   muN[1000];
float mueta[1000];
float muphi[1000];
float mudr[1000];
float muptrel[1000];
float discr_ssvHighEff[1000];
float discr_ssvHighPur[1000];
int   nsvtx[1000];
int   svtxntrk[1000];
float svtxdl[1000];
float svtxdls[1000];
float svtxm[1000];
float svtxpt[1000];

int nIPtrk[1000];//per jet
int nselIPtrk[1000];//per jet
int nIP;//per event
int ipJetIndex[10000];//per track
float ipPt[10000];
float ipProb0[10000];
float ipProb1[10000];
float ip2d[10000];
float ip2dSig[10000];
float ip3d[10000];
float ip3dSig[10000];
float ipDist2Jet[10000];
float ipDist2JetSig[10000];
float ipClosest2Jet[10000];
float trkChi2[10000];
float trkPt[10000];
float trkEta[10000];
float trkPhi[10000];

float pthat;                      // MC
float refpt[1000];                // MC
int   refparton_flavorForB[1000]; // MC

// hltanalysis/HltTree
int HLT_PAMu3_v1;
int HLT_PAMu7_v1;
int HLT_PAMu12_v1;
int HLT_PAMu3PFJet40_v1;

// hiEvtAnalyzer/HiTree
float vz;

// skimanalysis/HltTree
int pPAcollisionEventSelectionPA;
int pHBHENoiseFilter;

// New data storage vars
double nJtpt;
double nJteta;
double nJtphi;

double nTrackMax;

double nMupt;
double nRawpt;
double nMueta;
double nMuphi;
int    nMuN;
double nMudr;
double nMuptrel;

double nDiscr_ssvHighEff;
double nDiscr_ssvHighPur;

int nNIPtrk;
int nNselIPtrk;
int nNIP;
int nIPJetIndex[10000];
double nIPPt[10000];
double nIPProb0[10000];
//double nIPProb1[10000];
double nIP2d[10000];
double nIP2dsig[10000];
double nIP3d[10000];
double nIP3dsig[10000];
double nIPDist2Jet[10000];
//double nIPDist2JetSig[10000];
double nIPClosest2Jet[10000];
double nTrkChi2[10000];
double nTrkPt[10000];
double nTrkEta[10000];
double nTrkPhi[10000];

int    nNsvtx;
int    nSvtxntrk;
double nSvtxdl;
double nSvtxdls;
double nSvtxm;
double nSvtxpt;

//double nVz; //Event info
double nVz;
double nPthat; //MC
double nRefpt;                // MC
int    nRefparton_flavorForB; // MC

double nWeight;

string fileList; 

int result;
int dataType;

// Main functions
// Mode: 0-makeNTuple(), 1-mergeMCSamples()
// Type: 0-data, 1-QCD, 2-BJet, 3-CJet

int bTagNTuple_Original(int type)
{
  switch (type) 
    {
    case 0: fileList = dataFileList ; printf("\n you chose data"); weights_file = dataWeightsFile ; break ;
    case 1: fileList = QCDFileList  ; printf("\n you chose QCD") ; weights_file = QCDWeightsFile  ; break ;
    case 2: fileList = BJetFileList ; printf("\n you chose BJets"); weights_file = BJetWeightsFile ; break ;
    case 3: fileList = CJetFileList ; printf("\n you chose CJets"); weights_file = CJetWeightsFile ; break ;
    default:
      cerr << "Type must be from {0,1,2,3}" << endl;
      return -1;
    }

  result = makeNTuple(type);
  
  return result;
}

// Processes files and outputs an ntuple
int makeNTuple(int type)
{
  // Set global data/QCD/BJet/CJet switch
  dataType = type;  

  // Initialize output file
  TFile *outFile;

  switch (dataType) 
    {
    case 0: outFile = new TFile( Form("%s",dataOutFile) , "RECREATE" ); break;
    case 1: outFile = new TFile( Form("%s",QCDOutFile ) , "RECREATE" ); break;
    case 2: outFile = new TFile( Form("%s",BJetOutFile) , "RECREATE" ); break;
    case 3: outFile = new TFile( Form("%s",CJetOutFile) , "RECREATE" ); break;
    default:cerr<<"dataType not found"<<endl; return -1;
    }

  // New tree, new branches
  //
  TTree newTree("nt","nt");
  newBranches(&newTree);
  
  //outFile->cd();

  ifstream fileStream(fileList.c_str(), ifstream::in);
  string fileName;      
  fileStream >> fileName;
  
  //getline(fileStream, fileName);
  
  // For every file in file list, process trees
  while (!fileStream.eof()) 
    {
      // Open input file
      printf("\n Opening File: %s \n",fileName.c_str());
      TFile *inFile = TFile::Open( Form("%s",fileName.c_str() ) );
      
      // Open trees
      cout << "Opening Trees..." << endl;
      TTree *akPu3 = (TTree *)inFile->Get("akPu3PFJetAnalyzer/t");
      akPu3->AddFriend("hlt=hltanalysis/HltTree");
      akPu3->AddFriend("hiEvt=hiEvtAnalyzer/HiTree");
      akPu3->AddFriend("skim=skimanalysis/HltTree");
      akPu3->AddFriend("trk=ppTrack/trackTree");
      
      // Set branch addresses
      branchAddresses(akPu3);
      
      // Process every event
      int nEvents = akPu3->GetEntries();
      cout << nEvents << " events to loop over in " << fileName << endl;
      //int totNumTracks=0;
      //nEvents = 100;
      for (int i=0; i<nEvents; i++) 
	{
	  
	  if (i%10000 == 0 ) cout << "Processing Event " << i << endl;
	  
	  akPu3->GetEntry(i);

	  // Event Level Selection
	  if ((dataType == 0) && (0 || !pPAcollisionEventSelectionPA || !pHBHENoiseFilter || abs(vz)>15)) continue;
	  else if(abs(vz)>15)continue;// (dataType >= 1) 
	    	  
	  // Set weight
	  if (dataType == 0) nWeight = 1.0;
	  else nWeight = MCWeights(pthat);
	  
	  
	  //Event Info
	  nVz    = vz;
	  nNIP   = nIP;
	  nPthat = pthat; 
	  
	  int trackPosition=0;
	  //Jet Processing
	  for (int j=0; j<nref; j++) 
	    {
	      
	      trackPosition+=nselIPtrk[j];//at end of loop, this is number of tracks in our event
	      switch(dataType)
		{
		case 0: break;//no particle level selection for data
		case 1: break;//want all parton flavors from QCD file
		case 2:  if (fabs(refparton_flavorForB[j])!=5) continue; break;//want only bs from b file
		case 3:  if (fabs(refparton_flavorForB[j])!=4) continue; break;//want only cs from c file
		}
	      
	      ////jet level selection, might never apply it at this level...
	      //if(abs(jteta[j])<2.0||jtpt<40)continue;

	      //jet parameters
	      nJtpt  = jtpt[j];
	      nJteta = jteta[j];
	      nJtphi = jtphi[j];
	      nRawpt = rawpt[j];
	      if(dataType>=1) nRefpt                = refpt[j];
	      if(dataType>=1) nRefparton_flavorForB = refparton_flavorForB[j];
	      
	      //muons associated w/ jet
	      nMuN     = muN[j];
	      nMupt    = mupt[j];
	      nMueta   = mueta[j];
	      nMuphi   = muphi[j];
	      nMudr    = mudr[j];
	      nMuptrel = muptrel[j];
	      
	      //ssv discriminator values for jet
	      nDiscr_ssvHighEff = discr_ssvHighEff[j];
	      nDiscr_ssvHighPur = discr_ssvHighPur[j];
	     
	      //what is this variable?
	      nTrackMax = trackMax[j];
	      
	      //jet-track variables
	      nNIPtrk    =nIPtrk[j];
	      nNselIPtrk =nselIPtrk[j];
	      
	      //secondary vertex variables
	      nNsvtx    = nsvtx[j];
	      nSvtxntrk = svtxntrk[j];
	      nSvtxdl   = svtxdl[j];
	      nSvtxdls  = svtxdls[j];
	      nSvtxm    = svtxm[j];
	      nSvtxpt   = svtxpt[j];               
	      
	      //track based variables
	      if(doTracks)
		{
		  int counter=0;
		  for(int it = trackPosition-nselIPtrk[j];it<trackPosition;it++)
		    {

		      nIPJetIndex[counter]    = ipJetIndex[it];//this number reflects the index of the jet that the ip belongs to
		      nIPPt[counter]          = ipPt[it];
		      nIPProb0[counter]       = ipProb0[it];
		      //nIPProb1[counter]       = ipProb1[it];
		      nIP2d[counter]          = ip2d[it] ;
		      nIP2dsig[counter]       = ip2dSig[it];
		      nIP3d[counter]          = ip3d[it] ;
		      nIP3dsig[counter]       = ip3dSig[it];
		      nIPDist2Jet[counter]    = ipDist2Jet[it];
		      //nIPDist2JetSig[counter] = ipDist2JetSig[it];
		      nIPClosest2Jet[counter] = ipClosest2Jet[it];
		      nTrkChi2[counter]       =trkChi2[it];
		      nTrkPt[counter]         =trkPt[it];
		      nTrkEta[counter]        =trkEta[it];
		      nTrkPhi[counter]        =trkPhi[it];
		      counter++;	
		      
		    }
		}
	      
	    }//jetloop
	  
	  newTree.Fill();
	  
	}//eventloop
      

      // Cleanup
      cout << "closing " << fileName << endl;
      inFile->Close();

      fileStream >> fileName;
      //getline(fileStream, fileName);
    }
  
  // Write to output file
  outFile->cd();
  cout << "writing tree..." << newTree.GetName() <<endl;

  //kOverwrite overwrites backup/partially finished tree
  newTree.Write(newTree.GetName(), TObject::kOverwrite);
  
  // Cleanup
  cout << "cleaning up..." << endl;
  outFile->Close();
  fileStream.close();
  delete outFile;
  
  return 0;
}

//static basically means set once and then not again
// Return the corresponding weight for an event based on pthat
static double MCWeights(double MCPthat)
{
  static bool initialized = false;
  static double *weight;
  
  weight = new double[QCDBins+1];

  // Calculate the weight for each pthat bin in the QCD MC sample
  //if (!initialized) 
  if ( !(std::ifstream(weights_file.c_str())) )
    {
      if(!initialized) cout << "No weights_file found. Initializing weight function.\n";
      
      // Add QCD MC files to chain
      TChain *ch = new TChain("ak3PFJetAnalyzer/t");

      ifstream inStr(QCDFileList.c_str(), ifstream::in);
      string fileName;
      ofstream weightFile(weights_file.c_str(), ofstream::out);

      inStr >> fileName;
      while (!inStr.eof()) 
	{
          ch->Add(fileName.c_str());

          cout << "Added " << fileName << " to chain." << endl;
	  inStr >> fileName;
	}
      cout << "Done adding files to chain." << endl;
      

      // Count events across all files
      double pthatEntries[QCDBins+1];
      
      for (int i=0; i<QCDBins+1; i++) 
	{
	  pthatEntries[i] = ch->GetEntries( pthatCut[i].c_str() );
	  cout << "\tQCD pthatEntries with " << pthatCut[i] << ": " << pthatEntries[i] << endl;
	}

      // Modify event count for HF MC files
      if (dataType >= 2) 
	{
          cout << "Heavy flavor file. Updating weights." << endl;
          heavyJetWeights(pthatEntries);
	}
      

      // Calculate weights
      for (int i=0; i<QCDBins+1; i++) 
	{
	  //if there are no entries there's nothing to weigh
          if (pthatEntries[i]==0) weight[i] = 0.0;
	  else//there are entries to weight
	    {
	      if (i!=QCDBins) weight[i] = (xsection[i] - xsection[i+1])/pthatEntries[i];
	      else weight[i] = xsection[i]/pthatEntries[i];//so i dont run off the end of the xsec array, whats the xsec @ pthat=infty anyways? pretty sure it's zero.
	    }
	  
          cout << "weight[" << i << "] = " << weight[i] << endl;
          weightFile << weight[i] << endl;
	}

      // Cleanup
      inStr.close();
      delete ch;
      cout << "weight file closing..." << endl;
      weightFile.close();

      initialized = true;
      //cout << "MCWeights Initialized" << endl;
    }
  else
    {
      if(!initialized) 
	{
	  cout << "weights file detected. using previously calculated weights..." << endl;
	  cout << "using file: " << weights_file.c_str() << endl;
	}
      ifstream inWeights(weights_file.c_str(), ifstream::in);
      for(int i=0; i<QCDBins+1; i++) 
	{
	  inWeights >> weight[i];
	  if(!initialized)cout << "weight[" << i << "] = " << weight[i] << endl;
	}
      initialized = true;
    }

  int j=0;
  //say MCPthat = 53, then j = 2 after this loop, calls weight[2] = xsec[2] - xsec[3]/ pthatentries[2]
  //so weight[2] = \sigma_50 - \sigma_80 / entries[50,80]
  //largest j gets is j = QCD bins
  while (MCPthat>pthatBin[j] && j<QCDBins) j++;
  return weight[j];
}

static void heavyJetWeights(double *pthatEntries)
{
  // Add heavy flavor MC files to chain
  TChain *HFCh = new TChain("akPu3PFJetAnalyzer/t");
  TChain *HFCh_hiEvt = new TChain("hiEvtAnalyzer/HiTree");
  
  string HFfileList;
  int heavyFlavor;				 						 
  
  switch(dataType)
    {
    case 2: HFfileList=BJetFileList;heavyFlavor=5;break;
    case 3: HFfileList=CJetFileList;heavyFlavor=4;break;
    default: heavyFlavor=-1;break;
    }
  
  ifstream HFInStr(HFfileList.c_str(), ifstream::in);
  string HFFileName;
  
  HFInStr >> HFFileName;
  while (!HFInStr.eof())
    {
      HFCh->Add(HFFileName.c_str());
      HFCh_hiEvt->Add(HFFileName.c_str());
      cout << "Added " << HFFileName << " to chain." << endl;
      HFInStr >> HFFileName;
    }
  
  HFCh->AddFriend(HFCh_hiEvt);
  
  cout << "Done adding files to chain." << endl;
  
  // Count HF events across all files
  int HFPthatEntries[QCDBins+1];
  for (int i=0; i<QCDBins+1; i++) 
    {
      HFPthatEntries[i] = HFCh->GetEntries( pthatCut[i].c_str() );
      cout << "\tHFPthatEntries with " << pthatCut[i] << ": " << HFPthatEntries[i] << endl;
    }
  
  // Add QCD MC files to chain
  TChain *QCDCh = new TChain("akPu3PFJetAnalyzer/t");
  TChain *QCDCh_hiEvt = new TChain("hiEvtAnalyzer/HiTree");
  
  ifstream QCDInStr(QCDFileList.c_str(), ifstream::in);
  string QCDFileName;
  
  QCDInStr >> QCDFileName;
  while  (!QCDInStr.eof()) 
    {
      QCDCh->Add(QCDFileName.c_str());
      QCDCh_hiEvt->Add(QCDFileName.c_str());
      cout << "Added " << QCDFileName << " to chain." << endl;
      QCDInStr >> QCDFileName;
    }
  QCDCh->AddFriend(QCDCh_hiEvt);
  
  cout << "Done adding files to chain." << endl;
  
  
  // Calculate HF weight for each pthat bin
  double HFWeight[QCDBins+1];
  char HFJetsCut[200];
  
  TH1D *HFJetHist = new TH1D("HFJetHist", "HFJetHist", 1, 0, 10000);
  TH1D *QCDJetHist = new TH1D("QCDJetHist", "QCDJetHist", 1, 0, 10000);
  
  for (int i=0; i<QCDBins+1; i++) 
    {
      // Count (indirectly) number of b jets
      sprintf(HFJetsCut, "%s&&abs(jteta)<2&&abs(vz)<15&&refpt>0&&abs(refparton_flavorForB)==%d", pthatCut[i].c_str(), heavyFlavor);
      //sprintf(HFJetsCut, "%s&&abs(jteta)<2&&refpt>0&&abs(refparton_flavorForB)==%d", pthatCut[i].c_str(), heavyFlavor);
      //printf("%s\n",HFJetsCut);
      
      HFCh->Draw("jtpt>>HFJetHist", HFJetsCut, "goff");
      QCDCh->Draw("jtpt>>QCDJetHist", HFJetsCut, "goff");
      
      cout << "\tHF MC b/cJets: " << HFJetHist->Integral() << endl;
      cout << "\tQCD MC b/cJets: " << QCDJetHist->Integral() << endl;
      
      // Calculate b jets per **event
      double HFJetsPerEvent = (double)HFJetHist->Integral()/(double)HFPthatEntries[i];
      double QCDJetsPerEvent = (double)QCDJetHist->Integral()/(double)pthatEntries[i];
      
      cout << "\tHF MC b/cJets/Event: " << HFJetsPerEvent << endl;
      cout << "\tQCD MC b/cJets/Event: " << QCDJetsPerEvent << endl;
      
      // Check for NaN (IEEE method) and calculate HF weight
      if (HFJetsPerEvent != HFJetsPerEvent || QCDJetsPerEvent != QCDJetsPerEvent) HFWeight[i] = 0;
      else HFWeight[i] = HFJetsPerEvent/QCDJetsPerEvent;
      
      cout << "HFWeight for " << pthatCut[i] << ": " << HFWeight[i] << endl;
      
      HFJetHist->Reset();
      QCDJetHist->Reset();
      
    }

  // Weigh HF events by HF weight for each pthat bin, add to QCD events
  for (int i=0; i<QCDBins+1; i++) 
    {
      if(weightsMode == 1 ) pthatEntries[i] += HFWeight[i] * HFPthatEntries[i];  //schemeA
      else pthatEntries[i] = HFWeight[i] * HFPthatEntries[i]; //schemeB
      cout << "Effective pthat entries for pthat " << pthatCut[i] << ": " << pthatEntries[i] << endl;
    }

  // Cleanup
  delete HFCh;
  HFInStr.close();
  delete QCDCh;
  QCDInStr.close();
  delete HFJetHist;
  delete QCDJetHist;

  return;
}


// Create the branches for the new tree
static inline void newBranches(TTree *newTree) 
{
  //event specific
  newTree->Branch("weight", &nWeight, "weight/D");  
  newTree->Branch("vz", &nVz ,"vz/D");
  if(dataType>=1)newTree->Branch("pthat", &nPthat, "pthat/D"); // TEMPORARY
  
  //jet variables
  newTree->Branch("jtpt", &nJtpt, "jtpt/D");
  newTree->Branch("jteta", &nJteta, "jteta/D");
  newTree->Branch("jtphi", &nJtphi, "jtphi/D");
  if(dataType>=1)newTree->Branch("rawpt", &nRawpt, "rawpt/D");
  if(dataType>=1)newTree->Branch("refpt", &nRefpt, "refpt/D");
  if(dataType>=1)newTree->Branch("refparton_flavorForB", &nRefparton_flavorForB, "refparton_flavorForB/I");
  
  //muon variables
  newTree->Branch("mupt", &nMupt, "mupt/D");
  newTree->Branch("muN", &nMuN, "muN/I");
  newTree->Branch("mueta", &nMueta, "mueta/D");
  newTree->Branch("muphi", &nMuphi, "muphi/D");
  newTree->Branch("mudr", &nMudr, "mudr/D");
  newTree->Branch("muptrel", &nMuptrel, "muptrel/D");
  
  //ssv discriminator values
  newTree->Branch("discr_ssvHighEff", &nDiscr_ssvHighEff, "discr_ssvHighEff/D");
  newTree->Branch("discr_ssvHighPur", &nDiscr_ssvHighPur, "discr_ssvHighPur/D");
  
  //secondary vertex
  newTree->Branch("nsvtx", &nNsvtx, "nsvtx/I");
  newTree->Branch("svtxntrk", &nSvtxntrk, "svtxntrk/I");
  newTree->Branch("svtxdl", &nSvtxdl, "svtxdl/D");
  newTree->Branch("svtxdls", &nSvtxdls, "svtxdls/D");
  newTree->Branch("svtxm", &nSvtxm, "svtxm/D");
  newTree->Branch("svtxpt", &nSvtxpt, "svtxpt/D");
  
  //tracks
  newTree->Branch("trackMax", &nTrackMax, "trackMax/D");
  newTree->Branch("nIPtrk" ,&nNIPtrk  , "nIPtrk/I");
  newTree->Branch("nselIPtrk" ,&nNselIPtrk  , "nselIPtrk/I");
  newTree->Branch("nIP" ,&nNIP  , "nIP/I");
  if(doTracks)
    {
      newTree->Branch("ipJetIndex" ,&nIPJetIndex  , "ipJetIndex[nIP]/I");
      newTree->Branch("ipPt" ,&nIPPt  , "ipPt[nIP]/I");
      newTree->Branch("ipProb0" ,&nIPProb0  , "ipProb0[nIP]/I");
      //newTree->Branch("ipProb1" ,&nIPProb1  , "ipProb1[nIP]/I");
      newTree->Branch("ip2d" ,&nIP2d  ," ip2d[nIP]/D");
      newTree->Branch("ip2dSig",&nIP2dsig , "ip2dSig[nIP]/D");
      newTree->Branch("ip3d" ,&nIP3d  , "ip3d[nIP]/D");
      newTree->Branch("ip3dSig",&nIP3dsig , "ip3dSig[nIP]/D");
      newTree->Branch("ipDist2Jet",&nIPDist2Jet , "ipDist2Jet[nIP]/D");
      //newTree->Branch("ipDist2JetSig",&nIPDist2JetSig , "ipDist2JetSig[nIP]/D");
      newTree->Branch("ipCloset2Jet",&nIPClosest2Jet , "ipClosest2Jet[nIP]/D");
      //new track variables from ppTrack tree
      newTree->Branch("trkChi2",&nTrkChi2,"trkChi2[nIP]/D");
      newTree->Branch("trkPt"  ,&nTrkPt  ,"trkPt[nIP]/D");
      newTree->Branch("trkEta" ,&nTrkEta ,"trkEta[nIP]/D");
      newTree->Branch("trkPhi" ,&nTrkPhi ,"trkPhi[nIP]/D");
  }  

  //HLT
  newTree->Branch("HLT_PAMu3_v1", &HLT_PAMu3_v1, "HLT_PAMu3_v1/I");
  newTree->Branch("HLT_PAMu7_v1", &HLT_PAMu7_v1, "HLT_PAMu7_v1/I");
  newTree->Branch("HLT_PAMu12_v1", &HLT_PAMu12_v1, "HLT_PAMu12_v1/I");
  newTree->Branch("HLT_PAMu3PFJet40_v1", &HLT_PAMu3PFJet40_v1, "HLT_PAMu3PFJet40_v1/I");
  
  return;
}

// Set all the input branch addresses
static inline void branchAddresses(TTree *akPu3) 
{


  //EVENT INFO
  akPu3->SetBranchAddress("nref", &nref);
  akPu3->SetBranchAddress("vz", &vz);
  akPu3->SetBranchAddress("pPAcollisionEventSelectionPA", &pPAcollisionEventSelectionPA);
  akPu3->SetBranchAddress("pHBHENoiseFilter", &pHBHENoiseFilter);
  if(dataType>=1)akPu3->SetBranchAddress("pthat", &pthat);

  //JET
  akPu3->SetBranchAddress("jtpt" , &jtpt  );
  akPu3->SetBranchAddress("jteta", &jteta );
  akPu3->SetBranchAddress("jtphi", &jtphi );
  akPu3->SetBranchAddress("rawpt", &rawpt);
  if(dataType>=1)akPu3->SetBranchAddress("refpt", &refpt);
  if(dataType>=1)akPu3->SetBranchAddress("refparton_flavorForB", &refparton_flavorForB);
  
  //TRACK
  akPu3->SetBranchAddress("trackMax", &trackMax);
  akPu3->SetBranchAddress("nIPtrk"    ,&nIPtrk);
  akPu3->SetBranchAddress("nselIPtrk" ,&nselIPtrk);
  akPu3->SetBranchAddress("nIP"       ,&nIP);
  if(doTracks)
    {
      akPu3->SetBranchAddress("ipJetIndex",&ipJetIndex);
      akPu3->SetBranchAddress("ipProb0"   , &ipProb0    );
      //akPu3->SetBranchAddress("ipProb1"   , &ipProb1    );
      akPu3->SetBranchAddress("ipPt"   , &ipPt    );
      akPu3->SetBranchAddress("ip2d"   , &ip2d    );
      akPu3->SetBranchAddress("ip2dSig", &ip2dSig );
      akPu3->SetBranchAddress("ip3d"   , &ip3d    );
      akPu3->SetBranchAddress("ip3dSig", &ip3dSig );
      akPu3->SetBranchAddress("ipDist2Jet", &ipDist2Jet );
      //akPu3->SetBranchAddress("ipDist2JetSig", &ipDist2JetSig );
      akPu3->SetBranchAddress("ipClosest2Jet", &ipClosest2Jet );
      //new track variables from ppTree
      akPu3->SetBranchAddress("trkChi2",&trkChi2);
      akPu3->SetBranchAddress("trkPt"  ,&trkPt  );
      akPu3->SetBranchAddress("trkEta" ,&trkEta );
      akPu3->SetBranchAddress("trkPhi" ,&trkPhi );
    }

  //HLT
  akPu3->SetBranchAddress("HLT_PAMu3_v1", &HLT_PAMu3_v1);
  akPu3->SetBranchAddress("HLT_PAMu7_v1", &HLT_PAMu7_v1);
  akPu3->SetBranchAddress("HLT_PAMu12_v1", &HLT_PAMu12_v1);
  akPu3->SetBranchAddress("HLT_PAMu3PFJet40_v1", &HLT_PAMu3PFJet40_v1);
    
  //MUON
  akPu3->SetBranchAddress("mupt"   , &mupt   );
  akPu3->SetBranchAddress("muN"    , &muN    );
  akPu3->SetBranchAddress("mueta"  , &mueta  );
  akPu3->SetBranchAddress("muphi"  , &muphi  );
  akPu3->SetBranchAddress("mudr"   , &mudr   );
  akPu3->SetBranchAddress("muptrel", &muptrel);
  
  //DISCRIMINATORS
  akPu3->SetBranchAddress("discr_ssvHighEff", &discr_ssvHighEff);
  akPu3->SetBranchAddress("discr_ssvHighPur", &discr_ssvHighPur);
  
  //SVTX
  akPu3->SetBranchAddress("nsvtx"   , &nsvtx   );
  akPu3->SetBranchAddress("svtxntrk", &svtxntrk );
  akPu3->SetBranchAddress("svtxdl"  , &svtxdl   );
  akPu3->SetBranchAddress("svtxdls" , &svtxdls  );
  akPu3->SetBranchAddress("svtxm"   , &svtxm    );
  akPu3->SetBranchAddress("svtxpt"  , &svtxpt   );
  
  return;
}


//Leo's Old Notes about Jet Reweighting
/* 
  Heavy jet weighting
  Helper method for MCWeights when the MC is BJet or CJet
  Takes in QCD pthatEntries and augments them w/ weighted HF pthatEntries
  Kurt Jung: counting QCD B-Jets per event
  Kurt Jung: and then counting B-Jets per event in the HF MC
  Kurt Jung: and scaling back the HF such that the number of b-jets per event matches the QCD jet MC
  Kurt Jung: so you have high statistics b-jets without biasing the sample"
  
  ex.
  10 events in BJet MC, 20 jets
  10 events in QCD MC, 20 jets
  
  16 BJets in BJet MC
  2 BJets in QCD MC
  
  16/10 = 1.6 BJets/event in BJet MC
  2/10 = 0.2 BJets/event in QCD MC
  
  1.6/0.2 = 8
  
  BJet MC weight = xSec / (10 QCD Events + (10 BJet Events * 8))
  QCD MC weight = xSec / 10 QCD Events
  
  20 jets in BJet MC * xSec / 90 = 2/9 xSec
  20 jets in QCD MC * xSec / 10 = 2 xSec
  
  Total integral = (2 * 2/9) * xSec
  
  So the effect of the heavy jet weighting is to decrease
  the weight for the BJet MC jets -> better statistics/error bars
  and still have a similar spectrum
  
  It seems like it works, but I still donk't completely understand why we scale by BJets/Event.
*/

////////OLDERCODE
//
//
//
//// Merges the QCD, BJet, and CJet MC ntuples into a single file
//int mergeMCSamples()
//{
//
//  TFile *QCD  = TFile::Open( Form("%s",  QCDOutFile ) );
//  TFile *BJet = TFile::Open( Form("%s", BJetOutFile ) );
//  TFile *CJet = TFile::Open( Form("%s", CJetOutFile ) );
//
//  if (!QCD || !BJet || !CJet) 
//    {
//      cerr << "mergeMCSamples() requires QCD.root, QCD_2.root, BJet.root, and CJet.root to run." << endl;
//      return -1;
//    }
//
//  // Add MC files to chain
//  TChain *ch = new TChain("nt");
//  ch->Add(Form("%s",  QCDOutFile ) );
//  ch->Add(Form("%s", BJetOutFile ) );
//  ch->Add(Form("%s", CJetOutFile ) );
//
//  // Merge chain and output to new file
//  TFile *outFile = new TFile(Form("%s", mergedOutFile), "RECREATE");
//  outFile->cd();
//  ch->Merge(outFile, 0, "keep");
//
//  // Cleanup
//  QCD->Close();
//  BJet->Close();
//  CJet->Close();
//  outFile->Close();
//  delete ch;
//  delete outFile;
//
//  return 0;
//}
//
//int calculateWeights(int type)
//{
//  TTree newTree("nt","nt");
//  newBranches(&newTree);
//  MCWeights(pthat);
//  return 0;
//  
//}
//
//int forestStatistics(int type)
//{
//  int totalCount = 0;
//  int numFile = 0;
//  if (type <= 4 && type >= 1) fileList = QCDFileList; 
//  ifstream inStr(fileList.c_str(), ifstream::in);
//  string fileName;
//  
//  inStr >> fileName;
//  
//  TFile *file = TFile::Open(fileName.c_str());
// 
//  TTree *srctree = (TTree*)file->Get("hiEvtAnalyzer/HiTree");
//
//  while (!inStr.eof()) 
//    {
//      if (type == 0 )//data isn't generated in pt bins, just spit out the total number of entries
//	{
//	  numFile += 1;
//	  totalCount += srctree->GetEntries();
//	  if (numFile%50==0) 
//	    {
//	      cout << "Number of Events for file #"<<numFile << " : " << srctree->GetEntries() << endl;
//	      
//	    }
//	  
//	  inStr >> fileName;
//	  
//	  file = TFile::Open(fileName.c_str());
//	  srctree = (TTree*)file->Get("hiEvtAnalyzer/HiTree");
//	}
//      else//MC Case, spit out Entries by the file, which also corresponds to pt bins
//	{
//	  numFile += 1;
//	  totalCount += srctree->GetEntries();
//	  
//	  cout << "Number of Events for "<< fileName.c_str() << endl;
//	  cout << "is: "<< srctree->GetEntries() << endl;
//
//	  if (!inStr.eof()) cout <<"not the end of the file"<<endl;
//	  else
//	    {
//	      cout << "the end of the file is nigh!" <<endl;
//	      cout << fileName.c_str() << endl;
//	    }
//	  inStr >> fileName;
//	  if (!inStr.eof()) cout <<"not the end of the file"<<endl;
//	  else
//	    {
//	      cout << "the end of the file is nigh!" <<endl;
//	      cout << fileName.c_str() << endl;
//	    }
//	  file = TFile::Open(fileName.c_str());
//	  srctree = (TTree*)file->Get("hiEvtAnalyzer/HiTree");
//	}
//      
//     
//    }
//
//  cout << "Done counting data Events" << endl;
//  cout << "final count = " << totalCount << endl;
//
//  
//  return 0;
//}
