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
#include "TROOT.h"
#include "TStopwatch.h"
#include "TTree.h"
#include "TTreeCache.h"
#include "TChain.h"
#include "TH1D.h"

using namespace std;

// Function declarations
int makeNTuple(int type);
int MCCounts(int type);
int NTupleWeights(int type);
int NTupleTest(int type);//crude script, doesn't even use the type
//static void heavyJetWeights(double *pthatEntries);
static inline void newBranches(TTree *newTree);
static inline void branchAddresses(TTree *akPu3);

//int impactParameterExploration(int type);
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
const bool doTracks=true;
const bool doMostSignificantTracks=true;

const double pi = 3.1415926535897932384626433832795028841971693993751058209749445923078164062862089986280348253421170679;
// Macro settings/constants

const string fileListPath   = "filelists/";
//const string outFilePath    = "/mnt/hadoop/cms/store/user/ilaflott/Leos_Analysis_NTuples/";
const string outFilePath = "";
const string weightFilePath = "/mnt/hadoop/cms/store/user/ilaflott/Leos_Analysis_weight_info/";
const string NTuplePath = "/mnt/hadoop/cms/store/user/ilaflott/Leos_Analysis_NTuples/";
//const string NTuplePath = "";/*debug*/

const string dataFileList = "ppMuon2013A_runForest_filelist.txt";
const string QCDFileList  = "QCDJets_noVsJets_filelist.txt";
const string QCDFileList_1   = "QCDJets_noVsJets_filelist_1.txt";
const string QCDFileList_2   = "QCDJets_noVsJets_filelist_2.txt";
const string QCDFileList_3   = "QCDJets_noVsJets_filelist_3.txt";
const string QCDFileList_4   = "QCDJets_noVsJets_filelist_4.txt";
const string QCDFileList_5   = "QCDJets_noVsJets_filelist_5.txt";
const string QCDFileList_6   = "QCDJets_noVsJets_filelist_6.txt";
const string QCDFileList_7   = "QCDJets_noVsJets_filelist_7.txt";
const string QCDFileList_8   = "QCDJets_noVsJets_filelist_8.txt";
const string QCDFileList_9   = "QCDJets_noVsJets_filelist_9.txt";
const string QCDFileList_10  = "QCDJets_noVsJets_filelist_10.txt";
const string QCDFileList_11  = "QCDJets_noVsJets_filelist_11.txt";
const string BJetFileList = "BJets_filelist.txt";
const string CJetFileList = "CJets_filelist.txt";

//Weight Files
//the weights only change if i change the event selection, i.e. the weights i've already calculated should be fine for any NTuples i make in the future
const string dataWeightsFile = "data_weights.txt";
const string QCDWeightsFile  = "QCDJets_weights.txt";
const string BJetWeightsFile = "BJets_weights.txt";
const string CJetWeightsFile = "CJets_weights.txt";

//Weight Information Files,  #Events per pthat bin, #BJets per pthat bin, #Cjets per pthat bin
const string QCD_NEventsFile = "QCDJets_NEvents.txt";	  
const string QCD_NBJetsFile  = "QCDJets_NBJets.txt";
const string QCD_NCJetsFile  = "QCDJets_NCJets.txt";

const string C_NEventsFile = "CJets_NEvents.txt";	  
const string C_NCJetsFile  = "CJets_NCJets.txt";

const string B_NEventsFile = "BJets_NEvents.txt";	  
const string B_NBJetsFile  = "BJets_NBJets.txt";

//NTuple Files
const string dataNTuple   = "data_NTuple_8.27.15.root";
const string QCDNTuple_noWeights    = "QCDJets_NTuple_noWeights.root";
const string BJetNTuple_noWeights   =   "BJets_NTuple_noWeights.root";
const string CJetNTuple_noWeights   =   "CJets_NTuple_noWeights.root";
const string QCDNTuple_withWeights    = "QCDJets_NTuple_withWeights.root";
const string BJetNTuple_withWeights   =   "BJets_NTuple_withWeights.root";
const string CJetNTuple_withWeights   =   "CJets_NTuple_withWeights.root";

/////*debug*/
//const string dataNTuple   = "data_NTuple_TEST0.root";
//const string QCDNTuple    = "QCDJets_NTuple_TEST.root";
//const string BJetNTuple   = "BJets_NTuple_TEST.root";
//const string CJetNTuple   = "CJets_NTuple_TEST.root";

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
//JetAnalyzer
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
//float discr_csvMva[1000];
//float discr_csvSimple[1000];
//float discr_muByIp3[1000];
//float discr_muByPt[1000];
//float discr_prob[1000];
//float discr_probb[1000];
//float discr_tcHighEff[1000];
//float discr_tcHighPur[1000];

int   nsvtx[1000];
int   svtxntrk[1000];
float svtxdl[1000];
float svtxdls[1000];
float svtx2Ddl[1000];
float svtx2Ddls[1000];
float svtxm[1000];
float svtxpt[1000];
float svtxXPos[1000];
float svtxYPos[1000];
float svtxZPos[1000];

int nIP;//per event
int nIPtrk[1000];//per jet
int nselIPtrk[1000];//per jet
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
int ipNHitPixel[10000];
int ipNHitStrip[10000];

float trkChi2[10000];
float trkPt[10000];
float trkEta[10000];
float trkPhi[10000];
float trkDz1[10000];
float trkDxy1[10000];

float deltaRtrk2Jet[10000];//i compute this one

float pthat;                      // MC Only
float weight;                      // MC Only
float refpt[1000];                // MC Only
int   refparton_flavorForB[1000]; // MC Only

// hltanalysis/HltTree
int HLT_PAMu3_v1;
int HLT_PAMu7_v1;
int HLT_PAMu12_v1;
int HLT_PABTagMu_Jet20_Mu4_v1;
int HLT_PAMu3PFJet20_v1;
int HLT_PAMu3PFJet40_v1;
int HLT_PAMu7PFJet20_v1;

// hiEvtAnalyzer/HiTree
float vz ; 
int evt  ;
int run  ; 
//int lumi ;
// skimanalysis/HltTree
int pPAcollisionEventSelectionPA;
int pHBHENoiseFilter;
// New data storage vars
int nNref;
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
int nNIPtrk;
int nNselIPtrk;
int nNIP;
double n1stMost2dSigTrk;  
double n2ndMost2dSigTrk;  
double n3rdMost2dSigTrk;  
double n1stIP2dTrk;  
double n2ndIP2dTrk;  
double n3rdIP2dTrk;  
double n1stMost3dSigTrk;  
double n2ndMost3dSigTrk;  
double n3rdMost3dSigTrk;  
double n1stIP3dTrk;  
double n2ndIP3dTrk;  
double n3rdIP3dTrk;     
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
double nTrkDxy[10000];
double nTrkDz[10000];
double nDeltaEta[10000];
double nDeltaPhi[10000];
double nDeltaRtrk2Jet[10000];
int nIpNHitPixel[10000];
int nIpNHitStrip[10000];


int    nNsvtx;
int    nSvtxntrk;
double nSvtxdl;
double nSvtxdls;
double nSvtx2Ddl;
double nSvtx2Ddls;
double nSvtxm;
double nSvtxpt;
double nSvtxXPos;
double nSvtxYPos;
double nSvtxZPos;
double nSvtxEta;
double nSvtxPhi;
double nSvtxDeltaEta;
double nSvtxDeltaPhi;
double nSvtxDeltaR2Jet;
double nDiscr_ssvHighEff;
double nDiscr_ssvHighPur;

double nVz;  //Event Information
double nPthat;                // MC Only
double nRefpt;                // MC Only
int Evt;
//int Lumi;
int Run;
int    nRefparton_flavorForB; // MC Only
double nWeight;

string weights_file;
string fileList; 


int result;
int dataType;

// Main functions
// Mode: 0-makeNTuple(), 1-mergeMCSamples()
// Type: 0-data, 1-QCD, 2-BJet, 3-CJet
int bTagNTuple(int job, int type)
{
  dataType = type;  
  switch (type) 
    {
      //notice: string fileList is compatible with const char*, which has the + operator overloaded
    case 0 : fileList = fileListPath + dataFileList    ; cout << "you chose " << dataFileList   << endl   ;   break ;
    case 1 : fileList = fileListPath + QCDFileList     ; cout << "you chose " << QCDFileList    << endl   ;   break ;
    case 2 : fileList = fileListPath + BJetFileList    ; cout << "you chose " << BJetFileList   << endl   ;   break ;
    case 3 : fileList = fileListPath + CJetFileList    ; cout << "you chose " << CJetFileList   << endl   ;   break ;
      //QCD PIECES
    case 4 : fileList = fileListPath + QCDFileList_1   ; cout << "you chose " << QCDFileList_1  << endl   ;   break ;
    case 5 : fileList = fileListPath + QCDFileList_2   ; cout << "you chose " << QCDFileList_2  << endl   ;   break ;
    case 6 : fileList = fileListPath + QCDFileList_3   ; cout << "you chose " << QCDFileList_3  << endl   ;   break ;
    case 7 : fileList = fileListPath + QCDFileList_4   ; cout << "you chose " << QCDFileList_4  << endl   ;   break ;
    case 8 : fileList = fileListPath + QCDFileList_5   ; cout << "you chose " << QCDFileList_5  << endl   ;   break ;
    case 9 : fileList = fileListPath + QCDFileList_6   ; cout << "you chose " << QCDFileList_6  << endl   ;   break ;
    case 10: fileList = fileListPath + QCDFileList_7   ; cout << "you chose " << QCDFileList_7  << endl   ;   break ;
    case 11: fileList = fileListPath + QCDFileList_8   ; cout << "you chose " << QCDFileList_8  << endl   ;   break ;
    case 12: fileList = fileListPath + QCDFileList_9   ; cout << "you chose " << QCDFileList_9  << endl   ;   break ;
    case 13: fileList = fileListPath + QCDFileList_10  ; cout << "you chose " << QCDFileList_10 << endl   ;   break ;
    case 14: fileList = fileListPath + QCDFileList_11  ; cout << "you chose " << QCDFileList_11 << endl   ;   break ;
    default: cerr << "Type must be from {0,1,2,3}" << endl ; return -1 ;

    }
  switch (job)
    {
    case 0:cout << "You Chose makeNTuple" << endl; result = makeNTuple(type) ; break ; 
    case 1:cout << "You Chose MCCounts" << endl; result = MCCounts(type) ; break ; 
    case 2:cout << "You Chose Apply NTuple Weights" << endl; result = NTupleWeights(type) ; break ; 
    case 3:cout << "You Chose to Test NTuple" << endl; result = NTupleTest(type) ; break ; 
    default: cerr << "Job must be from {0,1,2,3}" << endl; return -1 ;
    }

  return result;
}

// Processes files and outputs an ntuple
int makeNTuple(int type)
{
  // Set global data/QCD/BJet/CJet switch
  // Initialize output file

  TFile *outFile;

  string outFileName;
  switch (dataType) 
    {
    case 0 : outFileName = outFilePath + dataNTuple ; cout << outFileName << endl ;  outFile = new TFile( Form( "%s" , outFileName.c_str() ) , "RECREATE" ); break ;
      //    case 1: outFileName = outFilePath + QCDNTuple  ; cout << outFileName << endl ;  outFile = new TFile( Form( "%s" , outFileName.c_str() ) , "RECREATE" ); break ;
    case 1 : cout << "consider submitting makeNTuple in 11 pieces using types 4-14..." << endl; return -1;
    case 2 : outFileName = outFilePath + BJetNTuple_noWeights ;   outFile = new TFile( Form( "%s" , outFileName.c_str() ) , "RECREATE" ); break ;
    case 3 : outFileName = outFilePath + CJetNTuple_noWeights ;   outFile = new TFile( Form( "%s" , outFileName.c_str() ) , "RECREATE" ); break ;
    case 4 : outFileName = outFilePath + "P_1"  + QCDNTuple_noWeights  ;   outFile = new TFile( Form( "%s" , outFileName.c_str() ) , "RECREATE" ); break ;
    case 5 : outFileName = outFilePath + "P_2"  + QCDNTuple_noWeights  ;   outFile = new TFile( Form( "%s" , outFileName.c_str() ) , "RECREATE" ); break ;
    case 6 : outFileName = outFilePath + "P_3"  + QCDNTuple_noWeights  ;   outFile = new TFile( Form( "%s" , outFileName.c_str() ) , "RECREATE" ); break ;
    case 7 : outFileName = outFilePath + "P_4"  + QCDNTuple_noWeights  ;   outFile = new TFile( Form( "%s" , outFileName.c_str() ) , "RECREATE" ); break ;
    case 8 : outFileName = outFilePath + "P_5"  + QCDNTuple_noWeights  ;   outFile = new TFile( Form( "%s" , outFileName.c_str() ) , "RECREATE" ); break ;
    case 9 : outFileName = outFilePath + "P_6"  + QCDNTuple_noWeights  ;   outFile = new TFile( Form( "%s" , outFileName.c_str() ) , "RECREATE" ); break ;
    case 10: outFileName = outFilePath + "P_7"  + QCDNTuple_noWeights  ;   outFile = new TFile( Form( "%s" , outFileName.c_str() ) , "RECREATE" ); break ;
    case 11: outFileName = outFilePath + "P_8"  + QCDNTuple_noWeights  ;   outFile = new TFile( Form( "%s" , outFileName.c_str() ) , "RECREATE" ); break ;
    case 12: outFileName = outFilePath + "P_9"  + QCDNTuple_noWeights  ;   outFile = new TFile( Form( "%s" , outFileName.c_str() ) , "RECREATE" ); break ;
    case 13: outFileName = outFilePath + "P_10" + QCDNTuple_noWeights  ;   outFile = new TFile( Form( "%s" , outFileName.c_str() ) , "RECREATE" ); break ;
    case 14: outFileName = outFilePath + "P_11" + QCDNTuple_noWeights  ;   outFile = new TFile( Form( "%s" , outFileName.c_str() ) , "RECREATE" ); break ;
    default:cerr<<"dataType not found"<<endl; return -1;
    }
  cout << "decalring new tree+branches" << endl;
  // Jet tree
  TTree * newTree = new TTree("nt","nt");
  newBranches(newTree);
  newTree->SetDirectory(0);
  // vtx tree for weighting
  TTree * vtxTree = new TTree("vertex","vertex");
  vtxTree->Branch("evtVz", &nVz ,"evtVz/D");									       
  vtxTree->SetDirectory(0);

  cout << "fileList is: " << fileList << endl;
  //grab filename
  ifstream fileStream(fileList.c_str(), ifstream::in);
  string fileName;      
  fileStream >> fileName;//this needed before loop so we don't miss the last file available
  int file_number = 0;

  // For every file in file list, process trees
  cout << "beginning file loop" << endl;

  //TStopwatch * clock = new TStopwatch(); /*debug*/
  //clock->Start();                        /*debug*/

  for(int kkk = 0 ; kkk < 1 ; kkk++)/*debug*/
    //while (!fileStream.eof()) 
    {
      // Open input file
      TFile *inFile = TFile::Open( Form("%s",fileName.c_str() ) );
      if(file_number%100==0 ) cout << "fileName is: " << fileName << endl;
      //if(file_number%100==0 ) cout << "Opening the " << file_number<< "th file" << endl;
      if(file_number%100==0 ) cout << "Opening the " << file_number<< "th file" << endl;/*debug*/
      file_number++;
      
      //Open trees
      //if(file_number%1000==0 )cout << "Opening Trees..." << endl;/*debug*/
      TTree *akPu3 = (TTree *)inFile->Get("akPu3PFJetAnalyzer/t");
      akPu3->AddFriend("hlt=hltanalysis/HltTree");
      akPu3->AddFriend("hiEvt=hiEvtAnalyzer/HiTree");
      akPu3->AddFriend("skim=skimanalysis/HltTree");
      akPu3->AddFriend("trk=ppTrack/trackTree");

      // Set branch addresses
      branchAddresses(akPu3);
      
      // Process every event
      int nEvents = akPu3->GetEntries();
      if(file_number%100==0 )cout << nEvents << " events to loop over in " << fileName << endl;
      
      nEvents = 10;/*debug*/
      for (int i=0; i<nEvents; i++) 
	{
	  //if (i%10000 == 0 && i != 0) cout << "Processing Event " << i << endl;
	  akPu3->GetEntry(i);

	  // Event Level Selection
	  //used to only apply  pPAcollisionEventSelectionPA+pHBHENoiseFileter to data.... 
	  //but they should apply to MC as well!
	  if ( !pPAcollisionEventSelectionPA || !pHBHENoiseFilter || abs(vz)>15)) continue;

	  //Event Info
	  nVz    = vz;
	  nNIP   = nIP;
	  nPthat = pthat; 
	 	       
	  Evt  =evt;   
	  //Lumi =lumi;  
	  Run  =run; 


	  int trackPosition=0;
	  //Jet Processing	  
	  for (int j=0; j<nref; j++) 
	    {
	      int counter = 0; //for counting the number of tracks in a jet after cuts
	      trackPosition+=nselIPtrk[j];//at end of loop, this is number of tracks in our event

	      switch(dataType)
		{
		case 0: break;//no particle level selection for data
		case 1: break;//want all parton flavors from QCD file
		case 2:  if (abs(refparton_flavorForB[j])!=5) continue; break;//want only bs from b file
		case 3:  if (abs(refparton_flavorForB[j])!=4) continue; break;//want only cs from c file
		}
	      
	      ////jet level selection, might never apply it at this level...
	      //if(abs(jteta[j])>2.0||jtpt<40)continue;

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
	      
	      //what is this variable?
	      nTrackMax = trackMax[j];
	      
	      //jet-track variables
	      nNIPtrk    = nIPtrk[j];
	      nNselIPtrk = nselIPtrk[j];
	      
	      //secondary vertex variables, not being filled correct
	      nNsvtx    = nsvtx[j];
	      nSvtxntrk = svtxntrk[j];
	      nSvtxdl   = svtxdl[j];
	      nSvtxdls  = svtxdls[j];
	      nSvtx2Ddl  = svtx2Ddl[j];
	      nSvtx2Ddls  = svtx2Ddls[j];
	      nSvtxm    = svtxm[j];
	      nSvtxpt   = svtxpt[j];           
	      nSvtxXPos = svtxXPos[j];
	      nSvtxYPos = svtxYPos[j];
	      nSvtxZPos = svtxZPos[j];

	      //compute secondary vertex eta, phi + deltaR to Jet
//	      cout << "nSvtxZPos = " << nSvtxZPos << endl;
//	      cout << "nSvtxYPos = " << nSvtxYPos << endl;
//	      cout << "nSvtxXPos = " << nSvtxXPos << endl;
//	      cout << "denominator = " << sqrt(nSvtxXPos*nSvtxXPos+nSvtxYPos*nSvtxYPos+nSvtxZPos*nSvtxZPos) << endl;
//	      cout << "inside fracion = " << nSvtxZPos/sqrt(nSvtxXPos*nSvtxXPos+nSvtxYPos*nSvtxYPos+nSvtxZPos*nSvtxZPos) << endl;
//	      cout << "acos = " << acos( nSvtxZPos/sqrt(nSvtxXPos*nSvtxXPos+nSvtxYPos*nSvtxYPos+nSvtxZPos*nSvtxZPos) ) << endl;
//	      cout << "tan = " << tan( 0.5*acos( nSvtxZPos/sqrt(nSvtxXPos*nSvtxXPos+nSvtxYPos*nSvtxYPos+nSvtxZPos*nSvtxZPos) ) )<< endl;//either 0 or nan, why?
//	      cout << "log = " << -log( tan( 0.5*acos( nSvtxZPos/sqrt(nSvtxXPos*nSvtxXPos+nSvtxYPos*nSvtxYPos+nSvtxZPos*nSvtxZPos) ) ) )<< endl;
	      
	      nSvtxEta = -log( tan( 0.5*acos( nSvtxZPos/sqrt(nSvtxXPos*nSvtxXPos+nSvtxYPos*nSvtxYPos+nSvtxZPos*nSvtxZPos) ) ) );
	      nSvtxPhi = acos(nSvtxXPos/sqrt(nSvtxXPos*nSvtxXPos+nSvtxYPos*nSvtxYPos));
	      if(nSvtxYPos<0) nSvtxPhi = -nSvtxPhi;
	      
	      //NaN check, IEEE method
	      if(nSvtxEta!=nSvtxEta || nSvtxPhi!=nSvtxPhi)
		{
		  //if NaN, stick in garbage values for later analysis of these events while not affecting the distributions I'll want to see
		  nSvtxEta=-19;
		  nSvtxPhi=-2*pi; // -pi<phi<pi
		  nSvtxDeltaPhi=-pi;// 0<deltaphi
		  nSvtxDeltaEta=-1; // 0<deltaEta
		  nSvtxDeltaR2Jet=-1; //0<deltaR2Jet
		}
	      else //no nans, we're fine
		{
		  nSvtxDeltaPhi = jtphi[j]-nSvtxPhi;
		  if(fabs(nSvtxDeltaPhi)>pi) nSvtxDeltaPhi = 2*pi-fabs(nSvtxDeltaPhi);
		  else nSvtxDeltaPhi = fabs(nSvtxDeltaPhi);	      
		  nSvtxDeltaEta = fabs(jteta[j]-nSvtxEta);
		  nSvtxDeltaR2Jet = sqrt(nSvtxDeltaEta*nSvtxDeltaEta+nSvtxDeltaPhi*nSvtxDeltaPhi);
		}
	      //	      cout << nSvtxEta << endl;
	      //discriminator values
	      nDiscr_ssvHighEff = discr_ssvHighEff[j];
	      nDiscr_ssvHighPur = discr_ssvHighPur[j];

	      //cout << "doing tracks for this jet" << endl;
	      //track based variables
	      if(doTracks)//tracks take awhile to run, may want to turn it off in the future?
		{
		  counter = 0;
		  n1stMost2dSigTrk = 0;
		  n2ndMost2dSigTrk = 0;
		  n3rdMost2dSigTrk = 0;
		  n1stMost3dSigTrk = 0;
		  n2ndMost3dSigTrk = 0;
		  n3rdMost3dSigTrk = 0;
		  n1stIP2dTrk	= 0 ;  
		  n2ndIP2dTrk	= 0 ;  
		  n3rdIP2dTrk	= 0 ;  
		  n1stIP3dTrk	= 0 ;  
		  n2ndIP3dTrk	= 0 ;  
		  n3rdIP3dTrk	= 0 ; 
		  
		  for(int it = trackPosition-nselIPtrk[j] ; it < trackPosition ; it++)
		    {
		      //basic track selection
		      //perhaps pixel/tracker hit selection was done at RECO step?
		      if( abs(trkDz1[it]) > 17     || 
			  abs(trkDxy1[it]) > 0.02  || 
			  trkPt[it] < 1            || 
			  trkChi2[it] > 5          ||
			  ipDist2Jet[it] < -0.07   || //shortest distance between track and jet axis
			  ipClosest2Jet[it] > 5.0  || // decay length
			  ipNHitPixel[it] < 2      ||
			  ipNHitStrip[it] < 8
			  ) continue;
		      
		      //"phi matching"
		      nDeltaPhi[counter] = jtphi[j] - trkPhi[it] ;
		      if( fabs(nDeltaPhi[counter]) > pi ) nDeltaPhi[counter] = 2*pi-fabs(nDeltaPhi[counter]) ;
		      else nDeltaPhi[counter] = fabs(nDeltaPhi[counter]) ;
		      
		      //deltaRcut
		      nDeltaEta[counter] = fabs(jteta[j] - trkEta[it]);
		      nDeltaRtrk2Jet[counter]=sqrt(nDeltaPhi[counter]*nDeltaPhi[counter]+nDeltaEta[counter]*nDeltaEta[counter]); 
		      if(nDeltaRtrk2Jet[counter]>0.3)continue;
		      
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
		      nTrkChi2[counter]       = trkChi2[it];
		      nTrkPt[counter]         = trkPt[it];
		      nTrkEta[counter]        = trkEta[it];
		      nTrkPhi[counter]        = trkPhi[it];
		      nTrkDz[counter]         = trkDz1[it];
		      nTrkDxy[counter]        = trkDxy1[it];
		      nIpNHitPixel[counter]   = ipNHitPixel[it];
		      //		      cout << ipNHitPixel[it] << endl;
		      nIpNHitStrip[counter]   = ipNHitStrip[it];
		      //		      cout << ipNHitStrip[it] << endl;
		      //		      nDeltaRtrk2Jet[counter] = deltaRtrk2Jet[it];
		      
		      counter++;
		      
		      //cout << "it = " << it << endl;
		      //cout << "trackPosition = " << trackPosition << endl ;
		      //cout << "counter = " << counter << endl;

		    }//track loop
		}//doTracks check
	      
	      //after we have done the tracks for the jet, NOW organize the tracks into top 3 most significant values
	      if(doMostSignificantTracks)
	      	{											
		  //if(counter>=3)cout<<"doing Most Significant Tracks"<<endl;						
	      	  for(int ij =0;ij<counter;ij++)							
	      	    {											
	      	      //first 2d values								 	
	      	      if (nIP2dsig[ij]>=n1stMost2dSigTrk)					 	
	      	      	{									 	
	      	      	  n3rdMost2dSigTrk = n2ndMost2dSigTrk;					 	
	      	      	  n2ndMost2dSigTrk = n1stMost2dSigTrk;					 	
	      	      	  n1stMost2dSigTrk = nIP2dsig[ij];					 	
	      		  		     								
	      	      	  n3rdIP2dTrk      = n2ndIP2dTrk;						
	      	      	  n2ndIP2dTrk      = n1stIP2dTrk;						
	      	      	  n1stIP2dTrk      = nIP2d[ij];						 	
	      	      	}									 	
	      	      else if (nIP2dsig[ij] <= n1stMost2dSigTrk && nIP2dsig[ij]>=n2ndMost2dSigTrk) 	
	      	      	{									 	
	      	      	  n3rdMost2dSigTrk=n2ndMost2dSigTrk;					 	
	      	      	  n2ndMost2dSigTrk=nIP2dsig[ij];					 	
	      	          										
	      	      	  n3rdIP2dTrk=n2ndIP2dTrk;						 	
	      	      	  n2ndIP2dTrk=nIP2d[ij];						 	
	      	      	}									 	
	      	      else if (nIP2dsig[ij] <= n2ndMost2dSigTrk && nIP2dsig[ij]>=n3rdMost2dSigTrk) 	
	      	      	{									 	
	      	      	  n3rdMost2dSigTrk=nIP2dsig[ij];					 	
	      	          										
	      	      	  n3rdIP2dTrk=nIP2d[ij];						 	
	      	      	}                                                                        	
	      	      											
	      	      //now 3d values									
	      	      if (nIP3dsig[ij]>n1stMost3dSigTrk)					 	
	      	      	{									 	
	      	      	  n3rdMost3dSigTrk=n2ndMost3dSigTrk;					 	
	      	      	  n2ndMost3dSigTrk=n1stMost3dSigTrk;					 	
	      	      	  n1stMost3dSigTrk=nIP3dsig[ij];					 	
	      		  										
	      	      	  n3rdIP3dTrk=n2ndIP3dTrk;						 	
	      	      	  n2ndIP3dTrk=n1stIP3dTrk;						 	
	      	      	  n1stIP3dTrk=nIP3d[ij];						 	
	      	      	}									 	
	      	      else if (nIP3dsig[ij] < n1stMost3dSigTrk && nIP3dsig[ij]>n2ndMost3dSigTrk) 	
	      	      	{									 	
	      	      	  n3rdMost3dSigTrk=n2ndMost3dSigTrk;					 	
	      	      	  n2ndMost3dSigTrk=nIP3dsig[ij];					 	
	      	                                                                                 	
	      	      	  n3rdIP3dTrk=n2ndIP3dTrk;						 	
	      	      	  n2ndIP3dTrk=nIP3d[ij];						 	
	      	      	}									 	
	      	      else if (nIP3dsig[ij] < n2ndMost3dSigTrk && nIP3dsig[ij]>n3rdMost3dSigTrk) 	
	      	      	{									 	
	      	      	  n3rdMost3dSigTrk=nIP3dsig[ij];					 	
	      	                                                                                 	
	      	      	  n3rdIP3dTrk=nIP3d[ij];						 	
	      	      	}   										
	      	      											
	      	    }//doMostSignificantTracks loop							
	      	}//doMostSignificanTracks Check                                                         
	      newTree->Fill();//note this means the event information gets filled in as many times as there are jets in the event to loop over
	                     //may be problematic for vz weighting?
	    }//jetloop
	  vtxTree->Fill();
	}//eventloop
      
      // Cleanup
      if(file_number%1000==0)cout << "closing " << fileName << endl;
      inFile->Close();
      //gROOT->GetListOfFiles()->Remove(inFile);
      fileStream >> fileName;
      
    }
  //clock->Stop();                                      /*debug*/
  //cout << "Real Time: "<<clock->GetRealTime()<<endl;	/*debug*/
  //cout << "CPU Time: "<<clock->GetCPUTime()<<endl;    /*debug*/
  
  // Write to output file
  outFile->cd();
  cout << "writing jet tree to file" << outFileName << endl;
  //kOverrwrite overwrites backup/partially finished tree
  newTree->Write(newTree->GetName(), TObject::kOverwrite);
  cout << "writing vtx tree to file" << outFileName << endl;
  vtxTree->Write(vtxTree->GetName(), TObject::kOverwrite);

  // Cleanup
  cout << "cleaning up..." << endl;
  outFile->Close();
  fileStream.close();
  delete outFile;
  
  return 0;
}

//static basically means set once and then not again
// Return the corresponding weight for an event based on pthat
//static double MCWeights(double MCPthat)
//int MCWeights(double MCPthat)'
int MCCounts(int type)
{
  bool doEventCountFile;
  bool doCJetCountFile ;
  bool doBJetCountFile ;
  bool loopOverFiles;
  string NEventsFile;
  string NbJetsFile;
  string NcJetsFile;
  

  switch(type)
    {
    case 1:  //QCD
      //NEventsFile = weightFilePath + QCD_NEventsFile ; 
      //NbJetsFile  = weightFilePath + QCD_NBJetsFile  ; 
      //NcJetsFile  = weightFilePath + QCD_NCJetsFile  ; 
      NEventsFile = QCD_NEventsFile ; 
      NbJetsFile  = QCD_NBJetsFile  ; 
      NcJetsFile  = QCD_NCJetsFile  ; 
      doEventCountFile = !(std::ifstream( (weightFilePath + NEventsFile).c_str() ) );
      doCJetCountFile  = !(std::ifstream( (weightFilePath + NbJetsFile ).c_str() ) );
      doBJetCountFile  = !(std::ifstream( (weightFilePath + NcJetsFile ).c_str() ) );
      loopOverFiles = doCJetCountFile || doBJetCountFile || doEventCountFile;
      break ;
    case 2:  //B
      NEventsFile =  B_NEventsFile ; 
      NbJetsFile  =  B_NBJetsFile  ; 
      NcJetsFile  = "";
      doEventCountFile = !(std::ifstream( (weightFilePath + NEventsFile).c_str() ) );
      doCJetCountFile  = false;
      doBJetCountFile  = !(std::ifstream( (weightFilePath + NcJetsFile ).c_str() ) );
      loopOverFiles =  doBJetCountFile || doEventCountFile;
      break ;
    case 3:  //C
      NEventsFile =  C_NEventsFile ; 
      NbJetsFile  = "";
      NcJetsFile  =  C_NCJetsFile  ; 
      doEventCountFile = !(std::ifstream( (weightFilePath + NEventsFile).c_str() ) );
      doCJetCountFile  = !(std::ifstream( (weightFilePath + NbJetsFile ).c_str() ) );
      doBJetCountFile  = false;
      loopOverFiles = doCJetCountFile || doEventCountFile;
      break ;
    default:cerr<<"dataType must be 1,2,3 to compute weights"<<endl; return -1;      
    }
  if(!loopOverFiles)
    {
      cout << "All Files Exist! Exiting." << endl;
      return -1;
    }

  double pthatEntries[QCDBins+1];
  double NCJets[QCDBins+1];
  double NBJets[QCDBins+1];

  //initialize the arrays, was getting NaN before. This may be a bit paranoid.
  for (int i=0; i<QCDBins+1; i++) 
    {
      pthatEntries[i]=0;
      NCJets[i]=0;
      NBJets[i]=0;
    }
    
  ifstream inStr(fileList.c_str(), ifstream::in);
  string fileName;
  inStr >> fileName;     
  
  //If you're not at the end of the list, and there's stuff to be calculated, loop over the files
  while ( !inStr.eof()) 
  //for(int kkk = 0;kkk < 50; kkk++)/*debug*/
    {
      //cout <<"test weight computation for 5 files" << endl;/*debug*/
      TFile *inFile=TFile::Open(Form("%s",fileName.c_str()));
      TTree *ch = (TTree *)inFile->Get("akPu3PFJetAnalyzer/t");
      ch->AddFriend("hiEvt=hiEvtAnalyzer/HiTree");
      ch->SetBranchAddress("pthat", &pthat);
      ch->SetBranchAddress("vz", &vz);
      ch->SetBranchAddress("jteta", &jteta);
      ch->SetBranchAddress("refpt", &refpt);
      ch->SetBranchAddress("refparton_flavorForB", &refparton_flavorForB);
      
      TH1D *CJetHist = new TH1D("CJetHist","CJetHist",1,0,10000);
      TH1D *BJetHist = new TH1D("BJetHist","BJetHist",1,0,10000);
      double NEntries = 0;
      double numCJets = 0;
      double numBJets = 0;
      for (int i=0; i<QCDBins+1; i++) 
	{
	  //count events for QCD Weights + HF Weights
	  if(doEventCountFile) 
	    {
	      NEntries = ch->GetEntries( pthatCut[i].c_str() );
	      pthatEntries[i] += NEntries;
	    }
	  //count heavy jets for HF Weights Later	  
	  if(doCJetCountFile)
	    {
	      ch->Draw("jtpt>>CJetHist",Form("%s&&abs(jteta)<2&&abs(vz)<15&&refpt>0&&abs(refparton_flavorForB)==%d",pthatCut[i].c_str(),4),"goff");//draw cjet histo
	      numCJets = (double)CJetHist->Integral();
	      NCJets[i] += numCJets;
	      if (NCJets[i]==NCJets[i]&&NCJets[i]!=NCJets[i]) cout <<"C NaN!!!"<<endl;/*debug*/
	      CJetHist->Reset();
	    }
	  if(doBJetCountFile )
	    {
	      ch->Draw("jtpt>>BJetHist",Form("%s&&abs(jteta)<2&&abs(vz)<15&&refpt>0&&abs(refparton_flavorForB)==%d",pthatCut[i].c_str(),5),"goff");//draw bjet histo
	      numBJets = (double)BJetHist->Integral();
	      NBJets[i] += numBJets;
	      BJetHist->Reset();
	    }
	  
	}//loop over bins
      inStr >> fileName;
      inFile->Close();
    }//loop over files
  
  inStr.close();        

  ofstream EventCountOutput( NEventsFile.c_str() , ofstream::out ); 
  ofstream CJetCountOutput(  NcJetsFile.c_str()  , ofstream::out );
  ofstream BJetCountOutput(  NbJetsFile.c_str()  , ofstream::out );
  
  //write out the information to a text file for later use
  for (int i=0; i<QCDBins+1; i++) 
    {
      cout << "data for " << pthatCut[i] << endl;
      if(doEventCountFile)
	{
	  cout << "pthatEntries[" << i << "] = " << pthatEntries[i] << endl;
	  EventCountOutput << pthatEntries[i] << endl;
	}
      if(doCJetCountFile)
	{
	  cout << "NCJets[" << i << "] = " << NCJets[i] << endl;
	  CJetCountOutput << NCJets[i]<<endl;
	}
      if(doBJetCountFile)
	{
	  cout << "NBJets[" << i << "] = " << NBJets[i] << endl;
	  BJetCountOutput << NBJets[i]<<endl;
	}
    }
  
  //Clean up
  EventCountOutput.close(); 
  CJetCountOutput.close();
  BJetCountOutput.close();
  
  return 0;
}

int NTupleWeights(int type)
{
  string NTupleFile     = "" ;
  string outNTupleFile  = "" ;
  string weightFile     = "" ;
  string outWeightFile  = "" ;
  string QCDHFNJetsFile = "" ;
  string HFNEventsFile  = "" ;
  string HFNJetsFile    = "" ;

  double weights[QCDBins+1]   ;
  double QCDNEvents[QCDBins+1];
  double QCDNHFJets[QCDBins+1];
  double HFNEvents[QCDBins+1] ;
  double NHFJets[QCDBins+1]   ;
  double HFWeights[QCDBins+1] ;
    
  switch(type)
    {
    case 1: 
      NTupleFile    = NTuplePath     + QCDNTuple_noWeights       ; 
      outNTupleFile = QCDNTuple_withWeights       ; 
      weightFile    = weightFilePath + QCDWeightsFile  ; 
      outWeightFile = QCDWeightsFile  ; 
      break;
    case 2: 
      NTupleFile     = NTuplePath     + BJetNTuple_noWeights      ;
      outNTupleFile  = BJetNTuple_withWeights       ;  
      weightFile     = weightFilePath + BJetWeightsFile ; 
      outWeightFile  = BJetWeightsFile ; 
      HFNEventsFile  = weightFilePath + B_NEventsFile   ; 
      HFNJetsFile    = weightFilePath + B_NBJetsFile    ; 
      QCDHFNJetsFile = weightFilePath + QCD_NBJetsFile  ; 
      break;
    case 3: 
      NTupleFile       = NTuplePath     + CJetNTuple_noWeights      ; 
      outNTupleFile    = CJetNTuple_withWeights       ;  
      weightFile       = weightFilePath + CJetWeightsFile ; 
      outWeightFile    = CJetWeightsFile ; 
      HFNEventsFile    = weightFilePath + C_NEventsFile   ; 
      HFNJetsFile      = weightFilePath + C_NCJetsFile    ; 
      QCDHFNJetsFile   = weightFilePath + QCD_NCJetsFile  ; 
      break;
    default:cerr<<"dataType must be 1,2,3 to apply weights (MC Only)"<<endl; return -1;      
    }

  bool allInfo;
  bool weightsExist = (std::ifstream(weightFile.c_str()));

  //compute weights
  if(type==1)//QCD
    {
      allInfo = (std::ifstream(QCD_NEventsFile.c_str())) || weightsExist;
      if(!allInfo){ cout << "not enough info to compute QCD weights" << endl; return -1; }
    
      if(weightsExist)//check for already computed weights.
	{
	  cout<<"QCDWeights exist, using them..."<<endl;
	  ifstream inWeights(weightFile.c_str(),ifstream::in);
	  for(int j=0;j<QCDBins+1;j++) inWeights >> weights[j];
	  inWeights.close();
	}
      else//weights not there, must compute
	{
	  cout<<"QCDWeights do not exist, computing them..."<<endl;
	  ifstream inQCDEvents(QCD_NEventsFile.c_str(),ifstream::in);
	  ofstream outWeights(outWeightFile.c_str(),ofstream::out);
	  for(int j=0;j<QCDBins+1;j++) 
	    {
	      inQCDEvents >> QCDNEvents[j];
	      if(j!=QCDBins)weights[j]=(xsection[j]-xsection[j+1])/(QCDNEvents[j]);
	      else weights[j]=xsection[j]/QCDNEvents[j];
	      outWeights << weights[j] << endl;
	    }							    
	  inQCDEvents.close();
	  outWeights.close();	  
	}
    }
  else//HF
    {
      allInfo = (std::ifstream(QCD_NEventsFile.c_str())&& std::ifstream(QCDHFNJetsFile.c_str()) &&std::ifstream(HFNEventsFile.c_str()) && std::ifstream(HFNJetsFile.c_str())) || weightsExist;
      if(!allInfo){ cout << "not enough info to compute HF weights" << endl; return -1; }
    
      if(weightsExist)//check for already computed weights.
	{
	  cout<<"HFWeights exist, using them..."<<endl;
	  ifstream inWeights(weightFile.c_str(),ifstream::in);
	  for(int j=0;j<QCDBins+1;j++) inWeights >> weights[j];
	  inWeights.close();
	}
      else//weights not there, must compute
	{
	  cout<<"HFWeights do not exist, computing them..."<<endl;
	  ifstream inQCDEvents(QCD_NEventsFile.c_str(),ifstream::in);
	  ifstream inQCDHFJets(QCDHFNJetsFile.c_str(),ifstream::in);
	  ifstream inHFEvents(HFNEventsFile.c_str(),ifstream::in);
	  ifstream inHFJets(HFNJetsFile.c_str(),ifstream::in);
	  
	  ofstream outWeights(outWeightFile.c_str(),ofstream::out);
	  
	  for(int j=0;j<QCDBins+1;j++) 
	    {
	      //grab all details needed for one bin
	      inQCDEvents >> QCDNEvents[j];
	      inQCDHFJets >> QCDNHFJets[j];
	      inHFEvents  >> HFNEvents[j];
	      inHFJets    >> NHFJets[j];
	      //heavy flavor event weight wrt to QCD event
	      HFWeights[j]= (NHFJets[j]/HFNEvents[j])/(QCDNHFJets[j]/QCDNEvents[j]);
	      //overall event weight wrt to distribution
	      if(j!=QCDBins)weights[j]=(xsection[j]-xsection[j+1])/(QCDNEvents[j]+HFWeights[j]*HFNEvents[j]);
	      else weights[j]=xsection[j]/(QCDNEvents[j]+HFWeights[j]*HFNEvents[j]);
	      outWeights << weights[j] << endl;
	    }
	  inQCDEvents.close();
	  inQCDHFJets.close();
	  inHFEvents.close();
	  inHFJets.close();
	  outWeights.close();
	}
    }

  cout<<"the weights are.."<<endl;
  for(int j=0;j<QCDBins+1;j++)cout<<"for " << pthatCut[j] << " weight is "<<weights[j]<<endl; 
  
//  a0->ls();
//  TTree *tinput = (TTree*)a0->Get("tinput");
//  a -> cd();
//  tinput->CloneTree()->Write();

  //Open Ntuple (filled once per jet w/ pthat info) for weighting
  //filled once per event w/ vz info for vz weighting
  //  TFile * NTuple = TFile::Open(Form("%s",NTupleFile.c_str()),"UPDATE");
  TFile * NTuple = TFile::Open(Form("%s",NTupleFile.c_str()),"");
  TTree* nt = (TTree * )NTuple->Get("nt");
  //  nt->SetBranchAddress("pthat",&nPthat);
  TFile * outNTuple = TFile::Open(Form("%s",outNTupleFile.c_str()),"RECREATE");
  outNTuple->cd();
  nt->CloneTree()->Write();
  TTree* newNt = (TTree *)outNTuple->Get("nt");

  //open new weight tree
  TTree* weightTree = new TTree("weightTree","weightTree");
  weightTree->Branch("Weight",&nWeight,"Weight/D");
  newNt->AddFriend(weightTree);
  
  int NEntries = newNt->GetEntries();
  cout << "NEntries = " << NEntries <<endl;
  
  //loop over jets
  for (int i = 0; i<NEntries; i++)
    //for (int i = 0; i<10; i++)
    {
      //grab a jet
      newNt->GetEntry(i);
      //nWeight=2.0;
      //figure out which pthat bin the jet belongs to
      if(i%1000000==0)cout << "weighting jet #" << i << endl;
      int j = 0;
      while(nPthat>pthatBin[j] && j<QCDBins)j++;
      nWeight = weights[j];
      //fill the weight tree once per jet
      weightTree->Fill();	  
    }
  
  //write the new info to the old file, hence "UPDATE" option
  outNTuple->Write();
  return 0;
}

// Create the branches for the new output tree
static inline void newBranches(TTree *newTree) 								       
{													       
  //event specific											       
  //newTree->Branch("weight", &nWeight, "weight/D");  							       
  newTree->Branch("vz", &nVz ,"vz/D");									       
  if(dataType>=1)newTree->Branch("pthat", &nPthat, "pthat/D");
  newTree->Branch("nref"  , &nNref ,"nref/I");								       
  newTree->Branch("evt"   , &Evt  ,"evt/I" );								       
  //  newTree->Branch("lumi"  , &Lumi ,"lumi/I");							       
  newTree->Branch("run"   , &Run  ,"run/I" );		
						       
  //jet variables											       
  newTree->Branch("jtpt", &nJtpt, "jtpt/D");								       
  newTree->Branch("jteta", &nJteta, "jteta/D");								       
  newTree->Branch("jtphi", &nJtphi, "jtphi/D");								       
  newTree->Branch("rawpt", &nRawpt, "rawpt/D");								       
  if(dataType>=1)newTree->Branch("refpt", &nRefpt, "refpt/D");						       
  if(dataType>=1)newTree->Branch("refparton_flavorForB", &nRefparton_flavorForB, "refparton_flavorForB/I");    
  													       
  //muon variables											       
  newTree->Branch("mupt", &nMupt, "mupt/D");								       
  newTree->Branch("muN", &nMuN, "muN/I");								       
  newTree->Branch("mueta", &nMueta, "mueta/D");								       
  newTree->Branch("muphi", &nMuphi, "muphi/D");								       
  newTree->Branch("mudr", &nMudr, "mudr/D");								       
  newTree->Branch("muptrel", &nMuptrel, "muptrel/D");//muon momentum, in the plane transverse to the jet axis  
                                                     //NOT muon transverse momentum projected onto the jet axis
  													       
  //discriminator values										       
  newTree->Branch("discr_ssvHighEff", &nDiscr_ssvHighEff, "discr_ssvHighEff/D");			       
  newTree->Branch("discr_ssvHighPur", &nDiscr_ssvHighPur, "discr_ssvHighPur/D");			       
													       
  //secondary vertex											       
  newTree->Branch("nsvtx", &nNsvtx, "nsvtx/I");								       
  newTree->Branch("svtxntrk", &nSvtxntrk, "svtxntrk/I");						       
  newTree->Branch("svtxdl", &nSvtxdl, "svtxdl/D");							       
  newTree->Branch("svtxdls", &nSvtxdls, "svtxdls/D");							       
  newTree->Branch("svtx2Ddl", &nSvtx2Ddl, "svtx2Ddl/D");						       
  newTree->Branch("svtx2Ddls", &nSvtx2Ddls, "svtx2Ddls/D");						       
  newTree->Branch("svtxm", &nSvtxm, "svtxm/D");								       
  newTree->Branch("svtxpt", &nSvtxpt, "svtxpt/D");							       
  newTree->Branch("svtxeta", &nSvtxEta, "svtxeta/D");							       
  newTree->Branch("svtxphi", &nSvtxPhi, "svtxphi/D");							       
  newTree->Branch("svtxDeltaR2Jet", &nSvtxDeltaR2Jet, "svtxDeltaR2Jet/D");				       
  newTree->Branch("svtxDeltaPhi2Jet", &nSvtxDeltaPhi, "svtxDeltaPhi2Jet/D");				       
  newTree->Branch("svtxDeltaEta2Jet", &nSvtxDeltaEta, "svtxDeltaEta2Jet/D");				       
													       
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
      newTree->Branch("ipNHitPixel" ,&nIpNHitPixel  , "ipNHitPixel[nIP]/I");				       
      newTree->Branch("ipNHitStrip" ,&nIpNHitStrip  , "ipNHitStrip[nIP]/I");				       
													       
      //new track variables from ppTrack tree								       
      newTree->Branch( "trkChi2" , &nTrkChi2, "trkChi2[nIP]/D" );					       
      newTree->Branch( "trkPt"   , &nTrkPt  , "trkPt[nIP]/D"   );					       
      newTree->Branch( "trkEta"  , &nTrkEta , "trkEta[nIP]/D"  );					       
      newTree->Branch( "trkPhi"  , &nTrkPhi , "trkPhi[nIP]/D"  );					       
      newTree->Branch( "trkDz1"  , &nTrkDz  , "trkDz1[nIP]/D"  );					       
      newTree->Branch( "trkDxy1" , &nTrkDz  , "trkDxy1[nIP]/D" );					       
      newTree->Branch( "deltaRtrk2Jet" , &nDeltaRtrk2Jet  , "deltaRtrk2Jet[nIP]/D" );			       
      newTree->Branch( "deltaPhitrk2Jet" , &nDeltaPhi  , "deltaPhitrk2Jet[nIP]/D" );			       
      newTree->Branch( "deltaEtatrk2Jet" , &nDeltaEta  , "deltaEtatrk2Jet[nIP]/D" );			       
													       
      if(doMostSignificantTracks)									       
	{												       
	  newTree->Branch("1stMost2dSigTrk",&n1stMost2dSigTrk ,"1stMost2dSigTrk/D");			       
	  newTree->Branch("2ndMost2dSigTrk",&n2ndMost2dSigTrk ,"2ndMost2dSigTrk/D");			       
	  newTree->Branch("3rdMost2dSigTrk",&n3rdMost2dSigTrk ,"3rdMost2dSigTrk/D");			       
	  newTree->Branch("1stIP2dTrk"     ,&n1stIP2dTrk      ,"1stIP2dTrk/D"     );			       
	  newTree->Branch("2ndIP2dTrk"     ,&n2ndIP2dTrk      ,"2ndIP2dTrk/D"     );			       
	  newTree->Branch("3rdIP2dTrk"     ,&n3rdIP2dTrk      ,"3rdIP2dTrk/D"     );			       
	  newTree->Branch("1stMost3dSigTrk",&n1stMost3dSigTrk ,"1stMost3dSigTrk/D");			       
	  newTree->Branch("2ndMost3dSigTrk",&n2ndMost3dSigTrk ,"2ndMost3dSigTrk/D");			       
	  newTree->Branch("3rdMost3dSigTrk",&n3rdMost3dSigTrk ,"3rdMost3dSigTrk/D");			       
	  newTree->Branch("1stIP3dTrk"     ,&n1stIP3dTrk      ,"1stIP3dTrk/D"     );			       
	  newTree->Branch("2ndIP3dTrk"     ,&n2ndIP3dTrk      ,"2ndIP3dTrk/D"     );			       
	  newTree->Branch("3rdIP3dTrk"     ,&n3rdIP3dTrk      ,"3rdIP3dTrk/D"     );			       
	}												       
  } 													       
  //HLT													       
  newTree->Branch("HLT_PAMu3_v1", &HLT_PAMu3_v1, "HLT_PAMu3_v1/I");					       
  newTree->Branch("HLT_PAMu7_v1", &HLT_PAMu7_v1, "HLT_PAMu7_v1/I");					       
  newTree->Branch("HLT_PAMu12_v1", &HLT_PAMu12_v1, "HLT_PAMu12_v1/I");					       
  newTree->Branch("HLT_PAMu3PFJet20_v1", &HLT_PAMu3PFJet20_v1, "HLT_PAMu3PFJet20_v1/I");		       
  newTree->Branch("HLT_PAMu3PFJet40_v1", &HLT_PAMu3PFJet40_v1, "HLT_PAMu3PFJet40_v1/I");		       
  newTree->Branch("HLT_PAMu7PFJet20_v1", &HLT_PAMu7PFJet20_v1, "HLT_PAMu7PFJet20_v1/I");		       
  newTree->Branch("HLT_PABTagMu_Jet20_Mu4_v1",&HLT_PABTagMu_Jet20_Mu4_v1,"HLT_PABTagMu_Jet20_Mu4_v1/I");                      
													       
  return;												       
}                                                                                                              

// Set all the input branch addresses for the input files
static inline void branchAddresses(TTree *akPu3) 
{
  //EVENT INFO
  akPu3->SetBranchAddress("nref", &nref);
  akPu3->SetBranchAddress("vz", &vz);
  akPu3->SetBranchAddress("evt", &evt);
  //  akPu3->SetBranchAddress("lumi", &lumi);
  akPu3->SetBranchAddress("run", &run);
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
      akPu3->SetBranchAddress("trkDz1" ,&trkDz1 );
      akPu3->SetBranchAddress("trkDxy1",&trkDxy1);
      akPu3->SetBranchAddress("ipNHitPixel",&ipNHitPixel);
      akPu3->SetBranchAddress("ipNHitStrip",&ipNHitStrip);

    }

  //HLT
  akPu3->SetBranchAddress("HLT_PAMu3_v1" , &HLT_PAMu3_v1);
  akPu3->SetBranchAddress("HLT_PAMu7_v1" , &HLT_PAMu7_v1);
  akPu3->SetBranchAddress("HLT_PAMu12_v1", &HLT_PAMu12_v1);
  akPu3->SetBranchAddress("HLT_PAMu3PFJet20_v1", &HLT_PAMu3PFJet20_v1);
  akPu3->SetBranchAddress("HLT_PAMu3PFJet40_v1", &HLT_PAMu3PFJet40_v1);
  akPu3->SetBranchAddress("HLT_PAMu7PFJet20_v1", &HLT_PAMu7PFJet20_v1);
  akPu3->SetBranchAddress("HLT_PABTagMu_Jet20_Mu4_v1", &HLT_PABTagMu_Jet20_Mu4_v1);
    
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
  akPu3->SetBranchAddress("svtx2Ddl" , &svtx2Ddl  );
  akPu3->SetBranchAddress("svtx2Ddls" , &svtx2Ddls  );
  akPu3->SetBranchAddress("svtxm"   , &svtxm    );
  akPu3->SetBranchAddress("svtxpt"  , &svtxpt   );
  akPu3->SetBranchAddress("svtxXPos"   , &svtxXPos    );
  akPu3->SetBranchAddress("svtxYPos"   , &svtxYPos    );
  akPu3->SetBranchAddress("svtxZPos"   , &svtxZPos    );
  return;
}



//a crude macro quickly written to make sure the weights and pthat of the jet are lined up appropriately
int NTupleTest(int type)
{
  int worthless = type;
  string file = NTuplePath + QCDNTuple;
  string theWeights = weightFilePath + QCDWeightsFile ;
  double weightArray[QCDBins+1];

  //read in weights here
  ifstream inWeights(theWeights.c_str(),ifstream::in);
  for(int j=0;j<QCDBins+1;j++) inWeights >> weightArray[j];
  inWeights.close();

  TFile *inFile = TFile::Open(Form("%s",file.c_str()));
  TTree *theTree = (TTree*)inFile->Get("weightTree");
  theTree->AddFriend("nt");
  theTree->SetBranchAddress("Weight",&nWeight);
  theTree->SetBranchAddress("pthat",&nPthat);
  int NEvents = theTree->GetEntries();
  for (int i =0; i<NEvents; i++)
    {
      theTree->GetEvent();//grab pthat and weight
      double actualWeight;
      int j = 0;
      while(nPthat>pthatBin[j] && j>QCDBins)j++;
      actualWeight = weightArray[j];
      if(actualWeight!=nWeight)cout<<"weight is incorrect!" << endl;
      //compare weight the tree gave us to what the weight SHOULD be
    }
  return 0;
}

