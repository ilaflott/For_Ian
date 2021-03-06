/*
Plot maker+formatter. Written for muon-tagged b-jet RpA + RAA
Inherited from: Leo Yu
HeavilyE dited: Ian Laflotte
*/

#include <iostream>
#include <stdio.h>
#include <string.h>

#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TLine.h"
#include "TAxis.h"
#include "TLegend.h"
#include "THStack.h"

void makeQAPlots(const char *, const char *);
//void makeHLTPlots();
static void formatQAPlots(const char *,int);
static void formatLeg(TLegend *);
static void formatHist(TH1 *, const char *, const char *);

//GLOBAL VARIABLES
const double pi = 3.1415926535897932384626433832795028841971693993751058209749445923078164062862089986280348253421170679;
const int n_types = 5; // data, MC, b, c, udsg (0,1,2,3,4...)
const float int_lumi = 4209000000;//inverse millibarns of data. according to lumiCalc2.py, golden lumimask for HLT_PAMu3_v1, 4.209 pb of data.

//file paths
const char *data_file_path = "/net/hisrv0001/home/ilaflott/pp_MC_2760GeV_bTag_forests_ntuples/NTuples/data_leo/";
const char *MC_file_path   = "/net/hisrv0001/home/ilaflott/pp_MC_2760GeV_bTag_forests_ntuples/NTuples/augmented_Samples/";
const char* QCD_file_path  = "/net/hisrv0001/home/ilaflott/pp_MC_2760GeV_bTag_forests_ntuples/NTuples/QCDMC_kurts/";
//const char* QCD_file_path  = "/net/hisrv0001/home/ilaflott/pp_MC_2760GeV_bTag_forests_ntuples/NTuples/HFMC_kurts/";
//const char* QCD_file_path  = "/net/hisrv0001/home/ilaflott/pp_MC_2760GeV_bTag_forests_ntuples/NTuples/HFMC_kurts/";
const char *hist_file_path = "/net/hisrv0001/home/ilaflott/pp_MC_2760GeV_bTag_forests_ntuples/Histograms/6.22.15_muTagbJetRpA_pp_QAplots/";
const char *pdf_file_path  = "/net/hisrv0001/home/ilaflott/pp_MC_2760GeV_bTag_forests_ntuples/Histograms/6.22.15_muTagbJetRpA_pp_QAplots/";

//filenames
const char *data_file_name = "data_6.1.15";
const char *MC_file_name   = "MC_HFaugmented_halfOfficial_6.1.15";
const char* QCD_file_name  = "QCD_6.1.15";


//plot formatting parameters
const int       color[]  = { kBlack, kGray+3, kRed-7, kGreen-6, kBlue-7};
const int    lineColor[] = { kBlack, kWhite, kRed-7, kGreen-6, kBlue-7};
const char *leg_label[]  = { "Data pp", "MC", "b", "c", "udsg" };
const char *leg_label_HLT[]  = { "Data pp (w/ HLT)", "MC (w/ HLT)", "b (w/ HLT)", "c (w/ HLT)", "udsg (w/ HLT)" };

//axis labels
const char *y_label   = "xsec (mb)";
const char *x_label[] =
  {
    "jtpt (GeV)" , "jteta" , "jtphi (rad)" ,         
    "mupt (Gev)" , "muptrel (GeV)" , "mueta" , "muphi (rad)" , "mudr" , 
    "discr_ssvHighEff" , "discr_ssvHighPur" ,
    "nsvtx" , "svtxntrk" , "svtxdl (cm)" , "svtxdls" ,
    "IP2d (cm)" , "IP2sSig" , "IP3d (cm)" , "IP3dSig" , "deltaRtrk2Jet",
    "IP2dSig 1st Trk"   ,  "IP2dSig 2nd Trk"   , "IP2dSig 3rd Trk"   ,
    "IP2d 1st Trk (cm)" ,  "IP2d 2nd Trk (cm)" , "IP2d 3rd Trk (cm)" ,
    "IP3dSig 1st Trk"   ,  "IP3dSig 2nd Trk"   , "IP3dSig 3rd Trk"   ,
    "IP3d 1st Trk (cm)" ,  "IP3d 2nd Trk (cm)" , "IP3d 3rd Trk (cm)" 
  };

//n_vars parameters
const int n_vars_low   = 3  ;//starting variable for formatting
const int n_vars_high  = 3 ;//ending variable for formatting
const int n_vars       = (n_vars_high - n_vars_low) + 1;//for formatting plots, reflects total number of plots being formatted
const int n_vars_TOTAL = 31;//for making plots, always make plots of all variables
//const int n_vars_TOTAL = 1;/*debug*/

//tree branch names
const char *var[] =
  {
    "jtpt",    "jteta",  "jtphi",                                  //jets, vars 0-2			  
    "mupt", "muptrel","mueta",    "muphi",   "mudr",               //muons, vars 3-7			  
    "discr_ssvHighEff", "discr_ssvHighPur",                        //discriminators, 8-9		  
    "nsvtx", "svtxntrk", "svtxdl", "svtxdls",                      //secondary vertex, vars 10-14	  
    "ip2d","ip2dSig","ip3d","ip3dSig","deltaRtrk2Jet",             //impact parameters+deltaR, vars 15-19 
    "1stMost2dSigTrk",    "2ndMost2dSigTrk",    "3rdMost2dSigTrk", //most significant tracks, vars 20-30  
    "1stIP2dTrk"     ,    "2ndIP2dTrk"     ,    "3rdIP2dTrk"     ,
    "1stMost3dSigTrk",    "2ndMost3dSigTrk",    "3rdMost3dSigTrk",
    "1stIP3dTrk"     ,    "2ndIP3dTrk"     ,    "3rdIP3dTrk"          
  };

//cuts
const char *particle_cut[] = 
  { 
    "1",         // Data
    "1",         // Total MC
    "abs(refparton_flavorForB)==5",//b
    "abs(refparton_flavorForB)==4",//c
    "(abs(refparton_flavorForB)==1 || abs(refparton_flavorForB)==2 || abs(refparton_flavorForB)==3 || abs(refparton_flavorForB)==21)"//udsg
  };
      
//plot parameters                                                                                     
const int   nbinsX[] = 
  {   
    25, 15, 15,             //jets, vars 0-2			  
    20, 25, 15, 15, 15,     //muons, vars 3-7			  
    6,   6, 		    //discriminators, 8-9		  
    5,  12,  8, 16, 	    //secondary vertex, vars 10-14	  
    20, 20, 20, 20, 7,	    //impact parameters+deltaR, vars 15-19 
    16, 16, 16,		    //most significant tracks, vars 20-31  
    8, 8, 8,
    16, 16, 16,
    8, 8, 8
  };
const double  lowX[] = 
  { 
      20,  -3,  -pi,            //jets, vars 0-2			  
       0,   0,   -3, -pi, 0,  	//muons, vars 3-7			  
       0,   0,       	  	//discriminators, 8-9		  
       0,   0,    0,   0, 	//secondary vertex, vars 10-14	  
      -0.1, -30, -0.1, -30, 0,  //impact parameters+deltaR, vars 15-19 
      -30,-30,-30,   		//most significant tracks, vars 20-31  
      -0.1,-0.1,-0.1,
      -30,-30,-30,   
      -0.1,-0.1,-0.1 
  };
const double highX[] = 
  {   
    270 ,  3  ,  pi,            //jets, vars 0-2			  
     80 ,  10 ,   3,  pi, 0.5,  //muons, vars 3-7			  
      6 ,   6 ,  		//discriminators, 8-9		  
      5 ,  12 ,   4, 240, 	//secondary vertex, vars 10-14	  
    0.1 ,  30 , 0.1,  30, 0.7,	//impact parameters+deltaR, vars 15-19 
    30  ,  30 , 30 ,   		//most significant tracks, vars 20-31  
    0.1 , 0.1 , 0.1,
    30  ,  30 , 30 ,   
    0.1 , 0.1 , 0.1 
  };
const bool  doLogy[] = 
  { 
    1, 1, 0,                    //jets, vars 0-2			  
    1, 1, 1, 0, 1,  		//muons, vars 3-7			  
    1, 1, 			//discriminators, 8-9		  
    1, 1, 1, 1, 		//secondary vertex, vars 10-14	  
    1,1,1,1, 1,		//impact parameters+deltaR, vars 15-19 
    1, 1, 1,			//most significant tracks, vars 20-31  
    1, 1, 1,
    1, 1, 1,
    1, 1, 1
  };

//const char *default_cut ="vz<15&&vz>-15&&jteta<2&&jteta>-2&&jtpt>40&&HLT_PAMu3_v1&&mupt!=0&&mupt/rawpt<0.95&&svtxdl>0.01&&svtxdl<2.5&&svtxdls>3.0&&svtxm<6.0";
const char *default_cut ="vz<15&&vz>-15&&jteta<2&&jteta>-2&&jtpt>40&&mupt!=0&&mupt/rawpt<0.95&&svtxm<6.0";

//const char *default_version ="HLT_PAMu3_v1_jtpt";
//const char *HLT_bit = "HLT_PAMu3_v1";
//const char *default_version ="HLT_PAMu7_v1_jtpt";
//const char *HLT_bit = "HLT_PAMu7_v1";
//const char *default_version ="HLT_PAMu12_v1_jtpt";
//const char *HLT_bit = "HLT_PAMu12_v1";
const char *default_version ="HLT_PAMu3PFJet40_v1_mupt";
const char *HLT_bit = "HLT_PAMu3PFJet40_v1";

const char *hist_file_name = "";

//MAIN FUNCTIONS
//stackOption == 0 -> Overlaid flavor curves
//stackOption == 1 -> stacked flavor curves
//NOTE fields which arent specified default to the values here. If one only wants to change one or two of the parameters,
//then they must be submitted in order and one after the other. Ex. you want to change the cuts, but because of the order
//of the arguments and the ambiguity inherent, one must also specify a version and an option. Just specifying one string as an input will
//be taken in as an input value to cutsVersion, even if what you really wanted to change is the cuts 
void bTagPlots_IanEdit_forHLTSelections( int option = 0 ,const char* cutsVersion = default_version, const char* cuts = default_cut , int stackOption = 0)
{
  printf("\nYour cuts are:\n %s\n",cuts);
  printf("\nYour version is:\n %s\n",cutsVersion);
  
  char outputFile[1000];
  sprintf(outputFile,"%sbTagPlots_pp_HFaugmented_halfOfficial_%s",hist_file_path,cutsVersion);
  
  switch(stackOption)
    {
    case 0:  printf("\nyou aren't stacking.\n"); break;
    case 1:  printf("\nyou are stacking.\n"); break;
    default: 
      printf("\nnot stacking, not NOT stacking... what are you doing?!\n"); 
      break;
    }
  switch(option)
    {
    case 0:  
      printf("\nMaking+formatting plots.\n");
      makeQAPlots(cuts,(const char*)outputFile);
      formatQAPlots((const char*)outputFile,stackOption); break;
    case 1:
      printf("Only making plots.\n");
      makeQAPlots(cuts,(const char*)outputFile); break;
    case 2:
      printf("Only formatting plots.\n");
      formatQAPlots((const char*)outputFile,stackOption); break;
    default:
      printf("ERROR: 0 for making+formatting plots, 1 for formating only, what are you doing?!\n");
      break;
    }
}

void makeQAPlots(const char* cuts, const char* outputFile)
{
  // Open files and trees
  printf("\nopening data ntuple : %s%s.root\n",data_file_path,data_file_name);
  TFile *data_file = TFile::Open(Form("%s%s.root",data_file_path,data_file_name));
  TTree *data_tree = (TTree *)data_file->Get("nt");
  data_tree->SetMakeClass(1);
  
  printf("opening MC ntuple : %s%s.root\n",MC_file_path,MC_file_name);
  TFile *MC_file = TFile::Open(Form("%s%s.root",MC_file_path,MC_file_name));
  TTree *MC_tree = (TTree *)MC_file->Get("nt");
  MC_tree->SetMakeClass(1);

  TFile *QCD_file = TFile::Open(Form("%s%s.root",QCD_file_path,QCD_file_name));
  TTree *QCD_tree = (TTree *)QCD_file->Get("nt");
  QCD_tree->SetMakeClass(1);

  // Output file
  TFile *out_file = new TFile(Form("%s.root",outputFile), "RECREATE");
  printf("\nthe output file is %s.root \n",outputFile);

  out_file->cd();
  
  // Declare histograms arrays
  TH1D     *hist[n_vars_TOTAL][n_types];      // Data and MC histograms
  TH1D     *hist_HLT[n_vars_TOTAL][n_types];      // Data and MC histograms
  TH1D     *QCDhist_HLT[n_vars_TOTAL][n_types];
  TH1D     *QCDhist[n_vars_TOTAL][n_types];
  //  TH1D     *ratio[n_vars][n_types]; //Data and/or stacked MC divided by HLT selected Data and/or stacked MC
  TH1D     *ratio[n_vars_TOTAL][n_types]; //Data and/or stacked MC divided by HLT selected Data and/or stacked MC
  double    integrals_HLT[n_vars_TOTAL][n_types]; // Integral of hists for stacking
  double    integrals[n_vars_TOTAL][n_types]; // Integral of hists for stacking
  
  double QCD_HFintegral = 0;
  double QCD_HLT_HFintegral = 0;
  int numEntries = 0;
  int numEntries_HLT = 0;

  //  for (int i_var = 0; i_var < n_vars_TOTAL; i_var++)
  for (int i_var = n_vars_low; i_var < n_vars_high+1; i_var++)
    {
      printf("\ni_var:  %d\n\nvariable:  %s\n\n", i_var, var[i_var]);

      //Loop to fill histograms for data and each MC jet flavor
      for (int i_type = 0; i_type < n_types; i_type++)
	{
	  // Initialize histogram
	  hist[i_var][i_type] = new TH1D( Form("hist_%d_%d",i_var,i_type), Form("%s",HLT_bit), nbinsX[i_var], lowX[i_var], highX[i_var]);
	  hist[i_var][i_type]->Sumw2();

	  hist_HLT[i_var][i_type] = new TH1D( Form("hist_HLT_%d_%d",i_var,i_type), Form("%s",HLT_bit), nbinsX[i_var], lowX[i_var], highX[i_var]);
	  hist_HLT[i_var][i_type]->Sumw2();

	  // Fill/draw histogram
	  if (i_type == 0) data_tree->Draw(Form("%s>>hist_%d_%d",var[i_var],i_var,i_type), Form("%s",cuts), "goff");
	  else /*i_type!=0*/ MC_tree->Draw(Form("%s>>hist_%d_%d",var[i_var],i_var,i_type), Form("weight*(%s&&%s)", particle_cut[i_type],cuts), "goff");

	  if (i_type == 0) data_tree->Draw(Form("%s>>hist_HLT_%d_%d",var[i_var],i_var,i_type), Form("%s&&%s",cuts,HLT_bit), "goff");
	  else /*i_type!=0*/ MC_tree->Draw(Form("%s>>hist_HLT_%d_%d",var[i_var],i_var,i_type), Form("weight*(%s&&%s&&%s)", particle_cut[i_type],cuts,HLT_bit), "goff");

	  integrals[i_var][i_type]=hist[i_var][i_type]->Integral();
	  integrals_HLT[i_var][i_type]=hist_HLT[i_var][i_type]->Integral();
	  
	  //Renomalize the B and C contributions, if no enriched files were used, the scale factor is 1.
//	  if(i_type == 2 || i_type == 3)
//	    {
//	      QCDhist[i_var][i_type] = new TH1D(Form("QCDhist_%d_%d",i_var,i_type), Form("QCDhist_%d_%d",i_var,i_type), nbinsX[i_var], lowX[i_var], highX[i_var]);
//	      QCDhist[i_var][i_type]->Sumw2();
//	      QCD_tree->Draw(Form("%s>>QCDhist_%d_%d",var[i_var],i_var,i_type), Form("weight*(%s&&%s)", particle_cut[i_type],cuts), "goff");
//	      QCD_HFintegral = QCDhist[i_var][i_type]->Integral();
//	      delete QCDhist[i_var][i_type];//gotta clean up, lest there exist memory leaks
//	      
//	      QCDhist_HLT[i_var][i_type] = new TH1D(Form("QCDhist_HLT_%d_%d",i_var,i_type), Form("QCDhist_HLT_%d_%d",i_var,i_type), nbinsX[i_var], lowX[i_var], highX[i_var]);
//	      QCDhist_HLT[i_var][i_type]->Sumw2();
//	      QCD_tree->Draw(Form("%s>>QCDhist_HLT_%d_%d",var[i_var],i_var,i_type), Form("weight*(%s&&%s&&%s)", particle_cut[i_type],cuts,HLT_bit), "goff");
//	      QCD_HLT_HFintegral = QCDhist_HLT[i_var][i_type]->Integral();
//	      delete QCDhist_HLT[i_var][i_type];//gotta clean up, lest there exist memory leaks
//
//	      double HFscale_factor = QCD_HFintegral/integrals[i_var][i_type];
//	      printf("QCD HF Integral for i_type == %d is %f\n", i_type, QCD_HFintegral);
//	      printf("HF Scale factor is: %f\n",HFscale_factor);
//	      hist[i_var][i_type]->Scale(HFscale_factor);
//
//	      double HFscale_HLT_factor = QCD_HLT_HFintegral/integrals_HLT[i_var][i_type];
//	      printf("QCD HLT-selected HF Integral for i_type == %d is %f\n", i_type, QCD_HLT_HFintegral);
//	      printf("HLT-selected HF Scale factor is: %f\n",HFscale_HLT_factor);
//	      hist_HLT[i_var][i_type]->Scale(HFscale_HLT_factor);
//	    }

	  if (i_type==0) hist[i_var][i_type]->Scale(1/int_lumi);//no need to scale MC by int_lumi, weighted by cross section
	  integrals[i_var][i_type]=hist[i_var][i_type]->Integral();
	  
	  if (i_type==0) hist_HLT[i_var][i_type]->Scale(1/int_lumi);//no need to scale MC by int_lumi, weighted by cross section
	  integrals_HLT[i_var][i_type]=hist_HLT[i_var][i_type]->Integral();
	}
      
      //second loop to renormalize MC distributions to data
//      double scale_factor = (integrals[i_var][0])/(integrals[i_var][2]+integrals[i_var][3]+integrals[i_var][4]);
//      double HLT_scale_factor = (integrals_HLT[i_var][0])/(integrals_HLT[i_var][2]+integrals_HLT[i_var][3]+integrals_HLT[i_var][4]);
//      for (int i_type = 2; i_type < n_types; i_type++)
//	{
//	  hist[i_var][i_type]->Scale(scale_factor);
//	  integrals[i_var][i_type]=hist[i_var][i_type]->Integral();
//	  numEntries = hist[i_var][i_type]->GetEntries(); 
//	  printf("%s, i_type = %i \nEntries = %i\n",  leg_label[i_type],i_type, numEntries);
//	  printf("Integral = %e\n\n", integrals[i_var][i_type]);
//	}
//      for (int i_type = 2; i_type < n_types; i_type++)
//	{
//	  hist_HLT[i_var][i_type]->Scale(HLT_scale_factor);
//	  integrals_HLT[i_var][i_type]=hist_HLT[i_var][i_type]->Integral();
//	  numEntries_HLT = hist_HLT[i_var][i_type]->GetEntries(); 
//	  printf("%s, i_type = %i \nEntries = %i\n",  leg_label[i_type],i_type, numEntries_HLT);
//	  printf("Integral = %e\n\n", integrals[i_var][i_type]);
//	}
      //we're going to remake the MC hist, so we de-register from current directory
      //if this is not done, we get potential memory leak warninigs
      TObject* old=gDirectory->GetList()->FindObject(Form("hist_%d_1",i_var));
      TObject* old_HLT=gDirectory->GetList()->FindObject(Form("hist_HLT_%d_1",i_var));
      gDirectory->GetList()->Remove(old);
      gDirectory->GetList()->Remove(old_HLT);
      
      //now remake the MC histogram by adding up appropriately scaled histogram types
      TH1D* newMCHist = new TH1D( Form("hist_%d_1",i_var), Form("hist_%d_1",i_var), nbinsX[i_var], lowX[i_var], highX[i_var]);
      newMCHist->Sumw2();
      newMCHist->Add(hist[i_var][2]);newMCHist->Add(hist[i_var][3]);newMCHist->Add(hist[i_var][4]);
      hist[i_var][1]=newMCHist;
      
      TH1D* newMCHist_HLT = new TH1D( Form("hist_HLT_%d_1",i_var), Form("hist_HLT_%d_1",i_var), nbinsX[i_var], lowX[i_var], highX[i_var]);
      newMCHist_HLT->Sumw2();
      newMCHist_HLT->Add(hist_HLT[i_var][2]);newMCHist_HLT->Add(hist_HLT[i_var][3]);newMCHist_HLT->Add(hist_HLT[i_var][4]);
      hist_HLT[i_var][1]=newMCHist_HLT;

      //Ratio plot
      for (int i_type=0; i_type<n_types; i_type++)
	{
	  ratio[i_var][n_types] = new TH1D(Form("ratio_%d_%d",i_var,i_type), Form("HLT Bit Cut/No HLT Bit Cut"), nbinsX[i_var], lowX[i_var], highX[i_var]);
	  ratio[i_var][n_types]->Sumw2();
	  ratio[i_var][n_types]->Divide(hist_HLT[i_var][i_type],hist[i_var][i_type],1,1,"b");
	  ratio[i_var][n_types]->Write();
	  delete ratio[i_var][n_types];
	}
      // Output to .root file
      for (int i_type=0; i_type<n_types; i_type++)
	{
	  hist[i_var][i_type]->Write();
	  delete hist[i_var][i_type];//gotta clean up, lest there exist memory leaks
	  hist_HLT[i_var][i_type]->Write();
	  delete hist_HLT[i_var][i_type];//gotta clean up, lest there exist memory leaks
	}
    }//variable loop

  //clean up
  data_file->Close();
  MC_file->Close();
  out_file->Close();
}

static void formatQAPlots(const char* input_file_name, int stackOption)
{
  //Set histogram style
  //gStyle->SetOptStat(1101);
  //gStyle->SetOptStat("irMe");
  gStyle->SetOptStat(0);
  gROOT->ForceStyle();
  
  // Open file
  //TFile *hist_file = TFile::Open(Form("%s%s_%s.root",pdf_file_path,hist_file_name,version));
  
  printf("opening %s.root\n", input_file_name);
  TFile *hist_file = TFile::Open(Form("%s.root",input_file_name));
  printf("Formatting %s.root\n", input_file_name);
  
  //const char* out_file_name;
  //sprintf(out_file_name, "%s.pdf", input_file_name);
  printf("out_file_name = %s.pdf\n", input_file_name);

  TCanvas *canv[n_vars_TOTAL];
  TLegend *leg[n_vars_TOTAL];
  TH1D    *hist[n_vars_TOTAL][n_types];
  TH1D    *hist_HLT[n_vars_TOTAL][n_types];
  TH1D    *ratio[n_vars_TOTAL][n_types];
  //  THStack *stacked_hist[n_vars];//declare n_vars stacked histograms
  TLine   *one[n_vars_TOTAL];
  
  TCanvas *temp_canv = new TCanvas("temp", "temp", 1200, 600);
  
  temp_canv->Print(Form("%s.pdf(",input_file_name));
  
  for (int i_var=n_vars_low; i_var<n_vars_high+1; i_var++)/*DEBUG*/
    //for (int i_var=0; i_var<n_vars; i_var++)
    {
      printf("Formatting %s\n", var[i_var]);

      canv[i_var] = new TCanvas(Form("canv_%d",i_var),Form("canv_%d",i_var),1200,600);
      canv[i_var]->Divide(2,1);
      
      if (doLogy[i_var]) canv[i_var]->cd(1)->SetLogy();

      //stacked_hist[i_var] = new THStack(Form("stacked_hist_%d",i_var),"Stacked MC");//create the stack
      
      for (int i_type=0; i_type<n_types; i_type++)
	{
	  hist[i_var][i_type] = (TH1D *)hist_file->Get(Form("hist_%d_%d",i_var,i_type));
	  formatHist(hist[i_var][i_type], x_label[i_var], y_label);
	  
	  hist_HLT[i_var][i_type] = (TH1D *)hist_file->Get(Form("hist_HLT_%d_%d",i_var,i_type));
	  formatHist(hist_HLT[i_var][i_type], x_label[i_var], y_label);
	  
	  hist[i_var][i_type]->SetMarkerColor(color[i_type]);
	  if(i_type==0)hist_HLT[i_var][i_type]->SetMarkerColor(kGray);
	  else hist_HLT[i_var][i_type]->SetMarkerColor(color[i_type]+6);
	  
	  hist[i_var][i_type]->SetLineColor(lineColor[i_type]);//if not stacked
	  if(i_type==0)hist_HLT[i_var][i_type]->SetLineColor(kGray+1);//if not stacked
	  else hist_HLT[i_var][i_type]->SetLineColor(lineColor[i_type]+6);//if not stacked
	  
	  //hist[i_var][i_type]->SetLineColor(kBlack);//if stacked
	  //hist[i_var][i_type]->SetLineWidth(1);
	  
	  if (i_type==0 )
	    {
	      //hist[i_var][i_type]->SetLineWidth(2);
	      hist[i_var][i_type]->SetMarkerStyle(8);
	      hist[i_var][i_type]->SetMarkerSize(1.1);

	      hist_HLT[i_var][i_type]->SetMarkerStyle(5);
	      hist_HLT[i_var][i_type]->SetMarkerSize(1.1);

	    }
	  else//i_type>0
	    {
	      stackOption==0;
	      if(stackOption==0)
		{
		}
	      else//stackOption==1
		{
		  hist[i_var][i_type]->SetFillStyle(1001);//solid fill
		  //hist[i_var][i_type]->SetFillColor(kGray+1);
		  hist[i_var][i_type]->SetFillColor(color[i_type]);
		  //hist[i_var][i_type]->SetFillStyle(3004 + i_type % 2); //various light stripes	      
		}
	    }



	  //hist[i_var][i_type]->GetYaxis()->SetMaximum(highY[i_var]);
	  //hist[i_var][i_type]->GetYaxis()->SetRangeUser(1,highY[i_var]);
	  
	  if (!doLogy[i_var]) hist[i_var][i_type]->SetMinimum(0);
	  // if (!doLogy[i_var]) hist_HLT[i_var][i_type]->SetMinimum(0);
	  
	  //if (i_type>=2&&stackOption==1) stacked_hist[i_var]->Add(hist[i_var][i_type]);
	    
	}
      
      canv[i_var]->cd(1);
      
      //trying to get the y-axis ranges right	
      //double yMax = 1.1*hist[i_var][0]->GetBinContent(hist[i_var][0]->GetMaximumBin());
      double yMax = 1.1*hist[i_var][0]->GetBinContent(hist[i_var][0]->GetMaximumBin());
      if(i_var==0) yMax = 1.0*hist[i_var][0]->GetBinContent(hist[i_var][0]->GetMaximumBin());

      if (doLogy[i_var]) hist[i_var][0]->GetYaxis()->SetRange(1/int_lumi, yMax);
      else hist[i_var][0]->GetYaxis()->SetRange(0 , yMax);
      
      printf("%s,i_var = %i, \nmax bin height = %f \n", var[i_var], i_var,  hist[i_var][0]->GetBinContent(hist[i_var][0]->GetMaximumBin()) );

      //hist[i_var][0]->Draw("SCAT E SAME");
      
      if(stackOption==0)
	{
	  hist[i_var][0]->Draw("SCAT E1");
	  hist[i_var][2]->Draw("E SAME");//b
	  hist[i_var][3]->Draw("E SAME");//c
	  hist[i_var][4]->Draw("E SAME");//udsg
	  hist_HLT[i_var][0]->Draw("SCAT SAME E1");
	  hist_HLT[i_var][2]->Draw("E SAME");//b
	  hist_HLT[i_var][3]->Draw("E SAME");//c
	  hist_HLT[i_var][4]->Draw("E SAME");//udsg
	}
      else//stackOption==1
	{
	  hist[i_var][1]->SetStats(0);
	  
	  hist[i_var][1]->Draw("E");//MC Total
	  //stacked_hist[i_var]->Draw("SAME HIST");	  
	  hist[i_var][0]->Draw("SAME SCAT E1");	  
	}

      leg[i_var] = new TLegend(0.77,0.73,0.95,0.93); // (xmin,ymin,xmax,ymax)
      formatLeg(leg[i_var]);

      //add entries to legend
      for (int i_type = 0; i_type < n_types; i_type++)
	{
	  if(i_type==0)//data
	    {    
	      leg[i_var]->AddEntry(hist[i_var][i_type], Form("%s",leg_label[i_type]), "lp");
	    }
	  else//i_type>=1
	    {
	      if(i_type==1)continue;//i.e. total MC
	      leg[i_var]->AddEntry(hist[i_var][i_type], Form("%s",leg_label[i_type]), "lp");
	    }
	}

      for (int i_type = 0; i_type < n_types; i_type++)
	{
	  if(i_type==0)//data
	    {    
	      leg[i_var]->AddEntry(hist_HLT[i_var][i_type], Form("%s",leg_label_HLT[i_type]), "lp");
	    }
	  else//i_type>=1
	    {
	      if(i_type==1)continue;//i.e. total MC
	      leg[i_var]->AddEntry(hist_HLT[i_var][i_type], Form("%s",leg_label_HLT[i_type]), "lp");
	    }
	}



      leg[i_var]->Draw();

      // Ratio plots
      ratio[i_var][0] = (TH1D *)hist_file->Get(Form("ratio_%d_0",i_var));
      
      ratio[i_var][0]->GetYaxis()->CenterTitle();

      formatHist(ratio[i_var][0], var[i_var], "");

      ratio[i_var][0]->SetMarkerColor(color[0]);
      ratio[i_var][0]->SetMarkerStyle(8);
      ratio[i_var][0]->SetMarkerSize(0.8);
      ratio[i_var][0]->SetMinimum(0.5);
      ratio[i_var][0]->SetMaximum(1.5);

      canv[i_var]->cd(2);
      
      //ratio[i_var]->Draw("E SAME"); // This has to come before the TLine for some reason...
      ratio[i_var][0]->Draw();
      
      one[i_var] = new TLine(lowX[i_var],1,highX[i_var],1);
      one[i_var]->SetLineColor(color[1]);
      one[i_var]->Draw("SAME");
      
      for(int i_type=2;i_type<n_types;i_type++)
	{
	  ratio[i_var][i_type] = (TH1D *)hist_file->Get(Form("ratio_%d_%d",i_var,i_type));
	  ratio[i_var][i_type]->SetLineColor(color[i_type]);
	  ratio[i_var][i_type]->Draw("SAME"); // Plot it again, over the line
	}
      //output pdf	
      canv[i_var]->Print(Form("%s.pdf",input_file_name));
      printf("i_var = %i and input_file_name=%s\n",i_var,input_file_name);
      
    }
  
  temp_canv->Print(Form("%s.pdf)",input_file_name));
  hist_file->Close();
}









static void formatLeg(TLegend *l)
{
    l->SetFillColor(0);
    l->SetTextSize(0.03);
    l->SetBorderSize(0.01);
    l->SetFillStyle(0);
}

static void formatHist(TH1 *h, const char *x_label1, const char *y_label1)
{
  //h->Sumw2();

  h->GetXaxis()->SetTitle(Form("%s",x_label1));
  //h->GetXaxis()->CenterTitle();

  h->GetYaxis()->SetTitle(Form("%s",y_label1));
  h->GetYaxis()->CenterTitle();
}
