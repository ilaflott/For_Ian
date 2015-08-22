/*
Plot maker+formatter. Written for muon-tagged b-jet RpA + RAA
Inherited from: Leo Yu
Heavily Edited: Ian Laflotte
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

void makeQAPlots(string, string);
//void makeHLTPlots();
static void formatQAPlots(string,int);
static void formatLeg(TLegend *);
static void formatHist(TH1 *, string, string);

//GLOBAL VARIABLES
const double pi = 3.1415926535897932384626433832795028841971693993751058209749445923078164062862089986280348253421170679;
const float int_lumi = 4209000000;//inverse millibarns of data. according to lumiCalc2.py, golden lumimask for HLT_PAMu3_v1, 4.209 pb of data.

//file paths
const string hist_file_path   = "/net/hisrv0001/home/ilaflott/Leos_Analysis/CMSSW_5_3_20_FOREST_PLOTS/src/For_Ian/5_Plot_Variables/Histograms/";
const string pdf_file_path    = "/net/hisrv0001/home/ilaflott/Leos_Analysis/CMSSW_5_3_20_FOREST_PLOTS/src/For_Ian/5_Plot_Variables/Histograms/";
const string NTuple_file_path = "/net/hisrv0001/home/ilaflott/Leos_Analysis/CMSSW_5_3_20_FOREST_PLOTS/src/For_Ian/4_Create_NTuples/good_NTuples/";
const string data_file_path   = NTuple_file_path;
const string MC_file_path     = NTuple_file_path;
const string QCD_file_path    = NTuple_file_path;
 
//filenames
const string data_file_name = "data_NTuple_7.23.15.root";
const string MC_file_name   = "TotalMCNTuple_WithWeights.root";
const string QCD_file_name  = "QCDJets_NTuple_8.3.15_withWeights.root";
//const string hist_file_name = "hist_Test";
//const string pdf_file_name = hist_file_name;

//cuts and naming, defualt values
//const string default_cut ="vz<15&&vz>-15&&jteta<2&&jteta>-2&&jtpt>40&&HLT_PAMu3_v1&&mupt!=0&&mupt/rawpt<0.95&&svtxdl>0.01&&svtxdl<2.5&&svtxdls>3.0&&svtxm<6.0";
//vz cut, jteta cut, jtpt cut, HLT cut, mupt!=0 cut, W cut for data
const string default_cut ="jteta<2&&jteta>-2&&jtpt>20";
const string default_version ="AllPlots_jteta2_jtpt20";

//n_vars parameters
const int n_types = 5; // data, MC, b, c, udsg (0,1,2,3,4...)
const int n_vars_TOTAL = 33;//for making plots, always make plots of all variables
const int n_vars_low   = 0  ;//starting variable
//const int n_vars_high  = 2  ;//ending variable 
const int n_vars_high  = n_vars_TOTAL - 1  ;//for all plots
const int n_vars       = (n_vars_high - n_vars_low) + 1;//reflects total number of plots being made/formattex


//plot formatting parameters
const int       color[]  = { kBlack, kGray+3, kRed-7, kGreen-6, kBlue-7};
const int    lineColor[] = { kBlack, kWhite, kRed-7, kGreen-6, kBlue-7};
const string leg_label[]  = { "Data pp", "MC", "b", "c", "udsg" };

//axis labels
const string y_label   = "xsec (mb)";
const string x_label[] =
  {
    "jtpt (GeV)" , "jteta" , "jtphi (rad)" ,         
    "mupt (Gev)" , "muptrel (GeV)" , "mueta" , "muphi (rad)" , "mudr" , 
    "discr_ssvHighEff" , "discr_ssvHighPur" ,
    "nsvtx" , "svtxntrk" , "svtx2Ddl (cm)" , "svtx2Ddls" , "svtxm (GeV)","svtxDeltaR2Jet",
    "IP2d (cm)" , "IP2sSig" , "IP3d (cm)" , "IP3dSig" , "deltaRtrk2Jet",
    "IP2dSig 1st Trk"   ,  "IP2dSig 2nd Trk"   , "IP2dSig 3rd Trk"   ,
    "IP2d 1st Trk (cm)" ,  "IP2d 2nd Trk (cm)" , "IP2d 3rd Trk (cm)" ,
    "IP3dSig 1st Trk"   ,  "IP3dSig 2nd Trk"   , "IP3dSig 3rd Trk"   ,
    "IP3d 1st Trk (cm)" ,  "IP3d 2nd Trk (cm)" , "IP3d 3rd Trk (cm)" 
  };

//tree branch names
const string var[] =
  {
    "jtpt",    "jteta",  "jtphi",                                                               //jets, vars 0-2		       	  
    "mupt", "muptrel","mueta",    "muphi",   "mudr",                                            //muons, vars 3-7		       	  
    "discr_ssvHighEff", "discr_ssvHighPur",                                                     //discriminators, 8-9		       
    "nsvtx", "svtxntrk", "svtx2Ddl", "svtx2Ddls", "svtxm","svtxDeltaR2Jet",                     //secondary vertex, vars 12-16	       
    "ip2d","ip2dSig","ip3d","ip3dSig","deltaRtrk2Jet",                                          //impact parameters+deltaR, vars 17-21 
    "1stMost2dSigTrk",    "2ndMost2dSigTrk",    "3rdMost2dSigTrk",                              //most significant tracks, vars 22-32  
    "1stIP2dTrk"     ,    "2ndIP2dTrk"     ,    "3rdIP2dTrk"     ,
    "1stMost3dSigTrk",    "2ndMost3dSigTrk",    "3rdMost3dSigTrk",
    "1stIP3dTrk"     ,    "2ndIP3dTrk"     ,    "3rdIP3dTrk"          
  };

//cuts
const string particle_cut[] = 
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
    54, 15, 15,              //jets, vars 0-2		       
    20, 25, 15, 15, 15,      //muons, vars 3-7		       
    12,   12, 		     //discriminators, 8-9		       
    5,  12,  8, 16, 14,10,    //secondary vertex, vars 12-16	       	  
    20, 20, 20, 20, 7,	     //impact parameters+deltaR, vars 17-21 
    16, 16, 16,		     //most significant tracks, vars 22-32  
    8, 8, 8,
    16, 16, 16,
    8, 8, 8
  };
const double  lowX[] = 
  { 
      0,  -3,  -pi,              //jets, vars 0-2		       
       0,   0,   -3, -pi, 0,  	 //muons, vars 3-7		        
       0,   0,       	      	 //discriminators, 8-9		       
      0,   0,    0,   0, 0,0, 	 //secondary vertex, vars 12-16	       
      -0.1, -30, -0.1, -30, 0, 	 //impact parameters+deltaR, vars 17-21 
      -30,-30,-30,   	       	 //most significant tracks, vars 22-32  
      -0.1,-0.1,-0.1,
      -30,-30,-30,   
      -0.1,-0.1,-0.1 
  };
const double highX[] = 
  {   
    270 ,  3  ,  pi,           	   //jets, vars 0-2		       
     80 ,  10 ,   3,  pi, 0.5, 	   //muons, vars 3-7		       
      6 ,   6 ,  	       	   //discriminators, 8-9		       
    5 ,  12 ,   4, 240, 7,5,  	   //secondary vertex, vars 12-16	       
    0.1 ,  30 , 0.1,  30, 0.7, 	   //impact parameters+deltaR, vars 17-21 
    30  ,  30 , 30 ,   	       	   //most significant tracks, vars 22-32  
    0.1 , 0.1 , 0.1,
    30  ,  30 , 30 ,   
    0.1 , 0.1 , 0.1 
  };
const bool  doLogy[] = 
  { 
    1, 1, 0,               //jets, vars 0-2		       
    1, 1, 1, 0, 1,  	   //muons, vars 3-7		       
    1, 1, 		   //discriminators, 8-9		       
    1, 1, 1, 1,1, 1, 	   //secondary vertex, vars 12-16	       
    1,1,1,1, 1,		   //impact parameters+deltaR, vars 17-21 
    1, 1, 1,		   //most significant tracks, vars 22-32  
    1, 1, 1,
    1, 1, 1,
    1, 1, 1
  };


//NOTE fields which arent specified default to the values here. If one only wants to change one or two of the parameters,
//then they must be submitted in order and one after the other. Ex. you want to change the cuts, but because of the order
//of the arguments and the ambiguity inherent, one must also specify a version and an option. Just specifying one string as an input will
//be taken in as an input value to cutsVersion, even if what you really wanted to change is the cuts 
int bTagPlots( int option = 0 , int stackOption = 1 , string cutsVersion = default_version , string cuts = default_cut  )
{
  string finalCuts = "";
  string finalVersion="";

  if(cutsVersion==default_version) finalVersion = default_version;
  else finalVersion = default_version + "_" + cutsVersion;

  const string outputFile = hist_file_path + "bTagPlots_pp_" + finalVersion ;//
  
  if(cuts==default_cut) finalCuts = default_cut ;
  else finalCuts = default_cut + "&&" + cuts;

  cout << "your final version is " << finalVersion << endl;
  cout << "your final cuts are " << finalCuts << endl;

  if (stackOption==0)       cout << "you aren't stacking" << endl ; 
  else if (stackOption==1)  cout << "you are stacking" << endl;
  else                      cout << "not stacking, not NOT stacking... what are you doing?!" <<endl;
 
  switch(option)
    {
    case 0:  
      cout << "Making plots..." << endl;
      makeQAPlots(finalCuts,outputFile);
      cout << "Formatting plots..." << endl;
      formatQAPlots(outputFile,stackOption); break;
    case 1:
      cout << "Only making plots..." << endl;
      makeQAPlots(finalCuts,outputFile); break;
    case 2:
      cout << "Only formatting plots..." << endl;;
      formatQAPlots(outputFile,stackOption); break;
    default:
      cout << "ERROR: 0 for making+formatting plots, 1 for formating only, what are you doing?!" << endl;
      return -1;
    }
  return 0;
}

void makeQAPlots(string cuts, string outputFileName)
{
  // Open files and trees
  cout << "opening data ntuple : " << data_file_path << data_file_name << endl;
  
  string theDataFile = data_file_path + data_file_name;
  TFile *data_file = TFile::Open(Form("%s",theDataFile.c_str()));
  TTree *data_tree = (TTree *)data_file->Get("nt");
  data_tree->SetMakeClass(1);
  
  cout << "opening MC ntuple : " << MC_file_path << MC_file_name << endl;
  
  string theMCFile = MC_file_path + MC_file_name;
  TFile *MC_file = TFile::Open(Form("%s",theMCFile.c_str()));
  TTree *MC_tree = (TTree *)MC_file->Get("nt");
  MC_tree->AddFriend("weightTree");
  MC_tree->SetMakeClass(1);

  cout << "opening QCD ntuple : " << QCD_file_path << QCD_file_name << endl;

  string theQCDFile=QCD_file_path+QCD_file_name;
  TFile *QCD_file = TFile::Open(Form("%s",theQCDFile.c_str()));
  TTree *QCD_tree = (TTree *)QCD_file->Get("nt");
  QCD_tree->AddFriend("weightTree");  
  QCD_tree->SetMakeClass(1);

  // Output file
  string theRootFile = outputFileName + ".root";
  cout << "Opening .root file for output: " << theRootFile << endl; 
  TFile *out_file = new TFile(Form("%s",theRootFile.c_str()), "RECREATE");
  out_file->cd();
  
  // Declare histograms arrays
  TH1D     *hist[n_vars][n_types];      // Data and MC histograms
  TH1D     *QCDhist[n_vars][n_types];
  TH1D     *ratio[n_vars];              // Data/MC ratio plots
  double    integrals[n_vars][n_types]; // Integral of hists for stacking
  
  double QCD_HFintegral = 0;
  int numEntries = 0;
  
  cout <<"creating histogram for each variable in range.."<< n_vars_low << " to " << n_vars_high << endl;
  for (int i_var = 0; i_var < n_vars; i_var++)
    {
      int variableIndex = i_var + n_vars_low;
      cout << "variable index:  " << variableIndex << endl;
      cout << "variable: " << var[variableIndex] << endl;

      //Loop to fill histograms for data and each MC jet flavor
      for (int i_type = 0; i_type < n_types; i_type++)
	{
	  // Initialize histogram
	  hist[i_var][i_type] = new TH1D( Form("hist_%d_%d",variableIndex,i_type), Form("hist_%d_%d",variableIndex,i_type), nbinsX[variableIndex], lowX[variableIndex], highX[variableIndex]);
	  hist[i_var][i_type]->Sumw2();

	  // Fill/draw histogram
	  if (i_type == 0) data_tree->Draw(Form("%s>>hist_%d_%d",var[variableIndex].c_str(),variableIndex,i_type), Form("%s&&mupt/rawpt<0.95",cuts.c_str()), "goff");
	  else /*i_type!=0*/ MC_tree->Draw(Form("%s>>hist_%d_%d",var[variableIndex].c_str(),variableIndex,i_type), Form("Weight*(%s&&%s)", particle_cut[i_type].c_str(),cuts.c_str()), "goff");
	    
	  integrals[i_var][i_type]=hist[i_var][i_type]->Integral();
	  
	  //Renomalize the B and C contributions, if no enriched files were used, the scale factor is 1.
	  if(i_type == 2 || i_type == 3)
	    {
	      cout << "renormalizing B and/or C to QCD distribution" << endl;
	      QCDhist[i_var][i_type] = new TH1D(Form("QCDhist_%d_%d",variableIndex,i_type), Form("QCDhist_%d_%d",variableIndex,i_type), nbinsX[variableIndex], lowX[variableIndex], highX[variableIndex]);
	      QCDhist[i_var][i_type]->Sumw2();
	      QCD_tree->Draw(Form("%s>>QCDhist_%d_%d",var[variableIndex].c_str(),variableIndex,i_type), Form("Weight*(%s&&%s)", particle_cut[i_type].c_str(),cuts.c_str()), "goff");
	      QCD_HFintegral = QCDhist[i_var][i_type]->Integral();
	      delete QCDhist[i_var][i_type];//gotta clean up, lest there exist memory leaks

	      double HFscale_factor = QCD_HFintegral/integrals[i_var][i_type];
	      cout << "QCD HF Integral for i_type == "<< i_type<< " is "  << QCD_HFintegral << endl;	      cout << "HF Scale factor is: " << HFscale_factor << endl;
	      
	      hist[i_var][i_type]->Scale(HFscale_factor);
	    }
	  if (i_type==0) hist[i_var][i_type]->Scale(1/int_lumi);//no need to scale MC by int_lumi, weighted by cross section
	  integrals[i_var][i_type]=hist[i_var][i_type]->Integral();
	  if(i_type==4)cout << "type loop Integral = " <<  integrals[i_var][i_type] << endl;
	}//type loop for filling
      
      //second loop to renormalize MC distributions to data
      double scale_factor = (integrals[i_var][0])/(integrals[i_var][2]+integrals[i_var][3]+integrals[i_var][4]);
      cout << "renormalizing MC distributions to data" << endl;
      for (int i_type = 2; i_type < n_types; i_type++)
	{
	  if(i_type==4)
	    {
	      cout << "i_type = 4" << endl;
	      cout << "scale_factor = " << scale_factor << endl;
	    }
	  hist[i_var][i_type]->Scale(scale_factor);
	  integrals[i_var][i_type]=hist[i_var][i_type]->Integral();
	  numEntries = hist[i_var][i_type]->GetEntries(); 
	  cout << "i_type = " << i_type << endl; 
	  if(i_type==4)cout << "Entries = " << numEntries << endl;
	  if(i_type==4)cout << "Integral = " <<  integrals[i_var][i_type] << endl;
	}

      //we're going to remake the MC hist, so we de-register from current directory
      //if this is not done, we get potential memory leak warninigs
      TObject* old=gDirectory->GetList()->FindObject(Form("hist_%d_1",variableIndex));
      gDirectory->GetList()->Remove(old);

      //now remake the MC histogram by adding up appropriately scaled histogram types
      TH1D* newMCHist = new TH1D( Form("hist_%d_1",variableIndex), Form("hist_%d_1",variableIndex), nbinsX[variableIndex], lowX[variableIndex], highX[variableIndex]);
      newMCHist->Sumw2();
      newMCHist->Add(hist[i_var][2]);
      newMCHist->Add(hist[i_var][3]);
      newMCHist->Add(hist[i_var][4]);
      
      hist[i_var][1]=newMCHist;
      
      //Ratio plot
      ratio[i_var] = new TH1D(Form("ratio_%d",variableIndex), Form("data/MC"), nbinsX[variableIndex], lowX[variableIndex], highX[variableIndex]);
      ratio[i_var]->Sumw2();
      ratio[i_var]->Divide(hist[i_var][0],hist[i_var][1],1,1,"b");
      
      // Output to .root file
      for (int i_type=0; i_type<n_types; i_type++)
	{
	  hist[i_var][i_type]->Write();
	  delete hist[i_var][i_type];//gotta clean up, lest there exist memory leaks
	}
      
      ratio[i_var]->Write();
      delete ratio[i_var];//gotta clean up, lest there exist memory leaks
    }//variable loop
  
  //clean up
  data_file->Close();
  MC_file->Close();
  out_file->Close();
}

static void formatQAPlots(string inputRootFileName, int stackOption)
{
  //Set histogram style
  //gStyle->SetOptStat(1101);
  gStyle->SetOptStat("irMe");
  gROOT->ForceStyle();
  
  string theRootFile = inputRootFileName + ".root";
  cout << "opening + formatting .root file: " <<  theRootFile << endl;
  TFile *hist_file = TFile::Open(Form("%s",theRootFile.c_str()));
  
  string thePDFFile = inputRootFileName + ".pdf";
  //cout << "writing out to .pdf file:  " <<  thePDFFile << endl;

  TCanvas *canv[n_vars];
  TLegend *leg[n_vars];
  TH1D    *hist[n_vars][n_types];
  TH1D    *ratio[n_vars];
  THStack *stacked_hist[n_vars];//declare n_vars stacked histograms
  TLine   *one[n_vars];
  
  TCanvas *temp_canv = new TCanvas("temp", "temp", 1200, 600);
  
  temp_canv->Print(Form("%s(",thePDFFile.c_str()));
  
  //for (int i_var=n_vars_low; i_var<n_vars_high+1; i_var++)/*DEBUG*/
  for (int i_var=0; i_var<n_vars; i_var++)
    {
      int variableIndex = i_var + n_vars_low;
      cout << "Formatting " <<  var[variableIndex] << endl;
      
      canv[i_var] = new TCanvas(Form("canv_%d",i_var),Form("canv_%d",i_var),1200,600);
      canv[i_var]->Divide(2,1);
      
      if (doLogy[variableIndex]) canv[i_var]->cd(1)->SetLogy();

      stacked_hist[i_var] = new THStack(Form("stacked_hist_%d",i_var),"Stacked MC");//create the stack
      

      for (int i_type=0; i_type<n_types; i_type++)
	{
	  hist[i_var][i_type] = (TH1D *)hist_file->Get(Form("hist_%d_%d",variableIndex,i_type));
	  formatHist(hist[i_var][i_type], x_label[variableIndex], y_label);
	  
	  hist[i_var][i_type]->SetMarkerColor(color[i_type]);
	  //hist[i_var][i_type]->SetLineColor(lineColor[i_type]);//if not stacked
	  hist[i_var][i_type]->SetLineColor(kBlack);//if stacked
	  //hist[i_var][i_type]->SetLineWidth(1);
	  
	  if (i_type==0 )
	    {
	      //hist[i_var][i_type]->SetLineWidth(2);
	      hist[i_var][i_type]->SetMarkerStyle(8);
	      hist[i_var][i_type]->SetMarkerSize(1.1);
	    }
	  else//i_type>0
	    {
	      if(stackOption==0) continue;
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
	  
	  if (!doLogy[variableIndex]) hist[i_var][i_type]->SetMinimum(0);
	  
	  if (i_type>=2&&stackOption==1) stacked_hist[i_var]->Add(hist[i_var][i_type]);
	    
	}
      
      canv[i_var]->cd(1);
      
      //trying to get the y-axis ranges right	
      //double yMax = 1.1*hist[i_var][0]->GetBinContent(hist[i_var][0]->GetMaximumBin());
      double yMax = 1.1*hist[i_var][0]->GetBinContent(hist[i_var][0]->GetMaximumBin());
      if(variableIndex==0) yMax = 1.0*hist[i_var][0]->GetBinContent(hist[i_var][0]->GetMaximumBin());

      if (doLogy[variableIndex]) hist[i_var][0]->GetYaxis()->SetRange(1/int_lumi, yMax);
      else hist[i_var][0]->GetYaxis()->SetRange(0 , yMax);
      
      cout << var[variableIndex] << ", variableIndex = " << variableIndex << endl;
      cout << "max bin height = " << hist[i_var][0]->GetBinContent( hist[i_var][0]->GetMaximumBin())  << endl;

      //hist[i_var][0]->Draw("SCAT E SAME");
      
      if(stackOption==0)
	{
	  hist[i_var][0]->Draw("SCAT E1");
	  hist[i_var][2]->Draw("HIST E SAME");//b
	  hist[i_var][3]->Draw("HIST E SAME");//c
	  hist[i_var][4]->Draw("HIST E SAME");//udsg
	}
      else//stackOption==1
	{
	  hist[i_var][1]->SetStats(0);
	  
	  hist[i_var][1]->Draw("E");//MC Total
	  stacked_hist[i_var]->Draw("SAME HIST");	  
	  hist[i_var][0]->Draw("SAME SCAT E1");	  
	}

      leg[i_var] = new TLegend(0.8,0.78,0.95,0.93); // (xmin,ymin,xmax,ymax)
      formatLeg(leg[i_var]);

      //add entries to legend
      for (int i_type = 0; i_type < n_types; i_type++)
	{
	  if(i_type==0)//data
	    {    
	      leg[i_var]->AddEntry(hist[i_var][i_type], Form("%s",leg_label[i_type].c_str()), "lp");
	    }
	  else//i_type>=1
	    {
	      if(i_type==1)continue;//i.e. total MC
	      leg[i_var]->AddEntry(hist[i_var][i_type], Form("%s",leg_label[i_type].c_str()), "f");
	    }
	}

      leg[i_var]->Draw();

      // Ratio plots
      ratio[i_var] = (TH1D *)hist_file->Get(Form("ratio_%d",i_var));
      
      ratio[i_var]->GetYaxis()->CenterTitle();

      formatHist(ratio[i_var], var[i_var], "data/MC");

      ratio[i_var]->SetMarkerColor(color[0]);
      ratio[i_var]->SetMarkerStyle(8);
      ratio[i_var]->SetMarkerSize(0.8);
      ratio[i_var]->SetMinimum(0.5);
      ratio[i_var]->SetMaximum(1.5);

      canv[i_var]->cd(2);
      
      //ratio[i_var]->Draw("E SAME"); // This has to come before the TLine for some reason...
      ratio[i_var]->Draw();
      
      one[i_var] = new TLine(lowX[variableIndex],1,highX[variableIndex],1);
      one[i_var]->SetLineColor(color[1]);
      one[i_var]->Draw("SAME");
      
      ratio[i_var]->SetLineColor(color[0]);
      ratio[i_var]->Draw("SAME"); // Plot it again, over the line
      
      //output pdf	
      canv[i_var]->Print(Form("%s",thePDFFile.c_str()));
      //cout << "i_var = " << i_var << " and input_file_name=" << input_file_name << endl;
      
    }
  cout << "writing out to .pdf file:  " <<  thePDFFile << endl;
  temp_canv->Print(Form("%s]",thePDFFile.c_str()));
  hist_file->Close();
}

////////////////////
//HELPER FUNCTIONS//
////////////////////

static void formatLeg(TLegend *l)
{
    l->SetFillColor(0);
    l->SetTextSize(0.03);
    l->SetBorderSize(0.01);
    l->SetFillStyle(0);
}

static void formatHist(TH1 *h, string xLabel, string yLabel)
{
  //h->Sumw2();

  h->GetXaxis()->SetTitle(Form("%s",xLabel.c_str()));
  //h->GetXaxis()->CenterTitle();

  h->GetYaxis()->SetTitle(Form("%s",yLabel.c_str()));
  //h->GetYaxis()->CenterTitle();
}
