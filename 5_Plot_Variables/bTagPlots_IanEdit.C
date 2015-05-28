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

void makePlots(const char *, const char *);

static void formatPlots(const char *,int);
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
const char *hist_file_path = "/net/hisrv0001/home/ilaflott/pp_MC_2760GeV_bTag_forests_ntuples/Histograms/5.22.15_muTagbJetRpA_pp_QAplots/";
const char *pdf_file_path  = "/net/hisrv0001/home/ilaflott/pp_MC_2760GeV_bTag_forests_ntuples/Histograms/5.22.15_muTagbJetRpA_pp_QAplots/";

//filenames
const char *data_file_name = "data_5.22.15";
const char *MC_file_name   = "MC_HFaugmented_halfOfficial_5.22.15";
const char* QCD_file_name  = "QCD_5.22.15";
const char *hist_file_name = "bTagPlots_pp";

//plot formatting parameters
const int       color[]  = { kBlack, kGray+3, kRed-7, kGreen-6, kBlue-7};
const int    lineColor[] = { kBlack, kWhite, kRed-7, kGreen-6, kBlue-7};
const char *leg_label[]  = { "Data pp", "MC", "b", "c", "udsg" };

//axis labels
const char *y_label   = "xsec (mb)";
const char *x_label[] =
  {
    "jtpt (GeV)",    "jteta",  "jtphi (rad)",         
    "mupt (Gev)", "muptrel (GeV)",    "muphi (rad)",   "mudr", "mueta",
    "discr_ssvHighEff", "discr_ssvHighPur",
    "nsvtx", "svtxntrk", "svtxdl (cm)",  "svtxdls",
    "ip2d (cm)","2dIP Sig","ip3d (cm)","3D IP Sig","deltaRtrk2Jet",
    "IP2dSig 1st Trk"  ,  "IP2dSig 2nd Trk"  , "IP2dSig 3rd Trk"  ,
    "IP2d 1st Trk (cm)",  "IP2d 2nd Trk (cm)", "IP2d 3rd Trk (cm)",
    "IP3dSig 1st Trk"  ,  "IP3dSig 2nd Trk"  , "IP3dSig 3rd Trk"  ,
    "IP3d 1st Trk (cm)",  "IP3d 2nd Trk (cm)", "IP3d 3rd Trk (cm)"
  };

//n_vars parameters
const int n_vars_low   = 15 ;//starting variable for formatting
const int n_vars_high  = 30 ;//ending variable for formatting
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
    0, 0, 0, 0, 1,		//impact parameters+deltaR, vars 15-19 
    1, 1, 1,			//most significant tracks, vars 20-31  
    1, 1, 1,
    1, 1, 1,
    1, 1, 1
  };

const char *default_cut ="vz<15&&vz>-15&&jteta<2&&jteta>-2&&jtpt>40&&HLT_PAMu3_v1&&mupt!=0&&mupt/rawpt<0.95&&svtxdl>0.01&&svtxdl<2.5&&svtxdls>3.0&&svtxm<6.5";
const char *default_version ="default_Cuts";

//MAIN FUNCTIONS
//stackOption == 0 -> Overlaid flavor curves
//stackOption == 1 -> stacked flavor curves
//NOTE fields which arent specified default to the values here. If one only wants to change one or two of the parameters,
//then they must be submitted in order and one after the other. Ex. you want to change the cuts, but because of the order
//of the arguments and the ambiguity inherent, one must also specify a version and an option. Just specifying one string as an input will
//be taken in as an input value to cutsVersion, even if what you really wanted to change is the cuts 
void bTagPlots_IanEdit( int option = 1 ,const char* cutsVersion = default_version, const char* cuts = default_cut , int stackOption = 1)
{
  
  printf("\nYour cuts are:\n %s\n",cuts);
  printf("\nYour version is:\n %s\n",cutsVersion);
  
  char outputFile[1000];
  sprintf(outputFile,"%s%s_HFaugmented_halfOfficial_%s",hist_file_path,hist_file_name,cutsVersion);
  
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
      makePlots(cuts,(const char*)outputFile);
      formatPlots((const char*)outputFile,stackOption); break;
    case 1:
      printf("Only making plots.\n");
      makePlots(cuts,(const char*)outputFile); break;
    case 2:
      printf("Only formatting plots.\n");
      formatPlots((const char*)outputFile,stackOption); break;
    default:
      printf("ERROR: 0 for making+formatting plots, 1 for formating only, what are you doing?!\n");
      break;
    }
}

void makePlots(const char* cuts, const char* outputFile)
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
  TH1D     *QCDhist[n_vars_TOTAL][n_types];
  TH1D     *ratio[n_vars_TOTAL];              // Data/MC ratio plots
  double    integrals[n_vars_TOTAL][n_types]; // Integral of hists for stacking
  
  double QCD_HFintegral = 0;
  int numEntries = 0;

  for (int i_var = 0; i_var < n_vars_TOTAL; i_var++)
    {
      printf("\ni_var:  %d\n\nvariable:  %s\n\n", i_var, var[i_var]);

      //Loop to fill histograms for data and each MC jet flavor
      for (int i_type = 0; i_type < n_types; i_type++)
	{
	  // Initialize histogram
	  hist[i_var][i_type] = new TH1D( Form("hist_%d_%d",i_var,i_type), Form("hist_%d_%d",i_var,i_type), nbinsX[i_var], lowX[i_var], highX[i_var]);
	  hist[i_var][i_type]->Sumw2();

	  // Fill/draw histogram
	  if (i_type == 0) data_tree->Draw(Form("%s>>hist_%d_%d",var[i_var],i_var,i_type), Form("%s",cuts), "goff");
	  else /*i_type!=0*/ MC_tree->Draw(Form("%s>>hist_%d_%d",var[i_var],i_var,i_type), Form("weight*(%s&&%s)", particle_cut[i_type],cuts), "goff");
	    
	  integrals[i_var][i_type]=hist[i_var][i_type]->Integral();
	  
	  //Renomalize the B and C contributions, if no enriched files were used, the scale factor is 1.
	  if(i_type == 2 || i_type == 3)
	    {
	      QCDhist[i_var][i_type] = new TH1D(Form("QCDhist_%d_%d",i_var,i_type), Form("QCDhist_%d_%d",i_var,i_type), nbinsX[i_var], lowX[i_var], highX[i_var]);
	      QCDhist[i_var][i_type]->Sumw2();
	      QCD_tree->Draw(Form("%s>>QCDhist_%d_%d",var[i_var],i_var,i_type), Form("weight*(%s&&%s)", particle_cut[i_type],cuts), "goff");
	      QCD_HFintegral = QCDhist[i_var][i_type]->Integral();
	      delete QCDhist[i_var][i_type];//gotta clean up, lest there exist memory leaks
	      double HFscale_factor = QCD_HFintegral/integrals[i_var][i_type];
	      printf("QCD HF Integral for i_type == %d is %f\n", i_type, QCD_HFintegral);
	      printf("HF Scale factor is: %f\n",HFscale_factor);
	      hist[i_var][i_type]->Scale(HFscale_factor);
	    }
	  if (i_type==0) hist[i_var][i_type]->Scale(1/int_lumi);//no need to scale MC by int_lumi, weighted by cross section
	  integrals[i_var][i_type]=hist[i_var][i_type]->Integral();
	}
      
      //second loop to renormalize MC distributions to data
      double scale_factor = (integrals[i_var][0])/(integrals[i_var][2]+integrals[i_var][3]+integrals[i_var][4]);
      for (int i_type = 2; i_type < n_types; i_type++)
	{
	  hist[i_var][i_type]->Scale(scale_factor);
	  integrals[i_var][i_type]=hist[i_var][i_type]->Integral();
	  numEntries = hist[i_var][i_type]->GetEntries(); 
	  printf("%s, i_type = %i \nEntries = %i\n",  leg_label[i_type],i_type, numEntries);
	  printf("Integral = %e\n\n", integrals[i_var][i_type]);
	}
      //we're going to remake the MC hist, so we de-register from current directory
      //if this is not done, we get potential memory leak warninigs
      TObject* old=gDirectory->GetList()->FindObject(Form("hist_%d_1",i_var));
      gDirectory->GetList()->Remove(old);
      //now remake the MC histogram by adding up appropriately scaled histogram types
      TH1D* newMCHist = new TH1D( Form("hist_%d_1",i_var), Form("hist_%d_1",i_var), nbinsX[i_var], lowX[i_var], highX[i_var]);
      newMCHist->Sumw2();
      newMCHist->Add(hist[i_var][2]);newMCHist->Add(hist[i_var][3]);newMCHist->Add(hist[i_var][4]);
      hist[i_var][1]=newMCHist;

      //Ratio plot
      ratio[i_var] = new TH1D(Form("ratio_%d",i_var), Form("data/MC"), nbinsX[i_var], lowX[i_var], highX[i_var]);
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

static void formatPlots(const char* input_file_name, int stackOption)
{
  //Set histogram style
  //gStyle->SetOptStat(1101);
  gStyle->SetOptStat("irMe");
  gROOT->ForceStyle();
  
  // Open file
  //TFile *hist_file = TFile::Open(Form("%s%s_%s.root",pdf_file_path,hist_file_name,version));
  
  printf("opening %s.root\n", input_file_name);
  TFile *hist_file = TFile::Open(Form("%s.root",input_file_name));
  printf("Formatting %s.root\n", input_file_name);
  
  //const char* out_file_name;
  //sprintf(out_file_name, "%s.pdf", input_file_name);
  printf("out_file_name = %s.pdf\n", input_file_name);

  TCanvas *canv[n_vars];
  TLegend *leg[n_vars];
  TH1D    *hist[n_vars][n_types];
  TH1D    *ratio[n_vars];
  THStack *stacked_hist[n_vars];//declare n_vars stacked histograms
  TLine   *one[n_vars];
  
  TCanvas *temp_canv = new TCanvas("temp", "temp", 1200, 600);
  
  temp_canv->Print(Form("%s.pdf(",input_file_name));
  
  for (int i_var=n_vars_low; i_var<n_vars_high+1; i_var++)/*DEBUG*/
    //for (int i_var=0; i_var<n_vars; i_var++)
    {
      printf("Formatting %s\n", var[i_var]);

      canv[i_var] = new TCanvas(Form("canv_%d",i_var),Form("canv_%d",i_var),1200,600);
      canv[i_var]->Divide(2,1);
      
      if (doLogy[i_var]) canv[i_var]->cd(1)->SetLogy();

      stacked_hist[i_var] = new THStack(Form("stacked_hist_%d",i_var),"Stacked MC");//create the stack
      

      for (int i_type=0; i_type<n_types; i_type++)
	{
	  hist[i_var][i_type] = (TH1D *)hist_file->Get(Form("hist_%d_%d",i_var,i_type));
	  formatHist(hist[i_var][i_type], x_label[i_var], y_label);
	  
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
	  
	  if (!doLogy[i_var]) hist[i_var][i_type]->SetMinimum(0);

	  if (i_type>=2&&stackOption==1) 
	    {
	      stacked_hist[i_var]->Add(hist[i_var][i_type]);
	    }
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
	      leg[i_var]->AddEntry(hist[i_var][i_type], Form("%s",leg_label[i_type]), "lp");
	    }
	  else//i_type>=1
	    {
	      if(i_type==1)continue;//i.e. total MC
	      leg[i_var]->AddEntry(hist[i_var][i_type], Form("%s",leg_label[i_type]), "f");
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
      
      one[i_var] = new TLine(lowX[i_var],1,highX[i_var],1);
      one[i_var]->SetLineColor(color[1]);
      one[i_var]->Draw("SAME");
      
      ratio[i_var]->SetLineColor(color[0]);
      ratio[i_var]->Draw("SAME"); // Plot it again, over the line
      
      //output pdf	
      canv[i_var]->Print(Form("%s.pdf",input_file_name));
      printf("i_var = %i and input_file_name=%s\n",i_var,input_file_name);
      
    }
  
  temp_canv->Print(Form("%s.pdf]",input_file_name));
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

static void formatHist(TH1 *h, const char *x_label1, const char *y_label1)
{
  //h->Sumw2();

  h->GetXaxis()->SetTitle(Form("%s",x_label1));
  //h->GetXaxis()->CenterTitle();

  h->GetYaxis()->SetTitle(Form("%s",y_label1));
  //h->GetYaxis()->CenterTitle();
}
