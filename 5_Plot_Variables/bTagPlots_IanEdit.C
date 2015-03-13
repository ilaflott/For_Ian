/*
Plot maker+formatter. Written for muon-tagged b-jet RpA + RAA
Inherited from: Leo Yu
Heavily Modified/Edited: Ian Laflotte
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

// File parameters
const char *data_file_path = "/net/hisrv0001/home/ilaflott/pp_MC_2760GeV_bTag_forests_ntuples/NTuples/data_leo/";
const char *data_file_name = "data_updated";

const char *MC_file_path = "/net/hisrv0001/home/ilaflott/pp_MC_2760GeV_bTag_forests_ntuples/NTuples/augmented_Samples/";
//const char *MC_file_path = "/net/hisrv0001/home/ilaflott/pp_MC_2760GeV_bTag_forests_ntuples/NTuples/kurts_QCDMC/"
//const char *MC_file_path = "/net/hisrv0001/home/ilaflott/pp_MC_2760GeV_bTag_forests_ntuples/NTuples/halfOfficial_HFMC/";
const char *MC_file_name = "MC_HFaugmented_halfOfficial_updated_schemeA";
//const char *MC_file_name = "MC_HFaugmented_halfOfficial_noCuts_schemeA";
//const char *MC_file_name = "BJet_halfOfficial_updated_schemeA";
//const char *MC_file_name = "CJet_halfOfficial_updated_schemeA";

const char* QCD_file_name = "QCD_updated";
const char* QCD_file_path = "/net/hisrv0001/home/ilaflott/pp_MC_2760GeV_bTag_forests_ntuples/NTuples/kurts_QCDMC/";

const char *hist_file_path = "/net/hisrv0001/home/ilaflott/pp_MC_2760GeV_bTag_forests_ntuples/Histograms/sevil_debug_plots/3.12.15_OvernightPlots/";
const char *hist_file_name = "bTagPlots_pp";

const char *pdf_file_path = "/net/hisrv0001/home/ilaflott/pp_MC_2760GeV_bTag_forests_ntuples/Histograms/sevil_debug_plots/3.12.15_OvernightPlots/";

// Histogram parameters
const int n_vars = 18; // Number of variables to plot
//const int n_vars = 4; /*debug*/
const int n_types = 5; // data, MC, b, c, udsg (0,1,2,3,4...)

const char *y_label[] =
  {
    "xsec (mb)", "xsec (mb)", "xsec (mb)", "xsec (mb)",
    "xsec (mb)", "xsec (mb)", "xsec (mb)", "xsec (mb)", "xsec (mb)",
    "xsec (mb)", "xsec (mb)", "xsec (mb)", "xsec (mb)", "xsec (mb)",
    "xsec (mb)", "xsec (mb)", "xsec (mb)", "xsec (mb)", "xsec (mb)"
  };

const char *x_label[] =
  {
    "jtpt (GeV)",    "jteta",  "jtphi (rad)",         
    "mupt (Gev)", "muptrel (GeV)",    "muphi (rad)",   "mudr", "mueta",
    "discr_ssvHighEff", "discr_ssvHighPur",
    "nsvtx", "svtxntrk", "svtxdl (cm)",  "svtxdls",
    "ip2d (cm)","ip2dSig","ip3d (cm)","ip3dSig"
  };

const char *var[] =
  {
    "jtpt",    "jteta",  "jtphi",                       //jets
    "mupt", "muptrel","mueta",    "muphi",   "mudr",     //muons
    "discr_ssvHighEff", "discr_ssvHighPur",             //discriminators
    "nsvtx", "svtxntrk", "svtxdl", "svtxdls",            //secondary vertex
    "ip2d","ip2dSig","ip3d","ip3dSig"                 //impact parameter
  };

const char *particle_cut[] = 
  { 
    "1",         // Data
    "1",         // Total MC
    "abs(refparton_flavorForB)==5",//b
    "abs(refparton_flavorForB)==4",//c
    "(abs(refparton_flavorForB)==1 || abs(refparton_flavorForB)==2 || abs(refparton_flavorForB)==3 || abs(refparton_flavorForB)==21)"//udsg
  };
                                                                                           
const int   nbinsX[] = {   25, 15,       15, /**/ 20, 25, 15,       15,   15,  /**/  6,  6, /**/  5, 12, 8,  16, /**/   40,   40,   40,   40 };
const double  lowX[] = {   20, -3, -3.14159, /**/  0,  0, -3, -3.14159,    0,  /**/  0,  0, /**/  0,  0, 0,   0, /**/ -0.1, -0.1, -0.1, -0.1 };
const double highX[] = {  270,  3,  3.14159, /**/ 80, 10,  3,  3.14159,  0.5,  /**/  6,  6, /**/  5, 12, 4, 240, /**/  0.1,  0.1,  0.1,  0.1 };
const bool  doLogy[] = {    1,  1,        0, /**/  1,  1,  1,        0,    1,  /**/  1,  1, /**/  1,  1, 1,   1, /**/    0,    0,    0,    0 };

const float int_lumi = 4209000000;//inverse millibarns of data. according to lumiCalc2.py, golden lumimask for HLT_PAMu3_v1, 4.209 pb of data.

// Selections
const char *event_cut = "HLT_PAMu3_v1";
//const char *event_cut = "HLT_PAMu7_v1";
//const char *event_cut = "HLT_PAMu12_v1";
//const char *event_cut = "HLT_PAMu3PFJet40_v1";

const char *default_cut = "vz<15&&vz>-15&&jteta<2&&jteta>-2&&jtpt>30&&HLT_PAMu3_v1&&mupt!=0&&mupt/rawpt<0.95&&svtxdl>0.01&&svtxdl<2.5&&svtxdls>3.0&&svtxm<6.5";
const char *default_version = "vz15_jteta2_jtpt40_HLTPAMu3v1_muCut_WCut_svtxCut";

const int       color[]  = { kBlack, kGray+3, kRed-7, kGreen-6, kBlue-7};
const int    lineColor[] = { kBlack, kWhite, kRed-7, kGreen-6, kBlue-7};
const char *leg_label[]  = { "Data pp", "MC", "b", "c", "udsg" };

//MAIN FUNCTIONS
//stackOption == 0 -> Overlaid flavor curves
//stackOption == 1 -> stacked flavor curves
//NOTE fields which arent specified default to the values here. If one only wants to change one or two of the parameters,
//then they must be submitted in order and one after the other. Ex. you want to change the cuts, but because of the order
//of the arguments and the ambiguity inherent, one must also specify a version and an option. Just specifying one string as an input will
// be taken in as an input value to cutsVersion, even if what you really wanted to change is the cuts 
void bTagPlots_IanEdit(const char* cutsVersion = default_version, int option = 0, const char* cuts = default_cut , int stackOption = 1)
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
  TH1D     *hist[n_vars][n_types];;      // Data and MC histograms
  TH1D     *QCDhist[n_vars][n_types];
  TH1D     *ratio[n_vars];              // Data/MC ratio plots
  double    integrals[n_vars][n_types]; // Integral of hists for stacking
  
  double QCD_HFintegral = 0;
  int numEntries = 0;
  //for (int i_var = 16; i_var < n_vars; i_var++)/*DEBUG*/
  for (int i_var = 0; i_var < n_vars; i_var++)
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

      //remake the MC histogram by adding up appropriately scaled histogram types
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
      //for (int i_type=0; i_type<2; i_type++)
	{
	  hist[i_var][i_type]->Write();
	}
      
      ratio[i_var]->Write();
    }

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
  //for (int i_var=16; i_var<n_vars; i_var++)/*DEBUG*/
  for (int i_var=0; i_var<n_vars; i_var++)
    {
      printf("Formatting %s\n", var[i_var]);

      canv[i_var] = new TCanvas(Form("canv_%d",i_var),Form("canv_%d",i_var),1200,600);
      canv[i_var]->Divide(2,1);
      
      if (doLogy[i_var]) canv[i_var]->cd(1)->SetLogy();

      stacked_hist[i_var] = new THStack(Form("stacked_hist_%d",i_var),"Stacked MC");//create the stack
      

      for (int i_type=0; i_type<n_types; i_type++)
	{
	  hist[i_var][i_type] = (TH1D *)hist_file->Get(Form("hist_%d_%d",i_var,i_type));
	  formatHist(hist[i_var][i_type], x_label[i_var], y_label[i_var]);
	  
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

///////////////////////////////////////
//OLDER makePlots() and formatPlots()//
///////////////////////////////////////
/*
void makePlots()
{
  //SETUP

  // Open files and trees
  cout << "opening input data files and trees" << endl;
  TFile *data_file = TFile::Open(Form("%s%s.root",data_file_path,data_file_name));
  TTree *data_tree = (TTree *)data_file->Get("nt");
  data_tree->SetMakeClass(1);

  cout << "opening input MC files and trees" << endl;
  TFile *MC_file = TFile::Open(Form("%s%s.root",MC_file_path,MC_file_name));
  TTree *MC_tree = (TTree *)MC_file->Get("nt");
  MC_tree->SetMakeClass(1);

  // Output file
  cout << "opening outputfile" << endl;
  TFile *out_file = new TFile(Form("%s%s_%s.root",hist_file_path,hist_file_name,version), "RECREATE");
  out_file->cd();

  // Declare histograms arrays
  TH1D     *hist[n_vars][n_types];      // Data and MC histograms
  TH1D     *ratio[n_vars];              // Data/MC ratio plots
  double    integrals[n_vars][n_types]; // Integral of hists for stacking

  //CREATE HISTOGRAMS

  printf("Opened files.\n");

    // For each variable:
    for (int i_var = 0; i_var < n_vars; i_var++)
      {
        printf("\ni_var:  %d\nvariable:  %s\n", i_var, var[i_var]);

        // Fill histograms for data and each MC jet flavor
        for (int i_type = 0; i_type < n_types; i_type++)
	  //for (int i_type = 0; i_type < n_1; i_type++)
	  {
	    
            // Initialize histogram
            hist[i_var][i_type] = new TH1D( Form("hist_%d_%d",i_var,i_type), Form("hist_%d_%d",i_var,i_type), nbinsX[i_var], lowX[i_var], highX[i_var]);
	    hist[i_var][i_type]->Sumw2();
	    
            // Fill/draw histogram
            if (i_type == 0) data_tree->Draw(Form("%s>>hist_%d_%d",var[i_var],i_var,i_type), Form("weight*(%s&&%s&&%s)",event_cut,mu_cut,svtx_cut), "goff");
	    else MC_tree->Draw(Form("%s>>hist_%d_%d",var[i_var],i_var,i_type), Form("weight*(%s&&%s&&%s&&%s)", event_cut, particle_cut[i_type],mu_cut,svtx_cut), "goff");
		      
            // Calc integral of hist
            integrals[i_var][i_type] = hist[i_var][i_type]->Integral();
            printf("\t %s Integral = %e\n", leg_label[i_type], integrals[i_var][i_type]);
	  }

        // Scale MC to data
        hist[i_var][1]->Scale(integrals[i_var][0]/integrals[i_var][1]);
        printf("\n\t MC scale factor = %e\n", integrals[i_var][0]/integrals[i_var][1]);

        // Scale the quark flavor hists for creating the stacked MC histogram w/
        // fractional flavor contributions
        double stacked_int = integrals[i_var][2] + integrals[i_var][3] + integrals[i_var][4]// + integrals[i_var][5];
        double data_int = integrals[i_var][0];
        for (int i_type=2; i_type<n_types; i_type++)
	  {
	    // See notes for how I got scale factor
	    double scale_factor = data_int/stacked_int;
	    hist[i_var][i_type]->Scale(scale_factor);
	    if (i_type==2) printf("\t stack scale factor = %e\n", scale_factor);
	  }

        // Ratio plot
        ratio[i_var] = new TH1D(Form("ratio_%d",i_var), Form("data/MC"), nbinsX[i_var], lowX[i_var], highX[i_var]);
	ratio[i_var]->Sumw2();
	if (i_var != 13) ratio[i_var]->Divide(hist[i_var][0],hist[i_var][1],1,1,"b");

        ///////////////////////////////////////////////////////////////////////////////////////////////
        //
        //          OUTPUT TO ROOT FIL
        //
        ///////////////////////////////////////////////////////////////////////////////////////////////

        // Output to .root file
        for (int i_type=0; i_type<n_types; i_type++)
	  {
	    hist[i_var][i_type]->Write();
	  }
        ratio[i_var]->Write();
    }

    ///////////////////////////////////////////////////////////////////////////////////////////////
    //
    //          CLEANUP
    //
    ///////////////////////////////////////////////////////////////////////////////////////////////

    data_file->Close();
    MC_file->Close();
    out_file->Close();
}

void formatPlots()
{
    ///////////////////////////////////////////////////////////////////////////////////////////////
    //
    //          SETUP
    //
    ///////////////////////////////////////////////////////////////////////////////////////////////

    // Set histogram style
    gStyle->SetOptStat(1101);
    gROOT->ForceStyle();
    
    // Open file
    TFile *hist_file = TFile::Open(Form("%s%s_%s.root",pdf_file_path,hist_file_name,version));
    printf("Formatting %s%s_%s.root\n", pdf_file_path, hist_file_name,version);

    // Output file
    //double version = -1;
    //char version[200] = "";
    //printf("Enter a version number ");
    //printf("Enter a version title ");
    //scanf("%s", &version);
    char out_file_name[200];
    sprintf(out_file_name, "%s%s_%s.pdf", pdf_file_path, hist_file_name, version);

    //const int       color[] = { data,MC,B,C,usdg/QCD };

    // const int       color[] = { kPink+6, kGray+3, kCyan+4, kCyan+3, kCyan+2, kCyan+1 };
    // const char *leg_label[] = { "Data pp", "MC", "b", "c", "uds", "gluon" };

    TCanvas *canv[n_vars];
    TLegend *leg[n_vars];
    TH1D    *hist[n_vars][n_types];
    TH1D    *ratio[n_vars];
    THStack *stacked_hist[n_vars];
    TLine   *one[n_vars];

    TCanvas *temp_canv = new TCanvas("temp", "temp", 1200, 600);
    temp_canv->Print(Form("%s(",out_file_name));

    ///////////////////////////////////////////////////////////////////////////////////////////////
    //
    //          FORMAT HISTOGRAMS
    //
    ///////////////////////////////////////////////////////////////////////////////////////////////

    for (int i_var=0; i_var<n_vars; i_var++)
      {
        printf("Formatting %s\n", var[i_var]);
	canv[i_var] = new TCanvas(Form("canv_%d",i_var),Form("canv_%d",i_var),1200,600);
        canv[i_var]->Divide(2,1);
	//canv[i_var]->Divide(1,2);

	if (doLogy[i_var]) canv[i_var]->cd(1)->SetLogy();

        //Histograms
        stacked_hist[i_var] = new THStack(Form("stacked_hist_%d",i_var), "Stacked MC");
	//stacked_hist->Sumw2();

        for (int i_type=0; i_type<n_types; i_type++)
	  {
	    hist[i_var][i_type] = (TH1D *)hist_file->Get(Form("hist_%d_%d",i_var,i_type));

	    formatHist(hist[i_var][i_type], x_label[i_var], y_label[i_var]);

	    hist[i_var][i_type]->SetMarkerColor(color[i_type]);
	    hist[i_var][i_type]->SetLineColor(lineColor[i_type]);
	    //hist[i_var][i_type]->SetLineColor(kBlack+3);
	    hist[i_var][i_type]->SetLineWidth(1);

	    if (i_type==0 )
	      {
		//hist[i_var][i_type]->SetLineWidth(2);
		hist[i_var][i_type]->SetMarkerStyle(8);
		hist[i_var][i_type]->SetMarkerSize(1.1);
	      }

	    //hist[i_var][i_type]->GetYaxis()->SetMaximum(highY[i_var]);
	    //hist[i_var][i_type]->GetYaxis()->SetRangeUser(1,highY[i_var]);

	    if (i_type>=1)
	      {
		hist[i_var][i_type]->SetFillColor(color[i_type]);
		//hist[i_var][i_type]->SetFillStyle(3004 + i_type % 2); //various light stripes
		hist[i_var][i_type]->SetFillStyle(1001);//solid fill
	      }

            if (!doLogy[i_var]) hist[i_var][i_type]->SetMinimum(0);
            if (i_type>=2) stacked_hist[i_var]->Add(hist[i_var][i_type]);
        }

        canv[i_var]->cd(1);

	////draw mc curve with error
        hist[i_var][1]->Draw("HIST E1");

	////draw stacked udsg+b+c curve with error
        //stacked_hist[i_var]->Draw("HIST E SAME");
	//stacked_hist[i_var]->Draw("E SAME");
	//stacked_hist[i_var]->Draw("HIST SAME");
	
	//draw data curve
        //hist[i_var][0]->Draw("HIST E SAME");
	hist[i_var][0]->Draw("E SAME");

        // Legend
        //leg[i_var] = new TLegend(0.7,0.75,0.85,0.90); // (xmin,ymin,xmax,ymax)
	leg[i_var] = new TLegend(0.5,0.78,0.65,0.93); // (xmin,ymin,xmax,ymax)

	formatLeg(leg[i_var]);

	//add entries to legend
	for (int i_type = 0; i_type < n_types; i_type++)
	  {
	    if(i_type==0)
	      {
		leg[i_var]->AddEntry(hist[i_var][i_type], Form("%s",leg_label[i_type]), "lp");
	      }
	    else
	      {
		//if (i_type == 1) continue;
		//leg[i_var]->AddEntry(hist[i_var][i_type], Form("%s",leg_label[i_type]), "f");
		//if plotting mcOnly
		if (i_type==1) leg[i_var]->AddEntry(hist[i_var][i_type], Form("%s",leg_label[i_type]), "f");
		else continue;
	      }
	  }
	leg[i_var]->Draw();

        // Ratio plots
        ratio[i_var] = (TH1D *)hist_file->Get(Form("ratio_%d",i_var));

        formatHist(ratio[i_var], x_label[i_var], "data/MC");

        ratio[i_var]->GetYaxis()->CenterTitle();
	ratio[i_var]->SetMarkerColor(color[0]);
	ratio[i_var]->SetMarkerStyle(8);
        ratio[i_var]->SetMarkerSize(0.8);
	ratio[i_var]->SetMinimum(0.5);
        ratio[i_var]->SetMaximum(1.5);	

        canv[i_var]->cd(2);

        ratio[i_var]->Draw();

	one[i_var] = new TLine(lowX[i_var],1,highX[i_var],1);
        one[i_var]->SetLineColor(color[1]);
        one[i_var]->Draw("SAME");

	ratio[i_var]->SetLineColor(color[0]);
        ratio[i_var]->Draw("SAME"); // Plot it again, over the line

        //output pdf
        canv[i_var]->Print(out_file_name);
    }

     //CLEAN UP
    temp_canv->Print(Form("%s]",out_file_name));
    hist_file->Close();
}


void fitDebugging()
{

  //TFile *_file3 = TFile::Open("/net/hisrv0001/home/ilaflott/pp_MC_2760GeV_bTag_forests_ntuples/NTuples/data_leo/data.root");
  //TFile *_file2 = TFile::Open("/net/hisrv0001/home/ilaflott/pp_MC_2760GeV_bTag_forests_ntuples/NTuples/halfOfficial_HFMC/BJet_halfOfficial_schemeB.root");
  //TFile *_file1 = TFile::Open("/net/hisrv0001/home/ilaflott/pp_MC_2760GeV_bTag_forests_ntuples/NTuples/halfOfficial_HFMC/CJet_halfOfficial_schemeB.root");
  //TFile *_file0 = TFile::Open("/net/hisrv0001/home/ilaflott/pp_MC_2760GeV_bTag_forests_ntuples/NTuples/kurts_QCDMC/QCD_sparecopy.root");

  TFile *_file0 = TFile::Open("/net/hisrv0001/home/ilaflott/pp_MC_2760GeV_bTag_forests_ntuples/NTuples/data_leo/data_noCuts.root");
  TFile *_file1 = TFile::Open("/net/hisrv0001/home/ilaflott/pp_MC_2760GeV_bTag_forests_ntuples/NTuples/augmented_Samples/MC_HFaugmented_halfOfficial_noCuts_schemeA.root"); //an hadd of the three MC files
  //TFile *_file2 = TFile::Open("/net/hisrv0001/home/ilaflott/pp_MC_2760GeV_bTag_forests_ntuples/NTuples/halfOfficial_HFMC/BJet_halfOfficial_noCuts_schemeA.root");
  //TFile *_file1 = TFile::Open("/net/hisrv0001/home/ilaflott/pp_MC_2760GeV_bTag_forests_ntuples/NTuples/halfOfficial_HFMC/CJet_halfOfficial_noCuts_schemeA.root");
  //TFile *_file0 = TFile::Open("/net/hisrv0001/home/ilaflott/pp_MC_2760GeV_bTag_forests_ntuples/NTuples/kurts_QCDMC/QCD_noCuts.root");

  TTree* nt0 = (TTree*)_file0->Get("nt");
  TTree* nt1 = (TTree*)_file1->Get("nt");

  TH1F* htest0 = new TH1F("htest0","",30,0,600);
  TH1F* htest1 = new TH1F("htest1","",30,0,600);

  htest0->Sumw2();
  htest1->Sumw2();

  TCanvas* c1=  new TCanvas("c1","",600,600);

  htest0->SetMarkerStyle(8);
  htest0->SetMarkerStyle(8);
  htest0->SetMarkerColor(kBlack);

  htest0->SetLineColor(kBlack);
  htest1->SetFillColor(kGray+1);

  htest1->SetFillStyle(1001);//solid fill

  nt0->Draw("jtpt>>htest0",Form(""),"goff");
  nt1->Draw("jtpt>>htest1",Form("weight"),"goff");


  htest0->Scale(1/int_lumi);

  c1->cd(1)->SetLogy();

  //printf("\n\t MC scale factor = %e\n", htest3->Integral()/MC_Integral);

  TLegend* Legend = new TLegend(.5,0.78,0.65,0.93);

  Legend->AddEntry(htest0,"data ntuple","lp");
  Legend->AddEntry(htest1,"MC Curve","l");
  //Legend->AddEntry(htest1,"c ntuple","l");
  //Legend->AddEntry(htest0,"udsg ntuple","l");

  htest0->GetXaxis()->SetTitle("jtpt (GeV), no Cuts");
  htest0->GetXaxis()->SetTitleSize(0.03);
  htest0->GetXaxis()->SetLabelSize(0.03);

  htest0->GetYaxis()->SetTitle("#sigma (mb)");
  htest0->GetYaxis()->SetTitleSize(0.03);
  htest0->GetYaxis()->SetLabelSize(0.02);

  htest0->Draw("SCAT E1");
  htest1->Draw("HIST SAME E");

  Legend->Draw("SAME");


}*/



