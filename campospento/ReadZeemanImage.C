//
// macro ReadZeemanImage.C per leggere i bidimensionali creati dal programma di controllo
// degli esperimenti Zeeman1 e Zeeman2 (files .zee)
// 
#include <Riostream.h>
#include <stdlib.h>
#include <TROOT.h>
#include <TSystem.h>
#include "TNtuple.h"
#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TH1D.h"
#include "TH2S.h"
#include "TF1.h"
#include <stdio.h>

using namespace std;

int ReadZeemanImage(const char *fname, int bkgfrom=-1, int bkgto=-1, int projYfrom=-1, int projYto=-1, char *hname=NULL) {
  // fname = file name (.zee)
  // bkgfrom, bkgto =  range in x (colonne) dove calcolare il background
  // projYfrom, projYto = range in x (colonne) dove proiettare lo spettro di segnale
  // hname = se c'e', crea un file root e un file testo della proiezione

  const int npixels=7926;
  const int maxcolumns=1024;

  if (gROOT->FindObject("hzee")) delete (TH2S*)gROOT->FindObject("hzee");
  if (gROOT->FindObject("hzee2")) delete (TH2S*)gROOT->FindObject("hzee2");
  if (gROOT->FindObject("czee")) delete (TCanvas*)gROOT->FindObject("czee");
  if (gROOT->FindObject("czee2")) delete (TCanvas*)gROOT->FindObject("czee2");
  if (gROOT->FindObject("cpx")) delete (TCanvas*)gROOT->FindObject("cpx");
  if (gROOT->FindObject("hzee_px")) delete (TH1D*)gROOT->FindObject("hzee_px");
  if (gROOT->FindObject("hzee_py")) delete (TH1D*)gROOT->FindObject("hzee_py");
  if (gROOT->FindObject("cpy")) delete (TCanvas*)gROOT->FindObject("cpy");


  FILE *pf = fopen(fname,"rb");
  if (!pf) {
    std::cout << "file not found! " << std::endl;
    return -1;
  }

  //short int (*image)[npixels] = new short int[maxcolumns][npixels];
  short int *image = new short int[npixels];

  int ncolumns=0;
  while ( fread(image, npixels*2, 1, pf)  ) {
    ncolumns++;
  }
  fclose(pf); pf=NULL;
  std::cout << "number of columns : " << ncolumns << std::endl;
  if (bkgto>ncolumns) bkgto=ncolumns;


  TH2S *hzee = new TH2S("hzee","Original Zeeman image",ncolumns,0,ncolumns,npixels,0,npixels);

  pf = fopen(fname,"rb");
  int icol=0;
  while ( fread(image, npixels*2, 1, pf)  ) {
    int sumpx=0;
    for (int i=0; i<npixels; i++) {
      hzee->SetBinContent(icol+1, i+1, image[i]);
      sumpx += image[i];
    }    
    //if (sumpx<10) printf("sumpx[%d]=%d\n",icol,sumpx);
    icol++;
  }
  fclose(pf);
  delete [] image;

  if (hname) hzee->SetName(hname); 
  hzee->SetStats(0);

  new TCanvas("czee");
  hzee->Draw("col2");

  TH1D *hpro=NULL;
  hpro = hzee->ProjectionX("_px");

  new TCanvas("cpx");
  hpro->Draw();


  if (bkgfrom<0) return 0;
  if (bkgfrom==0) bkgfrom=1; // first bin is 1

  TH2S *hzee2=NULL;
  TH1D *hb=NULL;

  if ( bkgto > bkgfrom ) {
    
    hzee2 = (TH2S*) hzee->Clone("hzee2");
    hzee2->SetTitle("Background subtracted Zeeman image");    

    std::cout << "subtracting background... " << std::endl;    
    for (int i=0; i<npixels; i++) {      
      // calculate bkg for row i
      float bkg = 0;
      for (int k=bkgfrom; k<=bkgto; k++) bkg += hzee->GetBinContent(k,i+1); // ok
      bkg /= (bkgto-bkgfrom+1);
      // subtract bkg from row i
      for (int j=0; j<ncolumns; j++) {
	short int content = hzee->GetBinContent(j+1,i+1) - bkg;
	hzee2->SetBinContent(j+1,i+1,content);
      }
    }
    new TCanvas("czee2");
    hzee2->Draw("col2");
  }

  if (projYfrom<0) return 0;
  
  hpro=NULL;
  TCanvas *cpy=NULL;
  if ( (projYto>projYfrom) ) { // project image on Y axis (pixel axis)
    if (hzee2) {
      hpro = hzee2->ProjectionY("_py",projYfrom, projYto);
    }
    cpy=new TCanvas("cpy");
    hpro->Draw();
  }

  if (hname && cpy) {
    TString fout=hname;
    fout += ".root";
    cpy->Print(fout.Data());
    TString fout2=hname;
    fout2 += ".txt";
    FILE *pf=fopen(fout2.Data(),"w");
    if (!pf) printf("output file open error!\n");
    else {
      fprintf(pf,"# Zeeman image prpjection from %d to %d\n",projYfrom,projYto);
      for (int i=0; i<hpro->GetNbinsX(); i++)
	fprintf(pf,"%.0f\n",hpro->GetBinContent(i+1));
      fclose(pf);
      printf("projection saved!\n");
    }
  }

  return 0;
}

