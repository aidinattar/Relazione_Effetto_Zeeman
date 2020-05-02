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
#include "TH1F.h"
#include "TF1.h"
#include <iostream>
#include <fstream>
#include <string>

using namespace std;

TH1F *ReadHistoFromTextFile(const char *fname, const char *histname=NULL, bool draw=1) {
  
  ifstream infile;
  string line;
  const int maxbins=4096*4;

  infile.open(fname, ios::in);
  if (!infile.is_open()) {
    std::cout << "file not found! " << std::endl;
    return NULL;
  }

  float *binc = new float[maxbins];
  int nbins=0;

  while ( getline(infile,line) ) {
    if (line[0] == '#') continue; // skip comments
    binc[nbins++] = atof(line.c_str()); // legge 1 float dalla riga
  }

  infile.close();

  std::cout << "creating new histo with " << nbins << " bins..." << std::endl;

  char hname[100];
  if (histname) strcpy(hname,histname);
  else strcpy(hname,"histo");

  if (gROOT->FindObject(hname)) gROOT->FindObject(hname)->Delete();

  TH1F *h = new TH1F(hname,hname,nbins,0,nbins);
  for (int i=1; i<=nbins; i++) {
    h->SetBinContent(i,binc[i-1]);
  }

  h->SetEntries(h->Integral());

  if (draw) h->Draw();

  delete [] binc;

  return h;
}
