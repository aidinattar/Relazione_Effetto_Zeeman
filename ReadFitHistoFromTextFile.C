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
#include <vector>

using namespace std;

void set_style(){
    // setup graphics
    //gStyle -> SetOptTitle   (   kFALSE   );
    gStyle -> SetOptStat    (      0     );
    gStyle -> SetLabelOffset( 0.01 , "x" );
    gStyle -> SetLabelOffset( 0.005, "y" );
    gStyle -> SetTitleOffset( 1.2  , "x" );
    gStyle -> SetTitleOffset( 1.1  , "y" );
}

void ReadFitHistoFromTextFile(const char *fname, const char *histname=NULL, bool draw=1) {

  cout << "*********************************************************" << endl
       << "-------Inserire il minimo e massimo per ogni range-------" << endl
       << "*********************************************************" << endl
       << endl << endl;

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

  delete [] binc;

  h->SetEntries(h->Integral());

  vector <TF1*> gauss;
  int min, max;

  cout << "massimo e minimo dei 3 picchi: ";
  cin  >>  min  >>  max;
  h->GetXaxis()->SetRange( min, max );

  for( int i = 0; i < 3; ++i ){
    cout << "min e max: ";
    cin  >>  min  >>  max;
    string name = "gauss" + to_string(i);
    gauss.push_back( new TF1 ( name.c_str(), "gaus", min, max ) );

    h->Fit( gauss[i], "R0" ); 
  }

  cout << endl << endl;

  double c1  = gauss[1]->GetParameter( 1 ) -
               gauss[0]->GetParameter( 1 );
  double ec1 = sqrt( pow( gauss[1]->GetParError( 1 ), 2 ) +
                     pow( gauss[0]->GetParError( 1 ), 2 ) );

  cout << "c1 = " << c1 << " +/- " << ec1 << endl;

  double c2  = gauss[2]->GetParameter( 1 ) -
               gauss[1]->GetParameter( 1 );
  double ec2 = sqrt( pow( gauss[2]->GetParError( 1 ), 2 ) +
                     pow( gauss[1]->GetParError( 1 ), 2 ) );

  cout << "c2 = " << c2 << " +/- " << ec2 << endl;

  double Dx  = ( c1 + c2 ) / 2;
  double eDx = sqrt( ec1 * ec1 + ec2 * ec2 ) / 2;

  cout << "Dx = " << Dx << " +/- " << eDx << endl;

  cout << endl << endl;

  vector <double>  FWHM;
  vector <double> eFWHM;
  for( int i = 0; i < 3; ++i ){
     FWHM.push_back( 2 * sqrt( 2 * log( 2 ) * gauss[i]->GetParameter( 2 )));
    eFWHM.push_back( 2 * sqrt( 2 * log( 2 ) * gauss[i]->GetParError(  2 )));
  }

  for( int i = 0; i < 3; ++i )
    cout << i << ". FWHM = " <<  FWHM[i] 
              << " +/- "     << eFWHM[i] << endl;
  
  TCanvas *canvas = new TCanvas( "canvas", "My ROOT Plots 2", 1280, 720 );
  canvas -> SetGrid(); //griglia
  set_style();

  h->Draw();

  for( int i = 0; i < 3; ++i )
    gauss[i]->Draw("SAME");

  cout << endl << endl
       << "----------------by aidin attar & Co----------------" 
       << endl << endl;

  return;
}