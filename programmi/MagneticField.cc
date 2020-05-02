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
#include <cmath>

using namespace std;

double  h = 6.62607015E-34; // Js
double  c = 299792458; // m/s
double μB = 9.27400949E-28; // J/G
double  λ = 585.3E-9; // m


void set_style(){
    // setup graphics
    gStyle -> SetOptTitle   (   kFALSE   );
    gStyle -> SetOptStat    (      0     );
    gStyle -> SetLabelOffset( 0.01 , "x" );
    gStyle -> SetLabelOffset( 0.005, "y" );
    gStyle -> SetTitleOffset( 1.2  , "x" );
    gStyle -> SetTitleOffset( 1.1  , "y" );
}


double mediapesata (vector <double> dati, vector <double> errori){
	double med, sumx=0, sum1e=0;
	for (int i=0; i<dati.size(); i++){
		sumx+=dati.at(i)/(errori.at(i)*errori.at(i));
		sum1e+=1/(errori.at(i)*errori.at(i));
		}
	med = sumx / sum1e;
	return med;
}


double err (vector <double> errori){
	double error, sum1e=0;
	for (int i=0; i<errori.size(); i++)
		sum1e+=1/(errori.at(i)*errori.at(i));
	error = sqrt ( 1 / sum1e ) ;
	return error;
}

void MagneticField(const char *fname, const char *histname=NULL, bool draw=1) {

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

  TH1F *hg = new TH1F(hname,hname,nbins,0,nbins);
  for (int i=1; i<=nbins; i++) {
    hg->SetBinContent(i,binc[i-1]);
  }

  delete [] binc;

  hg->SetEntries(hg->Integral());
  hg->Rebin( 4 );

  vector <TF1*> gauss;
  int min, max;

  cout << "massimo e minimo dei 3 picchi: ";
  cin  >>  min  >>  max;
  hg->GetXaxis()->SetRangeUser( min, max );

  hg->GetXaxis()->SetTitle( "CCD position" );
  hg->GetYaxis()->SetTitle( "counts" );
  hg->GetXaxis()->CenterTitle();
  hg->GetYaxis()->CenterTitle();
  
  int j = 0;
  double  par[18];
  double epar[18];

  double c1, c2, μ1, μ2, σ1, σ2;
  for( int i = 0; i < 3; ++i ){
    cout << "min e max: ";
    cin  >>  min  >>  max;
    string name = "gauss" + to_string(i);
    gauss.push_back( new TF1 ( name.c_str(), "gaus(0) + gaus(3)", min, max ) );

    cout << "valori iniziali IN ORDINE delle due gaussiane: ";
    cin  >> c1 >> μ1 >> σ1
         >> c2 >> μ2 >> σ2;
    gauss[i]->SetParameters( c1, μ1, σ1,
                             c2, μ2, σ2 );

    hg->Fit( gauss[i], "R0+" ); 

    gauss[i]->GetParameters( &par[j] );
    epar[j++] = gauss[i]->GetParError( 0 );
    epar[j++] = gauss[i]->GetParError( 1 );
    epar[j++] = gauss[i]->GetParError( 2 );
    epar[j++] = gauss[i]->GetParError( 3 );
    epar[j++] = gauss[i]->GetParError( 4 );
    epar[j++] = gauss[i]->GetParError( 5 );
  }

  cout << endl << endl;

  double F, eF; // nm/pixel
  cout << "inserire il fattore di conversione con errore: ";
  cin  >> F >> eF;

  double G;
  cout << "inserire il valore di campo magnetico in gauss: ";
  cin  >> G;

  cout << endl;

  double zee1  =       gauss[0]->GetParameter( 4 ) -
                       gauss[0]->GetParameter( 1 );
  double ezee1 = sqrt( pow( gauss[0]->GetParError( 4 ), 2 ) +
                       pow( gauss[0]->GetParError( 1 ), 2 ) );

  cout << "zee1 = " <<  zee1 
       <<  " +/- "  << ezee1 << endl;

  double zee2  = gauss[1]->GetParameter( 4 ) -
                 gauss[1]->GetParameter( 1 );
  double ezee2 = sqrt( pow( gauss[1]->GetParError( 4 ), 2 ) +
                       pow( gauss[1]->GetParError( 1 ), 2 ) );

  cout << "zee2 = " <<  zee2 
       <<  " +/- "  << ezee2 << endl;

  double zee3  = gauss[2]->GetParameter( 4 ) -
                 gauss[2]->GetParameter( 1 );
  double ezee3 = sqrt( pow( gauss[2]->GetParError( 4 ), 2 ) +
                       pow( gauss[2]->GetParError( 1 ), 2 ) );

  cout << "zee3 = " << zee3 << " +/- " << ezee3 << endl;


  double dzee  = ( zee1 + zee2 + zee3 ) / 3;
  double edzee = sqrt( ezee1 * ezee1 + 
                       ezee2 * ezee2 +
                       ezee3 * ezee3 ) / 3;

  cout << "dzee = " <<  dzee 
       <<  " +/- "  << edzee 
       << endl      << endl;

  double  dlzee = F * dzee;
  double edlzee = dlzee * sqrt( pow( eF / F, 2 ) + pow( edzee / dzee, 2 ));

  cout << "dlzee = " <<  dlzee 
       <<  " +/- "   << edlzee 
       << endl;
  cout << endl << endl;

  double  g = h * c * dlzee * pow( 10, -9 ) /
            ( 2 * μB * G * λ * λ );
  double eg = g * sqrt( 0.0001 + pow( edlzee / dlzee, 2 ));

  cout << "g =  " <<  g
       << " +/- " << eg << endl;

  TCanvas *canvas = new TCanvas( "canvas", "My ROOT Plots 2", 1280, 720 );
  canvas -> SetGrid(); //griglia
  set_style();

  gPad   -> SetTopMargin   ( 0.01  );
  gPad   -> SetBottomMargin( 0.09  );
  gPad   -> SetRightMargin ( 0.01  );
  gPad   -> SetLeftMargin  ( 0.075 );

  hg->Draw();

  for( int i = 0; i < 3; ++i )
    gauss[i]->Draw("SAME");


  vector <TPaveText*> legend;
  string name;
  int k = 0;

  legend.push_back( new TPaveText( .05, .75, .25, .95, "NDC" ));
  legend.push_back( new TPaveText( .35, .75, .55, .95, "NDC" ));
  legend.push_back( new TPaveText( .65, .75, .85, .95, "NDC" ));

  for( int i = 0; i < 3; ++i ){
    name = "c_{1} = " + to_string(  par[k] ) +
           " #pm "       + to_string( epar[k] );
    legend[i]->AddText( name.c_str() );
    k++;

    name = "#mu_{1} = "     + to_string(  par[k] ) +
           " #pm "       + to_string( epar[k] );
    legend[i]->AddText( name.c_str() );
    k++;

    name = "#sigma_{1} = "    + to_string(  par[k] ) +
           " #pm "       + to_string( epar[k] );
    legend[i]->AddText( name.c_str() );
    k++;

    name = "c_{2} = " + to_string(  par[k] ) +
           " #pm "       + to_string( epar[k] );
    legend[i]->AddText( name.c_str() );
    k++;

    name = "#mu_{2} = "     + to_string(  par[k] ) +
           " #pm "       + to_string( epar[k] );
    legend[i]->AddText( name.c_str() );
    k++;

    name = "#sigma_{2} = "    + to_string(  par[k] ) +
           " #pm "       + to_string( epar[k] );
    legend[i]->AddText( name.c_str() );
    k++;

    legend[i]->Draw();
  }

/*  vector <TLegend*> legend;
  string name;
  int k = 0;

  for( int i = 0; i < 3; ++i ){
    legend.push_back( new TLegend );

    name = "Constant1  " + to_string(  par[k] ) +
           " +/- "      + to_string( epar[k] );
    legend[i]->AddEntry(( TObject* )0, name.c_str(), "");
    k++;

    name = "Mean1  "     + to_string(  par[k] ) +
           " +/- "      + to_string( epar[k] );
    legend[i]->AddEntry(( TObject* )0, name.c_str(), "");
    k++;

    name = "Sigma1  "    + to_string(  par[k] ) +
           " +/- "      + to_string( epar[k] );
    legend[i]->AddEntry(( TObject* )0, name.c_str(), "");
    k++;

    name = "Constant2  " + to_string(  par[k] ) +
           " +/- "      + to_string( epar[k] );
    legend[i]->AddEntry(( TObject* )0, name.c_str(), "");
    k++;

    name = "Mean2  "     + to_string(  par[k] ) +
           " +/- "      + to_string( epar[k] );
    legend[i]->AddEntry(( TObject* )0, name.c_str(), "");
    k++;

    name = "Sigma2  "    + to_string(  par[k] ) +
           " +/- "      + to_string( epar[k] );
    legend[i]->AddEntry(( TObject* )0, name.c_str(), "");
    k++;

    legend[i]->Draw( "SAME" );
  }
*/

  cout << endl << endl
       << "----------------by aidin attar & Co----------------" 
       << endl << endl;

  return;
}