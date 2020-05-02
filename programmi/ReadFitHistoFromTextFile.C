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


const double lambda =              585.3E-9;
const double      d =                250E-3;
const double      n =      1.51198429503326;
const double     en =     0.001265109776924;
const double   dndl = -4.19418931599459E-05;
const double  edndl =  2.89106089571122E-07;


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
  h->Rebin( 6 );

  vector <TF1*> gauss;
  int min, max;

  cout << "massimo e minimo dei 3 picchi: ";
  cin  >>  min  >>  max;
  h->GetXaxis()->SetRangeUser( min, max );

  h->GetXaxis()->SetTitle( "CCD position" );
  h->GetYaxis()->SetTitle( "counts" );
  h->GetXaxis()->CenterTitle();
  h->GetYaxis()->CenterTitle();

  int j = 0;
  double  par[9];
  double epar[9];

  for( int i = 0; i < 3; ++i ){
    cout << "min e max: ";
    cin  >>  min  >>  max;
    string name = "gauss" + to_string(i);
    gauss.push_back( new TF1 ( name.c_str(), "gausn", min, max ) );

    h->Fit( gauss[i], "R0+" ); 
    gauss[i]->GetParameters( &par[j] );
    epar[j++] = gauss[i]->GetParError( 0 );
    epar[j++] = gauss[i]->GetParError( 1 );
    epar[j++] = gauss[i]->GetParError( 2 );

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

  double  Dx = ( c1 + c2 ) / 2;
  double eDx = sqrt( ec1 * ec1 + ec2 * ec2 ) / 2;

  cout << "Dx = " << Dx << " +/- " << eDx << endl;

  cout << endl << endl;

  vector <double>  FWHM;
  vector <double> eFWHM;
  for( int i = 0; i < 3; ++i ){
     FWHM.push_back( 2 * sqrt( 2 * log( 2 ) ) * gauss[i]->GetParameter( 2 ));
    eFWHM.push_back( 2 * sqrt( 2 * log( 2 ) ) * gauss[i]->GetParError(  2 ));
  }

  for( int i = 0; i < 3; ++i )
    cout << i << ". FWHM = " <<  FWHM[i] 
              << " +/- "     << eFWHM[i] << endl;

  double meaFWHM = mediapesata( FWHM, eFWHM );
  double rmsFWHM = err(               eFWHM );
  cout << endl << "FWHM = " << meaFWHM 
               <<  " +/- "  << rmsFWHM
               <<   endl;

  double i, ei;
  cout << "inserire l'angolo i e il suo errore: ";
  cin  >>  i  >>  ei;


//  GUARDA QUI
//  CONTROLLA CHE dlambda ABBIA L'ORDINE DI GRANDEZZA GIUSTO
//  E IMPLEMENTA IL CALCOLO DELL'ERRORE SOSTITUENDO IL VALORE
//  SCRITTO IN PRECEDENZA: IL CALCOLO DEVE ESSERE IN 
//  NANOMETRI, LE VARIABILI CON LE COSTANTI E GLI ERRORI SONO
//  IN CIMA. BUON LAVORO :)


  double  dlambda = lambda*lambda /      ( 2 * d )   * 
                    sqrt( n*n - pow( sin( i ), 2 ) ) /
                  ( n*n - pow( sin( i ), 2 )         -
                    n * lambda * dndl              );
  double edlambda = 0.000160081307397 ; // nm

  double  F = dlambda / Dx;
  double eF = F * sqrt( pow( edlambda / dlambda, 2 ) + pow( eDx / Dx, 2 ));
  cout <<   endl  << endl 
       <<  "F = " << F 
       << " +/- " << eF << endl;

  TCanvas *canvas = new TCanvas( "canvas", "My ROOT Plots 2", 1280, 720 );
  canvas -> SetGrid(); //griglia
  set_style();
  
  gPad   -> SetTopMargin   ( 0.01  );
  gPad   -> SetBottomMargin( 0.09  );
  gPad   -> SetRightMargin ( 0.01  );
  gPad   -> SetLeftMargin  ( 0.075 );

  h->Draw();

  for( int i = 0; i < 3; ++i )
    gauss[i]->Draw("SAME");

  vector <TPaveText*> legend;
  string name;
  int k = 0;

  legend.push_back( new TPaveText( .05, .80, .20, .95, "NDC" ));
  legend.push_back( new TPaveText( .35, .80, .50, .95, "NDC" ));
  legend.push_back( new TPaveText( .70, .80, .85, .95, "NDC" ));

  for( int i = 0; i < 3; ++i ){
    name = "c = "       + to_string(  par[k] ) +
           " #pm "      + to_string( epar[k] );
    legend[i]->AddText( name.c_str() );
    k++;

    name = "#mu = "     + to_string(  par[k] ) +
           " #pm "      + to_string( epar[k] );
    legend[i]->AddText( name.c_str() );
    k++;

    name = "#sigma = "  + to_string(  par[k] ) +
           " #pm "      + to_string( epar[k] );
    legend[i]->AddText( name.c_str() );
    k++;
    
    legend[i]->Draw( );
  }


  cout << endl << endl
       << "----------------by aidin attar & Co----------------" 
       << endl << endl;

  return;
}
