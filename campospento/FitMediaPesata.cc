#include "TROOT.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TGraph.h"

#include <iostream>
#include <fstream>
#include <string>
#include <cmath>


void set_style(){
    // setup graphics
    //gStyle -> SetOptTitle   (   kFALSE   );
    gStyle -> SetOptStat    (      0     );
    gStyle -> SetLabelOffset( 0.01 , "x" );
    gStyle -> SetLabelOffset( 0.005, "y" );
    gStyle -> SetTitleOffset( 1.2  , "x" );
    gStyle -> SetTitleOffset( 1.1  , "y" );
}


double MEA ( TGraphErrors* g ){
	double med, sumx=0, sum1e=0;
	for( int i = 0; i < g->GetN(); ++i ){
		sumx  += g->GetY()[i] / ( pow( g->GetEY()[i], 2 ));
		sum1e += 1            / ( pow( g->GetEY()[i], 2 ));
	}
	med = sumx / sum1e;
	return med;
}


double RMS ( TGraphErrors* g ){
	double error, sum1e=0;
	for( int i = 0; i < g->GetN(); ++i )
		sum1e += 1 / ( pow( g->GetEY()[i], 2 ));

	return sqrt( 1 / sum1e );
}


void FitMediaPesata( void ){

  TCanvas *c1 = new TCanvas("c1", "My ROOT Plots", 1280, 720);
  c1 -> SetGrid(); //griglia
  set_style();

  gPad   -> SetTopMargin   ( 0.01  );
  gPad   -> SetBottomMargin( 0.09  );
  gPad   -> SetRightMargin ( 0.01  );
  gPad   -> SetLeftMargin  ( 0.075 );

  string title, asx, asy;
  std::cout << "file con i dati sperimentali: ";
  std::getline( std::cin, title );
  TGraphErrors *g = new TGraphErrors( title.c_str() );

  std::cout << "asse x: ";
  std::getline( std::cin, asx );
  g -> GetXaxis() -> SetTitle( asx.c_str() );  
  g -> GetXaxis() -> CenterTitle();

  std::cout << "asse y: ";
  std::getline( std::cin, asy );
  g -> GetYaxis() -> SetTitle( asy.c_str() );
  g -> GetYaxis() -> CenterTitle();

  g -> SetMarkerStyle( 20 );

  // disegna il grafico a punti (p) in un nuovo frame (a)
  g->Draw("ap");

  //TF1 *f = new TF1 ( "Retta", "pol1", g -> GetMinimum(), g -> GetMaximum() );
  TF1 *f = new TF1 ( "Retta", "pol1", -99999, 99999 );
  f -> SetParNames( "offset", "slope" );

  f -> SetLineColor(kBlue);

  // R = Region   S = StoreResults
  TFitResultPtr r = g->Fit( f, "RS" ); 

  gStyle -> SetOptFit();

  double  A = f -> GetParameter( 0 );
  double eA = f ->  GetParError( 0 );
  double  B = f -> GetParameter( 1 );
  double eB = f ->  GetParError( 1 );

  std::cout << "A:\t "                               <<   A 
            << "\t +/- \t"                           <<  eA
            << "\t \t \t " 
            << round( fabs( eA / A * 10000 ) ) / 100 << "%"
            << std::endl;

  std::cout << "B:\t "                               <<   B
            << "\t +/- \t"                           <<  eB
            << "\t \t \t " 
            << round( fabs( eB / B * 10000 ) ) / 100 << "%"
            << std::endl << std::endl;


  // to access the covariance matrix
  TMatrixDSym  cov = r ->  GetCovarianceMatrix(); 
  // to retrieve the fit chi2
  double      chi2 = r ->                 Chi2(); 
  // to retrieve the fit rho
  double       rho = g -> GetCorrelationFactor(); 
  // to evaluate the fit t
  double         t = rho * sqrt( ( f -> GetNDF() ) 
                   / ( 1 - pow( rho, 2 ) ) );
  // to evaluate the chi compatibility
  double      comp = ( fabs( chi2 - f -> GetNDF() ) )
                   / sqrt( 2 * f -> GetNDF() );
  // to evaluate the fit maximum a posteriori
  double post = 0;
  for( unsigned int i = 0; i < g -> GetN(); ++i )
    post += pow( A + B * g -> GetX()[ i ] - g -> GetY()[ i ], 2 );
  post = sqrt( post / f -> GetNDF() );

  std::cout << "\nCOVARIANCE MATRIX :\n ";
  cov.Print();

  std::cout << "chi2:\t\t"  << chi2 << std::endl;
  std::cout << "rho:\t\t"   << rho  << std::endl;
  std::cout << "t: \t\t"    <<  t   << std::endl;
  std::cout << "comp:\t\t"  << comp << std::endl;
  std::cout << "sigma:\t\t" << post << std::endl;

  TF1 *h = new TF1 ( "h", "[0]", -99999, 99999 );
  h->SetParNames( "mean" );
  h->SetLineColor( kRed );

  double mea = MEA( g );
  double rms = RMS( g );

  h->SetParameter( 0, mea );
  h->Draw( "SAME" );

  std::cout << "mean:\t\t" << mea
            << " +/- "     << rms  << std::endl;

  TLegend *legend = new TLegend();
  legend->AddEntry( g, "experimental data", "lpf");
  legend->AddEntry( f, "linear fit",         "l" );
  legend->AddEntry( h, "mean",               "l" );

  legend->Draw( "SAME" );

  TCanvas *c2 = new TCanvas( "c2" , "My ROOT Plots", 1280, 720 );
  c2 -> SetGrid(); //griglia
  set_style();

  gPad   -> SetTopMargin   ( 0.01  );
  gPad   -> SetBottomMargin( 0.09  );
  gPad   -> SetRightMargin ( 0.01  );
  gPad   -> SetLeftMargin  ( 0.075 );

  // parte da una copia del grafico originale
  TGraphErrors *gr = new TGraphErrors( title.c_str() ); 
  for ( unsigned int i = 0; i < g -> GetN(); i++ ){
    // residuo
    double   res = g -> GetY()[ i ] - f -> Eval(g -> GetX()[ i ]); 
    gr -> SetPoint( i, g -> GetX()[ i ], res );
    // contributo delle Yi
    double eresy = g -> GetEY()[ i ];
    // contrib. Xi (approx. 1 ordine)
    double eresx = f -> Derivative( g -> GetX()[ i ] ) *g -> GetEX()[ i ]; 
    double eres = TMath::Sqrt( eresy * eresy + eresx * eresx );
    gr -> SetPointError( i, 0, eres );
  }

  //TF1 *fr = new TF1( "zero", "[0]", gr -> GetMinimum(), gr -> GetMaximum() );
  TF1 *fr = new TF1( "zero", "[0]", -999999, 999999 );

  fr      -> SetParameters( 0, 0 );
  fr      -> SetLineColor(  kRed );

  gr      ->       SetName( "gr" );
  gr      ->   SetMarkerColor( 4 );
  gr      ->  SetMarkerStyle( 20 );

  gr      ->  GetXaxis() -> SetTitle( asx.c_str() );
  gr      ->  GetYaxis() -> SetTitle( "Residuals" );
  gr      ->  GetXaxis() -> CenterTitle();
  gr      ->  GetYaxis() -> CenterTitle();

  TLegend *legendr = new TLegend();

  legendr -> AddEntry( gr, "Residuals", "lpf" );

  gr      -> Draw(  "AP"  );
  fr      -> Draw( "SAME" );
  legendr -> Draw( "SAME" );
}
