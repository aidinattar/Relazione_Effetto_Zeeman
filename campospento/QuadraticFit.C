#include<string>
#include<iostream>

int QuadraticFit(string filename, string title, string xname, string yname) {

  // crea l'oggetto GRAFICO CON ERRORI leggendo da file di testo
  TGraphErrors *g = new TGraphErrors(filename.c_str());
  // formattazione del grafico con comandi testuali
  g->SetTitle(title.c_str());
  g->GetXaxis()->SetTitle(xname.c_str());
  g->GetXaxis()->SetTitleSize(0.04);
  g->GetXaxis()->CenterTitle();
  g->GetYaxis()->SetTitle(yname.c_str());
  g->GetYaxis()->SetTitleSize(0.04);
  g->GetYaxis()->CenterTitle();

  TCanvas *c1 = new TCanvas("c1"); // apre la finestra grafica

  g->Draw("ap"); // disegna il grafico a punti (p) in un nuovo frame (a)
  g->SetMarkerStyle(20);

  TF1 *f1 = new TF1("f1","pol2",-10e15,10e15); // formula definita nel costruttore 

  // setta i valori iniziali dei parametri
  f1->SetParNames("a","b", "c");
  f1->SetLineColor(kBlue);

  TFitResultPtr r = g->Fit(f1,"RS"); // R = Region   S = StoreResults
  c1->SetGrid();
  gStyle->SetOptFit(); // cliccare sulla pad dei risultati con l'editor e provare a cambiare qualcosa...
  // gStyle->SetoptFit(1111); // stampa tutto il necessario

  double a = f1->GetParameter(0);
  double b = f1->GetParameter(1);
  double c = f1->GetParameter(2);

  // GRAFICO DEI RESIDUI
  TGraphErrors *gr = new TGraphErrors(filename.c_str()); // parte da una copia del grafico originale
  for (int i=0; i<g->GetN(); i++) {
    double res = g->GetY()[i] - f1->Eval(g->GetX()[i]); // residuo
    gr->SetPoint(i,g->GetX()[i],res);
    double eresy = g->GetEY()[i]; // contributo delle Yi
    double eresx = f1->Derivative(g->GetX()[i])*g->GetEX()[i]; // contrib. Xi (approx. 1 ordine)
    double eres = TMath::Sqrt(eresy*eresy+eresx*eresx);
    gr->SetPointError(i,0,eres);
  }
  gr->SetName("gr");
  // formattazione del grafico con comandi testuali
  gr->SetMarkerColor(4);
  gr->SetMarkerStyle(20);
  gr->SetTitle(" ");
  gr->GetXaxis()->SetTitle(xname.c_str());
  gr->GetXaxis()->SetTitleSize(0.04);
  gr->GetXaxis()->CenterTitle();
  gr->GetYaxis()->SetTitle("Residuals");
  gr->GetYaxis()->SetTitleSize(0.04);
  gr->GetYaxis()->CenterTitle();

  TCanvas *c2 = new TCanvas("c2");
  c2->SetGrid();
  gr->Draw("ap");

  TF1 *n = new TF1("n", "[0]", -1000000, 1000000);
  n->SetParameters(0);
  n->Draw("same");
  

  return 1;

}

