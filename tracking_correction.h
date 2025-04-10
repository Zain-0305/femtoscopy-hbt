#include "call_libraries.h"  // call libraries from ROOT and C++

/*
Check if the given pT and eta values are within the valid bounds
--> Arguments:
pt: transverse momentum of the track
eta: pseudorapidity of the track
*/
bool checkBounds(double pt, double eta){
  if( TMath::Abs(eta) > 2.4 ){return false;}  // eta cut
  if( pt < 0 || pt > 500 ){return false;}     // pT cut
  return true;
}

/*
Get the track correction weight based on the tracking efficiency
--> Arguments:
trkeff_file: the input file containing tracking efficiency histograms
centrality: centrality bin to be used for the correction
pT: transverse momentum of the track
eta: pseudorapidity of the track
*/
double getTrkCorrWeight(TFile *trkeff_file, int centrality, double pT, double eta){
  if( !checkBounds(pT, eta) ) return 0;  // If track is out of bounds, return 0
  
  double factor = 1.0;
  TH2 *eff_factor = nullptr;  // Placeholder for the 2D efficiency histogram

  // Select the corresponding efficiency histogram based on centrality range
  if(centrality <= 20){
    trkeff_file->GetObject("rTotalEff3D_0_10", eff_factor);  // Centrality 0-10%
  }else if (centrality > 20 && centrality <= 60){
    trkeff_file->GetObject("rTotalEff3D_10_30", eff_factor);  // Centrality 10-30%
  }else if (centrality > 60 && centrality <= 100){
    trkeff_file->GetObject("rTotalEff3D_30_50", eff_factor);  // Centrality 30-50%
  }else if (centrality > 100 && centrality <= 140){
    trkeff_file->GetObject("rTotalEff3D_50_70", eff_factor);  // Centrality 50-70%
  }else if (centrality > 140){
    trkeff_file->GetObject("rTotalEff3D_70_100", eff_factor);  // Centrality 70-100%
  }

  // Get the efficiency value from the histogram
  double eff = eff_factor->GetBinContent( eff_factor->GetXaxis()->FindBin(eta), eff_factor->GetYaxis()->FindBin(pT) );

  // Calculate the tracking efficiency correction factor (inverse of efficiency)
  factor = (1. / eff);  // Return inverse of the efficiency
  
  return factor;
}
