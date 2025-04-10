#include "call_libraries.h"  // call libraries from ROOT and C++
#include "MixEventsHBT.h"    // Mixing events header file

// define the bins
// qinv
const int nQBins = 200;   // number of qinv bins
const double minQ = 0.0;  // minimum qinv
const double maxQ = 2.0;  // maximum qinv

// q3D
const int nQBins3D = 100;   // number of q3D bins
const double minQ3D = 0.0;  // minimum q3D
const double maxQ3D = 2.0;  // maximum q3D

// kT (adjusted for finer granularity)
const int nKtBins = 9; // number of average transverse momentum bins
double KtBins[nKtBins + 1] = {0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 1.0, 1.5}; 

// centrality (adjusted for finer granularity)
const int nCentBins = 13; // number of centrality bins
double CentBins[nCentBins + 1] = {0.0, 5.0, 10.0, 15.0, 20.0, 25.0, 30.0, 40.0, 60.0, 80.0, 100.0, 120.0, 150.0, 200.0}; 
// Centrality bins: 0-5%, 5-10%, 10-15%, 15-20%, 20-25%, 25-30%, 30-40%, 40-60%, 60-80%, 80-100%, 100-120%, 120-150%, 150-200%

// Event histograms
TH1I *Nevents = new TH1I("Nevents", "Nevents", 10, 0, 10);
TH1D *centrality_beforefilters = new TH1D("centrality_beforefilters", "centrality_beforefilters", 150, 0.0, 300.0);
TH1D *centrality = new TH1D("centrality", "centrality", 150, 0.0, 300.0);
TH1D *vzhist_beforefilters = new TH1D("vzhist_beforefilters", "vzhist_beforefilters", 80, -20., 20.);
TH1D *vzhist = new TH1D("vzhist", "vzhist", 80, -20., 20.);
TH1D *multiplicity = new TH1D("multiplicity", "multiplicity", 400, 0.0, 4000.0);
TH2D *MultVSCent = new TH2D("MultVSCent", "MultVSCent", 200, 0.0, 4000.0, 100, 0.0, 200.0);
TH1I *NeventsAss = new TH1I("NeventsAss", "NeventsAss", 11, 0, 11);
TH1I *NeventsAssGEN = new TH1I("NeventsAssGEN", "NeventsAssGEN", 11, 0, 11);
TH1D *CheckNtrk = new TH1D("CheckNtrk", "CheckNtrk", 20000, 0, 20000);


//histograms before selection
TH1D *dxyoversigmadxy_beforeselection = new TH1D("dxyoversigmadxy_beforeselection", "dxyoversigmadxy_beforeselection", 100, -6.0, 6.0);
TH1D *dzoversigmadz_beforeselection = new TH1D("dzoversigmadz_beforeselection", "dzoversigmadz_beforeselection", 100, -6.0, 6.0);
TH1D *ptresolution_beforeselection = new TH1D("ptresolution_beforeselection", "ptresolution_beforeselection", 50, 0.0, 0.25);
TH1D *chi2overNDFoverNLayer_beforeselection = new TH1D("chi2overNDFoverNLayer_beforeselection", "chi2overNDFoverNLayer_beforeselection", 100, 0.0, 0.5);
TH1D *nhits_beforeselection = new TH1D("nhits_beforeselection", "nhits_beforeselection", 60, 0.0, 60.0);
TH1D *npixelhit_beforeselection = new TH1D("npixelhit_beforeselection", "npixelhit_beforeselection", 5, 0.0, 5.0);

//histograms after selections
TH1D *dxyoversigmadxy = new TH1D("dxyoversigmadxy", "dxyoversigmadxy", 100, -6.0, 6.0);
TH1D *dzoversigmadz = new TH1D("dzoversigmadz", "dzoversigmadz", 100, -6.0, 6.0);
TH1D *ptresolution = new TH1D("ptresolution", "ptresolution", 50, 0.0, 0.25);
TH1D *chi2overNDFoverNLayer = new TH1D("chi2overNDFoverNLayer", "chi2overNDFoverNLayer", 100, 0.0, 0.5);
TH1D *nhits = new TH1D("nhits", "nhits", 60, 0.0, 60.0);
TH1D *npixelhit = new TH1D("npixelhit", "npixelhit", 5, 0.0, 5.0);


// Track/Particle histograms
// Axis : 0 -> track pT, 1 -> trk eta, 2 -> trk phi, 3 -> trk charge, 4 -> centrality bin
int	bins_trk[5]      =   { 200   ,  24  ,   30				      , 3   , nCentBins};
double xmin_trk[5]   =   { 0.0   , -2.4 ,   -TMath::Pi()  		  , -1.5, CentBins[0]};
double xmax_trk[5]   =   { 50.0  ,  2.4 ,   TMath::Pi()  		  ,  1.5, CentBins[nCentBins]};

