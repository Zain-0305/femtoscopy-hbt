#include "call_libraries.h"  // Call libraries from ROOT and C++
#define coscut 0.99996
#define dptcut 0.04
#define pimass 0.1396

/*
Find Ntrk offline -> updated for all systems (and easy to update for future systems)
The Ntrk offline is a definition with specific cuts (we should not change it). The track systematics must be applied using the input_variables.h!
--> Arguments
size: track collection size per event
pt: track pT
eta: track eta
charge: track charge
hp: track high purity workflow
pterr: track pT uncertainty
dcaxy: track DCA in the transverse plane
dcaxyerr: track DCA in the transverse plane uncertainty
dcaz: track DCA in the longitudinal plane
dcazerr: track DCA in the longitudinal plane uncertainty
*/
int get_Ntrkoff(int size, float *pt, float *eta, int *charge, bool *hp, float *pterr, float *dcaxy, float *dcaxyerr,  float *dcaz, float *dcazerr) {
    int Ntrk_off = 0;
    for(int ii = 0; ii < size; ii++) { 
        // Apply cuts on track quality
        if(pt[ii] <= 0.3) continue;
        if(fabs(eta[ii]) > 2.4) continue; 
        if(fabs(charge[ii]) == 0) continue;
        if(hp[ii] == false) continue;
        if(fabs(pterr[ii]/pt[ii]) >= 0.1) continue;
        if(fabs(dcaxy[ii]/dcaxyerr[ii]) >= 3.0) continue;
        if(fabs(dcaz[ii]/dcazerr[ii]) >= 3.0) continue;
        Ntrk_off = Ntrk_off + 1;
    }
    return Ntrk_off;
}

/*
Calculate q invariant
--> Arguments
p1: particle 1 4-vector
p2: particle 2 4-vector
cos_cut: cut used in the cosine of angle between the particles
dpt_cut: cut used in the difference in pT of the particles
*/
bool splitcomb(ROOT::Math::PtEtaPhiMVector &vec1, ROOT::Math::PtEtaPhiMVector &vec2, double cos_cut, double dpt_cut) {
    bool issplit = false;
    double cosa = TMath::Abs(vec1.Px() * vec2.Px() + vec1.Py() * vec2.Py() + vec1.Pz() * vec2.Pz()) / (vec1.P() * vec2.P());
    double deltapt = fabs(vec1.Pt() - vec2.Pt());
    if ((cosa > cos_cut) && (deltapt < dpt_cut)) { 
        issplit = true;
    }
    return issplit;
}

/*
Calculate q invariant
--> Arguments
p1: particle 1 4-vector
p2: particle 2 4-vector
*/
float GetQ(ROOT::Math::PtEtaPhiMVector &p1, ROOT::Math::PtEtaPhiMVector &p2) {
    ROOT::Math::PtEtaPhiMVector Sum4V = p1 + p2;
    Double_t q = Sum4V.M2() - 4.0 * p1.mass() * p2.mass();
    return (q > 0 ? TMath::Sqrt(q) : -TMath::Sqrt(-q));
}

/*
Calculate q long in the LCMS
--> Arguments
p1: particle 1 4-vector
p2: particle 2 4-vector
*/
float GetQlongLCMS(ROOT::Math::PtEtaPhiMVector &p1, ROOT::Math::PtEtaPhiMVector &p2) {
    Double_t num = 2 * ((p1.Pz()) * (p2.E()) - (p2.Pz()) * (p1.E()));
    Double_t den = TMath::Sqrt((p1.E() + p2.E()) * (p1.E() + p2.E()) - (p1.Pz() + p2.Pz()) * (p1.Pz() + p2.Pz()));
    Double_t qlongLCMS = 0.0;
    if (den != 0) qlongLCMS = fabs(num / den);
    return qlongLCMS;
}

/*
Calculate q out in the LCMS
--> Arguments
p1: particle 1 4-vector
p2: particle 2 4-vector
*/
float GetQout(ROOT::Math::PtEtaPhiMVector &p1, ROOT::Math::PtEtaPhiMVector &p2) {
    TVector3 qT;
    qT.SetXYZ(p1.Px() - p2.Px(), p1.Py() - p2.Py(), 0.0);
    TVector3 kT;
    kT.SetXYZ((p1.Px() + p2.Px()) / 2.0, (p1.Py() + p2.Py()) / 2.0, 0.0);
    TVector3 qout;
    qout = qT.Dot(kT.Unit()) * kT.Unit();
    Double_t absValue = qout.Mag();
    return absValue; 
}

/*
Calculate q side in the LCMS
--> Arguments
p1: particle 1 4-vector
p2: particle 2 4-vector
*/
float GetQside(ROOT::Math::PtEtaPhiMVector &p1, ROOT::Math::PtEtaPhiMVector &p2) {
    TVector3 qT;
    qT.SetXYZ(p1.Px() - p2.Px(), p1.Py() - p2.Py(), 0.0);
    TVector3 kT;
    kT.SetXYZ((p1.Px() + p2.Px()) / 2.0, (p1.Py() + p2.Py()) / 2.0, 0.0);
    TVector3 qout;
    qout = qT.Dot(kT.Unit()) * kT.Unit();
    TVector3 qsid;
    qsid = qT - qout;
    Double_t absValue = qsid.Mag();
    return absValue;
}

/*
Invert Px, Py and Pz for reference sample
--> Arguments
p1: particle 1 4-vector
p2: particle 2 4-vector
*/
ROOT::Math::PtEtaPhiMVector InvertPVector(ROOT::Math::PtEtaPhiMVector &vec) {
    ROOT::Math::PtEtaPhiMVector ovec = vec;
    ovec.SetPxPyPzE(-vec.Px(), -vec.Py(), -vec.Pz(), vec.E());
    return ovec;
}

/*
Rotate X and Y for reference sample
--> Arguments
p1: particle 1 4-vector
p2: particle 2 4-vector
*/
ROOT::Math::PtEtaPhiMVector InvertXYVector(ROOT::Math::PtEtaPhiMVector &vec) {
    ROOT::Math::PtEtaPhiMVector ovec = vec;
    ovec.SetXYZT(-vec.X(), -vec.Y(), vec.Z(), vec.T());
    return ovec;
}

/*
Return the weight factor due to Coulomb repulsion [Gamow] same charge
--> Arguments
q: q invariant
systematic: systematic number
*/
const double CoulombSS(const double &q, int systematic) {
    const double alpha = 1. / 137.0;
    double x = 2.0 * TMath::Pi() * (alpha * pimass / q);
    double weight = 1.0;
    if (systematic == 9) weight = 1.15; // For syst. +15%
    if (systematic == 10) weight = 0.85; // For syst. -15%
    return weight * ((TMath::Exp(x) - 1.0) / x - 1.0) + 1.0;
}

/*
Return the weight factor due to Coulomb attraction [Gamow] opposite charge
--> Arguments
q: q invariant
systematic: systematic number
*/
const double CoulombOS(const double &q, int systematic) {
    const double alpha = 1. / 137.0;
    double x = 2.0 * TMath::Pi() * (alpha * pimass / q);
    double weight = 1.0;
    if (systematic == 9) weight = 1.15; // For syst. +15%
    if (systematic == 10) weight = 0.85; // For syst. -15%
    return weight * ((1. - TMath::Exp(-x)) / x - 1.0) + 1.0;
}

