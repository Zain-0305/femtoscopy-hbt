#ifndef MIXEVENTSHBT_H
#define MIXEVENTSHBT_H

#include <vector>
#include <cmath>
#include <iostream>
#include "TVector2.h"
#include "TH1.h"
#include "THnSparse.h"
#include "Math/Vector4D.h"

// Use this if you apply ΔpT and cosθ pair rejection
bool splitcomb(ROOT::Math::PtEtaPhiMVector p1, ROOT::Math::PtEtaPhiMVector p2, double coscut, double dptcut);

// Main mixing function for HBT and correlation analysis
void MixEvents(
    bool use_centrality,
    int centrality_or_ntrkoff_int,
    int nEvt_to_mix,
    std::vector<int> ev_centrality,
    std::vector<int> ev_multiplicity,
    std::vector<double> vtx_z_vec,
    float vzcut,
    std::vector<std::vector<ROOT::Math::PtEtaPhiMVector>> Track_Vector,
    std::vector<std::vector<int>> Track_Chg_Vector,
    std::vector<std::vector<double>> Track_Eff_Vector,
    THnSparseD* histo_SS,
    THnSparseD* histo_SS3D,
    THnSparseD* histo_OS,
    THnSparseD* histo_OS3D,
    bool docostdptcut,
    bool do_hbt3d,
    bool dogamovcorrection,
    int systematic,
    TH1I* NeventsAss
);

#endif // MIXEVENTSHBT_H
