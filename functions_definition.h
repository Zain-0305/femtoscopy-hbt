#ifndef FUNCTION_DEFINITIONS_H
#define FUNCTION_DEFINITIONS_H

#include "call_libraries.h"
#include <vector>
#include <utility>

// HBT Constants ==============================================================
constexpr double PI_MASS = 0.13957039;       // Pion mass (GeV/c²)
constexpr double COS_CUT = 0.99996;          // Cosine angle cut for split pairs
constexpr double DPT_CUT = 0.04;             // ΔpT cut for split pairs (GeV/c)
constexpr double ALPHA_EM = 1.0/137.035999;  // Fine structure constant

// HBT Quality Cuts ===========================================================
struct HBTTrackCut {
    float minPt = 0.15;       // GeV/c
    float maxEta = 2.4;       // Pseudorapidity
    int minHits = 11;         // Minimum tracker hits
    int minPixelHits = 2;     // Minimum pixel hits
    float maxDcaXY = 3.0;     // cm
    float maxDcaZ = 3.0;      // cm
    float maxChi2 = 5.0;      // χ²/ndof
    bool requireHighPurity = true;
};

// Core HBT Functions =========================================================

/**
 * Calculates invariant q for pion pairs
 * @param p1 First particle 4-vector
 * @param p2 Second particle 4-vector
 * @return q invariant (GeV/c)
 */
template<typename LorentzVec>
double CalculateQinv(const LorentzVec& p1, const LorentzVec& p2) {
    auto sum = p1 + p2;
    double q = sum.M2() - 4.0 * PI_MASS * PI_MASS;
    return (q > 0) ? std::sqrt(q) : -std::sqrt(-q);
}

/**
 * Calculates q_long in LCMS frame
 * @param p1 First particle 4-vector
 * @param p2 Second particle 4-vector
 * @return q_long (GeV/c)
 */
template<typename LorentzVec>
double CalculateQlongLCMS(const LorentzVec& p1, const LorentzVec& p2) {
    double num = 2.0 * (p1.Pz()*p2.E() - p2.Pz()*p1.E());
    double den = std::hypot(p1.E() + p2.E(), p1.Pz() + p2.Pz());
    return (den > 0) ? std::abs(num/den) : 0.0;
}

/**
 * Calculates q_out component
 * @param p1 First particle 4-vector
 * @param p2 Second particle 4-vector
 * @return q_out (GeV/c)
 */
template<typename LorentzVec>
double CalculateQout(const LorentzVec& p1, const LorentzVec& p2) {
    TVector3 qT(p1.Px()-p2.Px(), p1.Py()-p2.Py(), 0);
    TVector3 kT((p1.Px()+p2.Px())/2.0, (p1.Py()+p2.Py())/2.0, 0);
    return qT.Dot(kT.Unit());
}

/**
 * Calculates q_side component
 * @param p1 First particle 4-vector
 * @param p2 Second particle 4-vector
 * @return q_side (GeV/c)
 */
template<typename LorentzVec>
double CalculateQside(const LorentzVec& p1, const LorentzVec& p2) {
    TVector3 qT(p1.Px()-p2.Px(), p1.Py()-p2.Py(), 0);
    TVector3 kT((p1.Px()+p2.Px())/2.0, (p1.Py()+p2.Py())/2.0, 0);
    TVector3 qout = qT.Dot(kT.Unit()) * kT.Unit();
    return (qT - qout).Mag();
}

/**
 * Gamow factor for same-sign pion pairs
 * @param qinv q invariant (GeV/c)
 * @param syst Systematic variation (0=nominal, 1=+15%, 2=-15%)
 */
double CoulombSSWeight(double qinv, int syst = 0) {
    if (qinv < 1e-6) return 1.0;
    double x = 2.0 * M_PI * ALPHA_EM * PI_MASS / qinv;
    double weight = 1.0 + (syst == 1 ? 0.15 : (syst == 2 ? -0.15 : 0.0));
    return weight * ((std::exp(x) - 1.0)/x - 1.0) + 1.0;
}

/**
 * Gamow factor for opposite-sign pion pairs
 * @param qinv q invariant (GeV/c)
 * @param syst Systematic variation (0=nominal, 1=+15%, 2=-15%)
 */
double CoulombOSWeight(double qinv, int syst = 0) {
    if (qinv < 1e-6) return 1.0;
    double x = 2.0 * M_PI * ALPHA_EM * PI_MASS / qinv;
    double weight = 1.0 + (syst == 1 ? 0.15 : (syst == 2 ? -0.15 : 0.0));
    return weight * ((1.0 - std::exp(-x))/x - 1.0) + 1.0;
}

// Track Selection ============================================================

/**
 * Checks if track pair should be split (split/merge protection)
 * @param p1 First track 4-vector
 * @param p2 Second track 4-vector
 * @return true if pair should be split
 */
template<typename LorentzVec>
bool IsSplitPair(const LorentzVec& p1, const LorentzVec& p2) {
    double cosa = (p1.Px()*p2.Px() + p1.Py()*p2.Py() + p1.Pz()*p2.Pz()) / (p1.P() * p2.P());
    double dpt = std::abs(p1.Pt() - p2.Pt());
    return (cosa > COS_CUT) && (dpt < DPT_CUT);
}

/**
 * Counts tracks passing HBT analysis cuts
 * @param tracks Vector of track 4-vectors
 * @param cuts Track selection criteria
 * @return Number of good tracks
 */
template<typename LorentzVec>
int CountHBTTracks(const std::vector<LorentzVec>& tracks, const HBTTrackCut& cuts) {
    return std::count_if(tracks.begin(), tracks.end(), [&](const LorentzVec& trk) {
        return (trk.Pt() > cuts.minPt) && 
               (std::abs(trk.Eta()) < cuts.maxEta);
        // Add additional cuts as needed
    });
}

// Correlation Analysis =======================================================

/**
 * Fills HBT correlation histograms for an event
 * @param tracks Track 4-vectors
 * @param charges Track charges
 * @param weights Track weights
 * @param hSame Same-sign histograms
 * @param hOpp Opposite-sign histograms
 * @param cent Centrality bin
 * @param applyCoulomb Whether to apply Gamow correction
 * @param syst Systematic variation
 */
template<typename LorentzVec>
void AnalyzeHBTCorrelations(
    const std::vector<LorentzVec>& tracks,
    const std::vector<int>& charges,
    const std::vector<double>& weights,
    HBTHistograms& hSame,
    HBTHistograms& hOpp,
    int cent,
    bool applyCoulomb = true,
    int syst = 0) 
{
    const int nTracks = tracks.size();
    
    for (int i = 0; i < nTracks; ++i) {
        for (int j = i + 1; j < nTracks; ++j) {
            // Skip split pairs
            if (IsSplitPair(tracks[i], tracks[j])) continue;
            
            double weight = weights[i] * weights[j];
            double qinv = CalculateQinv(tracks[i], tracks[j]);
            bool isSameSign = (charges[i] * charges[j] > 0);
            
            // Apply Coulomb correction
            if (applyCoulomb) {
                weight *= isSameSign ? CoulombSSWeight(qinv, syst) 
                                     : CoulombOSWeight(qinv, syst);
            }
            
            // Fill appropriate histograms
            auto& h = isSameSign ? hSame : hOpp;
            h.FillCorrelation(tracks[i], tracks[j], cent, weight);
        }
    }
}

// Utility Functions ==========================================================

/**
 * Inverts momentum vector for mixed-event analysis
 * @param vec Input 4-vector
 * @return Inverted 4-vector
 */
template<typename LorentzVec>
LorentzVec InvertMomentum(const LorentzVec& vec) {
    LorentzVec inverted = vec;
    inverted.SetPxPyPzE(-vec.Px(), -vec.Py(), -vec.Pz(), vec.E());
    return inverted;
}

/**
 * Rotates XY components for background estimation
 * @param vec Input 4-vector
 * @return Rotated 4-vector
 */
template<typename LorentzVec>
LorentzVec RotateXY(const LorentzVec& vec) {
    LorentzVec rotated = vec;
    rotated.SetXYZT(-vec.X(), -vec.Y(), vec.Z(), vec.T());
    return rotated;
}

#endif // FUNCTION_DEFINITIONS_H
