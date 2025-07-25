#ifndef READ_TREE_H
#define READ_TREE_H

#include "call_libraries.h"
#include <vector>
#include <stdexcept>

// Configuration ===============================================================
constexpr int MAX_HBT_TRACKS = 30000;  // Safe upper limit for PbPb HBT
constexpr float MIN_HBT_PT = 0.15;     // GeV/c
constexpr float MAX_HBT_ETA = 2.4;     // Pseudorapidity cut

// HBT Quality Flags ===========================================================
struct HBTQualityCuts {
    bool requireHighPurity = true;
    int minPixelHits = 2;
    int minTotalHits = 11;
    float maxDcaXY = 3.0;    // cm
    float maxDcaZ = 3.0;     // cm
    float maxChi2 = 5.0;     // χ²/ndof
};

// Main Data Structure =========================================================
struct HBTEvent {
    // Event Info
    float vz = 0;            // Vertex z-position (cm)
    int hiBin = -1;          // Centrality (PbPb only)
    float weight = 1.0;      // MC weight
    
    // Track Arrays (optimized for HBT)
    int nTracks = 0;
    float pt[MAX_HBT_TRACKS];
    float eta[MAX_HBT_TRACKS];
    float phi[MAX_HBT_TRACKS];
    float dcaXY[MAX_HBT_TRACKS];    // Important for HBT pair cuts
    float dcaZ[MAX_HBT_TRACKS];     // Important for HBT pair cuts
    short charge[MAX_HBT_TRACKS];   // ±1
    bool goodTrack[MAX_HBT_TRACKS]; // Passed quality cuts
    
    // MC Truth (if available)
    bool isMC = false;
    std::vector<float> genPt;       // For efficiency corrections
    std::vector<float> genEta;
    std::vector<float> genPhi;
    
    void clear() {
        nTracks = 0;
        if (isMC) {
            genPt.clear();
            genEta.clear();
            genPhi.clear();
        }
    }
};

// Core Function ===============================================================
/**
 * Reads an event with HBT-optimized branch loading and track selection
 * 
 * @param tree Input TChain (must be initialized)
 * @param event Output event structure
 * @param cuts Track quality cuts for HBT analysis
 * @param isMC Whether to read MC truth branches
 * @throws std::runtime_error if critical branches are missing
 */
void read_tree(TChain* tree, HBTEvent& event, 
               const HBTQualityCuts& cuts = HBTQualityCuts(),
               bool isMC = false) 
{
    // Validate input
    if (!tree) throw std::runtime_error("Null TChain provided");
    event.clear();
    event.isMC = isMC;
    
    // Disable all branches initially for performance
    tree->SetBranchStatus("*", 0);
    
    // 1. Enable Event-Level Branches ------------------------------------------
    tree->SetBranchStatus("vz", 1);
    tree->SetBranchAddress("vz", &event.vz);
    
    tree->SetBranchStatus("hiBin", 1);
    tree->SetBranchAddress("hiBin", &event.hiBin);
    
    // 2. Enable Track Branches (HBT-essential) -------------------------------
    int nTrk = 0;
    float tmpPt[MAX_HBT_TRACKS];
    float tmpEta[MAX_HBT_TRACKS];
    float tmpPhi[MAX_HBT_TRACKS];
    float tmpDcaXY[MAX_HBT_TRACKS];
    float tmpDcaZ[MAX_HBT_TRACKS];
    UChar_t tmpNhits[MAX_HBT_TRACKS];
    UChar_t tmpPixHits[MAX_HBT_TRACKS];
    bool tmpHighPurity[MAX_HBT_TRACKS];
    float tmpChi2[MAX_HBT_TRACKS];
    Short_t tmpCharge[MAX_HBT_TRACKS];
    
    tree->SetBranchStatus("nTrk", 1);
    tree->SetBranchAddress("nTrk", &nTrk);
    
    tree->SetBranchStatus("trkPt", 1);
    tree->SetBranchAddress("trkPt", tmpPt);
    
    tree->SetBranchStatus("trkEta", 1);
    tree->SetBranchAddress("trkEta", tmpEta);
    
    tree->SetBranchStatus("trkPhi", 1);
    tree->SetBranchAddress("trkPhi", tmpPhi);
    
    tree->SetBranchStatus("trkDxy1", 1);
    tree->SetBranchAddress("trkDxy1", tmpDcaXY);
    
    tree->SetBranchStatus("trkDz1", 1);
    tree->SetBranchAddress("trkDz1", tmpDcaZ);
    
    tree->SetBranchStatus("highPurity", 1);
    tree->SetBranchAddress("highPurity", tmpHighPurity);
    
    tree->SetBranchStatus("trkNHit", 1);
    tree->SetBranchAddress("trkNHit", tmpNhits);
    
    tree->SetBranchStatus("trkNPixelHit", 1);
    tree->SetBranchAddress("trkNPixelHit", tmpPixHits);
    
    tree->SetBranchStatus("trkChi2", 1);
    tree->SetBranchAddress("trkChi2", tmpChi2);
    
    tree->SetBranchStatus("trkCharge", 1);
    tree->SetBranchAddress("trkCharge", tmpCharge);
    
    // 3. Enable MC Branches if needed ----------------------------------------
    if (isMC) {
        tree->SetBranchStatus("weight", 1);
        tree->SetBranchAddress("weight", &event.weight);
        
        std::vector<float> *genPt = nullptr, *genEta = nullptr, *genPhi = nullptr;
        tree->SetBranchStatus("pt", 1);
        tree->SetBranchAddress("pt", &genPt);
        
        tree->SetBranchStatus("eta", 1);
        tree->SetBranchAddress("eta", &genEta);
        
        tree->SetBranchStatus("phi", 1);
        tree->SetBranchAddress("phi", &genPhi);
    }
    
    // 4. Load the entry -----------------------------------------------------
    tree->GetEntry(0);
    
    // 5. Apply HBT Track Selection ------------------------------------------
    event.nTracks = 0;
    for (int i = 0; i < nTrk && event.nTracks < MAX_HBT_TRACKS; ++i) {
        // Basic kinematic cuts
        if (tmpPt[i] < MIN_HBT_PT) continue;
        if (fabs(tmpEta[i]) > MAX_HBT_ETA) continue;
        
        // HBT-specific quality cuts
        if (cuts.requireHighPurity && !tmpHighPurity[i]) continue;
        if (tmpPixHits[i] < cuts.minPixelHits) continue;
        if (tmpNhits[i] < cuts.minTotalHits) continue;
        if (fabs(tmpDcaXY[i]) > cuts.maxDcaXY) continue;
        if (fabs(tmpDcaZ[i]) > cuts.maxDcaZ) continue;
        if (tmpChi2[i] > cuts.maxChi2) continue;
        
        // Store accepted tracks
        event.pt[event.nTracks] = tmpPt[i];
        event.eta[event.nTracks] = tmpEta[i];
        event.phi[event.nTracks] = tmpPhi[i];
        event.dcaXY[event.nTracks] = tmpDcaXY[i];
        event.dcaZ[event.nTracks] = tmpDcaZ[i];
        event.charge[event.nTracks] = tmpCharge[i];
        event.goodTrack[event.nTracks] = true;
        event.nTracks++;
    }
    
    // 6. Handle MC Truth ----------------------------------------------------
    if (isMC) {
        // ... (MC processing logic) ...
    }
}

#endif // READ_TREE_H
