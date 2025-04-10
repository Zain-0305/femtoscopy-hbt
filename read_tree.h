#include "call_libraries.h"  // Call libraries from ROOT and C++

// Declare variables

// Event quantities
float vertexz;               // Event z vertex
int hiBin;                   // Event centrality (centrality bin)
float hfplus;                // Event HF + energy deposit
float hfminus;               // Event HF - energy deposit
int PrimaryVertexFilter;     // Primary vertex filter flag
int BeamScrapingFilter;      // Beam scraping filter flag
int HfCoincFilter;           // HF coincidence filter flag

// Reconstructed tracks
int ntrk;                    // Number of tracks
float trkpt[9999];           // Track pT
float trketa[9999];          // Track eta
float trkphi[9999];          // Track phi
float trkpterr[9999];        // Track pT error (uncertainty)
float trkdcaxy[9999];        // Track dxy impact parameter
float trkdcaz[9999];         // Track dz impact parameter
float trkdcaxyerr[9999];     // Track dxy error (uncertainty)
float trkdcazerr[9999];      // Track dz error (uncertainty)
float trkchi2[9999];         // Track reconstruction chi2
float pfEcal[9999];          // Particle flow energy deposit in ECAL
float pfHcal[9999];          // Particle flow energy deposit in HCAL
float trkmva[9999];          // Track MVA for each step
int trkalgo[9999];           // Track algorithm/step
unsigned char trkndof[9999]; // Track number of degrees of freedom in the fitting
int trkcharge[9999];         // Track charge
unsigned char trknhits[9999];// Number of hits in the tracker
unsigned char trknlayer[9999]; // Number of layers with measurement in the tracker
unsigned char trkpixhits[9999]; // Number of hits in the pixel detector
bool highpur[9999];          // Tracker steps MVA selection flag

// Generated tracks (MC-specific)
std::vector<float> *gen_trkpt = 0;    // Gen particle pT
std::vector<float> *gen_trketa = 0;   // Gen particle eta
std::vector<float> *gen_trkphi = 0;   // Gen particle phi
std::vector<int> *gen_trkchg = 0;     // Gen particle charge
std::vector<int> *gen_trkpid = 0;     // Gen particle PID

// All variables listed above are read in the function below
/*
Function to read the Forest/Skim tree.
Arguments -> Transfer quantities from trees to our variables.
tree: input TChain from jet_analyzer.C file
is_MC: true -> MC, false -> Data
*/
void read_tree(TChain *tree, bool is_MC){

    tree->SetBranchStatus("*", 0); // Disable all branches initially

    // Enable branches of interest (based on the variables defined above)
    // Event quantities
    tree->SetBranchStatus("vz", 1);
    tree->SetBranchAddress("vz", &vertexz);
    
    tree->SetBranchStatus("hiBin", 1); 
    tree->SetBranchAddress("hiBin", &hiBin);
    
    tree->SetBranchStatus("hiHFplus", 1); 
    tree->SetBranchAddress("hiHFplus", &hfplus);
    
    tree->SetBranchStatus("hiHFminus", 1); 
    tree->SetBranchAddress("hiHFminus", &hfminus);
    
    tree->SetBranchStatus("pPAprimaryVertexFilter", 1); 
    tree->SetBranchAddress("pPAprimaryVertexFilter", &PrimaryVertexFilter); 
    
    tree->SetBranchStatus("pBeamScrapingFilter", 1); 
    tree->SetBranchAddress("pBeamScrapingFilter", &BeamScrapingFilter);
    
    tree->SetBranchStatus("phfCoincFilter3", 1);
    tree->SetBranchAddress("phfCoincFilter3", &HfCoincFilter);

    // Track quantities
    tree->SetBranchStatus("nTrk", 1);
    tree->SetBranchStatus("trkPt", 1);
    tree->SetBranchStatus("trkEta", 1);
    tree->SetBranchStatus("trkPhi", 1);
    tree->SetBranchStatus("trkPtError", 1);
    tree->SetBranchStatus("trkDxy1", 1);
    tree->SetBranchStatus("trkDxyError1", 1);
    tree->SetBranchStatus("trkDz1", 1);
    tree->SetBranchStatus("trkDzError1", 1);
    tree->SetBranchStatus("trkChi2", 1);
    tree->SetBranchStatus("trkNdof", 1);
    tree->SetBranchStatus("trkCharge", 1);
    tree->SetBranchStatus("trkNHit", 1);
    tree->SetBranchStatus("trkNlayer", 1);
    tree->SetBranchStatus("trkNPixelHit", 1);
    tree->SetBranchStatus("highPurity", 1);
    tree->SetBranchStatus("pfEcal", 1);
    tree->SetBranchStatus("pfHcal", 1);

    tree->SetBranchAddress("nTrk", &ntrk);
    tree->SetBranchAddress("trkPt", &trkpt);
    tree->SetBranchAddress("trkEta", &trketa);
    tree->SetBranchAddress("trkPhi", &trkphi);
    tree->SetBranchAddress("trkPtError", &trkpterr);
    tree->SetBranchAddress("trkDxy1", &trkdcaxy);
    tree->SetBranchAddress("trkDxyError1", &trkdcaxyerr);
    tree->SetBranchAddress("trkDz1", &trkdcaz);
    tree->SetBranchAddress("trkDzError1", &trkdcazerr);
    tree->SetBranchAddress("trkChi2", &trkchi2);
    tree->SetBranchAddress("trkNdof", &trkndof);
    tree->SetBranchAddress("trkCharge", &trkcharge);
    tree->SetBranchAddress("trkNHit", &trknhits);
    tree->SetBranchAddress("trkNlayer", &trknlayer);
    tree->SetBranchAddress("trkNPixelHit", &trkpixhits);
    tree->SetBranchAddress("highPurity", &highpur);
    tree->SetBranchAddress("pfEcal", &pfEcal);
    tree->SetBranchAddress("pfHcal", &pfHcal);

    // Gen particle quantities (only for MC)
    if(is_MC){
        tree->SetBranchStatus("pt", 1);
        tree->SetBranchStatus("eta", 1);
        tree->SetBranchStatus("phi", 1);
        tree->SetBranchStatus("chg", 1);
        tree->SetBranchStatus("pdg", 1);

        tree->SetBranchAddress("pt", &gen_trkpt);
        tree->SetBranchAddress("eta", &gen_trketa);
        tree->SetBranchAddress("phi", &gen_trkphi);
        tree->SetBranchAddress("chg", &gen_trkchg);
        tree->SetBranchAddress("pdg", &gen_trkpid);
    }
}
