// correlation_XeXe.C - Unified HBT Femtoscopy Analysis for XeXe Collisions
// Maintains original function signature for compatibility

#include "call_libraries.h"
#include "hbt_analysis_utils.h"  // Contains common utilities
#include "track_corrections.h"   // Your improved tracking corrections

void correlation_XeXe(
    TString input_file,          // List of input files
    TString output_tag,          // Output file identifier
    int isMC = 0,                // 0=data, 1=MC
    int quick_test = 0,          // 0=full run, 1=test mode
    int mixing_mode = 1,         // 0=no mixing, 1=standard mixing
    int n_mix_events = 10,       // Number of events to mix
    int cent_mult_window = 5,    // Centrality/multiplicity window
    float vz_window = 2.0,       // Vertex z window (cm)
    int hbt3d = 1,               // 0=3D analysis, 1=1D only
    int coulomb_corr = 0,        // 0=no Coulomb, 1=apply correction
    int centrality_mode = 0,     // 0=centrality, 1=multiplicity
    int systematic = 0           // Systematic variation
) {
    // Start timing and logging
    TStopwatch timer;
    timer.Start();
    std::cout << "=== Starting HBT analysis for XeXe collisions ===" << std::endl;
    
    // ======================
    // 1. Configuration Setup
    // ======================
    const bool is_mc = (isMC == 1);
    const bool do_quick_test = (quick_test == 1);
    const bool do_mixing = (mixing_mode == 1);
    const bool do_3d = (hbt3d == 0);
    const bool do_coulomb = (coulomb_corr == 1);
    const bool use_cent = (centrality_mode == 0);
    
    // Systematic configuration
    TString syst_tag = GetSystematicTag(systematic);
    if (systematic == 9 || systematic == 10) do_coulomb = true;
    
    // ======================
    // 2. Initialize Components
    // ======================
    
    // a) Track correction setup
    TFile* eff_file = OpenEfficiencyFile(systematic);
    std::vector<TH2D*> eff_hists = LoadEfficiencyHists(eff_file);
    
    // b) Output file naming
    TDatime date;
    TString output_name = Form("%s_%s_Nmix%d_MixWin%d_VzWin%.1f_%d",
                              output_tag.Data(),
                              syst_tag.Data(),
                              n_mix_events,
                              cent_mult_window,
                              vz_window,
                              date.GetDate());
    
    // Replace special characters
    output_name.ReplaceAll(".", "p");
    
    // c) Create output file structure
    TFile output(output_name + ".root", "RECREATE");
    SetupOutputDirectories(output);
    
    // ======================
    // 3. Event Processing
    // ======================
    
    // a) Initialize chains
    TChain* event_chain = new TChain("hiEvtAnalyzer/HiTree");
    TChain* track_chain = new TChain("ppTrack/trackTree");
    TChain* skim_chain = new TChain("skimanalysis/HltTree");
    
    // b) Add input files
    AddInputFiles(input_file, {event_chain, track_chain, skim_chain});
    
    // c) Main event loop
    int n_events = event_chain->GetEntries();
    int n_processed = 0;
    
    for (int i = 0; i < n_events; ++i) {
        // Event selection
        if (!LoadEvent(event_chain, i)) continue;
        if (!PassEventCuts(systematic)) continue;
        
        // Track processing
        auto [tracks, weights, charges] = ProcessTracks(track_chain, eff_hists, systematic);
        
        // Pair analysis
        AnalyzePairs(tracks, weights, charges, do_3d, do_coulomb, systematic);
        
        // Store for mixing
        if (do_mixing) {
            StoreForMixing(tracks, weights, charges);
        }
        
        // Progress reporting
        if (++n_processed % 1000 == 0) {
            std::cout << "Processed " << n_processed << " events (" 
                      << 100.0*n_processed/n_events << "%)" << std::endl;
        }
        
        if (do_quick_test && n_processed >= 1000) break;
    }
    
    // ======================
    // 4. Event Mixing
    // ======================
    if (do_mixing) {
        std::cout << "Performing event mixing..." << std::endl;
        PerformMixing(n_mix_events, cent_mult_window, vz_window);
    }
    
    // ======================
    // 5. Finalization
    // ======================
    
    // Write results
    WriteResults(output);
    output.Close();
    
    // Cleanup
    delete event_chain;
    delete track_chain;
    delete skim_chain;
    
    // Timing information
    timer.Stop();
    std::cout << "=== Analysis completed ===" << std::endl;
    std::cout << "Processed " << n_processed << " events in " 
              << timer.RealTime() << " seconds" << std::endl;
    std::cout << "Output saved to: " << output_name << ".root" << std::endl;
}

// Helper function implementations
namespace {

TString GetSystematicTag(int syst) {
    static const std::map<int, TString> syst_map = {
        {0, "nominal"}, {1, "vznarrow"}, {2, "vzwide"},
        {3, "trktight"}, {4, "trkloose"}, {5, "centup"},
        {6, "centdown"}, {7, "nodupcut"}, {8, "nopixcut"},
        {9, "coulombp15"}, {10, "coulombm15"}
    };
    return syst_map.count(syst) ? syst_map.at(syst) : "unknown";
}

TFile* OpenEfficiencyFile(int syst) {
    TString eff_path = "efftables/XeXe_eff_table_94x_cent.root";
    if (syst == 1) eff_path = "efftables/XeXe_eff_narrow_table_94x_cent.root";
    if (syst == 2) eff_path = "efftables/XeXe_eff_wide_table_94x_cent.root";
    if (syst == 3) eff_path = "efftables/XeXe_eff_tight_table_94x_cent.root";
    if (syst == 4) eff_path = "efftables/XeXe_eff_loose_table_94x_cent.root";
    
    TFile* file = TFile::Open(eff_path);
    if (!file || file->IsZombie()) {
        throw std::runtime_error("Could not open efficiency file: " + eff_path);
    }
    return file;
}

} // anonymous namespace
