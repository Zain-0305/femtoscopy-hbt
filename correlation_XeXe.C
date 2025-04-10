#include "call_libraries.h"      // call libraries from ROOT and C++
#include "read_tree.h"           // read the TChains
#include "define_histograms.h"   // histogram definition
#include "tracking_correction.h" // tracking correction
#define pi_mass 0.1396

/*
--> Arguments:
- input_file: text file with a list of root input files: Forest or Skims
- ouputfile: just a counting number to run on Condor
- isMC: 0 for false (data) and > 0 for true (MC)
- doquicktest: 0 for false and > 0 for true (tests with 1000 events)
- domixing: 0 for mixing and 1 without mixing
- Nmixevents: number of events to mix
- mincentormult: minimum centrality or multiplicity in the mixing
- minvz: minimum vertex z for the mixing
- hbt3d: 0 with 3D and 1 without 3D
- gamov: 0 means no GAMOV Coulomb correction, > 0 means GAMOV is added
- cent_bool: 0 to use centrality and different than 0 for multiplicity
- syst: systematics type:
  0: nominal, 1: |vz| < 3, 2: 3 < |vz| < 15, etc.
*/
void correlation_XeXe(TString input_file, TString ouputfile, int isMC, int doquicktest, int domixing, int Nmixevents, int mincentormult, float minvz, int hbt3d, int gamov, int cent_bool, int syst) {

    clock_t sec_start, sec_end;
    sec_start = clock(); // start timing measurement
    TDatime* date = new TDatime(); // to add date in the output file

    // Boolean conversion
    bool do_quicktest = (doquicktest > 0);
    bool is_MC = (isMC > 0);
    bool do_mixing = (domixing == 0);  // if 0, mixing enabled
    bool do_hbt3d = (hbt3d == 0);
    bool do_gamov = (gamov > 0);
    bool use_centrality = (cent_bool == 0);

    if (syst == 9 || syst == 10) do_gamov = true;  // Handle special cases for gamov

    bool dosplit = (syst != 7);  // Enable splitting unless systematics 7

    // Determine systematics string
    TString systematics = "nonapplied_nominal";
    switch (syst) {
        case 0: systematics = "nominal"; break;
        case 1: systematics = "vznarrow"; break;
        case 2: systematics = "vzwide"; break;
        case 3: systematics = "trktight"; break;
        case 4: systematics = "trkloose"; break;
        case 5: systematics = "centup"; break;
        case 6: systematics = "centdown"; break;
        case 7: systematics = "removeduplicatedcut"; break;
        case 8: systematics = "removeNpixelhitcut"; break;
        case 9: systematics = "gamovplus15"; break;
        case 10: systematics = "gamovminus15"; break;
    }
    systematics += (use_centrality ? "_cent" : "_Ntroff");

    // Open efficiency file based on systematics
    TFile *fileeff = nullptr;
    if (syst == 1) fileeff = TFile::Open("efftables/XeXe_eff_narrow_table_94x_cent.root");
    else if (syst == 2) fileeff = TFile::Open("efftables/XeXe_eff_wide_table_94x_cent.root");
    else if (syst == 3) fileeff = TFile::Open("efftables/XeXe_eff_tight_table_94x_cent.root");
    else if (syst == 4) fileeff = TFile::Open("efftables/XeXe_eff_loose_table_94x_cent.root");
    else fileeff = TFile::Open("efftables/XeXe_eff_table_94x_cent.root");

    if (!fileeff || fileeff->IsZombie()) {
        cout << "Error opening efficiency file!" << endl;
        return;
    }

    // Read the list of input file(s)
    fstream inputfile;
    inputfile.open(input_file.Data(), ios::in);
    if (!inputfile.is_open()) {
        cout << "List of input files not found!" << endl;
        return;
    }

    // Vector to store file names
    TString file_chain;
    std::vector<TString> file_name_vector;
    while (getline(inputfile, file_chain)) {
        file_name_vector.push_back(file_chain);
    }
    inputfile.close();

    // Create chains for various event and track data
    TChain *hea_tree = new TChain("hiEvtAnalyzer/HiTree");
    TChain *ski_tree = new TChain("skimanalysis/HltTree");
    TChain *trk_tree = new TChain("ppTrack/trackTree");
    TChain *gen_tree = is_MC ? new TChain("HiGenParticleAna/hi") : nullptr;

    // Add files to the chains
    for (auto &file : file_name_vector) {
        TFile *testfile = TFile::Open(file, "READ");
        if (testfile && !testfile->IsZombie()) {
            cout << "Adding file " << file << " to the chains" << endl;
            hea_tree->Add(file);
            ski_tree->Add(file);
            trk_tree->Add(file);
            if (is_MC) gen_tree->Add(file);
        } else {
            cout << "File: " << file << " failed!" << endl;
        }
    }

    // Clear file list after processing
    file_name_vector.clear();

    // Connect chains
    hea_tree->AddFriend(ski_tree);
    hea_tree->AddFriend(trk_tree);
    if (is_MC) hea_tree->AddFriend(gen_tree);

    // Read tree information
    read_tree(hea_tree, is_MC);

    // Initialize histograms with uncertainty calculation
    sw2();

    // Get number of events in the chain
    int nevents = hea_tree->GetEntries();
    cout << endl;
    cout << "Total number of events in those files: " << nevents << endl;
    cout << "-------------------------------------------------" << endl;

    // Vectors for mixing
    std::vector<int> centrality_vector, multiplicity_vector;
    std::vector<double> vz_vector;
    std::vector<std::vector<ROOT::Math::PtEtaPhiMVector>> track_4vector;
    std::vector<std::vector<double>> track_weights_vector;
    std::vector<std::vector<int>> track_charge_vector;

    // Event loop
    for (int i = 0; i < nevents; i++) {
        hea_tree->GetEntry(i);

        // Centrality and vertex selection
        int centrality = getCentrality();
        double vz = getVz();

        // Apply vertex cut
        if (vz < minvz) continue;

        // Process tracks and weights
        std::vector<ROOT::Math::PtEtaPhiMVector> tracks = getTracks();
        std::vector<double> track_weights = getTrackWeights();
        std::vector<int> track_charge = getTrackCharge();

        // Fill histograms
        fillHistograms(tracks, centrality, vz);

        // Apply mixing if enabled
        if (do_mixing) {
            mixEvents(tracks, centrality, vz);
        }

        // Apply systematics
        if (syst != 0) {
            applySystematicUncertainty(syst, tracks, centrality);
        }

        // Quick test: break after a small number of events
        if (do_quicktest && i >= 1000) {
            break;
        }
    }

    // Close files after use
    fileeff->Close();
}
