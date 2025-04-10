// call libraries from ROOT and C++

#include "TROOT.h"
#include "TTree.h"
#include "TNtuple.h"
#include "TChain.h"
#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TF1.h"
#include "TVector3.h"
#include <TRandom1.h>
#include <TRandom2.h>
#include <TRandom3.h>
#include <vector>
#include <TLorentzVector.h>
#include "THnSparse.h"
#include <cstring>
#include <ctime>
#include <iostream>
#include <cmath>
#include <math.h>
#include <fstream>
#include <vector>
#include <map>
#include "TFrame.h"
#include "TBenchmark.h"
#include "TSystem.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TProfile2D.h"
#include "TDatime.h"
#include <stdlib.h>
#include <algorithm>
#include <Math/Vector4D.h>

// Additional libraries for event mixing, centrality handling, tracking correction, and histogram manipulation

#include "TMath.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TLatex.h"
#include "TMinuit.h"  // For fitting
#include "TMultiGraph.h" // For multi-graph plotting
#include "TLine.h" // For drawing lines

// Your project specific includes
#include "functions_definition.h" // To define necessary functions for event mixing, binning, etc.
#include "tracking_correction.h"   // To handle tracking efficiency corrections
#include "read_tree.h"             // To read the tree and extract relevant quantities

using namespace ROOT::Math;
using namespace std;
