# femtoscopy-hbt
 HBT correlation code
This code is going to be used to measure the HBT correlations in 1 and 3D for charged hadrons in XeXe collisions using 2017 dataset.

Instructions
The code was updated to work on LXPLUS9 in order to use more resources from the cluster in HTCondor. To login in LXPLUS9 use:

ssh username@lxplus9.cern.ch
Once logged into LXPLUS machines setup the CMSSW version as follows:

cmsrel CMSSW_13_0_5
cd CMSSW_13_0_5/src/
cmsenv
and work inside of the src folder. The code was created to generate histograms as function of centrality or multiplicity.

For centrality depentency we have the following bins: 0-10%, 10-30%, 30-50%, 50-70% and 70-100% (not a good idea to use given the higher EM contamination). This code takes a long time to run due the correlations between large number of tracks and also because of mixing. To get the centrality dependency, use
