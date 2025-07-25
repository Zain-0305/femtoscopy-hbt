#ifndef TRACK_CORRECTIONS_H
#define TRACK_CORRECTIONS_H

#include "call_libraries.h"
#include <stdexcept>

/**
 * @file track_corrections.h
 * @brief Track efficiency and fake rate corrections for HBT analysis
 * 
 * Combines functionality from tracking_correction.h and trk_efficiency_correction.h
 * with added safety checks and HBT-specific optimizations
 */

namespace HBT {
    namespace Corrections {

        // Constants ===========================================================
        constexpr double MAX_ETA = 2.4;       // Tracking acceptance
        constexpr double MIN_PT = 0.15;       // GeV/c (HBT typical cutoff)
        constexpr double MAX_PT = 500.0;      // GeV/c

        /** 
         * @brief Validate track kinematics before correction
         * @throws std::range_error if values are outside physical limits
         */
        inline void ValidateTrack(double pt, double eta) {
            if (std::abs(eta) > MAX_ETA) {
                throw std::range_error("Track eta " + std::to_string(eta) + 
                                      " exceeds acceptance");
            }
            if (pt < MIN_PT || pt > MAX_PT) {
                throw std::range_error("Track pT " + std::to_string(pt) + 
                                      " outside valid range");
            }
        }

        /**
         * @brief Get bin content with safety checks
         * @return Bin content or fallback value if bin is invalid
         */
        inline double SafeGetBinContent(const TH2* hist, double eta, double pt, 
                                       double fallback = 1.0) {
            int xbin = hist->GetXaxis()->FindBin(eta);
            int ybin = hist->GetYaxis()->FindBin(pt);
            
            // Check if bins are in valid range
            if (xbin < 1 || xbin > hist->GetNbinsX() || 
                ybin < 1 || ybin > hist->GetNbinsY()) {
                return fallback;
            }
            
            double val = hist->GetBinContent(xbin, ybin);
            return (val < 0.0001 || val > 0.9999) ? fallback : val;
        }

        /**
         * @brief Basic efficiency correction (simplified interface)
         * @param eff_map 2D histogram of efficiency vs (η,pT)
         */
        inline double EfficiencyCorrection(const TH2* eff_map, 
                                         double pt, double eta) {
            ValidateTrack(pt, eta);
            double eff = SafeGetBinContent(eff_map, eta, pt);
            return 1.0 / eff;
        }

        /**
         * @brief Complete track correction including fakes/seconadries/multiples
         * @param eff_map Efficiency map
         * @param fake_map Fake rate map
         * @param sec_map Secondary fraction map  
         * @param mul_map Multiple reconstruction map
         * @param mode Correction combination mode:
         *             0 = full correction (default)
         *             1 = efficiency+fakes only
         *             2 = efficiency only
         */
        inline double FullTrackCorrection(
            const TH2* eff_map, const TH2* fake_map,
            const TH2* sec_map = nullptr, const TH2* mul_map = nullptr,
            double pt = 0, double eta = 0, int mode = 0) 
        {
            ValidateTrack(pt, eta);
            
            // Get all components
            double eff = SafeGetBinContent(eff_map, eta, pt);
            double fake = SafeGetBinContent(fake_map, eta, pt, 0.0);
            
            // Early return for simplified modes
            if (mode == 2) return 1.0 / eff;          // Efficiency only
            if (mode == 1) return (1.0 - fake) / eff; // Efficiency + fakes
            
            // Full correction
            double sec = sec_map ? SafeGetBinContent(sec_map, eta, pt, 0.0) : 0.0;
            double mul = mul_map ? SafeGetBinContent(mul_map, eta, pt, 0.0) : 0.0;
            
            return (1.0 - fake) * (1.0 - sec) / eff / (1.0 + mul);
        }

    } // namespace Corrections
} // namespace HBT

#endif // TRACK_CORRECTIONS_H
