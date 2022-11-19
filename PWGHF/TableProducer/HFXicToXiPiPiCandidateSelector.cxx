// Copyright 2019-2020 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

/// \file HFXicToXiPiPiCandidateSelector.cxx
/// \brief Ξc± → Ξ∓ π± π± selection task
///
/// \author Jinjoo Seo <jseo@cern.ch>, Inha University

#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "PWGHF/DataModel/HFSecondaryVertex.h"
#include "PWGHF/DataModel/HFCandidateSelectionTables.h"
#include "Common/Core/TrackSelectorPID.h"

#include <TPDGCode.h>
#include <TDatabasePDG.h>

using namespace o2;
using namespace o2::framework;
using namespace o2::aod::hf_cand_casc_prong3;
using namespace o2::analysis::hf_cuts_xic_toxipipi;

/// Struct for applying Ξc selection cuts
struct HFXicToXiPiPiCandidateSelector {
  Produces<aod::HFSelXicCandidate> hfSelXicCandidate;

  Configurable<double> d_pTCandMin{"d_pTCandMin", 0., "Lower bound of candidate pT"};
  Configurable<double> d_pTCandMax{"d_pTCandMax", 36., "Upper bound of candidate pT"};
  Configurable<bool> d_FilterPID{"d_FilterPID", true, "Bool to use or not the PID based on nSigma cut at filtering level"};
  // TPC
  Configurable<bool> b_requireTPC{"b_requireTPC", true, "Flag to require a positive Number of found clusters in TPC"};
  Configurable<double> d_pidTPCMinpT{"d_pidTPCMinpT", 0.15, "Lower bound of track pT for TPC PID"};
  Configurable<double> d_pidTPCMaxpT{"d_pidTPCMaxpT", 20, "Upper bound of track pT for TPC PID"};
  Configurable<double> d_nSigmaTPC{"d_nSigmaTPC", 3., "Nsigma cut on TPC only"};
  Configurable<double> d_nSigmaTPCCombined{"d_nSigmaTPCCombined", 5., "Nsigma cut on TPC combined with TOF"};
  // Configurable<double> d_TPCNClsFindablePIDCut{"d_TPCNClsFindablePIDCut", 70., "Lower bound of TPC findable clusters for good PID"};
  //  TOF
  Configurable<double> d_pidTOFMinpT{"d_pidTOFMinpT", 0.15, "Lower bound of track pT for TOF PID"};
  Configurable<double> d_pidTOFMaxpT{"d_pidTOFMaxpT", 20, "Upper bound of track pT for TOF PID"};
  Configurable<double> d_nSigmaTOF{"d_nSigmaTOF", 3., "Nsigma cut on TOF only"};
  Configurable<double> d_nSigmaTOFCombined{"d_nSigmaTOFCombined", 5., "Nsigma cut on TOF combined with TPC"};

  // topological cuts
  Configurable<std::vector<double>> pTBins{"pTBins", std::vector<double>{hf_cuts_xic_toxipipi::pTBins_v}, "pT bin limits"};
  Configurable<LabeledArray<double>> cuts{"Xic_to_Xi_pi_pi_cuts", {hf_cuts_xic_toxipipi::cuts[0], npTBins, nCutVars, pTBinLabels, cutVarLabels}, "Xic candidate selection per pT bin"};

  /// Conjugate-independent topological cuts
  /// \param candidate is candidate
  /// \return true if candidate passes all cuts
  template <typename T>
  bool selectionTopol(const T& candidate)
  {
    auto candpT = candidate.pt();

    int pTBin = findBin(pTBins, candpT);
    if (pTBin == -1) {
      return false;
    }

    // Topology selection for V0
    //if (candidate.v0cosPA() < cuts->get(pTBin, "V0 CosPA")) {
    //  return false;
    //}

    if (candidate.dcaV0daughters() < cuts->get(pTBin, "DCA V0 Daughters")) {
      return false;
    }

    //if (candidate.dcav0topv() < cuts->get(pTBin, "DCA V0 To PV")) {
    //  return false;
    //}

    if (candidate.dcanegtopv() < cuts->get(pTBin, "DCA Neg To PV")) {
      return false;
    }

    if (candidate.dcapostopv() < cuts->get(pTBin, "DCA Pos To PV")) {
      return false;
    }

    if (candidate.v0radius() < cuts->get(pTBin, "v0radius")) {
      return false;
    }

    if (TMath::Abs(candidate.mLambda() - RecoDecay::getMassPDG(kLambda0)) > cuts->get(pTBin, "v0masswindow")) {
      return false;
    }

    //if (candidate.casccosPA() < cuts->get(pTBin, "Casc CosPA")) {
    //  return false;
    //}

    if (candidate.dcacascdaughters() < cuts->get(pTBin, "DCA Casc Daughters")) {
      return false;
    }

    if (candidate.dcabachtopv() < cuts->get(pTBin, "DCA Bach To PV")) {
      return false;
    }

    if (candidate.cascradius() < cuts->get(pTBin, "cascradius")) {
      return false;
    }

    if (TMath::Abs(candidate.mXi() - RecoDecay::getMassPDG(kXiMinus)) > cuts->get(pTBin, "cascmasswindow")) {
      return false;
    }

    // check that the candidate pT is within the analysis range
    if (candpT < d_pTCandMin || candpT >= d_pTCandMax) {
      return false;
    }

    // cosine of pointing angle
    if (candidate.cpa() <= cuts->get(pTBin, "cos pointing angle")) {
      return false;
    }

    // candidate chi2PCA
    if (candidate.chi2PCA() > cuts->get(pTBin, "Chi2PCA")) {
      return false;
    }

    if (candidate.decayLength() <= cuts->get(pTBin, "decay length")) {
      return false;
    }
    return true;
  }

  /// Conjugate-dependent topological cuts
  /// \param candidate is candidate
  /// \param cascXi is the Cascade with the Xi hypothesis
  /// \param trackPion1 is the track with the pion hypothesis
  /// \param trackPion2 is the track with the pion hypothesis
  /// \return true if candidate passes all cuts for the given Conjugate
  template <typename T1, typename T2, typename T3>
  bool selectionTopolConjugate(const T1& candidate, const T2& cascXi, const T3& trackPion1, const T3& trackPion2)
  {

    auto candpT = candidate.pt();
    int pTBin = findBin(pTBins, candpT);
    if (pTBin == -1) {
      return false;
    }

    // cut on daughter pT
    if (candidate.ptProng2() < cuts->get(pTBin, "pT Xi") || trackPion1.pt() < cuts->get(pTBin, "pT Pi1") || trackPion2.pt() < cuts->get(pTBin, "pT Pi2")) {
      return false;
    }

    if (trackPion1.globalIndex() == candidate.index0Id()) {
      if (std::abs(InvMassXic(candidate) - RecoDecay::getMassPDG(pdg::Code::kXiCPlus)) > cuts->get(pTBin, "m")) {
        return false;
      }
    } else {
      if (std::abs(InvMassXicbar(candidate) - RecoDecay::getMassPDG(pdg::Code::kXiCPlus)) > cuts->get(pTBin, "m")) {
        return false;
      }
    }

    return true;
  }

  using TrksPID = soa::Join<aod::BigTracksPID, aod::pidBayesPi, aod::pidBayesKa, aod::pidBayesPr, aod::pidBayes>;

  void process(aod::HfCandCascProng3 const& candidates, TrksPID const&)
  {
    TrackSelectorPID selectorPion1(kPiPlus);
    selectorPion1.setRangePtTPC(d_pidTPCMinpT, d_pidTPCMaxpT);
    selectorPion1.setRangeNSigmaTPC(-d_nSigmaTPC, d_nSigmaTPC);
    selectorPion1.setRangeNSigmaTPCCondTOF(-d_nSigmaTPCCombined, d_nSigmaTPCCombined);
    selectorPion1.setRangePtTOF(d_pidTOFMinpT, d_pidTOFMaxpT);
    selectorPion1.setRangeNSigmaTOF(-d_nSigmaTOF, d_nSigmaTOF);
    selectorPion1.setRangeNSigmaTOFCondTPC(-d_nSigmaTOFCombined, d_nSigmaTOFCombined);
    // selectorPion1.setRangePtBayes(d_pidBayesMinpT, d_pidBayesMaxpT);

    TrackSelectorPID selectorPion2(selectorPion1);

    // looping over 3-prong candidates
    for (const auto& candidate : candidates) {

      // final selection flag: 0 - rejected, 1 - accepted
      auto statusXic = 0;
      auto statusXicbar = 0;

      if (!(candidate.hfflag() & 1 << DecayType::XicToXiPiPi)) {
        hfSelXicCandidate(0, 0, statusXic, statusXicbar);
        continue;
      }

      const auto& trackPos1 = candidate.index0_as<TrksPID>();  // positive daughter
      const auto& trackPos2 = candidate.index1_as<TrksPID>();  // positive daughter
      const auto& trackCasc = candidate.cascade_as<TrksPID>(); // Cascade daughter (positive for the antiparticles)

      // conjugate-independent topological selectione
      if (!selectionTopol(candidate)) {
        hfSelXicCandidate(0, 0, statusXic, statusXicbar);
        continue;
      }

      // conjugate-dependent topological selection for Xic

      bool topolXic = selectionTopolConjugate(candidate, trackCasc, trackPos1, trackPos2);
      bool topolXicbar = selectionTopolConjugate(candidate, trackCasc, trackPos2, trackPos1);

      if (!topolXic && !topolXicbar) {
        hfSelXicCandidate(0, 0, statusXic, statusXicbar);
        continue;
      }

      auto pidXic = -1;
      auto pidXicbar = -1;

      if (!d_FilterPID) {
        // PID non applied
        pidXic = 1;
        pidXicbar = 1;
      } else {
        // track-level PID selection
        int pidTrackPos1Pion = selectorPion1.getStatusTrackPIDAll(trackPos1);
        int pidTrackPos2Pion = selectorPion2.getStatusTrackPIDAll(trackPos2);

        if (pidTrackPos1Pion == TrackSelectorPID::Status::PIDAccepted &&
            pidTrackPos2Pion == TrackSelectorPID::Status::PIDAccepted) {
          pidXic = 1; // accept Xic
        } else if (pidTrackPos1Pion == TrackSelectorPID::Status::PIDRejected ||
                   pidTrackPos2Pion == TrackSelectorPID::Status::PIDRejected) {
          pidXic = 0; // exclude Xic
        }
        if (pidTrackPos2Pion == TrackSelectorPID::Status::PIDAccepted &&
            pidTrackPos1Pion == TrackSelectorPID::Status::PIDAccepted) {
          pidXicbar = 1; // accept Xicbar
        } else if (pidTrackPos2Pion == TrackSelectorPID::Status::PIDRejected ||
                   pidTrackPos1Pion == TrackSelectorPID::Status::PIDRejected) {
          pidXicbar = 0; // exclude Xicbar
        }
      }

      if (pidXic == 0 && pidXicbar == 0) {
        hfSelXicCandidate(0, 0, statusXic, statusXicbar);
        continue;
      }

      if ((pidXic == -1 || pidXic == 1) && topolXic) {
        statusXic = 1; // identified as Xic
      }
      if ((pidXicbar == -1 || pidXicbar == 1) && topolXicbar) {
        statusXicbar = 1; // identified as Xicbar
      }

      hfSelXicCandidate(0, 0, statusXic, statusXicbar);
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<HFXicToXiPiPiCandidateSelector>(cfgc, TaskName{"hf-xic-toxipipi-candidate-selector"})};
}
