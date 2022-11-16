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

/// \file HFCandidateCreatorCascade3Prong.cxx
/// \brief Reconstruction of Xic 3-prong cascade decay candidates
/// \note Inspired from HFCandidateCreator3Prong
///
/// \author Jinjoo Seo <jseo@cern.ch>, Inha University

#include "Framework/AnalysisTask.h"
#include "DetectorsVertexing/DCAFitterN.h"
#include "PWGHF/DataModel/HFSecondaryVertex.h"
#include "PWGHF/Utils/utilsBfieldCCDB.h"
#include "Common/Core/trackUtilities.h"
#include "ReconstructionDataFormats/DCA.h"
#include "ReconstructionDataFormats/V0.h"
#include "PWGHF/DataModel/HFCandidateSelectionTables.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::aod::hf_cand;
using namespace o2::aod::hf_cand_casc_prong3;


void customize(std::vector<o2::framework::ConfigParamSpec>& workflowOptions)
{
  ConfigParamSpec optionDoMC{"doMC", VariantType::Bool, true, {"Perform MC matching."}};
  workflowOptions.push_back(optionDoMC);
}

#include "Framework/runDataProcessing.h"

/// Reconstruction of Xic 3-prong decay candidates
struct HFCandidateCreatorCascade3Prong {
  Produces<aod::HfCandCascProng3Base> rowCandidateBase;

  //Configurable<double> magneticField{"d_bz", 5., "magnetic field"};
  Configurable<bool> b_propdca{"b_propdca", true, "create tracks version propagated to PCA"};
  Configurable<double> d_maxr{"d_maxr", 200., "reject PCA's above this radius"};
  Configurable<double> d_maxdzini{"d_maxdzini", 4., "reject (if>0) PCA candidate if tracks DZ exceeds threshold"};
  Configurable<double> d_minparamchange{"d_minparamchange", 1.e-3, "stop iterations if largest change of any X is smaller than this"};
  Configurable<double> d_minrelchi2change{"d_minrelchi2change", 0.9, "stop iterations is chi2/chi2old > this"};
  Configurable<bool> b_dovalplots{"b_dovalplots", true, "do validation plots"};

  OutputObj<TH1F> hMass3{TH1F("hMass3", "3-prong candidates;inv. mass (#Xi #pi #pi) (GeV/#it{c}^{2});entries", 500, 2.3, 2.7)};
  OutputObj<TH1F> hCovPVXX{TH1F("hCovPVXX", "3-prong candidates;XX element of cov. matrix of prim. vtx position (cm^{2});entries", 100, 0., 1.e-4)};
  OutputObj<TH1F> hCovSVXX{TH1F("hCovSVXX", "3-prong candidates;XX element of cov. matrix of sec. vtx position (cm^{2});entries", 100, 0., 0.2)};
  OutputObj<TH2F> hDcaXYProngs{TH2F("hDcaXYProngs", "DCAxy of 3-prong candidates;#it{p}_{T} (GeV/#it{c};#it{d}_{xy}) (#mum);entries", 100, 0., 20., 200, -500., 500.)};
  OutputObj<TH2F> hDcaZProngs{TH2F("hDcaZProngs", "DCAz of 3-prong candidates;#it{p}_{T} (GeV/#it{c};#it{d}_{z}) (#mum);entries", 100, 0., 20., 200, -500., 500.)};

  /// magnetic field setting from CCDB
  Configurable<bool> isRun2{"isRun2", false, "enable Run 2 or Run 3 GRP objects for magnetic field"};
  Configurable<std::string> ccdbUrl{"ccdburl", "http://alice-ccdb.cern.ch", "url of the ccdb repository"};
  Configurable<std::string> ccdbPathLut{"ccdbPathLut", "GLO/Param/MatLUT", "Path for LUT parametrization"};
  Configurable<std::string> ccdbPathGeo{"ccdbPathGeo", "GLO/Config/GeometryAligned", "Path of the geometry file"};
  Configurable<std::string> ccdbPathGrp{"ccdbPathGrp", "GLO/GRP/GRP", "Path of the grp file (Run 2)"};
  Configurable<std::string> ccdbPathGrpMag{"ccdbPathGrpMag", "GLO/Config/GRPMagField", "CCDB path of the GRPMagField object (Run 3)"};
  Service<o2::ccdb::BasicCCDBManager> ccdb;
  o2::base::MatLayerCylSet* lut;
  o2::base::Propagator::MatCorrType matCorr = o2::base::Propagator::MatCorrType::USEMatCorrLUT;
  int runNumber;
  
  float toMicrometers = 10000.; // from cm to µm

  double massPi = RecoDecay::getMassPDG(kPiPlus);
  double massXi = RecoDecay::getMassPDG(kXiMinus);
  double massXiPiPi{0.};
  double magneticField = 0.;

  void init(InitContext const&)
  {
    ccdb->setURL(ccdbUrl);
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking();
    lut = o2::base::MatLayerCylSet::rectifyPtrFromFile(ccdb->get<o2::base::MatLayerCylSet>(ccdbPathLut));
    if (!o2::base::GeometryManager::isGeometryLoaded()) {
      ccdb->get<TGeoManager>(ccdbPathGeo);
    }
    runNumber = 0;
  }

  void process(aod::Collisions const& collisions,
               aod::HfCascade3Prongs const& rowsTrackIndexProng3,
               aod::V0sLinked const&,
               aod::V0Datas const&,
               aod::CascDataExt const& cascades,
               aod::BigTracks const& tracks)
  {
    // 3-prong vertex fitter
    o2::vertexing::DCAFitterN<3> df;
    //df.setBz(magneticField);
    df.setPropagateToPCA(b_propdca);
    df.setMaxR(d_maxr);
    df.setMaxDZIni(d_maxdzini);
    df.setMinParamChange(d_minparamchange);
    df.setMinRelChi2Change(d_minrelchi2change);
    df.setUseAbsDCA(true);

    // loop over triplets of track indices
    for (const auto& rowTrackIndexProng3 : rowsTrackIndexProng3) {

            //for (auto& casc : cascades) {
      auto casc = rowTrackIndexProng3.cascade_as<aod::CascDataExt>();
      if (!casc.has_v0()) {
        continue;
      }
      if (!casc.v0_as<aod::V0sLinked>().has_v0Data()) {
        continue;
      }

      const auto& v0 = casc.v0_as<aod::V0sLinked>().v0Data();
      auto track0 = rowTrackIndexProng3.index0_as<aod::BigTracksExtended>();
      auto track1 = rowTrackIndexProng3.index1_as<aod::BigTracksExtended>();

      /// Set the magnetic field from ccdb.
      /// The static instance of the propagator was already modified in the HFTrackIndexSkimCreator,
      /// but this is not true when running on Run2 data/MC already converted into AO2Ds.
      auto bc = track0.collision().bc_as<aod::BCsWithTimestamps>();
      if (runNumber != bc.runNumber()) {
        LOG(info) << ">>>>>>>>>>>> Current run number: " << runNumber;
        initCCDB(bc, runNumber, ccdb, isRun2 ? ccdbPathGrp : ccdbPathGrpMag, lut, isRun2);
        magneticField = o2::base::Propagator::Instance()->getNominalBz();
        LOG(info) << ">>>>>>>>>>>> Magnetic field: " << magneticField;
        // df.setBz(magneticField); /// put it outside the 'if'! Otherwise we have a difference wrt bz Configurable (< 1 permille) in Run2 conv. data
        // df.print();
      }
      df.setBz(magneticField);

      auto trackV0DaughPos = v0.posTrack_as<aod::BigTracksExtended>();
      auto trackV0DaughNeg = v0.negTrack_as<aod::BigTracksExtended>();
      auto bachTrackCascade = casc.bachelor_as<aod::BigTracksExtended>();

      auto trackParVar0 = getTrackParCov(track0);
      auto trackParVar1 = getTrackParCov(track1);
      auto collision = track0.collision();

      auto charge = -1; 
      if (bachTrackCascade.signed1Pt() > 0) {
        charge = +1; 
      } 


      o2::vertexing::DCAFitterN<2> fitterV0, fitterCasc;
      fitterV0.setBz(magneticField);
      fitterV0.setPropagateToPCA(b_propdca);
      fitterV0.setMaxR(d_maxr);
      fitterV0.setMaxDZIni(d_maxdzini);
      fitterV0.setMinParamChange(d_minparamchange);
      fitterV0.setMinRelChi2Change(d_minrelchi2change);
      fitterV0.setUseAbsDCA(true);

      fitterCasc.setBz(magneticField);
      fitterCasc.setPropagateToPCA(b_propdca);
      fitterCasc.setMaxR(d_maxr);
      fitterCasc.setMaxDZIni(d_maxdzini);
      fitterCasc.setMinParamChange(d_minparamchange);
      fitterCasc.setMinRelChi2Change(d_minrelchi2change);
      fitterCasc.setUseAbsDCA(true);
            
            std::array<float, 3> posCascade = {0.};
            std::array<float, 3> pvecpos = {0.};
            std::array<float, 3> pvecneg = {0.};
            std::array<float, 3> pvecBach = {0.};
            
            auto trackParCovV0DaughPos = getTrackParCov(trackV0DaughPos);
            trackParCovV0DaughPos.propagateTo(v0.posX(), magneticField); // propagate the track to the X closest to the V0 vertex
            auto trackParCovV0DaughNeg = getTrackParCov(trackV0DaughNeg);
            trackParCovV0DaughNeg.propagateTo(v0.negX(), magneticField); // propagate the track to the X closest to the V0 vertex
            
            // Acquire basic tracks
            auto pTrack = getTrackParCov(trackV0DaughPos);
            auto nTrack = getTrackParCov(trackV0DaughNeg);
            auto bTrack = getTrackParCov(bachTrackCascade);
            
            int nCandV0 = fitterV0.process(pTrack, nTrack);
            if (nCandV0 == 0) {
              continue;
            }
            fitterV0.propagateTracksToVertex();
            const auto& v0vtx = fitterV0.getPCACandidate();
            
            // Covariance matrix calculation
            std::array<float, 21> cov0 = {0};
            std::array<float, 21> cov1 = {0};
            std::array<float, 21> covV0 = {0};

            const int momInd[6] = {9, 13, 14, 18, 19, 20}; // cov matrix elements for momentum component
            fitterV0.getTrack(0).getPxPyPzGlo(pvecpos);
            fitterV0.getTrack(1).getPxPyPzGlo(pvecneg);
            fitterV0.getTrack(0).getCovXYZPxPyPzGlo(cov0);
            fitterV0.getTrack(1).getCovXYZPxPyPzGlo(cov1);
            
            const std::array<float, 3> posV0 = {(float)v0vtx[0], (float)v0vtx[1], (float)v0vtx[2]};
            const std::array<float, 3> pvecV0 = {pvecpos[0] + pvecneg[0], pvecpos[1] + pvecneg[1], pvecpos[2] + pvecneg[2]};
            
            for (int i = 0; i < 6; i++) {
              int j = momInd[i]; 
              covV0[j] = cov0[j] + cov1[j];
            }
            auto covVtxV0 = fitterV0.calcPCACovMatrix();
            covV0[0] = covVtxV0(0, 0);
            covV0[1] = covVtxV0(1, 0);
            covV0[2] = covVtxV0(1, 1);
            covV0[3] = covVtxV0(2, 0);
            covV0[4] = covVtxV0(2, 1);
            covV0[5] = covVtxV0(2, 2);
            
            // we build the neutral track to then build the cascade
            auto trackV0 = o2::track::TrackParCov(posV0, pvecV0, covV0, 0);
            trackV0.setQ2Pt(0); // No bending, please
            
            int nCandCasc = fitterCasc.process(trackV0, bTrack);
            if (nCandCasc == 0) {
              continue;
            }
            fitterCasc.propagateTracksToVertex();
            const auto& cascvtx = fitterCasc.getPCACandidate(); 
            posCascade = {(float)cascvtx[0], (float)cascvtx[1], (float)cascvtx[2]};
            
            fitterCasc.getTrack(1).getPxPyPzGlo(pvecBach);
            std::array<float, 3> pvecCascade = {pvecV0[0] + pvecBach[0], pvecV0[1] + pvecBach[1], pvecV0[2] + pvecBach[2]};
            
            // Covariance matrix calculation for cascade
            std::array<float, 21> cov2 = {0};
            std::array<float, 21> cov3 = {0};
            std::array<float, 21> covCascade = {0};
            
            fitterCasc.getTrack(0).getCovXYZPxPyPzGlo(cov2);
            fitterCasc.getTrack(1).getCovXYZPxPyPzGlo(cov3);
            for (int i = 0; i < 6; i++) {
              int j = momInd[i];
              covCascade[j] = cov2[j] + cov3[j];
            }
            auto covVtxCascade = fitterCasc.calcPCACovMatrix();
            covCascade[0] = covVtxCascade(0, 0);
            covCascade[1] = covVtxCascade(1, 0);
            covCascade[2] = covVtxCascade(1, 1);
            covCascade[3] = covVtxCascade(2, 0);
            covCascade[4] = covVtxCascade(2, 1);
            covCascade[5] = covVtxCascade(2, 2);
            
            auto trackCascade = o2::track::TrackParCov(posCascade, pvecCascade, covCascade, 0);
            trackCascade.setQ2Pt(0); // FIXME not sure if this is needed




      if (df.process(trackParVar0, trackParVar1, trackCascade) == 0) {
        continue;
      }

      const auto& secondaryVertex = df.getPCACandidate();
      auto chi2PCA = df.getChi2AtPCACandidate();
      auto covMatrixPCA = df.calcPCACovMatrixFlat();
      hCovSVXX->Fill(covMatrixPCA[0]); // FIXME: Calculation of errorDecayLength(XY) gives wrong values without this line.
      trackParVar0 = df.getTrack(0);
      trackParVar1 = df.getTrack(1);
      auto trackParVar2 = df.getTrack(2);

      // get track momenta
      array<float, 3> pvec0;
      array<float, 3> pvec1;
      array<float, 3> pvec2;
      trackParVar0.getPxPyPzGlo(pvec0);
      trackParVar1.getPxPyPzGlo(pvec1);
      trackParVar2.getPxPyPzGlo(pvec2);

      // get track impact parameters
      // This modifies track momenta!
      auto primaryVertex = getPrimaryVertex(collision);
      auto covMatrixPV = primaryVertex.getCov();
      // TODO PVrefit
      hCovPVXX->Fill(covMatrixPV[0]);
      o2::dataformats::DCA impactParameter0;
      o2::dataformats::DCA impactParameter1;
      o2::dataformats::DCA impactParameter2;
      trackParVar0.propagateToDCA(primaryVertex, magneticField, &impactParameter0);
      trackParVar1.propagateToDCA(primaryVertex, magneticField, &impactParameter1);
      trackParVar2.propagateToDCA(primaryVertex, magneticField, &impactParameter2);
      hDcaXYProngs->Fill(track0.pt(), impactParameter0.getY() * toMicrometers);
      hDcaXYProngs->Fill(track1.pt(), impactParameter1.getY() * toMicrometers);
      hDcaXYProngs->Fill(casc.pt(), impactParameter2.getY() * toMicrometers);
      hDcaZProngs->Fill(track0.pt(), impactParameter0.getZ() * toMicrometers);
      hDcaZProngs->Fill(track1.pt(), impactParameter1.getZ() * toMicrometers);
      hDcaZProngs->Fill(casc.pt(), impactParameter2.getZ() * toMicrometers);

      // get uncertainty of the decay length
      double phi, theta;
      getPointDirection(array{primaryVertex.getX(), primaryVertex.getY(), primaryVertex.getZ()}, secondaryVertex, phi, theta);
      auto errorDecayLength = std::sqrt(getRotatedCovMatrixXX(covMatrixPV, phi, theta) + getRotatedCovMatrixXX(covMatrixPCA, phi, theta));
      auto errorDecayLengthXY = std::sqrt(getRotatedCovMatrixXX(covMatrixPV, phi, 0.) + getRotatedCovMatrixXX(covMatrixPCA, phi, 0.));

      int hfFlag = 1 << DecayType::XicToXiPiPi;

      // fill candidate table rows
      rowCandidateBase(collision.globalIndex(),
                       primaryVertex.getX(), primaryVertex.getY(), primaryVertex.getZ(),
                       secondaryVertex[0], secondaryVertex[1], secondaryVertex[2],
                       errorDecayLength, errorDecayLengthXY,
                       chi2PCA,
                       pvec0[0], pvec0[1], pvec0[2],
                       pvec1[0], pvec1[1], pvec1[2],
                       pvec2[0], pvec2[1], pvec2[2],
                       impactParameter0.getY(), impactParameter1.getY(), impactParameter2.getY(),
                       std::sqrt(impactParameter0.getSigmaY2()), std::sqrt(impactParameter1.getSigmaY2()), std::sqrt(impactParameter2.getSigmaY2()),
                       rowTrackIndexProng3.cascadeId(), rowTrackIndexProng3.index0Id(), rowTrackIndexProng3.index1Id(),
                       charge,
                       fitterV0.getChi2AtPCACandidate(), fitterCasc.getChi2AtPCACandidate(),
                       trackV0DaughPos.dcaXY(),
                       trackV0DaughNeg.dcaXY(),
                       bachTrackCascade.dcaXY(),
                       hfFlag);

      // fill histograms
      if (b_dovalplots) {
        // calculate invariant mass
        auto arrayMomenta = array{pvec1, pvec0, pvec1};
        massXiPiPi = RecoDecay::m(std::move(arrayMomenta), array{massXi, massPi, massPi});
        hMass3->Fill(massXiPiPi);
      }
      //}
    }
  }
};

/// Extends the base table with expression columns.
struct HFCandidateCreatorCascade3ProngExpressions {
  Produces<aod::HfCandCascProng3MCRec> rowMCMatchRec;
  Produces<aod::HfCandCascProng3MCGen> rowMCMatchGen;

  Spawns<aod::HfCandCascProng3Ext> rowCandidateProng3;
  void init(InitContext const&) {}

  /// Performs MC matching.
  void processMC(aod::HfCascade3Prongs const& candidates,
                 aod::BigTracksMC const& tracks,
                 aod::CascDataExt const&,
                 aod::McParticles const& particlesMC)
  {
    rowCandidateProng3->bindExternalIndices(&tracks);
    int indexRec = -1;
    int8_t sign = 0;
    int8_t flag = 0;
    int8_t origin = 0;

    // Match reconstructed candidates.n
    // Spawned table can be used directly
    for (auto& candidate : candidates) {
      // Printf("New rec. candidate");
      flag = 0;
      origin = 0;
      //auto XiTrack = candidate.cascade_as();
      auto arrayDaughters = array{candidate.cascade_as<aod::BigTracksMC>(), candidate.index0_as<aod::BigTracksMC>(), candidate.index1_as<aod::BigTracksMC>()};

      // Ξc± → Ξ∓ π± π±
      // Printf("Checking Ξc± → Ξ∓ π± π±");
      indexRec = RecoDecay::getMatchedMCRec(particlesMC, arrayDaughters, pdg::Code::kXiCPlus, array{-kXiMinus, +kPiPlus, +kPiPlus}, true, &sign);
      if (indexRec > -1) {
        flag = sign * (1 << DecayType::XicToXiPiPi);
      }

      // Check whether the particle is non-prompt (from a b quark).
      if (flag != 0) {
        auto particle = particlesMC.rawIteratorAt(indexRec);
        origin = RecoDecay::getCharmHadronOrigin(particlesMC, particle);
      }

      rowMCMatchRec(flag, origin);
    }

    // Match generated particles.
    for (auto& particle : particlesMC) {
      // Printf("New gen. candidate");
      flag = 0;
      origin = 0;

      // Ξc± → p± K∓ π±
      // Printf("Checking Ξc± → p± K∓ π±");
      if (RecoDecay::isMatchedMCGen(particlesMC, particle, pdg::Code::kXiCPlus, array{-kXiMinus, +kPiPlus, +kPiPlus}, true, &sign)) {
        flag = sign * (1 << DecayType::XicToXiPiPi);
      }

      // Check whether the particle is non-prompt (from a b quark).
      if (flag != 0) {
        origin = RecoDecay::getCharmHadronOrigin(particlesMC, particle);
      }

      rowMCMatchGen(flag, origin);
    }
  }

  PROCESS_SWITCH(HFCandidateCreatorCascade3ProngExpressions, processMC, "Process MC", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<HFCandidateCreatorCascade3Prong>(cfgc, TaskName{"hf-cand-creator-cascade-3prong"})
    adaptAnalysisTask<HFCandidateCreatorCascade3ProngExpressions>(cfgc, TaskName{"hf-cand-creator-cascade-3prong-expressions"})};
}
