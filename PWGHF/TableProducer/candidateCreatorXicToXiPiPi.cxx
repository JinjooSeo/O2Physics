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

/// \file candidateCreatorXicToXiPiPi.cxx
/// \brief Reconstruction of heavy-flavour cascade 3-prong decay candidates
/// \note Extended from candidateCreator2Prong.cxx candidateCreator3Prong.cxx candidateCreatorToXiPi.cxx
///
/// \author Jinjoo Seo <jseo@cern.ch>, Heidelberg University

#include <KFParticleBase.h>
#include <KFParticle.h>
#include <KFPTrack.h>
#include <KFPVertex.h>
#include <KFVertex.h>

#include <TPDGCode.h>

#include "CommonConstants/PhysicsConstants.h"
#include "DCAFitter/DCAFitterN.h"
#include "Framework/AnalysisTask.h"
#include "Framework/runDataProcessing.h"
#include "ReconstructionDataFormats/DCA.h"

#include "Common/Core/trackUtilities.h"
#include "Tools/KFparticle/KFUtilities.h"

#include "PWGHF/DataModel/CandidateReconstructionTables.h"
#include "PWGHF/Utils/utilsBfieldCCDB.h"

using namespace o2;
using namespace o2::analysis;
using namespace o2::aod::hf_cand_casc_lf;
using namespace o2::framework;

/// Reconstruction of heavy-flavour 3-prong decay candidates
struct HfCandidateCreatorCasc3Prong {
  Produces<aod::HfCandCascProng3Base> rowCandidateBase;

  // vertexing
  Configurable<bool> propagateToPCA{"propagateToPCA", true, "create tracks version propagated to PCA"};
  Configurable<bool> useAbsDCA{"useAbsDCA", false, "Minimise abs. distance rather than chi2"};
  Configurable<bool> useWeightedFinalPCA{"useWeightedFinalPCA", false, "Recalculate vertex position using track covariances, effective only if useAbsDCA is true"};
  Configurable<double> maxR{"maxR", 200., "reject PCA's above this radius"};
  Configurable<double> maxDZIni{"maxDZIni", 4., "reject (if>0) PCA candidate if tracks DZ exceeds threshold"};
  Configurable<double> minParamChange{"minParamChange", 1.e-3, "stop iterations if largest change of any X is smaller than this"};
  Configurable<double> minRelChi2Change{"minRelChi2Change", 0.9, "stop iterations is chi2/chi2old > this"};
  Configurable<double> maxChi2{"maxChi2", 100., "discard vertices with chi2/Nprongs > this (or sum{DCAi^2}/Nprongs for abs. distance minimization)"};
  Configurable<bool> rejDiffCollTrack{"rejDiffCollTrack", true, "Reject tracks coming from different collisions"};
  Configurable<bool> fillHistograms{"fillHistograms", true, "do validation plots"};
  
  // magnetic field setting from CCDB
  Configurable<bool> isRun2{"isRun2", false, "enable Run 2 or Run 3 GRP objects for magnetic field"};
  Configurable<std::string> ccdbUrl{"ccdbUrl", "http://alice-ccdb.cern.ch", "url of the ccdb repository"};
  Configurable<std::string> ccdbPathLut{"ccdbPathLut", "GLO/Param/MatLUT", "Path for LUT parametrization"};
  Configurable<std::string> ccdbPathGrp{"ccdbPathGrp", "GLO/GRP/GRP", "Path of the grp file (Run 2)"};
  Configurable<std::string> ccdbPathGrpMag{"ccdbPathGrpMag", "GLO/Config/GRPMagField", "CCDB path of the GRPMagField object (Run 3)"};
  
  // cascade cuts
  Configurable<bool> doCascadePreselection{"doCascadePreselection", true, "Use invariant mass and dcaXY cuts to preselect cascade candidates"};
  Configurable<double> massToleranceCascade{"massToleranceCascade", 0.01, "Invariant mass tolerance for cascade"};
  Configurable<float> dcaXYToPVCascadeMax{"dcaXYToPVCascadeMax", 3, "Max cascade DCA to PV in xy plane"};

  Service<o2::ccdb::BasicCCDBManager> ccdb;
  o2::base::MatLayerCylSet* lut;
  o2::base::Propagator::MatCorrType matCorr = o2::base::Propagator::MatCorrType::USEMatCorrLUT;

  int runNumber{0};
  float toMicrometers = 10000.; // from cm to µm
  double massPi{0.};
  double massXi{0.};
  double massXiPiPi{0.};
  double bz{0.};

  OutputObj<TH1F> hMass3{TH1F("hMass3", "3-prong candidates;inv. mass (#Xi #pi #pi) (GeV/#it{c}^{2});entries", 500, 2.3, 2.7)};
  OutputObj<TH1F> hCovPVXX{TH1F("hCovPVXX", "3-prong candidates;XX element of cov. matrix of prim. vtx. position (cm^{2});entries", 100, 0., 1.e-4)};
  OutputObj<TH1F> hCovSVXX{TH1F("hCovSVXX", "3-prong candidates;XX element of cov. matrix of sec. vtx. position (cm^{2});entries", 100, 0., 0.2)};
  OutputObj<TH1F> hCovPVYY{TH1F("hCovPVYY", "3-prong candidates;YY element of cov. matrix of prim. vtx. position (cm^{2});entries", 100, 0., 1.e-4)};
  OutputObj<TH1F> hCovSVYY{TH1F("hCovSVYY", "3-prong candidates;YY element of cov. matrix of sec. vtx. position (cm^{2});entries", 100, 0., 0.2)};
  OutputObj<TH1F> hCovPVXZ{TH1F("hCovPVXZ", "3-prong candidates;XZ element of cov. matrix of prim. vtx. position (cm^{2});entries", 100, -1.e-4, 1.e-4)};
  OutputObj<TH1F> hCovSVXZ{TH1F("hCovSVXZ", "3-prong candidates;XZ element of cov. matrix of sec. vtx. position (cm^{2});entries", 100, -1.e-4, 0.2)};
  OutputObj<TH1F> hCovPVZZ{TH1F("hCovPVZZ", "3-prong candidates;ZZ element of cov. matrix of prim. vtx. position (cm^{2});entries", 100, 0., 1.e-4)};
  OutputObj<TH1F> hCovSVZZ{TH1F("hCovSVZZ", "3-prong candidates;ZZ element of cov. matrix of sec. vtx. position (cm^{2});entries", 100, 0., 0.2)};
  OutputObj<TH2F> hDcaXYProngs{TH2F("hDcaXYProngs", "DCAxy of 3-prong candidates;#it{p}_{T} (GeV/#it{c};#it{d}_{xy}) (#mum);entries", 100, 0., 20., 200, -500., 500.)};
  OutputObj<TH2F> hDcaZProngs{TH2F("hDcaZProngs", "DCAz of 3-prong candidates;#it{p}_{T} (GeV/#it{c};#it{d}_{z}) (#mum);entries", 100, 0., 20., 200, -500., 500.)};
  OutputObj<TH1F> hVertexerType{TH1F("hVertexerType", "Use KF or DCAFitterN;Vertexer type;entries", 2, -0.5, 1.5)}; // See o2::aod::hf_cand::VertexerType

  void init(InitContext const&)
  {
    std::array<bool, 2> doprocessDF{doprocessPvRefitWithDCAFitterN, doprocessNoPvRefitWithDCAFitterN};
    std::array<bool, 2> doprocessKF{doprocessPvRefitWithKFParticle, doprocessNoPvRefitWithKFParticle};
    if ((std::accumulate(doprocessDF.begin(), doprocessDF.end(), 0) + std::accumulate(doprocessKF.begin(), doprocessKF.end(), 0)) != 1) {
      LOGP(fatal, "Only one process function can be enabled at a time.");
    }
    if (std::accumulate(doprocessDF.begin(), doprocessDF.end(), 0) == 1) {
      hVertexerType->Fill(aod::hf_cand::VertexerType::DCAFitter);
    }
    if (std::accumulate(doprocessKF.begin(), doprocessKF.end(), 0) == 1) {
      hVertexerType->Fill(aod::hf_cand::VertexerType::KfParticle);
    }

    massPi = MassPiPlus;
    massXi = MassXiMinus;
    ccdb->setURL(ccdbUrl);
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking();
    lut = o2::base::MatLayerCylSet::rectifyPtrFromFile(ccdb->get<o2::base::MatLayerCylSet>(ccdbPathLut));
    runNumber = 0;
  }

  // Will be modified
  using CascFull = soa::Join<aod::CascDatas, aod::CascCovs>;
  using V0Full = soa::Join<aod::V0Datas, aod::V0Covs>;

  template <bool doPvRefit, typename CandType, typename TTracks>
  void runCreatorCasc3ProngWithDCAFitterN(aod::Collisions const& collisions,
                                          CandType const& rowsTrackIndexCasc3Prong,
                                          aod::V0sLinked const& v0Linked,
                                          V0Full const& v0s,
                                          CascFull const& cascs,
                                          TTracks const& tracks,
                                          aod::BCsWithTimestamps const& bcWithTimeStamps)
  {
    // 3-prong vertex fitter
    o2::vertexing::DCAFitterN<3> df;
    // df.setBz(bz);
    df.setPropagateToPCA(propagateToPCA);
    df.setMaxR(maxR);
    df.setMaxDZIni(maxDZIni);
    df.setMinParamChange(minParamChange);
    df.setMinRelChi2Change(minRelChi2Change);
    df.setUseAbsDCA(useAbsDCA);
    df.setWeightedFinalPCA(useWeightedFinalPCA);

    // loop over triplets of track indices
    for (const auto& rowTrackIndexCasc3Prong : rowsTrackIndexCasc3Prong) {
      auto track0 = rowTrackIndexCasc3Prong prong0_as<TTracks>();
      auto track1 = rowTrackIndexCasc3Prong prong1_as<TTracks>();
      auto casc = rowTrackIndexCasc3Prong.cascade_as<CascFull>();

      // preselect cascade candidates
      if (doCascadePreselection) {
        if (std::abs(casc.dcaXYCascToPV()) > dcaXYToPVCascadeMax) {
          continue;
        }
        if (std::abs(casc.mXi() - massXiFromPDG) > massToleranceCascade) {
          continue;
        }
      }

      //----------------accessing particles in the decay chain-------------
      // cascade daughter - charged particle
      auto trackXiDauCharged = casc.bachelor_as<TTracks>(); // pion <- xi track from TTracks table

      if (!casc.v0_as<aod::V0sLinked>().has_v0Data()) { // check that V0 data are stored
          continue;
        }
      auto v0 = casc.v0_as<aod::V0sLinked>();
      auto v0Element = v0.v0Data_as<V0Full>(); // V0 element from LF table containing V0 info
      // V0 positive daughter
      auto trackV0Dau0 = v0Element.posTrack_as<TTracks>(); // p <- V0 track (positive track) from TTracks table
      // V0 negative daughter
      auto trackV0Dau1 = v0Element.negTrack_as<TTracks>(); // pion <- V0 track (negative track) from TTracks table

      // check that particles come from the same collision
      if (rejDiffCollTrack) {
        if (trackV0Dau0.collisionId() != trackV0Dau1.collisionId()) {
          continue;
        }
        if (trackXiDauCharged.collisionId() != trackV0Dau0.collisionId()) {
          continue;
        }
      }

      //--------------------------reconstruct V0 track---------------------------
      // pseudorapidity
      double pseudorapV0PosDau = trackV0Dau0.eta();
      double pseudorapV0NegDau = trackV0Dau1.eta();

      // pion & p <- V0 tracks
      auto trackParCovV0Dau0 = getTrackParCov(trackV0Dau0);
      auto trackParCovV0Dau1 = getTrackParCov(trackV0Dau1);

      // info from LF table
      std::array<float, 3> pVecV0 = {casc.pxlambda(), casc.pylambda(), casc.pzlambda()}; // pVec stands for vector containing the 3-momentum components
      std::array<float, 3> vertexV0 = {casc.xlambda(), casc.ylambda(), casc.zlambda()};
      std::array<float, 3> pVecV0Dau0 = {casc.pxpos(), casc.pypos(), casc.pzpos()};
      std::array<float, 3> pVecV0Dau1 = {casc.pxneg(), casc.pyneg(), casc.pzneg()};

      //-----------------------------reconstruct cascade track-----------------------------
      // pseudorapidity
      double pseudorapPiFromCas = trackXiDauCharged.eta();

      // pion <- casc track to be processed with DCAfitter
      auto trackParCovXiDauCharged = getTrackParCov(trackXiDauCharged);

      // info from LF table
      std::array<float, 3> vertexCasc = {casc.x(), casc.y(), casc.z()};
      std::array<float, 3> pVecCasc = {casc.px(), casc.py(), casc.pz()};
      std::array<float, 21> covCasc = {0.};
      constexpr int MomInd[6] = {9, 13, 14, 18, 19, 20}; // cov matrix elements for momentum component
      for (int i = 0; i < 6; i++) {
        covCasc[MomInd[i]] = casc.momentumCovMat()[i];
        covCasc[i] = casc.positionCovMat()[i];
      }

      // create cascade track
      o2::track::TrackParCov trackCasc;
      if (trackXiDauCharged.sign() > 0) {
        trackCasc = o2::track::TrackParCov(vertexCasc, pVecCasc, covCasc, 1, true);
      } else if (trackXiDauCharged.sign() < 0) {
        trackCasc = o2::track::TrackParCov(vertexCasc, pVecCasc, covCasc, -1, true);
      } else {
        continue;
      }
      
      trackCasc.setAbsCharge(1);
      trackCasc.setPID(o2::track::PID::XiMinus);

      std::array<float, 3> pVecPionFromCasc = {casc.pxbach(), casc.pybach(), casc.pzbach()};

      auto trackParVar0 = getTrackParCov(track0);
      auto trackParVar1 = getTrackParCov(track1);
      auto collision = rowTrackIndexCasc3Prong.collision();

      /// Set the magnetic field from ccdb.
      /// The static instance of the propagator was already modified in the HFTrackIndexSkimCreator,
      /// but this is not true when running on Run2 data/MC already converted into AO2Ds.
      auto bc = collision.template bc_as<aod::BCsWithTimestamps>();
      if (runNumber != bc.runNumber()) {
        LOG(info) << ">>>>>>>>>>>> Current run number: " << runNumber;
        initCCDB(bc, runNumber, ccdb, isRun2 ? ccdbPathGrp : ccdbPathGrpMag, lut, isRun2);
        bz = o2::base::Propagator::Instance()->getNominalBz();
        LOG(info) << ">>>>>>>>>>>> Magnetic field: " << bz;
        // df.setBz(bz); /// put it outside the 'if'! Otherwise we have a difference wrt bz Configurable (< 1 permille) in Run2 conv. data
        // df.print();
      }
      df.setBz(bz);

      // reconstruct the 3-prong secondary vertex
      if (df.process(trackParVar0, trackParVar1, trackCasc) == 0) {
        continue;
      }
  
      const auto& secondaryVertex = df.getPCACandidate();
      auto chi2PCA = df.getChi2AtPCACandidate();
      auto covMatrixPCA = df.calcPCACovMatrixFlat();
      hCovSVXX->Fill(covMatrixPCA[0]); // FIXME: Calculation of errorDecayLength(XY) gives wrong values without this line.
      hCovSVYY->Fill(covMatrixPCA[2]);
      hCovSVXZ->Fill(covMatrixPCA[3]);
      hCovSVZZ->Fill(covMatrixPCA[5]);
      trackParVar0 = df.getTrack(0);
      trackParVar1 = df.getTrack(1);
      trackCasc = df.getTrack(2);

      // get track momenta
      std::array<float, 3> pvec0;
      std::array<float, 3> pvec1;
      std::array<float, 3> pvec2;
      trackParVar0.getPxPyPzGlo(pvec0);
      trackParVar1.getPxPyPzGlo(pvec1);
      trackCasc.getPxPyPzGlo(pvec2);

      // get track impact parameters
      // This modifies track momenta!
      auto primaryVertex = getPrimaryVertex(collision);
      auto covMatrixPV = primaryVertex.getCov();
      if constexpr (doPvRefit) {
        /// use PV refit
        /// Using it in the rowCandidateBase all dynamic columns shall take it into account
        // coordinates
        primaryVertex.setX(rowTrackIndexCasc3Prong.pvRefitX());
        primaryVertex.setY(rowTrackIndexCasc3Prong.pvRefitY());
        primaryVertex.setZ(rowTrackIndexCasc3Prong.pvRefitZ());
        // covariance matrix
        primaryVertex.setSigmaX2(rowTrackIndexCasc3Prong.pvRefitSigmaX2());
        primaryVertex.setSigmaXY(rowTrackIndexCasc3Prong.pvRefitSigmaXY());
        primaryVertex.setSigmaY2(rowTrackIndexCasc3Prong.pvRefitSigmaY2());
        primaryVertex.setSigmaXZ(rowTrackIndexCasc3Prong.pvRefitSigmaXZ());
        primaryVertex.setSigmaYZ(rowTrackIndexCasc3Prong.pvRefitSigmaYZ());
        primaryVertex.setSigmaZ2(rowTrackIndexCasc3Prong.pvRefitSigmaZ2());
        covMatrixPV = primaryVertex.getCov();
      }
      hCovPVXX->Fill(covMatrixPV[0]);
      hCovPVYY->Fill(covMatrixPV[2]);
      hCovPVXZ->Fill(covMatrixPV[3]);
      hCovPVZZ->Fill(covMatrixPV[5]);
      o2::dataformats::DCA impactParameter0;
      o2::dataformats::DCA impactParameter1;
      o2::dataformats::DCA impactParameterCasc;
      trackParVar0.propagateToDCA(primaryVertex, bz, &impactParameter0);
      trackParVar1.propagateToDCA(primaryVertex, bz, &impactParameter1);
      trackCasc.propagateToDCA(primaryVertex, bz, &impactParameter2);
      hDcaXYProngs->Fill(track0.pt(), impactParameter0.getY() * toMicrometers);
      hDcaXYProngs->Fill(track1.pt(), impactParameter1.getY() * toMicrometers);
      hDcaXYProngs->Fill(trackCasc.pt(), impactParameterCasc.getY() * toMicrometers);
      hDcaZProngs->Fill(track0.pt(), impactParameter0.getZ() * toMicrometers);
      hDcaZProngs->Fill(track1.pt(), impactParameter1.getZ() * toMicrometers);
      hDcaZProngs->Fill(trackCasc.pt(), impactParameterCasc.getZ() * toMicrometers);

      // get uncertainty of the decay length
      double phi, theta;
      getPointDirection(std::array{primaryVertex.getX(), primaryVertex.getY(), primaryVertex.getZ()}, secondaryVertex, phi, theta);
      auto errorDecayLength = std::sqrt(getRotatedCovMatrixXX(covMatrixPV, phi, theta) + getRotatedCovMatrixXX(covMatrixPCA, phi, theta));
      auto errorDecayLengthXY = std::sqrt(getRotatedCovMatrixXX(covMatrixPV, phi, 0.) + getRotatedCovMatrixXX(covMatrixPCA, phi, 0.));

      // fill candidate table rows !TODO
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
                       rowTrackIndexProng3.index0Id(), rowTrackIndexProng3.index1Id(), rowTrackIndexProng3.cascadeId(),
                       charge, posCascade[0], posCascade[1], posCascade[2], posV0[0], posV0[1], posV0[2],
                       pvecpos[0], pvecpos[1], pvecpos[2],
                       pvecneg[0], pvecneg[1], pvecneg[2],
                       pvecBach[0], pvecBach[1], pvecBach[2],
                       fitterV0.getChi2AtPCACandidate(), fitterCasc.getChi2AtPCACandidate(),
                       trackV0DaughPos.dcaXY(),
                       trackV0DaughNeg.dcaXY(),
                       bachTrackCascade.dcaXY(),
                       rowTrackIndexCasc3Prong.hfflag());

      // fill histograms
      if (fillHistograms) {
        // calculate invariant mass
        auto arrayMomenta = std::array{pvec0, pvec1, pvec2};
        massXiPiPi = RecoDecay::m(std::move(arrayMomenta), std::array{massXi, massPi, massPi});
        hMass3->Fill(massXiPiPi);
      }
    }
  }

  template <bool doPvRefit, typename CandType, typename TTracks>
  void runCreatorCasc3ProngWithKFParticle(aod::Collisions const& collisions,
                                          CandType const& rowsTrackIndexCasc3Prong,
                                          aod::V0sLinked const&,
                                          V0Full const&,
                                          CascFull const&,
                                          TTracks const& tracks,
                                          aod::BCsWithTimestamps const& bcWithTimeStamps)
  {

     // loop over triplets of track indices
    for (const auto& rowTrackIndexCasc3Prong : rowsTrackIndexCasc3Prong) {
      auto track0 = rowTrackIndexCasc3Prong prong0_as<TTracks>();
      auto track1 = rowTrackIndexCasc3Prong prong1_as<TTracks>();
      auto casc = rowTrackIndexCasc3Prong.cascade_as<CascFull>();

      // preselect cascade candidates
      if (doCascadePreselection) {
        if (std::abs(casc.dcaXYCascToPV()) > dcaXYToPVCascadeMax) {
          continue;
        }
        if (std::abs(casc.mXi() - massXiFromPDG) > massToleranceCascade) {
          continue;
        }
      }

      //----------------accessing particles in the decay chain-------------
      // cascade daughter - charged particle
      auto trackXiDauCharged = casc.bachelor_as<TTracks>(); // pion <- xi track from TTracks table

      if (!casc.v0_as<aod::V0sLinked>().has_v0Data()) { // check that V0 data are stored
          continue;
        }
      auto v0 = casc.v0_as<aod::V0sLinked>();
      auto v0Element = v0.v0Data_as<V0Full>(); // V0 element from LF table containing V0 info
      // V0 positive daughter
      auto trackV0Dau0 = v0Element.posTrack_as<TTracks>(); // p <- V0 track (positive track) from TTracks table
      // V0 negative daughter
      auto trackV0Dau1 = v0Element.negTrack_as<TTracks>(); // pion <- V0 track (negative track) from TTracks table

      // check that particles come from the same collision
      if (rejDiffCollTrack) {
        if (trackV0Dau0.collisionId() != trackV0Dau1.collisionId()) {
          continue;
        }
        if (trackXiDauCharged.collisionId() != trackV0Dau0.collisionId()) {
          continue;
        }
      }

      //--------------------------reconstruct V0 track---------------------------
      // pseudorapidity
      double pseudorapV0PosDau = trackV0Dau0.eta();
      double pseudorapV0NegDau = trackV0Dau1.eta();

      // pion & p <- V0 tracks
      auto trackParCovV0Dau0 = getTrackParCov(trackV0Dau0);
      auto trackParCovV0Dau1 = getTrackParCov(trackV0Dau1);

      // info from LF table
      std::array<float, 3> pVecV0 = {casc.pxlambda(), casc.pylambda(), casc.pzlambda()}; // pVec stands for vector containing the 3-momentum components
      std::array<float, 3> vertexV0 = {casc.xlambda(), casc.ylambda(), casc.zlambda()};
      std::array<float, 3> pVecV0Dau0 = {casc.pxpos(), casc.pypos(), casc.pzpos()};
      std::array<float, 3> pVecV0Dau1 = {casc.pxneg(), casc.pyneg(), casc.pzneg()};

      //-----------------------------reconstruct cascade track-----------------------------
      // pseudorapidity
      double pseudorapPiFromCas = trackXiDauCharged.eta();

      // pion <- casc track to be processed with DCAfitter
      auto trackParCovXiDauCharged = getTrackParCov(trackXiDauCharged);

      // info from LF table
      std::array<float, 3> vertexCasc = {casc.x(), casc.y(), casc.z()};
      std::array<float, 3> pVecCasc = {casc.px(), casc.py(), casc.pz()};
      std::array<float, 21> covCasc = {0.};
      constexpr int MomInd[6] = {9, 13, 14, 18, 19, 20}; // cov matrix elements for momentum component
      for (int i = 0; i < 6; i++) {
        covCasc[MomInd[i]] = casc.momentumCovMat()[i];
        covCasc[i] = casc.positionCovMat()[i];
      }

      // create cascade track
      o2::track::TrackParCov trackCasc;
      if (trackXiDauCharged.sign() > 0) {
        trackCasc = o2::track::TrackParCov(vertexCasc, pVecCasc, covCasc, 1, true);
      } else if (trackXiDauCharged.sign() < 0) {
        trackCasc = o2::track::TrackParCov(vertexCasc, pVecCasc, covCasc, -1, true);
      } else {
        continue;
      }
      
      trackCasc.setAbsCharge(1);
      trackCasc.setPID(o2::track::PID::XiMinus);

      std::array<float, 3> pVecPionFromCasc = {casc.pxbach(), casc.pybach(), casc.pzbach()};

      auto trackParVar0 = getTrackParCov(track0);
      auto trackParVar1 = getTrackParCov(track1);
      auto collision = rowTrackIndexCasc3Prong.collision();

      /// Set the magnetic field from ccdb.
      /// The static instance of the propagator was already modified in the HFTrackIndexSkimCreator,
      /// but this is not true when running on Run2 data/MC already converted into AO2Ds.
      auto bc = collision.template bc_as<aod::BCsWithTimestamps>();
      if (runNumber != bc.runNumber()) {
        LOG(info) << ">>>>>>>>>>>> Current run number: " << runNumber;
        initCCDB(bc, runNumber, ccdb, isRun2 ? ccdbPathGrp : ccdbPathGrpMag, lut, isRun2);
        bz = o2::base::Propagator::Instance()->getNominalBz();
        LOG(info) << ">>>>>>>>>>>> Magnetic field: " << bz;
        // df.setBz(bz); /// put it outside the 'if'! Otherwise we have a difference wrt bz Configurable (< 1 permille) in Run2 conv. data
        // df.print();
      }
      float covMatrixPV[6];

      KFParticle::SetField(bz);
      KFPVertex kfpVertex = createKFPVertexFromCollision(collision);

      if constexpr (doPvRefit) {
        /// use PV refit
        /// Using it in the rowCandidateBase all dynamic columns shall take it into account
        // coordinates
        kfpVertex.SetXYZ(rowTrackIndexCasc3Prong.pvRefitX(), rowTrackIndexCasc3Prong.pvRefitY(), rowTrackIndexCasc3Prong.pvRefitZ());
        // covariance matrix
        kfpVertex.SetCovarianceMatrix(rowTrackIndexCasc3Prong.pvRefitSigmaX2(), rowTrackIndexCasc3Prong.pvRefitSigmaXY(), rowTrackIndexCasc3Prong.pvRefitSigmaY2(), rowTrackIndexCasc3Prong.pvRefitSigmaXZ(), rowTrackIndexCasc3Prong.pvRefitSigmaYZ(), rowTrackIndexCasc3Prong.pvRefitSigmaZ2());
      }
      kfpVertex.GetCovarianceMatrix(covMatrixPV);
      KFParticle KFPV(kfpVertex);
      hCovPVXX->Fill(covMatrixPV[0]);
      hCovPVYY->Fill(covMatrixPV[2]);
      hCovPVXZ->Fill(covMatrixPV[3]);
      hCovPVZZ->Fill(covMatrixPV[5]);

      KFPTrack kfpTrack0 = createKFPTrackFromTrack(track0);
      KFPTrack kfpTrack1 = createKFPTrackFromTrack(track1);
      KFPTrack kfpCasc = createKFPTrackFromTrack(trackCasc);

      // fill candidate table rows
      rowCandidateBase(collision.globalIndex(),
                       rowTrackIndexProng3.index0Id(), rowTrackIndexProng3.index1Id(), rowTrackIndexProng3.cascadeId(),
                       rowTrackIndexCasc3Prong.hfflag());
                     
      // fill histograms
      if (fillHistograms) {
        //hMass3->Fill(massXiPiPi);
      }
    }
  }

  void processPvRefitWithDCAFitterN(aod::Collisions const& collisions,
                                    soa::Join<aod::HfCascades, aod::HfPvRefit3Prong> const& rowsTrackIndexProng3, // HfPvRefit3Prong work?
                                    aod::V0sLinked const& v0Linked,
                                    V0Full const& v0s,
                                    CascFull const& cascs,
                                    od::TracksWCov const& tracks,
                                    aod::BCsWithTimestamps const& bcWithTimeStamps)
  {
    runCreatorCasc3ProngWithDCAFitterN<true>(collisions, rowsTrackIndexCasc3Prong, v0Linked, v0s, cascs, tracks, bcWithTimeStamps);
  }

  PROCESS_SWITCH(HfCandidateCreatorCasc3Prong, processPvRefit, "Run candidate creator with PV refit", false);

  void processNoPvRefitWithDCAFitterN(aod::Collisions const& collisions,
                                    soa::Join<aod::HfCascades, aod::HfPvRefit3Prong> const& rowsTrackIndexProng3, // HfPvRefit3Prong work?
                                    aod::V0sLinked const& v0Linked,
                                    V0Full const& v0s,
                                    CascFull const& cascs,
                                    od::TracksWCov const& tracks,
                                    aod::BCsWithTimestamps const& bcWithTimeStamps)
  {
    runCreatorCasc3ProngWithDCAFitterN<false>(collisions, rowsTrackIndexCasc3Prong, v0Linked, v0s, cascs, tracks, bcWithTimeStamps);
  }

  PROCESS_SWITCH(HfCandidateCreatorCasc3Prong, processNoPvRefit, "Run candidate creator without PV refit", true);
  
  void processPvRefitWithKFParticle(aod::Collisions const& collisions,
                                    aod::HfCascades const& rowsTrackIndexProng3, 
                                    aod::V0sLinked const& v0Linked,
                                    V0Full const& v0s,
                                    CascFull const& cascs,
                                    od::TracksWCov const& tracks,
                                    aod::BCsWithTimestamps const& bcWithTimeStamps)
  {
    runCreatorCasc3ProngWithKFParticle<true>(collisions, rowsTrackIndexCasc3Prong, v0Linked, v0s, cascs, tracks, bcWithTimeStamps);
  }
  PROCESS_SWITCH(HfCandidateCreatorCasc3Prong, processPvRefitWithKFParticle, "Run candidate creator with PV refit", false);

  void processNoPvRefitWithKFParticle(aod::Collisions const& collisions,
                                    aod::HfCascades const& rowsTrackIndexProng3,
                                    aod::V0sLinked const& v0Linked,
                                    V0Full const& v0s,
                                    CascFull const& cascs,
                                    od::TracksWCov const& tracks,
                                    aod::BCsWithTimestamps const& bcWithTimeStamps)
  {
    runCreatorCasc3ProngWithKFParticle<false>(collisions, rowsTrackIndexCasc3Prong, v0Linked, v0s, cascs, tracks, bcWithTimeStamps);
  }
  PROCESS_SWITCH(HfCandidateCreatorCasc3Prong, processNoPvRefitWithKFParticle, "Run candidate creator without PV refit", false);
};

/// Performs MC matching.
struct HfCandidateCreatorXicToXiPiPiMc {
  // TBA after finaliz candidate table
  // Spawns<aod::HfCandCasc3ProngExt> rowTrackIndexProng3;
  Produces<aod::HfXicToXiPiPiMCRec> rowMCMatchRec;
  Produces<aod::HfXicToXiPiPiMCGen> rowMCMatchGen;

  Configurable<bool> matchXicMc{"matchXicMc", true, "Do MC matching for XicPlus"};

  void init(InitContext const&) {}

  void processDoNoMc(aod::Collisions::iterator const& collision)
  {
    // dummy process function - should not be required in the future
  }
  PROCESS_SWITCH(HfCandidateCreatorToXiPiMc, processDoNoMc, "Do not run MC process function", true);

  void processMc(aod::HfCandXicToXiPiPi const& candidates,
                 aod::TracksWMc const& tracks,
                 aod::McParticles const& mcParticles)
  {
    int indexRec = -1;
    int8_t sign = -9;
    int8_t flag = -9;
    // int8_t origin = 0; //to be used for prompt/non prompt
    int8_t debug = 0;
    int8_t debugGenCharmBar = 0;
    int8_t debugGenXi = 0;
    int8_t debugGenLambda = 0;

    int pdgCodeXicPlus = pdg::Code::kXiCPlus;    // 4232
    int pdgCodeXiMinus = kXiMinus;            // 3312
    int pdgCodeLambda = kLambda0;             // 3122
    int pdgCodePiPlus = kPiPlus;              // 211
    int pdgCodePiMinus = kPiMinus;            // -211
    int pdgCodeProton = kProton;              // 2212

    // Match reconstructed candidates.
    for (const auto& candidate : candidates) {
      flag = 0;
      // origin = 0;
      debug = 0;
      auto arrayDaughters = std::array{candidate.prong0_as<aod::TracksWMc>(), // pi <- Xic
                                       candidate.prong1_as<aod::TracksWMc>(), // pi <- Xic
                                       candidate.bachelor_as<aod::TracksWMc>(),  // pi <- cascade
                                       candidate.posTrack_as<aod::TracksWMc>(),  // p <- lambda
                                       candidate.negTrack_as<aod::TracksWMc>()}; // pi <- lambda
      auto arrayDaughtersCasc = std::array{candidate.bachelor_as<aod::TracksWMc>(),
                                           candidate.posTrack_as<aod::TracksWMc>(),
                                           candidate.negTrack_as<aod::TracksWMc>()};
      auto arrayDaughtersV0 = std::array{candidate.posTrack_as<aod::TracksWMc>(),
                                         candidate.negTrack_as<aod::TracksWMc>()};

      // Xic matching
      if (matchXicMc) {
        // Xic → pi pi pi pi p
        indexRec = RecoDecay::getMatchedMCRec(mcParticles, arrayDaughters, pdgCodeXicPlus, std::array{pdgCodePiPlus, pdgCodePiPlus, pdgCodePiMinus, pdgCodeProton, pdgCodePiMinus}, true, &sign, 3);
        if (indexRec == -1) {
          debug = 1;
        }
        if (indexRec > -1) {
          // Xi- → pi pi p
          indexRec = RecoDecay::getMatchedMCRec(mcParticles, arrayDaughtersCasc, pdgCodeXiMinus, std::array{pdgCodePiMinus, pdgCodeProton, pdgCodePiMinus}, true, &sign, 2);
          if (indexRec == -1) {
            debug = 2;
          }
          if (indexRec > -1) {
            // Lambda → p pi
            indexRec = RecoDecay::getMatchedMCRec(mcParticles, arrayDaughtersV0, pdgCodeLambda, std::array{pdgCodeProton, pdgCodePiMinus}, true, &sign, 1);
            if (indexRec == -1) {
              debug = 3;
            }
            if (indexRec > -1) {
              flag = sign * (1 << aod::hf_cand_casc_lf_3prong::DecayType::XicplusToXiPiPi);
            }
          }
        }
      }

      if (debug == 2 || debug == 3) {
        LOGF(info, "WARNING: Charm baryon decays in the expected final state but the condition on the intermediate states are not fulfilled");
      }
      rowMCMatchRec(flag, debug);

    } // close loop over candidates

    // Match generated particles.
    for (const auto& particle : mcParticles) {
      flag = -9;
      sign = -9;
      debugGenCharmBar = 0;
      debugGenXi = 0;
      debugGenLambda = 0;
      // origin = 0;
      if (matchXicMc) {
        //  Xic → Xi pi pi
        if (RecoDecay::isMatchedMCGen(mcParticles, particle, pdgCodeXicPlus, std::array{pdgCodeXiMinus, pdgCodePiPlus, pdgCodePiPlus}, true, &sign)) {
          debugGenCharmBar = 1;
          // Xi- -> Lambda pi
          auto cascMC = mcParticles.rawIteratorAt(particle.daughtersIds().front());
          if (RecoDecay::isMatchedMCGen(mcParticles, cascMC, pdgCodeXiMinus, std::array{pdgCodeLambda, pdgCodePiMinus}, true)) {
            debugGenXi = 1;
            // Lambda -> p pi
            auto v0MC = mcParticles.rawIteratorAt(cascMC.daughtersIds().front());
            if (RecoDecay::isMatchedMCGen(mcParticles, v0MC, pdgCodeLambda, std::array{pdgCodeProton, pdgCodePiMinus}, true)) {
              debugGenLambda = 1;
              flag = sign * (1 << aod::hf_cand_casc_lf_3prong::DecayType::XicplusToXiPiPi);
            }
          }
        }
      }

      // rowMCMatchGen(flag, origin);
      rowMCMatchGen(flag, debugGenCharmBar, debugGenXi, debugGenLambda);
    }
  } // close process
  PROCESS_SWITCH(HfCandidateCreatorXicToXiPiPiMc, processMc, "Process MC", false);
}; // close struct

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<HfCandidateCreatorXicToXiPiPi>(cfgc),
    adaptAnalysisTask<HfCandidateCreatorXicToXiPiPiMc>(cfgc)};
}
