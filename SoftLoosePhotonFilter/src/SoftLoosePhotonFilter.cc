#include "../interface/SoftLoosePhotonFilter.h"

// system include files
#include <memory>
#include <cmath>

// user include files
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/Utilities/interface/Exception.h"

#include "DataFormats/Common/interface/Handle.h"

#include "DataFormats/EgammaCandidates/interface/PhotonFwd.h"
#include "DataFormats/EgammaCandidates/interface/Photon.h"
#include "DataFormats/EgammaCandidates/interface/ConversionFwd.h"

#include "DataFormats/BeamSpot/interface/BeamSpot.h"

//#include "DataFormats/ParticleFlowCandidate/interface/PFCandidateFwd.h"

//#include "DataFormats/VertexReco/interface/VertexFwd.h"

#include "RecoEgamma/EgammaTools/interface/ConversionTools.h"

//#include "EGamma/EGammaAnalysisTools/interface/PFIsolationEstimator.h"

#include "TString.h"
#include "TObjArray.h"

SoftLoosePhotonFilter::SoftLoosePhotonFilter(const edm::ParameterSet& iConfig) :
  photonCollectionTag_(iConfig.getUntrackedParameter<edm::InputTag>("photonCollectionTag")),
  electronCollectionTag_(iConfig.getUntrackedParameter<edm::InputTag>("electronCollectionTag")),
  conversionCollectionTag_(iConfig.getUntrackedParameter<edm::InputTag>("conversionCollectionTag")),
  beamSpotTag_(iConfig.getUntrackedParameter<edm::InputTag>("beamSpotTag")),
  //  pfCandidateCollectionTag_(iConfig.getUntrackedParameter<edm::InputTag>("pfCandidateCollectionTag")),
  //  vertexCollectionTag_(iConfig.getUntrackedParameter<edm::InputTag>("vertexCollectionTag")),
  //  isolator_(),
  ptThreshold_(iConfig.getUntrackedParameter<double>("ptThreshold")),
  doCut_(nCutType, false)
{
  std::vector<int> cuts(iConfig.getUntrackedParameter<std::vector<int> >("cuts"));
  for(unsigned iC(0); iC != cuts.size(); ++iC)
    if(cuts[iC] < nCutType) doCut_[cuts[iC]] = true;

//   isolator_.initializePhotonIsolation(true);
//   isolator_.setConeSize(0.3);
}


SoftLoosePhotonFilter::~SoftLoosePhotonFilter()
{
}

bool
SoftLoosePhotonFilter::filter(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  edm::Handle<reco::PhotonCollection> photons;
  if(!iEvent.getByLabel(photonCollectionTag_, photons))
    throw cms::Exception("ProductNotFound") << photonCollectionTag_;

  edm::Handle<reco::GsfElectronCollection> electrons;
  edm::Handle<reco::ConversionCollection> conversions;
  edm::Handle<reco::BeamSpot> beamSpot;
//   edm::Handle<reco::PFCandidateCollection> pfCandidates;
//   edm::Handle<reco::VertexCollection> vertices;
  if(!iEvent.getByLabel(electronCollectionTag_, electrons))
    throw cms::Exception("ProductNotFound") << electronCollectionTag_;
  if(!iEvent.getByLabel(conversionCollectionTag_, conversions))
    throw cms::Exception("ProductNotFound") << conversionCollectionTag_;
  if(!iEvent.getByLabel(beamSpotTag_, beamSpot)) 
    throw cms::Exception("ProductNotFound") << beamSpotTag_;
//     if(!iEvent.getByLabel(pfCandidateCollectionTag_, pfCandidates)) return false;
//     if(!iEvent.getByLabel(vertexCollectionTag_, vertices)) return false;

  reco::PhotonCollection::const_iterator phItr(photons->begin());
  reco::PhotonCollection::const_iterator phEnd(photons->end());
  for(; phItr != phEnd; ++phItr){
    if(phItr->pt() < ptThreshold_) continue;

    bool isEB(std::abs(phItr->eta()) < 1.4442);
    bool isEE(std::abs(phItr->eta()) > 1.566);
    if(!isEB && !isEE) continue;

    if(doCut_[HOverELoose] && phItr->hadTowOverEm() > 0.1) continue;
    if(doCut_[ElectronVeto] && ConversionTools::hasMatchedPromptElectron(phItr->superCluster(), electrons, conversions, beamSpot->position())) continue;
    if(doCut_[ShowerShapeLoose] && phItr->sigmaIetaIeta() > (isEB ? 0.014 : 0.035)) continue;

    break;
  }

  return phItr != phEnd;
}

// ------------ method called once each job just before starting event loop  ------------
void 
SoftLoosePhotonFilter::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
SoftLoosePhotonFilter::endJob() {
}

// ------------ method called when starting to processes a run  ------------
bool 
SoftLoosePhotonFilter::beginRun(edm::Run&, edm::EventSetup const&)
{ 
  return true;
}

// ------------ method called when ending the processing of a run  ------------
bool 
SoftLoosePhotonFilter::endRun(edm::Run&, edm::EventSetup const&)
{
  return true;
}

// ------------ method called when starting to processes a luminosity block  ------------
bool 
SoftLoosePhotonFilter::beginLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&)
{
  return true;
}

// ------------ method called when ending the processing of a luminosity block  ------------
bool 
SoftLoosePhotonFilter::endLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&)
{
  return true;
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
SoftLoosePhotonFilter::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}
//define this as a plug-in
DEFINE_FWK_MODULE(SoftLoosePhotonFilter);
