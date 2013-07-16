// -*- C++ -*-
//
// Package:    MuonIdSimple
// Class:      MuonIdSimple
// 
/**\class MuonIdSimple MuonIdSimple.cc SusyAnalysis/MuonIdSimple/src/MuonIdSimple.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Yutaro Iiyama,512 1-005,+41227670489,
//         Created:  Sun May 20 23:42:00 CEST 2012
// $Id: MuonIdSimple.cc,v 1.1 2012/09/19 11:32:05 yiiyama Exp $
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "FWCore/Utilities/interface/Exception.h"
#include "FWCore/Utilities/interface/InputTag.h"

#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/Common/interface/View.h"
#include "DataFormats/Common/interface/RefToBaseVector.h"
#include "DataFormats/Common/interface/Association.h"

#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"

#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonSelectors.h"

#include "DataFormats/HepMCCandidate/interface/GenParticle.h"

#include "DataFormats/HLTReco/interface/TriggerEvent.h"

#include "TrackingTools/IPTools/interface/IPTools.h"

#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"

#include "SusyAnalysis/PhotonIdSimple/interface/TrigObjMatchFinder.h"

//
// class declaration
//

class MuonIdSimple : public edm::EDProducer {
   public:
      explicit MuonIdSimple(const edm::ParameterSet&);
      ~MuonIdSimple();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

   private:
      virtual void beginJob() ;
      virtual void produce(edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;
      
      virtual void beginRun(edm::Run&, edm::EventSetup const&);
      virtual void endRun(edm::Run&, edm::EventSetup const&);
      virtual void beginLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&);
      virtual void endLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&);

  enum WP {
    none,
    loose,
    tight,
    nWP
  };

      // ----------member data ---------------------------

  edm::InputTag sourceTag_;
  edm::InputTag mcMatchTag_;
  edm::InputTag trigEventTag_;
//   edm::InputTag beamSpotTag_;
  edm::InputTag vertexTag_;
  edm::InputTag selectedMuonsTag_;

  unsigned wp_;

  double minPt_;
  double matchDR_;
  double maxIso_;

  bool barrelOnly_;
  bool clone_;
  bool mcMatch_;
  bool preselection_;
  unsigned maxOutputSize_;

  TrigObjMatchFinder* trigMatch_;
};

//
// constants, enums and typedefs
//


//
// static data member definitions
//

//
// constructors and destructor
//
MuonIdSimple::MuonIdSimple(const edm::ParameterSet& iConfig) :
  sourceTag_(iConfig.getParameter<edm::InputTag>("sourceTag")),
  mcMatchTag_(),
  trigEventTag_(),
  //  beamSpotTag_(),
  vertexTag_(),
  selectedMuonsTag_(),
  wp_(nWP),
  minPt_(10.),
  matchDR_(0.),
  maxIso_(-1.),
  barrelOnly_(false),
  clone_(false),
  mcMatch_(false),
  preselection_(false),
  maxOutputSize_(-1),
  trigMatch_(0)
{
//   if(iConfig.existsAs<edm::InputTag>("beamSpotTag"))
//     beamSpotTag_ = iConfig.getParameter<edm::InputTag>("beamSpotTag");
  if(iConfig.existsAs<edm::InputTag>("vertexTag"))
    vertexTag_ = iConfig.getParameter<edm::InputTag>("vertexTag");
  if(iConfig.existsAs<edm::InputTag>("mcMatchTag"))
    mcMatchTag_ = iConfig.getParameter<edm::InputTag>("mcMatchTag");
  if(iConfig.existsAs<bool>("barrelOnly"))
    barrelOnly_ = iConfig.getParameter<bool>("barrelOnly");
  if(iConfig.existsAs<bool>("mcMatch"))
    mcMatch_ = iConfig.getParameter<bool>("mcMatch");
  if(iConfig.existsAs<bool>("clone"))
    clone_ = iConfig.getParameter<bool>("clone");
  if(iConfig.existsAs<double>("minPt"))
    minPt_ = iConfig.getParameter<double>("minPt");
  if(iConfig.existsAs<double>("maxIso"))
    maxIso_ = iConfig.getParameter<double>("maxIso");
  if(iConfig.existsAs<int>("maxOutputSize"))
    maxOutputSize_ = iConfig.getParameter<int>("maxOutputSize");

  std::string wp(iConfig.getParameter<std::string>("wp"));
  if(wp == "loose") wp_ = loose;
  else if(wp == "tight") wp_ = tight;
  else if(wp == "none") wp_ = none;
  else throw cms::Exception("NotImplemented") << wp;

  if(iConfig.existsAs<std::vector<edm::InputTag> >("trigFilterTag")){
    trigEventTag_ = iConfig.getParameter<edm::InputTag>("trigEventTag");
    trigMatch_ = new TrigObjMatchFinder(iConfig.getParameter<std::vector<edm::InputTag> >("trigFilterTag"));
    matchDR_ = iConfig.getParameter<double>("matchDR");
  }

  if(iConfig.existsAs<edm::InputTag>("selectedMuonsTag")){
    selectedMuonsTag_ = iConfig.getParameter<edm::InputTag>("selectedMuonsTag");
    preselection_ = true;
  }

  if(clone_){
    produces<reco::MuonCollection>();
  }
  else
    produces<reco::MuonRefVector>();
}


MuonIdSimple::~MuonIdSimple()
{
  delete trigMatch_;
}


//
// member functions
//
// ------------ method called to produce the data  ------------
void
MuonIdSimple::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;

   Handle<reco::MuonCollection> srcHndl;
   if(!iEvent.getByLabel(sourceTag_, srcHndl))
     throw cms::Exception("ProductNotFound") << sourceTag_;

//    // beam spot
//    reco::BeamSpot const* beamSpot(0);
   // vertices
   Handle<reco::VertexCollection> vtxHndl;
   ESHandle<TransientTrackBuilder> trackBuilder;

   if(wp_ == tight){
//      Handle<reco::BeamSpot> bsHndl;
//      iEvent.getByLabel(beamSpotTag_, bsHndl);
//      beamSpot = bsHndl.product();

     iEvent.getByLabel(vertexTag_, vtxHndl);
     if(vtxHndl->size() == 0) return;
     iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder", trackBuilder);
   }

   // mc matching
   Association<std::vector<reco::GenParticle> > const* mcMatches(0);
   if(mcMatch_){
     Handle<Association<std::vector<reco::GenParticle> > > matchHndl;
     if(!iEvent.getByLabel(mcMatchTag_, matchHndl))
       throw cms::Exception("ProductNotFound") << mcMatchTag_;

     mcMatches = matchHndl.product();
   }

   RefToBaseVector<reco::Muon> const* preselected(0);
   if(preselection_){
     Handle<RefToBaseVector<reco::Muon> > refHndl;
     if(!iEvent.getByLabel(selectedMuonsTag_, refHndl))
       throw cms::Exception("ProductNotFound") << selectedMuonsTag_;

     preselected = refHndl.product();
   }

   if(trigMatch_){
     Handle<trigger::TriggerEvent> teHndl;
     if(!iEvent.getByLabel(trigEventTag_, teHndl))
       throw cms::Exception("ProductNotFound") << trigEventTag_;

     trigMatch_->init(*teHndl);
   }

   reco::MuonRefVector refs;

   unsigned nMu(srcHndl->size());
     //   loop on electrons
   for(unsigned iMu(0); iMu < nMu; iMu++){

     reco::MuonRef const mu(srcHndl, iMu);

     if(mu->pt() < minPt_) continue;

     if(barrelOnly_ && std::abs(mu->eta()) > 1.4442) continue;

     if(maxIso_ > 0.){
       reco::MuonPFIsolation const& pfIso(mu->pfIsolationR04());
       double relIso((pfIso.sumChargedHadronPt + std::max(0., pfIso.sumNeutralHadronEt + pfIso.sumPhotonEt - 0.5 * pfIso.sumPUPt)) / mu->pt());
       if(wp_ == loose && relIso > 0.2) continue;
       if(wp_ == tight && relIso > 0.12) continue;
     }

     if(preselection_){
       bool isInSelection(false);
       for(RefToBaseVector<reco::Muon>::const_iterator psItr(preselected->begin()); psItr != preselected->end(); ++psItr)
         if(psItr->get() == mu.get()) isInSelection = true;

       if(!isInSelection) continue;
     }

     if(trigMatch_ && !trigMatch_->match(*mu, matchDR_)) continue;

     if(wp_ != none && !mu->isPFMuon()) continue;

     if(wp_ == loose && !mu->isGlobalMuon() && !mu->isTrackerMuon()) continue;
     else if(wp_ == tight){
       if(!mu->isGlobalMuon()) continue;
       if(mu->numberOfMatchedStations() <= 1) continue;
       if(mu->innerTrack().isNull() || mu->globalTrack().isNull()) continue;
       if(mu->globalTrack()->normalizedChi2() > 10. ||
          mu->globalTrack()->hitPattern().numberOfValidMuonHits() == 0 ||
          mu->innerTrack()->hitPattern().trackerLayersWithMeasurement() <= 5 ||
          mu->innerTrack()->hitPattern().numberOfValidPixelHits() == 0) continue;

       reco::TransientTrack tt(trackBuilder->build(mu->innerTrack()));

       std::pair<bool, Measurement1D> absB(IPTools::absoluteTransverseImpactParameter(tt, vtxHndl->at(0)));
       if(!absB.first) continue;
       if(absB.second.value() > 0.2) continue;

       std::pair<bool, Measurement1D> absD(IPTools::absoluteImpactParameter3D(tt, vtxHndl->at(0)));
       if(!absD.first) continue;
       double dZ(std::sqrt(absD.second.value() * absD.second.value() - absB.second.value() * absB.second.value()));
       if(dZ > 0.5) continue;
     }

     if(mcMatch_ && (*mcMatches)[mu].isNull()) continue;

     refs.push_back(mu);
   }

   if(clone_){
     std::auto_ptr<reco::MuonCollection> output(new reco::MuonCollection);
     unsigned iOut(0);
     for(reco::MuonRefVector::const_iterator muItr(refs.begin()); muItr != refs.end() && iOut != maxOutputSize_; ++muItr, iOut++)
       output->push_back(**muItr);
     iEvent.put(output);
   }
   else{
     std::auto_ptr<reco::MuonRefVector> output(new reco::MuonRefVector);
     unsigned iOut(0);
     for(reco::MuonRefVector::const_iterator muItr(refs.begin()); muItr != refs.end() && iOut != maxOutputSize_; ++muItr, iOut++)
       output->push_back(*muItr);
     iEvent.put(output);
   }

}

// ------------ method called once each job just before starting event loop  ------------
void 
MuonIdSimple::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
MuonIdSimple::endJob() {
}

// ------------ method called when starting to processes a run  ------------
void 
MuonIdSimple::beginRun(edm::Run&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a run  ------------
void 
MuonIdSimple::endRun(edm::Run&, edm::EventSetup const&)
{
}

// ------------ method called when starting to processes a luminosity block  ------------
void 
MuonIdSimple::beginLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a luminosity block  ------------
void 
MuonIdSimple::endLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&)
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
MuonIdSimple::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(MuonIdSimple);
