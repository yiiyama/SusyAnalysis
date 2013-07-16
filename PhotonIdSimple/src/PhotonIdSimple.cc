// -*- C++ -*-
//
// Package:    PhotonIdSimple
// Class:      PhotonIdSimple
// 
/**\class PhotonIdSimple PhotonIdSimple.cc SusyAnalysis/PhotonIdSimple/src/PhotonIdSimple.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Yutaro Iiyama,512 1-005,+41227670489,
//         Created:  Sun May 20 23:42:00 CEST 2012
// $Id: PhotonIdSimple.cc,v 1.4 2012/09/19 11:32:05 yiiyama Exp $
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
#include "DataFormats/Common/interface/RefToPtr.h"
#include "DataFormats/Common/interface/Association.h"

#include "DataFormats/EgammaCandidates/interface/PhotonFwd.h"
#include "DataFormats/EgammaCandidates/interface/Photon.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectronFwd.h"

#include "DataFormats/HepMCCandidate/interface/GenParticle.h"

#include "DataFormats/HLTReco/interface/TriggerEvent.h"

#include "SusyAnalysis/PhotonIdSimple/interface/TrigObjMatchFinder.h"

#include "DataFormats/Math/interface/deltaR.h"

//
// class declaration
//

class PhotonIdSimple : public edm::EDProducer {
   public:
      explicit PhotonIdSimple(const edm::ParameterSet&);
      ~PhotonIdSimple();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

   private:
      virtual void beginJob() ;
      virtual void produce(edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;
      
      virtual void beginRun(edm::Run&, edm::EventSetup const&);
      virtual void endRun(edm::Run&, edm::EventSetup const&);
      virtual void beginLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&);
      virtual void endLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&);

      // ----------member data ---------------------------

  edm::InputTag sourceTag_;
  edm::InputTag rhoSourceTag_;
  edm::InputTag mcMatchTag_;
  edm::InputTag trigEventTag_;
  std::vector<edm::InputTag> cleaningTags_;
  std::vector<edm::InputTag> bypassTags_;
  bool fiducialCut_;
  bool baseCut_;
  bool pixelSeedCut_;
  bool caloIdCut_;
  bool isoCut_;
  bool mcCut_;
  bool trigMatchCut_;
  bool trigBypassCut_;
  bool cleaningCut_;

  bool PUCorrection_;
  double ecalAeff_;
  double hcalAeff_;
  double trkAeff_;
  double matchDR_;
  double minPt_;
  double cleaningDR_;

  bool clone_;
  bool useSC_;
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
PhotonIdSimple::PhotonIdSimple(const edm::ParameterSet& iConfig) :
  sourceTag_(iConfig.getParameter<edm::InputTag>("sourceTag")),
  rhoSourceTag_(),
  mcMatchTag_(),
  trigEventTag_(),
  cleaningTags_(),
  bypassTags_(),
  fiducialCut_(false),
  baseCut_(false),
  pixelSeedCut_(false),
  caloIdCut_(false),
  isoCut_(false),
  mcCut_(false),
  trigMatchCut_(false),
  trigBypassCut_(false),
  cleaningCut_(false),
  PUCorrection_(true),
  ecalAeff_(0.),
  hcalAeff_(0.),
  trkAeff_(0.),
  matchDR_(0.),
  minPt_(20.),
  cleaningDR_(0.),
  clone_(false),
  useSC_(false),
  maxOutputSize_(-1),
  trigMatch_(0)
{
  std::vector<std::string> filters(iConfig.getParameter<std::vector<std::string> >("filters"));
  for(unsigned iF(0); iF < filters.size(); iF++){
    if(filters[iF] == "Fiducial") fiducialCut_ = true;
    if(filters[iF] == "Base") baseCut_ = true;
    if(filters[iF] == "PixelSeed") pixelSeedCut_ = true;
    if(filters[iF] == "CaloId") caloIdCut_ = true;
    if(filters[iF] == "Iso") isoCut_ = true;
    if(filters[iF] == "GenMatch") mcCut_ = true;
    if(filters[iF] == "TrigMatch") trigMatchCut_ = true;
    if(filters[iF] == "TrigBypass") trigBypassCut_ = true;
    if(filters[iF] == "Cleaning") cleaningCut_ = true;
  }

  PUCorrection_ = isoCut_;

  if(iConfig.existsAs<bool>("clone")) clone_ = iConfig.getParameter<bool>("clone");

  if(iConfig.existsAs<bool>("useSC")) useSC_ = iConfig.getParameter<bool>("useSC");

  if(iConfig.existsAs<int>("maxOutputSize")) maxOutputSize_ = iConfig.getParameter<int>("maxOutputSize");

  if(iConfig.existsAs<bool>("PUCorrection")) PUCorrection_ = iConfig.getParameter<bool>("PUCorrection");

  if(iConfig.existsAs<double>("minPt")) minPt_ = iConfig.getParameter<double>("minPt");

  if(trigBypassCut_ && iConfig.existsAs<std::vector<edm::InputTag> >("bypassTags")) bypassTags_ = iConfig.getParameter<std::vector<edm::InputTag> >("bypassTags");

  if(mcCut_){
    mcMatchTag_ = iConfig.getParameter<edm::InputTag>("mcMatchTag");
  }

  if(PUCorrection_){
    rhoSourceTag_ = iConfig.getParameter<edm::InputTag>("rhoSourceTag");
    edm::ParameterSet const& areas(iConfig.getParameterSet("effectiveAreas"));
    ecalAeff_ = areas.getParameter<double>("ecal");
    hcalAeff_ = areas.getParameter<double>("hcal");
    trkAeff_ = areas.getParameter<double>("trk");
  }

  if(trigMatchCut_ || trigBypassCut_){
    trigEventTag_ = iConfig.getParameter<edm::InputTag>("trigEventTag");
    trigMatch_ = new TrigObjMatchFinder(iConfig.getParameter<std::vector<edm::InputTag> >("trigFilterTag"));
    matchDR_ = iConfig.getParameter<double>("matchDR");
  }

  if(cleaningCut_){
    cleaningTags_ = iConfig.getParameter<std::vector<edm::InputTag> >("cleaningTags");
    cleaningDR_ = iConfig.getParameter<double>("cleaningDR");
  }

  if(clone_)
    produces<reco::PhotonCollection>();
  else
    produces<edm::RefToBaseVector<reco::Photon> >();
}


PhotonIdSimple::~PhotonIdSimple()
{
  delete trigMatch_;
}


//
// member functions
//

// ------------ method called to produce the data  ------------
void
PhotonIdSimple::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;

   RefToBaseVector<reco::Photon> const* srcRefs(0);

   Handle<View<reco::Photon> > srcHndl;
   if(iEvent.getByLabel(sourceTag_, srcHndl))
     srcRefs = &(srcHndl->refVector());
   else{
     Handle<RefToBaseVector<reco::Photon> > refSrcHndl;
     if(iEvent.getByLabel(sourceTag_, refSrcHndl))
       srcRefs = refSrcHndl.product();
     else
       throw cms::Exception("ProductNotFound") << "Source";
   }

   if(!srcRefs)
     throw cms::Exception("ProductNotFound") << "Source";

   double rho(0.);
   if(isoCut_ && PUCorrection_){
     Handle<double> rhoHndl;
     if(!iEvent.getByLabel(rhoSourceTag_, rhoHndl))
       throw cms::Exception("ProductNotFound") << rhoSourceTag_;

     rho = *(rhoHndl.product());
   }

   Association<std::vector<reco::GenParticle> > const* mcMatches(0);
   if(mcCut_){
     Handle<Association<std::vector<reco::GenParticle> > > matchHndl;
     if(!iEvent.getByLabel(mcMatchTag_, matchHndl))
       throw cms::Exception("ProductNotFound") << mcMatchTag_;

     mcMatches = matchHndl.product();
   }

   if(trigMatchCut_ || trigBypassCut_){
     Handle<trigger::TriggerEvent> teHndl;
     if(!iEvent.getByLabel(trigEventTag_, teHndl))
       throw cms::Exception("ProductNotFound") << trigEventTag_;

     trigMatch_->init(*teHndl);
   }

   std::vector<RefToBaseVector<reco::Candidate> const*> candidatesToAvoid;
   if(cleaningCut_){
     for(unsigned iC(0); iC < cleaningTags_.size(); iC++){
       Handle<View<reco::Candidate> > cHndl;
       if(iEvent.getByLabel(cleaningTags_[iC], cHndl))
         candidatesToAvoid.push_back(&(cHndl->refVector()));
       else{
         Handle<RefToBaseVector<reco::Candidate> > cHndl2;
         if(iEvent.getByLabel(cleaningTags_[iC], cHndl2)){
           candidatesToAvoid.push_back(cHndl2.product());
         }
         else
           throw cms::Exception("ProductNotFound") << cleaningTags_[iC];
       }
     }
   }

   std::vector<View<reco::Candidate> const*> extCollections;
   for(unsigned iB(0); iB < bypassTags_.size(); iB++){
     Handle<View<reco::Candidate> > colHndl;
     if(iEvent.getByLabel(bypassTags_[iB], colHndl))
       extCollections.push_back(colHndl.product());
     else
       throw cms::Exception("ProductNotFound") << bypassTags_[iB];
   }

   RefToBaseVector<reco::Photon> refs;

   for(RefToBaseVector<reco::Photon>::const_iterator phItr(srcRefs->begin()); phItr != srcRefs->end(); ++phItr){

     reco::Photon const& photon(*(phItr->get()));

     float et(useSC_ ? photon.superCluster()->rawEnergy() * std::sin(photon.superCluster()->position().theta()) : photon.et());
     if(et < minPt_) continue;

     if(trigBypassCut_){
       bool matched(false);
       for(RefToBaseVector<reco::Photon>::const_iterator it(srcRefs->begin()); it != srcRefs->end(); ++it){
         if(it == phItr) continue;
         matched |= trigMatch_->match(*(it->get()), matchDR_);
         if(matched) break;
       }
       if(!matched){
         for(unsigned iB(0); iB < extCollections.size(); iB++){
           for(View<reco::Candidate>::const_iterator cItr(extCollections[iB]->begin()); cItr != extCollections[iB]->end(); ++cItr){
             if(reco::deltaR(photon, *cItr) < matchDR_) continue;
             matched |= trigMatch_->match(*cItr, matchDR_);
             if(matched) break;
           }
           if(matched) break;
         }
       }
       if(!matched) continue;
     }

     if(trigMatchCut_)
       if(!trigMatch_->match(photon, matchDR_)) continue;

     if(fiducialCut_){
       float eta(useSC_ ? photon.superCluster()->eta() : photon.eta());
       if(std::abs(eta) > 1.4442) continue;
     }

     if(baseCut_){
       if(photon.r9() > 1.) continue;
       if(photon.hadronicOverEm() > 0.05) continue;
     }

     if(pixelSeedCut_){
       if(photon.hasPixelSeed()) continue;
     }

     if(caloIdCut_){
       if(photon.sigmaIetaIeta() > 0.011) continue;
     }

     if(isoCut_){
       float ecalIso(photon.ecalRecHitSumEtConeDR03() - ecalAeff_ * rho);
       float hcalIso(photon.hcalTowerSumEtConeDR03() - hcalAeff_ * rho);
       float trkIso(photon.trkSumPtHollowConeDR03() - trkAeff_ * rho);
       if(ecalIso + hcalIso + trkIso > 6.) continue;
     }

     if(mcCut_)
       if((*mcMatches)[*phItr].isNull()) continue;

     if(cleaningCut_){
       bool matches(false);
       for(unsigned iC(0); iC < candidatesToAvoid.size(); iC++){
         RefToBaseVector<reco::Candidate> const* refvec(candidatesToAvoid[iC]);
         for(RefToBaseVector<reco::Candidate>::const_iterator cItr(refvec->begin()); cItr != refvec->end(); ++cItr){
           if(reco::deltaR(photon, **cItr) < cleaningDR_){
             matches = true;
             break;
           }
         }
         if(matches) break;
       }
       if(matches) continue;
     }

     refs.push_back(*phItr);
   }

   if(maxOutputSize_ < refs.size()){
     RefToBaseVector<reco::Photon> refsTmp;
     refsTmp.swap(refs);
     for(unsigned iR(0); iR < maxOutputSize_; iR++)
       refs.push_back(refsTmp.at(iR));
   }

   if(clone_){
     std::auto_ptr<reco::PhotonCollection> output(new reco::PhotonCollection);
     for(unsigned iR(0); iR < refs.size(); iR++)
       output->push_back(*(refs[iR]));
     iEvent.put(output);
   }
   else{
     std::auto_ptr<RefToBaseVector<reco::Photon> > output(new RefToBaseVector<reco::Photon>);
     for(unsigned iR(0); iR < refs.size(); iR++)
       output->push_back(refs[iR]);
     iEvent.put(output);
   }
}

// ------------ method called once each job just before starting event loop  ------------
void 
PhotonIdSimple::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
PhotonIdSimple::endJob() {
}

// ------------ method called when starting to processes a run  ------------
void 
PhotonIdSimple::beginRun(edm::Run&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a run  ------------
void 
PhotonIdSimple::endRun(edm::Run&, edm::EventSetup const&)
{
}

// ------------ method called when starting to processes a luminosity block  ------------
void 
PhotonIdSimple::beginLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a luminosity block  ------------
void 
PhotonIdSimple::endLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&)
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
PhotonIdSimple::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(PhotonIdSimple);
