// -*- C++ -*-
//
// Package:    ElectronIdSimple
// Class:      ElectronIdSimple
// 
/**\class ElectronIdSimple ElectronIdSimple.cc SusyAnalysis/ElectronIdSimple/src/ElectronIdSimple.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Yutaro Iiyama,512 1-005,+41227670489,
//         Created:  Sun May 20 23:42:00 CEST 2012
// $Id: ElectronIdSimple.cc,v 1.4 2012/09/19 11:32:05 yiiyama Exp $
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

#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectronFwd.h"
#include "DataFormats/EgammaCandidates/interface/Conversion.h"

#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"

#include "FWCore/Utilities/interface/Exception.h"
#include "FWCore/Utilities/interface/InputTag.h"

#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/Common/interface/ValueMap.h"
#include "DataFormats/Common/interface/View.h"
#include "DataFormats/Common/interface/RefToBaseVector.h"
#include "DataFormats/Common/interface/Association.h"

#include "DataFormats/HepMCCandidate/interface/GenParticle.h"

#include "EGamma/EGammaAnalysisTools/interface/EGammaCutBasedEleId.h"

#include "SusyAnalysis/PhotonIdSimple/interface/TrigObjMatchFinder.h"

//
// class declaration
//

class ElectronIdSimple : public edm::EDProducer {
   public:
      explicit ElectronIdSimple(const edm::ParameterSet&);
      ~ElectronIdSimple();

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
    veto = EgammaCutBasedEleId::VETO,
    loose = EgammaCutBasedEleId::LOOSE,
    medium = EgammaCutBasedEleId::MEDIUM,
    tight = EgammaCutBasedEleId::TIGHT,
    nWP
  };

      // ----------member data ---------------------------

  edm::InputTag sourceTag_;
  edm::InputTag conversionTag_;
  edm::InputTag rhoTag_;
  edm::InputTag beamSpotTag_;
  edm::InputTag vertexTag_;
  std::vector<edm::InputTag> isoValTags_;
  edm::InputTag mcMatchTag_;
  edm::InputTag selectedElectronsTag_;
  edm::InputTag trigEventTag_;
  unsigned wp_;

  bool barrelOnly_;
  bool mcMatch_;
  bool clone_;
  bool preselection_;

  unsigned maxOutputSize_;
  double matchDR_;
  double minPt_;

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
ElectronIdSimple::ElectronIdSimple(const edm::ParameterSet& iConfig) :
  sourceTag_(iConfig.getParameter<edm::InputTag>("sourceTag")),
  conversionTag_(),
  rhoTag_(),
  beamSpotTag_(),
  vertexTag_(),
  isoValTags_(),
  mcMatchTag_(),
  selectedElectronsTag_(),
  trigEventTag_(),
  barrelOnly_(false),
  mcMatch_(false),
  clone_(false),
  preselection_(false),
  maxOutputSize_(-1),
  matchDR_(0.),
  minPt_(36.),
  trigMatch_(0)
{
  if(iConfig.existsAs<edm::InputTag>("conversionTag"))
    conversionTag_ = iConfig.getParameter<edm::InputTag>("conversionTag");
  if(iConfig.existsAs<edm::InputTag>("rhoTag"))
    rhoTag_ = iConfig.getParameter<edm::InputTag>("rhoTag");
  if(iConfig.existsAs<edm::InputTag>("beamSpotTag"))
    beamSpotTag_ = iConfig.getParameter<edm::InputTag>("beamSpotTag");
  if(iConfig.existsAs<edm::InputTag>("vertexTag"))
    vertexTag_ = iConfig.getParameter<edm::InputTag>("vertexTag");
  if(iConfig.existsAs<std::vector<edm::InputTag> >("isoValTags"))
    isoValTags_ =iConfig.getParameter<std::vector<edm::InputTag> >("isoValTags");
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
  if(iConfig.existsAs<int>("maxOutputSize"))
    maxOutputSize_ = iConfig.getParameter<int>("maxOutputSize");

  std::string wp(iConfig.getParameter<std::string>("wp"));
  if(wp == "veto") wp_ = veto;
  else if(wp == "loose") wp_ = loose;
  else if(wp == "medium") wp_ = medium;
  else if(wp == "tight") wp_ = tight;
  else if(wp == "none") wp_ = none;
  else throw cms::Exception("NotImplemented") << wp;

  if(iConfig.existsAs<edm::InputTag>("selectedElectronsTag")){
    selectedElectronsTag_ = iConfig.getParameter<edm::InputTag>("selectedElectronsTag");
    preselection_ = true;
  }

  if(iConfig.existsAs<std::vector<edm::InputTag> >("trigFilterTag")){
    trigEventTag_ = iConfig.getParameter<edm::InputTag>("trigEventTag");
    trigMatch_ = new TrigObjMatchFinder(iConfig.getParameter<std::vector<edm::InputTag> >("trigFilterTag"));
    matchDR_ = iConfig.getParameter<double>("matchDR");
  }

  if(clone_){
    produces<reco::GsfElectronCollection>();
  }
  else
    produces<reco::GsfElectronRefVector>();
}


ElectronIdSimple::~ElectronIdSimple()
{
  delete trigMatch_;
}


//
// member functions
//

// ------------ method called to produce the data  ------------
void
ElectronIdSimple::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;

   Handle<reco::GsfElectronCollection> srcHndl;
   if(!iEvent.getByLabel(sourceTag_, srcHndl))
     throw cms::Exception("ProductNotFound") << sourceTag_;

//    reco::GsfElectronRefVector* srcRefs(0);

//    Handle<reco::GsfElectronCollection> srcHndl;
//    if(iEvent.getByLabel(sourceTag_, srcHndl)){
//      srcRefs = new reco::GsfElectronRefVector;
//      for(unsigned iS(0); iS < srcHndl->size(); iS++)
//        srcRefs->push_back(Ref<reco::GsfElectronCollection>(srcHndl.product(), iS));
//    }
//    else{
//      Handle<reco::GsfElectronRefVector> refSrcHndl;
//      if(iEvent.getByLabel(sourceTag_, refSrcHndl))
//        srcRefs = new reco::GsfElectronRefVector(*refSrcHndl.product());
//      else{
//        Handle<RefToBaseVector<reco::GsfElectron> > baseRefSrcHndl;
//        if(iEvent.getByLabel(sourceTag_, baseRefSrcHndl)){
//          srcRefs = new reco::GsfElectronRefVector;
//          for(RefToBaseVector<reco::GsfElectron>::const_iterator rItr(baseRefSrcHndl->begin()); rItr != baseRefSrcHndl->end(); ++rItr)
//            srcRefs->push_back(rItr->castTo<Ref<reco::GsfElectronCollection> >());
//        }
//        else
//          throw cms::Exception("ProductNotFound") << sourceTag_;
//      }
//    }

   // conversions
   Handle<reco::ConversionCollection> cnvHndl;
   // iso deposits
   ValueMap<double> const* chargedHadIso(0);
   ValueMap<double> const* emIso(0);
   ValueMap<double> const* neutralHadIso(0);
   // beam spot
   reco::BeamSpot const* beamSpot(0);
   // vertices
   Handle<reco::VertexCollection> vtxHndl;
   // rho for isolation
   double rho;

   if(wp_ != none){
     Handle<ValueMap<double> > chIsoHndl, emIsoHndl, nhIsoHndl;
     if(!iEvent.getByLabel(isoValTags_[0], chIsoHndl))
       throw cms::Exception("ProductNotFound") << isoValTags_[0];
     if(!iEvent.getByLabel(isoValTags_[1], emIsoHndl))
       throw cms::Exception("ProductNotFound") << isoValTags_[1];
     if(!iEvent.getByLabel(isoValTags_[2], nhIsoHndl))
       throw cms::Exception("ProductNotFound") << isoValTags_[2];
     chargedHadIso = chIsoHndl.product();
     emIso = emIsoHndl.product();
     neutralHadIso = nhIsoHndl.product();

     if(!iEvent.getByLabel(conversionTag_, cnvHndl))
       throw cms::Exception("ProductNotFound") << conversionTag_;

     Handle<reco::BeamSpot> bsHndl;
     iEvent.getByLabel(beamSpotTag_, bsHndl);
     beamSpot = bsHndl.product();

     iEvent.getByLabel(vertexTag_, vtxHndl);

     Handle<double> rhoHndl;
     iEvent.getByLabel(rhoTag_, rhoHndl);
     rho = *(rhoHndl.product());
   }

   // mc matching
   Association<std::vector<reco::GenParticle> > const* mcMatches(0);
   if(mcMatch_){
     Handle<Association<std::vector<reco::GenParticle> > > matchHndl;
     if(!iEvent.getByLabel(mcMatchTag_, matchHndl))
       throw cms::Exception("ProductNotFound") << mcMatchTag_;

     mcMatches = matchHndl.product();
   }

   RefToBaseVector<reco::GsfElectron> const* preselected(0);
   if(preselection_){
     Handle<RefToBaseVector<reco::GsfElectron> > refHndl;
     if(!iEvent.getByLabel(selectedElectronsTag_, refHndl))
       throw cms::Exception("ProductNotFound") << selectedElectronsTag_;

     preselected = refHndl.product();
   }

   if(trigMatch_){
     Handle<trigger::TriggerEvent> teHndl;
     if(!iEvent.getByLabel(trigEventTag_, teHndl))
       throw cms::Exception("ProductNotFound") << trigEventTag_;

     trigMatch_->init(*teHndl);
   }

   reco::GsfElectronRefVector refs;

   unsigned nE(srcHndl->size());
     //   loop on electrons
   for(unsigned iE(0); iE < nE; iE++){
       //   for(reco::GsfElectronRefVector::const_iterator eItr(srcRefs->begin()); eItr != srcRefs->end(); ++eItr){

     // get reference to electron
     reco::GsfElectronRef const ele(srcHndl, iE);

     if(preselection_){
       bool isInSelection(false);
       for(RefToBaseVector<reco::GsfElectron>::const_iterator psItr(preselected->begin()); psItr != preselected->end(); ++psItr)
         if(psItr->get() == ele.get()) isInSelection = true;

       if(!isInSelection) continue;
     }

     if(ele->pt() < minPt_) continue;

     if(barrelOnly_ && std::abs(ele->eta()) > 1.4442) continue;

     if(trigMatch_ && !trigMatch_->match(*ele, matchDR_)) continue;

     // working points
     bool pass(true);

     if(wp_ != none){
       double chIsoVal(0.), emIsoVal(0.), nhIsoVal(0.);

       chIsoVal = (*chargedHadIso)[ele];
       emIsoVal = (*emIso)[ele];
       nhIsoVal = (*neutralHadIso)[ele];

       pass = EgammaCutBasedEleId::PassWP(EgammaCutBasedEleId::WorkingPoint(wp_), ele, cnvHndl, *beamSpot, vtxHndl, chIsoVal, emIsoVal, nhIsoVal, rho);
     }

     if(mcMatch_ && pass)
       pass = (*mcMatches)[ele].isNonnull();

     if(pass)
       refs.push_back(ele);
   }

   if(clone_){
     std::auto_ptr<reco::GsfElectronCollection> output(new reco::GsfElectronCollection);
     unsigned iOut(0);
     for(reco::GsfElectronRefVector::const_iterator eItr(refs.begin()); eItr != refs.end() && iOut != maxOutputSize_; ++eItr, iOut++)
       output->push_back(**eItr);
     iEvent.put(output);
   }
   else{
     std::auto_ptr<reco::GsfElectronRefVector> output(new reco::GsfElectronRefVector);
     unsigned iOut(0);
     for(reco::GsfElectronRefVector::const_iterator eItr(refs.begin()); eItr != refs.end() && iOut != maxOutputSize_; ++eItr, iOut++)
       output->push_back(*eItr);
     iEvent.put(output);
   }

   //   delete srcRefs;
}

// ------------ method called once each job just before starting event loop  ------------
void 
ElectronIdSimple::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
ElectronIdSimple::endJob() {
}

// ------------ method called when starting to processes a run  ------------
void 
ElectronIdSimple::beginRun(edm::Run&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a run  ------------
void 
ElectronIdSimple::endRun(edm::Run&, edm::EventSetup const&)
{
}

// ------------ method called when starting to processes a luminosity block  ------------
void 
ElectronIdSimple::beginLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a luminosity block  ------------
void 
ElectronIdSimple::endLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&)
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
ElectronIdSimple::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(ElectronIdSimple);
