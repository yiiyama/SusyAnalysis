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

#include <limits>

class PhotonIdSimple : public edm::EDProducer {
public:
  explicit PhotonIdSimple(const edm::ParameterSet&);
  ~PhotonIdSimple();

private:
  enum WorkingPoint {
    kLoose,
    kMedium,
    kTight,
    nWorkingPoints
  };

  enum Cut {
    kFiducial,
    kElectronVeto,
    kHOverE,
    kShowerShape,
    kCHIsolation,
    kNHIsolation,
    kPhIsolation,
    kGenMatch,
    kTriggerMatch,
    kTriggerBypass,
    kCleaning,
    nCuts
  };

  virtual void produce(edm::Event&, const edm::EventSetup&);
      
  edm::InputTag sourceTag_;
  edm::InputTag rhoSourceTag_;
  edm::InputTag genMatchTag_;
  edm::InputTag trigEventTag_;
  edm::InputTag conversionTag_;
  edm::InputTag beamSpotTag_;
  edm::InputTag electronTag_;
  std::vector<edm::InputTag> cleaningTags_;
  std::vector<edm::InputTag> bypassTags_;

  WorkingPoint wp_;
  bool cuts_[nCuts];

  bool pfBasedIso_;
  bool PUCorrection_;
  double effArea_[3];
  double matchDR_;
  double minPt_;
  double cleaningDR_;

  bool clone_;
  bool usePixelVeto_;
  bool useSC_;
  unsigned maxOutputSize_;

  TrigObjMatchFinder* trigMatch_;
};

PhotonIdSimple::PhotonIdSimple(const edm::ParameterSet& _ps) :
  sourceTag_(_ps.getParameter<edm::InputTag>("sourceTag")),
  rhoSourceTag_(),
  genMatchTag_(),
  trigEventTag_(),
  conversionTag_(),
  beamSpotTag_(),
  electronTag_(),
  cleaningTags_(),
  bypassTags_(),
  wp_(nWorkingPoints),
  pfBasedIso_(true),
  doPUCorrection_(true),
  matchDR_(0.),
  minPt_(20.),
  cleaningDR_(0.),
  clone_(false),
  usePixelVeto_(false),
  useSC_(false),
  maxOutputSize_(-1),
  trigMatch_(0)
{
  std::string workingPoint(_ps.getParameter<std::string>("workingPoint"));
  if(workingPoint == "Loose") wp_ = kLoose;
  else if(workingPoint == "Medium") wp_ = kMedium;
  else if(workingPoint == "Tight") wp_ = kTight;
  else throw cms::Exception("InvalidConfiguration") << "Undefined working point " + workingPoint;

  std::fill_n(cuts_, nCuts, false);
  std::vector<std::string> filters(_ps.getParameter<std::vector<std::string> >("filters"));
  for(unsigned iF(0); iF < filters.size(); iF++){
    if(filters[iF] == "Fiducial") cuts_[kFiducial] = true;
    else if(filters[iF] == "ElectronVeto") cuts_[kElectronVeto] = true;
    else if(filters[iF] == "HOverE") cuts_[kHOverE] = true;
    else if(filters[iF] == "ShowerShape") cuts_[kShowerShape] = true;
    else if(filters[iF] == "CHIsolation") cuts_[kCHIsolation] = true;
    else if(filters[iF] == "NHIsolation") cuts_[kNHIsolation] = true;
    else if(filters[iF] == "PhIsolation") cuts_[kPhIsolation] = true;
    else if(filters[iF] == "GenMatch") cuts_[kGenMatch] = true;
    else if(filters[iF] == "TriggerMatch") cuts_[kTriggerMatch] = true;
    else if(filters[iF] == "TriggerBypass") cuts_[kTriggerBypass] = true;
    else if(filters[iF] == "Cleaning") cuts_[kCleaning] = true;
    else throw cms::Exception("InvalidConfiguration") << "Undefined cut " + filters[iF];
  }

  if(_ps.existsAs<bool>("clone")) clone_ = _ps.getParameter<bool>("clone");

  if(_ps.existsAs<bool>("usePixelVeto")) usePixelVeto_ = _ps.getParameter<bool>("usePixelVeto");

  if(_ps.existsAs<bool>("useSC")) useSC_ = _ps.getParameter<bool>("useSC");

  if(_ps.existsAs<int>("maxOutputSize")) maxOutputSize_ = _ps.getParameter<int>("maxOutputSize");

  if(_ps.existsAs<bool>("pfBasedIso")) pfBasedIso_ = _ps.getParameter<bool>("pfBasedIso");

  // PU correction on isolation is default ON
  if(_ps.existsAs<bool>("doPUCorrection")) doPUCorrection_ = _ps.getParameter<bool>("doPUCorrection");

  if(_ps.existsAs<double>("minPt")) minPt_ = _ps.getParameter<double>("minPt");

  if(cuts_[kElectronVeto] && !usePixelVeto_){
    conversionTag_ = _ps.getParameter<edm::InputTag>("conversionTag");
    electronTag_ = _ps.getParameter<edm::InputTag>("electronTag");
    beamSpotTag_ = _ps.getParameter<edm::InputTag>("beamSpotTag");
  }

  if(cuts_[kTriggerBypass]) bypassTags_ = _ps.getParameter<std::vector<edm::InputTag> >("bypassTags");

  if(cuts_[kGenMatch]) genMatchTag_ = _ps.getParameter<edm::InputTag>("genMatchTag");

  std::fill_n(effArea_, 3, 0.);
  if(doPUCorrection_){
    rhoSourceTag_ = _ps.getParameter<edm::InputTag>("rhoSourceTag");
    edm::ParameterSet const& areas(_ps.getParameterSet("effectiveAreas"));
    if(pfBasedIso_){
      effArea_[0] = areas.getParameter<double>("chargedHadron");
      effArea_[1] = areas.getParameter<double>("neutralHadron");
      effArea_[2] = areas.getParameter<double>("photon");
    }
    else{
      throw cms::Exception("NotSupported");
      effArea_[0] = areas.getParameter<double>("tracker");
      effArea_[1] = areas.getParameter<double>("hcal");
      effArea_[2] = areas.getParameter<double>("ecal");
    }
  }

  if(cuts_[kTriggerMatch] || cuts_[kTriggerBypass]){
    trigEventTag_ = _ps.getParameter<edm::InputTag>("trigEventTag");
    trigMatch_ = new TrigObjMatchFinder(_ps.getParameter<std::vector<edm::InputTag> >("trigFilterTag"));
    matchDR_ = _ps.getParameter<double>("matchDR");
  }

  if(cuts_[kCleaning]){
    cleaningTags_ = _ps.getParameter<std::vector<edm::InputTag> >("cleaningTags");
    cleaningDR_ = _ps.getParameter<double>("cleaningDR");
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

void
PhotonIdSimple::produce(edm::Event& _event, const edm::EventSetup&)
{
  edm::RefToBaseVector<reco::Photon> const* srcRefs(0);

  edm::Handle<View<reco::Photon> > srcHndl;
  if(_event.getByLabel(sourceTag_, srcHndl))
    srcRefs = &(srcHndl->refVector());
  else{
    edm::Handle<edm::RefToBaseVector<reco::Photon> > refSrcHndl;
    if(_event.getByLabel(sourceTag_, refSrcHndl))
      srcRefs = refSrcHndl.product();
    else
      throw cms::Exception("ProductNotFound") << "Source";
  }

  if(!srcRefs)
    throw cms::Exception("ProductNotFound") << "Source";

  double rho(0.);
  if(doPUCorrection_){
    Handle<double> rhoHndl;
    if(!_event.getByLabel(rhoSourceTag_, rhoHndl))
      throw cms::Exception("ProductNotFound") << rhoSourceTag_;

    rho = *(rhoHndl.product());
  }

  edm::Association<std::vector<reco::GenParticle> > const* genMatches(0);
  if(cuts_[kGenMatch]){
    edm::Handle<edm::Association<std::vector<reco::GenParticle> > > matchHndl;
    if(!_event.getByLabel(genMatchTag_, matchHndl))
      throw cms::Exception("ProductNotFound") << genMatchTag_;

    genMatches = matchHndl.product();
  }

  if(cuts_[kTriggerMatch] || cuts_[kTriggerBypass]){
    edm::Handle<trigger::TriggerEvent> teHndl;
    if(!_event.getByLabel(trigEventTag_, teHndl))
      throw cms::Exception("ProductNotFound") << trigEventTag_;

    trigMatch_->init(*teHndl);
  }

  std::vector<edm::RefToBaseVector<reco::Candidate> const*> candidatesToAvoid;
  if(cuts_[kCleaning]){
    for(unsigned iC(0); iC < cleaningTags_.size(); iC++){
      edm::Handle<View<reco::Candidate> > cHndl;
      if(_event.getByLabel(cleaningTags_[iC], cHndl))
        candidatesToAvoid.push_back(&(cHndl->refVector()));
      else{
        edm::Handle<RefToBaseVector<reco::Candidate> > cHndl2;
        if(_event.getByLabel(cleaningTags_[iC], cHndl2)){
          candidatesToAvoid.push_back(cHndl2.product());
        }
        else
          throw cms::Exception("ProductNotFound") << cleaningTags_[iC];
      }
    }
  }

  std::vector<edm::View<reco::Candidate> const*> extCollections;
  for(unsigned iB(0); iB < bypassTags_.size(); iB++){
    edm::Handle<edm::View<reco::Candidate> > colHndl;
    if(_event.getByLabel(bypassTags_[iB], colHndl))
      extCollections.push_back(colHndl.product());
    else
      throw cms::Exception("ProductNotFound") << bypassTags_[iB];
  }

  edm::Handle<reco::ConversionCollection> conversionHndl;
  edm::Handle<reco::GsfElectronCollection> electonHndl;
  reco::BeamSpot const* beamSpot(0);
  if(cuts_[kElectronVeto] && !usePixelVeto_){
    _event.getByLabel(conversionTag_, conversionHndl);
    _event.getByLabel(electronTag_, electronHndl);
    edm::Handle<reco::BeamSpot> beamSpotHndl;
    _event.getByLabel(beamSpotTag_, beamSpotHndl);
    beamSpot = *beamSpotHndl;
  }

  double cutValues[nWorkingPoints][nCuts][2];

  for(unsigned iW(0); iW != nWorkingPoints; ++iW)
    cutValues[iW][kHOverE][kBarrel] = cutValues[iW][kHOverE][kEndcap] = 0.05;

  cutValues[kLoose][kShowerShape][kBarrel] = 0.012;
  cutValues[kMedium][kShowerShape][kBarrel] = 0.011;
  cutValues[kTight][kShowerShape][kBarrel] = 0.011;

  cutValues[kLoose][kShowerShape][kEndcap] = 0.034;
  cutValues[kMedium][kShowerShape][kEndcap] = 0.033;
  cutValues[kTight][kShowerShape][kEndcap] = 0.031;

  if(pfBasedIso_){
    cutValues[kLoose][kCHIsolation][kBarrel] = 2.6;
    cutValues[kMedium][kCHIsolation][kBarrel] = 1.5;
    cutValues[kTight][kCHIsolation][kBarrel] = 0.7;
    cutValues[kLoose][kNHIsolation][kBarrel] = 3.5;
    cutValues[kMedium][kNHIsolation][kBarrel] = 1.0;
    cutValues[kTight][kNHIsolation][kBarrel] = 0.4;
    cutValues[kLoose][kPhIsolation][kBarrel] = 1.3;
    cutValues[kMedium][kPhIsolation][kBarrel] = 0.7;
    cutValues[kTight][kPhIsolation][kBarrel] = 0.5;

    cutValues[kLoose][kCHIsolation][kEndcap] = 2.3;
    cutValues[kMedium][kCHIsolation][kEndcap] = 1.2;
    cutValues[kTight][kCHIsolation][kEndcap] = 0.5;
    cutValues[kLoose][kNHIsolation][kEndcap] = 2.9;
    cutValues[kMedium][kNHIsolation][kEndcap] = 1.5;
    cutValues[kTight][kNHIsolation][kEndcap] = 1.5;
    cutValues[kLoose][kPhIsolation][kEndcap] = std::numeric_limits<double>::max();
    cutValues[kMedium][kPhIsolation][kEndcap] = 1.0;
    cutValues[kTight][kPhIsolation][kEndcap] = 1.0;
  }

  edm::RefToBaseVector<reco::Photon> refs;

  // Loop over photons
  for(edm::RefToBaseVector<reco::Photon>::const_iterator phItr(srcRefs->begin()); phItr != srcRefs->end(); ++phItr){

    reco::Photon const& photon(*(phItr->get()));

    double et(useSC_ ? photon.superCluster()->rawEnergy() * std::sin(photon.superCluster()->position().theta()) : photon.et());
    if(et < minPt_) continue;

    if(cuts_[kTriggerBypass]){
      edm::RefToBaseVector<reco::Photon>::const_iterator it(srcRefs->begin());
      for(; it != srcRefs->end(); ++it){
        if(it == phItr) continue;
        if(trigMatch_->match(*(it->get()), matchDR_));
      }
      if(it == srcRefs->end()){ // no matching reco photon was found for the given trigger filter
        // is there something offline that matches the trigger object? (e.g. electrons)
        unsigned iB(0);
        for(; iB < extCollections.size(); iB++){
          edm::View<reco::Candidate>::const_iterator cItr(extCollections[iB]->begin());
          for(; cItr != extCollections[iB]->end(); ++cItr){
            if(reco::deltaR(photon, *cItr) < matchDR_) continue;
            if(trigMatch_->match(*cItr, matchDR_)) break;
          }
          if(cItr != extCollections[iB]->end()) break;
        }
        if(iB == extCollections.size()) continue; // no offline object matched the trigger filter
      }
    }

    if(cuts_[kTriggerMatch] && !trigMatch_->match(photon, matchDR_)) continue;

    if(cuts_[kElectronVeto]){
      if(usePixelVeto_){
        if(photon.hasPixelSeed()) continue;
      }
      else{
        ConversionTools::hasMatchedPromptElectron(photon.superCluster(), electronHndl, conversionHndl, beamSpot) continue;
      }
    }

    double absEta(std::abs(photon.superCluster()->eta()));

    enum Subdetector {
      kBarrel,
      kEndcap,
      kGap
    };

    Subdetector subdet(kGap);
    if(absEta < 1.4442) subdet = kBarrel; // barrel
    else if(absEta < 1.566) subdet = kGap;
    else if(absEta < 2.5) subdet = kEndcap; // endcap

    if(cuts_[kFiducial] && subdet == kGap) continue;

    if(cuts_[kHOverE] && photon.hadTowOverEm() > cutValues[wp_][kHOverE][subdet]) continue;

    if(cuts_[kShowerShape] && photon.sigmaIetaIeta() > cutValues[wp_][kShowerShape][subdet]) continue;

    // TODO HERE PFISOLATOR

    if(cuts_[kCHIsolation]){
      float ecalIso(photon.ecalRecHitSumEtConeDR03() - ecalAeff_ * rho);
      float hcalIso(photon.hcalTowerSumEtConeDR03() - hcalAeff_ * rho);
      float trkIso(photon.trkSumPtHollowConeDR03() - trkAeff_ * rho);
      if(ecalIso + hcalIso + trkIso > 6.) continue;
    }

    if(mcCut_)
      if((*genMatches)[*phItr].isNull()) continue;

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
    _event.put(output);
  }
  else{
    std::auto_ptr<RefToBaseVector<reco::Photon> > output(new RefToBaseVector<reco::Photon>);
    for(unsigned iR(0); iR < refs.size(); iR++)
      output->push_back(refs[iR]);
    _event.put(output);
  }
}

//define this as a plug-in
DEFINE_FWK_MODULE(PhotonIdSimple);


