#ifndef TrigObjMatchFinder_H
#define TrigObjMatchFinder_H

#include <vector>

#include "FWCore/Utilities/interface/InputTag.h"

#include "DataFormats/HLTReco/interface/TriggerEvent.h"
#include "DataFormats/HLTReco/interface/TriggerObject.h"

#include "DataFormats/RecoCandidate/interface/RecoCandidate.h"

class TrigObjMatchFinder {
 public:
  TrigObjMatchFinder(std::vector<edm::InputTag> const& _tags) : filterTags_(_tags), filterProducts_() {}
  void init(trigger::TriggerEvent const&);
  bool match(reco::Candidate const&, double);

 private:
  std::vector<edm::InputTag> filterTags_;
  std::vector<trigger::TriggerObjectCollection> filterProducts_;
};

#endif
