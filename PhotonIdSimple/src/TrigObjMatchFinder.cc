#include "SusyAnalysis/PhotonIdSimple/interface/TrigObjMatchFinder.h"

#include "DataFormats/Math/interface/deltaR.h"

void
TrigObjMatchFinder::init(trigger::TriggerEvent const& trigEvent)
{
  filterProducts_.clear();
  filterProducts_.resize(filterTags_.size(), trigger::TriggerObjectCollection());

  trigger::TriggerObjectCollection const& toCollection(trigEvent.getObjects());

  for(unsigned iT(0); iT < filterTags_.size(); iT++){
    edm::InputTag& filterTag(filterTags_[iT]);

    for(unsigned iF(0); iF < trigEvent.sizeFilters(); iF++){

      if(filterTag == trigEvent.filterTag(iF)){

        trigger::Keys const& keys(trigEvent.filterKeys(iF));
        for(unsigned iK(0); iK < keys.size(); iK++)
          filterProducts_[iT].push_back(toCollection.at(keys[iK]));

        break;

      }
    
    }

  }
}

bool
TrigObjMatchFinder::match(reco::Candidate const& _cand, double _dR)
{
  for(unsigned iT(0); iT < filterProducts_.size(); iT++){
    trigger::TriggerObjectCollection& collection(filterProducts_[iT]);

    bool matched(false);
    for(unsigned iO(0); iO < collection.size(); iO++){
      matched |= reco::deltaR(collection[iO], _cand) < _dR;
      if(matched) break;
    }

    if(!matched) return false;
  }

  return true;
}
