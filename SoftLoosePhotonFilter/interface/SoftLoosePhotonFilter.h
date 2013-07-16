#ifndef SoftLoosePhotonFilter_h
#define SoftLoosePhotonFilter_h

// -*- C++ -*-
//
// Package:    SoftLoosePhotonFilter
// Class:      SoftLoosePhotonFilter
// 
/**\class SoftLoosePhotonFilter SoftLoosePhotonFilter.cc SusyAnalysis/SoftLoosePhotonFilter/src/SoftLoosePhotonFilter.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Yutaro Iiyama,512 1-005,+41227670489,
//         Created:  Fri Nov  2 15:34:30 CET 2012
// $Id: SoftLoosePhotonFilter.h,v 1.3 2013/06/07 13:48:36 yiiyama Exp $
//
//

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDFilter.h"
#include "FWCore/Framework/interface/Event.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "FWCore/Utilities/interface/InputTag.h"

//
// class declaration
//

class SoftLoosePhotonFilter : public edm::EDFilter {
   public:
      explicit SoftLoosePhotonFilter(const edm::ParameterSet&);
      ~SoftLoosePhotonFilter();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

   private:
      virtual void beginJob() ;
      virtual bool filter(edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;
      
      virtual bool beginRun(edm::Run&, edm::EventSetup const&);
      virtual bool endRun(edm::Run&, edm::EventSetup const&);
      virtual bool beginLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&);
      virtual bool endLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&);

      // ----------member data ---------------------------
  enum CutType {
    HOverELoose,
    ShowerShapeLoose,
    IsoLoose,
    ElectronVeto,
    nCutType
  };

  edm::InputTag photonCollectionTag_;
  edm::InputTag electronCollectionTag_;
  edm::InputTag conversionCollectionTag_;
  edm::InputTag beamSpotTag_;
//   edm::InputTag pfCandidateCollectionTag_;
//   edm::InputTag vertexCollectionTag_;
  //  PFIsolationEstimator isolator_;
  double ptThreshold_;
  std::vector<bool> doCut_;
};

#endif
