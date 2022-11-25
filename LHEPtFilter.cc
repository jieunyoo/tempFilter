
//
// Package:    LHEPtFilter
// Class:      LHEPtFilter
//
/* 

 Description: Filter to select events with pT in a given range.
 (Based on LHEGenericFilter)

     
*/
//

// system include files
#include <memory>
#include <iostream>
#include <set>

// user include files
#include "Math/Vector4D.h"
#include "Math/Vector4Dfwd.h"

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/global/EDFilter.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "SimDataFormats/GeneratorProducts/interface/LHEEventProduct.h"

//
// class declaration
//

class LHEPtFilter : public edm::global::EDFilter<> {
public:
  explicit LHEPtFilter(const edm::ParameterSet&);
  ~LHEPtFilter() override;

  bool filter(edm::StreamID, edm::Event&, const edm::EventSetup&) const override;

//private: //change to protected so derived classes can access this
protected:
  // ----------member data ---------------------------

  edm::EDGetTokenT<LHEEventProduct> src_;
  std::vector<int> pdgIdVec_;
  std::set<int> pdgIds_;  // Set of PDG Ids to include
  double ptMin_;          // number of particles required to pass filter
  double ptMax_;          // number of particles required to pass filter

  //bool scanIntermediateRes_; ///ADDED THIS NEW PARAMETER
};

using namespace edm;
using namespace std;

LHEPtFilter::LHEPtFilter(const edm::ParameterSet& iConfig)
    : pdgIdVec_(iConfig.getParameter<std::vector<int>>("selectedPdgIds")),
      ptMin_(iConfig.getParameter<double>("ptMin")),
      ptMax_(iConfig.getParameter<double>("ptMax")) {
  //here do whatever other initialization is needed
  src_ = consumes<LHEEventProduct>(iConfig.getParameter<edm::InputTag>("src"));
  pdgIds_ = std::set<int>(pdgIdVec_.begin(), pdgIdVec_.end());
}

LHEPtFilter::~LHEPtFilter() {
  // do anything here that needs to be done at destruction time
  // (e.g. close files, deallocate resources etc.)
}

// ------------ method called to skim the data  ------------
bool LHEPtFilter::filter(edm::StreamID, edm::Event& iEvent, const edm::EventSetup& iSetup) const {
  edm::Handle<LHEEventProduct> EvtHandle;
  iEvent.getByToken(src_, EvtHandle);

  std::vector<lhef::HEPEUP::FiveVector> lheParticles = EvtHandle->hepeup().PUP;
  std::vector<ROOT::Math::PxPyPzEVector> cands;

  for (unsigned int i = 0; i < lheParticles.size(); ++i) {
    //std::cout << "lheParticles.size()" << lheParticles.size() << std::endl;
    //std::cout << "EvtHandle->hepeup().ISTUP[i]" << EvtHandle->hepeup().ISTUP[i] << std::endl;

    if (EvtHandle->hepeup().ISTUP[i] != 2) {  // keep only intermediate resonances
    //if (EvtHandle->hepeup().ISTUP[i] != 1) {  // keep only outgoing particles
      continue;
    }
    int pdgId = EvtHandle->hepeup().IDUP[i];
    //cout << "pdgId " << pdgId << endl;
    if (pdgIds_.count(pdgId)) {
      cands.push_back(
          ROOT::Math::PxPyPzEVector(lheParticles[i][0], lheParticles[i][1], lheParticles[i][2], lheParticles[i][3]));
    }
  }
  double vpt_ = -1;
  if (!cands.empty()) {
    ROOT::Math::PxPyPzEVector tot = cands.at(0);
    for (unsigned icand = 1; icand < cands.size(); ++icand) {
      tot += cands.at(icand);
    }
    vpt_ = tot.pt();
  }
  if ((ptMax_ < 0. || vpt_ <= ptMax_) && vpt_ > ptMin_) {
    return true;
  } else {
    return false;
  }
};


//***********************************changed this***********************************
//#derived class - i want this new class to derive from the exisint LHEPtFilter one and just create a new function that changes
//one line in the below filter function to filter for intermediate particles
class LHEIntermediateResonancePtFilter : public LHEPtFilter {

public:
  bool filterIntermediateParticle(edm::StreamID, edm::Event&, const edm::EventSetup&) const;
//bool LHEIntermediateResonancePtFilter::filterIntermediateParticle(edm::StreamID, edm::Event& iEvent, const edm::EventSetup& iSetup) const {

bool filterIntermediateParticle(edm::StreamID, edm::Event& iEvent, const edm::EventSetup& iSetup)  {
  edm::Handle<LHEEventProduct> EvtHandle;
  iEvent.getByToken(src_, EvtHandle);

  std::vector<lhef::HEPEUP::FiveVector> lheParticles = EvtHandle->hepeup().PUP;
  std::vector<ROOT::Math::PxPyPzEVector> cands;

  for (unsigned int i = 0; i < lheParticles.size(); ++i) {
    //std::cout << "lheParticles.size()" << lheParticles.size() << std::endl;
    //std::cout << "EvtHandle->hepeup().ISTUP[i]" << EvtHandle->hepeup().ISTUP[i] << std::endl;

    if (EvtHandle->hepeup().ISTUP[i] != 2) {  // CHANGED THIS: want a cut only on intermediate resonances
      continue;
    }
    int pdgId = EvtHandle->hepeup().IDUP[i];
    //cout << "pdgId " << pdgId << endl;
    if (pdgIds_.count(pdgId)) {
      cands.push_back(
          ROOT::Math::PxPyPzEVector(lheParticles[i][0], lheParticles[i][1], lheParticles[i][2], lheParticles[i][3]));
    }
  }
  double vpt_ = -1;
  if (!cands.empty()) {
    ROOT::Math::PxPyPzEVector tot = cands.at(0);
    for (unsigned icand = 1; icand < cands.size(); ++icand) {
      tot += cands.at(icand);
    }
    vpt_ = tot.pt();
  }
  if ((ptMax_ < 0. || vpt_ <= ptMax_) && vpt_ > ptMin_) {
    return true;
  } else {
    return false;
  }
}



}



//define this as a plug-in
DEFINE_FWK_MODULE(LHEPtFilter);

//try adding this module as well
//DEFINE_FWK_MODULE(LHEIntermediateResonancePtFilter);
