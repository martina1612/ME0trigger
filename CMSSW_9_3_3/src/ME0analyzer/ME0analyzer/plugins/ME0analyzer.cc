// -*- C++ -*-
//
// Package:    ME0analyzer/ME0analyzer
// Class:      ME0analyzer
// 
/**\class ME0analyzer ME0analyzer.cc ME0analyzer/ME0analyzer/plugins/ME0analyzer.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Martina Ressegotti
//         Created:  Wed, 24 Jan 2018 13:41:06 GMT
//
//


// system include files
#include <memory>

#include <iterator>
#include <stdlib.h>
#include <stdio.h>
#include <vector>
#include <map>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "FWCore/Utilities/interface/EDGetToken.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "DataFormats/GEMDigi/interface/ME0TriggerDigi.h"
#include "Geometry/GEMGeometry/interface/ME0Geometry.h"
#include "DataFormats/MuonData/interface/MuonDigiCollection.h"
#include "DataFormats/GEMDigi/interface/ME0Digi.h"
#include "Geometry/Records/interface/MuonGeometryRecord.h"
#include "Geometry/GEMGeometry/interface/ME0EtaPartitionSpecs.h"
#include "Geometry/CommonTopologies/interface/StripTopology.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/GEMDigi/interface/ME0PadDigi.h"
#include "DataFormats/GEMDigi/interface/ME0PadDigiCluster.h"
#include "DataFormats/GEMDigi/interface/ME0PadDigiClusterCollection.h"
//
// class declaration
//

// If the analyzer does not use TFileService, please remove
// the template argument to the base class so the class inherits
// from  edm::one::EDAnalyzer<> and also remove the line from
// constructor "usesResource("TFileService");"
// This will improve performance in multithreaded jobs.

using namespace edm;
using namespace std;
using namespace reco;

class ME0analyzer : public edm::one::EDAnalyzer<edm::one::SharedResources>  {
   public:
      explicit ME0analyzer(const edm::ParameterSet&);
      ~ME0analyzer();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


   private:
      virtual void beginJob() override;
      virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
      virtual void endJob() override;
      edm::EDGetTokenT<MuonDigiCollection<ME0DetId,ME0Digi>> ME0DigiToken_;
      edm::EDGetTokenT<MuonDigiCollection<ME0DetId,ME0PadDigi>> ME0PadDigiToken_;
      edm::EDGetTokenT<MuonDigiCollection<ME0DetId,ME0PadDigiCluster>> ME0PadDigiClusterToken_;
      edm::EDGetTokenT<GenParticleCollection> genToken_;

      // ----------member data ---------------------------
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
ME0analyzer::ME0analyzer(const edm::ParameterSet& iConfig):
    ME0DigiToken_(consumes<MuonDigiCollection<ME0DetId,ME0Digi>>(iConfig.getParameter<edm::InputTag>("me0DigiToken"))),
    ME0PadDigiToken_(consumes<MuonDigiCollection<ME0DetId,ME0PadDigi>>(iConfig.getParameter<edm::InputTag>("me0PadDigiToken"))),
    ME0PadDigiClusterToken_(consumes<MuonDigiCollection<ME0DetId,ME0PadDigiCluster>>(iConfig.getParameter<edm::InputTag>("me0PadDigiClusterToken"))),
    genToken_(consumes<GenParticleCollection>(iConfig.getParameter<InputTag>("genParticles")))

{
   //now do what ever initialization is needed
   usesResource("TFileService");
}


ME0analyzer::~ME0analyzer()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
ME0analyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;
 
   Handle<MuonDigiCollection<ME0DetId,ME0Digi>> me0digiH;
   Handle<MuonDigiCollection<ME0DetId,ME0PadDigi>> me0padDigiH;
   Handle<MuonDigiCollection<ME0DetId,ME0PadDigiCluster>> me0padDigiClusterH;
   iEvent.getByToken(ME0DigiToken_, me0digiH);
   iEvent.getByToken(ME0PadDigiToken_, me0padDigiH);
   iEvent.getByToken(ME0PadDigiClusterToken_, me0padDigiClusterH);

   vector<ME0DetId> me0detid_v1;
   vector<ME0DetId> me0detid_v2;
   vector<ME0DetId> me0detid_v3;

//------------------- DIGI --------------------------------------  
   for ( DigiContainerIterator<ME0DetId,ME0Digi> it = me0digiH->begin(); it != me0digiH->end(); ++it )
   {
    ME0DetId me0id = (*it).first;
    cout << me0id << endl;
    MuonDigiCollection<ME0DetId,ME0Digi>::const_iterator itr;
    for ( itr =((*it).second).first; itr!= ((*it).second).second; itr++ )
     {
      int strip = (*itr).strip();
      int bx    = (*itr).bx();
      cout << "strip: " << strip << "  bx: " << bx << endl;
     }
    me0detid_v1.push_back(me0id);
   }

//------------------- PAD DIGI --------------------------------------  
   for ( DigiContainerIterator<ME0DetId,ME0PadDigi> it = me0padDigiH->begin(); it != me0padDigiH->end(); ++it )
   {
    ME0DetId me0id = (*it).first;
    cout << me0id << endl;
    MuonDigiCollection<ME0DetId,ME0PadDigi>::const_iterator itr;
    for ( itr =((*it).second).first; itr!= ((*it).second).second; itr++ )
     {
      int pad = (*itr).pad();
      int bx    = (*itr).bx();
      cout << "pad: " << pad << "  bx: " << bx << endl;
     }
    me0detid_v2.push_back(me0id);
   }

//------------------- PAD DIGI CLUSTER ------------------------------  
   for ( DigiContainerIterator<ME0DetId,ME0PadDigiCluster> it = me0padDigiClusterH->begin(); it != me0padDigiClusterH->end(); ++it )
   {
    ME0DetId me0id = (*it).first;
    cout << me0id << endl;
    MuonDigiCollection<ME0DetId,ME0PadDigiCluster>::const_iterator itr;
    for ( itr =((*it).second).first; itr!= ((*it).second).second; itr++ )
     {
      vector<uint16_t> pads = (*itr).pads();
      int bx    = (*itr).bx();
      cout << "pads: " ;
      for ( unsigned int i=0; i<pads.size(); i++)  cout << pads[i] << " " ;
     }
    me0detid_v3.push_back(me0id);
   }

//Verify that the ME0DetId list is the same in the three cases above
if( equal(me0detid_v1.begin(), me0detid_v1.end(), me0detid_v2.begin()) )
  if( equal(me0detid_v1.begin(), me0detid_v1.end(), me0detid_v3.begin()) )
    if ( me0detid_v1.size() == me0detid_v2.size() && me0detid_v1.size() == me0detid_v3.size() )
      cout << "The three classes have the same list of ME0DetId" << endl;
    

//------------------ GEN PARTICLES ------------------------------
     Handle<GenParticleCollection> genParticles;
     iEvent.getByToken(genToken_, genParticles);

     for(size_t i = 0; i < genParticles->size(); ++ i) 
     {
       cout << "GenParticle index: " << i << endl;
       const GenParticle & p = (*genParticles)[i];
       //const Candidate * mom =( p.mother());
       cout  <<"id: "<<p.pdgId()<<"\tstatus: "<<p.status()
       <<"\npt: "<<p.pt()<<"\teta: "<<p.eta()
       <<"\nphi: "<<p.phi()<<"\tmass: "<<p.mass()
       <<"\nvx: "<<p.vx()<<"\tvy "<<p.vy()
       <<"\tvz: "<<p.vz() << "\n" << endl;
      // << "\nmother: "<<p.mother()<<"\tMOTHER id: "<<mom->pdgId()
      // <<"\tMOTHER pt:"<<mom->pt()<<endl<<endl;
     }

}


// ------------ method called once each job just before starting event loop  ------------
void 
ME0analyzer::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
ME0analyzer::endJob() 
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
ME0analyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(ME0analyzer);
