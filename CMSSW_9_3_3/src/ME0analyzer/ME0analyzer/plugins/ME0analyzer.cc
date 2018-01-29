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

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
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

#include <TTree.h>
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

      //variables for efficiency calculation
      float totmu_pos = 0;	float totmu_neg = 0;	//nr. of GenParticle muons in positive/negative endcap
      float meas_any_pos = 0;  	float meas_any_neg = 0; //measured, any nr. of layers fired accepted
      float meas_1L_pos = 0; 	float meas_1L_neg = 0;  //measured, 1 layer fired accepted
      float meas_2L_pos = 0;  	float meas_2L_neg = 0;  //measured, 2 layers fired accepted
      float meas_3L_pos = 0;  	float meas_3L_neg = 0;  //measured, 3 layers fired accepted
      float meas_4L_pos = 0;  	float meas_4L_neg = 0;  //measured, 4 layers fired accepted
      float meas_5L_pos = 0;  	float meas_5L_neg = 0;  //measured, 5 layers fired accepted
      float meas_6L_pos = 0;  	float meas_6L_neg = 0;  //measured, 6 layers fired accepted

      float eff_any_pos = 0;  float eff_any_neg = 0;  float eff_any_posneg = 0;	
      float eff_1L_pos  = 0;  float eff_1L_neg  = 0;  float eff_1L_posneg  = 0;
      float eff_2L_pos  = 0;  float eff_2L_neg  = 0;  float eff_2L_posneg  = 0;
      float eff_3L_pos  = 0;  float eff_3L_neg  = 0;  float eff_3L_posneg  = 0;
      float eff_4L_pos  = 0;  float eff_4L_neg  = 0;  float eff_4L_posneg  = 0;
      float eff_5L_pos  = 0;  float eff_5L_neg  = 0;  float eff_5L_posneg  = 0;
      float eff_6L_pos  = 0;  float eff_6L_neg  = 0;  float eff_6L_posneg  = 0;

      float erreff_any_pos  = 0;   float erreff_any_neg  = 0;	float erreff_any_posneg  = 0;
      float erreff_1L_pos  =  0;   float erreff_1L_neg  =  0;	float erreff_1L_posneg  =  0;
      float erreff_2L_pos  =  0;   float erreff_2L_neg  =  0;	float erreff_2L_posneg  =  0;
      float erreff_3L_pos  =  0;   float erreff_3L_neg  =  0;	float erreff_3L_posneg  =  0;
      float erreff_4L_pos  =  0;   float erreff_4L_neg  =  0;	float erreff_4L_posneg  =  0;
      float erreff_5L_pos  =  0;   float erreff_5L_neg  =  0;	float erreff_5L_posneg  =  0;
      float erreff_6L_pos  =  0;   float erreff_6L_neg  =  0;	float erreff_6L_posneg  =  0;

      float nL_neg, nL_pos;

      bool v = 1 ; //verbose





      // ----------member data ---------------------------
      TTree *tr;
      TTree *tr2;
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
   

   v = iConfig.getParameter<bool>("verbose");

   Service<TFileService> fs;
//   h_n3mu = fs->make<TH1F>("n3mu", "", 10, 0, 10);
   tr2 = fs->make<TTree>("Run", "");
   //efficiencies
   tr2->Branch("eff_any_pos",	&eff_any_pos,	"eff_any_pos");	tr2->Branch("eff_any_neg",	&eff_any_neg,	"eff_any_neg");
   tr2->Branch("eff_1L_pos",	&eff_1L_pos,	"eff_1L_pos");	tr2->Branch("eff_1L_neg",	&eff_1L_neg,	"eff_1L_neg");
   tr2->Branch("eff_2L_pos",	&eff_2L_pos,	"eff_2L_pos");	tr2->Branch("eff_2L_neg",	&eff_2L_neg,	"eff_2L_neg");
   tr2->Branch("eff_3L_pos",	&eff_3L_pos,	"eff_3L_pos");	tr2->Branch("eff_3L_neg",	&eff_3L_neg,	"eff_3L_neg");
   tr2->Branch("eff_4L_pos",	&eff_4L_pos,	"eff_4L_pos");	tr2->Branch("eff_4L_neg",	&eff_4L_neg,	"eff_4L_neg");
   tr2->Branch("eff_5L_pos",	&eff_5L_pos,	"eff_5L_pos");	tr2->Branch("eff_5L_neg",	&eff_5L_neg,	"eff_5L_neg");
   tr2->Branch("eff_6L_pos",	&eff_6L_pos,	"eff_6L_pos");	tr2->Branch("eff_6L_neg",	&eff_6L_neg,	"eff_6L_neg");
   //efficiencies errors
   tr2->Branch("erreff_any_pos",&erreff_any_pos,"erreff_any_pos");	tr2->Branch("erreff_any_neg",&erreff_any_neg,"erreff_any_neg");
   tr2->Branch("erreff_1L_pos",	&erreff_1L_pos,	"erreff_1L_pos");	tr2->Branch("erreff_1L_neg",	&erreff_1L_neg,	"erreff_1L_neg");
   tr2->Branch("erreff_2L_pos",	&erreff_2L_pos,	"erreff_2L_pos");	tr2->Branch("erreff_2L_neg",	&erreff_2L_neg,	"erreff_2L_neg");
   tr2->Branch("erreff_3L_pos",	&erreff_3L_pos,	"erreff_3L_pos");	tr2->Branch("erreff_3L_neg",	&erreff_3L_neg,	"erreff_3L_neg");
   tr2->Branch("erreff_4L_pos",	&erreff_4L_pos,	"erreff_4L_pos");	tr2->Branch("erreff_4L_neg",	&erreff_4L_neg,	"erreff_4L_neg");
   tr2->Branch("erreff_5L_pos",	&erreff_5L_pos,	"erreff_5L_pos");	tr2->Branch("erreff_5L_neg",	&erreff_5L_neg,	"erreff_5L_neg");
   tr2->Branch("erreff_6L_pos",	&erreff_6L_pos,	"erreff_6L_pos");	tr2->Branch("erreff_6L_neg",	&erreff_6L_neg,	"erreff_6L_neg");

   //measured hits on >= nr. layers
   tr2->Branch("totmu_pos",	&totmu_pos,	"totmu_pos");
   tr2->Branch("totmu_neg",	&totmu_neg,	"totmu_neg");
   tr2->Branch("meas_1L_pos",	&meas_1L_pos,	"meas_1L_pos");
   tr2->Branch("meas_2L_pos",	&meas_2L_pos,	"meas_2L_pos");
   tr2->Branch("meas_3L_pos",	&meas_3L_pos,	"meas_3L_pos");
   tr2->Branch("meas_4L_pos",	&meas_4L_pos,	"meas_4L_pos");
   tr2->Branch("meas_5L_pos",	&meas_5L_pos,	"meas_5L_pos");
   tr2->Branch("meas_6L_pos",	&meas_6L_pos,	"meas_6L_pos");
   tr2->Branch("meas_1L_neg",	&meas_1L_neg,	"meas_1L_neg");
   tr2->Branch("meas_2L_neg",	&meas_2L_neg,	"meas_2L_neg");
   tr2->Branch("meas_3L_neg",	&meas_3L_neg,	"meas_3L_neg");
   tr2->Branch("meas_4L_neg",	&meas_4L_neg,	"meas_4L_neg");
   tr2->Branch("meas_5L_neg",	&meas_5L_neg,	"meas_5L_neg");
   tr2->Branch("meas_6L_neg",	&meas_6L_neg,	"meas_6L_neg");

   tr = fs->make<TTree>("Event", "");
   tr->Branch("nL_neg"     ,	&nL_neg     , 	"nL_neg"     );
   tr->Branch("nL_pos"     ,	&nL_pos     , 	"nL_pos"     );
}

ME0analyzer::~ME0analyzer()
{
 
   tr2->Fill();
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

   vector<ME0DetId> 	me0detid_v1;
   vector<ME0DetId> 	me0detid_v2;
   vector<ME0DetId> 	me0detid_v3;
   vector<ME0DetId> 	me0detid_pos; 
   vector<ME0DetId> 	me0detid_neg;
   vector<unsigned int> posGenPart;
   vector<unsigned int> negGenPart;
   vector<int>		layers_pos;
   vector<int>		layers_neg;



//------------------ GEN PARTICLES ------------------------------
     Handle<GenParticleCollection> genParticles;
     iEvent.getByToken(genToken_, genParticles);

     for(size_t i = 0; i < genParticles->size(); ++ i) 
     {
       if (v) cout << "GenParticle index: " << i << endl;
       const GenParticle & p = (*genParticles)[i];
       //const Candidate * mom =( p.mother());
       if (v) 
         {
         cout  <<"id: "<<p.pdgId()<<"\tstatus: "<<p.status()
         <<"\npt: "<<p.pt()<<"\teta: "<<p.eta()
         <<"\npx: "<<p.px()<<"\tpy: "<<p.py()<<"\tpz: "<<p.pz()
         <<"\nphi: "<<p.phi()<<"\tmass: "<<p.mass()
         <<"\nvx: "<<p.vx()<<"\tvy "<<p.vy()
         <<"\tvz: "<<p.vz() << "\n" << endl;
         }
      // << "\nmother: "<<p.mother()<<"\tMOTHER id: "<<mom->pdgId()
      // <<"\tMOTHER pt:"<<mom->pt()<<endl<<endl;

      if ( fabs(p.eta()) < 2.0 || fabs(p.eta()) > 2.8 ) return;

      if ( p.pz()>0 ) { posGenPart.push_back(i); totmu_pos++; }
      else            { negGenPart.push_back(i); totmu_neg++; }
     }


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
    cout << "region " << me0id.region() << "layer " << me0id.layer() << endl;
    if ( me0id.region() > 0 )  	{ me0detid_pos.push_back( me0id ) ;
    			  	  layers_pos.push_back( me0id.layer() ) ;
			       	}
    else			{ me0detid_neg.push_back( me0id ) ;
    				  layers_neg.push_back( me0id.layer() ) ;
			       	}
   }


//-------------------- EFFICIENCY ------------------------------------
unsigned int nL = layers_pos.size();
nL_pos = nL;
if (v) cout << "nL_pos = " << nL_pos << endl;
if ( nL_pos>0 )	
    {
    meas_any_pos++;
    if (nL >= 1) 	meas_1L_pos++;
    if (nL >= 2) 	meas_2L_pos++;
    if (nL >= 3) 	meas_3L_pos++;
    if (nL >= 4) 	meas_4L_pos++;
    if (nL >= 5) 	meas_5L_pos++;
    if (nL >= 6) 	meas_6L_pos++;
    }

nL = layers_neg.size();
nL_neg = nL;
if (v) cout << "nL_neg = " << nL_neg << endl;
if ( nL_neg>0 )	
    {
    meas_any_neg++;
    if (nL >= 1) 	meas_1L_neg++;
    if (nL >= 2) 	meas_2L_neg++;
    if (nL >= 3) 	meas_3L_neg++;
    if (nL >= 4) 	meas_4L_neg++;
    if (nL >= 5) 	meas_5L_neg++;
    if (nL >= 6) 	meas_6L_neg++;
    }

if (v) {

  cout << "meas_1L_pos " << meas_1L_pos << endl;
  cout << "meas_2L_pos " << meas_2L_pos << endl;
  cout << "meas_3L_pos " << meas_3L_pos << endl;
  cout << "meas_4L_pos " << meas_4L_pos << endl;
  cout << "meas_5L_pos " << meas_5L_pos << endl;
  cout << "meas_6L_pos " << meas_6L_pos << endl;
  cout << "meas_1L_neg " << meas_1L_neg << endl;
  cout << "meas_2L_neg " << meas_2L_neg << endl;
  cout << "meas_3L_neg " << meas_3L_neg << endl;
  cout << "meas_4L_neg " << meas_4L_neg << endl;
  cout << "meas_5L_neg " << meas_5L_neg << endl;
  cout << "meas_6L_neg " << meas_6L_neg << endl;

}

//efficiencies
eff_any_pos = meas_any_pos/totmu_pos;		eff_any_neg = meas_any_neg/totmu_neg;
eff_1L_pos  = meas_1L_pos/totmu_pos;		eff_1L_neg  = meas_1L_neg/totmu_neg;
eff_2L_pos  = meas_2L_pos/totmu_pos;		eff_2L_neg  = meas_2L_neg/totmu_neg;
eff_3L_pos  = meas_3L_pos/totmu_pos;		eff_3L_neg  = meas_3L_neg/totmu_neg;
eff_4L_pos  = meas_4L_pos/totmu_pos;		eff_4L_neg  = meas_4L_neg/totmu_neg;
eff_5L_pos  = meas_5L_pos/totmu_pos;		eff_5L_neg  = meas_5L_neg/totmu_neg;
eff_6L_pos  = meas_6L_pos/totmu_pos;		eff_6L_neg  = meas_6L_neg/totmu_neg;

eff_any_posneg = ( meas_any_pos + meas_any_neg	) / ( totmu_pos	+ totmu_neg ) ;	
eff_1L_posneg  = ( meas_1L_pos  + meas_1L_neg	) / ( totmu_pos	+ totmu_neg ) ;	
eff_2L_posneg  = ( meas_2L_pos  + meas_2L_neg   ) / ( totmu_pos	+ totmu_neg ) ;	
eff_3L_posneg  = ( meas_3L_pos  + meas_3L_neg  	) / ( totmu_pos	+ totmu_neg ) ;	
eff_4L_posneg  = ( meas_4L_pos  + meas_4L_neg  	) / ( totmu_pos	+ totmu_neg ) ;	
eff_5L_posneg  = ( meas_5L_pos  + meas_5L_neg  	) / ( totmu_pos	+ totmu_neg ) ;	
eff_6L_posneg  = ( meas_6L_pos  + meas_6L_neg  	) / ( totmu_pos	+ totmu_neg ) ;	

//efficiencies errors
erreff_any_pos  = sqrt ( eff_any_pos*(1-eff_any_pos)/totmu_pos ) ;	erreff_any_neg  = sqrt ( eff_any_neg*(1-eff_any_neg)/totmu_neg ) ;
erreff_1L_pos  = sqrt ( eff_1L_pos*(1-eff_1L_pos)/totmu_pos ) ;		erreff_1L_neg  = sqrt ( eff_1L_neg*(1-eff_1L_neg)/totmu_neg ) ;
erreff_2L_pos  = sqrt ( eff_2L_pos*(1-eff_2L_pos)/totmu_pos ) ;		erreff_2L_neg  = sqrt ( eff_2L_neg*(1-eff_2L_neg)/totmu_neg ) ;
erreff_3L_pos  = sqrt ( eff_3L_pos*(1-eff_3L_pos)/totmu_pos ) ;		erreff_3L_neg  = sqrt ( eff_3L_neg*(1-eff_3L_neg)/totmu_neg ) ;
erreff_4L_pos  = sqrt ( eff_4L_pos*(1-eff_4L_pos)/totmu_pos ) ;		erreff_4L_neg  = sqrt ( eff_4L_neg*(1-eff_4L_neg)/totmu_neg ) ;
erreff_5L_pos  = sqrt ( eff_5L_pos*(1-eff_5L_pos)/totmu_pos ) ;		erreff_5L_neg  = sqrt ( eff_5L_neg*(1-eff_5L_neg)/totmu_neg ) ;
erreff_6L_pos  = sqrt ( eff_6L_pos*(1-eff_6L_pos)/totmu_pos ) ;		erreff_6L_neg  = sqrt ( eff_6L_neg*(1-eff_6L_neg)/totmu_neg ) ;

erreff_any_posneg  = sqrt ( eff_any_posneg*(1-eff_any_posneg) /	( totmu_pos	+ totmu_neg ) ) ;	
erreff_1L_posneg  = sqrt ( eff_1L_posneg*(1-eff_1L_posneg) /	( totmu_pos	+ totmu_neg ) ) ;		
erreff_2L_posneg  = sqrt ( eff_2L_posneg*(1-eff_2L_posneg) /	( totmu_pos	+ totmu_neg ) ) ;		
erreff_3L_posneg  = sqrt ( eff_3L_posneg*(1-eff_3L_posneg) /	( totmu_pos	+ totmu_neg ) ) ;		
erreff_4L_posneg  = sqrt ( eff_4L_posneg*(1-eff_4L_posneg) /	( totmu_pos	+ totmu_neg ) ) ;		
erreff_5L_posneg  = sqrt ( eff_5L_posneg*(1-eff_5L_posneg) /	( totmu_pos	+ totmu_neg ) ) ;		
erreff_6L_posneg  = sqrt ( eff_6L_posneg*(1-eff_6L_posneg) /	( totmu_pos	+ totmu_neg ) ) ;		

if (v) {

cout << "eff_any_pos "	<< eff_any_pos	<< ",\teff_any_neg " <<	eff_any_neg <<	",\teff_any_posneg " <<	eff_any_posneg << endl ;
cout << "eff_1L_pos  "	<< eff_1L_pos	<< ",\teff_1L_neg  " <<	eff_1L_neg  <<	",\teff_1L_posneg  " <<	eff_1L_posneg  << endl ;
cout << "eff_2L_pos  "	<< eff_2L_pos	<< ",\teff_2L_neg  " <<	eff_2L_neg  <<	",\teff_2L_posneg  " <<	eff_2L_posneg  << endl ;
cout << "eff_3L_pos  "	<< eff_3L_pos	<< ",\teff_3L_neg  " <<	eff_3L_neg  <<	",\teff_3L_posneg  " <<	eff_3L_posneg  << endl ;
cout << "eff_4L_pos  "	<< eff_4L_pos	<< ",\teff_4L_neg  " <<	eff_4L_neg  <<	",\teff_4L_posneg  " <<	eff_4L_posneg  << endl ;
cout << "eff_5L_pos  "	<< eff_5L_pos	<< ",\teff_5L_neg  " <<	eff_5L_neg  <<	",\teff_5L_posneg  " <<	eff_5L_posneg  << endl ;
cout << "eff_6L_pos  "	<< eff_6L_pos	<< ",\teff_6L_neg  " <<	eff_6L_neg  <<	",\teff_6L_posneg  " <<	eff_6L_posneg  << endl ;





}

//------------------- TEST ----------------------

//pair<MuonDigiCollection<ME0DetId,ME0Digi>::const_iterator,MuonDigiCollection<ME0DetId,ME0Digi>::const_iterator> digiByDet_it;
//digiByDet_it = MuonDigiCollection<ME0DetId,ME0Digi>::get(me0detid_v1[0]);

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
      //int bx    = (*itr).bx();
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
    


   tr->Fill();
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
