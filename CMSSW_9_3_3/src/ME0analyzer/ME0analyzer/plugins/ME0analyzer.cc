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
#include <bitset>
#include <algorithm>

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
#include <TH1F.h>
#include <TGraphErrors.h>
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

      bool v = 1 ; //verbose initialization

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

      Float_t 	eff_x_pos[6] = { 1,2,3,4,5,6 }; 	Float_t eff_y_pos[6]; 	Float_t erreff_x_pos[6] = {}; 	 Float_t erreff_y_pos[6];
      Float_t 	eff_x_neg[6] = { 1,2,3,4,5,6 }; 	Float_t eff_y_neg[6]; 	Float_t erreff_x_neg[6] = {}; 	 Float_t erreff_y_neg[6];

      float bx_pos, bx_neg;
      float nL_pos, nL_neg;

      short int		patternL_pos		= 0;
      short int		patternL_neg		= 0;
      short int		patternD_pos		= 0;
      short int		patternD_neg		= 0;
      float		patternLtree_pos	= 0;
      float		patternLtree_neg	= 0;
      float		patternDtree_pos	= 0;
      float		patternDtree_neg	= 0;
      vector<float>	patternDtree_vec_pos;
      vector<float>	patternDtree_vec_neg;

      //histograms GenParticles
      TH1F * h_p_pos 	;
      TH1F * h_p_neg 	;
      TH1F * h_pt_pos 	;
      TH1F * h_pt_neg 	;
      TH1F * h_px_pos 	;
      TH1F * h_px_neg 	;
      TH1F * h_py_pos 	;
      TH1F * h_py_neg 	;
      TH1F * h_pz_pos 	;
      TH1F * h_pz_neg 	;
      TH1F * h_eta_pos 	;
      TH1F * h_eta_neg 	;
      TH1F * h_vx_pos 	;
      TH1F * h_vx_neg 	;
      TH1F * h_vy_pos 	;
      TH1F * h_vy_neg 	;
      TH1F * h_vz_pos 	;
      TH1F * h_vz_neg 	;
      TH1F * h_id_pos 	;
      TH1F * h_id_neg 	;

      //histograms ME0
      TH1F * h_nStr_pos		;
      TH1F * h_nStr_neg		;
      TH1F * h_holeSize_pos	;
      TH1F * h_holeSize_neg	;
      TH1F * h_digiPattern_pos	;
      TH1F * h_digiPattern_neg	;

      //TGraph efficiency
      TGraphErrors * g_effVsnL_pos;
      TGraphErrors * g_effVsnL_neg;

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
   Service<TFileService> fs;

   v = iConfig.getParameter<bool>("verbose");

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
   tr->Branch("bx_pos"     ,	&bx_pos     , 	"bx_pos"     );
   tr->Branch("bx_neg"     ,	&bx_neg     , 	"bx_neg"     );
   tr->Branch("nL_pos"     ,	&nL_pos     , 	"nL_pos"     );
   tr->Branch("nL_neg"     ,	&nL_neg     , 	"nL_neg"     );
   tr->Branch("patternLtree_pos" ,	&patternLtree_pos , "patternL_pos" );
   tr->Branch("patternLtree_neg" ,	&patternLtree_neg , "patternL_neg" );
//   tr->Branch("patternDtree_vec_pos" ,	&patternDtree_vec_pos , "patternD_vec_pos" );
//   tr->Branch("patternDtree_vec_neg" ,	&patternDtree_vec_neg , "patternD_vec_neg" );
  

   //histograms GenParticles
   h_p_pos 	= fs->make<TH1F>("h_p_pos","h_p_pos", 10000, 0, 1000);
   h_p_neg 	= fs->make<TH1F>("h_p_neg","h_p_neg", 10000, 0, 1000);
   h_pt_pos 	= fs->make<TH1F>("h_pt_pos","h_pt_pos", 10000, 0, 1000);
   h_pt_neg 	= fs->make<TH1F>("h_pt_neg","h_pt_neg", 10000, 0, 1000);
   h_px_pos 	= fs->make<TH1F>("h_px_pos","h_px_pos", 10000, 0, 1000);
   h_px_neg 	= fs->make<TH1F>("h_px_neg","h_px_neg", 10000, 0, 1000);
   h_py_pos 	= fs->make<TH1F>("h_py_pos","h_py_pos", 10000, 0, 1000);
   h_py_neg 	= fs->make<TH1F>("h_py_neg","h_py_neg", 10000, 0, 1000);
   h_pz_pos 	= fs->make<TH1F>("h_pz_pos","h_pz_pos", 10000, 0, 1000);
   h_pz_neg 	= fs->make<TH1F>("h_pz_neg","h_pz_neg", 10000, 0, 1000);
   h_eta_pos 	= fs->make<TH1F>("h_eta_pos","h_eta_pos", 500, -5, +5);
   h_eta_neg 	= fs->make<TH1F>("h_eta_neg","h_eta_neg", 500, -5, +5);
   h_vx_pos 	= fs->make<TH1F>("h_vx_pos","h_vx_pos", 10000, -500, 500);
   h_vx_neg 	= fs->make<TH1F>("h_vx_neg","h_vx_neg", 10000, -500, 500);
   h_vy_pos 	= fs->make<TH1F>("h_vy_pos","h_vy_pos", 10000, -500, 500);
   h_vy_neg 	= fs->make<TH1F>("h_vy_neg","h_vy_neg", 10000, -500, 500);
   h_vz_pos 	= fs->make<TH1F>("h_vz_pos","h_vz_pos", 10000, -800, 800);
   h_vz_neg 	= fs->make<TH1F>("h_vz_neg","h_vz_neg", 10000, -800, 800);
   h_id_pos 	= fs->make<TH1F>("h_id_pos","h_id_pos", 2001, -1000, 1000);
   h_id_neg 	= fs->make<TH1F>("h_id_neg","h_id_neg", 2001, -1000, 1000);

   //histograms ME0
   h_nStr_pos 	= fs->make<TH1F>("h_nStr_pos","h_nStr_pos", 301, 0, 300);
   h_nStr_neg 	= fs->make<TH1F>("h_nStr_neg","h_nStr_neg", 301, 0, 300);
   h_holeSize_pos = fs->make<TH1F>("h_holeSize_pos","h_holeSize_pos", 301, 0, 300);
   h_holeSize_neg = fs->make<TH1F>("h_holeSize_neg","h_holeSize_neg", 301, 0, 300);
   h_digiPattern_pos = fs->make<TH1F>("h_digiPattern_pos","h_digiPattern_pos", 301, 0, 300);
   h_digiPattern_neg = fs->make<TH1F>("h_digiPattern_neg","h_digiPattern_neg", 301, 0, 300);
   
   //TGraph efficiency
   //g_effVsnL_pos = fs->make<TGraphErrors>(6,eff_x_pos, eff_y_pos, erreff_x_pos, erreff_y_pos);
   //g_effVsnL_neg = fs->make<TGraphErrors>(6,eff_x_neg, eff_y_neg, erreff_x_neg, erreff_y_neg);

   //patterns
   //map<short int,int> firedLayers_pos;
   //map<short int,int> firedLayers_neg;

}

ME0analyzer::~ME0analyzer()
{
   //arrays for efficiency TGraphs
   eff_y_pos[0] = eff_1L_pos ; 	eff_y_pos[1] = eff_2L_pos ; 	eff_y_pos[2] = eff_3L_pos ; 
   eff_y_pos[3] = eff_4L_pos ; 	eff_y_pos[4] = eff_5L_pos ; 	eff_y_pos[5] = eff_6L_pos ;
   erreff_y_pos[0] = erreff_1L_pos ; erreff_y_pos[1] = erreff_2L_pos ; erreff_y_pos[2] = erreff_3L_pos ; 
   erreff_y_pos[3] = erreff_4L_pos ; erreff_y_pos[4] = erreff_5L_pos ; erreff_y_pos[5] = erreff_6L_pos ;
   
   eff_y_neg[0] = eff_1L_neg ; 	eff_y_neg[1] = eff_2L_neg ; 	eff_y_neg[2] = eff_3L_neg ; 
   eff_y_neg[3] = eff_4L_neg ; 	eff_y_neg[4] = eff_5L_neg ; 	eff_y_neg[5] = eff_6L_neg ;
   erreff_y_neg[0] = erreff_1L_neg ; erreff_y_neg[1] = erreff_2L_neg ; erreff_y_neg[2] = erreff_3L_neg ; 
   erreff_y_neg[3] = erreff_4L_neg ; erreff_y_neg[4] = erreff_5L_neg ; erreff_y_neg[5] = erreff_6L_neg ;

   if (v)
        {
   	cout << "Positive efficiency:" << endl;
	for ( unsigned int i = 0; i<6; i++)
	    {
	    cout << eff_y_pos[i] << " +- " << erreff_y_pos[i] << endl;
	    }
   	cout << "Negative efficiency:" << endl;
	for ( unsigned int i = 0; i<6; i++)
	    {
	    cout << eff_y_neg[i] << " +- " << erreff_y_neg[i] << endl;
	    }
	}
   
   //TGraph efficiency
   usesResource("TFileService");
   Service<TFileService> fs;
   g_effVsnL_pos = fs->make<TGraphErrors>(6,eff_x_pos, eff_y_pos, erreff_x_pos, erreff_y_pos);
   g_effVsnL_neg = fs->make<TGraphErrors>(6,eff_x_neg, eff_y_neg, erreff_x_neg, erreff_y_neg);

   g_effVsnL_pos->SetTitle("g_effVsnL_pos");
   g_effVsnL_neg->SetTitle("g_effVsnL_neg");
   g_effVsnL_pos->SetName("g_effVsnL_pos");
   g_effVsnL_neg->SetName("g_effVsnL_neg");

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

   vector<ME0DetId> 	me0detid_v1;		me0detid_v1.clear();
   vector<ME0DetId> 	me0detid_v2;            me0detid_v2.clear();
   vector<ME0DetId> 	me0detid_v3;            me0detid_v3.clear();
   vector<ME0DetId> 	me0detid_pos;           me0detid_pos.clear();
   vector<ME0DetId> 	me0detid_neg;           me0detid_neg.clear();
   vector<unsigned int> posGenPart;             posGenPart.clear();
   vector<unsigned int> negGenPart;             negGenPart.clear();
   vector<int>		layers_pos;             layers_pos.clear();
   vector<int>		layers_neg;             layers_neg.clear();

   patternDtree_vec_neg.clear();
   patternDtree_vec_pos.clear();

   map<ME0DetId,vector<int>> 			strips;		strips.clear();
   map<ME0DetId,vector<int>> 			pads;		pads.clear();
   map<ME0DetId,vector<vector<uint16_t>>> 	clusters;	clusters.clear();
   map<ME0DetId,short int> 			stripsPat;	stripsPat.clear();

   nL_pos	= 0;
   nL_neg	= 0;
   patternL_pos	= 0;
   patternL_neg	= 0;
   patternD_pos	= 0;
   patternD_neg	= 0;
   patternLtree_pos	= 0;
   patternLtree_neg	= 0;
   patternDtree_pos	= 0;
   patternDtree_neg	= 0;
   bx_pos		= 99;
   bx_neg		= 99;
   int nStr_pos	= 99;
   int nStr_neg	= 99;

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
      if ( fabs(p.pt()) < 0.5 ) return;

      if ( p.pz()>0 ) {
      		      posGenPart.push_back(i); totmu_pos++;
		      
		      h_p_pos 	->Fill(p.p());
                      h_pt_pos 	->Fill(p.pt());
                      h_px_pos 	->Fill(p.px());
                      h_py_pos 	->Fill(p.py());
                      h_pz_pos 	->Fill(p.pz());
                      h_eta_pos	->Fill(p.eta());
                      h_vx_pos 	->Fill(p.vx());
                      h_vy_pos 	->Fill(p.vy());
                      h_vz_pos 	->Fill(p.vz());
                      h_id_pos 	->Fill(p.pdgId());
		      }
      else            { 
      		      negGenPart.push_back(i); totmu_neg++; 
		      
		      h_p_neg 	->Fill(p.p());
                      h_pt_neg 	->Fill(p.pt());
                      h_px_neg 	->Fill(p.px());
                      h_py_neg 	->Fill(p.py());
                      h_pz_neg 	->Fill(p.pz());
                      h_eta_neg	->Fill(p.eta());
                      h_vx_neg 	->Fill(p.vx());
                      h_vy_neg 	->Fill(p.vy());
                      h_vz_neg 	->Fill(p.vz());
                      h_id_neg 	->Fill(p.pdgId());
		      }
     }


//------------------- DIGI --------------------------------------  
   for ( DigiContainerIterator<ME0DetId,ME0Digi> it = me0digiH->begin(); it != me0digiH->end(); ++it )
   {
    ME0DetId me0id = (*it).first;
    //cout << me0id << endl;
    MuonDigiCollection<ME0DetId,ME0Digi>::const_iterator itr;
    int bx = 99;
    for ( itr =((*it).second).first; itr!= ((*it).second).second; itr++ )
     {
      int strip = (*itr).strip();
      bx  = (*itr).bx();
      strips[me0id].push_back(strip);
      if (v)	cout << "strip: " << strip << "  bx: " << bx << endl;
     }
    
    me0detid_v1.push_back(me0id);
    cout << "region " << me0id.region() << "layer " << me0id.layer() << endl;
    if ( me0id.region() > 0 )  	{ 
      				  bx_pos = bx;
				  int layer = me0id.layer();
    				  if ( ! (std::find(me0detid_pos.begin(), me0detid_pos.end(),me0id)!=me0detid_pos.end()) )
    				      me0detid_pos.push_back( me0id ) ;
				  if ( ! (std::find(layers_pos.begin(), layers_pos.end(),layer)!=layers_pos.end()) )
    			  	     { layers_pos.push_back( layer ) ;
				       if (v) cout << "new layers_pos added : " << layer << endl;
				     }
	 			  if (layer==1) patternL_pos = 1  | patternL_pos;
	 			  if (layer==2) patternL_pos = 2  | patternL_pos;
	 			  if (layer==3) patternL_pos = 4  | patternL_pos;
	 			  if (layer==4) patternL_pos = 8  | patternL_pos;
	 			  if (layer==5) patternL_pos = 16 | patternL_pos;
	 			  if (layer==6) patternL_pos = 32 | patternL_pos;
			       	}
    else			{ 
      				  bx_neg = bx;
				  int layer = me0id.layer();
    				  if ( ! (std::find(me0detid_neg.begin(), me0detid_neg.end(),me0id)!=me0detid_neg.end()) )
    				      me0detid_neg.push_back( me0id ) ;
				  if ( ! (std::find(layers_neg.begin(), layers_neg.end(),layer)!=layers_neg.end()) )
    			  	     { layers_neg.push_back( layer ) ;
				       if (v) cout << "new layers_neg added : " << layer << endl;
				     }
    
    				  //me0detid_neg.push_back( me0id ) ;
    				  //layers_neg.push_back( me0id.layer() ) ;
				  //int layer = me0id.layer();
				  //if (v) cout << "new layers_neg added : " << layer << endl;
	 			  if (layer==1) patternL_neg = 1  | patternL_neg;
	 			  if (layer==2) patternL_neg = 2  | patternL_neg;
	 			  if (layer==3) patternL_neg = 4  | patternL_neg;
	 			  if (layer==4) patternL_neg = 8  | patternL_neg;
	 			  if (layer==5) patternL_neg = 16 | patternL_neg;
	 			  if (layer==6) patternL_neg = 32 | patternL_neg;

			       	}

   }

   patternLtree_pos = patternL_pos;	patternLtree_neg = patternL_neg;
   //make sure that strips and pads vectors in maps are sorted ascending
   for (auto itmap = strips.begin() ; itmap != strips.end() ; ++itmap)	
      std::sort( (itmap->second).begin() , (itmap->second).end() );
   for (auto itmap = pads.begin() ; itmap != pads.end() ; ++itmap)	
      std::sort( (itmap->second).begin() , (itmap->second).end() );
   //(above lines probably not necessary)

   if (v) { std::bitset<6> x(patternL_pos);
      		std::cout << "patternL_pos : " << x ;
		std::cout << "\t = " << patternL_pos;
          }
   if (v) { std::bitset<6> x(patternL_neg);
       		std::cout << "patternL_neg : " << x ;
		std::cout << "\t = " << patternL_neg;
          }
   if (v) {
	  cout << "\n\nPositive ME0DetId:" << endl;
	  for (unsigned int i=0; i<me0detid_pos.size(); i++)	
	     {
	     cout << me0detid_pos[i] << endl;
	     cout << "strips: " ;
	     for (auto it=strips[me0detid_pos[i]].begin(); it!=strips[me0detid_pos[i]].end(); ++it)	cout << " " << *it << " " ;
	     cout << endl;
	     cout << "pads: " ;
	     for (auto it=pads[me0detid_pos[i]].begin(); it!=pads[me0detid_pos[i]].end(); ++it)	cout << " " << *it << " " ;
	     cout << endl;
	     }
	  cout << "\n\nNegative ME0DetId:" << endl;
	  //for (unsigned int i=0; i<me0detid_neg.size(); i++)	cout << me0detid_neg[i] << endl;
	  for (unsigned int i=0; i<me0detid_neg.size(); i++)	
	     {
	     cout << me0detid_neg[i] << endl;
	     cout << "strips: " ;
	     for (auto it=strips[me0detid_neg[i]].begin(); it!=strips[me0detid_neg[i]].end(); ++it)	cout << " " << *it << " " ;
	     cout << endl;
	     cout << "pads: " ;
	     for (auto it=pads[me0detid_neg[i]].begin(); it!=pads[me0detid_neg[i]].end(); ++it)	cout << " " << *it << " " ;
	     cout << endl;
	     }
          }

int patternD_pos = 0;
int patternD_neg = 0;
int n = 0;
//-------------------- Fill nStrips, holeSize histograms -----------------------------------
   if (v)	cout << "\n\nPositive ME0DetId:" << endl;
   for (unsigned int i=0; i<me0detid_pos.size(); i++)	
      {
      patternD_pos = 0;
      n = 0;
      nStr_pos = strips[me0detid_pos[i]].size();
      h_nStr_pos->Fill(nStr_pos);
      if (v) { 
      	     cout << me0detid_pos[i] << endl;
             cout << "holesize: " ;
	     }
      patternD_pos = patternD_pos | static_cast<int>(pow(2,n));
      for (auto it= (strips[me0detid_pos[i]].rbegin()+1); it!= strips[me0detid_pos[i]].rend(); ++it)
         {
	 int hole = (*(it-1)-*(it)-1);
	 if (v)		cout << " ********************" << hole << " " ;
	 h_holeSize_pos->Fill(hole);
	 n = n + (hole+1);
	 patternD_pos = patternD_pos | static_cast<int>(pow(2,n));
         }
	 
      stripsPat[me0detid_pos[i]] = patternD_pos;
      patternDtree_pos = patternD_pos; 
      patternDtree_vec_pos.push_back(patternDtree_pos);
      h_digiPattern_pos->Fill(patternDtree_pos);

      if (v) { std::bitset<20> x(patternD_pos);
      		std::cout << "patternD_pos : " << x ;
		std::cout << "\t = " << patternD_pos;
              }

      cout << endl; 
      }

   if (v)	cout << "\n\nNegative ME0DetId:" << endl;
   for (unsigned int i=0; i<me0detid_neg.size(); i++)	
      {
      patternD_neg = 0;
      n = 0;
      nStr_neg = strips[me0detid_neg[i]].size();
      h_nStr_neg->Fill(nStr_neg);
      if (v) { 
      	     cout << me0detid_neg[i] << endl;
             cout << "holesize: " ;
	     }
      patternD_neg = patternD_neg | static_cast<int>(pow(2,n));
      for (auto it= (strips[me0detid_neg[i]].rbegin()+1); it!= strips[me0detid_neg[i]].rend(); ++it)
         {
	 int hole = (*(it-1)-*(it)-1);
	 if (v)		cout << " ********************" << hole << " " ;
	 h_holeSize_neg->Fill(hole);
	 n = n + (hole+1);
	 patternD_neg = patternD_neg | static_cast<int>(pow(2,n));
         }

      stripsPat[me0detid_neg[i]] = patternD_neg;
      patternDtree_neg = patternD_neg;
      patternDtree_vec_neg.push_back(patternDtree_neg);
      h_digiPattern_neg->Fill(patternDtree_neg);

      if (v) { std::bitset<20> x(patternD_neg);
      		std::cout << "patternD_neg : " << x ;
		std::cout << "\t = " << patternD_neg;
             }

      cout << endl; 
      }
//   for (unsigned int i=0; i<me0detid_neg.size(); i++)	
//      {
//      nStr_neg = strips[me0detid_neg[i]].size();
//      h_nStr_neg->Fill(nStr_neg);
//      if (v) { 
//      	     cout << me0detid_neg[i] << endl;
//             cout << "holesize: " ;
//	     }
//      for (auto it=strips[me0detid_neg[i]].rbegin(); it!= (strips[me0detid_neg[i]].rend()-1); ++it)
//         {
//	 int hole = ((*it)-*(it+1)-1);
//	 if (v)		cout << " ********************" << hole << " " ;
//	 h_holeSize_neg->Fill(hole);
//         }
//      cout << endl; 
//      }


//-----------------------Fill LAYER PATTERNS map --------------------
//if ( firedLayers_pos.count(patternL_pos) ) 	firedLayers_pos[patternL_pos]++;
//else					  	firedLayers_pos[patternL_pos] = 1;
//if ( firedLayers_neg.count(patternL_neg) ) 	firedLayers_neg[patternL_neg]++;
//else					  	firedLayers_neg[patternL_neg] = 1;

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
    if (v)	cout << me0id << endl;
    MuonDigiCollection<ME0DetId,ME0PadDigi>::const_iterator itr;
    for ( itr =((*it).second).first; itr!= ((*it).second).second; itr++ )
     {
      int pad = (*itr).pad();
      int bx    = (*itr).bx();
      pads[me0id].push_back(pad);
      if (v)	cout << "pad: " << pad << "  bx: " << bx << endl;
     }
    me0detid_v2.push_back(me0id);
   }

//------------------- PAD DIGI CLUSTER ------------------------------  
   for ( DigiContainerIterator<ME0DetId,ME0PadDigiCluster> it = me0padDigiClusterH->begin(); it != me0padDigiClusterH->end(); ++it )
   {
    ME0DetId me0id = (*it).first;
    if (v)	cout << me0id << endl;
    MuonDigiCollection<ME0DetId,ME0PadDigiCluster>::const_iterator itr;
    for ( itr =((*it).second).first; itr!= ((*it).second).second; itr++ )
     {
      vector<uint16_t> pads = (*itr).pads();
      //int bx    = (*itr).bx();
      clusters[me0id].push_back(pads);
      if (v)	cout << "pads: " ;
      for ( unsigned int i=0; i<pads.size(); i++)  cout << pads[i] << " " ;
     }
    me0detid_v3.push_back(me0id);
   }

    
   tr->Fill();

//-------------------- CHECKS ----------------------
if ( nL_pos > 6 )	{ cout << "ERROR! More than 6 fired layers found on pos endcap." << endl; return; }
if ( nL_neg > 6 )	{ cout << "ERROR! More than 6 fired layers found on neg endcap." << endl; return; }
//Verify that the ME0DetId list is the same in the three cases above
if( equal(me0detid_v1.begin(), me0detid_v1.end(), me0detid_v2.begin()) )
  if( equal(me0detid_v1.begin(), me0detid_v1.end(), me0detid_v3.begin()) )
    if ( me0detid_v1.size() == me0detid_v2.size() && me0detid_v1.size() == me0detid_v3.size() )
      { cout << "The three classes have the same list of ME0DetId" << endl; return; }

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
