// -*- C++ -*-
//
// Package:    Analyzer/DemoAnalyzer
// Class:      DemoAnalyzer
// 
/**\class DemoAnalyzer DemoAnalyzer.cc Analyzer/DemoAnalyzer/plugins/DemoAnalyzer.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Georgios Krintiras
//         Created:  Fri, 13 Mar 2015 11:02:38 GMT
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "DQMServices/Core/interface/DQMEDAnalyzer.h"
#include "DQMServices/Core/interface/MonitorElement.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"


#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/TrackerCommon/interface/TrackerTopology.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "SimDataFormats/TrackingHit/interface/PSimHit.h"
#include "SimDataFormats/Track/interface/SimTrack.h"
#include "SimDataFormats/Track/interface/SimTrackContainer.h"
#include "SimGeneral/HepPDTRecord/interface/PdtEntry.h"

#include "DataFormats/Math/interface/Vector3D.h"
#include "Geometry/Records/interface/IdealGeometryRecord.h"
#include "Geometry/Records/interface/TrackerDigiGeometryRecord.h"
#include "Geometry/TrackerGeometryBuilder/interface/TrackerGeometry.h"

#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "TH2.h"
#include "TMath.h"
#include "TString.h"
#include <list>
#include <map>
#include <vector>
#include <set>
#include <sstream>
#include <stdlib.h> 
#include <string>

//
// class declaration
// migrate to:https://twiki.cern.ch/twiki/bin/viewauth/CMS/ThreadedDQM
//https://github.com/cirkovic/my-cmssw/blob/master/DQMOffline/RecoB/plugins/BTagPerformanceAnalyzerOnData.cc

class DemoAnalyzer : public DQMEDAnalyzer {
   public:
      explicit DemoAnalyzer(const edm::ParameterSet&);
      ~DemoAnalyzer();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


   private:
  
      virtual void dqmBeginRun(edm::Run const& ,edm::EventSetup const& ) override;
      virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
      void bookHistograms(DQMStore::IBooker &, edm::Run const &, edm::EventSetup const &) override;
  
      // ----------member data ---------------------------
      bool processFastSim_;
      void bookHistosPerParticle( const std::string& particles,  DQMStore::IBooker & ibooker );


  
      void bookEnergyLossesRelatedInfo2D( std::vector<MonitorElement*>&, DQMStore::IBooker & ibooker,
					  int nBinsX, float rangeX,  int nBinsY, float rangeY, 
					  const TString &varX, const TString &varY, const TString &det, unsigned int nHistos );

      void bookEnergyLosses1D( std::vector<MonitorElement*>&, DQMStore::IBooker & ibooker, int nBins, float range, 
			       const TString &var, const TString &det, unsigned int nHistos );

      std::vector<edm::InputTag> simHitsTag_;
      std::string particles_;
      std::vector<std::string> bins_E;
      std::vector<std::string> bins_p;

      std::vector<float> bins1_;
  
  
  /*
      edm::InputTag simHitsTag_PXB_Lowtof_;
      edm::InputTag simHitsTag_PXB_Hightof_;
      edm::InputTag simHitsTag_PXF_Lowtof_;
      edm::InputTag simHitsTag_PXF_Hightof_;
      edm::InputTag simHitsTag_TEC_Lowtof_;
      edm::InputTag simHitsTag_TEC_Hightof_;
      edm::InputTag simHitsTag_TIB_Lowtof_;
      edm::InputTag simHitsTag_TIB_Hightof_;
      edm::InputTag simHitsTag_TID_Lowtof_;
      edm::InputTag simHitsTag_TID_Hightof_;
      edm::InputTag simHitsTag_TOB_Lowtof_;
      edm::InputTag simHitsTag_TOB_Hightof_;
  */

      //edm::Service<TFileService> fs_;
  /*
      //Funtions
      void bookHistosPerDetector();
      void bookEnergyLosses1D( std::vector<TH1F*>&, int nBins, float range, const TString &det, unsigned int nHistos );
  
      void bookEnergyLossesRelatedInfo2D( std::vector<TH2F*>&, 
					  int nBinsX, float rangeX,  int nBinsY, float rangeY, 
					  const TString &var, const TString &det, unsigned int nHistos );
  
      void bookEnergyLossesRelatedInfo2D( std::vector<TH2F*>&,
					  int nBinsX, float rangeX,  int nBinsY, float rangeY,
					  const TString &var, const TString &det );
  */
      //Helper Functions
      void identifyToken(std::set<char> & , const std::string& );
      void removeWhiteSpaces( std::string &strg );

      std::vector<std::vector<MonitorElement*> > histos_protons_p_dedx_;
      typedef std::vector<std::vector<MonitorElement*> > MyClassSetVector;
      typedef std::map<std::string,  MyClassSetVector > MyClassSetMap;
      MyClassSetMap my_map1;
      MyClassSetMap my_map2;


  //typedef typename std::vector<T>::iterator iterator;

      std::vector<std::vector<MonitorElement*> > histos_protons_dedx_;
  
      std::vector<std::vector<MonitorElement*> > histos_pions_p_dedx_;
      std::vector<std::vector<MonitorElement*> > histos_pions_dedx_;


      //  Detectors
      std::map<std::string, int> map_subdet_nlayers_;

     
  /*
      // Histograms      
      //Generator Level
      TH1F* h_Genpart_Mother_Pt_;
      TH1F* h_Genpart_Daughter_Pt_;
      
      // Simulation Level 
      ////SimTrack Class
      TH1F* h_Track_P_;
      TH1F* h_Track_Pt_;
      TH1F* h_Track_Phi_;
      TH1F* h_Track_Eta_;
      TH1I* h_Track_Multiplicity_;
      TH2F* h_Track_ID_Track_Pt_;
      ////PSimHit Class
      TH1F* h_Hit_Pixels_P_AtEntry_;
      TH1F* h_Hit_Strips_P_AtEntry_;
      TH1I* h_Hit_PXBPXF_PerLayer_Multiplicity_;
      TH1I* h_Hit_PixelDets_Multiplicity_;
      TH1I* h_Hit_StripDets_Multiplicity_;
      TH2F* h_Hit_PixelDets_Map_;
      TH2F* h_Hit_StripDets_Map_;
      TH2F* h_Hit_PixelDets_dedx_P_AtEntry_;
      TH2F* h_Hit_StripDets_dedx_P_AtEntry_;
      ////PSimHit against SimTrack Class
      TH2F* h_Hit_PixelDets_dedx_Track_P_;
      TH2F* h_Hit_StripDets_dedx_Track_P_;
      TH2F* h_Hit_PixelDets_Multiplicity_Track_P_;
      TH2F* h_Hit_StripDets_Multiplicity_Track_P_;

      //  Detectors
      std::map<std::string, int> map_subdet_nlayers_;
      ////Pixel Barrel - 3 different detectors
      static const unsigned int nHistos_PXB_ = 3;
      std::vector<TH1F*> histos_PXB_dedx_;
      std::vector<TH2F*> histos_PDG_PXB_dedx_;

      std::vector<std::vector<TH2F*> > histos_pions_p_dedx_;
      std::vector<TH2F*> histos_PDG_PXB_dx_;
      std::vector<TH2F*> histos_p_PXB_dedx_;
      ////Pixel Endcap - 2 different detectors
      static const unsigned int nHistos_PXF_ = 2;
      std::vector<TH1F*> histos_PXF_dedx_;
      std::vector<TH2F*> histos_PDG_PXF_dedx_;
      std::vector<TH2F*> histos_PDG_PXF_dx_;
      std::vector<TH2F*> histos_p_PXF_dedx_;
      ////TIB - 4 different detectors
      static const unsigned int nHistos_TIB_ = 4;
      std::vector<TH1F*> histos_TIB_dedx_;
      std::vector<TH2F*> histos_PDG_TIB_dedx_;
      std::vector<TH2F*> histos_PDG_TIB_dx_;
      std::vector<TH2F*> histos_p_TIB_dedx_;
      ////TOB - 6 different detectors
      static const unsigned int nHistos_TOB_ = 6;
      std::vector<TH1F*> histos_TOB_dedx_;
      std::vector<TH2F*> histos_PDG_TOB_dedx_;
      std::vector<TH2F*> histos_PDG_TOB_dx_;
      std::vector<TH2F*> histos_p_TOB_dedx_;
      ////TID - 3 different detectors
      static const unsigned int nHistos_TID_ = 3;
      std::vector<TH1F*> histos_TID_dedx_;
      std::vector<TH2F*> histos_PDG_TID_dedx_;
      std::vector<TH2F*> histos_PDG_TID_dx_;
      std::vector<TH2F*> histos_p_TID_dedx_;
      ////TEC - 9 different detectors
      static const unsigned int nHistos_TEC_ = 7;
      std::vector<TH1F*> histos_TEC_dedx_;
      std::vector<TH2F*> histos_PDG_TEC_dedx_;
      std::vector<TH2F*> histos_PDG_TEC_dx_;
      std::vector<TH2F*> histos_p_TEC_dedx_;
  */
      // PSimHits
      std::vector<edm::InputTag> trackerContainers_;
      //
    
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
DemoAnalyzer::DemoAnalyzer(const edm::ParameterSet& iConfig):

  simHitsTag_(iConfig.getParameter<std::vector<edm::InputTag> >("SimHitTagLabels")),
  particles_(),
  bins_E(),
  bins_p()
  //fNbEvtToPrint(iConfig.getUntrackedParameter<std::vector<float> >("NbEventsToPrint"))
  

  //fs_()


{
    //now do what ever initialization is needed
  /*
    //Generator Level
    h_Genpart_Mother_Pt_ = fs_->make<TH1F>("Genpart_Mother_Pt", ";p_{T}", 200, 0., 12.); //ok
    h_Genpart_Daughter_Pt_ = fs_->make<TH1F>("Genpart_Daughter_Pt", ";p_{T}", 200, 0., 12.); //ok
  
    //Simulation Level 
    ////SimTrack Class
    h_Track_P_ = fs_->make<TH1F>("SimTrack_P", ";P", 200, 0., 100.); //ok
    h_Track_Pt_ = fs_->make<TH1F>("SimTrack_Pt", ";#p_{T}", 200, 0., 100.); //ok
    h_Track_Phi_ = fs_->make<TH1F>("SimTrack_Phi", ";#phi", 200, -TMath::Pi(), TMath::Pi()); //ok
    h_Track_Eta_ = fs_->make<TH1F>("SimTrack_Eta", ";#eta", 200, -2.5, 2.5); //ok
    h_Track_Multiplicity_ = fs_->make<TH1I>("SimTrack_Multiplicity", ";N_{Tracks}", 29, 1, 30); //ok
    h_Track_ID_Track_Pt_ = fs_->make<TH2F>("Track_ID_Track_Pt", ";P;PDG", 2400, 0., 100., 2*2212, 0., 2213.); //ok
    
    ////PSimHit Class
    h_Hit_Pixels_P_AtEntry_ = fs_->make<TH1F>("PSimHit_Pixels_P_AtEntry", ";P", 200, 0., 100.); //ok
    h_Hit_Strips_P_AtEntry_ = fs_->make<TH1F>("PSimHit_Strips_P_AtEntry", ";P", 200, 0., 100.); //ok
    
    h_Hit_PXBPXF_PerLayer_Multiplicity_ = fs_->make<TH1I>("PSimHit_Multiplicity_PXBPXF", ";Pixel Layer", 5, 1, 6); //ok
    h_Hit_PixelDets_Multiplicity_ = fs_->make<TH1I>("PSimHit_AllPixelLayers", ";PSimHits_{Pixels}", 29, 1, 30); //ok
    h_Hit_StripDets_Multiplicity_ = fs_->make<TH1I>("PSimHit_AllStripLayers", ";PSimHits_{Strips}", 29, 1, 30); //ok
    
    h_Hit_PixelDets_Map_ = fs_->make<TH2F>("PSimHit_AllPixelLayers_Map", ";z;y", 240, -300., 300., 240, -140., 140.); //ok
    h_Hit_StripDets_Map_ = fs_->make<TH2F>("PSimHit_AllStripLayers_Map", ";z;y", 240, -300., 300., 240, -140., 140.); //ok
    
    h_Hit_PixelDets_dedx_P_AtEntry_ = fs_->make<TH2F>("PSimHit_AllPixelLayers_dedx_P_AtEntry", ";P;dedx", 240, 0., 100., 240, 0., .01); //ok
    h_Hit_StripDets_dedx_P_AtEntry_ = fs_->make<TH2F>("PSimHit_AllStripLayersd_dedx_P_AtEntry", ";P;dedx", 240, 0., 100., 240, 0., .01); //ok
    
    
    ////PSimHit against SimTrack Class
    h_Hit_PixelDets_Multiplicity_Track_P_ = fs_->make<TH2F>("PSimHit_AllPixelLayers_Multiplicity_SimTrack_P", ";PSimHits_{Pixels};P",
							    29, 1, 30, 640, 0., 100.); //ok
    h_Hit_StripDets_Multiplicity_Track_P_ = fs_->make<TH2F>("PSimHit_AllStripLayers_Multiplicity_SimTrack_P", ";PSimHits_{Strips};P",
							    29, 1, 30, 640, 0., 100.); //ok
    
    h_Hit_PixelDets_dedx_Track_P_ = fs_->make<TH2F>("PSimHit_AllPixelLayers_dedx_Track_P", ";P;dedx", 640, 0., 100., 240, 0., .01); //ok
    h_Hit_StripDets_dedx_Track_P_ = fs_->make<TH2F>("PSimHit_AllStripLayers_dedx_Track_P", ";P;dedx", 640, 0., 100., 240, 0., .01); //ok
    
    bookHistosPerDetector();
    map_subdet_nlayers_.insert(std::pair<std::string,int>("PXB",3));
    map_subdet_nlayers_.insert(std::pair<std::string,int>("PXF",2));

    bookHistosPerParticle("protons:pions");
  */
    processFastSim_ = iConfig.getUntrackedParameter<bool>("processFastSim");
    if (processFastSim_)
      trackerContainers_.push_back(simHitsTag_[0]);
    
    else {
      for (unsigned int i_subDet=0; i_subDet<simHitsTag_.size(); ++i_subDet){
        trackerContainers_.push_back(simHitsTag_[i_subDet]);
	/*
	trackerContainers_.push_back(simHitsTag_PXB_Hightof_);
	trackerContainers_.push_back(simHitsTag_PXF_Lowtof_);
	trackerContainers_.push_back(simHitsTag_PXF_Hightof_);
	trackerContainers_.push_back(simHitsTag_TEC_Lowtof_);
	trackerContainers_.push_back(simHitsTag_TEC_Hightof_);
	trackerContainers_.push_back(simHitsTag_TIB_Lowtof_);
	trackerContainers_.push_back(simHitsTag_TIB_Hightof_);
	trackerContainers_.push_back(simHitsTag_TID_Lowtof_);
	trackerContainers_.push_back(simHitsTag_TID_Hightof_);
	trackerContainers_.push_back(simHitsTag_TOB_Lowtof_);
	trackerContainers_.push_back(simHitsTag_TOB_Hightof_);
	*/
      }

    }

    particles_ = iConfig.getUntrackedParameter <std::string > ("particles");
    bins_E = iConfig.getUntrackedParameter <std::vector <std::string> > ("bins_E", std::vector<std::string>());
    bins_p = iConfig.getUntrackedParameter <std::vector <std::string> > ("bins_p", std::vector<std::string>());

    //if (!bins_p.empty()) 
    if (iConfig.exists("bins_p")) {
      sort(bins_p.begin(),bins_p.end());	
	
      for (unsigned int it=0; it!=bins_p.size(); ++it)
	bins1_.push_back(std::atof(bins_p[it].data()));
	//std::cout << " !!!!!!!!!!!!!!!!!!!!!!!!!!!!!! " << bins_[it]<< std::endl;
      sort(bins1_.begin(),bins1_.end());	
    }
    
    /*
    std::list<float>::iterator it;

      std::cout << ' ' << *it;
    */

}

DemoAnalyzer::~DemoAnalyzer()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
DemoAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
#ifdef rrDEBUG
  std::cout << "Famos analysis" << std::endl;
#endif
  // get event and run number
#ifdef rrDEBUG
  int t_Run   = iEvent.id().run();
  int t_Event = iEvent.id().event();
  std::cout
    << " #################################### Run " << t_Run 
    << " Event "                                    << t_Event << " #################################### " 
    << std::endl;
#endif

     using namespace edm;

    //generator level particles
    Handle<reco::GenParticleCollection> genParticles;
    iEvent.getByLabel("genParticles", genParticles);
    /*
    for(unsigned int iPart = 0; iPart < genParticles->size(); ++ iPart) 
    {
        const reco::GenParticle & mPart = (*genParticles)[iPart];
	h_Genpart_Mother_Pt_->Fill(mPart.pt());
	size_t n = mPart.numberOfDaughters();
	for(unsigned int jPart = 0; jPart < n; ++ jPart) 
        {
	    const reco::Candidate * dPart = mPart.daughter(jPart);
	    h_Genpart_Daughter_Pt_->Fill(dPart->pt());
	}
    }
   
    ESHandle< TrackerGeometry > TG;
    iSetup.get<TrackerDigiGeometryRecord>().get(TG);
    const TrackerGeometry *theTrackerGeometry = TG.product ();
    */

    
    // Retrieve tracker topology from geometry
    ESHandle< TrackerTopology > tTopoHandle;
    iSetup.get<IdealGeometryRecord>().get(tTopoHandle);
    const TrackerTopology* const tTopo = tTopoHandle.product();
    

    //Retrieve the Monte Carlo truth Particles(SimTracks)
    Handle< SimTrackContainer > simTracksHandle;
    iEvent.getByLabel(trackerContainers_[0].label(),simTracksHandle); 
    const SimTrackContainer simTracks = *(simTracksHandle.product());

    std::vector<float>::iterator it;

    for ( SimTrackContainer::const_iterator  simTrack = simTracks.begin(); simTrack != simTracks.end(); ++simTrack) 
    {
      /*
        h_Track_P_->Fill((*simTrack).momentum().P());
	h_Track_Pt_->Fill((*simTrack).momentum().Pt());
	h_Track_Phi_->Fill((*simTrack).momentum().Phi());
	h_Track_Eta_->Fill((*simTrack).momentum().Eta());
	h_Track_ID_Track_Pt_->Fill((*simTrack).momentum().Pt(), fabs((*simTrack).type()));
      */
	unsigned int i_SimHit_Multiplicity = 0, i_SimHit_Pixels_Multiplicity = 0;

	for (unsigned i=0; i<trackerContainers_.size(); ++i) 
	{
           Handle<std::vector<PSimHit> > simHitCollection;
	   iEvent.getByLabel(trackerContainers_[i], simHitCollection);
	   const std::vector<PSimHit>& simHits = *simHitCollection.product();
	   if(simHitCollection.isValid()) 
	   {
	       for( std::vector<PSimHit>::const_iterator Hits=simHits.begin(); Hits!=simHits.end(); ++Hits) {
	           if((*Hits).trackId() == (*simTrack).trackId()) 
		   {
		       i_SimHit_Multiplicity++;
	     
		       std::set<unsigned int> detIds;
		       unsigned int detId = (*Hits).detUnitId();
		       unsigned int isub  = DetId(detId).subdetId();
		       
		       
		       float f_dx = TMath::Power(
						 TMath::Power( (*Hits).entryPoint().x() - (*Hits).exitPoint().x(), 2) + 
						 TMath::Power( (*Hits).entryPoint().y() - (*Hits).exitPoint().y(), 2) + 
						 TMath::Power( (*Hits).entryPoint().z() - (*Hits).exitPoint().z(), 2), 
						 1/2.);
		       
		       
		       //float f_dx = TMath::Power(TMath::Power( (*Hits).entryPoint().z() - (*Hits).exitPoint().z(), 2), 
		       //1/2.);
		       //const DetId& theDetId = (*Hits).detUnitId();
		       //const GeomDet* theGeomDet = theTrackerGeometry->idToDet(theDetId);
		       //const GlobalPoint& globalPoint = theGeomDet->toGlobal (localPoint); 
		       //std::cout<< globalPoint.x() <<std::endl;
		       float momentumAtEntry = (*Hits).momentumAtEntry().mag();
		       std::cout<<momentumAtEntry<<" prin "<<bins1_[1]-bins1_[0]<<" "<<abs((*Hits).particleType())<<std::endl;
		       //if( momentumAtEntry<TMath::Abs((bins1_[0])) ||  momentumAtEntry>TMath::Abs((bins1_[bins1_.size()-1])) ) continue;
		       bins1_.push_back(momentumAtEntry);
		       sort(bins1_.begin(),bins1_.end());	
		       it=find(bins1_.begin(),bins1_.end(), momentumAtEntry);
		       auto pos = std::distance(bins1_.begin(), it);
		       
		       if (isub == static_cast<int>(PixelSubdetector::PixelBarrel) || isub == static_cast<int>(PixelSubdetector::PixelEndcap)) 
		       {
			 //std::cout<<std::to_string((*Hits).momentumAtEntry().mag()).data()<<std::endl;
			 

			 //sort(bins_.begin(),bins_.end());
		       /*
		           h_Hit_PXBPXF_PerLayer_Multiplicity_->Fill(tTopo->pxbLayer(detId));
			   h_Hit_Pixels_P_AtEntry_->Fill((*Hits).momentumAtEntry().mag());
			   h_Hit_PixelDets_dedx_P_AtEntry_->Fill((*Hits).momentumAtEntry().mag(), (*Hits).energyLoss()/f_dx);
			   h_Hit_PixelDets_dedx_Track_P_->Fill((*simTrack).momentum().P(), (*Hits).energyLoss()/f_dx);
			   h_Hit_PixelDets_Map_->Fill(theGeomDet->toGlobal((*Hits).localPosition()).z(),
						      theGeomDet->toGlobal((*Hits).localPosition()).y());
			   histos_PXB_dedx_[tTopo->pxbLayer(detId)-1]->Fill((*Hits).energyLoss()/f_dx);
			   histos_PDG_PXB_dedx_[tTopo->pxbLayer(detId)-1]->Fill(abs((*Hits).particleType()), (*Hits).energyLoss()/f_dx);
			   histos_PDG_PXB_dx_[tTopo->pxbLayer(detId)-1]->Fill(abs((*Hits).particleType()), f_dx);
			   histos_p_PXB_dedx_[tTopo->pxbLayer(detId)-1]->Fill((*Hits).pabs(), (*Hits).energyLoss()/f_dx);
			 */
			   i_SimHit_Pixels_Multiplicity++;
			   if (abs((*Hits).particleType())==2212)
			   {  
			     //histos_protons_p_dedx_[0][tTopo->pxbLayer(detId)-1]->Fill((*Hits).pabs(), (*Hits).energyLoss()/f_dx);
			     if (pos!=0 && pos<int(bins1_.size()-1)){
			       std::cout<<momentumAtEntry<<" meta "<<pos<<" meta "<<bins1_.size()<<" "<<(*Hits).particleType()<<std::endl;
			       //histos_protons_dedx_[0][pos-1]->Fill((*Hits).energyLoss()/f_dx);
			       typedef std::map<std::string,  std::vector<std::vector<MonitorElement*> > >::iterator iter;
			
				for( iter iterator = my_map1.begin(); iterator != my_map1.end(); ++iterator ) {
				  MyClassSetVector foo = iterator->second;
				  std::cout<<" "<<iterator->first<<std::endl;    
				  //for (unsigned i =0; i<foo.size(); i++) {
				  foo[0][pos-1]->Fill((*Hits).pabs(), (*Hits).energyLoss()/f_dx);
				  //std::cout <<  " " << foo[0][0]->Fill(0) << std::endl;
				  //}
				}
			       
			     }
			     histos_protons_p_dedx_[0][0]->Fill((*Hits).pabs(), (*Hits).energyLoss()/f_dx);
			     
			     typedef std::map<std::string,  std::vector<std::vector<MonitorElement*> > >::iterator iter;
    			     for( iter iterator = my_map2.begin(); iterator != my_map2.end(); ++iterator ) {
			       MyClassSetVector foo = iterator->second;
			       std::cout<<" "<<iterator->first<<std::endl;    
			       //for (unsigned i =0; i<foo.size(); i++) {
				 foo[0][0]->Fill((*Hits).pabs(), (*Hits).energyLoss()/f_dx);
				 //std::cout <<  " " << foo[0][0]->Fill(0) << std::endl;
				 //}
			     }
			     

			   }
			   else if (abs((*Hits).particleType())==211)
			   {
			     //histos_pions_p_dedx_[0][tTopo->pxbLayer(detId)-1]->Fill((*Hits).pabs(), (*Hits).energyLoss()/f_dx);
			     if (pos!=0 && pos<int(bins1_.size()-1)){
			       std::cout<<momentumAtEntry<<" meta "<<pos<<" meta "<<bins1_.size()<<" "<<(*Hits).particleType()<<std::endl;
			       histos_pions_dedx_[0][pos-1]->Fill((*Hits).energyLoss()/f_dx);
			     }
			     histos_pions_p_dedx_[0][0]->Fill((*Hits).pabs(), (*Hits).energyLoss()/f_dx);

			   }
			   			   
		       }

		       else if (isub == StripSubdetector::TIB || isub == StripSubdetector::TOB || 
				isub == StripSubdetector::TID || isub == StripSubdetector::TEC )
			 {
			   //if(momentumAtEntry<(bins1_[1]-bins2_[0])) continue;

			    if (abs((*Hits).particleType())==2212)
			   {  
			     if (isub == StripSubdetector::TEC && tTopo->tecRing(detId)>=5) {
			       if (pos!=0 && pos<int(bins1_.size()-1)){
				 std::cout<<momentumAtEntry<<" meta "<<pos<<" meta "<<bins1_.size()<<" "<<(*Hits).particleType()<<std::endl;
				 //histos_protons_dedx_[2][pos-1]->Fill((*Hits).energyLoss()/f_dx);
				 typedef std::map<std::string,  std::vector<std::vector<MonitorElement*> > >::iterator iter;
				 
				 for( iter iterator = my_map1.begin(); iterator != my_map1.end(); ++iterator ) {
				   MyClassSetVector foo = iterator->second;
				   std::cout<<" "<<iterator->first<<std::endl;    
				   //for (unsigned i =0; i<foo.size(); i++) {
				   foo[2][pos-1]->Fill((*Hits).pabs(), (*Hits).energyLoss()/f_dx);
				   //std::cout <<  " " << foo[0][0]->Fill(0) << std::endl;
				   //}
				}
			       }
			       histos_protons_p_dedx_[2][0]->Fill((*Hits).pabs(), (*Hits).energyLoss()/f_dx);
			     }
			     
			     else {
			       if (pos!=0 && pos<int(bins1_.size()-1)){
				 std::cout<<momentumAtEntry<<" meta "<<pos<<" meta "<<bins1_.size()<<" "<<(*Hits).particleType()<<std::endl;
				 //histos_protons_dedx_[1][pos-1]->Fill((*Hits).energyLoss()/f_dx);
				  typedef std::map<std::string,  std::vector<std::vector<MonitorElement*> > >::iterator iter;
				  for( iter iterator = my_map1.begin(); iterator != my_map1.end(); ++iterator ) {
				  MyClassSetVector foo = iterator->second;
				  std::string particle = iterator->first;
				  std::cout<<" "<<iterator->first<<std::endl;    
				  //for (unsigned i =0; i<foo.size(); i++) {
				  if (particle.compare("protons")==0)
				    foo[1][pos-1]->Fill((*Hits).pabs(), (*Hits).energyLoss()/f_dx);
				  //std::cout <<  " " << foo[0][0]->Fill(0) << std::endl;
				  //}
				  }
			       histos_protons_p_dedx_[1][0]->Fill((*Hits).pabs(), (*Hits).energyLoss()/f_dx);
			       }
			     }
			   }
			   else if (abs((*Hits).particleType())==211)
			   {
			     //histos_pions_p_dedx_[0][tTopo->pxbLayer(detId)-1]->Fill((*Hits).pabs(), (*Hits).energyLoss()/f_dx);
			     if (isub == StripSubdetector::TEC && tTopo->tecRing(detId)>=5){
			       if (pos!=0 && pos<int(bins1_.size()-1)){
				 std::cout<<momentumAtEntry<<" meta "<<pos<<" meta "<<bins1_.size()<<" "<<(*Hits).particleType()<<std::endl;
				 //histos_pions_dedx_[2][pos-1]->Fill((*Hits).energyLoss()/f_dx);
				  typedef std::map<std::string,  std::vector<std::vector<MonitorElement*> > >::iterator iter;
				  for( iter iterator = my_map1.begin(); iterator != my_map1.end(); ++iterator ) {
				  MyClassSetVector foo = iterator->second;
				  std::cout<<" "<<iterator->first<<std::endl;    
				  //for (unsigned i =0; i<foo.size(); i++) {
				  foo[2][pos-1]->Fill((*Hits).pabs(), (*Hits).energyLoss()/f_dx);
				  //std::cout <<  " " << foo[0][0]->Fill(0) << std::endl;
				  //}
				  }
			       }
				  histos_pions_p_dedx_[2][0]->Fill((*Hits).pabs(), (*Hits).energyLoss()/f_dx);
			     }
			     else {
			       if (pos!=0 && pos<int(bins1_.size()-1)){
				 std::cout<<momentumAtEntry<<" meta "<<pos<<" meta "<<bins1_.size()<<" "<<(*Hits).particleType()<<std::endl;
				 histos_pions_dedx_[1][pos-1]->Fill((*Hits).energyLoss()/f_dx);
			       }
			       histos_pions_p_dedx_[1][0]->Fill((*Hits).pabs(), (*Hits).energyLoss()/f_dx);
			     }
			   }
			 }
		   
		   bins1_.erase(std::remove(bins1_.begin(), bins1_.end(), momentumAtEntry), bins1_.end());
		   }
	       }
	   }
	}
    }
    /*
		       else if (isub == static_cast<int>(PixelSubdetector::PixelEndcap)) 
		       {
			 
			   h_Hit_PXBPXF_PerLayer_Multiplicity_->Fill(3 + tTopo->pxfDisk(detId));
			   h_Hit_Pixels_P_AtEntry_->Fill((*Hits).momentumAtEntry().mag());
			   h_Hit_PixelDets_dedx_P_AtEntry_->Fill((*Hits).momentumAtEntry().mag(), (*Hits).energyLoss()/f_dx);
			   h_Hit_PixelDets_dedx_Track_P_->Fill((*simTrack).momentum().P(), (*Hits).energyLoss()/f_dx);
			   h_Hit_PixelDets_Map_->Fill(theGeomDet->toGlobal((*Hits).localPosition()).z(),
						      theGeomDet->toGlobal((*Hits).localPosition()).y());
			   histos_PXF_dedx_[tTopo->pxfDisk(detId)-1]->Fill((*Hits).energyLoss()/f_dx);
			   histos_PDG_PXF_dedx_[tTopo->pxfDisk(detId)-1]->Fill(abs((*Hits).particleType()), (*Hits).energyLoss()/f_dx);
			   histos_PDG_PXF_dx_[tTopo->pxfDisk(detId)-1]->Fill(abs((*Hits).particleType()), f_dx);
			   histos_p_PXF_dedx_[tTopo->pxfDisk(detId)-1]->Fill((*Hits).pabs(), (*Hits).energyLoss()/f_dx);
			 
			   if (abs((*Hits).particleType())==2212)
			   {
                               histos_protons_p_dedx_[1][tTopo->pxfDisk(detId)-1]->Fill((*Hits).pabs(), (*Hits).energyLoss()/f_dx);
			   }
                           else if (abs((*Hits).particleType())==211)
			   {
			       histos_pions_p_dedx_[1][tTopo->pxfDisk(detId)-1]->Fill((*Hits).pabs(), (*Hits).energyLoss()/f_dx);
			   }

			   i_SimHit_Pixels_Multiplicity++;
		       }
		   }
	       }
	   }
	}
    }

		       else if (isub == StripSubdetector::TIB) 
		       {
		           h_Hit_Strips_P_AtEntry_->Fill((*Hits).momentumAtEntry().mag());
			   h_Hit_StripDets_dedx_P_AtEntry_->Fill((*Hits).momentumAtEntry().mag(), (*Hits).energyLoss()/f_dx);
			   h_Hit_StripDets_dedx_Track_P_->Fill((*simTrack).momentum().P(), (*Hits).energyLoss()/f_dx);
			   h_Hit_StripDets_Map_->Fill(theGeomDet->toGlobal((*Hits).localPosition()).z(),
						      theGeomDet->toGlobal((*Hits).localPosition()).y());
			   histos_TIB_dedx_[tTopo->tibLayer(detId)-1]->Fill((*Hits).energyLoss()/f_dx);
			   histos_PDG_TIB_dedx_[tTopo->tibLayer(detId)-1]->Fill(abs((*Hits).particleType()), (*Hits).energyLoss()/f_dx);
			   histos_PDG_TIB_dx_[tTopo->tibLayer(detId)-1]->Fill(abs((*Hits).particleType()), f_dx);
			   histos_p_TIB_dedx_[tTopo->tibLayer(detId)-1]->Fill((*Hits).pabs(), (*Hits).energyLoss()/f_dx);
		       }
		       else if (isub == StripSubdetector::TOB) 
		       {
			   h_Hit_Strips_P_AtEntry_->Fill((*Hits).momentumAtEntry().mag());
			   h_Hit_StripDets_dedx_P_AtEntry_->Fill((*Hits).momentumAtEntry().mag(), (*Hits).energyLoss()/f_dx);
			   h_Hit_StripDets_dedx_Track_P_->Fill((*simTrack).momentum().P(), (*Hits).energyLoss()/f_dx);
			   h_Hit_StripDets_Map_->Fill(theGeomDet->toGlobal((*Hits).localPosition()).z(),
						      theGeomDet->toGlobal((*Hits).localPosition()).y());
			   histos_TOB_dedx_[tTopo->tobLayer(detId)-1]->Fill((*Hits).energyLoss()/f_dx); 
			   histos_PDG_TOB_dedx_[tTopo->tobLayer(detId)-1]->Fill(abs((*Hits).particleType()), (*Hits).energyLoss()/f_dx);
			   histos_PDG_TOB_dx_[tTopo->tobLayer(detId)-1]->Fill(abs((*Hits).particleType()), f_dx);
			   histos_p_TOB_dedx_[tTopo->tobLayer(detId)-1]->Fill((*Hits).pabs(), (*Hits).energyLoss()/f_dx);
		       }
		       else if (isub == StripSubdetector::TID) 
		       {
			   h_Hit_Strips_P_AtEntry_->Fill((*Hits).momentumAtEntry().mag());
			   h_Hit_StripDets_dedx_P_AtEntry_->Fill((*Hits).momentumAtEntry().mag(), (*Hits).energyLoss()/f_dx);
			   h_Hit_StripDets_dedx_Track_P_->Fill((*simTrack).momentum().P(), (*Hits).energyLoss()/f_dx);
			   h_Hit_StripDets_Map_->Fill(theGeomDet->toGlobal((*Hits).localPosition()).z(),
						      theGeomDet->toGlobal((*Hits).localPosition()).y());
			   histos_TID_dedx_[tTopo->tidRing(detId)-1]->Fill((*Hits).energyLoss()/f_dx);
			   histos_PDG_TID_dedx_[tTopo->tidRing(detId)-1]->Fill(abs((*Hits).particleType()), (*Hits).energyLoss()/f_dx);
			   histos_PDG_TID_dx_[tTopo->tidRing(detId)-1]->Fill(abs((*Hits).particleType()), f_dx);
			   histos_p_TID_dedx_[tTopo->tidRing(detId)-1]->Fill((*Hits).pabs(), (*Hits).energyLoss()/f_dx);
		       }
		       else if (isub == StripSubdetector::TEC) 
		       {
			   h_Hit_Strips_P_AtEntry_->Fill((*Hits).momentumAtEntry().mag());
			   h_Hit_StripDets_dedx_P_AtEntry_->Fill((*Hits).momentumAtEntry().mag(), (*Hits).energyLoss()/f_dx);
			   h_Hit_StripDets_dedx_Track_P_->Fill((*simTrack).momentum().P(), (*Hits).energyLoss()/f_dx);
			   h_Hit_StripDets_Map_->Fill(theGeomDet->toGlobal((*Hits).localPosition()).z(),
						      theGeomDet->toGlobal((*Hits).localPosition()).y());
			   histos_TEC_dedx_[tTopo->tecRing(detId)-1]->Fill((*Hits).energyLoss()/f_dx);
			   histos_PDG_TEC_dedx_[tTopo->tecRing(detId)-1]->Fill(abs((*Hits).particleType()), (*Hits).energyLoss()/f_dx);
			   histos_PDG_TEC_dx_[tTopo->tecRing(detId)-1]->Fill(abs((*Hits).particleType()), f_dx);
			   histos_p_TEC_dedx_[tTopo->tecRing(detId)-1]->Fill((*Hits).pabs(), (*Hits).energyLoss()/f_dx);
		       }
		   }// if((*Hits).trackId() == (*simTrack).trackId())
	       }// Hits!=simHits.end()
	   }// simHitCollection.isValid()
	} // trackerContainers
	h_Hit_PixelDets_Multiplicity_->Fill(i_SimHit_Pixels_Multiplicity) ;
	h_Hit_StripDets_Multiplicity_->Fill(i_SimHit_Multiplicity-i_SimHit_Pixels_Multiplicity);
	h_Hit_PixelDets_Multiplicity_Track_P_->Fill(i_SimHit_Pixels_Multiplicity, (*simTrack).momentum().P());
	h_Hit_StripDets_Multiplicity_Track_P_->Fill(i_SimHit_Multiplicity-i_SimHit_Pixels_Multiplicity, (*simTrack).momentum().P());
    }// simTrack != simTracks.end()
    h_Track_Multiplicity_->Fill(simTracks.size());
     */    
#ifdef THIS_IS_AN_EVENT_EXAMPLE
    Handle<ExampleData> pIn;
    iEvent.getByLabel("example",pIn);
#endif
    
#ifdef THIS_IS_AN_EVENTSETUP_EXAMPLE
    ESHandle<SetupData> pSetup;
    iSetup.get<SetupRecord>().get(pSetup);
#endif
}

/*
// ------------ method called once each job just before starting event loop  ------------
void 
DemoAnalyzer::beginJob()
{
}
*/
/*
// ------------ method called once each job just after ending the event loop  ------------
void 
DemoAnalyzer::endJob() 
{
}
*/
// ------------ method called when starting to processes a run  ------------

void 
DemoAnalyzer::dqmBeginRun(edm::Run const&, edm::EventSetup const&)
{
}


// ------------ method called when ending the processing of a run  ------------
/*
void 
DemoAnalyzer::endRun(edm::Run const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when starting to processes a luminosity block  ------------
/*
void 
DemoAnalyzer::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when ending the processing of a luminosity block  ------------
/*
void 
DemoAnalyzer::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
DemoAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) 
{
    //The following says we do not know what parameters are allowed so do no validation
    // Please change this to state exactly what you do use, even if it is no parameters
    edm::ParameterSetDescription desc;
    desc.setUnknown();
    descriptions.addDefault(desc);
}

/*
void 
DemoAnalyzer::bookHistosPerDetector()
{
    //EnergyLosses1D
    int    i_nbins = 200;
    float  f_range = 1e-2;
    //EnergyLosses2D-PDG binning
    int    i_nbins_PDG = 2*2212;
    float  f_range_PDG = 2213;
    //EnergyLosses2D-momentum(p) binning
    int    i_nbins_p = 50;
    float  f_range_p = 60;
    //EnergyLosses2D-specific energy loss(dedx) binning
    int    i_nbins_dedx = 200;    
    float  f_range_dedx = 1e-2; //1e-3 (eloss)
    //EnergyLosses2D-thickness(dx) binning
    int    i_nbins_dx = 200; 
    float  f_range_dx = 1e-1;
    
    //Detectors
    //Pixel Barrel - 3 different detector
    bookEnergyLosses1D( histos_PXB_dedx_,  i_nbins, f_range, "PXB", nHistos_PXB_);
    bookEnergyLossesRelatedInfo2D( histos_PDG_PXB_dedx_, i_nbins_PDG, f_range_PDG, i_nbins_dedx, f_range_dedx,"dEdx", "PDG_PXB", nHistos_PXB_ );
    bookEnergyLossesRelatedInfo2D( histos_PDG_PXB_dx_, i_nbins_PDG, f_range_PDG, i_nbins_dx, f_range_dx, "dx", "PDG_PXB", nHistos_PXB_ );
    bookEnergyLossesRelatedInfo2D( histos_p_PXB_dedx_, i_nbins_p, f_range_p, i_nbins_dedx, f_range_dedx, "dEdx", "p_PXB", nHistos_PXB_ );
    //Pixel Endcap - 2 different detector
    bookEnergyLosses1D(histos_PXF_dedx_,  i_nbins, f_range, "PXF", nHistos_PXF_);
    bookEnergyLossesRelatedInfo2D( histos_PDG_PXF_dedx_, i_nbins_PDG, f_range_PDG, i_nbins_dedx, f_range_dedx,"dEdx", "PDG_PXF", nHistos_PXF_ );
    bookEnergyLossesRelatedInfo2D( histos_PDG_PXF_dx_, i_nbins_PDG, f_range_PDG, i_nbins_dx, f_range_dx,"dx", "PDG_PXF", nHistos_PXF_ );
    bookEnergyLossesRelatedInfo2D( histos_p_PXF_dedx_, i_nbins_p, f_range_p, i_nbins_dedx, f_range_dedx, "dEdx", "p_PXF", nHistos_PXF_ );
    //TIB - 4 different detectors
    bookEnergyLosses1D( histos_TIB_dedx_, i_nbins, f_range, "TIB", nHistos_TIB_ );
    bookEnergyLossesRelatedInfo2D( histos_PDG_TIB_dedx_, i_nbins_PDG, f_range_PDG, i_nbins_dedx, f_range_dedx,"dEdx", "PDG_TIB", nHistos_TIB_ );
    bookEnergyLossesRelatedInfo2D( histos_PDG_TIB_dx_, i_nbins_PDG, f_range_PDG, i_nbins_dx, f_range_dx,"dx", "PDG_TIB", nHistos_TIB_ );
    bookEnergyLossesRelatedInfo2D( histos_p_TIB_dedx_, i_nbins_p, f_range_p, i_nbins_dedx, f_range_dedx, "dEdx", "p_TIB", nHistos_TIB_ );
    //TOB - 6 different detectors
    bookEnergyLosses1D( histos_TOB_dedx_, i_nbins, f_range, "TOB", nHistos_TOB_ );
    bookEnergyLossesRelatedInfo2D( histos_PDG_TOB_dedx_, i_nbins_PDG, f_range_PDG, i_nbins_dedx, f_range_dedx, "dEdx", "PDG_TOB", nHistos_TOB_ );
    bookEnergyLossesRelatedInfo2D( histos_PDG_TOB_dx_, i_nbins_PDG, f_range_PDG, i_nbins_dx, f_range_dx, "dx", "PDG_TOB", nHistos_TOB_ );
    bookEnergyLossesRelatedInfo2D( histos_p_TOB_dedx_, i_nbins_p, f_range_p, i_nbins_dedx, f_range_dedx, "dEdx", "p_TOB", nHistos_TOB_ );
    //TID - 3 different detectors
    bookEnergyLosses1D( histos_TID_dedx_, i_nbins, f_range, "TID", nHistos_TID_ );
    bookEnergyLossesRelatedInfo2D( histos_PDG_TID_dedx_, i_nbins_PDG, f_range_PDG, i_nbins_dedx, f_range_dedx,"dEdx", "PDG_TID", nHistos_TID_ );
    bookEnergyLossesRelatedInfo2D( histos_PDG_TID_dx_, i_nbins_PDG, f_range_PDG, i_nbins_dx, f_range_dx,"dx", "PDG_TID", nHistos_TID_ );
    bookEnergyLossesRelatedInfo2D( histos_p_TID_dedx_, i_nbins_p, f_range_p, i_nbins_dedx, f_range_dedx, "dEdx", "p_TID", nHistos_TID_ );
    //TEC - 7 different detectors
    bookEnergyLosses1D( histos_TEC_dedx_, i_nbins, f_range, "TEC", nHistos_TEC_ );
    bookEnergyLossesRelatedInfo2D( histos_PDG_TEC_dedx_, i_nbins_PDG, f_range_PDG, i_nbins_dedx, f_range_dedx,"dEdx", "PDG_TEC", nHistos_TEC_ );
    bookEnergyLossesRelatedInfo2D( histos_PDG_TEC_dx_, i_nbins_PDG, f_range_PDG, i_nbins_dx, f_range_dx,"dx", "PDG_TEC", nHistos_TEC_ );
    bookEnergyLossesRelatedInfo2D( histos_p_TEC_dedx_, i_nbins_p, f_range_p, i_nbins_dedx, f_range_dedx, "dEdx", "p_TEC", nHistos_TEC_ );

}
*/

void 
DemoAnalyzer::bookHistosPerParticle(const std::string& str_particles,  DQMStore::IBooker & ibooker)
{


    //EnergyLosses2D-momentum(p) binning
    int    i_nbins_p = 50;
    float  f_range_p = 60;
    //EnergyLosses2D-specific energy loss(dedx) binning
    int    i_nbins_dedx = 200;    
    float  f_range_dedx = 1e-2; //1e-3 (eloss)

    /*    
    //EnergyLosses2D-PDG binning
    int    i_nbins_PDG = 2*2212;
    float  f_range_PDG = 2213;
    //EnergyLosses2D-momentum(p) binning
    int    i_nbins_p = 50;
    float  f_range_p = 60;
    //EnergyLosses2D-specific energy loss(dedx) binning
    int    i_nbins_dedx = 200;    
    float  f_range_dedx = 1e-2; //1e-3 (eloss)
    //EnergyLosses2D-thickness(dx) binning
    int    i_nbins_dx = 200; 
    float  f_range_dx = 1e-1;
    */
    //
    std::vector<std::string> vctr_particles;
    std::set<char> token;
    identifyToken(token, str_particles);
    std::istringstream particles_toDraw(str_particles);
    std::string s_particle;
    while (std::getline(particles_toDraw, s_particle, (*token.begin()))) {
      removeWhiteSpaces(s_particle);
      vctr_particles.push_back(s_particle);
    }
  
    //    std::cout<<vctr_particles.size()<<std::endl;
    //
    MyClassSetVector histos_particles_p_dedx;
    //MyClassSetMap my_map;
    

    for (unsigned int i=0; i<vctr_particles.size(); ++i) 
    {
      //bool b_initialize_EnergyLossesRelatedInfo2D = false;
	//bool b_store_protons_EnergyLossesRelatedInfo2D= false;
	//bool b_store_pions_EnergyLossesRelatedInfo2D= false;
        for ( auto& x: map_subdet_nlayers_)//=map_subdet_nlayers_.begin(); it!=map_subdet_nlayers_.end(); ++it)
	{
	    
	    std::vector<MonitorElement*> histos_PXB_dedx;
	    std::vector<MonitorElement*> histos_PDG_PXF_dedx;

	    std::string str_subdet; 
	    str_subdet.append(x.first);

	    //if(!b_initialize_EnergyLossesRelatedInfo2D) {
	    bookEnergyLossesRelatedInfo2D( histos_PDG_PXF_dedx, ibooker, i_nbins_p, f_range_p, i_nbins_dedx, f_range_dedx, "p", "dEdx", Form("%s_%s", str_subdet.data(), vctr_particles[i].data()), 1 );
	      // }

	    bookEnergyLosses1D( histos_PXB_dedx, ibooker, i_nbins_dedx, f_range_dedx, "dEdx", Form("%s_%s", str_subdet.data(), vctr_particles[i].data()), x.second );
	    my_map1[vctr_particles[i]].push_back(histos_PXB_dedx);
	    if (vctr_particles[i].compare("protons")==0)
	    {
	      //if(!b_initialize_EnergyLossesRelatedInfo2D){
	      histos_protons_p_dedx_.push_back(histos_PDG_PXF_dedx);
	      histos_particles_p_dedx.push_back(histos_PDG_PXF_dedx);
	      
	      my_map2[vctr_particles[i]].push_back(histos_PDG_PXF_dedx);
		//b_initialize_EnergyLossesRelatedInfo2D = true;
		//	      }
	      histos_protons_dedx_.push_back(histos_PXB_dedx);

	    }
	    else if (vctr_particles[i].compare("pions")==0)
	    {
	      //if(!b_initialize_EnergyLossesRelatedInfo2D){
		histos_pions_p_dedx_.push_back(histos_PDG_PXF_dedx);
		//b_initialize_EnergyLossesRelatedInfo2D = true;
		// }
	      histos_pions_dedx_.push_back(histos_PXB_dedx);
	    }
	}
    }
    std::cout<<my_map1.size()<<std::endl;

    /*It works!
    typedef std::map<std::string,  std::vector<std::vector<MonitorElement*> > >::iterator iter;
    
    for( iter iterator = my_map.begin(); iterator != my_map.end(); ++iterator ) {
      MyClassSetVector foo = iterator->second;
      
      //std::string Key = iter.first;
      //iterator->first = key
      std::cout<<" "<<iterator->first<<std::endl;    
      
      for (unsigned i =0; i<foo.size(); i++) {
	foo[0][0]->Fill(0);
	//std::cout <<  " " << foo[0][0]->Fill(0) << std::endl;
      }
    }
    */
    
    //my_map["protons"]
}


void 
DemoAnalyzer::bookEnergyLosses1D( std::vector<MonitorElement*>& histos_det_dedx, DQMStore::IBooker & ibooker, int nBins, float range, 
				  const TString &var, const TString &det, unsigned int nHistos )
{
    if (!bins_p.empty()) {

    for(unsigned int iHist = 0; iHist < nHistos; iHist++) 
      histos_det_dedx.push_back( ibooker.book1D( Form( "SimHit_%s_%s_%u", var.Data(), det.Data(), iHist+1 ) ,
						 Form( "(%s,%s);%s;", bins_p[iHist].data(), bins_p[iHist+1].data(), var.Data() ) ,
						 nBins , 0. , range ) );
    }
}


void 
DemoAnalyzer::bookEnergyLossesRelatedInfo2D( std::vector<MonitorElement*>& histos_PDG_det_dedx, DQMStore::IBooker & ibooker,
					     int nBinsX, float rangeX, int nBinsY, float rangeY, 
					     const TString &varX, const TString &varY, const TString &det, unsigned int nHistos )
{
    //TString x_Axis_Label  = ( (TObjString*) (det.Tokenize("_")->At(0)))->GetString().Data();
    for(unsigned int iHist = 0; iHist < nHistos; iHist++) 
      histos_PDG_det_dedx.push_back( ibooker.book2D( Form( "SimHit_%s_%s_%s_%u", varY.Data(), varX.Data(), det.Data(), iHist+1 ) ,
						     Form( ";%s;%s", varX.Data(), varY.Data() ) ,
						     nBinsX , 0. , rangeX, nBinsY , 0.0 ,  rangeY ) );

}


/*
void 
DemoAnalyzer::bookEnergyLossesRelatedInfo2D( std::vector<TH2F*>& histos_PDG_det_dedx, 
					     int nBinsX, float rangeX, int nBinsY, float rangeY, 
					     const TString &var, const TString &det )
{
    TString x_Axis_Label  = ( (TObjString*) (det.Tokenize("_")->At(0)))->GetString().Data();
    int iHist = 0;
    histos_PDG_det_dedx.push_back( fs_->make<TH2F>(Form( "hist_%s_%u_%s" , det.Data() , iHist+1, var.Data() ) ,
						   Form( "Sim_Hit_%s_%s_%u;%s;" , var.Data(), det.Data() , iHist+1, x_Axis_Label.Data() ) ,
						   nBinsX , 0. , rangeX, nBinsY , 0.0 ,  rangeY) );
    
}

*/

void DemoAnalyzer::bookHistograms(DQMStore::IBooker & ibooker, edm::Run const & iRun,
				  edm::EventSetup const & iSetup)
{

     ibooker.cd();
     ibooker.setCurrentFolder("testMaterialEffects");
     //map_subdet_nlayers_.insert(std::pair<std::string,int>("PXF",2));
     //map_subdet_nlayers_.insert(std::pair<std::string,int>("PXB",3));
     map_subdet_nlayers_.insert(std::pair<std::string,int>("Pixels",bins1_.size()-1));
     map_subdet_nlayers_.insert(std::pair<std::string,int>("StripsNoTEC5to7",bins_p.size()-1));
     map_subdet_nlayers_.insert(std::pair<std::string,int>("StripsOnlyTEC5to7",bins_p.size()-1));


     //bookHistosPerParticle("protons:pions", ibooker);
     if (!particles_.empty()) {
     bookHistosPerParticle(particles_, ibooker);
     }
          
}

void 
DemoAnalyzer::identifyToken (std::set<char> & token, const std::string& tokenize ) 
{
    //Function to automatically identify token used by the user in the string literal 'tokenize'.
    //Token could be everything other a number, i.e. everything with ASCII code outside the [48,57] range
    //other than a letter, i.e. everything with ASCII code outside the [65,90] and [97,122] range
    //and other than an underscore, i.e. with ASCII code different than 95
    for(unsigned int count =0; count < tokenize.size(); count++)
      if (( (int) tokenize[count]>=48 && (int) tokenize[count]<=57 ) ||
	  ( (int) tokenize[count]>=65 && (int) tokenize[count]<=90 ) || 
	  ( (int) tokenize[count]>=97 && (int) tokenize[count]<=122 ) ||
	  (int) tokenize[count]==95 ) 
      
	continue;
      else
	if ( (int) tokenize[count]!=32 )
	  token.insert(tokenize[count]);

    if (token.size()>1) throw "Please define one delimeter at most!";
}

void 
DemoAnalyzer::removeWhiteSpaces( std::string &strg ) 
{  
    //Function to remove white spaces that accidentally have been added by the user in the string literal
    int i;
    i = strg.find("  ");
  
    if ( i > -1 )
    {
        removeWhiteSpaces (  strg.replace( i,2,"" ) ); //recursive calling of the function itself
    }
    else if ( i == -1 )
    {
        i = strg.find(" ");
	if ( i > -1 )
	{
	    strg.erase(0,1);
	}
    }
}

//define this as a plug-in
DEFINE_FWK_MODULE(DemoAnalyzer);
