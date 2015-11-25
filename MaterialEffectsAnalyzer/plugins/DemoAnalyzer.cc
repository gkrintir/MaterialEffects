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
#include "FWCore/Framework/interface/EDFilter.h"
#include "DQMServices/Core/interface/MonitorElement.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "CommonTools/Utils/interface/StringCutObjectSelector.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/TrackerCommon/interface/TrackerTopology.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
#include "SimDataFormats/TrackingHit/interface/PSimHit.h"
#include "SimDataFormats/Track/interface/SimTrack.h"
#include "SimDataFormats/Track/interface/SimTrackContainer.h"
#include "HepPDT/defs.h"
#include "HepPDT/TableBuilder.hh"
#include "HepPDT/ParticleDataTable.hh"
#include "SimGeneral/HepPDTRecord/interface/PdtEntry.h"
#include "SimGeneral/HepPDTRecord/interface/ParticleDataTable.h"

#include "DataFormats/Math/interface/Vector3D.h"
#include "Geometry/Records/interface/IdealGeometryRecord.h"
#include "Geometry/Records/interface/TrackerDigiGeometryRecord.h"
#include "Geometry/TrackerGeometryBuilder/interface/TrackerGeometry.h"

#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "TH2.h"
#include "TMath.h"
#include "TString.h"
#include <algorithm>
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
//https://github.com/cms-tau-pog/HighPtTau_539/blob/master/PhysicsTools/PatAlgos/plugins/PATGenCandsFromSimTracksProducer.cc


class DemoAnalyzer : public DQMEDAnalyzer {
   public:
      explicit DemoAnalyzer(const edm::ParameterSet&);
      ~DemoAnalyzer();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

      struct formating1D
      {
	  std::string title;
	  std::string name;
	  std::string labelx;
	  std::string labely;
	  std::vector<double>  rangex;
      };
  
      struct formating2D
      {
          std::string title;
          std::string name;
          std::string labelx;
          std::string labely;
          std::vector<double> rangex;
          std::vector<double> rangey;
      };


   private:
  
      virtual void dqmBeginRun(edm::Run const& ,edm::EventSetup const& ) override;
      virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
      void bookHistograms(DQMStore::IBooker &, edm::Run const &, edm::EventSetup const &) override;
  
      // ----------member data ---------------------------
      template <typename T> 
      void bookHistosPerParticle(std::vector<T> & pdg, DQMStore::IBooker & ibooker);

  
      void bookEnergyLossesRelatedInfo2D( std::vector<MonitorElement*>&, DQMStore::IBooker & ibooker,
					  const std::string &det, const std::string &pdgName, int iHisto );

      void bookEnergyLossesRelatedInfo1D( std::vector<MonitorElement*>&, DQMStore::IBooker & ibooker, 
					  const std::string &det, const std::string &pdgName, int iHisto );
  
      void bookEnergyLossesRelatedInfo1D( std::vector<MonitorElement*>&, DQMStore::IBooker & ibooker,
					  int iHisto );

  
      typedef std::list<std::vector<MonitorElement*> > MyClassSetVector;
      typedef std::map<std::string,  MyClassSetVector > MyClassSetMap;
      MyClassSetMap my_map1;
      MyClassSetMap my_map2;
  

      void defaultHistoFormating1D(int iHisto, bool publishedParticleFilterWarning);
      void defaultHistoFormating2D(int iHisto, bool publishedParticleFilterWarnings);
      void fillHistoInfo1D(MyClassSetMap::const_iterator iterator, int subdetPosBegin, int subdetPosEnd, int bin, float dE, float dx);
      void fillHistoInfo2D(MyClassSetMap::const_iterator iterator, int subdetPosBegin, int subdetPosEnd, float p, float dE, float dx);
      std::vector<int> findSubDet (int detId, int subdetId, const TrackerTopology* const tTopo);
      bool passesFilterSelection(edm::SimTrackContainer::const_iterator simTrack);
      void parseConfiguration(const edm::ParameterSet& iConfig);
      void parseHistoFormating1D(const std::vector<edm::ParameterSet>& selectEventPSets);
      void parseHistoFormating2D(const std::vector<edm::ParameterSet>& selectEventPSets);
      void publishParticleFilterWarnings();
      void publishHistoFormating1DWarnings(const std::vector<edm::ParameterSet>& selectEventPSets);
      void publishHistoFormating2DWarnings(const std::vector<edm::ParameterSet>& selectEventPSets);
  
      std::vector<edm::InputTag> simHitsTag_;
      std::vector<double> bins_test_;
      bool usePBins; 
      bool useEBins;    
      edm::ParameterSet particleSelection_;
      std::vector<PdtEntry> pdts_;   // these are needed before we get the EventSetup
      std::vector<int> pdgIdsToSelect_;
      bool selectAntiParticle_;

      edm::ESHandle<HepPDT::ParticleDataTable> fPDGTable ;

      //Helper Functions
      void identifyToken(std::set<char> & , const std::string& );
      void removeWhiteSpaces( std::string &strg );

      edm::ParameterSet particleFilter_;
      double etaMin_;
      double etaMax_;
      double pTMin_;
      double pTMax_;
      double pMin_;
      double pMax_;
      double EMin_;
      double EMax_;
      std::vector<int> pdgIdsToFilter_;
      bool filterAntiParticle_;

      std::vector<formating1D> selectedProcessPaths_;
      std::vector<formating2D> selectedProcessPaths2_;
      // Detectors
      std::map<std::string, int> map_subdet_nlayers_;

     
      // PSimHits
      std::vector<edm::InputTag> trackerContainers_;

      edm::Service<TFileService> fs_;
      TH1F* h_Track_Eta_;
      TH1F* h_SimHit_Eta_;
      TH1F* h_SimHit_Eta_CM;
  
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

  simHitsTag_(iConfig.getParameter<std::vector<edm::InputTag> >("SimHitTags")),
  fs_()
  
{
  //initConfiguration();
    parseConfiguration(iConfig);
    h_Track_Eta_ = fs_->make<TH1F>("SimTrack_Eta", ";#eta", 200, -2.5, 2.5); //ok
    h_SimHit_Eta_ = fs_->make<TH1F>("SimHit_Eta", ";#eta", 200, -2.5, 2.5); //ok
    h_SimHit_Eta_CM = fs_->make<TH1F>("SimHit_Eta_CM", ";#eta", 200, -2.5, 2.5); //ok

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

    // Retrieve tracker topology from geometry
    ESHandle< TrackerTopology > tTopoHandle;
    iSetup.get<IdealGeometryRecord>().get(tTopoHandle);
    const TrackerTopology* const tTopo = tTopoHandle.product();

    ESHandle< TrackerGeometry > TG;
    iSetup.get<TrackerDigiGeometryRecord>().get(TG);
    const TrackerGeometry *theTrackerGeometry = TG.product ();

    //Retrieve the Monte Carlo truth Particles(SimTracks)
    Handle< SimTrackContainer > simTracksHandle;
    iEvent.getByLabel(trackerContainers_[0].label(),simTracksHandle); 
    const SimTrackContainer simTracks = *(simTracksHandle.product());

    std::vector<double>::iterator it;
    typedef std::map<std::string,  std::list<std::vector<MonitorElement*> > >::iterator iter;

    for ( SimTrackContainer::const_iterator  simTrack = simTracks.begin(); simTrack != simTracks.end(); ++simTrack) 
    {
        //unsigned int i_SimHit_Multiplicity = 0, i_SimHit_Pixels_Multiplicity = 0;
	for (unsigned i=0; i<trackerContainers_.size(); ++i) 
	{
           Handle<std::vector<PSimHit> > simHitCollection;
	   iEvent.getByLabel(trackerContainers_[i], simHitCollection);
	   const std::vector<PSimHit>& simHits = *simHitCollection.product();
	   if(simHitCollection.isValid()) 
	   {
	       for( std::vector<PSimHit>::const_iterator Hits=simHits.begin(); Hits!=simHits.end(); ++Hits) 
	       {
		   int simHitID = (*Hits).trackId();
		   int simTrackID = (*simTrack).trackId();
		   int simHitparticleID = (*Hits).particleType();
		   if(simHitID == simTrackID && passesFilterSelection(simTrack))// &&  
		   {
		       //i_SimHit_Multiplicity++;
	     
		       std::set<unsigned int> detIds;
		       unsigned int detId = (*Hits).detUnitId();
		       unsigned int subdetId = DetId(detId).subdetId();
		       const GeomDet* theGeomDet = theTrackerGeometry->idToDet((*Hits).detUnitId());
		       h_Track_Eta_->Fill((*simTrack).momentum().Eta());
		       h_SimHit_Eta_->Fill( 1/2.*std::log((theGeomDet->toGlobal((*Hits).momentumAtEntry()).mag()+theGeomDet->toGlobal((*Hits).momentumAtEntry()).z())/ (theGeomDet->toGlobal((*Hits).momentumAtEntry()).mag()-theGeomDet->toGlobal((*Hits).momentumAtEntry()).z())) ); 
		       h_SimHit_Eta_CM->Fill( 1/2.*std::log(((*Hits).momentumAtEntry().mag()+(*Hits).momentumAtEntry().z())/ ((*Hits).momentumAtEntry().mag()-(*Hits).momentumAtEntry().z())));
		  		       
		       float f_dx = TMath::Power(
						 TMath::Power( (*Hits).entryPoint().x() - (*Hits).exitPoint().x(), 2) + 
						 TMath::Power( (*Hits).entryPoint().y() - (*Hits).exitPoint().y(), 2) + 
						 TMath::Power( (*Hits).entryPoint().z() - (*Hits).exitPoint().z(), 2), 
						 1/2.);
		       
		       float momentumAtEntry = 0;
		       if (usePBins)
		       {
			   momentumAtEntry = (*Hits).momentumAtEntry().mag();
			   bins_test_.push_back(momentumAtEntry);
		       }
		       else
		       {
			   momentumAtEntry = (*Hits).momentumAtEntry().mag();
			   const HepPDT::ParticleData* PData = fPDGTable->particle(HepPDT::ParticleID(abs(simHitparticleID))) ;
			   float mass = PData->mass().value() ;
			   momentumAtEntry+= sqrt(pow(momentumAtEntry,2) + pow(mass,2));
			   bins_test_.push_back(momentumAtEntry);
		       }
		       
		       sort(bins_test_.begin(),bins_test_.end());	
		       it=find(bins_test_.begin(),bins_test_.end(), momentumAtEntry);

		       int pos = std::distance(bins_test_.begin(), it);
		       for( iter iterator = my_map2.begin(); iterator != my_map2.end(); ++iterator )
		       {
			   std::vector<int> subdetPosBeginAndEnd = findSubDet(detId, subdetId, tTopo );
			   fillHistoInfo2D(iterator, subdetPosBeginAndEnd[0], subdetPosBeginAndEnd[1], momentumAtEntry, (*Hits).energyLoss(), f_dx);
		       }
 
		       for( iter iterator = my_map1.begin(); iterator != my_map1.end(); ++iterator )
		       {
		           if (pos!=0 && pos<int(bins_test_.size()-1))
			   {
			       std::vector<int> subdetPosBeginAndEnd = findSubDet(detId, subdetId, tTopo );
			       fillHistoInfo1D(iterator, subdetPosBeginAndEnd[0], subdetPosBeginAndEnd[1], pos, (*Hits).energyLoss(), f_dx);
			   }
			   bins_test_.erase(std::remove(bins_test_.begin(), bins_test_.end(), momentumAtEntry), bins_test_.end());

		       }
		   }
	       }
	   }
	}
    }
    
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
DemoAnalyzer::dqmBeginRun(edm::Run const& r, edm::EventSetup const& es)
{
  es.getData( fPDGTable ) ;
  
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

template <typename T> 
void 
DemoAnalyzer::bookHistosPerParticle(std::vector<T> & pdgNames, DQMStore::IBooker & ibooker)
{
  
  //MyClassSetVector histos_particles_p_dedx;

    std::cout<< " size  parseHistoFormating1 "<< selectedProcessPaths_.size()<<" "<<bins_test_.size()<<std::endl;

    for (unsigned int iPart=0; iPart<pdgNames.size(); ++iPart) 
    {
        for ( auto& x: map_subdet_nlayers_)
	{
	    
	    std::vector<MonitorElement*> histos_PXB_dedx;
	    std::vector<MonitorElement*> histos_PDG_PXF_dedx;

	    std::string str_subdet; 
	    str_subdet.append(x.first);

	    int nHistos = x.second;

	    for(unsigned int iHisto = 0; iHisto < nHistos; iHisto++) 
	    {   
	        std::cout<< " edo1 !!!!!!!!!!!!!! "<<" iHisto " <<iHisto<< " det" << x.first.data()<<std::endl;
		
		bookEnergyLossesRelatedInfo1D( histos_PXB_dedx, ibooker, str_subdet, pdgNames[iPart], iHisto );
		my_map1[pdgNames[iPart]].push_back(histos_PXB_dedx);
		if (iHisto == 0)
		{
		    bookEnergyLossesRelatedInfo2D( histos_PDG_PXF_dedx, ibooker, str_subdet, pdgNames[iPart], iHisto );
		    my_map2[pdgNames[iPart]].push_back(histos_PDG_PXF_dedx);
		}
	    }
	}
    }

}


void 
DemoAnalyzer::bookEnergyLossesRelatedInfo1D( std::vector<MonitorElement*>& histos_det_dedx, DQMStore::IBooker & ibooker, 
					     const std::string &det, const std::string &pdgName, int iHisto )
{
    std::string name(selectedProcessPaths_[iHisto].name.data());
    std::string title(selectedProcessPaths_[iHisto].title.data());
    std::string labelx(selectedProcessPaths_[iHisto].labelx.data());
    std::string labely(selectedProcessPaths_[iHisto].labely.data());
    int nBins = selectedProcessPaths_[iHisto].rangex[0];
    double lowBinX = selectedProcessPaths_[iHisto].rangex[1];
    double highBinX = selectedProcessPaths_[iHisto].rangex[2];
    std::string gevOverC = usePBins ? "GeV/c" : "GeV/c^{2}";
    if (name.compare("default")==0 && title.compare("default")==0) 
    {
         histos_det_dedx.push_back( ibooker.book1D( Form( "SimHit_%s_%s_%s_%u", labelx.data(), det.data(), pdgName.data(), iHisto+1 ) ,
						    Form( "(%.2f,%.2f) %s;%s;%s", bins_test_[iHisto] , bins_test_[iHisto+1], gevOverC.data(),
							  labelx.data(), labely.data() ) , 
						    nBins , lowBinX, highBinX ) );
    }
    else
    {
	 histos_det_dedx.push_back( ibooker.book1D( Form( "SimHit_%s_%s_%s_%u", name.data(), det.data(),  pdgName.data(), iHisto+1), 
						    Form( "(%.2f,%.2f) %s_%s;%s;%s", bins_test_[iHisto] , bins_test_[iHisto+1], gevOverC.data(),
							  title.data(), labelx.data(), labely.data() ) ,
						    nBins , lowBinX, highBinX ) );
    }

}
/*
void
DemoAnalyzer::bookEnergyLossesRelatedInfo1D( std::vector<MonitorElement*>& histos_det_dedx, DQMStore::IBooker & ibooker,
                                             int iHisto )
{
  
  std::string name(selectedProcessPaths_[iHisto].name.data());
  std::string title(selectedProcessPaths_[iHisto].title.data());
  std::string labelx(selectedProcessPaths_[iHisto].labelx.data());
  std::string labely(selectedProcessPaths_[iHisto].labely.data());
  int nBins = selectedProcessPaths_[iHisto].rangex[0];
  double lowBinX = selectedProcessPaths_[iHisto].rangex[1];
  double highBinX = selectedProcessPaths_[iHisto].rangex[2];
  histos_det_dedx.push_back( ibooker.book1D( Form( "%s_%s_%u", det.data, name.data(), iHisto+1), 
					     Form( "(%.2f,%.2f)_%s;%s;%s", bins_test_[iHisto] , bins_test_[iHisto+1], 
						   title.data(), labelx.data(), labely.data() ) ,
					     nBins , lowBinX, highBinX ) );
  //histos_det_dedx.push_back( ibooker.book1D( name.data(), title.data(), 100 , 1.,2. ) );
}
*/

void 
DemoAnalyzer::bookEnergyLossesRelatedInfo2D( std::vector<MonitorElement*>& histos_det_dedx, DQMStore::IBooker & ibooker,
                                             const std::string &det, const std::string &pdgName, int iHisto)

{
    std::string name(selectedProcessPaths2_[iHisto].name.data());
    std::string title(selectedProcessPaths2_[iHisto].title.data());
    std::string labelx(selectedProcessPaths2_[iHisto].labelx.data());
    std::string labely(selectedProcessPaths2_[iHisto].labely.data());
    int nBinsX = selectedProcessPaths2_[iHisto].rangex[0];
    double lowBinX = selectedProcessPaths2_[iHisto].rangex[1];
    double highBinX = selectedProcessPaths2_[iHisto].rangex[2];
    int nBinsY = selectedProcessPaths2_[iHisto].rangey[0];
    double lowBinY = selectedProcessPaths2_[iHisto].rangey[1];
    double highBinY = selectedProcessPaths2_[iHisto].rangey[2];
    std::string gevOverC = usePBins ? "GeV/c" : "GeV/c^{2}";
    if (name.compare("default")==0 && title.compare("default")==0)
    {
      std::cout<< " edo2 !!!!!!!!!!!!!! "<<" iHisto " <<iHisto<< " det" <<" "<< labelx.data()<< " "<<  det.data() <<" "<< pdgName.data()<<std::endl;

      histos_det_dedx.push_back( ibooker.book2D( Form( "SimHit_%s_vs_%s_%s_%s_%u", labelx.data(), labely.data(), det.data(), pdgName.data(), iHisto+1) ,
						 Form( ";%s (%s);%s", labelx.data(), gevOverC.data(), labely.data() ) ,
						 nBinsX , lowBinX, highBinX,  nBinsY , lowBinY, highBinY) );
    }
    else
    {
        std::cout<< " edo3 !!!!!!!!!!!!!! "<<" iHisto " <<iHisto<< " det" << std::endl;
	histos_det_dedx.push_back( ibooker.book2D( Form( "SimHit_%s_%s_%s", name.data(), det.data(),  pdgName.data() ),
						   Form( "%s;%s;%s", title.data(), labelx.data(), labely.data() ) ,
						   nBinsX , lowBinX, highBinX, nBinsY , lowBinY, highBinY ) );
    }

}


void DemoAnalyzer::bookHistograms(DQMStore::IBooker & ibooker, edm::Run const & iRun,
				  edm::EventSetup const & iSetup)
{

     ibooker.cd();
     ibooker.setCurrentFolder("testMaterialEffects");
     
     std::vector<std::string> pdgNames; // these are the ones we really use
     if (!pdts_.empty()) 
     {
         for (std::vector<PdtEntry>::iterator itp = pdts_.begin(), edp = pdts_.end(); itp != edp; ++itp) 
	 {
	     itp->setup(iSetup); // decode string->pdgId and vice-versa
	     pdgNames.push_back(itp->name());
	     pdgIdsToSelect_.push_back(std::abs(itp->pdgId()));
	 }
     }
     else
     {
         pdgNames.push_back("all");
     }
     
     sort(pdgIdsToSelect_.begin(),pdgIdsToSelect_.end());
     publishParticleFilterWarnings();
     bookHistosPerParticle(pdgNames, ibooker);
	  
}


void
DemoAnalyzer::defaultHistoFormating1D(int iHisto, bool publishedParticleFilterWarnings=false)
{
    //========
    //dynamic title and name (will be substistuted in bookEnergyLossesRelatedInfo1D( std::vector<MonitorElement*>& , DQMStore::IBooker & , const std::string &, const std::string &, int ) ) 
    //========
    if(!publishedParticleFilterWarnings)
    {
	selectedProcessPaths_[iHisto].title= "default"; 
	selectedProcessPaths_[iHisto].name = "default";
	selectedProcessPaths_[iHisto].labelx = "dEdx";
	selectedProcessPaths_[iHisto].labely = "";
    }
    selectedProcessPaths_[iHisto].rangex.push_back(200); //nBins
    selectedProcessPaths_[iHisto].rangex.push_back(0.); //xMin
    selectedProcessPaths_[iHisto].rangex.push_back(0.01); //xMax
        
}


void
DemoAnalyzer::defaultHistoFormating2D(int iHisto, bool publishedParticleFilterWarnings=false)
{
    //========
    //dynamic title and name (will be substistuted in bookEnergyLossesRelatedInfo1D( std::vector<MonitorElement*>& , DQMStore::IBooker  & , const std::string &, const std::string &, int ) )
    //========
    selectedProcessPaths2_[iHisto].title= "default";
    selectedProcessPaths2_[iHisto].name = "default";
    selectedProcessPaths2_[iHisto].labelx = "p";
    selectedProcessPaths2_[iHisto].labely = "dEdx";
    selectedProcessPaths2_[iHisto].rangex.push_back(50); //nBinsX
    selectedProcessPaths2_[iHisto].rangex.push_back(0.); //xMin
    selectedProcessPaths2_[iHisto].rangex.push_back(60); //xMax
    selectedProcessPaths2_[iHisto].rangey.push_back(200); //nBinsY
    selectedProcessPaths2_[iHisto].rangey.push_back(0.); //yMin
    selectedProcessPaths2_[iHisto].rangey.push_back(0.01); //yMax
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

/*
void 
DemoAnalyzer::fillHistoInfo1D(MyClassSetMap::const_iterator iterator, int subdetPosBegin, int subdetPosEnd, int pos, 
			      int simHitparticleID, float dE, float dx)
{

    MyClassSetVector foo = iterator->second;
			 
    if (std::find(pdgIdsToSelect_.begin(), pdgIdsToSelect_.end(), simHitparticleID) != pdgIdsToSelect_.end())
      {
	int detct = 0;
	for (std::list<std::vector<MonitorElement*> >::iterator it1 = foo.begin(); it1 != foo.end(); ++it1){
	  
	  int count =0 ;
	  printf("element out : %d %d %.2f ",count, detct, dx);
	  std::vector<MonitorElement*>::iterator it2;
	  for (it2 = (*it1).begin(); it2 != (*it1).end(); ++ it2){
	    if ((pos-1) == count && (detct/int(bins_test_.size()-2))>=subdetPosBegin && (detct/int(bins_test_.size()-2))<subdetPosEnd)
	    //if ((pos-1) == count &&  (detct/int(bins_test_.size()-2))<subdetPosEnd)
	      {
	      //std::printf("element: %s \n",(*it2)->getName());
		std::cout<<(*it2)->getName()<<" "<<count<<" "<<detct<<" "<<int(bins_test_.size()-1) <<" "<<subdetPosEnd<<std::endl;
	      (*it2)->Fill(dE/dx);
	    }
	    //std::cout<< pos<< it- pos<int(bins_test_.size()-1))
	    count++;
	  }
	  detct++;
	  
	}
      }
}
*/


void 
DemoAnalyzer::fillHistoInfo1D(MyClassSetMap::const_iterator iterator,  int subdetPosBegin, int subdetPosEnd, int pos, float dE, float dx)
{

    MyClassSetVector foo = iterator->second;
    int detct = 0;
    for (std::list<std::vector<MonitorElement*> >::iterator it1 = foo.begin(); it1 != foo.end(); ++it1)
    {
      int count =0 ;
      //printf("element out : %d %d %.2f ",count, detct, dx);
      std::vector<MonitorElement*>::iterator it2;
      for (it2 = (*it1).begin(); it2 != (*it1).end(); ++ it2)
      {
	  if ((pos-1) == count && (detct/int(bins_test_.size()-2))>=subdetPosBegin && (detct/int(bins_test_.size()-2))<subdetPosEnd)
	  {
	      //std::printf("element: %s \n",(*it2)->getName());
	    //std::cout<<(*it2)->getName()<<" "<<count<<std::endl;
	      (*it2)->Fill(dE/dx);
	  }
	  //std::cout<< pos<< it- pos<int(bins_test_.size()-1))
	  count++;
      }
      detct++;
    }
}

void
DemoAnalyzer::fillHistoInfo2D(MyClassSetMap::const_iterator iterator,  int subdetPosBegin, int subdetPosEnd, float p, float dE, float dx)
{

  MyClassSetVector foo = iterator->second;
  int detct = 0;
  for (std::list<std::vector<MonitorElement*> >::iterator it1 = foo.begin(); it1 != foo.end(); ++it1)
    {
      int count =0 ;
      std::vector<MonitorElement*>::iterator it2;
      for (it2 = (*it1).begin(); it2 != (*it1).end(); ++ it2)
	{
          if (detct/1>=subdetPosBegin && detct/1<subdetPosEnd)
	    {
              //std::printf("element: %s \n",(*it2)->getName());
	    printf("element out : %d %d %.2f \n",count, detct, dx);
	  std::cout<<(*it2)->getName()<<" "<<detct<<" "<< detct/int(map_subdet_nlayers_.size()) <<std::endl;
	    (*it2)->Fill(p, dE/dx);
	    }
          //std::cout<< pos<< it- pos<int(bins_test_.size()-1))
          count++;
	}
      detct++;
    }
}

std::vector<int> 
DemoAnalyzer::findSubDet (int detId, int subdetId, const TrackerTopology* const tTopo)
{
    std::vector<int> subdetPosBeginAndEnd; 
    if (subdetId == static_cast<int>(PixelSubdetector::PixelBarrel) || subdetId == static_cast<int>(PixelSubdetector::PixelEndcap)) 
    {
	subdetPosBeginAndEnd.push_back(0);
	subdetPosBeginAndEnd.push_back(1);
    }
    else if (subdetId == StripSubdetector::TIB || subdetId == StripSubdetector::TOB || 
	     subdetId == StripSubdetector::TID || subdetId == StripSubdetector::TEC )
    {
	if (subdetId == StripSubdetector::TEC && tTopo->tecRing(detId)>=5) 
	{
	    subdetPosBeginAndEnd.push_back(2);
	    subdetPosBeginAndEnd.push_back(3);
	}
	else 
	{
	    subdetPosBeginAndEnd.push_back(1);
	    subdetPosBeginAndEnd.push_back(2);	    
	}
    }
    
    return subdetPosBeginAndEnd;
}    

bool
DemoAnalyzer::passesFilterSelection(edm::SimTrackContainer::const_iterator simTrack)
{

   if ((*simTrack).momentum().Eta() < etaMin_)
   {
       return false;
   }
   if ((*simTrack).momentum().Eta() > etaMax_)
   {
       return false;
   }
   if ((*simTrack).momentum().Pt() < pTMin_)
   {
       return false;
   }
   if ((*simTrack).momentum().Pt() > pTMax_)
   {
       return false;
   }
   if ((*simTrack).momentum().P() < pMin_)
   {
       return false;
   }
   if ((*simTrack).momentum().E() < EMin_)
   {
       return false;
   }
   int particleIDToFilter = filterAntiParticle_ ? (*simTrack).type() : fabs((*simTrack).type());
   if (std::find(pdgIdsToFilter_.begin(), pdgIdsToFilter_.end(), particleIDToFilter) != pdgIdsToFilter_.end())
   {
       return false;
   }
   if (!pdts_.empty()) 
   {
       int particleIDToSelect = selectAntiParticle_ ? (*simTrack).type() : fabs((*simTrack).type());
       if (std::find(pdgIdsToSelect_.begin(), pdgIdsToSelect_.end(), particleIDToSelect) != pdgIdsToSelect_.end())
       {
	   return true; 
       }
       else
       {
	   return false;
       }
   }

   return true;

}

void
DemoAnalyzer::parseConfiguration(const edm::ParameterSet& iConfig)
{
    for (unsigned int i_subDet=0; i_subDet<simHitsTag_.size(); ++i_subDet)
    {
	trackerContainers_.push_back(simHitsTag_[i_subDet]);
    }

    // Possibly allow a list of particle types
    if (iConfig.exists("TestParticleSelection"))
    {
    	particleSelection_ = iConfig.getParameter<edm::ParameterSet>("TestParticleSelection");
	pdts_ = particleSelection_.getUntrackedParameter<std::vector<PdtEntry> >("particleTypes");
	selectAntiParticle_ =  particleSelection_.getUntrackedParameter<bool>("selectAntiParticle");
    }
    
    if (iConfig.exists("bins_p") && iConfig.exists("bins_E"))
    {	
        edm::LogWarning("TooManyBins") << "Only momentum bins will be used";
	bins_test_ = iConfig.getUntrackedParameter <std::vector <double> >("bins_p");
	usePBins = true;useEBins = false;
    }
    else if (iConfig.exists("bins_p") && !iConfig.exists("bins_E"))
    {
	bins_test_ = iConfig.getUntrackedParameter <std::vector <double> >("bins_p");
	usePBins = true;useEBins = false;
    }
    else if (!iConfig.exists("bins_p") && iConfig.exists("bins_E"))
    {
        bins_test_ = iConfig.getUntrackedParameter <std::vector <double> >("bins_E");
	usePBins = false;useEBins = true;
    }


    sort(bins_test_.begin(),bins_test_.end());
    
    if (iConfig.exists("TestParticleFilter"))
    {
    	particleFilter_ = iConfig.getParameter<edm::ParameterSet>("TestParticleFilter");
	etaMin_ = particleFilter_.getParameter<double>("etaMin");
        etaMax_ = particleFilter_.getParameter<double>("etaMax");
	pTMin_ = particleFilter_.getParameter<double>("pTMin");
	pTMax_ = particleFilter_.getParameter<double>("pTMax");
	pMin_ = particleFilter_.getParameter<double>("pMin");
	pMax_ = particleFilter_.getParameter<double>("pMax");
	EMin_ = particleFilter_.getParameter<double>("EMin");
	EMax_ = particleFilter_.getParameter<double>("EMax");
	pdgIdsToFilter_ = particleFilter_.getParameter<std::vector<int> >("pdgIdsToFilter");
	filterAntiParticle_ =  particleFilter_.getParameter<bool>("filterAntiParticle");
	sort(pdgIdsToFilter_.begin(),pdgIdsToFilter_.end());
    }

    
    map_subdet_nlayers_.insert(std::pair<std::string,int>("Pixels",bins_test_.size()-1));
    map_subdet_nlayers_.insert(std::pair<std::string,int>("StripsNoTEC5to7",bins_test_.size()-1));
    map_subdet_nlayers_.insert(std::pair<std::string,int>("StripsOnlyTEC5to7",bins_test_.size()-1));
    

    if (iConfig.exists("formating1D"))
    {
        const std::vector<edm::ParameterSet>& selectEventPSets = iConfig.getParameter<std::vector<edm::ParameterSet>>("formating1D");
	parseHistoFormating1D(selectEventPSets);
    }
    if (iConfig.exists("formating2D"))
    {
        const std::vector<edm::ParameterSet>& selectEventPSets = iConfig.getParameter<std::vector<edm::ParameterSet>>("formating2D");
	parseHistoFormating2D(selectEventPSets);
    }

}

void
DemoAnalyzer::parseHistoFormating1D(const std::vector<edm::ParameterSet>& selectEventPSets)
{
    for (unsigned int iHisto=0; iHisto< bins_test_.size()-1; ++iHisto)
    {   

        selectedProcessPaths_.push_back(formating1D());
        if(iHisto<selectEventPSets.size()) //Read parameters given in the configuration file 
	{
	    selectedProcessPaths_[iHisto].title=selectEventPSets[iHisto].getParameter<std::string>("title");
	    selectedProcessPaths_[iHisto].name=selectEventPSets[iHisto].getParameter<std::string>("name");
	    selectedProcessPaths_[iHisto].labelx=selectEventPSets[iHisto].getUntrackedParameter<std::string>("labelx");
	    selectedProcessPaths_[iHisto].labely=selectEventPSets[iHisto].getUntrackedParameter<std::string>("labely");
	    selectedProcessPaths_[iHisto].rangex=selectEventPSets[iHisto].getUntrackedParameter<std::vector<double> >("rangex");
	}
	else //Default parameters 
	{
	    defaultHistoFormating1D(iHisto);
	}
	    
    }
    publishHistoFormating1DWarnings(selectEventPSets);
    
}

void
DemoAnalyzer::parseHistoFormating2D(const std::vector<edm::ParameterSet>& selectEventPSets)
{

    for (unsigned int iHisto=0; iHisto< map_subdet_nlayers_.size(); ++iHisto)
    {
        selectedProcessPaths2_.push_back(formating2D());
        if(iHisto<selectEventPSets.size()) //Read parameters given in the configuration file
        {
	    selectedProcessPaths2_[iHisto].title=selectEventPSets[iHisto].getParameter<std::string>("title");
	    selectedProcessPaths2_[iHisto].name=selectEventPSets[iHisto].getParameter<std::string>("name");
	    selectedProcessPaths2_[iHisto].labelx=selectEventPSets[iHisto].getUntrackedParameter<std::string>("labelx");
	    selectedProcessPaths2_[iHisto].labely=selectEventPSets[iHisto].getUntrackedParameter<std::string>("labely");
	    selectedProcessPaths2_[iHisto].rangex=selectEventPSets[iHisto].getUntrackedParameter<std::vector<double> >("rangex");
	    selectedProcessPaths2_[iHisto].rangey=selectEventPSets[iHisto].getUntrackedParameter<std::vector<double> >("rangey");
	}
	else //Default parameters
	{
	    defaultHistoFormating2D(iHisto);
	}
    }
    
}

void
DemoAnalyzer::publishHistoFormating1DWarnings(const std::vector<edm::ParameterSet>&  selectEventPSets)
{
    if(selectEventPSets.size()>bins_test_.size()+1)
    {
	edm::LogWarning("TooMany1DHistograms") << "Number of reserved 1D Histos will be reduced to number of ! minus one";
    }
    for (unsigned int iHisto=0; iHisto< bins_test_.size()-1; ++iHisto)
    {  
        if(iHisto<selectEventPSets.size())
	{
	    if (selectedProcessPaths_[iHisto].rangex.size()!=3)
	      {
		edm::LogWarning("NotProperBinConfiguration1DHistograms") << "Wrong bin configuration in the reserved 1D Histos; will be reduced to default ones";

	    selectedProcessPaths_[iHisto].rangex=selectEventPSets[iHisto].getUntrackedParameter<std::vector<double> >("rangex");
	      }
	}
	else //Default parameters 
	{
	    defaultHistoFormating1D(iHisto, true);
	}
    }
}

/*
void
DemoAnalyzer::publishHistoFormating2DWarnings(const std::vector<edm::ParameterSet>&  selectEventPSets)
{
    if(selectEventPSets.size()>map_subdet_nlayers_.size())   
    {
        edm::LogWarning("TooMany2DHistograms") << "Number of reserved 2D Histos will be reduced to number of !";
    }
}
*/
void
DemoAnalyzer::publishParticleFilterWarnings()
{
    if (etaMin_== etaMax_)
    {
        edm::LogWarning("IdenticalFilterValues") << "Unique (min=max) pseudo-rapidity filter value will be used";
    }
    if (pTMin_== pTMax_)
    {
	if (pTMin_< std::numeric_limits<double>::epsilon())
	{
	    edm::LogWarning("DisabledFilterValues") << "No pT filter value will be applied";
	    pTMin_ = -1.;
	    pTMax_ = std::numeric_limits<double>::max();
	}
	else
	{
	    edm::LogWarning("IdenticalFilterValues") << "Unique (min=max) pT filter value will be used with machine tolerance";
	    pTMin_ = pTMin_ - std::numeric_limits<double>::epsilon();
	    pTMax_ = pTMax_ + std::numeric_limits<double>::epsilon();
	}
    }
    if (pMin_== pMax_)
    {
        if (pMin_< std::numeric_limits<double>::epsilon())
	{
	    edm::LogWarning("DisabledFilterValues") << "No p (momentum) filter value will be applied";
	    pMin_ = -1.;
            pMax_ = std::numeric_limits<double>::max();
	}
        else
	{
	    edm::LogWarning("IdenticalFilterValues") << "Unique p (momentum) filter value will be used";
	    pMin_ = pMin_ - std::numeric_limits<double>::epsilon();
            pMax_ = pMax_ + std::numeric_limits<double>::epsilon();
	}
    }
    if (EMin_== EMax_)
    {
        if (EMin_< std::numeric_limits<double>::epsilon())
        {
	    edm::LogWarning("DisabledFilterValues") << "No E (energy) filter value will be applied";
	    EMin_ = -1.;
            EMax_ = std::numeric_limits<double>::max();
	}
        else
        {
	    edm::LogWarning("IdenticalFilterValues") << "Unique E (energy) filter value will be used";
	    EMin_ = pMin_ - std::numeric_limits<double>::epsilon();
            EMax_ = pMax_ + std::numeric_limits<double>::epsilon();
	}
    }
    
    std::set<int> intersect;
    std::set<int>::iterator it;
    std::set_intersection(pdgIdsToSelect_.begin(), pdgIdsToSelect_.end(), pdgIdsToFilter_.begin(), pdgIdsToFilter_.end(), std::inserter(intersect,intersect.begin()));
    if (intersect.size()!=0)
    {
	edm::LogWarning("IdenticalToFilterAndToSelectPDGValues") << "No PDGfilter can be applied for identical PDGselect values!";
	for (it=intersect.begin(); it!=intersect.end(); ++it)
	{
	    pdgIdsToFilter_.erase(std::remove(pdgIdsToFilter_.begin(), pdgIdsToFilter_.end(), *it), pdgIdsToFilter_.end());
	}
    }
    
}

//define this as a plug-in
DEFINE_FWK_MODULE(DemoAnalyzer);
