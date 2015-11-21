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
					  int nBinsX, float rangeX,  int nBinsY, float rangeY, 
					  const TString &varX, const TString &varY, const TString &det );

      void bookEnergyLossesRelatedInfo1D( std::vector<MonitorElement*>&, DQMStore::IBooker & ibooker, 
					  const std::string &det, const std::string &pdgName, int iHisto );
  
      void bookEnergyLossesRelatedInfo1D( std::vector<MonitorElement*>&, DQMStore::IBooker & ibooker,
					  int iHisto );


      void defaultHistoFormating1D(int iHisto);
      void defaultHistoFormating2D(int iHisto);
      bool passesFilterSelection(edm::SimTrackContainer::const_iterator simTrack);
      void parseConfiguration(const edm::ParameterSet& iConfig);
      void parseHistoFormating1D(const std::vector<edm::ParameterSet>& selectEventPSets);
      void parseHistoFormating2D(const std::vector<edm::ParameterSet>& selectEventPSets);
      void publishParticleFilterWarnings();
      void publishHistoFormating1DWarnings(const std::vector<edm::ParameterSet>& selectEventPSets);
      void publishHistoFormating2DWarnings(const std::vector<edm::ParameterSet>& selectEventPSets);
  
      std::vector<edm::InputTag> simHitsTag_;
      std::vector<double> bins_test_;

      std::vector<PdtEntry> pdts_;   // these are needed before we get the EventSetup
      std::vector<int> pdgIdsToSelect_;
      edm::ESHandle<HepPDT::ParticleDataTable> fPDGTable ;

      //Helper Functions
      void identifyToken(std::set<char> & , const std::string& );
      void removeWhiteSpaces( std::string &strg );

      typedef std::list<std::vector<MonitorElement*> > MyClassSetVector;
      typedef std::map<std::string,  MyClassSetVector > MyClassSetMap;
      MyClassSetMap my_map1;
      MyClassSetMap my_map2;
  
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
	unsigned int i_SimHit_Multiplicity = 0, i_SimHit_Pixels_Multiplicity = 0;
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
		       i_SimHit_Multiplicity++;
	     
		       std::set<unsigned int> detIds;
		       unsigned int detId = (*Hits).detUnitId();
		       unsigned int isub  = DetId(detId).subdetId();
		       const GeomDet* theGeomDet = theTrackerGeometry->idToDet((*Hits).detUnitId());
		       h_Track_Eta_->Fill((*simTrack).momentum().Eta());
		       h_SimHit_Eta_->Fill( 1/2.*std::log((theGeomDet->toGlobal((*Hits).momentumAtEntry()).mag()+theGeomDet->toGlobal((*Hits).momentumAtEntry()).z())/ (theGeomDet->toGlobal((*Hits).momentumAtEntry()).mag()-theGeomDet->toGlobal((*Hits).momentumAtEntry()).z())) ); 
		       h_SimHit_Eta_CM->Fill( 1/2.*std::log(((*Hits).momentumAtEntry().mag()+(*Hits).momentumAtEntry().z())/ ((*Hits).momentumAtEntry().mag()-(*Hits).momentumAtEntry().z())));
		       // It works!
		       //const HepPDT::ParticleData* 
		       //	 PData = fPDGTable->particle(HepPDT::ParticleID(abs(simHitparticleID))) ;
		       //double mass   = PData->mass().value() ;
		       
		       float f_dx = TMath::Power(
						 TMath::Power( (*Hits).entryPoint().x() - (*Hits).exitPoint().x(), 2) + 
						 TMath::Power( (*Hits).entryPoint().y() - (*Hits).exitPoint().y(), 2) + 
						 TMath::Power( (*Hits).entryPoint().z() - (*Hits).exitPoint().z(), 2), 
						 1/2.);
		       
		       float momentumAtEntry = (*Hits).momentumAtEntry().mag();
		       //std::cout<< f_dx<<std::endl;

		       bins_test_.push_back(momentumAtEntry);
		       sort(bins_test_.begin(),bins_test_.end());	
		       it=find(bins_test_.begin(),bins_test_.end(), momentumAtEntry);
		       auto pos = std::distance(bins_test_.begin(), it);
		       
		       for( iter iterator = my_map1.begin(); iterator != my_map1.end(); ++iterator )
		       {
			 MyClassSetVector foo = iterator->second;
			 std::cout<< " foo !! "<<std::endl; 
			 if (isub == static_cast<int>(PixelSubdetector::PixelBarrel) || isub == static_cast<int>(PixelSubdetector::PixelEndcap)) 
			   {
			   i_SimHit_Pixels_Multiplicity++;
			   if (pos!=0 && pos<int(bins_test_.size()-1))
			   {
			     /*
			       for( iter iterator = my_map1.begin(); iterator != my_map1.end(); ++iterator ) 
			       {
				   std::cout<<" bug1 "<<iterator->first<<" xliara1 "<<pdgIdsToSelect_.size()<< std::endl;   
				   if (pdgIdsToSelect_.size()!=0)
				   {
				       if (std::find(pdgIdsToSelect_.begin(), pdgIdsToSelect_.end(), simHitparticleID) != pdgIdsToSelect_.end())
				       {
					 MyClassSetVector foo = iterator->second;
				         std::cout<<" bug1 "<<iterator->first<<" xliara2 "<<momentumAtEntry<<" pos "<< pos << std::endl;
					 foo[0][pos-1]->Fill((*Hits).energyLoss()/f_dx);
					 std::cout<<" bug1 "<<iterator->first<<" xliara3 "<<pdgIdsToSelect_.size()<<" pos1 "<< pos-1 << std::endl;
				       }
				   }
				   else
				   {
				        //foo[0][pos-1]->Fill((*Hits).energyLoss()/f_dx);
				   }
				   //std::cout <<  " bug1 " << foo[0][pos-1]->getBinContent(1) << std::endl;
			       }
			     */if (std::find(pdgIdsToSelect_.begin(), pdgIdsToSelect_.end(), simHitparticleID) != pdgIdsToSelect_.end())
			       {
				 int detct = 0;
				 for (std::list<std::vector<MonitorElement*> >::iterator it1 = foo.begin(); it1 != foo.end(); ++it1){

				   int count =0 ;
				   printf("element out : %d %d %.2f ",count, detct, f_dx);
				   std::vector<MonitorElement*>::iterator it2;
				   for (it2 = (*it1).begin(); it2 != (*it1).end(); ++ it2){
				     if ((pos-1) == count && detct/3<1){
				       //std::printf("element: %s \n",(*it2)->getName());
				       std::cout<<(*it2)->getName()<<" "<<count<<std::endl;
				       (*it2)->Fill((*Hits).energyLoss()/f_dx);
				     }
				     //std::cout<< pos<< it- pos<int(bins_test_.size()-1))
				     count++;
				   }
				   detct++;
				  
				 }
			       }
				 //foo[0][0]->Fill((*Hits).energyLoss()/f_dx);
			   }
			  
			   /*
			   for( iter iterator = my_map2.begin(); iterator != my_map2.end(); ++iterator ) 
			   {
			       MyClassSetVector foo = iterator->second;
			       std::cout<<" bug2 "<<iterator->first<<std::endl;    
			       if (pdgIdsToSelect_.size()!=0)
			       {
				   if (std::find(pdgIdsToSelect_.begin(), pdgIdsToSelect_.end(), simHitparticleID) != pdgIdsToSelect_.end())
				   {
				     //foo[0][0]->Fill((*Hits).pabs(), (*Hits).energyLoss()/f_dx);
				   }
			       }
			       else
			       {
				 //foo[0][0]->Fill((*Hits).pabs(), (*Hits).energyLoss()/f_dx);
			       }
			       //std::cout <<  " bug2 " << foo[0][0]->getBinContent(1) << std::endl;
			   }
			   */
		       }
		       else if (isub == StripSubdetector::TIB || isub == StripSubdetector::TOB || 
				isub == StripSubdetector::TID || isub == StripSubdetector::TEC )
		       {
			   if (isub == StripSubdetector::TEC && tTopo->tecRing(detId)>=5) 
			   {
			       if (pos!=0 && pos<int(bins_test_.size()-1))
			       {
				 /*
				 //std::cout<<momentumAtEntry<<" meta "<<pos<<" meta "<<bins_test_.size()<<" "<<(*Hits).particleType()<<std::endl;
				   for( iter iterator = my_map1.begin(); iterator != my_map1.end(); ++iterator ) 
				   {
				       MyClassSetVector foo = iterator->second;
				       std::cout<<" bug3 "<<iterator->first<<std::endl;    
				       if (pdgIdsToSelect_.size()!=0)
				       {
					   if (std::find(pdgIdsToSelect_.begin(), pdgIdsToSelect_.end(), simHitparticleID) != pdgIdsToSelect_.end())
					   {
					     //foo[2][pos-1]->Fill((*Hits).energyLoss()/f_dx);
					   }
				       }
				       else
				       {
					 //foo[2][pos-1]->Fill((*Hits).energyLoss()/f_dx);
				       }
				       //std::cout <<  " bug3 " << foo[2][pos-1]->getBinContent(1) << std::endl;
				   }
				 */
				 if (std::find(pdgIdsToSelect_.begin(), pdgIdsToSelect_.end(), simHitparticleID) != pdgIdsToSelect_.end())
				 {
				   				 
				   //foo[2][0]->Fill((*Hits).energyLoss()/f_dx);
				   /*
				   for (std::list<std::vector<MonitorElement*> >::iterator it1 = foo.begin(); it1 != foo.end(); ++it1){
				     int count =0 ;
				     printf("element out : %d %.2f ",count, f_dx);
				     std::vector<MonitorElement*>::iterator it2;
				     for (it2 = (*it1).begin(); it2 != (*it1).end(); ++ it2){
				       if ((pos-1) == count){
					 std::printf("element: %d %.2f \n",count, f_dx);
					 (*it2)->Fill((*Hits).energyLoss()/f_dx);
				       }
				       //std::cout<< pos<< it- pos<int(bins_test_.size()-1))
				       count++;
				     }
				   }
				   */
				 }
			       }
			   }
			   
			   else 
			   {
			       if (pos!=0 && pos<int(bins_test_.size()-1))
			       {
				 /*
				   for( iter iterator = my_map1.begin(); iterator != my_map1.end(); ++iterator ) 
				   {
				       MyClassSetVector foo = iterator->second;
				       std::cout<<" bug4 "<<iterator->first<<std::endl;    
				       if (pdgIdsToSelect_.size()!=0)
				       {
					   if (std::find(pdgIdsToSelect_.begin(), pdgIdsToSelect_.end(), simHitparticleID) != pdgIdsToSelect_.end())
					   {
					     //foo[1][pos-1]->Fill((*Hits).energyLoss()/f_dx);
					   }
				       }
				       else
				       {
					 //foo[1][pos-1]->Fill((*Hits).energyLoss()/f_dx);
				       }
				       //std::cout <<  " bug4 " << foo[1][pos-1]->getBinContent(1) << std::endl;
				   }
				 */
				 if (pos==2)
				   std::cout <<  " bug4!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! " <<  std::endl;
				 //foo[1][pos-1]->Fill((*Hits).energyLoss()/f_dx);
			       }
			   }
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
		
	      //if (selectedProcessPaths_[iHisto].title.compare("default")==0 && selectedProcessPaths_[iHisto].name.compare("default")==0) 
	      //    {
		        //bookEnergyLossesRelatedInfo2D( histos_PDG_PXF_dedx, ibooker, i_nbins_p, f_range_p, i_nbins_dedx, f_range_dedx, "p", "dEdx", Form("%s_%s", str_subdet.data(), pdg[i].data()));
		      std::cout<< " edo1 !!!!!!!!!!!!!! "<<" iHisto " <<iHisto<< " det" << x.first.data()<<std::endl;
					    
		        bookEnergyLossesRelatedInfo1D( histos_PXB_dedx, ibooker, str_subdet, pdgNames[iPart], iHisto );
			my_map1[pdgNames[iPart]].push_back(histos_PXB_dedx);
			//my_map2[pdg[i]].push_back(histos_PDG_PXF_dedx);
			//    }
			//else
			//                    {
			//		        std::cout<< " edo2 !!!!!!!!!!!!!! "<<std::endl;
		        //bookEnergyLossesRelatedInfo2D( histos_PDG_PXF_dedx, ibooker, i_nbins_p, f_range_p, i_nbins_dedx, f_range_dedx, "p", "dEdx", Form("%s_%s", str_subdet.data(), pdg[i].data())); 
			//bookEnergyLossesRelatedInfo1D( histos_PXB_dedx, ibooker, iHisto );
			//my_map1[pdgNames[iPart]].push_back(histos_PXB_dedx);
			//my_map2[pdg[i]].push_back(histos_PDG_PXF_dedx); 
			//		    }
		    //}
		    //else
		    //{
		    //bookEnergyLossesRelatedInfo1D( histos_PXB_dedx, ibooker, str_subdet, pdgNames[iPart], iHisto );
		    //my_map1[pdgNames[iPart]].push_back(histos_PXB_dedx);
		    //}
	    }
	}
    }
    //std::cout<<my_map1.size()<<std::endl;

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
     if (name.compare("default")==0 && title.compare("default")==0) 
     {
         histos_det_dedx.push_back( ibooker.book1D( Form( "SimHit_%s_%s_%s_%u", labelx.data(), det.data(), pdgName.data(), iHisto+1 ) ,
						    Form( "(%.2f,%.2f);%s;%s", bins_test_[iHisto] , bins_test_[iHisto+1], labelx.data(),
							  labely.data() ) , 
						    nBins , lowBinX, highBinX ) );
     }
     else
     {
	 histos_det_dedx.push_back( ibooker.book1D( Form( "SimHit_%s_%s_%s_%u", name.data(), det.data(),  pdgName.data(), iHisto+1), 
						    Form( "(%.2f,%.2f)_%s;%s;%s", bins_test_[iHisto] , bins_test_[iHisto+1], 
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
DemoAnalyzer::bookEnergyLossesRelatedInfo2D( std::vector<MonitorElement*>& histos_PDG_det_dedx, DQMStore::IBooker & ibooker,
					     int nBinsX, float rangeX, int nBinsY, float rangeY, 
					     const TString &varX, const TString &varY, const TString &det )
{
    //TString x_Axis_Label  = ( (TObjString*) (det.Tokenize("_")->At(0)))->GetString().Data();
    histos_PDG_det_dedx.push_back( ibooker.book2D( Form( "SimHit_%s_%s_%s", varY.Data(), varX.Data(), det.Data() ) ,
						   Form( ";%s;%s", varX.Data(), varY.Data() ) ,
						   nBinsX , 0. , rangeX, nBinsY , 0.0 ,  rangeY ) );

}


void DemoAnalyzer::bookHistograms(DQMStore::IBooker & ibooker, edm::Run const & iRun,
				  edm::EventSetup const & iSetup)
{

     ibooker.cd();
     ibooker.setCurrentFolder("testMaterialEffects");
     
     map_subdet_nlayers_.insert(std::pair<std::string,int>("Pixels",bins_test_.size()-1));//bins_test_.size()-1));
     map_subdet_nlayers_.insert(std::pair<std::string,int>("StripsNoTEC5to7",bins_test_.size()-1));
     map_subdet_nlayers_.insert(std::pair<std::string,int>("StripsOnlyTEC5to7",bins_test_.size()-1));
     

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
DemoAnalyzer::defaultHistoFormating1D(int iHisto)
{
    //========
    //dynamic title and name (will be substistuted in bookEnergyLossesRelatedInfo1D( std::vector<MonitorElement*>& , DQMStore::IBooker & ibooker, const std::string &, const std::string &, int iHisto) ) 
    //========
    //std::cout<< " size  parseHistoFormatingd "<< iHisto<<" "<<selectedProcessPaths_.size()<<" "<<bins_test_.size()<<std::endl;

    selectedProcessPaths_[iHisto].title= "default"; 
    selectedProcessPaths_[iHisto].name = "default";
    selectedProcessPaths_[iHisto].labelx = "dEdx";
    selectedProcessPaths_[iHisto].labely = "";
    selectedProcessPaths_[iHisto].rangex.push_back(200); //nBins
    selectedProcessPaths_[iHisto].rangex.push_back(0.); //xMin
    selectedProcessPaths_[iHisto].rangex.push_back(0.01); //xMax
    //std::cout<< " size  parseHistoFormatingd "<< iHisto<<" "<<selectedProcessPaths_.size()<<" "<<bins_test_.size()<<std::endl;


}

/*
void
DemoAnalyzer::defaultHistoFormating2D(int iHisto)
{

    selectedProcessPaths2_.labelx = "p";
    selectedProcessPaths2_.labely = "dEdx";
    selectedProcessPaths2_.rangex.push_back(50); //nBins
    selectedProcessPaths2_.rangex.push_back(0.); //xMin
    selectedProcessPaths2_.rangex.push_back(60); //xMax
    selectedProcessPaths2_.rangey.push_back(200); //nBins
    selectedProcessPaths2_.rangey.push_back(0.); //xMin
    selectedProcessPaths2_.rangey.push_back(0.01); //xMax
}
*/

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
   if (std::find(pdgIdsToFilter_.begin(), pdgIdsToFilter_.end(), fabs((*simTrack).type())) != pdgIdsToFilter_.end())
   {
       return false;
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
    if (iConfig.exists("particleTypes")) 
    {
        pdts_ = iConfig.getUntrackedParameter<std::vector<PdtEntry> >("particleTypes");
    }
    
    if (iConfig.exists("bins_p") && iConfig.exists("bins_E"))
    {	
	edm::LogWarning("TooManyBins") << "Only momentum bins will be used";
	bins_test_ = iConfig.getUntrackedParameter <std::vector <double> >("bins_p");
    }
    else if (iConfig.exists("bins_p") && !iConfig.exists("bins_E"))
    {
	bins_test_ = iConfig.getUntrackedParameter <std::vector <double> >("bins_p");
    }
    else if (!iConfig.exists("bins_p") && iConfig.exists("bins_E"))
    {
        bins_test_ = iConfig.getUntrackedParameter <std::vector <double> >("bins_E");
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
	sort(pdgIdsToFilter_.begin(),pdgIdsToFilter_.end());
    }

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
    publishHistoFormating1DWarnings(selectEventPSets);
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
	    std::cout<< " size  parseHistoFormating "<< iHisto<<" "<<selectedProcessPaths_.size()<<" "<<bins_test_.size()<<std::endl;

	    defaultHistoFormating1D(iHisto);
	}
	    
    }
    
}

void
DemoAnalyzer::parseHistoFormating2D(const std::vector<edm::ParameterSet>& selectEventPSets)
{
    for (unsigned int iHisto=0; iHisto< selectEventPSets.size(); ++iHisto)
    {
        selectedProcessPaths2_.push_back(formating2D());
        selectedProcessPaths2_[iHisto].title=selectEventPSets[iHisto].getParameter<std::string>("title");
	selectedProcessPaths2_[iHisto].name=selectEventPSets[iHisto].getParameter<std::string>("name");
	selectedProcessPaths2_[iHisto].labelx=selectEventPSets[iHisto].getUntrackedParameter<std::string>("labelx");
	selectedProcessPaths2_[iHisto].labely=selectEventPSets[iHisto].getUntrackedParameter<std::string>("labely");
	selectedProcessPaths2_[iHisto].rangex=selectEventPSets[iHisto].getUntrackedParameter<std::vector<double> >("rangex");
	selectedProcessPaths2_[iHisto].rangey=selectEventPSets[iHisto].getUntrackedParameter<std::vector<double> >("rangey");
    }

}


void
DemoAnalyzer::publishHistoFormating1DWarnings(const std::vector<edm::ParameterSet>&  selectEventPSets)
{
    if(selectEventPSets.size()>bins_test_.size()+1)
    {
	edm::LogWarning("TooMany1DHistograms") << "Number of reserved 1D Histos will be reduced to number of ! minus one";
    }
}


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
