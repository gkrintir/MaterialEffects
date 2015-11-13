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
	  std::vector<double> rangex;
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
      void bookHistosPerParticle(std::vector<T> & vec,  DQMStore::IBooker & ibooker);

  
      void bookEnergyLossesRelatedInfo2D( std::vector<MonitorElement*>&, DQMStore::IBooker & ibooker,
					  int nBinsX, float rangeX,  int nBinsY, float rangeY, 
					  const TString &varX, const TString &varY, const TString &det );

      void bookEnergyLossesRelatedInfo1D( std::vector<MonitorElement*>&, DQMStore::IBooker & ibooker, int nBins, float range, 
					  const TString &var, const TString &det, unsigned int nHistos );
      
      bool passesFilterSelection(std::vector<PSimHit>::const_iterator Hits);
      void parseHistoFormating1D(const std::vector<edm::ParameterSet>& selectEventPSets);
      void parseHistoFormating2D(const std::vector<edm::ParameterSet>& selectEventPSets);
   
      std::vector<edm::InputTag> simHitsTag_;
      std::string particles_;
      std::vector<std::string> bins_E;
      std::vector<std::string> bins_p;
    
      std::vector<double> bins_test_;
      std::vector<float> bins1_;
      std::vector<PdtEntry> pdts_;   // these are needed before we get the EventSetup
      std::vector<int> pdgIds;
      edm::ESHandle<HepPDT::ParticleDataTable> fPDGTable ;

      //Helper Functions
      void identifyToken(std::set<char> & , const std::string& );
      void removeWhiteSpaces( std::string &strg );

      typedef std::vector<std::vector<MonitorElement*> > MyClassSetVector;
      typedef std::map<std::string,  MyClassSetVector > MyClassSetMap;
      MyClassSetMap my_map1;
      MyClassSetMap my_map2;
  
      typedef StringCutObjectSelector<reco::GenParticle> StrFilter;
      std::auto_ptr<StrFilter> filter_;

      edm::ParameterSet particleFilter_;
      double etaMax_;
      double pTMin_;
      double EMin_;
      formating1D selectedProcessPaths_;
      formating2D selectedProcessPaths2_;
      // Detectors
      std::map<std::string, int> map_subdet_nlayers_;

     
      // PSimHits
      std::vector<edm::InputTag> trackerContainers_;
    
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
  bins_p(),
  selectedProcessPaths_()
  //bins_test(iConfig.getUntrackedParameter<std::vector<double> >("NbEventsToPrint",0.))
  //bins_test(2,static_cast<float>(0.))

  //fs_()


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
    // Possibly allow a string cut
    if (iConfig.existsAs<std::string>("filter")) 
    {
        std::string filter = iConfig.getParameter<std::string>("filter");
        if (!filter.empty())
	{
 	    filter_ = std::auto_ptr<StrFilter>(new StrFilter(filter));
        }
    }

    particles_ = iConfig.getUntrackedParameter <std::string > ("particles");
    bins_E = iConfig.getUntrackedParameter <std::vector <std::string> > ("bins_E", std::vector<std::string>());
    bins_p = iConfig.getUntrackedParameter <std::vector <std::string> > ("bins_p", std::vector<std::string>());

    bins_test_ = iConfig.getUntrackedParameter< std::vector<double> >("DaughterIDs");
    std::cout << " !!!!!!!!!!!!!!!!!!!!!!!!!!!!!! " << bins_test_[0]<< std::endl;
    sort(bins_test_.begin(),bins_test_.end());
    
    //if (!bins_p.empty()) 
    if (iConfig.exists("bins_p")) {
      sort(bins_p.begin(),bins_p.end());	
	
      for (unsigned int it=0; it!=bins_p.size(); ++it)
	bins1_.push_back(std::atof(bins_p[it].data()));
	//std::cout << " !!!!!!!!!!!!!!!!!!!!!!!!!!!!!! " << bins_[it]<< std::endl;
      sort(bins1_.begin(),bins1_.end());	
    }
    
    if (iConfig.exists("TestParticleFilter"))
    {
    	particleFilter_ = iConfig.getParameter<edm::ParameterSet>("TestParticleFilter");
        etaMax_ = particleFilter_.getParameter<double>("etaMax");
	pTMin_ = particleFilter_.getParameter<double>("pTMin");
	EMin_ = particleFilter_.getParameter<double>("EMin");
	std::cout<<" eta!!!! "<<EMin_<<std::endl;
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

    //Retrieve the Monte Carlo truth Particles(SimTracks)
    Handle< SimTrackContainer > simTracksHandle;
    iEvent.getByLabel(trackerContainers_[0].label(),simTracksHandle); 
    const SimTrackContainer simTracks = *(simTracksHandle.product());

    std::vector<double>::iterator it;
    typedef std::map<std::string,  std::vector<std::vector<MonitorElement*> > >::iterator iter;

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
		   if(simHitID == simTrackID && passesFilterSelection(Hits))// &&  
		   {
		       i_SimHit_Multiplicity++;
	     
		       std::set<unsigned int> detIds;
		       unsigned int detId = (*Hits).detUnitId();
		       unsigned int isub  = DetId(detId).subdetId();
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
		       
		       bins_test_.push_back(momentumAtEntry);
		       sort(bins_test_.begin(),bins_test_.end());	
		       it=find(bins_test_.begin(),bins_test_.end(), momentumAtEntry);
		       auto pos = std::distance(bins_test_.begin(), it);
		       
		       if (isub == static_cast<int>(PixelSubdetector::PixelBarrel) || isub == static_cast<int>(PixelSubdetector::PixelEndcap)) 
		       {
			   i_SimHit_Pixels_Multiplicity++;
			   if (pos!=0 && pos<int(bins_test_.size()-1))
			   {
			       for( iter iterator = my_map1.begin(); iterator != my_map1.end(); ++iterator ) 
			       {
				   MyClassSetVector foo = iterator->second;
				   std::cout<<" bug1 "<<iterator->first<<" xliara "<<pdgIds.size()<<std::endl;   
				   if (pdgIds.size()!=0)
				   {
				       if (std::find(pdgIds.begin(), pdgIds.end(), simHitparticleID) != pdgIds.end())
				       {
					   foo[0][pos-1]->Fill((*Hits).energyLoss()/f_dx);
				       }
				   }
				   else
				   {
				       foo[0][pos-1]->Fill((*Hits).energyLoss()/f_dx);
				   }
				   std::cout <<  " bug1 " << foo[0][pos-1]->getBinContent(1) << std::endl;
			       }
			       
			   }
			   for( iter iterator = my_map2.begin(); iterator != my_map2.end(); ++iterator ) 
			   {
			       MyClassSetVector foo = iterator->second;
			       std::cout<<" bug2 "<<iterator->first<<std::endl;    
			       if (pdgIds.size()!=0)
			       {
				   if (std::find(pdgIds.begin(), pdgIds.end(), simHitparticleID) != pdgIds.end())
				   {
				       foo[0][0]->Fill((*Hits).pabs(), (*Hits).energyLoss()/f_dx);
				   }
			       }
			       else
			       {
				   foo[0][0]->Fill((*Hits).pabs(), (*Hits).energyLoss()/f_dx);
			       }
			       std::cout <<  " bug2 " << foo[0][0]->getBinContent(1) << std::endl;
			     }
			   
		       }
		       else if (isub == StripSubdetector::TIB || isub == StripSubdetector::TOB || 
				isub == StripSubdetector::TID || isub == StripSubdetector::TEC )
		       {
			   if (isub == StripSubdetector::TEC && tTopo->tecRing(detId)>=5) 
			   {
			       if (pos!=0 && pos<int(bins_test_.size()-1))
			       {
				   std::cout<<momentumAtEntry<<" meta "<<pos<<" meta "<<bins_test_.size()<<" "<<(*Hits).particleType()<<std::endl;
				   for( iter iterator = my_map1.begin(); iterator != my_map1.end(); ++iterator ) 
				   {
				       MyClassSetVector foo = iterator->second;
				       std::cout<<" bug3 "<<iterator->first<<std::endl;    
				       if (pdgIds.size()!=0)
				       {
					   if (std::find(pdgIds.begin(), pdgIds.end(), simHitparticleID) != pdgIds.end())
					   {
					       foo[2][pos-1]->Fill((*Hits).energyLoss()/f_dx);
					   }
				       }
				       else
				       {
					   foo[2][pos-1]->Fill((*Hits).energyLoss()/f_dx);
				       }
				       std::cout <<  " bug3 " << foo[2][pos-1]->getBinContent(1) << std::endl;
				   }
			       }
			   }
			   
			   else 
			   {
			       if (pos!=0 && pos<int(bins_test_.size()-1))
			       {
				   for( iter iterator = my_map1.begin(); iterator != my_map1.end(); ++iterator ) 
				   {
				       MyClassSetVector foo = iterator->second;
				       std::cout<<" bug4 "<<iterator->first<<std::endl;    
				       if (pdgIds.size()!=0)
				       {
					   if (std::find(pdgIds.begin(), pdgIds.end(), simHitparticleID) != pdgIds.end())
					   {
					       foo[1][pos-1]->Fill((*Hits).energyLoss()/f_dx);
					   }
				       }
				       else
				       {
					     foo[1][pos-1]->Fill((*Hits).energyLoss()/f_dx);
				       }
				       std::cout <<  " bug4 " << foo[1][pos-1]->getBinContent(1) << std::endl;
				   }
			       }
			   }
		       }
		       
		       bins_test_.erase(std::remove(bins_test_.begin(), bins_test_.end(), momentumAtEntry), bins_test_.end());
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
DemoAnalyzer::bookHistosPerParticle(std::vector<T> & vec,  DQMStore::IBooker & ibooker)
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
    /*
    std::vector<std::string> vctr_particles;
    std::set<char> token;
    identifyToken(token, str_particles);
    std::istringstream particles_toDraw(str_particles);
    std::string s_particle;
    while (std::getline(particles_toDraw, s_particle, (*token.begin()))) {
      removeWhiteSpaces(s_particle);
      vctr_particles.push_back(s_particle);
    }
    */
    //    std::cout<<vctr_particles.size()<<std::endl;
    //
    MyClassSetVector histos_particles_p_dedx;
    //MyClassSetMap my_map;
    

    for (unsigned int i=0; i<vec.size(); ++i) 
      {
        for ( auto& x: map_subdet_nlayers_)//=map_subdet_nlayers_.begin(); it!=map_subdet_nlayers_.end(); ++it)
	{
	    
	    std::vector<MonitorElement*> histos_PXB_dedx;
	    std::vector<MonitorElement*> histos_PDG_PXF_dedx;

	    std::string str_subdet; 
	    str_subdet.append(x.first);
	    std::cout << " !!!!!!!!!!!!!!!!!!!!!!!!!!!!!! " << selectedProcessPaths_.name.data()<< std::endl;

	    bookEnergyLossesRelatedInfo2D( histos_PDG_PXF_dedx, ibooker, i_nbins_p, f_range_p, i_nbins_dedx, f_range_dedx, "p", "dEdx", Form("%s_%s", str_subdet.data(), vec[i].data()));
	    bookEnergyLossesRelatedInfo1D( histos_PXB_dedx, ibooker, i_nbins_dedx, f_range_dedx, "dEdx", Form("%s_%s", str_subdet.data(), vec[i].data()), x.second );
	    my_map1[vec[i]].push_back(histos_PXB_dedx);
	    my_map2[vec[i]].push_back(histos_PDG_PXF_dedx);
	}
    }
    std::cout<<my_map1.size()<<std::endl;

}


void 
DemoAnalyzer::bookEnergyLossesRelatedInfo1D( std::vector<MonitorElement*>& histos_det_dedx, DQMStore::IBooker & ibooker, int nBins, float range, 
					     const TString &var, const TString &det, unsigned int nHistos )
{
    if (!bins_p.empty()) 
    {

        for(unsigned int iHist = 0; iHist < nHistos; iHist++) 
	{
	    histos_det_dedx.push_back( ibooker.book1D( Form( "SimHit_%s_%s_%u", var.Data(), det.Data(), iHist+1 ) ,
						       Form( "(%s,%s);%s;",  bins_p[iHist].data() , bins_p[iHist+1].data(), var.Data() ) ,
						       nBins , 0. , range ) );
	}
    }
}


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
     
     map_subdet_nlayers_.insert(std::pair<std::string,int>("Pixels",bins_test_.size()-1));
     map_subdet_nlayers_.insert(std::pair<std::string,int>("StripsNoTEC5to7",bins_p.size()-1));
     map_subdet_nlayers_.insert(std::pair<std::string,int>("StripsOnlyTEC5to7",bins_p.size()-1));
     
     std::vector<std::string> pdgNames; // these are the ones we really use
     if (!pdts_.empty()) 
     {
         for (std::vector<PdtEntry>::iterator itp = pdts_.begin(), edp = pdts_.end(); itp != edp; ++itp) 
	 {
	     itp->setup(iSetup); // decode string->pdgId and vice-versa
	     pdgNames.push_back(itp->name());
	     pdgIds.push_back(std::abs(itp->pdgId()));
	 }
     }
     else
     {
         pdgNames.push_back("all");
     }

     bookHistosPerParticle(pdgNames, ibooker);
	  
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

bool
DemoAnalyzer::passesFilterSelection(std::vector<PSimHit>::const_iterator Hits)
{
  float f_dx = TMath::Power(
			    TMath::Power( (*Hits).entryPoint().x() - (*Hits).exitPoint().x(), 2) + 
			    TMath::Power( (*Hits).entryPoint().y() - (*Hits).exitPoint().y(), 2) + 
			    TMath::Power( (*Hits).entryPoint().z() - (*Hits).exitPoint().z(), 2), 
			    1/2.);
  
  if ( (*Hits).energyLoss()/f_dx<0.0028 )
    return false;
  else
    return true;
}

void
DemoAnalyzer::parseHistoFormating1D(const std::vector<edm::ParameterSet>& selectEventPSets)
{
    for (unsigned int iname=0; iname< selectEventPSets.size(); ++iname)
    {
        selectedProcessPaths_.title=selectEventPSets[iname].getParameter<std::string>("title");
	selectedProcessPaths_.name=selectEventPSets[iname].getParameter<std::string>("name");
	selectedProcessPaths_.labelx=selectEventPSets[iname].getUntrackedParameter<std::string>("labelx");
	selectedProcessPaths_.labely=selectEventPSets[iname].getUntrackedParameter<std::string>("labely");
	selectedProcessPaths_.rangex=selectEventPSets[iname].getUntrackedParameter<std::vector<double> >("rangex");
	std::cout<< " bins !!!!!!!! "<<selectedProcessPaths_.rangex[1]<<std::endl;
    }
}

void
DemoAnalyzer::parseHistoFormating2D(const std::vector<edm::ParameterSet>& selectEventPSets)
{
    for (unsigned int iname=0; iname< selectEventPSets.size(); ++iname)
    {
        selectedProcessPaths2_.title=selectEventPSets[iname].getParameter<std::string>("title");
	selectedProcessPaths2_.name=selectEventPSets[iname].getParameter<std::string>("name");
	selectedProcessPaths2_.labelx=selectEventPSets[iname].getUntrackedParameter<std::string>("labelx");
	selectedProcessPaths2_.labely=selectEventPSets[iname].getUntrackedParameter<std::string>("labely");
	selectedProcessPaths2_.rangex=selectEventPSets[iname].getUntrackedParameter<std::vector<double> >("rangex");
	selectedProcessPaths2_.rangey=selectEventPSets[iname].getUntrackedParameter<std::vector<double> >("rangey");
    }
}


//define this as a plug-in
DEFINE_FWK_MODULE(DemoAnalyzer);
