#include <EventLoop/Job.h>
#include <EventLoop/StatusCode.h>
#include <EventLoop/Worker.h>
#include <ZanalysispPb/Zreco.h>


// ASG status code check
#include <AsgTools/MessageCheck.h>
#include <TSystem.h>
#include <TFile.h>
#include <iostream>
// Infrastructure includes:
#include "xAODRootAccess/Init.h"
#include "xAODRootAccess/TEvent.h"
#include "xAODRootAccess/tools/Message.h"
// EDM includes:
#include "xAODEventInfo/EventInfo.h"
#include "xAODTracking/TrackParticleContainer.h"
#include "xAODTracking/TrackingPrimitives.h"
#include "xAODTracking/VertexContainer.h"
#include "xAODMuon/MuonContainer.h"
#include "xAODMuon/MuonSegmentContainer.h"
#include "xAODBPhys/BPhysHelper.h"
#include "xAODTrigMuon/L2StandAloneMuonContainer.h"
#include "xAODTrigMuon/L2CombinedMuonContainer.h"
#include "xAODTrigger/MuonRoIContainer.h"
#include "xAODTruth/TruthParticleContainer.h"
#include "xAODTruth/TruthVertexContainer.h"
#include "xAODHIEvent/HIEventShapeContainer.h"
#include "TrigMuonMatching/TrigMuonMatching.h"
#include "PATInterfaces/CorrectionCode.h" // to check the return correction code status of tools
#include "xAODCore/ShallowAuxContainer.h"
#include "xAODCore/ShallowCopy.h"
#include "TrigMuonMatching/TrigMuonMatching.h"
#include <TLorentzVector.h>
//to use the fcal tool
#include "xAODHIEvent/HIEventShapeContainer.h"
#include "xAODForward/ZdcModuleContainer.h"
#include "HIEventUtils/HICentralityTool.h"
#include "CxxUtils/make_unique.h"
#include <memory>
#include <string>
#include "TError.h"
//calorimeter
#include "xAODCaloEvent/CaloClusterContainer.h"

#include "ElectronPhotonSelectorTools/AsgElectronPhotonIsEMSelectorConfigHelper.h"
#include "xAODEgamma/EgammaxAODHelpers.h"

#include "GoodRunsLists/GoodRunsListSelectionTool.h"
#include "xAODHIEvent/HIEventShape.h"
#include "egTriggerHelpers/egTriggerHelpers.h"

using std::string;
using std::unique_ptr;
using CxxUtils::make_unique;
#define CHECK( CONTEXT, EXP )                                  \
{                                                              \
   if (!EXP.isSuccess())                                       \
   Fatal( CONTEXT, "Failed to execute: %s. Exiting.", #EXP);   \
 }

 #define CP_TOOL_CHECK( CONTEXT, EXP )                                                                  \
 {                                                                                                    \
  CP::CorrectionCode correctionCode = EXP;                                                             \
  if (correctionCode==CP::CorrectionCode::Ok)                                                          \
     ;                                                                                                  \
	  else if (correctionCode==CP::CorrectionCode::Error)                                                  \
	    Error( CONTEXT, Form("Error when calling %s", #EXP) );                                             \
	  else if (correctionCode==CP::CorrectionCode::OutOfValidityRange)                                     \
	    Warning( CONTEXT, Form("Out-of-range return code when calling %s", #EXP) );                        \
	  else                                                                                                 \
	    Warning( CONTEXT, Form("Unknown correction code %d when calling %s",(int) correctionCode, #EXP));  \
 }


// this is needed to distribute the algorithm to the workers
ClassImp(Zreco)



Zreco :: Zreco ()
{

 //Nothing is importnat is here

 }



EL::StatusCode Zreco :: setupJob (EL::Job& job)
{

 job.useXAOD ();
 ANA_CHECK(xAOD::Init()); // call before opening first file
 return EL::StatusCode::SUCCESS;
}



EL::StatusCode Zreco :: histInitialize ()
{
    TFile *outputFile = wk()->getOutputFile (outputName);
    tree = new TTree ("tree", "tree");
    tree->SetDirectory (outputFile);


    tree2 = new TTree ("tree2", "tree2");
    tree2->SetDirectory (outputFile);

    tree3 = new TTree ("tree3", "tree3");
    tree3->SetDirectory (outputFile);


    //Run and MC
    tree->Branch("RunNumber",   	            &m_RunNumber);
    tree->Branch("EventNumber", 	            &m_EventNumber);
    tree->Branch("LumiBlock",   	            &m_LumiBlock);
    tree->Branch("mcChannelNumber",             &m_mcChannelNumber);
    tree->Branch("FCal_Et",   		            &m_FCal_Et);

    //Lets start from defining the HLT triggers
    tree->Branch("HLT_mu4",                      &m_HLT_mu4); //1
    tree->Branch("HLT_mu6",                      &m_HLT_mu6); //2
    tree->Branch("HLT_mu8",                      &m_HLT_mu8); //3
    tree->Branch("HLT_mu10",                     &m_HLT_mu10); //4
    tree->Branch("HLT_mu15",   		  	         &m_HLT_mu15); //5
    tree->Branch("HLT_mu15_L1MU10",              &m_HLT_mu15_L1MU10);//6
    tree->Branch("HLT_mu15_L1MU6",               &m_HLT_mu15_L1MU6);//7
    tree->Branch("HLT_mu10_L1MU6",	             &m_HLT_mu10_L1MU6);//8
    tree->Branch("HLT_mb_sptrk",				 &m_HLT_mb_sptrk);

    tree->Branch("HLT_e15_lhloose",		         &m_HLT_e15_lhloose);

    tree->Branch("averageIntPerXing",            &m_averageIntPerXing);//average interactions per bunch crossing
    tree->Branch("actualIntPerXing",			 &m_actualIntPerXing);

    tree->Branch("bcID",                       &m_bcid); //number of muons per event
    tree->Branch("is_good_bcid",				&is_good_bcid);


    tree->Branch("Record",						 &m_record);

    //Lets make the L1 triggers
    tree->Branch("L1_MU4",                       &m_L1_MU4);  //1
    tree->Branch("L1_MU6",                       &m_L1_MU6);  //2
    tree->Branch("L1_MU10",                      &m_L1_MU10); //3
    tree->Branch("L1_MU15",                      &m_L1_MU15); //4

    //Muon kinematics
    tree->Branch("Muon_n",                       &m_Muon_n);
    tree->Branch("Muon_pt",   	     		     &m_Muon_pt);
    tree->Branch("Muon_pt_corr",                 &m_Muon_ptcorr);
    tree->Branch("Muon_eta",   	     		     &m_Muon_eta);
    tree->Branch("Muon_phi",   	     		     &m_Muon_phi);
    tree->Branch("Muon_charge",      		     &m_Muon_charge);
    tree->Branch("Muon_quality",     		     &m_Muon_quality);
    tree->Branch("Muon_type",        		     &m_Muon_type);
    tree->Branch("Muon_passedIDCutsMed",    	 &m_Muon_passedIDCutsMed);
    tree->Branch("Muon_passedIDCutsTight",    	 &m_Muon_passedIDCutsTight);
    tree->Branch("Muon_eLoss",    		         &m_Muon_ELoss);
    tree->Branch("Muon_isTrack",                 &m_Muon_isTrack);
    tree->Branch("isCombined",     		         &M_isCombined);
    tree->Branch("isSegTag",     		         &M_isSegTag);
    tree->Branch("isCaloTag",     		         &M_isCaloTag);
    tree->Branch("isStandAl",                    &M_isStandAl);
    tree->Branch("Muon_isMedium",                &m_Muon_isMedium);
    tree->Branch("Muon_isTight",                 &m_Muon_isTight);

    //ID track muon
    tree->Branch("Muon_id_pt",      		     &m_Muon_id_pt);
    tree->Branch("Muon_id_phi",      		     &m_Muon_id_phi);
    tree->Branch("Muon_id_theta",      		     &m_Muon_id_theta);
    tree->Branch("Muon_id_charge",      	     &m_Muon_id_charge);
    tree->Branch("Muon_id_eta",                  &m_Muon_id_eta);


    //ME muon track
    tree->Branch("Me_pt",      		             &m_me_pt);
    tree->Branch("Me_eta",      		         &m_me_eta);
    tree->Branch("Me_phi",      		         &m_me_phi);
    tree->Branch("Me_theta",      	             &m_me_theta);
    tree->Branch("Me_z0",                        &m_me_z0);
    tree->Branch("Me_d0",                        &m_me_d0);
    tree->Branch("Me_ed0",                       &m_me_ed0);
    tree->Branch("Me_charge",                    &m_me_charge);


    //Triggers which we can study the efficincy which were atched to the tracks
    //using the trig matching tool  --HLT
    tree->Branch("Muon_MatchdR_mu4",    	     &m_Muon_MatchdR_mu4); //1
    tree->Branch("Muon_MatchdR_mu6",    	     &m_Muon_MatchdR_mu6); //2
    tree->Branch("Muon_MatchdR_mu8",    	     &m_Muon_MatchdR_mu8); //3
    tree->Branch("Muon_MatchdR_mu10",    	     &m_Muon_MatchdR_mu10);//4
    tree->Branch("Muon_MatchdR_mu15",  		     &m_Muon_MatchdR_mu15);//5
    tree->Branch("Muon_MatchdR_mu15_MU10",  	 &m_Muon_MatchdR_mu15_MU10);//6
    tree->Branch("Muon_MatchdR_mu10_MU6",  	     &m_Muon_MatchdR_mu10_MU6);//7
    tree->Branch("Muon_MatchdR_mu15_MU6",  	     &m_Muon_MatchdR_mu15_MU6);//8
    tree->Branch("El_MatchdR_mu15",  		    &m_El_MatchdR_mu15);//5

    //Trigger matching for L1
    tree->Branch("Muon_MatchdR_L1MU4",           &m_Muon_MatchdR_L1MU4); //1
    tree->Branch("Muon_MatchdR_L1MU6",           &m_Muon_MatchdR_L1MU6); //2
    tree->Branch("Muon_MatchdR_L1MU10",          &m_Muon_MatchdR_L1MU10);//3
    tree->Branch("Muon_MatchdR_L1MU15",          &m_Muon_MatchdR_L1MU15);//4

    //Muon track container ID
    tree->Branch("Muon_Track_d0",      	      	 &m_Muon_Track_d0);
    tree->Branch("Muon_Track_z0",      	      	 &m_Muon_Track_z0);
    tree->Branch("Muon_Track_nPixelHits",        &m_Muon_Track_nPixelHits);
    tree->Branch("Muon_Track_nPixelDeadSensors", &m_Muon_Track_nPixelDeadSensors);
    tree->Branch("Muon_Track_nSCTHits",          &m_Muon_Track_nSCTHits);
    tree->Branch("Muon_Track_nSCTDeadSensors",   &m_Muon_Track_nSCTDeadSensors);
    tree->Branch("Muon_Track_nTRTHits",          &m_Muon_Track_nTRTHits);
    tree->Branch("Muon_Track_nTRTOutliers",      &m_Muon_Track_nTRTOutliers);
    tree->Branch("Muon_Track_nPixelHoles",       &m_Muon_Track_nPixelHoles);
    tree->Branch("Muon_Track_nSCTHoles",         &m_Muon_Track_nSCTHoles);
    tree->Branch("sumEt",                        &m_sumEt);
    tree->Branch("percentile",                   &percentile);
    tree->Branch("Track_number_pileup",		     &trackn);
    tree->Branch("PileUp_vertices",				     &pile_up_vertices);

    //Vertex position
    tree->Branch("Vertex_PV_x",                     &m_vertex_PV_x);
    tree->Branch("Vertex_PV_y",                     &m_vertex_PV_y);
    tree->Branch("Vertex_PV_z",                     &m_vertex_PV_z);
    tree->Branch("Vertex_PU_x",                     &m_vertex_PU_x);
    tree->Branch("Vertex_PU_y",                     &m_vertex_PU_y);
    tree->Branch("Vertex_PU_z",                     &m_vertex_PU_z);
    tree->Branch("PrimaryVertex",				 &m_PV);
    tree->Branch("El_pt",   	     		     &m_El_pt);
    tree->Branch("El_eta",   	     		     &m_El_eta);
    tree->Branch("El_phi",   	     		     &m_El_phi);
    tree->Branch("El_charge",      		     &m_El_charge);
    tree->Branch("El_isMedium",                &m_El_isMedium);
    tree->Branch("El_isTight",                 &m_El_isTight);
    tree->Branch("El_isLoose",                 &m_El_isLoose);


    tree2->Branch("RunNumber",   	            &m_RunNumber);
    tree2->Branch("EventNumber", 	            &m_EventNumber);
    tree2->Branch("LumiBlock",   	            &m_LumiBlock);
    tree2->Branch("mcChannelNumber",            &m_mcChannelNumber);
    tree2->Branch("sumEt",                      &m_sumEt);
    tree2->Branch("percentile",                 &percentile);
    tree2->Branch("Track_number_pileup",		&trackn);
    tree2->Branch("PileUp_vertices",			&pile_up_vertices);
    tree2->Branch("HLT_e15_lhloose",		    &m_HLT_e15_lhloose);
    tree2->Branch("El_MatchdR_mu15",  		    &m_El_MatchdR_mu15);//
    tree2->Branch("averageIntPerXing",          &m_averageIntPerXing);//average interactions per bunch crossing
    tree2->Branch("bcID",                       &m_bcid); //number of muons per event
    tree2->Branch("is_good_bcid",				&is_good_bcid);
    tree2->Branch("El_pt",   	     		     &m_El_pt);
    tree2->Branch("El_eta",   	     		     &m_El_eta);
    tree2->Branch("El_phi",   	     		     &m_El_phi);
    tree2->Branch("El_charge",      		     &m_El_charge);
    tree2->Branch("El_isMedium",                &m_El_isMedium);
    tree2->Branch("El_isTight",                 &m_El_isTight);
    tree2->Branch("El_isLoose",                 &m_El_isLoose);

    tree3->Branch("RunNumber",   	            &m_RunNumber);
    tree3->Branch("EventNumber", 	            &m_EventNumber);
    tree3->Branch("LumiBlock",   	            &m_LumiBlock);
    tree3->Branch("mcChannelNumber",            &m_mcChannelNumber);
    tree3->Branch("sumEt",                      &m_sumEt);
    tree3->Branch("percentile",                 &percentile);
    tree3->Branch("Track_number_pileup",		&trackn);
    tree3->Branch("PileUp_vertices",			&pile_up_vertices);
    tree3->Branch("HLT_g35_loose",				&m_HLT_g35_loose);

    tree3->Branch("averageIntPerXing",          &m_averageIntPerXing);//average interactions per bunch crossing
    tree3->Branch("bcID",                       &m_bcid); //number of muons per event
    tree3->Branch("is_good_bcid",				&is_good_bcid);

  return EL::StatusCode::SUCCESS;
}



EL::StatusCode Zreco :: fileExecute ()
{
  // Here you do everything that needs to be done exactly once for every
  // single file, e.g. collect a list of all lumi-blocks processed
  return EL::StatusCode::SUCCESS;
}



EL::StatusCode Zreco :: changeInput (bool firstFile)
{
  // Here you do everything you need to do when we change input files,
  // e.g. resetting branch addresses on trees.  If you are using
  // D3PDReader or a similar service this method is not needed.
  return EL::StatusCode::SUCCESS;
}



EL::StatusCode Zreco :: initialize ()
{

    xAOD::TEvent* event = wk()->xaodEvent();
    // as a check, let's see the number of events in our xAOD
    Info("initialize()", "Number of events = %lli", event->getEntries() ); // print long long int
    m_eventCounter = 0;

    //for the fcal tool // ZDCAnalysisTool
    m_zdcTools = new ZDC::ZdcAnalysisTool("ZdcAnalysisTool");
    ANA_CHECK(m_zdcTools->initializeTool());
    // HIPileupTool
    m_hiPileup = new HI::HIPileupTool("PileupTool");
    ANA_CHECK(m_hiPileup->initialize());
    //for the centralicy
    centTool =	new HI::HICentralityTool("CentralityTool");
     ANA_CHECK(	centTool->setProperty("RunSpecies","pPb2016"));
    ANA_CHECK(	centTool->initialize()	);
    // GRL
    m_grl = new GoodRunsListSelectionTool("GoodRunsListSelectionTool");
    const char* GRLFilePath = "$ROOTCOREBIN/data/ZanalysispPb/data16_hip8TeV.periodAllYear_DetStatus-v86-pro20-19_DQDefects-00-02-04_PHYS_HeavyIonP_All_Good.xml";
    const char* fullGRLFilePath = gSystem->ExpandPathName (GRLFilePath);
    std::vector<std::string> vecStringGRL;
    vecStringGRL.push_back(fullGRLFilePath);
    ANA_CHECK(m_grl->setProperty( "GoodRunsListVec", vecStringGRL));
    ANA_CHECK(m_grl->setProperty("PassThrough", false)); // if true (default) will ignore result of GRL and will just pass all events
    ANA_CHECK(m_grl->initialize());
    // TDT -trigger decision tool from my code
    configTool = new TrigConf::xAODConfigTool("xAODConfigTool");
    ANA_CHECK(configTool->initialize());
	ToolHandle<TrigConf::ITrigConfigTool> configHandle(configTool);
    //ANA_CHECK(configHandle->initialize())
    trigDecTool = new Trig::TrigDecisionTool("TrigDecisionTool");
    ANA_CHECK(trigDecTool->setProperty("ConfigTool",configHandle));
    ANA_CHECK(trigDecTool->setProperty("OutputLevel", MSG::INFO));
    ANA_CHECK(trigDecTool->setProperty("TrigDecisionKey","xTrigDecision"));
    ANA_CHECK(trigDecTool->initialize());
    // TriggerMatchingTool - the new one where the matching to the muon and track at low pt is good
    m_tmt = new Trig::MatchingTool("MyMatchingTool");
    ANA_CHECK(m_tmt->initialize());
    // MuonTriggerMatching Tool
    m_trigMuonMatching = new Trig::TrigMuonMatching("TrigMuonMatching");
    ToolHandle<Trig::TrigDecisionTool> m_trigDec(trigDecTool);
    ANA_CHECK(m_trigMuonMatching->setProperty("TriggerTool",m_trigDec));
    ANA_CHECK(m_trigMuonMatching->initialize());
    // MuonSelectionTool
    m_muonSelection1 = new CP::MuonSelectionTool("MuonSelection1");
    m_muonSelection1->msg().setLevel(MSG::ERROR);
    ANA_CHECK(m_muonSelection1->setProperty("MaxEta",2.5));
    ANA_CHECK(m_muonSelection1->setProperty("MuQuality",0));
	m_muonSelection2 = new CP::MuonSelectionTool("MuonSelection2");
    m_muonSelection2->msg().setLevel(MSG::ERROR);
    ANA_CHECK(m_muonSelection2->setProperty("MaxEta",2.5));
    ANA_CHECK(m_muonSelection2->setProperty("MuQuality",1));
    m_muonSelection3 = new CP::MuonSelectionTool("MuonSelection2");
    m_muonSelection3->msg().setLevel(MSG::ERROR);
    ANA_CHECK(m_muonSelection3->setProperty("MaxEta",2.5));
    ANA_CHECK(m_muonSelection3->setProperty("MuQuality",2));
    //for the calo tagged muons
    // MuonSelectionTool
    m_muonSelection = new CP::MuonSelectionTool("MuonSelection");
    m_muonSelection->msg().setLevel(MSG::ERROR);
    //ana check
    ANA_CHECK(m_muonSelection->initialize());
    ANA_CHECK(m_muonSelection1->initialize());
	ANA_CHECK(m_muonSelection2->initialize());
	ANA_CHECK(m_muonSelection3->initialize());
    // MuonCalibrationAndSmearingTool
    m_muonCorr = new CP::MuonCalibrationAndSmearingTool( "muonCorrectionTool" );
    m_muonCorr->msg().setLevel( MSG::INFO );
    //m_muonCorr.setProperty( "Release", "PreRecs" );
    ANA_CHECK(m_muonCorr->initialize());


    //Electron selection  tool Tight,   medium, loose
    //LH selector
    std::string confDir = "ElectronPhotonSelectorTools/offline/mc15_20160512/";  //mc15_20150712

    m_electronLHLooseSelector = new AsgElectronLikelihoodTool( "m_electronLHLooseSelector" );
    CHECK("SetupEgammaTools()", m_electronLHLooseSelector->setProperty("primaryVertexContainer","PrimaryVertices"));
    CHECK("SetupEgammaTools()", m_electronLHLooseSelector->setProperty("ConfigFile",confDir+"ElectronLikelihoodLooseOfflineConfig2016_Smooth.conf"));
    CHECK ("SetupEgammaTools()", m_electronLHLooseSelector->initialize());

    m_electronLHMediumSelector = new AsgElectronLikelihoodTool( "m_electronLHMediumSelector" );
    CHECK("SetupEgammaTools()", m_electronLHMediumSelector->setProperty("primaryVertexContainer","PrimaryVertices"));
    CHECK("SetupEgammaTools()", m_electronLHMediumSelector->setProperty("ConfigFile",confDir+"ElectronLikelihoodMediumOfflineConfig2016_Smooth.conf"));
    CHECK ("SetupEgammaTools()", m_electronLHMediumSelector->initialize());

    m_electronLHTightSelector = new AsgElectronLikelihoodTool( "m_electronLHTightSelector" );
    CHECK("SetupEgammaTools()", m_electronLHTightSelector->setProperty("primaryVertexContainer","PrimaryVertices"));
    CHECK("SetupEgammaTools()", m_electronLHTightSelector->setProperty("ConfigFile",confDir+"ElectronLikelihoodTightOfflineConfig2016_Smooth.conf"));
    CHECK ("SetupEgammaTools()", m_electronLHTightSelector->initialize());

      //matching
    m_match_tool = new Trig::TrigEgammaMatchingTool("TrigEgammaMatching");
    CHECK("SetupEgammaTools()", m_match_tool->initialize());

    return EL::StatusCode::SUCCESS;


}



EL::StatusCode Zreco :: execute ()
{

    xAOD::TEvent* event = wk()->xaodEvent();
    if( (m_eventCounter % 500) ==0 ) Info("execute()", "Event number = %i", m_eventCounter );
    m_eventCounter++;

    // Event information
    const xAOD::EventInfo* eventInfo = 0;
    ANA_CHECK(event->retrieve( eventInfo, "EventInfo"));
    bool isMC = false;
    if (eventInfo->eventType( xAOD::EventInfo::IS_SIMULATION ) ){isMC = true;}


    //Filling in the event info
    m_EventNumber = eventInfo->eventNumber();
    m_RunNumber   = eventInfo->runNumber();
    m_LumiBlock   = eventInfo->lumiBlock();
    m_bcid        = eventInfo->bcid();
    m_averageIntPerXing = eventInfo->averageInteractionsPerCrossing();
    m_actualIntPerXing  = eventInfo->actualInteractionsPerCrossing();

    clearVector();
    //GRL cut
    cout<<"================================================"<<endl;
    cout<<"---------------Checking the GRL ...-------------"<<endl;
    if (!isMC)
       { // it's data!
      	 if (!m_grl->passRunLB(*eventInfo))
     	   {
        	    return EL::StatusCode::SUCCESS; // go to next event
      	   }
        }// end if not MC


     //DAQ errors
    cout<<"=============================================================="<<endl;
    cout<<"----------------Checking the Data quality ...-----------------"<<endl;
     if(!isMC)
        {
          if((eventInfo->errorState(xAOD::EventInfo::LAr)==xAOD::EventInfo::Error ) || (eventInfo->errorState(xAOD::EventInfo::Tile)==xAOD::EventInfo::Error ) ||
             (eventInfo->errorState(xAOD::EventInfo::SCT)==xAOD::EventInfo::Error ) || (eventInfo->isEventFlagBitSet(xAOD::EventInfo::Core, 18) ) )
             return EL::StatusCode::SUCCESS;

        }


     //Checking the BCID number
    is_good_bcid = is_this_good_bcid(m_RunNumber,m_bcid);


     cout<<"=============================================================="<<endl;
    cout<<"------------------------Looking at PV ...---------------------"<<endl;
    //Vertex requirement
    const xAOD::VertexContainer * vertices = 0;
    if ( !event->retrieve( vertices, "PrimaryVertices" ).isSuccess() )
      {
          Error("execute()", "Failed to retrieve VertexContainer container. Exiting." );

      }

    if(vertices->size()<2) return EL::StatusCode::SUCCESS;
    xAOD::VertexContainer::const_iterator vtx_itr = vertices->begin();
    xAOD::VertexContainer::const_iterator vtx_end = vertices->end();
    // find primary vertex
    const xAOD::Vertex* primaryVertex = 0;
    for(;vtx_itr!=vtx_end;++vtx_itr)
       {
          if((*vtx_itr)->vertexType()==xAOD::VxType::PriVtx)
            {
                primaryVertex = (*vtx_itr);

                m_vertex_PV_x .  push_back(primaryVertex->x());
                m_vertex_PV_y.   push_back(primaryVertex->y());
                m_vertex_PV_z.   push_back(primaryVertex->z());

               break;
            }
       }

    cout<<"=============================================================="<<endl;
    cout<<"----------------Checking the pile up vertices ...-----------------"<<endl;
    vtx_itr = vertices->begin();
    const xAOD::Vertex* pileup = 0;
    pile_up_vertices=0;
    trackn = 0;
    for(;vtx_itr!=vtx_end;++vtx_itr)
       {
          if((*vtx_itr)->vertexType()==xAOD::VxType::PileUp)
            {
                pileup = (*vtx_itr);
                if((int)pileup->nTrackParticles()>trackn) trackn = pileup->nTrackParticles();
                if(trackn>6)pile_up_vertices=1;


                m_vertex_PU_x.   push_back(pileup->x());
                m_vertex_PU_y.   push_back(pileup->y());
                m_vertex_PU_z.   push_back(pileup->z());

            }
       }//if the number of tracks>6 then this is a pile up vertex


    //Primary vertex requirement errors
    cout<<"======================================================================="<<endl;
    cout<<"------------------=========Cutting on the PV ...========-----------------"<<endl;
    if (!primaryVertex) return EL::StatusCode::SUCCESS;
    if(fabs(primaryVertex->z())>150.) return EL::StatusCode::SUCCESS;
    m_PV=0;
    if(primaryVertex)m_PV=1;



    //the event triggers
    m_HLT_mu4 			    = false;
    m_HLT_mu6 			    = false;
    m_HLT_mu8 			    = false;
    m_HLT_mu15 			    = false;
    m_HLT_mu15_L1MU10 		= false;
    m_HLT_mu10_L1MU6 		= false;
    m_HLT_mu15_L1MU6 		= false;
    m_L1_MU4                = false;
    m_L1_MU6                = false;
    m_L1_MU10               = false;
    m_L1_MU15               = false;
    m_record                = false;
    m_HLT_mb_sptrk          = false;
    m_HLT_e15_lhloose       = false;
    m_HLT_g35_loose         = false;


    if (trigDecTool->isPassed( "HLT_mu4"   	     	    )) { m_HLT_mu4 			    = true;}
    if (trigDecTool->isPassed( "HLT_mu6"   	       	    )) { m_HLT_mu6 			    = true;}
    if (trigDecTool->isPassed( "HLT_mu8"   	       	    )) { m_HLT_mu8 			    = true;}

    if (trigDecTool->isPassed( "HLT_mu15"   	       	)) { m_HLT_mu15 		    = true; }
    if (trigDecTool->isPassed( "HLT_mu15_L1MU10"        )) { m_HLT_mu15_L1MU10 		= true; }
    if (trigDecTool->isPassed( "HLT_mu10_L1MU6"         )) { m_HLT_mu10_L1MU6 		= true; }
    if (trigDecTool->isPassed( "HLT_mu15_L1MU6"         )) { m_HLT_mu15_L1MU6 		= true; }

    if (trigDecTool->isPassed( "L1_MU4"   	     	    )) { m_L1_MU4 			    = true; }
    if (trigDecTool->isPassed( "L1_MU6"   	     	    )) { m_L1_MU6 			    = true; }
    if (trigDecTool->isPassed( "L1_MU10"   	     	    )) { m_L1_MU10 			    = true; }
    if (trigDecTool->isPassed( "L1_MU15"   	     	    )) { m_L1_MU15			    = true; }
    if (trigDecTool->isPassed( "HLT_e15_lhloose"        )) { m_HLT_e15_lhloose = true; }
    if (trigDecTool->isPassed( "HLT_g35_loose"          )) { m_HLT_g35_loose        = true; }
    if (trigDecTool->isPassed( "HLT_mb_sptrk"           )) {m_HLT_mb_sptrk 		    = true; }



    m_sumEt = 0;
    m_sumEt = centTool->getCentralityEstimator();
    percentile = 0;
    percentile = centTool->getCentralityPercentile();


    cout<<"=================="<<endl;
    cout<<"Creating mstore"<<endl;
    cout<<"=================="<<endl;
    xAOD::TStore* m_store = wk()->xaodStore();
    //Open and work with muon container
    const xAOD::MuonContainer* muonsnotcorr = 0;
    const xAOD::MuonContainer* m_newMuons = 0;
    ANA_CHECK(event->retrieve( muonsnotcorr, "Muons" ));
    int m_Muon_number = 0;
    // create a shallow copy of the muons container
    std::pair< xAOD::MuonContainer*, xAOD::ShallowAuxContainer* > muons_corr = xAOD::shallowCopyContainer( *muonsnotcorr );
    xAOD::MuonContainer::iterator begin = (muons_corr.first)->begin();
    xAOD::MuonContainer::iterator end = (muons_corr.first)->end();

    for(; begin!=end; ++begin)
      {
	    // recalibrate 4-momentum (pt and E are both rescaled, but not caloCluster) of shallow-copy
	    float initial_pt = (*begin)->pt();
	    //if (m_config.recalibrateElectrons) {
	      Info("CalibrateMuons()", "  original muon pt = %.5f GeV", ((*begin)->pt() / 1000.));
	      //m_egammaCalibrationAndSmearingTool->setRandomSeed(m_eventInfo->eventNumber()+100*i++);
	      CP_TOOL_CHECK("CalibrateMuons()", m_muonCorr->applyCorrection(**begin));

	    //}
	    float final_pt = (*begin)->pt();
	    float correction = final_pt/initial_pt;
	    Info("CalibrateMuons()", "  corrected muon pt = %.5f GeV", ((*begin)->pt() / 1000.));
	    (*begin)->auxdata< float >( "ptcorr" ) = correction;

	  }

	 // store shallow-copies of original muons as container CalibMuons in transient event store
	if( !m_store->record(muons_corr.first, "CalibMuons") || !m_store->record(muons_corr.second, "CalibMuonsAux.") )
	{
	    Fatal("CalibrateMuons()", "Failed to store CalibMuons. Exiting.");
	    return EL::StatusCode::FAILURE;
	    //return;
	  }

    //retrieve container of calibrated muons
    if ( !m_store->retrieve( m_newMuons, "CalibMuons" ).isSuccess() )
       { // retrieve arguments: container type, container key
         Fatal("SelectMuons()", "Failed to retrieve CalibMuons container. Exiting." );
         return EL::StatusCode::FAILURE;
        }



    const xAOD::ElectronContainer* m_electrons = 0;
    ANA_CHECK(event->retrieve( m_electrons, "Electrons" ));

    std::pair< xAOD::ElectronContainer*, xAOD::ShallowAuxContainer*> electrons_corr = xAOD::shallowCopyContainer(*m_electrons);
	xAOD::ElectronContainer::iterator begin_el = electrons_corr.first->begin();
	xAOD::ElectronContainer::iterator end_el = electrons_corr.first->end();






   	//Trakcs
   	const xAOD::TrackParticleContainer* Recotracks = 0;
    ANA_CHECK(event->retrieve( Recotracks, "InDetTrackParticles" ));
    //1.  Lets make the muon loop for the MUON VECTORS in order to det the info about the muons and their tracks in general
    //also I have added here the loop corresponding to the ME tracks in order to study efficiency of the ID/////////////////

    xAOD::MuonContainer::const_iterator begin_new = m_newMuons->begin();
    xAOD::MuonContainer::const_iterator end_new = m_newMuons->end();
    bool is_a_good_muon = false;
    for(;begin_new!=end_new;++begin_new)

       {
         bool passTightID  = false;
	     bool passMediumID = false;
	     m_Muon_number++;


	     //the cuts on muons
	     if((*begin_new)->pt()/1000.<20)continue;
	     if(abs((*begin_new)->eta())>2.5)continue;
	     if(m_trigMuonMatching->minDelR((*begin_new),"HLT_mu15",1)>0.01)continue;

	     if (m_muonSelection1->accept(*begin_new)) passTightID = true; //based on eta quality and ID cuts - tight muons
	     if (m_muonSelection2->accept(*begin_new)) passMediumID = true; //medium muons
	     //muon type
	     if(passMediumID==false)continue;

	     bool isCombined   = false; if((*begin_new)->muonType() == xAOD::Muon_v1::Combined) isCombined=true;
	     bool isSegTag     = false; if((*begin_new)->muonType() == xAOD::Muon_v1::SegmentTagged) isSegTag=true;
	     bool isCaloTag    = false; if((*begin_new)->muonType() == xAOD::Muon_v1::CaloTagged) isCaloTag=true;
	     bool isStandAl    = false; if((*begin_new)->muonType() == xAOD::Muon_v1::MuonStandAlone) isStandAl=true;
	     M_isCombined              .push_back(isCombined);
	     M_isSegTag                .push_back(isSegTag);
	     M_isCaloTag               .push_back(isCaloTag);
	     M_isStandAl               .push_back(isStandAl);
         //Fill in the muon reco info
	     m_Muon_number++;
	     m_Muon_pt	                .push_back( (*begin_new)->pt()/1000.);
	     m_Muon_ptcorr              .push_back( (*begin_new)->auxdata< float >( "ptcorr" ));
	     m_Muon_eta	                .push_back( (*begin_new)->eta() );
	     m_Muon_phi	                .push_back( (*begin_new)->phi() );
	     m_Muon_charge	            .push_back( (*begin_new)->charge() );
	     //m_Muon_quality	            .push_back( m_muonSelection2->getQuality(*Muon) );
	     m_Muon_type	            .push_back( (*begin_new)->muonType() ); //Combined of sigment tag and so on.
	     m_Muon_ELoss               .push_back( (*begin_new)->floatParameter(xAOD::Muon::EnergyLoss) );
         m_Muon_isTight             .push_back( passTightID );
	     m_Muon_isMedium            .push_back( passMediumID );
	     m_Muon_MatchdR_mu4         .push_back(m_trigMuonMatching->minDelR((*begin_new),"HLT_mu4",1));
         m_Muon_MatchdR_mu6         .push_back(m_trigMuonMatching->minDelR((*begin_new),"HLT_mu6",1));
         m_Muon_MatchdR_mu8         .push_back(m_trigMuonMatching->minDelR((*begin_new),"HLT_mu8",1));
         m_Muon_MatchdR_mu15        .push_back(m_trigMuonMatching->minDelR((*begin_new),"HLT_mu15",1));
         m_Muon_MatchdR_mu15_MU10   .push_back(m_trigMuonMatching->minDelR((*begin_new),"HLT_mu15_L1MU10",1));
         m_Muon_MatchdR_mu10_MU6    .push_back(m_trigMuonMatching->minDelR((*begin_new),"HLT_mu10_L1MU6",1));
         m_Muon_MatchdR_mu15_MU6    .push_back(m_trigMuonMatching->minDelR((*begin_new),"HLT_mu15_L1MU6",1));



         is_a_good_muon = true;

/*
         //ID track information
         const xAOD::TrackParticle* Track = (*begin_new)->trackParticle( xAOD::Muon::InnerDetectorTrackParticle );
         if (Track)
            {
               m_Muon_passedIDCutsMed        .push_back( m_muonSelection2->passedIDCuts(*Track) );
               m_Muon_passedIDCutsTight      .push_back( m_muonSelection1->passedIDCuts(*Track) );
               m_Muon_Track_d0               .push_back( Track->d0()    );
               m_Muon_Track_z0               .push_back( Track->z0()    );
               m_Muon_id_theta               .push_back( Track->theta() );
               m_Muon_id_phi                 .push_back( Track->phi()   );
               m_Muon_id_pt                  .push_back( Track->pt()    );
               m_Muon_id_charge              .push_back( Track->charge());
               m_Muon_id_eta                 .push_back( Track->eta()   );


               int IdTrkValue=0;

               //Lets find the number of holes in the detector and so on left from muons
               IdTrkValue=Track->auxdata< unsigned char >("numberOfPixelHits");
               m_Muon_Track_nPixelHits            .push_back(IdTrkValue);
               IdTrkValue=Track->auxdata< unsigned char >("numberOfPixelDeadSensors");
               m_Muon_Track_nPixelDeadSensors     .push_back(IdTrkValue);
               IdTrkValue=Track->auxdata< unsigned char >("numberOfSCTHits");
               m_Muon_Track_nSCTHits              .push_back(IdTrkValue);
               IdTrkValue=Track->auxdata< unsigned char >("numberOfSCTDeadSensors");
               m_Muon_Track_nSCTDeadSensors       .push_back(IdTrkValue);
               IdTrkValue=Track->auxdata< unsigned char >("numberOfTRTHits");
               m_Muon_Track_nTRTHits              .push_back(IdTrkValue);
               IdTrkValue=Track->auxdata< unsigned char >("numberOfTRTOutliers");
               m_Muon_Track_nTRTOutliers          .push_back(IdTrkValue);
               IdTrkValue=Track->auxdata< unsigned char >("numberOfPixelHoles");
               m_Muon_Track_nPixelHoles           .push_back(IdTrkValue);
               IdTrkValue=Track->auxdata< unsigned char >("numberOfSCTHoles");
               m_Muon_Track_nSCTHoles             .push_back(IdTrkValue);
            } //End of track if
         else
            {
               m_Muon_passedIDCutsMed        .push_back( -2222 );
               m_Muon_passedIDCutsTight      .push_back( -2222 );
               m_Muon_Track_d0               .push_back( -2222 );
               m_Muon_Track_z0               .push_back( -2222 );
               m_Muon_Track_nPixelHits       .push_back( -2222 );
               m_Muon_Track_nPixelDeadSensors.push_back( -2222 );
               m_Muon_Track_nSCTHits         .push_back( -2222 );
               m_Muon_Track_nSCTDeadSensors  .push_back( -2222 );
               m_Muon_Track_nTRTHits         .push_back( -2222 );
               m_Muon_Track_nTRTOutliers     .push_back( -2222 );
               m_Muon_Track_nPixelHoles      .push_back( -2222 );
               m_Muon_Track_nSCTHoles        .push_back( -2222 );
               m_Muon_id_theta               .push_back( -2222 );
               m_Muon_id_phi                 .push_back( -2222 );
               m_Muon_id_pt                  .push_back( -2222 );
               m_Muon_id_charge              .push_back( -2222 );
               m_Muon_id_eta                 .push_back( -2222 );
            }//End of track else

          ////////////////////////////////////////////////////////////////////////////////////////
         //Lets do the ID efficiency///////////////////////////////////////////////////////////////
         /////////////////////////////////////////////////////////////////////////////////////////
         const xAOD::TrackParticle* MeTrack = (*begin_new)->trackParticle( xAOD::Muon::ExtrapolatedMuonSpectrometerTrackParticle );
         if(MeTrack)
            {

              m_me_pt                       .push_back(MeTrack->pt()*0.001);
              m_me_eta                      .push_back(MeTrack->eta());
              m_me_phi                      .push_back(MeTrack->phi());
              m_me_charge                   .push_back(MeTrack->charge());
              m_me_d0                       .push_back(MeTrack->d0());
              m_me_z0                       .push_back(MeTrack->z0());
              m_me_theta                    .push_back(MeTrack->theta());
              m_me_ed0                      .push_back(sqrt(MeTrack->definingParametersCovMatrix()(0,0)));

              int m_isTrack=0;
              ElementLink< xAOD::TrackParticleContainer > trkmulink = (*begin_new)->auxdata<ElementLink< xAOD::TrackParticleContainer > >("inDetTrackParticleLink");
              if(trkmulink.isValid())
               {
                 int keep = 1;
                 //Lets do some track cuts
                 if ((*trkmulink)->pt()*0.001<0.4) keep=0; //Pt cut
				 //eta cut
                 float eta=(*trkmulink)->eta();
	             if(fabs(eta)>2.5) keep=0; //eta cut
				 //theta cut
                 float theta=(*trkmulink)->theta();
	             if(theta==0.0) keep=0;
                 //hits
	             int nPixhits   = (*trkmulink)->auxdata< unsigned char >("numberOfPixelHits") +   (*trkmulink)->auxdata< unsigned char >("numberOfPixelDeadSensors") ;
                 int nSCThits   = (*trkmulink)->auxdata< unsigned char >("numberOfSCTHits")   +   (*trkmulink)->auxdata< unsigned char >("numberOfSCTDeadSensors");
	             int nPixholes  = (*trkmulink)->auxdata< unsigned char >("numberOfPixelHoles");
	             //int nSCTholes  = (*trkmulink)->auxdata< unsigned char >("numberOfSCTHoles");
                 int nIBLHits   = (*trkmulink)->auxdata< unsigned char >("numberOfInnermostPixelLayerHits");
                 int nBLHits    = (*trkmulink)->auxdata< unsigned char >("numberOfNextToInnermostPixelLayerHits");
                 int expIBLHits = (*trkmulink)->auxdata< unsigned char >("expectInnermostPixelLayerHit");
                 int expBLHits  = (*trkmulink)->auxdata< unsigned char >("expectNextToInnermostPixelLayerHit");
				 //cuts
	             if((nPixhits+nSCThits)<8)keep=0;
	             if((nPixhits+nSCThits)<10 && fabs(eta)>1.65)keep=0;
	             //if((nPixholes+nSCTholes)>=2)continue;
	             if( nPixholes >=1)keep=0;
	             if( (nIBLHits + nBLHits) <1 && expIBLHits && expBLHits)keep=0;
                 if(keep) m_isTrack=1;
                }//end if statement
               m_Muon_isTrack. push_back(m_isTrack);
            }//end of if me statement
*/
        }// End of the muon loop

       // m_Muon_n = m_Muon_number;
        //fill the tree
        //if(m_HLT_mu15 || m_HLT_mu15_L1MU10 || m_HLT_mu15_L1MU6)tree->Fill();



        begin_el = (electrons_corr.first)->begin();
        bool is_a_good_electron = false;
        for(;begin_el!=end_el;++begin_el)
       {

         xAOD::Electron* el1 = (*begin_el);

         //Quality
         bool passTight  = false;
	     bool passMedium = false;
	     bool passLoose  = false;

	     if (m_electronLHTightSelector ->accept(el1) )   passTight  = true; //based on eta quality and ID cuts - tight muons
	     if (m_electronLHMediumSelector->accept(el1) )   passMedium = true; //medium muons
	     if (m_electronLHLooseSelector ->accept(el1) )   passLoose  = true; //medium muons

	     if(passMedium  = false) continue;
	     if(el1->pt()/1000.<20)continue;
	     if(abs(el1->eta())>2.47)continue;
	     if(m_match_tool->matchHLT(el1,"HLT_e15_lhloose")==0)continue;


         m_El_pt	              .push_back(el1->pt()/1000.);
	     m_El_eta	              .push_back(el1->eta() );
	     m_El_phi	              .push_back(el1->phi() );
	     m_El_charge	          .push_back(el1->charge() );
         m_El_isTight             .push_back( passTight );
	     m_El_isMedium            .push_back( passMedium );
	     m_El_isLoose             .push_back( passLoose );
	     //trigger match
	     m_El_MatchdR_mu15         .push_back(m_match_tool->matchHLT(el1,"HLT_e15_lhloose"));





          is_a_good_electron = true;
        }

        if(!is_a_good_electron && !is_a_good_muon) return EL::StatusCode::SUCCESS;
        //if (m_HLT_e15_lhloose)tree2->Fill();
        tree->Fill();


delete electrons_corr.first;
delete electrons_corr.second;





    return EL::StatusCode::SUCCESS;
}



//Clear vector before each events
void Zreco :: clearVector ()
{


    m_Muon_MatchdR_mu4        	.clear();
    m_Muon_MatchdR_mu6        	.clear();
    m_Muon_MatchdR_mu8        	.clear();
    m_Muon_MatchdR_mu10       	.clear();
    m_Muon_MatchdR_mu15       	.clear();
    m_Muon_MatchdR_mu15_MU10   	.clear();
    m_Muon_MatchdR_mu10_MU6   	  .clear();
	m_Muon_MatchdR_mu15_MU6   	  .clear();


    m_Muon_MatchdR_L1MU4          .clear();
	m_Muon_MatchdR_L1MU6          .clear();
	m_Muon_MatchdR_L1MU10         .clear();
	m_Muon_MatchdR_L1MU15         .clear();

    isgap.clear();

    m_me_pt                       .clear();
    m_me_eta                      .clear();
    m_me_phi                      .clear();
    m_me_charge                   .clear();
    m_me_d0                       .clear();
    m_me_z0                       .clear();
    m_me_theta                    .clear();
    m_me_ed0                      .clear();


    m_Muon_pt 			.clear();
    m_Muon_ptcorr       .clear();
    m_Muon_eta			.clear();
    m_Muon_phi			.clear();
    m_Muon_charge		.clear();
    m_Muon_quality		.clear();
    m_Muon_type			.clear();


    m_Muon_passedIDCutsMed 	.clear();
    m_Muon_passedIDCutsTight 	.clear();
    m_Muon_ELoss 		.clear();
    m_Muon_reducedChi2    	.clear();
    m_Muon_rho    		.clear();
    m_Muon_me_pt    		.clear();
    m_Muon_me_phi    		.clear();
    m_Muon_me_theta    		.clear();
    m_Muon_me_charge        .clear();
    m_Muon_me_eta           .clear();
    m_Muon_id_pt    		.clear();
    m_Muon_id_phi    		.clear();
    m_Muon_id_theta    		.clear();
    m_Muon_id_charge        .clear();
    m_Muon_id_eta           .clear();
    m_Muon_qOverPsigma  	.clear();
    m_Muon_qOverPsignif  	.clear();
   // m_Muon_is_truth             .clear();


    m_Muon_isTight  		.clear();
    m_Muon_isMedium  		.clear();
    m_Muon_isTrack          .clear();



    m_Muon_Track_d0               .clear();
    m_Muon_Track_z0               .clear();
    m_Muon_Track_nPixelHits       .clear();
    m_Muon_Track_nPixelDeadSensors.clear();
    m_Muon_Track_nSCTHits         .clear();
    m_Muon_Track_nSCTDeadSensors  .clear();
    m_Muon_Track_nTRTHits         .clear();
    m_Muon_Track_nTRTOutliers     .clear();
    m_Muon_Track_nPixelHoles      .clear();
    m_Muon_Track_nSCTHoles        .clear();


    m_vertex_PV_x                    .clear();
    m_vertex_PV_y                    .clear();
    m_vertex_PV_z                    .clear();


    m_vertex_PU_x                    .clear();
    m_vertex_PU_y                    .clear();
    m_vertex_PU_z                    .clear();


    M_isCombined              .clear();
	M_isSegTag                .clear();
	M_isCaloTag               .clear();
	M_isStandAl               .clear();

	m_El_MatchdR_mu15       	.clear();
	m_El_isTight  		.clear();
    m_El_isMedium  		.clear();
    m_El_pt    		.clear();
    m_El_phi    		.clear();
    m_El_charge        .clear();
    m_El_eta           .clear();

    m_El_isLoose. clear();






}

bool Zreco :: selectCluster( float eta, float cl_cell_sig, int cl_cell_sig_samp)
{
  float cut[] = {0      ,4.7426, 5.11018,5.07498,5.0969, 5.10695,5.04098,5.07106,4.98087,5.11647,
                 5.08988,5.16267,5.17202,5.23803,5.25314,5.29551,5.35092,5.40863,5.44375,5.38075,
                 5.25022,5.37933,5.25459,5.37719,5.25169,5.73985,5.79174,5.79266,5.79588,5.7963,
                 5.81949,5.82273,5.85658,5.85442,5.84779,5.77679,5.83323,5.84524,5.84439,5.84488,
                 5.84744,5.84683,5.84524,5.84594,5.84656,5.84639,5.84461,5.84515,5.84206,5.8396,
                 5.84497,5.84801,5.84608,5.84608,5.84783,5.84726,5.84844,5.8477, 5.84796,5.84757,
                 5.84822,5.84814,5.84617,5.83451,5.77658,5.84309,5.85496,5.85761,5.82555,5.82206,
                 5.78982,5.78482,5.7778, 5.78327,5.74898,5.25459,5.37503,5.25459,5.37283,5.25169,
                 5.37862,5.44473,5.41041,5.34498,5.29551,5.25602,5.2283, 5.17428,5.14504,5.09342,
                 5.12256,4.98721,5.07106,5.02642,5.10031,5.11018,5.05447,5.10031,4.7426 ,0};

  //eta-dependent cut on the signficance of the cell in the cluster with the maximum significance
  for(int i=0; i<100; i++)
  {
    if( eta > (-5 + 0.1*i) && eta <= (-5 + 0.1*(i+1)) )
    {
      //std::cout << "selectCluster: cell_sig = " << cl_cell_sig << std::endl;
      //skip the cluster if its cell significance is not above the threshold
      if(cl_cell_sig < cut[i])
      {
        //std::cout << "cluster rejected: cell_sig = " << cl_cell_sig << std::endl;
        return false;
      }

      //if the most significant cell is in the range of layers corresponding to the tile, reject it
      else if(cl_cell_sig_samp>=CaloSampling::TileBar0&&cl_cell_sig_samp<=CaloSampling::TileExt2)
        return false;

      else
        return true;
    }//end of if
  }//end of for loop

  return false;
}


EL::StatusCode Zreco :: postExecute ()
{
  // Here you do everything that needs to be done after the main event
  // processing.  This is typically very rare, particularly in user
  // code.  It is mainly used in implementing the NTupleSvc.
  return EL::StatusCode::SUCCESS;
}



EL::StatusCode Zreco :: finalize ()
{

    if (m_grl)                  { delete m_grl;                 m_grl = 0;                       }
    if (configTool)             { delete configTool;            configTool = 0;                 }
    if (trigDecTool)            { delete trigDecTool;           trigDecTool = 0;                }
    if (m_muonSelection)        { delete m_muonSelection;       m_muonSelection = 0;            }
    if (m_muonSelection1)       { delete m_muonSelection1;      m_muonSelection1 = 0;           }
    if (m_muonSelection2)       { delete m_muonSelection2;      m_muonSelection2 = 0;           }
    if (m_muonSelection3)       { delete m_muonSelection3;      m_muonSelection3 = 0;           }
    if (m_muonCorr)             { delete m_muonCorr;            m_muonCorr = 0;                 }
    if (m_trigMuonMatching) 	{ delete m_trigMuonMatching; 	m_trigMuonMatching = 0; 	    }
    if (m_tmt) 			        { delete m_tmt; 		        m_tmt = 0; 			            }
    if (centTool)				{ delete centTool;				centTool =0;                    }

     if(m_electronLHTightSelector)            {delete m_electronLHTightSelector;     		m_electronLHTightSelector=0; }
    if(m_electronLHMediumSelector)           {delete m_electronLHMediumSelector; 			m_electronLHMediumSelector=0; }
    if(m_electronLHLooseSelector)            {delete m_electronLHLooseSelector;  			m_electronLHLooseSelector=0; }
    // if(m_egammaCalibrationAndSmearingTool)  {delete m_egammaCalibrationAndSmearingTool; 	m_egammaCalibrationAndSmearingTool=0; }
      if(m_match_tool)                         {delete m_match_tool; 							m_match_tool=0; }


  return EL::StatusCode::SUCCESS;
}



EL::StatusCode Zreco :: histFinalize ()
{
  // This method is the mirror image of histInitialize(), meaning it
  // gets called after the last event has been processed on the worker
  // node and allows you to finish up any objects you created in
  // histInitialize() before they are written to disk.  This is
  // actually fairly rare, since this happens separately for each
  // worker node.  Most of the time you want to do your
  // post-processing on the submission node after all your histogram
  // outputs have been merged.  This is different from finalize() in
  // that it gets called on all worker nodes regardless of whether
  // they processed input events.
  return EL::StatusCode::SUCCESS;
}



int Zreco :: is_this_good_bcid(int run, int bcid)
{
    switch(run)
          {
    		case 313063: return 0;
			case 313067: return(bcid==1 || bcid==177 || bcid==892 || bcid==1068 || bcid==1783 || bcid==1959);
			case 313100: return(bcid==1 || bcid==177 || bcid==353 || bcid==892 || bcid==1068 || bcid==1244 || bcid==1420 || bcid==1596 || bcid==1783 || bcid==1959 || bcid==2135 || bcid==2311 || bcid==2487 || bcid==2674 || bcid==2850 || bcid==3286);
			case 313107: return(bcid==1 || bcid==177 || bcid==353 || bcid==529 || bcid==705 || bcid==892 || bcid==1068 || bcid==1420 || bcid==1783 || bcid==1959 || bcid==2135 || bcid==2311 || bcid==2487 || bcid==2663 || bcid==2839 || bcid==3015 || bcid==3183 || bcid==3283);
			case 313136: return(bcid==1 || bcid==177 || bcid==353 || bcid==529 || bcid==705 || bcid==892 || bcid==1068 || bcid==1420 || bcid==1783 || bcid==1959 || bcid==2135 || bcid==2311 || bcid==2487 || bcid==2663 || bcid==2839 || bcid==3015 || bcid==3183 || bcid==3283);
			case 313187: return(bcid==1 || bcid==177 || bcid==353 || bcid==529 || bcid==705 || bcid==892 || bcid==1068 || bcid==1420 || bcid==1783 || bcid==1959 || bcid==2135 || bcid==2311 || bcid==2487 || bcid==2663 || bcid==2839 || bcid==3015 || bcid==3183 || bcid==3283);
			case 313259: return(bcid==1 || bcid==177 || bcid==353 || bcid==529 || bcid==705 || bcid==2490 || bcid==2677 || bcid==2821 || bcid==2853 || bcid==2985 || bcid==3029 || bcid==3109 || bcid==3209);
			case 313285: return(bcid==1 || bcid==177 || bcid==353 || bcid==529 || bcid==705 || bcid==892 || bcid==1068 || bcid==1420 || bcid==1783 || bcid==1959 || bcid==2135 || bcid==2311 || bcid==2487 || bcid==2663 || bcid==2839 || bcid==3015 || bcid==3183 || bcid==3283);
			case 313295: return(bcid==1 || bcid==177 || bcid==353 || bcid==529 || bcid==705 || bcid==892 || bcid==1068 || bcid==1420 || bcid==1783 || bcid==1959 || bcid==2135 || bcid==2311 || bcid==2487 || bcid==2663 || bcid==2839 || bcid==3015 || bcid==3183 || bcid==3283);
			case 313333: return(bcid==82 || bcid==258 || bcid==434 || bcid==610 || bcid==793 || bcid==1691 || bcid==1871 || bcid==2047 || bcid==2223 || bcid==2399 || bcid==2586 || bcid==2758 || bcid==2934 || bcid==3110 || bcid==3286);
			case 313435: return(bcid==1 || bcid==177 || bcid==353 || bcid==892 || bcid==1068 || bcid==1244 || bcid==1420 || bcid==1596 || bcid==1783 || bcid==1959 || bcid==2135 || bcid==2311 || bcid==2487 || bcid==2674 || bcid==2850 || bcid==3286);
			case 313572: return 0;
			case 313574: return(bcid==1 || bcid==177 || bcid==892 || bcid==1068 || bcid==1783 || bcid==1959);
			case 313575: return(bcid==73 || bcid==785 || bcid==959 || bcid==1139 || bcid==1315 || bcid==1491 || bcid==1680 || bcid==1852 || bcid==2031 || bcid==2207 || bcid==2387 || bcid==2571 || bcid==2922 || bcid==3098 || bcid==3274);
			case 313603: return(bcid==73 || bcid==785 || bcid==959 || bcid==1139 || bcid==1315 || bcid==1491 || bcid==1680 || bcid==1852 || bcid==2031 || bcid==2207 || bcid==2387 || bcid==2571 || bcid==2922 || bcid==3098 || bcid==3274);
			case 313629: return(bcid==73 || bcid==785 || bcid==959 || bcid==1139 || bcid==1315 || bcid==1491 || bcid==1680 || bcid==1852 || bcid==2031 || bcid==2207 || bcid==2387 || bcid==2571 || bcid==2922 || bcid==3098 || bcid==3274);
			case 313630: return(bcid==73 || bcid==785 || bcid==959 || bcid==1139 || bcid==1315 || bcid==1491 || bcid==1680 || bcid==1852 || bcid==2031 || bcid==2207 || bcid==2387 || bcid==2571 || bcid==2922 || bcid==3098 || bcid==3274);
			case 313688: return(bcid==73 || bcid==785 || bcid==959 || bcid==1139 || bcid==1315 || bcid==1491 || bcid==1680 || bcid==1852 || bcid==2031 || bcid==2207 || bcid==2387 || bcid==2571 || bcid==2922 || bcid==3098 || bcid==3274);
			case 313695: return(bcid==73 || bcid==785 || bcid==959 || bcid==1139 || bcid==1315 || bcid==1491 || bcid==1680 || bcid==1852 || bcid==2031 || bcid==2207 || bcid==2387 || bcid==2571 || bcid==2922 || bcid==3098 || bcid==3274);
			case 313833: return(bcid==73 || bcid==785 || bcid==959 || bcid==1139 || bcid==1315 || bcid==1491 || bcid==1680 || bcid==1852 || bcid==2031 || bcid==2207 || bcid==2387 || bcid==2571 || bcid==2922 || bcid==3098 || bcid==3274);
			case 313878: return(bcid==75 || bcid==251 || bcid==427 || bcid==603 || bcid==780 || bcid==958 || bcid==1134 || bcid==1310 || bcid==1486 || bcid==1678 || bcid==1856 || bcid==2032 || bcid==2208 || bcid==2384 || bcid==2576 || bcid==2756 || bcid==2932 || bcid==3108 || bcid==3284);
			case 313929: return(bcid==75 || bcid==251 || bcid==427 || bcid==603 || bcid==780 || bcid==958 || bcid==1134 || bcid==1310 || bcid==1486 || bcid==1678 || bcid==1856 || bcid==2032 || bcid==2208 || bcid==2384 || bcid==2576 || bcid==2756 || bcid==2932 || bcid==3108 || bcid==3284);
			case 313935: return(bcid==75 || bcid==251 || bcid==427 || bcid==603 || bcid==780 || bcid==958 || bcid==1134 || bcid==1310 || bcid==1486 || bcid==1678 || bcid==1856 || bcid==2032 || bcid==2208 || bcid==2384 || bcid==2576 || bcid==2756 || bcid==2932 || bcid==3108 || bcid==3284);
			case 313984: return(bcid==75 || bcid==251 || bcid==427 || bcid==603 || bcid==780 || bcid==958 || bcid==1134 || bcid==1310 || bcid==1486 || bcid==1678 || bcid==1856 || bcid==2032 || bcid==2208 || bcid==2384 || bcid==2576 || bcid==2756 || bcid==2932 || bcid==3108 || bcid==3284);
			case 314014: return(bcid==75 || bcid==251 || bcid==427 || bcid==603 || bcid==780 || bcid==958 || bcid==1134 || bcid==1310 || bcid==1486 || bcid==1678 || bcid==1856 || bcid==2032 || bcid==2208 || bcid==2384 || bcid==2576 || bcid==2756 || bcid==2932 || bcid==3108 || bcid==3284);
			case 314077: return(bcid==75 || bcid==251 || bcid==427 || bcid==603 || bcid==780 || bcid==958 || bcid==1134 || bcid==1310 || bcid==1486 || bcid==1678 || bcid==1856 || bcid==2032 || bcid==2208 || bcid==2384 || bcid==2576 || bcid==2756 || bcid==2932 || bcid==3108 || bcid==3284);
			case 314105: return(bcid==75 || bcid==251 || bcid==427 || bcid==603 || bcid==780 || bcid==958 || bcid==1134 || bcid==1310 || bcid==1486 || bcid==1678 || bcid==1856 || bcid==2032 || bcid==2208 || bcid==2384 || bcid==2576 || bcid==2756 || bcid==2932 || bcid==3108 || bcid==3284);
			case 314112: return(bcid==75 || bcid==251 || bcid==427 || bcid==603 || bcid==780 || bcid==958 || bcid==1134 || bcid==1310 || bcid==1486 || bcid==1678 || bcid==1856 || bcid==2032 || bcid==2208 || bcid==2384 || bcid==2576 || bcid==2756 || bcid==2932 || bcid==3108 || bcid==3284);
			case 314157: return(bcid==75 || bcid==251 || bcid==427 || bcid==603 || bcid==780 || bcid==958 || bcid==1134 || bcid==1310 || bcid==1486 || bcid==1678 || bcid==1856 || bcid==2032 || bcid==2208 || bcid==2384 || bcid==2576 || bcid==2756 || bcid==2932 || bcid==3108 || bcid==3284);
	        case 314170: return(bcid==75 || bcid==251 || bcid==427 || bcid==603 || bcid==780 || bcid==958 || bcid==1134 || bcid==1310 || bcid==1486 || bcid==1678 || bcid==1856 || bcid==2032 || bcid==2208 || bcid==2384 || bcid==2576 || bcid==2756 || bcid==2932 || bcid==3108 || bcid==3284);

            default: return -1;

           }







}



