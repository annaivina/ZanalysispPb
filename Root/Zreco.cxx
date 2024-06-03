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

    ////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////
    //TREE//////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////
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
    tree2->Branch("El_MatchdR_mu15",  		    &m_El_MatchdR_mu15);//5


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

    ////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////
    //TREE2/////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////

    //Second tree
    //Track containers
    /*
    tree2->Branch("Track_pt",		             &m_Track_pt);
    tree2->Branch("Track_eta",		             &m_Track_eta);
    tree2->Branch("Track_phi",		             &m_Track_phi);
    tree2->Branch("Track_charge", 	             &m_Track_charge);
    tree2->Branch("Track_d0",              	     &m_Track_d0);
    tree2->Branch("Track_ed0",              	 &m_Track_ed0);
    tree2->Branch("Track_ez0",              	 &m_Track_ez0);
    tree2->Branch("Track_z0",         	         &m_Track_z0);
    tree2->Branch("Track_ez0sin",				 &m_Track_ez0sin);
    tree2->Branch("Track_z0_calc",         	     &m_Track_z0_calc);
    tree2->Branch("Track_theta",                 &m_Track_theta);
    tree2->Branch("Track_nPixelHits",            &m_Track_nPixelHits);
    tree2->Branch("Track_nPixelDeadSensors",     &m_Track_nPixelDeadSensors);
    tree2->Branch("Track_nSCTHits",              &m_Track_nSCTHits);
    tree2->Branch("Track_nSCTDeadSensors",       &m_Track_nSCTDeadSensors);
    tree2->Branch("Track_nIBLHits",       		 &m_Track_nIBLHits);
    tree2->Branch("Track_nBLHits",       		 &m_Track_nBLHits);
    tree2->Branch("Track_expIBLHits",            &m_Track_expIBLHits);
    tree2->Branch("Track_expBLHits",             &m_Track_expBLHits);
    tree2->Branch("Track_nTRTHits",              &m_Track_nTRTHits);
    tree2->Branch("Track_nTRTOutliers",          &m_Track_nTRTOutliers);
    tree2->Branch("Track_nPixelHoles",           &m_Track_nPixelHoles);
    tree2->Branch("Track_nSCTHoles",             &m_Track_nSCTHoles);
    */

/*
    //For the calo tagged muons - probe
    tree2->Branch("Muon_isCalotagged",     	     &m_isCaloTag_c);
    tree2->Branch("Muon_ispassedCalo",     		 &m_passed_calo_c);
    tree2->Branch("Muon_ispassedIDcuts",     	 &m_passed_IDcuts_c);
    tree2->Branch("Muon_calo_pt",     	         &m_caloPt);
    tree2->Branch("Muon_calo_ptcorr",     	     &m_caloPtcorr);
    tree2->Branch("Muon_calo_eta",     	         &m_caloEta);
    tree2->Branch("Muon_calo_phi",     	         &m_caloPhi);
    tree2->Branch("Muon_calo_charge",     	     &m_caloCharge);
    tree2->Branch("Muon_calo_id_d0",			 &m_Muonidcalo_d0);
    tree2->Branch("Muon_calo_id_d0sin",			 &m_Muonidcalo_d0sig);
    tree2->Branch("Muon_calo_id_theta",			 &m_Muonidcalo_theta);
    tree2->Branch("Muon_calo_id_z0",			 &m_Muonidcalo_z0);
    tree2->Branch("Muon_calo_id_z0sig",			 &m_Muonidcalo_z0sig);
    tree2->Branch("Muon_calo_id_ez0sin",		 &m_Muonidcalo_ez0sin);
    tree2->Branch("Muon_calo_id_phi",			 &m_Muonidcalo_phi);
    tree2->Branch("Muon_calo_id_pt",			 &m_Muonidcalo_pt);
    tree2->Branch("Muon_calo_id_charge",		 &m_Muonidcalo_charge);
    tree2->Branch("Muon_calo_id_eta",		     &m_Muonidcalo_eta);
    tree2->Branch("Muon_calo_nPixelHit",         &m_Muonidcalo_nPixelHits);
    tree2->Branch("Muon_calo_nPixelDeadSensors", &m_Muonidcalo_nPixelDeadSensors);
    tree2->Branch("Muon_calo_nSCTHits",          &m_Muonidcalo_nSCTHits);
    tree2->Branch("Muon_calo_nSCTDeadSensors",   &m_Muonidcalo_nSCTDeadSensors);
    tree2->Branch("Muon_calo_nPixelHoles",       &m_Muonidcalo_nPixelHoles);
    tree2->Branch("Muon_calo_nSCTHoles",         &m_Muonidcalo_nSCTHoles);
    tree2->Branch("Muon_calo_nTRTHits",          &m_Muonidcalo_nTRTHits);
    tree2->Branch("Muon_calo_nTRTOutliers",      &m_Muonidcalo_nTRTOutliers);
    tree2->Branch("Muon_calo_nIBLHits",          &m_Muonidcalo_nIBLHits);
    tree2->Branch("Muon_calo_nBLHits",           &m_Muonidcalo_nBLHits);
    tree2->Branch("Muon_calo_expIBLHits",        &m_Muonidcalo_expIBLHits);
    tree2->Branch("Muon_calo_expBLHits",         &m_Muonidcalo_expBLHits);
    tree2->Branch("Muon_isMedium_macthed",		 &m_passMediumID_c);
    tree2->Branch("Muon_isTight_matched",		 &m_passTightID_c);
    tree2->Branch("dR_medium",					 &deltaR_med);
    tree2->Branch("dR_tight",					 &deltaR_tight);
    tree2->Branch("Track_nPixelHits_calo",       &nPixelHits_calo);
    tree2->Branch("Track_nPixelDeadSensors_calo",&nPixelDeadSensors_calo);
    tree2->Branch("Track_nSCTHits_calo",         &nSCTHits_calo);
    tree2->Branch("Track_nSCTDeadSensors_calo",  &nSCTDeadSensors_calo);
    tree2->Branch("Track_nPixelHoles_calo",      &nPixelHoles_calo);
    tree2->Branch("Track_nSCTHoles_calo",        &nSCTHoles_calo);
    tree2->Branch("Track_nTRTHits_calo",         &nTRTHits_calo);
    tree2->Branch("Track_nTRTOutliers_calo",     &nTRTOutliers_calo);
    tree2->Branch("Track_nIBLHits_calo",         &nIBLHits_calo);
    tree2->Branch("Track_nBLHits_calo",          &nBLHits_calo);
    tree2->Branch("Track_expIBLHits_calo",       &expIBLHits_calo);
    tree2->Branch("Track_expBLHits_calo",        &expBLHits_calo);
    tree2->Branch("Track_calo_passIDcuts",		 &m_calotrack_passIDcuts);


    //muon container - tag
    tree2->Branch("Muon_Pt",   	     		     &m_muonPt);
    tree2->Branch("Muon_Ptcorr",   	     		 &m_muonPtcorr);
    tree2->Branch("Muon_Eta",   	     		 &m_muonEta);
    tree2->Branch("Muon_Phi",   	     		 &m_muonPhi);
    tree2->Branch("Muon_Charge",      		     &m_muonCharge);
    tree2->Branch("Muon_Quality",     		     &m_muonQuality);
    tree2->Branch("Muon_Type",        		     &m_muonType);
    tree2->Branch("Muon_ELoss",    		         &m_muonELoss);
    tree2->Branch("Muon_isCombined_2",     		 &m_isCombined_2);
    tree2->Branch("Muon_isSegTag_2",     		 &m_isSegTag_2);
    tree2->Branch("Muon_isCaloTag_2",     		 &m_isCaloTag_2);
    tree2->Branch("Muon_isStandAl_2",            &m_isStandAl_2);
    tree2->Branch("Muon_MatchdR_mu15",           &m_MatchdR_mu15);
    tree2->Branch("Muon_isMedium",				 &m_passMediumID);
    tree2->Branch("Muon_isTight",				 &m_passTightID);
    tree2->Branch("Muon_passIDcuts",			 &m_passIDcuts);
    tree2->Branch("Muon_id_d0",				     &m_Muonid_d0);
    tree2->Branch("Muon_id_d0sin",				 &m_Muonid_d0sig);
    tree2->Branch("Muon_id_theta",				 &m_Muonid_theta);
    tree2->Branch("Muon_id_z0",				 	 &m_Muonid_z0);
    tree2->Branch("Muon_id_z0sig",				 &m_Muonid_z0sig);
    tree2->Branch("Muon_id_ez0sin",				 &m_Muonid_ez0sin);
    tree2->Branch("Muon_id_phi",				 &m_Muonid_phi);
    tree2->Branch("Muon_id_pt",				     &m_Muonid_pt);
    tree2->Branch("Muon_id_charge",				 &m_Muonid_charge);
    tree2->Branch("Muon_id_eta",				 &m_Muonid_eta);
    tree2->Branch("Track_passIDcuts",			 &m_track_passIDcuts);
    //the di-muon kinematics
    tree2->Branch("Pair_pt",                      &pairpt);
    tree2->Branch("Pair_m",                       &pairm);
    tree2->Branch("Pair_y",                       &pairy);
    tree2->Branch("Pair_eta",                     &paireta);
    tree2->Branch("Pair_phi",                     &pairphi);
    tree2->Branch("PrimaryVertex",				  &m_PV);
*/
    /*
    //matching
    tree2->Branch("Muon_isMuon_loose",            &m_isMuon_loose);
    tree2->Branch("Muon_isMuon_medium",           &m_isMuon_medium);
    tree2->Branch("Muon_isMuon_tight",            &m_isMuon_tight);
    */
/*
    //fcal and pileup
    tree2->Branch("sumEt",                        &m_sumEt);
    tree2->Branch("percentile",                   &percentile);
    tree2->Branch("Track_number_pileup",		  &trackn);
    tree2->Branch("PileUp_vertices",			  &pile_up_vertices);
    ////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////
    //TREE3/////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////
    tree3->Branch("Muon_Pt_3",   	     		 &m_muonPt_3);
    tree3->Branch("Muon_Eta_3",   	     		 &m_muonEta_3);
    tree3->Branch("Muon_Phi_3",   	     		 &m_muonPhi_3);
    tree3->Branch("Muon_Charge_3",      		 &m_muonCharge_3);
    tree3->Branch("Muon_Quality_3",     		 &m_muonQuality_3);
    tree3->Branch("Muon_Type_3",        		 &m_muonType_3);
    //tree3->Branch("MuonELoss_3",    		     &m_muonELoss_3);
    tree3->Branch("Muon_isCombined_3",     		 &m_isCombined_3);
    tree3->Branch("Muon_isSegTag_3",     		 &m_isSegTag_3);
    tree3->Branch("Muon_isCaloTag_3",     		 &m_isCaloTag_3);
    tree3->Branch("Muon_isStandAl_3",            &m_isStandAl_3);
    tree3->Branch("Muon_passTightID_3",          &m_passTightID_3);
	tree3->Branch("Muon_passMediumID_3",         &m_passMediumID_3);
	tree3->Branch("Muon_passLooseID_3",          &m_passLooseID_3);
	tree3->Branch("Muon_passIDcuts_3",			 &m_passIDcuts_3);
	tree3->Branch("Muon_id_d0_3",				 &m_Muonid_d0_3);
    tree3->Branch("Muon_id_d0sin_3",		     &m_Muonid_d0sig_3);
    tree3->Branch("Muon_id_theta_3",			 &m_Muonid_theta_3);
    tree3->Branch("Muon_id_z0_3",				 &m_Muonid_z0_3);
    tree3->Branch("Muon_id_z0sig_3",		     &m_Muonid_z0sig_3);
    tree3->Branch("Muon_id_ez0sin_3",			 &m_Muonid_ez0sin_3);
    tree3->Branch("Muon_id_phi_3",				 &m_Muonid_phi_3);
    tree3->Branch("Muon_id_pt_3",				 &m_Muonid_pt_3);
    tree3->Branch("Muon_id_charge_3",		     &m_Muonid_charge_3);
    tree3->Branch("Muon_id_eta_3",				 &m_Muonid_eta_3);
    tree3->Branch("Muon_id_track_passIDcuts",    &m_Muonidtrack_passIDcuts);

    tree3->Branch("Muon_isCaloTag_4",     		 &m_isCaloTag_4);
    tree3->Branch("Muon_isCombined_4",     		 &m_isCombined_4);
    tree3->Branch("Muon_passTightID_4",          &m_passTightID_4);
	tree3->Branch("Muon_passMediumID_4",         &m_passMediumID_4);
	tree3->Branch("Muon_passLooseID_4",          &m_passLooseID_4);
    tree3->Branch("passed",                      &passed);
    tree3->Branch("ME_Track_pt",                 &m_ME_Track_pt);
    tree3->Branch("ME_Track_eta",                &m_ME_Track_eta);
    tree3->Branch("ME_Track_phi",                &m_ME_Track_phi);
    tree3->Branch("ME_Track_charge",             &m_ME_Track_charge);
    tree3->Branch("ME_Track_d0",                 &m_ME_Track_d0);
    tree3->Branch("ME_Track_z0",                 &m_ME_Track_z0);
    tree3->Branch("ME_Track_theta",              &m_ME_Track_theta);
    tree3->Branch("ME_Track_ed0",                &m_ME_Track_ed0);

    tree3->Branch("Pair_pt_3",                    &pairpt_3);
    tree3->Branch("Pair_m_3",                     &pairm_3);
    tree3->Branch("Pair_y_3",                     &pairy_3);
    tree3->Branch("Pair_phi_3",                   &pairphi_3);
    tree3->Branch("Pair_eta_3",                   &paireta_3);
    tree3->Branch("PrimaryVertex",				  &m_PV);
	tree3->Branch("istrack",                      &m_istrack);
*/
	/*
	tree3->Branch("dR",                           &m_dR);
	tree3->Branch("d0",                           &d0_trk);
	tree3->Branch("ed0",                          &ed0_trk);
	tree3->Branch("ed0",                          &ez0_trk);
	tree3->Branch("z0",                           &z0_trk);
	tree3->Branch("theta",                        &theta_trk);
	tree3->Branch("pt",                           &pt_trk);
	tree3->Branch("eta",                          &eta_trk);
	tree3->Branch("phi",                          &phi_trk);
	*/
/*
    tree3->Branch("sumEt",                        &m_sumEt);
    tree3->Branch("percentile",                   &percentile);
    tree3->Branch("Track_number_pileup",		  &trackn);
    tree3->Branch("PileUp_vertices",		      &pile_up_vertices);

    tree3->Branch("Muon_iscalo",				  &is_calo);
    tree3->Branch("Muon_isID",					  &is_ID);
    tree3->Branch("Muon_passcalo",				  &pass_calo);
    tree3->Branch("dR_match",					  &dR);


    //Using the track link
	//tree3->Branch("dR_old",                      &m_dR_old);
	//ree3->Branch("istrack_old",                 &m_istrack_old);

	//tree3->Branch("d0_old",                      &d0_old);
	//tree3->Branch("ed0_old",                     &ed0_old);
	//tree3->Branch("z0_old",                      &z0_old);
	//tree3->Branch("theta_old",                   &theta_old);


	tree3->Branch("Muon2_pt",   	     		 &m_muon2_pt);
    tree3->Branch("Muon2_eta",   	     		 &m_muon2_eta);
    tree3->Branch("Muon2_phi",   	     		 &m_muon2_phi);
    tree3->Branch("Muon2_charge",      		     &m_muon2_charge);
    tree3->Branch("MatchdR_mu15_3",              &MatchdR_mu15_3);
*/
    /*
    tree3->Branch("nPixelHits_3",                &nPixelHits_3);
    tree3->Branch("nPixelDeadSensors_3",         &nPixelDeadSensors_3);
    tree3->Branch("nSCTHits_3",                  &nSCTHits_3);
    tree3->Branch("nSCTDeadSensors_3",           &nSCTDeadSensors_3);
    tree3->Branch("nPixelHoles_3",               &nPixelHoles_3);
    tree3->Branch("nSCTHoles_3",                 &nSCTHoles_3);
    tree3->Branch("nTRTHits_3",                  &nTRTHits_3);
    tree3->Branch("nTRTOutliers_3",              &nTRTOutliers_3);
    tree3->Branch("nIBLHits_3",                  &nIBLHits_3);
    tree3->Branch("nBLHits_3",                   &nBLHits_3);
    tree3->Branch("expIBLHits_3",                &expIBLHits_3);
    tree3->Branch("expBLHits_3",                 &expBLHits_3);
    */


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


/*

    myPhTight = new AsgPhotonIsEMSelector("myPhTight");
    myPhTight->setProperty("ConfigFile", Settings_PhotonTightConfig);
    myPhTight->setProperty("isEMMask",static_cast<unsigned int> (egammaPID::PhotonTight) );
    myPhTight->setProperty("PIDName",static_cast<int> (egammaPID::IsEMTight) );
    myPhTight->initialize();


    myTrigMedium = new AsgPhotonIsEMSelector("myTrigMedium");
    myTrigMedium->setProperty("ConfigFile",Settings_HLTPhotonMediumConfig);
    myTrigMedium->setProperty("isEMMask",static_cast<unsigned int> (egammaPID::PhotonMedium) );
    myTrigMedium->setProperty("PIDName",static_cast<int> (egammaPID::IsEMMedium) );
    myTrigMedium->setProperty("ForceConvertedPhotonPID", true);
    myTrigMedium->initialize();
*/

/*
    //photon selector tool
    m_photon_tight = new AsgPhotonIsEMSelector ( "PhotonTightIsEMSelector" );
    m_photon_tight->setProperty("isEMMask",egammaPID::PhotonTight);
    m_photon_tight->setProperty("ConfigFile","ElectronPhotonSelectorTools/offline/mc15_20150712/PhotonIsEMTightSelectorCutDefs.conf");
    if(!m_photon_tight->initialize().isSuccess()) {
       Fatal("MyFunction", "Failed to initialize PhotonTightIsEMSelector");}
*/
/*
    //photon selector tool
    m_photon_loose = new AsgPhotonIsEMSelector ( "PhotonTightIsEMSelector" );
    m_photon_loose->setProperty("isEMMask",egammaPID::PhotonLoose);
    m_photon_loose->setProperty("ConfigFile","ElectronPhotonSelectorTools/offline/mc15_20150408/PhotonIsEMLooseSelectorCutDefs.conf");
    if(!m_photon_loose->initialize().isSuccess()) {
       Fatal("MyFunction", "Failed to initialize PhotonTightIsEMSelector");}
*/
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



/*
    cout<<"=============================================================="<<endl;
    cout<<"----------------CLooking at the eta gap ...-----------------"<<endl;
    //calo cluster container
    const xAOD::CaloClusterContainer *Calocluster = 0;
    ANA_CHECK(event->retrieve( Calocluster, "CaloCalTopoClusters" ));
    //xAOD::CaloClusterContainer::const_iterator cal_beg = cluster->begin();
    //xAOD::CaloClusterContainer::const_iterator cal_end = cluster->end();
    for (const auto cluster : *Calocluster)
        {

            gap_decision=false;
			float clustEta = cluster->eta();
			float clustCellSig = cluster->getMomentValue( xAOD::CaloCluster::MomentType::CELL_SIGNIFICANCE );
			int clustCellSigSamp = cluster->getMomentValue( xAOD::CaloCluster::MomentType::CELL_SIG_SAMPLING );

            gap_decision = selectCluster(clustEta,clustCellSig,clustCellSigSamp); isgap.push_back(gap_decision);
		}
*/

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


   //if (!m_record) return EL::StatusCode::SUCCESS;
    //cleaning the vectors for new event



    //const xAOD::ZdcModuleContainer* zdcMod = 0;
    //ANA_CHECK(event->retrieve( zdcMod, "ZdcModules"));
    //const xAOD::HIEventShapeContainer* hiev = 0;
    //ANA_CHECK(event->retrieve( hiev, "HIEventShape"));
    // ZDC
    //m_zdcTools->reprocessZdc();
    // is Pileup
    //m_is_pileup = m_hiPileup->is_pileup( *hiev, *zdcMod); // SAVE pileup Decision HERE 0 = NO pileup, 1 = pileup
    //cout<<"Im here"<<endl;
    ///////////////////////////////
    //SUM ET AND CENTRALITY TOOL
    ////////////////////////////////

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
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //1.  Lets make the muon loop for the MUON VECTORS in order to det the info about the muons and their tracks in general
    //also I have added here the loop corresponding to the ME tracks in order to study efficiency of the ID/////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
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


/*
        const xAOD::PhotonContainer* photons = 0;
        ANA_CHECK(event->retrieve(photons, "PhotonCollection"));
        xAOD::PhotonContainer::const_iterator begin_ph = photons->begin();
	    xAOD::PhotonContainer::const_iterator end_ph = photons->end();
	     for(;begin_ph!=end_ph;++begin_ph)
       {
            xAOD::Photon* ph = (*begin_ph);

             if(ph->caloCluster()->pt()<20)continue;


        unsigned long int _m_isEM_Tight = 999999;
	    _m_isEM_Tight = myPhTight->IsemValue();
	    bool _b_IsNewTight = myPhTight->accept(begin_ph);
*/


/*
        ///////////////////////////////////////////////////////////////////////////////////////////////////////
        //Lets make the tree corresponding to the pairs of muons needed for the reconstruction efficiency of (ID&MS)/all ID - Probability to reconstruct Medium muon
        ///////////////////////////////////////////////////////////////////////////////////////////////////////
        xAOD::MuonContainer::const_iterator begin_new1 = m_newMuons->begin();
        xAOD::MuonContainer::const_iterator end_new1 = m_newMuons->end();
        for(;begin_new1!=end_new1;++begin_new1)
           {
             cout<<"Im in the first loop in tree2 "<<endl;
             m_passTightID  = false; if (m_muonSelection1->accept(*begin_new1)) m_passTightID = true; //based on eta quality and ID cuts - tight muons
	         m_passMediumID = false; if (m_muonSelection2->accept(*begin_new1)) m_passMediumID = true; //medium muons
	         m_passIDcuts   = false; if (m_muonSelection2->passedIDCuts(**begin_new1)) m_passIDcuts = true;
             //Lets do some cuts:
             if(m_HLT_mu15 == false)continue; //the event shall pass the muon15 trigger
             bool Muon_MatchdR_mu15      = m_trigMuonMatching->match((*begin_new1),"HLT_mu15",0.1);
             if (Muon_MatchdR_mu15 == false)continue; //the muon shall be matched to the trigger
             m_MatchdR_mu15          = m_trigMuonMatching->minDelR((*begin_new1),"HLT_mu15",1);
	         //if (!istag) continue;
	         if((*begin_new1)->pt()*0.001<17)continue; //muon pT shall be more then the 17 GeV - Plato for HLT_mu15
	         if(abs((*begin_new1)->eta())>2.4)continue;//the trigger points till 2.4
             //Muon kinematic
	         m_muonPt       = (*begin_new1)->pt()*0.001;
	         m_muonPtcorr   = (*begin_new1)->auxdata< float >( "ptcorr" );
	         m_muonEta      = (*begin_new1)->eta();
	         m_muonPhi      = (*begin_new1)->phi();
	         m_muonCharge   = (*begin_new1)->charge();
	         m_muonQuality  = (*begin_new1)->quality();
	         m_muonType     = (*begin_new1)->muonType();
	         m_isCombined_2   = false; if((*begin_new1)->muonType() == xAOD::Muon_v1::Combined) m_isCombined_2=true;
	         m_isSegTag_2     = false; if((*begin_new1)->muonType() == xAOD::Muon_v1::SegmentTagged) m_isSegTag_2=true;
	         m_isCaloTag_2    = false; if((*begin_new1)->muonType() == xAOD::Muon_v1::CaloTagged) m_isCaloTag_2=true;
	         m_isStandAl_2    = false; if((*begin_new1)->muonType() == xAOD::Muon_v1::MuonStandAlone) m_isStandAl_2=true;
	         //m_muonEloss    = muon->floatParameter(xAOD::muon::EnergyLoss);

	        //Muon ID track information
            const xAOD::TrackParticle* Track = (*begin_new1)->trackParticle( xAOD::Muon::InnerDetectorTrackParticle );
             if (Track)
            {

               m_Muonid_d0      =  Track->d0()    ;
               m_Muonid_d0sig   =  sqrt(Track->definingParametersCovMatrix()(0,0));
               m_Muonid_z0      =  Track->z0()    ;
               m_Muonid_z0sig   =  sqrt(Track->definingParametersCovMatrix()(1,1));
               m_Muonid_theta   =  Track->theta() ;
               m_Muonid_phi     =  Track->phi()   ;
               m_Muonid_pt      =  Track->pt()*0.001;
               m_Muonid_charge  =  Track->charge();
               m_Muonid_eta     =  Track->eta()   ;
               m_Muonid_ez0sin  = sqrt((Track->definingParametersCovMatrix()(1,1))*pow(sin(m_Muonid_theta),2) +
				                            (Track->definingParametersCovMatrix()(3,3))*pow(m_Muonid_z0*cos(m_Muonid_theta),2)+
				                             2*(Track->definingParametersCovMatrix()(1,3))*fabs(sin(m_Muonid_theta)*m_Muonid_z0*cos(m_Muonid_theta)) );


				nPixelHits_calo        = Track->auxdata< unsigned char >("numberOfPixelHits");
	            nPixelDeadSensors_calo = Track->auxdata< unsigned char >("numberOfPixelDeadSensors");
	            nSCTHits_calo          = Track->auxdata< unsigned char >("numberOfSCTHits");
	            nSCTDeadSensors_calo   = Track->auxdata< unsigned char >("numberOfSCTDeadSensors");
	            nPixelHoles_calo       = Track->auxdata< unsigned char >("numberOfPixelHoles");
	            nSCTHoles_calo         = Track->auxdata< unsigned char >("numberOfSCTHoles");
	            nTRTHits_calo          = Track->auxdata< unsigned char >("numberOfTRTHits");
	            nTRTOutliers_calo      = Track->auxdata< unsigned char >("numberOfTRTOutliers");
	            nIBLHits_calo          = Track->auxdata< unsigned char >("numberOfInnermostPixelLayerHits");
	            nBLHits_calo           = Track->auxdata< unsigned char >("numberOfNextToInnermostPixelLayerHits");
	            expIBLHits_calo        = Track->auxdata< unsigned char >("expectInnermostPixelLayerHit");
	            expBLHits_calo         = Track->auxdata< unsigned char >("expectNextToInnermostPixelLayerHit");

	            m_track_passIDcuts   = false; if (m_muonSelection->passedIDCuts(*Track)) m_track_passIDcuts = true;

            }
*/
/*
             //Probe track
             cout<<"The pt in the tree2 loop is "<<m_muonPt<<endl;
             cout<<"The corrected pt in tree2 loop is "<<m_muonPtcorr<<endl;
             for (const auto& trk : *Recotracks)
                {

                  //Lets do some track cuts
                  //Pt cut
                  if (trk->pt()*0.001<0.4) continue;
                  //eta cut
                  float eta=trk->eta();
	              if(fabs(eta)>2.5) continue;
                  //theta cut
                  float theta=trk->theta();//to make sure that there are no tracks pointin gdirectcly to the beam pipe where rapidity is infinite
	              if(theta==0.0) continue;
                  //Hits and holes
	              //int nPixhits   = trk->auxdata< unsigned char >("numberOfPixelHits") +   trk->auxdata< unsigned char >("numberOfPixelDeadSensors") ;
                  //int nSCThits   = trk->auxdata< unsigned char >("numberOfSCTHits")   +  trk->auxdata< unsigned char >("numberOfSCTDeadSensors");
	              //int nPixholes  = trk->auxdata< unsigned char >("numberOfPixelHoles");
	              //int nSCTholes  = trk->auxdata< unsigned char >("numberOfSCTHoles");
                  //int nIBLHits   = trk->auxdata< unsigned char >("numberOfInnermostPixelLayerHits");
                  ///int nBLHits    = trk->auxdata< unsigned char >("numberOfNextToInnermostPixelLayerHits");
                  //int expIBLHits = trk->auxdata< unsigned char >("expectInnermostPixelLayerHit");
                  //int expBLHits  = trk->auxdata< unsigned char >("expectNextToInnermostPixelLayerHit");
                  //these are tight track cuts
	              //if((nPixhits+nSCThits)<8)continue;
	              //if((nPixhits+nSCThits)<10 && fabs(eta)>1.65)continue;
	              //if((nPixholes+nSCTholes)>=2)continue;
	              //if( nPixholes >=1)continue;
	              //if( (nIBLHits + nBLHits) <1 && expIBLHits && expBLHits)continue;
	              //Letsfill in all the information related to the tracks
	              //hits
	              m_Track_nPixelHits        = trk->auxdata< unsigned char >("numberOfPixelHits");
	              m_Track_nPixelDeadSensors = trk->auxdata< unsigned char >("numberOfPixelDeadSensors");
	              m_Track_nSCTHits          = trk->auxdata< unsigned char >("numberOfSCTHits");
	              m_Track_nSCTDeadSensors   = trk->auxdata< unsigned char >("numberOfSCTDeadSensors");
	              m_Track_nPixelHoles       = trk->auxdata< unsigned char >("numberOfPixelHoles");
	              m_Track_nSCTHoles         = trk->auxdata< unsigned char >("numberOfSCTHoles");
	              m_Track_nTRTHits          = trk->auxdata< unsigned char >("numberOfTRTHits");
	              m_Track_nTRTOutliers      = trk->auxdata< unsigned char >("numberOfTRTOutliers");
	              m_Track_nIBLHits          = trk->auxdata< unsigned char >("numberOfInnermostPixelLayerHits");
	              m_Track_nBLHits           = trk->auxdata< unsigned char >("numberOfNextToInnermostPixelLayerHits");
	              m_Track_expIBLHits        = trk->auxdata< unsigned char >("expectInnermostPixelLayerHit");
	              m_Track_expBLHits         = trk->auxdata< unsigned char >("expectNextToInnermostPixelLayerHit");
                  //kinematics
                  m_Track_pt                = trk->pt()*0.001;
                  m_Track_eta               = trk->eta();
                  m_Track_phi               = trk->phi();
                  m_Track_theta             = trk->theta();
                  m_Track_charge            = trk->charge();
                  m_Track_d0                = trk->d0();
                  m_Track_z0_calc           = (trk->z0()+trk->vz()-m_vertex_z)*sin(m_Track_theta);//with respect to collision point//with respect to the vertex
                  m_Track_z0                = trk->z0();//with respect to beam line
                  m_Track_ed0               = sqrt(trk->definingParametersCovMatrix()(0,0));//d0 significance
                  m_Track_ez0               = sqrt(trk->definingParametersCovMatrix()(1,1));//z0 significance
                  m_Track_ez0sin            = sqrt((trk->definingParametersCovMatrix()(1,1))*pow(sin(m_Track_theta),2) +
				                            (trk->definingParametersCovMatrix()(3,3))*pow(m_Track_z0*cos(m_Track_theta),2)+
				                             2*(trk->definingParametersCovMatrix()(1,3))*fabs(sin(m_Track_theta)*m_Track_z0*cos(m_Track_theta)) );


*/           //for the calo tagged muons

/*
             xAOD::MuonContainer::const_iterator begin_new2 = m_newMuons->begin();
             xAOD::MuonContainer::const_iterator end_new2 = m_newMuons->end();
             cout<<"Im in the second loop in tree2 "<<endl;
             for(;begin_new2!=end_new2;++begin_new2)
                {

                    //check if this is a calo tagged muon
                    m_isCaloTag_c     = false; if((*begin_new2)->isAuthor(xAOD::Muon::CaloTag))m_isCaloTag_c=true; //to require the author you will peak the muons which are reconstructed by the calo algorithms - thats what you need
                    //in the muon type - the muons are already put in the hierarchy where the higher the hierarchy the better the muon - so calo tagged muons are left with almost nothing - the works types
	                m_passed_calo_c   = false; if(m_muonSelection->passedCaloTagQuality(**begin_new2))   m_passed_calo_c = true;//To study Calo tagged
                    m_passed_IDcuts_c = false; if(m_muonSelection->passedIDCuts(**begin_new2))           m_passed_IDcuts_c = true;


                    m_caloPt       = (*begin_new2)->pt()*0.001;
	                m_caloPtcorr   = (*begin_new2)->auxdata< float >( "ptcorr" );
	                m_caloEta      = (*begin_new2)->eta();
	                m_caloPhi      = (*begin_new2)->phi();
	                m_caloCharge   = (*begin_new2)->charge();


                    //tracks from calo
                    const xAOD::TrackParticle* Calo_Track = (*begin_new2)->trackParticle( xAOD::Muon::InnerDetectorTrackParticle );
                    if (Calo_Track)
                       {


                            m_Muonidcalo_d0      =  Calo_Track->d0()    ;
                            m_Muonidcalo_d0sig   =  sqrt(Calo_Track->definingParametersCovMatrix()(0,0));
               				m_Muonidcalo_z0      =  Calo_Track->z0()    ;
               				m_Muonidcalo_z0sig   =  sqrt(Calo_Track->definingParametersCovMatrix()(1,1));
               				m_Muonidcalo_theta   =  Calo_Track->theta() ;
               				m_Muonidcalo_phi     =  Calo_Track->phi()   ;
               				m_Muonidcalo_pt      =  Calo_Track->pt()*0.001;
               				m_Muonidcalo_charge  =  Calo_Track->charge();
               				m_Muonidcalo_eta     =  Calo_Track->eta()   ;
               				m_Muonidcalo_ez0sin  = sqrt((Calo_Track->definingParametersCovMatrix()(1,1))*pow(sin(m_Muonidcalo_theta),2) +
				                                   (Calo_Track->definingParametersCovMatrix()(3,3))*pow(m_Muonidcalo_z0*cos(m_Muonidcalo_theta),2)+
				                                    2*(Calo_Track->definingParametersCovMatrix()(1,3))*fabs(sin(m_Muonidcalo_theta)*m_Muonidcalo_z0*cos(m_Muonidcalo_theta)) );




				            m_Muonidcalo_nPixelHits        = Calo_Track->auxdata< unsigned char >("numberOfPixelHits");
	            			m_Muonidcalo_nPixelDeadSensors = Calo_Track->auxdata< unsigned char >("numberOfPixelDeadSensors");
	            			m_Muonidcalo_nSCTHits          = Calo_Track->auxdata< unsigned char >("numberOfSCTHits");
	            			m_Muonidcalo_nSCTDeadSensors   = Calo_Track->auxdata< unsigned char >("numberOfSCTDeadSensors");
	            			m_Muonidcalo_nPixelHoles       = Calo_Track->auxdata< unsigned char >("numberOfPixelHoles");
	            			m_Muonidcalo_nSCTHoles         = Calo_Track->auxdata< unsigned char >("numberOfSCTHoles");
	            			m_Muonidcalo_nTRTHits          = Calo_Track->auxdata< unsigned char >("numberOfTRTHits");
	            			m_Muonidcalo_nTRTOutliers      = Calo_Track->auxdata< unsigned char >("numberOfTRTOutliers");
	            			m_Muonidcalo_nIBLHits          = Calo_Track->auxdata< unsigned char >("numberOfInnermostPixelLayerHits");
	            			m_Muonidcalo_nBLHits           = Calo_Track->auxdata< unsigned char >("numberOfNextToInnermostPixelLayerHits");
	            			m_Muonidcalo_expIBLHits        = Calo_Track->auxdata< unsigned char >("expectInnermostPixelLayerHit");
	            			m_Muonidcalo_expBLHits         = Calo_Track->auxdata< unsigned char >("expectNextToInnermostPixelLayerHit");

	            			m_calotrack_passIDcuts   = false; if (m_muonSelection->passedIDCuts(*Calo_Track)) m_calotrack_passIDcuts = true;
                        }



                     //Make the pair of a tag muon and a tight track
                  	TLorentzVector mom1,mom2;
 				  	mom1.SetPtEtaPhiM(m_muonPt,m_muonEta,m_muonPhi,0.106);
		 	      	mom2.SetPtEtaPhiM(m_caloPt,m_caloEta,m_caloPhi,0.106);
		 	      	TLorentzVector mom = mom1+mom2;
                  	//mass cut
		 	      	if (mom.M() < 60|| mom.M() > 120) continue; //Z mass
		 	      	pairpt   = mom.Pt();
		 	      	pairy    = mom.Rapidity();
		 	      	paireta  = mom.PseudoRapidity();
		 	      	pairphi  = mom.Phi();
		 	      	pairm    = mom.M();



                  	xAOD::MuonContainer::const_iterator begin_calo = m_newMuons->begin();
                  	xAOD::MuonContainer::const_iterator end_calo = m_newMuons->end();
                  	for(;begin_calo!=end_calo;++begin_calo)
                     	{   cout<<"Checking the matching in tree2 "<<endl;
                         	m_passMediumID_c = false;    if (m_muonSelection2->accept(*begin_calo)) m_passMediumID_c = true;
                         	m_passTightID_c =  false;    if (m_muonSelection1->accept(*begin_calo)) m_passTightID_c = true;

                            deltaR_med  = +99999.9;
                            deltaR_tight  = +99999.9;
                         	if(m_passMediumID_c)
                           		{
                                	float phi=(*begin_calo)->phi();
                                	float deltaphi = fabs (m_caloPhi - phi);
                                	if (deltaphi > TMath::Pi()) deltaphi = 2*TMath::Pi() - deltaphi;
                                	double dr = sqrt( pow(deltaphi,2) + pow ((m_caloEta - (*begin_calo)->eta()),2) ); //to match
                                	if(dr<deltaR_med) deltaR_med=dr;

                           		}

                        	if(m_passTightID_c)
                           		{
                                	float phi=(*begin_calo)->phi();
                                	float deltaphi = fabs (m_caloPhi - phi);
                                	if (deltaphi > TMath::Pi()) deltaphi = 2*TMath::Pi() - deltaphi;
                                	double dr = sqrt( pow(deltaphi,2) + pow ((m_caloEta - (*begin_calo)->eta()),2) ); //to match
                                	if(dr<deltaR_tight) deltaR_tight=dr;
                           		}


                      	}//end of 3d loop to match the calo
                 	cout<<"Filling in the second tree "<<endl;
                 	tree2->Fill();
                }//end of the second loop
            }//end of the first loop

*/
/*

                  //Make the pair of a tag muon and a tight track
                  TLorentzVector mom1,mom2;
 				  mom1.SetPtEtaPhiM(m_muonPt,m_muonEta,m_muonPhi,0.106);
		 	      mom2.SetPtEtaPhiM(m_Track_pt,m_Track_eta,m_Track_phi,0.106);
		 	      TLorentzVector mom = mom1+mom2;
                  //mass cut
		 	      if (mom.M() < 60|| mom.M() > 120) continue; //Z mass
		 	      pairpt   = mom.Pt();
		 	      pairy    = mom.Rapidity();
		 	      paireta  = mom.PseudoRapidity();
		 	      pairphi  = mom.Phi();
		 	      pairm    = mom.M();


                  m_isMuon_medium = false;
                  m_isMuon_loose  = false;
                  m_isMuon_tight  = false;
                  //for (auto MuonZ: *(m_newMuons->first))
                  xAOD::MuonContainer::const_iterator begin_new2 = m_newMuons->begin();
                  xAOD::MuonContainer::const_iterator end_new2 = m_newMuons->end();
                  for(;begin_new2!=end_new2;++begin_new2)
                        {
                        // check muon quality here
                        ElementLink< xAOD::TrackParticleContainer > trkmulink = (*begin_new2)->auxdata<ElementLink< xAOD::TrackParticleContainer > >("inDetTrackParticleLink");
                        if(trkmulink.isValid())
                            {
                              if(TMath::Abs((*trkmulink)->pt() - trk->pt())<1e-6 && TMath::Abs((*trkmulink)->eta() - trk->eta())<1e-6 && TMath::Abs((*trkmulink)->phi() - trk->phi())<1e-6)
                                {
                                  if(m_muonSelection1->accept((*begin_new2)))m_isMuon_tight  = true;
                                  if(m_muonSelection2->accept((*begin_new2)))m_isMuon_medium = true;
                                  if(m_muonSelection3->accept((*begin_new2)))m_isMuon_loose  = true;
                                  break;
                                }
                            }
                        }
*/
                // cout<<"Filling in the second tree "<<endl;
               //  tree2->Fill();

              //  }// end track loop
         //   }// End of the muon loop
/*
        ////////////////////////////////////////////////////////////////////////////////////////////////////
        //Lets make the tree corresponding to the pairs of muons needed for the reconstruction efficiency ID - (MS && ID)/all MS
        ////////////////////////////////////////////////////////////////////////////////////////////////////
        xAOD::MuonContainer::const_iterator begin_new3 = m_newMuons->begin();
        xAOD::MuonContainer::const_iterator end_new3 = m_newMuons->end();
        for(;begin_new3!=end_new3;++begin_new3)
           {
           cout<<"Im in the first loop in tree3 "<<endl;
             //Lets do some cuts:
             m_passTightID_3  = false; if (m_muonSelection1->accept((*begin_new3))) m_passTightID_3 = true;
	         m_passMediumID_3 = false; if (m_muonSelection2->accept((*begin_new3))) m_passMediumID_3 = true; //medium muons
	         m_passIDcuts_3   = false; if (m_muonSelection2->passedIDCuts(**begin_new3)) m_passIDcuts_3 = true;
             //event mast pass the HLT_mu15 trrigger
             if (m_HLT_mu15 == false)continue;
             //1st muon mathed to the trigger
             bool Muon_MatchdR_mu15 = m_trigMuonMatching->match((*begin_new3),"HLT_mu15",0.1);
             if (Muon_MatchdR_mu15 == 0)continue;
             //record trrigger decision
             MatchdR_mu15_3          = m_trigMuonMatching->minDelR((*begin_new3),"HLT_mu15",1);
             //Cuts on the tag muon
	         if((*begin_new3)->pt()*0.001<17)continue;
	         if(abs((*begin_new3)->eta())>2.4)continue;

	         m_muonPt_3       = (*begin_new3)->pt()*0.001;
	         m_muonPtcorr_3   = (*begin_new3)->auxdata< float >( "ptcorr" );
	         m_muonEta_3      = (*begin_new3)->eta();
	         m_muonPhi_3      = (*begin_new3)->phi();
	         m_muonCharge_3   = (*begin_new3)->charge();
	         m_muonQuality_3  = (*begin_new3)->quality();
	         m_muonType_3     = (*begin_new3)->muonType();
	         m_isCombined_3   = false; if((*begin_new3)->muonType() == xAOD::Muon_v1::Combined) m_isCombined_3=true;
	         m_isSegTag_3     = false; if((*begin_new3)->muonType() == xAOD::Muon_v1::SegmentTagged) m_isSegTag_3=true;
	         m_isCaloTag_3    = false; if((*begin_new3)->muonType() == xAOD::Muon_v1::CaloTagged) m_isCaloTag_3=true;
	         m_isStandAl_3    = false; if((*begin_new3)->muonType() == xAOD::Muon_v1::MuonStandAlone) m_isStandAl_3=true;
             //Muon ID track information
             const xAOD::TrackParticle* Track = (*begin_new3)->trackParticle( xAOD::Muon::InnerDetectorTrackParticle );
             if (Track)
              {

               m_Muonid_d0_3      =  Track->d0()    ;
               m_Muonid_d0sig_3   =  sqrt(Track->definingParametersCovMatrix()(0,0));
               m_Muonid_z0_3      =  Track->z0()    ;
               m_Muonid_z0sig_3   =  sqrt(Track->definingParametersCovMatrix()(1,1));
               m_Muonid_theta_3   =  Track->theta() ;
               m_Muonid_phi_3     =  Track->phi()   ;
               m_Muonid_pt_3      =  Track->pt()*0.001;
               m_Muonid_charge_3  =  Track->charge();
               m_Muonid_eta_3     =  Track->eta()   ;
               m_Muonid_ez0sin_3  = sqrt((Track->definingParametersCovMatrix()(1,1))*pow(sin(m_Muonid_theta_3),2) +
				                            (Track->definingParametersCovMatrix()(3,3))*pow(m_Muonid_z0_3*cos(m_Muonid_theta_3),2)+
				                             2*(Track->definingParametersCovMatrix()(1,3))*fabs(sin(m_Muonid_theta_3)*m_Muonid_z0_3*cos(m_Muonid_theta_3)) );

			   m_Muonidtrack_passIDcuts   = false; if (m_muonSelection->passedIDCuts(*Track)) m_Muonidtrack_passIDcuts = true;
              }

             xAOD::MuonContainer::const_iterator begin_new4 = m_newMuons->begin();
             xAOD::MuonContainer::const_iterator end_new4 = m_newMuons->end();
             for(;begin_new4!=end_new4;++begin_new4)
             //for   (auto muon2: *(m_newMuons.first))
                {
                    cout<<"Im in the second loop in tree3 "<<endl;
	               	//passed         = false; for the calo tagged muons (form the validation)
	               	m_passTightID_4  = false; if (m_muonSelection1->accept((*begin_new4))) m_passTightID_4 = true; //based on eta quality and ID cuts - tight muons
	               	m_passMediumID_4 = false; if (m_muonSelection2->accept((*begin_new4))) m_passMediumID_4 = true; //medium muons
	               	m_passLooseID_4  = false; if (m_muonSelection3->accept((*begin_new4))) m_passLooseID_4 = true; //loose muons
                   	m_isCaloTag_4    = false; if((*begin_new4)->muonType() == xAOD::Muon_v1::CaloTagged) m_isCaloTag_4=true;//To study the calo tagged
                   	m_isCombined_4   = false; if((*begin_new4)->muonType() == xAOD::Muon_v1::Combined) m_isCombined_4=true;
	               	//if(m_muonSelection3->passedCaloTagQuality(*muon2) && muon2->isAuthor(xAOD::Muon::CaloTag) ) passed = true;//To study Calo tagged
    			   	//muon kinematics
	               	m_muon2_pt       = (*begin_new4)->pt()*0.001;
	               	m_muon2_ptcorr   = (*begin_new4)->auxdata< float >( "ptcorr" );
	               	m_muon2_eta      = (*begin_new4)->eta();
	               	m_muon2_phi      = (*begin_new4)->phi();
	               	m_muon2_charge   = (*begin_new4)->charge();
                   	//MS tracks
                   	//const xAOD::TrackParticle* MeTracks = (*begin_new4)->trackParticle( xAOD::Muon::ExtrapolatedMuonSpectrometerTrackParticle );//For the ID effic
                   	const xAOD::TrackParticle* MeTracks = (*begin_new4)->trackParticle( xAOD::Muon::MuonSpectrometerTrackParticle);//For calotagged muons
                   	//if(MeTracks) //for the Calotag muons
                    if(!MeTracks)continue;//for the extrapolatedmuonspectrometer track
                    m_ME_Track_pt =MeTracks->pt()*0.001;
                   	m_ME_Track_eta =MeTracks->eta();
                    m_ME_Track_phi =MeTracks->phi();
                    m_ME_Track_charge =MeTracks->charge();
                    m_ME_Track_d0 =MeTracks->d0();
                    m_ME_Track_z0 =MeTracks->z0();
                    m_ME_Track_theta =MeTracks->theta();
                    m_ME_Track_ed0= sqrt(MeTracks->definingParametersCovMatrix()(0,0));
                    //Define the lorentz vector
            		TLorentzVector mom1,mom2;
					mom1.SetPtEtaPhiM(m_muonPt_3,m_muonEta_3,m_muonPhi_3,0.106);
		 	       	mom2.SetPtEtaPhiM(m_ME_Track_pt,m_ME_Track_eta,m_ME_Track_phi,0.106);
		 	       	TLorentzVector mom3 = mom1+mom2;
		 	       	//mass cut
		 	       	if (mom3.M() < 60|| mom3.M() > 120) continue; //Z mass
		 	      	pairpt_3  = mom3.Pt();
		 	     	pairy_3   = mom3.Rapidity();
		 	     	paireta_3 = mom3.PseudoRapidity();
		 	     	pairphi_3 = mom3.Phi();
		 	     	pairm_3   = mom3.M();
*/
/*
                   	m_istrack      = 0; //using the something else gods knows what is the difference between them, but clearlt the pts are slightly different
	               	const  xAOD::TrackParticle* trkmulink = (*begin_new4)->trackParticle ( xAOD::Muon::InnerDetectorTrackParticle);
                   	if(trkmulink)
                     {  bool keep = true;
                        //Lets do some track cuts
                        //pt cut
                        if (trkmulink->pt()*0.001<0.4) keep=false; //Pt cut
  						//eta cut
                        float eta=trkmulink->eta();
	                    if(fabs(eta)>2.5) keep=false; //eta cut
						//theta cut - no tracks along the beam line
						float theta=trkmulink->theta();
	              		if(theta==0.0) keep=false;
 						//hits
                     	nPixelHits_3        = trkmulink->auxdata< unsigned char >("numberOfPixelHits");
	                    nPixelDeadSensors_3 = trkmulink->auxdata< unsigned char >("numberOfPixelDeadSensors");
	                    nSCTHits_3         = trkmulink->auxdata< unsigned char >("numberOfSCTHits");
	                    nSCTDeadSensors_3   = trkmulink->auxdata< unsigned char >("numberOfSCTDeadSensors");
	              		nPixelHoles_3       = trkmulink->auxdata< unsigned char >("numberOfPixelHoles");
	              		nSCTHoles_3         = trkmulink->auxdata< unsigned char >("numberOfSCTHoles");
	              		nTRTHits_3          = trkmulink->auxdata< unsigned char >("numberOfTRTHits");
	              		nTRTOutliers_3      = trkmulink->auxdata< unsigned char >("numberOfTRTOutliers");
	              		nIBLHits_3          = trkmulink->auxdata< unsigned char >("numberOfInnermostPixelLayerHits");
	              		nBLHits_3           = trkmulink->auxdata< unsigned char >("numberOfNextToInnermostPixelLayerHits");
	              		expIBLHits_3        = trkmulink->auxdata< unsigned char >("expectInnermostPixelLayerHit");
	              		expBLHits_3         = trkmulink->auxdata< unsigned char >("expectNextToInnermostPixelLayerHit");

	              		//if((nPixhits+nSCThits)<8)keep=false;
	              		//if((nPixhits+nSCThits)<10 && fabs(eta)>1.65)keep=false;
	              		//if((nPixholes+nSCTholes)>=2)continue;
	              		//if( nPixholes >=1)keep=false;
	                    ///if( (nIBLHits + nBLHits) <1 && expIBLHits && expBLHits)keep=false;
                        if(keep) m_istrack=1; //Here is my matched track
                        float phi=trkmulink->phi();
                        float deltaphi = fabs (m_ME_Track_phi - phi);
                        if (deltaphi > TMath::Pi()) deltaphi = 2*TMath::Pi() - deltaphi;
                        m_dR=sqrt( pow(deltaphi,2) + pow ((m_ME_Track_eta - eta),2) ); //to match in the code
	                    //some variables to store
						d0_trk     =trkmulink->d0();
	              		ed0_trk    =sqrt(trkmulink->definingParametersCovMatrix()(0,0));
	              		ez0_trk    =sqrt(trkmulink->definingParametersCovMatrix()(1,1));
                        z0_trk    =trkmulink->z0();
                        theta_trk  =trkmulink->theta();
                        pt_trk     =trkmulink->pt();
                        eta_trk	   =trkmulink->eta();
                        phi_trk    =trkmulink->phi();
                      }//end of track link if statement
                    cout<<"Filling in the third tree "<<endl;
*/                  //Checking for calo tagged muons
/*
                     xAOD::MuonContainer::const_iterator begin_cal = m_newMuons->begin();
                     xAOD::MuonContainer::const_iterator end_cal = m_newMuons->end();
                     is_calo = false;
                     is_ID   = false;
                     pass_calo = false;
                     dR      = 999999.9;
                     for(;begin_cal!=end_cal;++begin_cal)
                        {
                             cout<<"Checking the matching in tree3 "<<endl;
                             if(m_muonSelection->passedCaloTagQuality(**begin_cal)) pass_calo = true;
                             if((*begin_cal)->isAuthor(xAOD::Muon::CaloTag)) is_calo = true;
                             //if(m_muonSelection->passedCaloTagQuality(**begin_cal) && (*begin_cal)->muonType() == xAOD::Muon_v1::CaloTagged) is_calo=true;//To study the calo tagged
                             if(m_muonSelection->passedIDCuts(**begin_cal)) is_ID = true;

                             float phi=(*begin_cal)->phi();
                             float deltaphi = fabs (m_ME_Track_phi - phi);
                             if (deltaphi > TMath::Pi()) deltaphi = 2*TMath::Pi() - deltaphi;
                             double dr = sqrt( pow(deltaphi,2) + pow ((m_ME_Track_eta - (*begin_cal)->eta()),2) );
                             if(dr<dR) dR=dr;

                        }//the matching loop

                    tree3->Fill();

                }// end muon2 loop


            }//end of the first loop
*/

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





/*
                   if(trkmulinks.isValid())
                     {          bool keeps = true;
                                //Lets do some track cuts
                                if ((*trkmulinks)->pt()*0.001<10) keeps=false; //Pt cut

                                float etas=(*trkmulinks)->eta();
	                            if(fabs(etas)>2.5) keeps=false; //eta cut

	                            float phis=(*trkmulinks)->phi();


                                float thetas=(*trkmulinks)->theta();
	              				if(thetas==0.0) keeps=false;


                  				//float z0=(trk->z0())*sin(theta);
	              				int nPixhitss   = (*trkmulinks)->auxdata< unsigned char >("numberOfPixelHits") +   (*trkmulinks)->auxdata< unsigned char >("numberOfPixelDeadSensors") ;
                  				int nSCThitss   = (*trkmulinks)->auxdata< unsigned char >("numberOfSCTHits")   +   (*trkmulinks)->auxdata< unsigned char >("numberOfSCTDeadSensors");
	              				int nPixholess  = (*trkmulinks)->auxdata< unsigned char >("numberOfPixelHoles");
	              				//int nSCTholes  = (*trkmulink)->auxdata< unsigned char >("numberOfSCTHoles");
                  				int nIBLHitss   = (*trkmulinks)->auxdata< unsigned char >("numberOfInnermostPixelLayerHits");
                  				int nBLHitss    = (*trkmulinks)->auxdata< unsigned char >("numberOfNextToInnermostPixelLayerHits");
                  				int expIBLHitss = (*trkmulinks)->auxdata< unsigned char >("expectInnermostPixelLayerHit");
                     		    int expBLHitss  = (*trkmulinks)->auxdata< unsigned char >("expectNextToInnermostPixelLayerHit");

	              				if((nPixhitss+nSCThitss)<8)keeps=false;
	              				if((nPixhitss+nSCThitss)<10 && fabs(etas)>1.65)keeps=false;
	              				//if((nPixholes+nSCTholes)>=2)continue;
	              				if( nPixholess >=1)keeps=false;
	                            if( (nIBLHitss + nBLHitss) <1 && expIBLHitss && expBLHitss)keeps=false;


	                            if(keeps) m_istrack_old=1; //Here is my matched track

                                float deltaphis = fabs (m_ME_Track_phi - phis);
                                if (deltaphis > TMath::Pi()) deltaphis = 2*TMath::Pi() - deltaphis;
                                m_dR_old=sqrt( pow(deltaphis,2) + pow ((m_ME_Track_eta - etas),2) );

                                d0_old=(*trkmulinks)->d0();
	              				ed0_old=sqrt((*trkmulinks)->definingParametersCovMatrix()(0,0));
                                z0_old    =(*trkmulinks)->z0();
                                theta_old =(*trkmulinks)->theta();



                     }

*/


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



