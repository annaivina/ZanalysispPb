#ifndef ZanalysispPb_Zreco_H
#define ZanalysispPb_Zreco_H

#include <EventLoop/Algorithm.h>
#include "GoodRunsLists/GoodRunsListSelectionTool.h"
#include "TrigConfxAOD/xAODConfigTool.h"
#include "TrigDecisionTool/TrigDecisionTool.h"
#include "TrigMuonMatching/TrigMuonMatching.h"
#include "MuonSelectorTools/MuonSelectionTool.h"
#include "MuonMomentumCorrections/MuonCalibrationAndSmearingTool.h"
#include "TriggerMatchingTool/MatchingTool.h"
#include "TrigMuonMatching/TrigMuonMatching.h"

//to use the fcal tool
#include "HIEventUtils/HIPileupTool.h"
#include "ZdcAnalysis/ZdcAnalysisTool.h"
#include "HIEventUtils/HICentralityTool.h"

#include <TTree.h>
#include <vector>
#include <string>


// electron selectors
#include "ElectronPhotonSelectorTools/AsgElectronIsEMSelector.h"
#include "ElectronPhotonSelectorTools/AsgElectronLikelihoodTool.h"
#include "ElectronPhotonSelectorTools/AsgForwardElectronIsEMSelector.h"

#include "GoodRunsLists/GoodRunsListSelectionTool.h"
#include "TrigConfxAOD/xAODConfigTool.h"

//calibration
#include "ElectronPhotonFourMomentumCorrection/EgammaCalibrationAndSmearingTool.h"

//electrons
#include "xAODEgamma/ElectronContainer.h"
#include "xAODEgamma/ElectronAuxContainer.h"
#include "xAODEgamma/EgammaxAODHelpers.h"


//matching
#include "TrigEgammaMatchingTool/TrigEgammaMatchingTool.h"

//calibration
#include "ElectronPhotonFourMomentumCorrection/EgammaCalibrationAndSmearingTool.h"


using namespace std;



class Zreco : public EL::Algorithm
{
  // put your configuration variables here as public variables.
  // that way they can be set directly from CINT and python.
public:
  // float cutValue;



  // variables that don't get filled at submission time should be
  // protected from being send from the submission node to the worker
  // node (done by the //!)
public:

    int m_eventCounter;                     //!
    std::string outputName;
    TTree *tree;                            //!
    TTree *tree2;                           //!
    TTree *tree3;                           //!


    //for tree
    int  m_RunNumber; 						//!
    int  m_EventNumber; 					//!
    int  m_LumiBlock; 						//!
    int  m_mcChannelNumber; 				//!
    double  m_averageIntPerXing;            //!
    double  m_actualIntPerXing;            //!
    double  m_FCal_Et; 						//!
    int m_is_pileup;						//!

    int m_bcid;								//!
    int is_good_bcid;						//!

     // RoI besed single muon chains
    bool  m_HLT_mu4; 			            //!
    bool  m_HLT_mu4_nomucomb; 				//!
    bool  m_HLT_mu4_bJpsi_Trkloose; 		//!
    bool  m_HLT_mu6; 						//!
    bool  m_HLT_mu6_nomucomb; 				//!
    bool  m_HLT_mu6_bJpsi_Trkloose; 		//!
    bool  m_HLT_mu8; 						//!
    bool  m_HLT_mu8_bJpsi_Trkloose; 		//!
    bool  m_HLT_mu10; 						//!
    bool  m_HLT_mu15; 						//!
    bool  m_HLT_mu15_L1MU10; 				//!
    bool  m_HLT_mu10_L1MU6; 				//!
    bool  m_HLT_mu15_L1MU6; 				//!
    bool  m_L1_MU4;							//!
    bool  m_L1_MU6;							//!
    bool  m_L1_MU10;						//!
    bool  m_L1_MU15;						//!
    bool  m_HLT_2mu4; 						//!
    bool  m_HLT_2mu4_nomucomb; 				//!
    bool  m_HLT_2mu6; 						//!
    bool  m_HLT_2mu6_nomucomb; 				//!
    float m_vertex_x;						//!
    float m_vertex_y;						//!
    float m_vertex_z;						//!

    bool  m_HLT_mb_sptrk;        			//!
    bool  m_record;				 			//!


    //tree2
    float m_muonPt;    						//!
	float m_muonEta;  						//!
	float m_muonPhi;  						//!
	float m_muonCharge;  					//!
	float m_muonQuality; 					//!
    float m_muonType;    					//!
	float m_muonEloss;  					//!

	bool m_isCombined_2;					//!
	bool m_isSegTag_2; 						//!
	bool m_isCaloTag_2;						//!
	bool m_isStandAl_2;						//!
	int m_Track_nPixelHits;       			//!
    int m_Track_nPixelDeadSensors; 			//!
    int m_Track_nSCTHits;          			//!
    int m_Track_nSCTDeadSensors;   			//!
    int m_Track_nTRTHits;          			//!
    int m_Track_nTRTOutliers;      			//!
    int m_Track_nPixelHoles;       			//!
    int m_Track_nSCTHoles;         			//!
    int m_Track_nIBLHits;					//!
    int m_Track_nBLHits;					//!
    int m_Track_expIBLHits;					//!
    int m_Track_expBLHits;					//!

    float m_Track_pt; 						//!
    float m_Track_eta; 						//!
    float m_Track_phi; 						//!
    float m_Track_charge;					//!
    float m_Track_d0;						//!
    float m_Track_z0 ;						//!
    float m_Track_theta;					//!
    float m_Track_ez0;						//!
    float m_Track_ed0;						//!
    float m_Track_z0_calc; 					//!


	bool m_isCombined_3;					//!
	bool m_isSegTag_3; 						//!
	bool m_isCaloTag_3;						//!
	bool m_isCaloTag_4;						//!
	bool m_isStandAl_3;						//!

    bool passed;                            //!



    int nPixelHits_3;                       //!
    int nPixelDeadSensors_3;                //!
    int nSCTHits_3;                         //!
    int nSCTDeadSensors_3;                  //!
    int nPixelHoles_3;                      //!
    int nSCTHoles_3;                        //!
    int nTRTHits_3;                         //!
    int nTRTOutliers_3;                     //!
    int nIBLHits_3;                         //!
    int nBLHits_3;                          //!
    int expIBLHits_3;                       //!
    int expBLHits_3;                        //!

    bool   m_isMuon_medium;					//!
    bool   m_isMuon_loose;					//!
    bool   m_isMuon_tight;					//!

    double m_MatchdR_mu15_MU6;				//!
    double m_MatchdR_mu15;					//!
    double MatchdR_mu15_3;               //!


    bool m_passTightID;					//!
	bool m_passMediumID;					//!


    int   m_isMuon;							//!
    int   m_istrack;						//!
    int   m_istrack_old;					//!

    float pairdR;       					//!
    float pairpt;							//!
    float pairy;							//!
    float pairm;							//!
    float pairphi;							//!
    float paireta;							//!

    float pairdR_3;							//!
    float pairpt_3;							//!
    float pairy_3;							//!
    float paireta_3;					    //!
    float pairphi_3;					    //!
    float pairm_3;							//!

    float m_ME_Track_pt;					//!
    float m_ME_Track_eta;					//!
    float m_ME_Track_phi;					//!
    float m_ME_Track_charge;				//!
    float m_ME_Track_d0;					//!
    float m_ME_Track_z0;					//!
    float m_ME_Track_theta;					//!
    float m_ME_Track_ed0;					//!

    float d0_trk;								//!
	float ed0_trk;								//!
	float z0_trk;								//!
	float theta_trk;							//!
	float ez0_trk;							//!
	float d0_old;							//!
	float ed0_old;							//!
	float z0_old;							//!
	float theta_old;						//!
	float pt_trk;                               //!
	float phi_trk;							//!
	float eta_trk;							//!
    float m_dR;								//!
    float m_dR_old;							//!




    float m_muonPt_3;						//!
	float m_muonEta_3;						//!
	float m_muonPhi_3;						//!
	float m_muonCharge_3;					//!
	float m_muonQuality_3;					//!
	float m_muonType_3;						//!

    float m_muon2_pt;         				//!
    float m_muon2_eta;						//!
    float m_muon2_phi;						//!
    float m_muon2_charge;					//!



    int m_Muon_n;							//!

    float m_muonELoss;						//!

    bool m_passTightID_3;						//!
	bool m_passMediumID_3;					//!
	bool m_passLooseID_3;					//!

	bool m_passTightID_4;						//!
	bool m_passMediumID_4;					//!
	bool m_passLooseID_4;					//!
	bool m_isCombined_4;					//!




    float m_muonPtcorr; //!
    float m_muonPtcorr_3;//!
    float m_muon2_ptcorr;//!


    float m_Muonid_d0;  					//!
    float m_Muonid_d0sig;					//!
    float m_Muonid_theta;					//!
    float m_Muonid_z0;						//!
    float m_Muonid_phi;						//!
    float m_Muonid_pt;						//!
    float m_Muonid_charge;					//!
    float m_Muonid_eta;						//!


    float m_Muonid_d0_3;  					//!
    float m_Muonid_d0sig_3;					//!
    float m_Muonid_theta_3;					//!
    float m_Muonid_z0_3;						//!
    float m_Muonid_phi_3;						//!
    float m_Muonid_pt_3;						//!
    float m_Muonid_charge_3;					//!
    float m_Muonid_eta_3;						//!

    float m_Muonid_z0sig;						//!
    float m_Muonid_ez0sin;						//!
    float m_Track_ez0sin;						//!

    float m_Muonid_z0sig_3;						//!
    float m_Muonid_ez0sin_3;						//!
    double m_sumEt;								//!
    double percentile;							//!
    int trackn;									//!
    double pile_up_vertices;					//!
    int m_PV;									//!

    bool gap_decision;							//!
    std::vector<bool> isgap;					//!


    //for the calo tagged muons
    bool  m_isCaloTag_c;						//!
    bool  m_isCaloTag_c1;						//!
    bool  m_passed_calo_c;						//!
    bool  m_passed_IDcuts_c;					//!
    float m_caloPt;								//!
    float m_caloPtcorr;							//!
    float m_caloPhi;							//!
    float m_caloCharge;							//!
    float m_caloEta;							//!
    float m_Muonidcalo_d0;     					//!
    float m_Muonidcalo_d0sig;					//!
    float m_Muonidcalo_z0;						//!
    float m_Muonidcalo_z0sig;					//!
    float m_Muonidcalo_theta;					//!
    float m_Muonidcalo_phi;						//!
    float m_Muonidcalo_pt;						//!
    float m_Muonidcalo_charge;					//!
    float m_Muonidcalo_eta;						//!
    float m_Muonidcalo_ez0sin;					//!
    bool m_passMediumID_c;						//!
    bool m_passTightID_c;						//!
    float deltaR_med;							//!
    float deltaR_tight;							//!
    bool is_calo;							    //!
    bool is_ID;									//!
    float dR;									//!

    int nPixelHits_calo;						//!
    int nPixelDeadSensors_calo;					//!
    int nSCTHits_calo;							//!
    int nSCTDeadSensors_calo;					//!
    int nPixelHoles_calo;						//!
    int nSCTHoles_calo;							//!
    int nTRTHits_calo;							//!
    int nTRTOutliers_calo;						//!
    int nIBLHits_calo;							//!
    int nBLHits_calo;							//!
    int expIBLHits_calo;						//!
    int expBLHits_calo;							//!

    bool m_passIDcuts;							//!
    bool m_passIDcuts_3;						//!


    int  m_Muonidcalo_nPixelHits;				//!
	int m_Muonidcalo_nPixelDeadSensors ;		//!
	int m_Muonidcalo_nSCTHits;					//!
	int m_Muonidcalo_nSCTDeadSensors;			//!
	int m_Muonidcalo_nPixelHoles;				//!
	int m_Muonidcalo_nSCTHoles;					//!
	int m_Muonidcalo_nTRTHits;					//!
	int m_Muonidcalo_nTRTOutliers;				//!
	int m_Muonidcalo_nIBLHits;					//!
	int m_Muonidcalo_nBLHits;					//!
	int m_Muonidcalo_expIBLHits;				//!
	int m_Muonidcalo_expBLHits;					//!

	bool pass_calo;								//!
	bool m_calotrack_passIDcuts;				//!
	bool m_track_passIDcuts;					//!
	bool m_Muonidtrack_passIDcuts;				//!



	bool m_debug;                           //!
    bool m_HLT_e15_lhloose;		//!
    //bool m_HLT_e15_lhloose;//!
    bool m_HLT_g35_loose;			//!

    std::vector<double> m_Muon_MatchdR_mu4; 	//!
    std::vector<double> m_Muon_MatchdR_mu4nomucomb; //!
    std::vector<double> m_Muon_MatchdR_mu4trk; 	//!
    std::vector<double> m_Muon_MatchdR_mu6; 	//!
    std::vector<double> m_Muon_MatchdR_mu6trk; 	//!
    std::vector<double> m_Muon_MatchdR_mu6nomucomb; //!
    std::vector<double> m_Muon_MatchdR_mu8; 	//!
    std::vector<double> m_Muon_MatchdR_mu8trk; 	//!
    std::vector<double> m_Muon_MatchdR_mu10; 	//!
    std::vector<double> m_Muon_MatchdR_mu15; 	//!
    std::vector<double> m_Muon_MatchdR_mu15_MU10; 	//!

    std::vector<double> m_Muon_MatchdR_mu10_MU6; 	//!
    std::vector<double> m_Muon_MatchdR_mu15_MU6; 	//!
    std::vector<int>    m_Muon_isTrack;              //!


    std::vector<double> m_vertex_PV_x;//!
    std::vector<double> m_vertex_PV_y;//!
    std::vector<double> m_vertex_PV_z;//!

    std::vector<double> m_vertex_PU_x;//!
    std::vector<double> m_vertex_PU_y;//!
    std::vector<double> m_vertex_PU_z;//!


    std::vector<double> m_Muon_MatchdR_L1MU4;  //!
    std::vector<double> m_Muon_MatchdR_L1MU6;  //!

    std::vector<double> m_Muon_MatchdR_L1MU10;  //!
    std::vector<double> m_Muon_MatchdR_L1MU15;  //!

   /*
    std::vector<double> m_Muon_MatchdR_L2CB; 	//!
    std::vector<double> m_Muon_Matched_RoiEta; 	//!
    std::vector<double> m_Muon_Matched_RoiPhi; 	//!
    std::vector<int>    m_Muon_Matched_RoiWord; //!
    std::vector<int>    m_Muon_Matched_RoiThr; 	//!
    std::vector<int>    m_Muon_Matched_RoiNum; 	//!
    std::vector<int>    m_Muon_Matched_RoiMoreMu; //!
    */

    std::vector<double> m_Muon_pt; 		//!
    std::vector<float>  m_Muon_ptcorr; //!
    std::vector<double> m_Muon_eta; 		//!
    std::vector<double> m_Muon_phi; 		//!
    std::vector<double> m_Muon_ELoss;    	//!
    std::vector<int>    m_Muon_charge; 		//!
    std::vector<int>    m_Muon_quality;         //!
    std::vector<bool>   m_Muon_isTight;         //!
    std::vector<bool>   m_Muon_isMedium;        //!

    std::vector<bool>   M_isCombined;        //!
    std::vector<bool>   M_isSegTag;        //!
    std::vector<bool>   M_isCaloTag;        //!
    std::vector<bool>   M_isStandAl;        //!

    std::vector<int>    m_Muon_type;            //!
    std::vector<int>    m_Muon_passedIDCutsMed;    //!
    std::vector<int>    m_Muon_passedIDCutsTight;    //!


    std::vector<double> m_Muon_reducedChi2;     //!
    std::vector<double> m_Muon_rho;    	        //!
    std::vector<double> m_Muon_me_pt;           //!
    std::vector<double> m_Muon_me_phi;          //!
    std::vector<double> m_Muon_me_theta;        //!
    std::vector<double> m_Muon_me_charge;       //!
    std::vector<double> m_Muon_me_eta;          //!
    std::vector<double> m_Muon_id_pt;           //!
    std::vector<double> m_Muon_id_phi;          //!
    std::vector<double> m_Muon_id_theta;        //!
    std::vector<double> m_Muon_id_charge;       //!
    std::vector<double> m_Muon_id_eta;          //!
    std::vector<double> m_Muon_qOverPsigma;     //!
    std::vector<double> m_Muon_qOverPsignif;    //!
    //std::vector<int>    m_Muon_is_truth;        //!

    std::vector<double> m_me_pt;//!
    std::vector<double> m_me_eta;//!
    std::vector<double> m_me_phi ;//!
    std::vector<double> m_me_charge;//!
    std::vector<double> m_me_d0;//!
    std::vector<double> m_me_z0;//!
    std::vector<double> m_me_theta;//!
    std::vector<double> m_me_ed0;//!


    std::vector<int>    m_Muon_MS_nPrecisionLayers;	//!
    std::vector<int>    m_Muon_MS_nPrecisionHoleLayers;	//!
    std::vector<int>    m_Muon_MS_nPhiLayers;    	//!
    std::vector<int>    m_Muon_MS_nPhiHoleLayers;    	//!
    std::vector<int>    m_Muon_MS_nTrigEtaLayers;	//!
    std::vector<int>    m_Muon_MS_nTrigEtaHoleLayers;	//!
    std::vector<int>    m_Muon_MS_nGoodPrecisionLayers;	//!

    std::vector<double> m_Muon_Track_d0;       		//!
    std::vector<double> m_Muon_Track_z0;       		//!
    //std::vector<double> m_Muon_Track_chiSquared;   	//!
    //std::vector<double> m_Muon_Track_numberDoF;   	//!
    std::vector<int>    m_Muon_Track_nPixelHits;        //!
    std::vector<int>    m_Muon_Track_nPixelDeadSensors; //!
    std::vector<int>    m_Muon_Track_nSCTHits;          //!
    std::vector<int>    m_Muon_Track_nSCTDeadSensors;   //!
    std::vector<int>    m_Muon_Track_nTRTHits;          //!
    std::vector<int>    m_Muon_Track_nTRTOutliers;      //!
    std::vector<int>    m_Muon_Track_nPixelHoles;       //!
    std::vector<int>    m_Muon_Track_nSCTHoles;         //!


    /*
    std::vector<double> m_Track_pt; 		   //!
    std::vector<double> m_Track_eta;		   //!
    std::vector<double> m_Track_theta;		   //!
    std::vector<double> m_Track_phi;		   //!
    std::vector<int>    m_Track_charge;  	   //!
    std::vector<double> m_Track_d0;          	   //!
    std::vector<double> m_Track_z0;          	   //!
    std::vector<int>    m_Track_nPixelHits;        //!
    std::vector<int>    m_Track_nPixelDeadSensors; //!
    std::vector<int>    m_Track_nSCTHits;          //!
    std::vector<int>    m_Track_nSCTDeadSensors;   //!
    std::vector<int>    m_Track_nTRTHits;          //!
    std::vector<int>    m_Track_nTRTOutliers;      //!
    std::vector<int>    m_Track_nPixelHoles;       //!
    std::vector<int>    m_Track_nSCTHoles;         //!

    */

    std::vector<double> m_TruthMuon_pt; 	//!
    std::vector<double> m_TruthMuon_eta; 	//!
    std::vector<double> m_TruthMuon_phi; 	//!
    std::vector<int>    m_TruthMuon_charge;     //!
    std::vector<int>    m_TruthMuon_pdgId;      //!
    std::vector<int>    m_TruthMuon_barcode;    //!
    std::vector<int>    m_TruthMuon_status;     //!
    std::vector<int>    m_TruthMuon_matched;    //!



    //electrons
    std::vector<bool> m_El_MatchdR_mu15; 	//!
    std::vector<double> m_El_pt; 		//!
    std::vector<double> m_El_eta; 		//!
    std::vector<double> m_El_phi; 		//!
    std::vector<int>    m_El_charge; 		//!
    std::vector<bool>   m_El_isTight;         //!
    std::vector<bool>   m_El_isMedium;        //!
    std::vector<bool>   m_El_isLoose;       //!



    //for new calibrated muons
    //xAOD::TStore* m_store;		//!
    //xAOD::MuonContainer* m_newMuons; //!

    // CP tools
    GoodRunsListSelectionTool 			*m_grl; 			//!
    TrigConf::xAODConfigTool  			*configTool; 		//!
    Trig::TrigDecisionTool    			*trigDecTool; 		//!
    Trig::TrigMuonMatching    			*m_trigMuonMatching; 	//!
    //Trig::MatchingTool	     		*m_matchtool;             //!
    CP::MuonCalibrationAndSmearingTool 	*m_muonCorr;     //!
    CP::MuonSelectionTool         		*m_muonSelection; 	//!
     CP::MuonSelectionTool         		*m_muonSelection1; 	//!
    CP::MuonSelectionTool         		*m_muonSelection2; 	//!
     CP::MuonSelectionTool        		*m_muonSelection3; 	//!
     Trig::MatchingTool 	      		*m_tmt; 			//!
     HI::HIPileupTool                   *m_hiPileup;    //!
     HI::HICentralityTool         		*centTool;//!
     ZDC::ZdcAnalysisTool 				*m_zdcTools; //!


      Trig::TrigEgammaMatchingTool    *m_match_tool;   	//!
    CP::EgammaCalibrationAndSmearingTool* m_egammaCalibrationAndSmearingTool; //!


    AsgElectronLikelihoodTool* m_electronLHLooseSelector; //!
    AsgElectronLikelihoodTool* m_electronLHMediumSelector; //!
    AsgElectronLikelihoodTool* m_electronLHTightSelector; //!

   // AsgPhotonIsEMSelector* m_photon_tight; //!
  //  AsgPhotonIsEMSelector* m_photon_loose; //!


    //std::vector<std::string> trigList; //!
    // this is a standard constructor
    Zreco ();

    // these are the functions inherited from Algorithm
    virtual EL::StatusCode setupJob (EL::Job& job);
    virtual EL::StatusCode fileExecute ();
    virtual EL::StatusCode histInitialize ();
    virtual EL::StatusCode changeInput (bool firstFile);
    virtual EL::StatusCode initialize ();
    virtual EL::StatusCode execute ();
    virtual EL::StatusCode postExecute ();
    virtual EL::StatusCode finalize ();
    virtual EL::StatusCode histFinalize ();
    void clearVector();
    bool selectCluster(float,float,int);
    int is_this_good_bcid(int run, int bcid);
  // this is needed to distribute the algorithm to the workers
  ClassDef(Zreco, 1);
};

#endif
