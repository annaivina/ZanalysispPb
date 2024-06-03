#include "xAODRootAccess/Init.h"
#include "SampleHandler/SampleHandler.h"
#include "SampleHandler/ScanDir.h"
#include "SampleHandler/ToolsDiscovery.h"
#include "EventLoop/Job.h"
#include "EventLoop/OutputStream.h"
#include "EventLoopAlgs/NTupleSvc.h"
#include "EventLoop/DirectDriver.h"
#include "SampleHandler/DiskListLocal.h"
#include <iostream>
#include <TSystem.h>
#include "EventLoopGrid/PrunDriver.h"
#include "EventLoopGrid/GridDriver.h"

#include "ZanalysispPb/Zreco.h"

int main( int argc, char* argv[] )

{
	char submitDir[100]="submitDir";
	int vers=1;
	int runnumber = 313063;
	int mnumber   = 1736;
	int fnumber   = 774;

	if( argc > 1 )
		{
			vers=atoi(argv[1]);
		}

    if( argc > 2)
	   {
	        runnumber=atoi(argv[2]);
	        if(runnumber>313575) {fnumber = 781; mnumber = 1741;}

	   }

	sprintf(submitDir,"submitDirVal_v%d%d", vers,runnumber);
	//sprintf(submitDir,"submitDirVal_v%d", vers);

	xAOD::Init().ignore();
	SH::SampleHandler sh;

	char input[2000];
	char outputfile[2000];
	char outdataset[2000];
	switch(vers)
	{
        case 2:    sprintf(input,"/srv01/cgrp/users/annai/annai/2018_ZmumuAnalysis/pPb2016/data16_hip8TeV.00%d.physics_Main.merge.DAOD_HION5.f%d_m%d_p2961",runnumber,fnumber,mnumber);
		           sprintf(outputfile,"Reco_efficiency");
				   break;

		case 5:  sprintf(input,"/srv01/cgrp/users/annai/annai/2018_ZmumuAnalysis/pPb2016/data16_hip8TeV.00313063.physics_Main.merge.DAOD_HION5.f774_m1736_p2961");
		         sprintf(outputfile,"try");break;

		case 10:  sprintf(input,"/srv01/cgrp/users/annai/annai/2018_ZmumuAnalysis/min_bias_pPb_hion2/data16_hip8TeV.00313285.physics_Main.recon.AOD.f774_m1736");
		         sprintf(outputfile,"fcal_cent_run1");break;

		case 7:  sprintf(input,"/mnt/lustre/cgrp/atlas_hi/hion6/data16_hip8TeV.00313136.physics_Main.merge.DAOD_HION6.f774_m1736_p3052");
		         sprintf(outputfile,"test_hion6");break;


		case 20:  sprintf(input,"/srv01/cgrp/users/annai/annai/2018_ZmumuAnalysis/min_bias_pPb_hion2/test");
		          sprintf(outputfile,"test");break;







		//case 1:  sprintf(input,"/afs/cern.ch/user/a/aivina/aivina/Zmumuanalysis_pPb2016/data");
		       //  sprintf(outputfile,"Derivation_test");break;
		//case 3:   sprintf(input,"/afs/cern.ch/user/a/aivina/aivina/Zmumuanalysis_pPb2016/data/data16_hip8TeV.00314014.physics_Main.merge.DAOD_HION5.f781_m1741_p2961");
		       //   sprintf(outputfile,"Derivation_test_v3");break;

		//case 4:   sprintf(input,"/afs/cern.ch/work/a/aivina/Zmumuanalysis_pPb2016/data/data16_hip8TeV.00313629.physics_Main.merge.DAOD_HION5.f781_m1741_p2961");
		       //   sprintf(outputfile,"Derivation_test_v4");break;

        ///code testing
       // case 6:   sprintf(input,"/mnt/Lustre/cgrp/atlas_hi/annai/2018_ZmumuAnalysis/pPb2016/data16_hip8TeV.00313063.physics_Main.merge.DAOD_HION5.f774_m1736_p2961");
		       //   sprintf(outputfile,"test-data");break;

		//case 7:   sprintf(input,"/mnt/Lustre/cgrp/atlas_hi/annai/2018_ZmumuAnalysis/pPb2016/data16_hip8TeV.00313067.physics_Main.merge.DAOD_HION5.f774_m1736_p2961");
		        //  sprintf(outputfile,"test-data");break;

		//case 8:   sprintf(input,"/mnt/Lustre/cgrp/atlas_hi/annai/2018_ZmumuAnalysis/pPb2016/data16_hip8TeV.00313100.physics_Main.merge.DAOD_HION5.f774_m1736_p2961");
		        //  sprintf(outputfile,"test-data");break;

		//case 9:   sprintf(input,"/mnt/Lustre/cgrp/atlas_hi/annai/2018_ZmumuAnalysis/pPb2016/data16_hip8TeV.00313107.physics_Main.merge.DAOD_HION5.f774_m1736_p2961");
		       //   sprintf(outputfile,"test-data");break;

		//case 10:   sprintf(input,"/mnt/Lustre/cgrp/atlas_hi/annai/2018_ZmumuAnalysis/pPb2016/data16_hip8TeV.00313136.physics_Main.merge.DAOD_HION5.f774_m1736_p2961");
		        //  sprintf(outputfile,"test-data");break;

		//case 11:   sprintf(input,"/srv01/cgrp/users/annai/annai/2018_ZmumuAnalysis/data/mc15_pPb8TeV.361107.PowhegPythia8EvtGen_AZNLOCTEQ6L1_Zmumu.merge.AOD.e5366_s3164_r9441_r9006");
		        //  sprintf(outputfile,"test-mc");break;


         //Period B
        //version 6 is from the 3d of november //v10 -with mu4_trkloose trigger
		case 401: sprintf(input,"data16_hip8TeV.00313063.physics_Main.recon.AOD.f774_m1736");
		          sprintf(outputfile,"run00313063_reco_v2_01");
				  sprintf(outdataset,"user.aivina.00313063PeriodB_reco_v2_01");
				  break;


		case 402: sprintf(input,"data16_hip8TeV.00313067.physics_Main.recon.AOD.f774_m1736");
		          sprintf(outputfile,"run00313067_reco_v2");
				  sprintf(outdataset,"user.aivina.00313067PeriodB_reco_v2");
				  break;

		case 403: sprintf(input,"data16_hip8TeV.00313100.physics_Main.recon.AOD.f774_m1736");
		          sprintf(outputfile,"run00313100_reco_v2");
				  sprintf(outdataset,"user.aivina.00313100PeriodB_reco_v2");
				  break;

	    case 404: sprintf(input,"data16_hip8TeV.00313107.physics_Main.recon.AOD.f774_m1736");
		          sprintf(outputfile,"run00313107_reco_v2");
				  sprintf(outdataset,"user.aivina.00313107PeriodB_reco_v2");
				  break;

		case 405: sprintf(input,"data16_hip8TeV.00313136.physics_Main.recon.AOD.f774_m1736");
		          sprintf(outputfile,"run00313136_reco_v2");
				  sprintf(outdataset,"user.aivina.00313136PeriodB_reco_v2");
				  break;

		case 406: sprintf(input,"data16_hip8TeV.00313295.physics_Main.recon.AOD.f774_m1736");
		          sprintf(outputfile,"run00313295_reco_v2");
				  sprintf(outdataset,"user.aivina.00313295PeriodB_reco_v2");
				  break;

		case 407: sprintf(input,"data16_hip8TeV.00313187.physics_Main.recon.AOD.f774_m1736");
                  sprintf(outputfile,"run00313187_reco_v2");
				  sprintf(outdataset,"user.aivina.00313187PeriodB_reco_v2");
				  break;


        case 408: sprintf(input,"data16_hip8TeV.00313259.physics_Main.recon.AOD.f774_m1736");
                  sprintf(outputfile,"run00313259_reco_v2");
				  sprintf(outdataset,"user.aivina.00313259PeriodB_reco_v2");
				  break;

        case 409: sprintf(input,"data16_hip8TeV.00313285.physics_Main.recon.AOD.f774_m1736");
                  sprintf(outputfile,"run00313285_reco_v2");
				  sprintf(outdataset,"user.aivina.00313285PeriodB_reco_v2");
				  break;

        case 410: sprintf(input,"data16_hip8TeV.00313333.physics_Main.recon.AOD.f774_m1736");
                  sprintf(outputfile,"run00313333_reco_v2");
				  sprintf(outdataset,"user.aivina.00313333PeriodB_reco_v2");
				  break;

        case 411: sprintf(input,"data16_hip8TeV.00313435.physics_Main.recon.AOD.f774_m1736");
                  sprintf(outputfile,"run00313435_reco_v2");
				  sprintf(outdataset,"user.aivina.00313435PeriodB_reco_v2");
				  break;

        //Period C
	    case 601: sprintf(input,"data16_hip8TeV.00313572.physics_Main.recon.AOD.f774_m1736");
	              sprintf(outputfile,"run00313572_reco_v2");
	              sprintf(outdataset,"user.aivina.00313572PeriodC_reco_v2");
                  break;

        case 602: sprintf(input,"data16_hip8TeV.00313574.physics_Main.recon.AOD.f774_m1736");
	              sprintf(outputfile,"run00313574_reco_v2");
	              sprintf(outdataset,"user.aivina.00313574PeriodC_reco_v2");
                  break;

        case 603: sprintf(input,"data16_hip8TeV.00313575.physics_Main.recon.AOD.f774_m1736");
	              sprintf(outputfile,"run00313575_reco_v2");
	              sprintf(outdataset,"user.aivina.00313575PeriodC_reco_v2");
                  break;

        case 604: sprintf(input,"data16_hip8TeV.00313603.physics_Main.recon.AOD.f781_m1741");
	              sprintf(outputfile,"run00313603_reco_v2");
	              sprintf(outdataset,"user.aivina.00313603PeriodC_reco_v2");
                  break;

        case 605: sprintf(input,"data16_hip8TeV.00313629.physics_Main.recon.AOD.f781_m1741");
	              sprintf(outputfile,"run00313629_reco_v2");
	              sprintf(outdataset,"user.aivina.00313629PeriodC_reco_v2");
                  break;

        case 606: sprintf(input,"data16_hip8TeV.00313630.physics_Main.recon.AOD.f781_m1741");
	              sprintf(outputfile,"run00313630_reco_v2");
	              sprintf(outdataset,"user.aivina.00313630PeriodC_reco_v2");
                  break;

        case 607: sprintf(input,"data16_hip8TeV.00313688.physics_Main.recon.AOD.f781_m1741");
	              sprintf(outputfile,"run00313688_reco_v2");
	              sprintf(outdataset,"user.aivina.00313688PeriodC_reco_v2");
                  break;

        case 608: sprintf(input,"data16_hip8TeV.00313695.physics_Main.recon.AOD.f781_m1741");
	              sprintf(outputfile,"run00313695_reco_v2");
	              sprintf(outdataset,"user.aivina.00313695PeriodC_reco_v2");
                  break;

        case 609: sprintf(input,"data16_hip8TeV.00313833.physics_Main.recon.AOD.f781_m1741");
	              sprintf(outputfile,"run00313833_reco_v2");
	              sprintf(outdataset,"user.aivina.00313833PeriodC_reco_v2");
                  break;

        case 610: sprintf(input,"data16_hip8TeV.00313878.physics_Main.recon.AOD.f781_m1741");
	              sprintf(outputfile,"run00313878_reco_v2");
	              sprintf(outdataset,"user.aivina.00313878PeriodC_reco_v2");
                  break;

        case 611: sprintf(input,"data16_hip8TeV.00313929.physics_Main.recon.AOD.f781_m1741");
	              sprintf(outputfile,"run00313929_reco_v2");
	              sprintf(outdataset,"user.aivina.00313929PeriodC_reco_v2");
                  break;

        case 612: sprintf(input,"data16_hip8TeV.00313935.physics_Main.recon.AOD.f781_m1741");
	              sprintf(outputfile,"run00313935_reco_v2");
	              sprintf(outdataset,"user.aivina.00313935PeriodC_reco_v2");
                  break;

        case 613: sprintf(input,"data16_hip8TeV.00313984.physics_Main.recon.AOD.f781_m1741");
	              sprintf(outputfile,"run00313984_reco_v2");
	              sprintf(outdataset,"user.aivina.00313984PeriodC_reco_v2");
                  break;

        case 614: sprintf(input,"data16_hip8TeV.00314014.physics_Main.recon.AOD.f781_m1741");
	              sprintf(outputfile,"run00314014_reco_v2");
	              sprintf(outdataset,"user.aivina.00314014PeriodC_reco_v2");
                  break;

        case 615: sprintf(input,"data16_hip8TeV.00314077.physics_Main.recon.AOD.f781_m1741");
	              sprintf(outputfile,"run00314077_reco_v2");
	              sprintf(outdataset,"user.aivina.00314077PeriodC_reco_v2");
                  break;

        case 616: sprintf(input,"data16_hip8TeV.00314105.physics_Main.recon.AOD.f781_m1741");
	              sprintf(outputfile,"run00314105_reco_v2");
	              sprintf(outdataset,"user.aivina.00314105PeriodC_reco_v2");
                  break;

        case 617: sprintf(input,"data16_hip8TeV.00314112.physics_Main.recon.AOD.f781_m1741");
	              sprintf(outputfile,"run00314112_reco_v2");
	              sprintf(outdataset,"user.aivina.00314112PeriodC_reco_v2");
                  break;

        case 618: sprintf(input,"data16_hip8TeV.00314157.physics_Main.recon.AOD.f781_m1741");
	              sprintf(outputfile,"run00314157_reco_v2");
	              sprintf(outdataset,"user.aivina.00314157PeriodC_reco_v2");
                  break;

        case 619: sprintf(input,"data16_hip8TeV.00314170.physics_Main.recon.AOD.f781_m1741");
	              sprintf(outputfile,"run00314170_reco_v2");
	              sprintf(outdataset,"user.aivina.00314170PeriodC_reco_v2");
                  break;

        default: std::cout<<"Wrong version"<<std::endl;
		}

	const char* inputFilePath = gSystem->ExpandPathName (input);
	if(vers<100)
	{
	  SH::ScanDir().filePattern("*AOD*").scan(sh,inputFilePath);
	  sh.setMetaString( "nc_tree", "CollectionTree" );
	  sh.print();
    }
	else
	{

 	SH::scanDQ2 (sh, inputFilePath);
	sh.setMetaString( "nc_grid_filter", "*AOD*");
	}

  // Create an EventLoop job:
  EL::Job job;
  job.sampleHandler( sh );
  job.options()->setDouble (EL::Job::optRemoveSubmitDir, 1);
  job.options()->setString(EL::Job::optXaodAccessMode, EL::Job::optXaodAccessMode_athena);

  // define an output and an ntuple associated to that output
  EL::OutputStream output  (outputfile);
  job.outputAdd (output);
  EL::NTupleSvc *ntuple = new EL::NTupleSvc (outputfile);
  job.algsAdd (ntuple);

  // Add our analysis to the job:
  Zreco *alg = new Zreco();
  job.algsAdd( alg );



  alg->outputName = outputfile; // give the name of the output to our algorithm
  // Run the job using the local/direct driver:
  if(vers<100)
  {
    EL::DirectDriver driver;
	driver.submit( job, submitDir );
	}

  else
  {
  EL::PrunDriver driver;
  driver.options()->setString("nc_outputSampleName",outdataset);
  driver.options()->setString("nc_cmtConfig", "x86_64-slc6-gcc49-opt");
  driver.options()->setDouble(EL::Job::optGridMergeOutput, 1);
  driver.options()->setDouble("nc_nFilesPerJob", 10); //so the site will not crush
 driver.submitOnly( job, submitDir );
 }



//The place where I have put the vers means that this one u can submit on the grid (read info)

  return 0;
}
