#ifdef _MSC_VER
	#define _CRTDBG_MAP_ALLOC
	#define _CRTDBG_MAP_ALLOC_NEW
	#include <stdlib.h>	
	#include <crtdbg.h>
	//#define VLD_AGGREGATE_DUPLICATES
	//#define VLD_MAX_DATA_DUMP 0
	//#include <vld.h>
#endif

#include <SC_Lib.h>

int main(int argc, char *argv[]) {
	char audioFile[sclib::bufferSize], sceneFile[sclib::bufferSize], segmentFile[sclib::bufferSize], iniFile[sclib::bufferSize];
	double videoFrameRate;
	bool res = false;
	int method;
	SC_TweakableParameters* pTweak;
	SC_MainTasks mainTasks;
	
	setbuf(stdout, NULL); //get instant printf()s
	setbuf(stderr, NULL); //get instant printf()s
	std::ios::sync_with_stdio(); //syncronize couts and printf()s
#ifdef _MSC_VER
	_CrtSetDbgFlag ( _CRTDBG_ALLOC_MEM_DF | _CRTDBG_LEAK_CHECK_DF ); //get memory leaks reported in vc
#endif
	printf("One moment, I'm initializing..."); //reynolds-clustering-small
	sprintf(audioFile, "%s", ((argc > 1) ? argv[1] : "C:/Daten/TIMIT/timit-all.crp"));
	sprintf(segmentFile, "%s", ((argc > 2) ? argv[2] : "-"));
	sprintf(sceneFile, "%s", ((argc > 3) ? argv[3] : "-")); 
	sprintf(iniFile, "%s", ((argc > 4) ? argv[4] : "C:/Users/stdm/Documents/Projekte/audioseg/data/time_ZHAW.ini"));
	videoFrameRate = (argc > 5) ? atof(argv[5]) : 25.0; //<--- A T T E N T I O N  ! ! !
	method = (argc > 6) ? atoi(argv[6]) : -1;
	pTweak = new SC_TweakableParameters(iniFile);
	for (int i = 7; i < argc; i++) { //extract and set commandline arguments
		if (argv[i][0] == '-') {
			char *key=NULL, *value=NULL;
			if (sclib::extractKeyValue(&argv[i][1], key, value) == true) {
				if (pTweak->setByName(key, value) == true) {
					printf("\n  Argument '%s=%s' set", key, value);
				} else {
					printf("\n  Error in argument '%s': Key-value-pair '%s=%s' is not valid", argv[i], key, value);
				}
				MFree_1D(key);
				MFree_1D(value);
			}
		}
	}
	printf("\nUsed parameters:\n  audioFile=%s\n  segmentFile=%s\n  sceneFile=%s\n  iniFile=%s\n  videoFrameRate=%f\n  method=%d\n", audioFile, segmentFile, sceneFile, iniFile, videoFrameRate, method);
	printf("done!\n");

	if (!sclib::pathExists(pTweak->debug.debugDir)) {printf("\n\nPath for storing results '%s' not accesible! Giving up...\n", pTweak->debug.debugDir); exit(-1);}
	if (!sclib::fileExists(audioFile) && !strncmp(audioFile, "-", sclib::bufferSize) && !strncmp(audioFile, "", sclib::bufferSize)) {printf("\n\nAudiofile '%s' not found! Giving up...\n", audioFile); exit(-1);}
	if (!sclib::fileExists(segmentFile) && !strncmp(segmentFile, "-", sclib::bufferSize) && !strncmp(segmentFile, "", sclib::bufferSize)) {printf("\n\nGround truth for audiosegements '%s' not found! Giving up...\n", segmentFile); exit(-1);}
	if (!sclib::fileExists(sceneFile) && !strncmp(segmentFile, "-", sclib::bufferSize) && !strncmp(segmentFile, "", sclib::bufferSize)) {printf("\n\nGround truth for scenes '%s' not found! Creating pseudo scene boundarys instead!", sceneFile);}
 
	mainTasks.initParameters(pTweak, audioFile, segmentFile, sceneFile, videoFrameRate);

	switch (method) {
		case -1: //the case for changes in-code

			//build time model for 40 speakers (to later cluster)
			//mainTasks.trainTimeWorldModelAndNormMatrix("C:/Users/stdm/Documents/Projekte/audioseg/data/test/timit-train.mfccPitch");
			mainTasks.createTIMITcontextVectors("TIMIT_mfccPitch_tracjectories_all_speakers.csv");

			//trainTimeModels(): -> needs adaptation/renaming
			//1. build norm file for mfcc&pitch on complete TIMIT train
			//(actually not: 2. build time world model in complete TIMIT train (no 2nd time length normalization required!))
			//createContextVectorMatrix(): -> needs adaptation (e.g. transposition) / c&p
			//(3. extract all trajectories of TIMIT test clustering 40)
			//scivo(): -> adaptation needed?
			//4. build time models for each utterance in TIMIT test clustering 40
			//5. perform clustering using these models

			//mainTasks.selfSimilarityAnalysis(3000, 500, 0, sclib::metaFeature_segDis_BIC);

			//mainTasks.wav2timeModel("d:/data/code/sprecherklassifikation/data/test");

			/*
			pTweak->modelHandler.maxSpeakerModelOrder = 32;
			pTweak->modelHandler.orderGuessMode = sclib::guessHeuristic;
			pTweak->debug.debugMode |= sclib::dbSpeakerModels|sclib::dbWav; //|sclib::dbSpeechTrainData;
			pTweak->modelHandler.speakerModelFeature = sclib::featurePitch;
			pTweak->featureMfcc.method = sclib::modeSClib;
			pTweak->featureLpc.LPCorder = 12;
			pTweak->general.firstScene = 3;
			pTweak->general.lastScene = 3;
			pTweak->segmentationHandler.audioTypeMode = sclib::algorithm_nothing;
			pTweak->segmentationHandler.changeDetectorMode = sclib::algorithm_nothing;
			pTweak->segmentationHandler.silenceDetectorMode = sclib::algorithm_nothing;
			pTweak->segmentationHandler.vUvDetectorMode = sclib::algorithm_nothing;
			res = mainTasks.scivo();
			//*/
			
			/*
			pTweak->modelHandler.foregroundModelType = sclib::mtGMM_new;
			pTweak->modelHandler.orderGuessMode = sclib::guessHeuristic;
			pTweak->modelHandler.msPerGaussian = 10;
			pTweak->modelHandler.maxSpeakerModelOrder = 32;
			pTweak->speakerClusterer.distanceMeasure = sclib::dmGLR;
			res = mainTasks.buildReynoldsSpeakerModels(0.0);
			//res = mainTasks.testReynoldsSpeakerModels("../data/test/reynoldsModelList.txt", 0.0);
			//*/

			/*
			pTweak->debug.useDebugPrefix = false;
			pTweak->featureMfcc.method = sclib::modeSClib;
			res = mainTasks.kottiBICtest("../../../Rohdaten/conTIMIT/train_bicpenaltyfactor.crp", "../../../Rohdaten/conTIMIT/train_utterancedurationmodeling.crp", "../../../Rohdaten/conTIMIT/train_fs_change.crp", "../../../Rohdaten/conTIMIT/train_fs_nochange.crp");
			//*/

			/*	
			//pTweak->modelTime.trajectoryStep = 5;
			pTweak->modelTime.removeTiming = false;
			pTweak->modelHandler.speakerModelFeature = sclib::featureSDP;
			pTweak->featureSdp.pictureSize = 23;
			pTweak->featureSdp.tau = 11;
			res = mainTasks.createContextVectorMatrix("../data/test/matrix.dat", "../data/test/matrix.gt");
			//res = mainTasks.testContextVectorMatrix("../data/test/matrix.dat", "../data/test/matrix.gt");
			//*/

			break;
		case 0: //starting from here are the cases for command-line parameter settings
			res = mainTasks.trainAudioTypeLZL();
			break;
		case 1: 
			res = mainTasks.enhancer();
			break;
		case 2: 
			res = mainTasks.modeller(argv[7], argv[8], argv[9], atoi(argv[10]), atoi(argv[11]), atoi(argv[12]), atoi(argv[13]), atoi(argv[14]));
			break;
		case 3: 
			res = mainTasks.combineModels(argv[7], argv[8], argv[9], atoi(argv[10]));
			break;
		case 4: 
			res = mainTasks.statistician();
			break;
		//case 5: 
		//	res = mainTasks.testSDP();
		//	break;
		case 6: 
			res = mainTasks.vuvDetection();
			break;
		case 7: 
			res = mainTasks.wesley(atoi(argv[7]), atoi(argv[8]), argv[9]);
			break;
		case 8: 
			res = mainTasks.scivo();
			break;
		case 9: 
			res = mainTasks.pitchTest();
			break;
		//case 10: 
		//	res = mainTasks.buildAllANNresidualSpeakerModels();
		//	break;
		//case 11: 
		//	res = mainTasks.testAllANNresidualSpeakerModels(argv[7]);
		//	break;
		case 12:
			res = mainTasks.trainSilenceLZL();
			break;
		case 13:
			res = mainTasks.testAudioType();
			break;
		case 14:
			res = mainTasks.synthesisTest();
			break;
		case 15:
			res = mainTasks.javaTest();
			break;
		case 16:
			res = mainTasks.kottiBICtest(argv[7], argv[8], argv[9], argv[10]);
			break;
		case 17:
			res = mainTasks.buildReynoldsSpeakerModels(atof(argv[7]));
			break;
		case 18:
			res = mainTasks.testReynoldsSpeakerModels(argv[7], atof(argv[8]));
			break;
		case 19: 
			res = mainTasks.HHTtest(atof(argv[7]));
			break;
		case 20: 
			res = mainTasks.synthesizeCorpus();
			break;
		case 21: 
			res = mainTasks.wav2gmm(argv[7]);
			break;
		case 22: 
			res = mainTasks.timeModelTest();
			break;
		case 23:
			res = mainTasks.splicer(atoi(argv[7]), atoi(argv[8]), atof(argv[9]));
			break;
		case 24:
			res = mainTasks.featureDisplacemetTest(atoi(argv[7]), sclib::atob(argv[8]));
			break;
		case 25:
			res = mainTasks.createTemplates(atoi(argv[7]));
			break;
		case 26:
			res = mainTasks.createContextVectorMatrix(argv[7], argv[8]);
			break;
		case 27:
			res = mainTasks.testContextVectorMatrix(argv[7], argv[8]);
			break;
		//case 28:
		//	res = mainTasks.selfSimilarityAnalysis(atoi(argv[7]), atoi(argv[8]), atoi(argv[9]), atoi(argv[10]));
		//	break;
		case 29:
			res = mainTasks.gmmVisualizationTest();
			break;
		case 30: 
			res = mainTasks.videoFeatureExtractor(argv[7], atoi(argv[8]));
			break;
		default:
			printf("\nNo valid task choosen!\n");
			res = false;
			break;
	}

	MFree_0D(pTweak);

	printf("\nThe task returned: %d", res);
	if (method < 0) {
		printf("\nPress the anykey...");
		getchar();
	}

#ifdef _MSC_VER
	_CrtDumpMemoryLeaks();
#endif

	return (res == true) ? 0 : -1;
}
