/**************************************************************************/
/*    All tweakable parameters in all the	SC_* lib                        */
/*																																				*/
/*    Author  : Thilo Stadelmann																					*/
/*    Date    : 22.03.2004																								*/
/**************************************************************************/

#include <string.h>
#include <stdio.h>
#include "SC_TweakableParameters.h"
#include "SC_Aux.h"
#include "SC_Ini.h"
#include <SV_Error.h>

//====================================================================================================================
//	Constructors for SC_FeaturePar* classes
//====================================================================================================================
SC_TweakableParameters::SC_FeaturePar::SC_FeaturePar() {
	this->frameSize = 0;
	this->frameStep = 0;
	this->lowCut = 0.0;
	this->highCut = 0.0;
	this->preEmphasizeFactor = 1.0;
	this->sampleRate = 0.0;
	this->Next = NULL;
	this->client = sclib::algorithm_nothing;
	this->featureType = sclib::featureNoFeature;
}

SC_TweakableParameters::SC_FeatureBandPeriodicityPar::SC_FeatureBandPeriodicityPar() : SC_TweakableParameters::SC_FeaturePar() {
	this->featureType = sclib::featureBandPeriodicity;
}

SC_TweakableParameters::SC_FeatureBrightnessBandwidthPar::SC_FeatureBrightnessBandwidthPar() : SC_TweakableParameters::SC_FeaturePar() {
	this->featureType = sclib::featureBrightnessBandwidth;
}

SC_TweakableParameters::SC_FeatureFbEPar::SC_FeatureFbEPar() : SC_TweakableParameters::SC_FeaturePar() {
	this->featureType = sclib::featureFbE;
}

SC_TweakableParameters::SC_FeatureSTEPar::SC_FeatureSTEPar() : SC_TweakableParameters::SC_FeaturePar() {
	this->featureType = sclib::featureSTE;
}

SC_TweakableParameters::SC_FeatureMFCCPar::SC_FeatureMFCCPar() : SC_TweakableParameters::SC_FeaturePar() {
	this->featureType = sclib::featureMFCC;
}

SC_TweakableParameters::SC_FeatureNFRPar::SC_FeatureNFRPar() : SC_TweakableParameters::SC_FeaturePar() {
	this->featureType = sclib::featureNFR;
}

SC_TweakableParameters::SC_FeatureSpectrumPar::SC_FeatureSpectrumPar() : SC_TweakableParameters::SC_FeaturePar() {
	this->featureType = sclib::featureSpectrum;
}

SC_TweakableParameters::SC_FeatureSpectrumFluxPar::SC_FeatureSpectrumFluxPar() : SC_TweakableParameters::SC_FeaturePar() {
	this->featureType = sclib::featureSpectrumFlux;
}

SC_TweakableParameters::SC_FeatureSubBandPowerPar::SC_FeatureSubBandPowerPar() : SC_TweakableParameters::SC_FeaturePar() {
	this->featureType = sclib::featureSubbandPower;
}

SC_TweakableParameters::SC_FeatureZCRPar::SC_FeatureZCRPar() : SC_TweakableParameters::SC_FeaturePar() {
	this->featureType = sclib::featureZCR;
}

SC_TweakableParameters::SC_FeatureLPCPar::SC_FeatureLPCPar() : SC_TweakableParameters::SC_FeaturePar() {
	this->featureType = sclib::featureLPC;
}

SC_TweakableParameters::SC_FeatureLPCresidualPar::SC_FeatureLPCresidualPar() : SC_TweakableParameters::SC_FeaturePar() {
	this->featureType = sclib::featureLPCresidual;
}

SC_TweakableParameters::SC_FeatureLSPPar::SC_FeatureLSPPar() : SC_TweakableParameters::SC_FeaturePar() {
	this->featureType = sclib::featureLSP;
}

SC_TweakableParameters::SC_FeaturePitchPar::SC_FeaturePitchPar() : SC_TweakableParameters::SC_FeaturePar() {
	this->featureType = sclib::featurePitch;
}

SC_TweakableParameters::SC_FeatureFormantPar::SC_FeatureFormantPar() : SC_TweakableParameters::SC_FeaturePar() {
	this->featureType = sclib::featureFormant;
}

SC_TweakableParameters::SC_FeatureSDPPar::SC_FeatureSDPPar() : SC_TweakableParameters::SC_FeaturePar() {
	this->featureType = sclib::featureSDP;
}

SC_TweakableParameters::SC_FeatureSamplesPar::SC_FeatureSamplesPar() : SC_TweakableParameters::SC_FeaturePar() {
	this->featureType = sclib::featureSamples;
}

//====================================================================================================================
//	The constructor
//====================================================================================================================
SC_TweakableParameters::SC_TweakableParameters(const char *iniFileName, bool verbose) {
  char *key = NULL, *value = NULL;
	SC_Ini *pIni = NULL;

	//initialize the mapping for setByName(); lowercase here because it is needed by Mediana (grr) and to be more fault-tolerant
	this->parameters.map["classifieradaboost.maxweakclassifiers"] = 1;
	this->parameters.map["classifieradaboost.weakclassifiertype"] = 2;
	this->parameters.map["classifierdecisionstump.splitmultiway"] = 3;
	this->parameters.map["classifierml.maxorder"] = 4;
	this->parameters.map["classifierml.minorder"] = 5;
	this->parameters.map["classifierml.modeltype"] = 6;
	this->parameters.map["classifiersvm.c"] = 7;
	this->parameters.map["classifiersvm.cache_size"] = 8;
	this->parameters.map["classifiersvm.coef0"] = 9;
	this->parameters.map["classifiersvm.cvfolds"] = 10;
	this->parameters.map["classifiersvm.cvmaxdatasetsize"] = 11;
	this->parameters.map["classifiersvm.cvcoarsecmin"] = 12;
	this->parameters.map["classifiersvm.cvcoarsecmax"] = 13; 
	this->parameters.map["classifiersvm.cvcoarsecstep"] = 14; 
	this->parameters.map["classifiersvm.cvcoarsegammamin"] = 15; 
	this->parameters.map["classifiersvm.cvcoarsegammamax"] = 16;
	this->parameters.map["classifiersvm.cvcoarsegammastep"] = 17;
	this->parameters.map["classifiersvm.cvfinecstep"] = 18;
	this->parameters.map["classifiersvm.cvfinecradius"] = 19;
	this->parameters.map["classifiersvm.cvfinegammastep"] = 20;
	this->parameters.map["classifiersvm.cvfinegammaradius"] = 21;
	this->parameters.map["classifiersvm.degree"] = 22;
	this->parameters.map["classifiersvm.docv"] = 23;
	this->parameters.map["classifiersvm.eps"] = 24; 
	this->parameters.map["classifiersvm.gamma"] = 25;
	this->parameters.map["classifiersvm.kernel_type"] = 26;
	this->parameters.map["classifiersvm.nu"] = 27;
	this->parameters.map["classifiersvm.p"] = 28;
	this->parameters.map["classifiersvm.probability"] = 29;
	this->parameters.map["classifiersvm.shrinking"] = 30;
	this->parameters.map["classifiersvm.svm_type"] = 31;
	this->parameters.map["classifiersvm.oneclassgammasearchmaxiterations"] = 32;
	this->parameters.map["classifiersvm.oneclassgammasearchrepeats"] = 33;
	this->parameters.map["debug.debugdir"] = 34;
	this->parameters.map["cluster.mergemode"] = 35; 
	this->parameters.map["debug.usedebugprefix"] = 36;
	this->parameters.map["debug.debugmode"] = 37;
	this->parameters.map["distancemeasure.bicpenaltyfactor"] = 38;
	//tata, one number missed...
	this->parameters.map["speakerclusterer.distancemeasure"] = 40;
	this->parameters.map["speakerclusterer.globalcriterion"] = 41;
	this->parameters.map["distancemeasure.mergemode"] = 42;
	this->parameters.map["speakerclusterer.terminationcriterion"] = 43;
	this->parameters.map["distancemeasure.wcdpenaltyfactor"] = 44;
	this->parameters.map["emd.debuglevel"] = 45; 
	this->parameters.map["emd.maxiterations"] = 46;
	this->parameters.map["emd.maxsigsize"] = 47;
	this->parameters.map["enhancement.doenhancement"] = 48;
	this->parameters.map["enhancement.minnoiseduration"] = 49;
	this->parameters.map["enhancement.noisemodelupdaterate"] = 50;
	this->parameters.map["enhancement.speechmodelfile"] = 51;
	this->parameters.map["featurebandperiodicity.fftsize"] = 52;
	this->parameters.map["featurebandperiodicity.framesize"] = 53;
	this->parameters.map["featurebandperiodicity.framestep"] = 54;
	this->parameters.map["featurebandperiodicity.highcut"] = 55;
	this->parameters.map["featurebandperiodicity.lowcut"] = 56;
	this->parameters.map["featurebandperiodicity.preemphasizefactor"] = 57;
	this->parameters.map["featurebandperiodicity.window"] = 58;
	this->parameters.map["featurebrightnessbandwidth.fftsize"] = 59;
	this->parameters.map["featurebrightnessbandwidth.framesize"] = 60;
	this->parameters.map["featurebrightnessbandwidth.framestep"] = 61;
	this->parameters.map["featurebrightnessbandwidth.highcut"] = 62;
	this->parameters.map["featurebrightnessbandwidth.lowcut"] = 63;
	this->parameters.map["featurebrightnessbandwidth.preemphasizefactor"] = 64;
	this->parameters.map["featurebrightnessbandwidth.window"] = 65;
	this->parameters.map["featurefbe.adddeltadeltas"] = 66;
	this->parameters.map["featurefbe.adddeltas"] = 67;
	this->parameters.map["featurefbe.cmn"] = 68;
	this->parameters.map["featurefbe.denergy"] = 69;
	this->parameters.map["featurefbe.fftsize"] = 70;
	this->parameters.map["featurefbe.filterbanksize"] = 71;
	this->parameters.map["featurefbe.framesize"] = 72;
	this->parameters.map["featurefbe.framestep"] = 73;
	this->parameters.map["featurefbe.frequencyscale"] = 74;
	this->parameters.map["featurefbe.highcut"] = 75;
	this->parameters.map["featurefbe.lowcut"] = 76;
	this->parameters.map["featurefbe.mfcccoeffselection"] = 77;
	this->parameters.map["featurefbe.mfccorder"] = 78;
	this->parameters.map["featurefbe.preemphasizefactor"] = 79;
	this->parameters.map["featurefbe.smoothing"] = 80;
	this->parameters.map["featurefbe.window"] = 81;
	this->parameters.map["featurefbe.minfilterbankfrequency"] = 82;
	this->parameters.map["featurefbe.maxfilterbankfrequency"] = 83;
	this->parameters.map["featurefbe.resulttype"] = 84;
	this->parameters.map["featuremfcc.adddeltadeltas"] = 85;
	this->parameters.map["featuremfcc.adddeltas"] = 86;
	this->parameters.map["featuremfcc.cmn"] = 87;
	this->parameters.map["featuremfcc.denergy"] = 88;
	this->parameters.map["featuremfcc.fftsize"] = 89;
	this->parameters.map["featuremfcc.filterbanksize"] = 90;
	this->parameters.map["featuremfcc.framesize"] = 91;
	this->parameters.map["featuremfcc.framestep"] = 92;
	this->parameters.map["featuremfcc.highcut"] = 93;
	this->parameters.map["featuremfcc.lowcut"] = 94;
	this->parameters.map["featuremfcc.coeffselection"] = 95;
	this->parameters.map["featuremfcc.mfccorder"] = 96;
	this->parameters.map["featuremfcc.preemphasizefactor"] = 97;
	this->parameters.map["featuremfcc.window"] = 98;
	this->parameters.map["featuremfcc.method"] = 99;
	this->parameters.map["featuremfcc.sclib_frequencyscale"] = 100;
	this->parameters.map["featuremfcc.sclib_minfilterbankfrequency"] = 101;
	this->parameters.map["featuremfcc.sclib_maxfilterbankfrequency"] = 102;
	this->parameters.map["featuremfcc.sclib_smoothing"] = 103;
	this->parameters.map["featurenfr.fftsize"] = 104;
	this->parameters.map["featurenfr.framesize"] = 105;
	this->parameters.map["featurenfr.framestep"] = 106;
	this->parameters.map["featurenfr.highcut"] = 107;
	this->parameters.map["featurenfr.lowcut"] = 108;
	this->parameters.map["featurenfr.nfrthreshold"] = 109;
	this->parameters.map["featurenfr.preemphasizefactor"] = 110;
	this->parameters.map["featurenfr.window"] = 111;
	this->parameters.map["featurespectrum.fftsize"] = 112;
	this->parameters.map["featurespectrum.framesize"] = 113;
	this->parameters.map["featurespectrum.framestep"] = 114;
	this->parameters.map["featurespectrum.highcut"] = 115;
	this->parameters.map["featurespectrum.lowcut"] = 116;
	this->parameters.map["featurespectrum.preemphasizefactor"] = 117;
	this->parameters.map["featurespectrum.window"] = 118;
	this->parameters.map["featurespectrumflux.fftsize"] = 119;
	this->parameters.map["featurespectrumflux.framesize"] = 120;
	this->parameters.map["featurespectrumflux.framestep"] = 121;
	this->parameters.map["featurespectrumflux.highcut"] = 122;
	this->parameters.map["featurespectrumflux.lowcut"] = 123;
	this->parameters.map["featurespectrumflux.preemphasizefactor"] = 124;
	this->parameters.map["featurespectrumflux.window"] = 125;
	this->parameters.map["featuresubbandpower.fftsize"] = 126;
	this->parameters.map["featuresubbandpower.framesize"] = 127;
	this->parameters.map["featuresubbandpower.framestep"] = 128;
	this->parameters.map["featuresubbandpower.highcut"] = 129;
	this->parameters.map["featuresubbandpower.lowcut"] = 130;
	this->parameters.map["featuresubbandpower.preemphasizefactor"] = 131;
	this->parameters.map["featuresubbandpower.window"] = 132;
	this->parameters.map["featurezcr.framesize"] = 133;
	this->parameters.map["featurezcr.framestep"] = 134;
	this->parameters.map["featurezcr.highcut"] = 135;
	this->parameters.map["featurezcr.lowcut"] = 136;
	this->parameters.map["featurezcr.preemphasizefactor"] = 137;
	this->parameters.map["featurezcr.usechebyshev"] = 138;
	this->parameters.map["featurezcr.scaleresult"] = 139;
	this->parameters.map["featureste.framesize"] = 140;
	this->parameters.map["featureste.framestep"] = 141;
	this->parameters.map["featureste.highcut"] = 142;
	this->parameters.map["featureste.lowcut"] = 143;
	this->parameters.map["featureste.preemphasizefactor"] = 144;
	this->parameters.map["featureste.usebutterworth"] = 145;
	this->parameters.map["featureste.scaleresult"] = 146;
	//this->parameters.map["featurecepstralpeak.fftsize"] = 147;
	//this->parameters.map["featurecepstralpeak.framesize"] = 148;
	//this->parameters.map["featurecepstralpeak.framestep"] = 149;
	//this->parameters.map["featurecepstralpeak.highcut"] = 150;
	//this->parameters.map["featurecepstralpeak.lowcut"] = 151;
	//this->parameters.map["featurecepstralpeak.preemphasizefactor"] = 152;
	//this->parameters.map["featurecepstralpeak.lowerlimitofpickingrange"] = 153;
	//this->parameters.map["featurecepstralpeak.upperlimitofpickingrange"] = 154;
	//this->parameters.map["featurecepstralpeak.linearweightingsize"] = 155;
	//this->parameters.map["featurecepstralpeak.window"] = 156;
	//this->parameters.map["featurecepstralpeak.usebutterworth"] = 157;
	//this->parameters.map["featurewaveletenergydistribution.framesize"] = 158;
	//this->parameters.map["featurewaveletenergydistribution.framestep"] = 159;
	//this->parameters.map["featurewaveletenergydistribution.highcut"] = 160;
	//this->parameters.map["featurewaveletenergydistribution.lowcut"] = 161;
	//this->parameters.map["featurewaveletenergydistribution.preemphasizefactor"] = 162;
	//this->parameters.map["featurewaveletenergydistribution.level"] = 163;
	this->parameters.map["featurelpc.lowcut"] = 164;
	this->parameters.map["featurelpc.highcut"] = 165;
	this->parameters.map["featurelpc.framesize"] = 166;
	this->parameters.map["featurelpc.framestep"] = 167;
	this->parameters.map["featurelpc.preemphasizefactor"] = 168;
	this->parameters.map["featurelpc.lpcorder"] = 169;
	this->parameters.map["featurelpcresidual.lowcut"] = 170;
	this->parameters.map["featurelpcresidual.highcut"] = 171;
	this->parameters.map["featurelpcresidual.framesize"] = 172;
	this->parameters.map["featurelpcresidual.framestep"] = 173;
	this->parameters.map["featurelpcresidual.preemphasizefactor"] = 174;
	this->parameters.map["featurelpcresidual.order"] = 175;
	this->parameters.map["featurelsp.lowcut"] = 176;
	this->parameters.map["featurelsp.highcut"] = 177;
	this->parameters.map["featurelsp.framesize"] = 178;
	this->parameters.map["featurelsp.framestep"] = 179;
	this->parameters.map["featurelsp.preemphasizefactor"] = 180;
	this->parameters.map["featurelsp.lpcorder"] = 181;
	this->parameters.map["featurelsp.minseparation"] = 182;
	this->parameters.map["featurelsp.bisections"] = 183;
	this->parameters.map["featurelsp.delta"] = 184;
	this->parameters.map["featurelsp.method"] = 185;
	this->parameters.map["featurelsp.maxloops"] = 186;
	this->parameters.map["featurepitch.framesize"] = 187;
	this->parameters.map["featurepitch.framestep"] = 188;
	this->parameters.map["featurepitch.method"] = 189;
	this->parameters.map["featurepitch.esps_cand_thresh"] = 190;
	this->parameters.map["featurepitch.esps_lag_weight"] = 191;
	this->parameters.map["featurepitch.esps_freq_weight"] = 192;
	this->parameters.map["featurepitch.esps_trans_cost"] = 193;
	this->parameters.map["featurepitch.esps_trans_amp"] = 194;
	this->parameters.map["featurepitch.esps_trans_spec"] = 195;
	this->parameters.map["featurepitch.esps_voice_bias"] = 196;
	this->parameters.map["featurepitch.esps_double_cost"] = 197;
	this->parameters.map["featurepitch.esps_min_f0"] = 198;
	this->parameters.map["featurepitch.esps_max_f0"] = 199;
	this->parameters.map["featurepitch.esps_n_cands"] = 200;
	this->parameters.map["featurepitch.esps_wind_dur"] = 201;
	this->parameters.map["featuresdp.color"] = 202;
	this->parameters.map["featuresdp.framesize"] = 203; 
	this->parameters.map["featuresdp.framestep"] = 204;
	this->parameters.map["featuresdp.lag"] = 205;
	this->parameters.map["featuresdp.m"] = 206;
	this->parameters.map["featuresdp.n"] = 207;
	this->parameters.map["featuresdp.picturesize"] = 208;
	this->parameters.map["featuresdp.tau"] = 209;
	this->parameters.map["featuresdp.preemphasizefactor"] = 210;
	this->parameters.map["general.firstscene"] = 211;
	this->parameters.map["general.lastscene"] = 212;
	this->parameters.map["general.sceneselection"] = 213;
	this->parameters.map["general.shortspeechthreshold"] = 214;
	this->parameters.map["general.pausesilencethreshold"] = 215;
	this->parameters.map["general.preclusteringresultsprefix"] = 216;
	this->parameters.map["general.featureprefix"] = 217;
	this->parameters.map["groundtruth.internalframesize"] = 218;
	this->parameters.map["groundtruth.pseudoscenelength"] = 219;
	this->parameters.map["groundtruth.videoframemachineoffset"] = 220;
	this->parameters.map["groundtruth.storeprobabilityinformation"] = 221;
	this->parameters.map["mixturemodelbgmm.emthreshold"] = 222;
	this->parameters.map["mixturemodelbgmm.maxemiterations"] = 223;
	this->parameters.map["mixturemodelbgmm.variancelimit"] = 224;
	this->parameters.map["mixturemodelbgmm.weightlimit"] = 225;
	this->parameters.map["mixturemodelgmm.emthreshold"] = 226;
	this->parameters.map["mixturemodelgmm.maxemiterations"] = 227;
	this->parameters.map["mixturemodelgmm.kmeansiterations"] = 228;
	this->parameters.map["mixturemodelgmm.variancelimit"] = 229;
	this->parameters.map["mixturemodelgmm.weightlimit"] = 230;
	this->parameters.map["mixturemodelmixmax.bgmodelcombination"] = 231;
	this->parameters.map["mixturemodelmixmax.emthreshold"] = 232;
	this->parameters.map["mixturemodelmixmax.kmeansiterations"] = 233; 
	this->parameters.map["mixturemodelmixmax.maxemiterations"] = 234; 
	this->parameters.map["mixturemodelmixmax.maxerfloops"] = 235; 
	this->parameters.map["mixturemodelmixmax.noisecorruptiontype"] = 236; 
	this->parameters.map["mixturemodelmixmax.variancelimit"] = 237;
	this->parameters.map["mixturemodelmixmax.weightlimit"] = 238;
	this->parameters.map["mixturemodelgmmubm.adaptmeans"] = 239;
	this->parameters.map["mixturemodelgmmubm.adaptvariances"] = 240;
	this->parameters.map["mixturemodelgmmubm.adaptweights"] = 241;
	this->parameters.map["mixturemodelgmmubm.relevancefactor"] = 242;
	this->parameters.map["mixturemodelgmmubm.scoringmethod"] = 243;
	this->parameters.map["mixturemodelgmmubm.topcmixtures"] = 244;
	this->parameters.map["mixturemodelgmmubm.ubmfilename"] = 245;
	this->parameters.map["mixturemodelgmmubm.variancelimit"] = 246; 
	this->parameters.map["mixturemodelmix2max.emthreshold"] = 247;
	this->parameters.map["mixturemodelmix2max.kmeansiterations"] = 248;
	this->parameters.map["mixturemodelmix2max.maxemiterations"] = 249;
	this->parameters.map["mixturemodelmix2max.maxerfloops"] = 250;
	this->parameters.map["mixturemodelmix2max.variancelimit"] = 251;
	this->parameters.map["mixturemodelmix2max.weightlimit"] = 252;
	this->parameters.map["mixturemodelmix2max.bgmodelcombination"] = 253;
	this->parameters.map["mixturemodelmix2maxex.bgmodelcombination"] = 254;
	this->parameters.map["mixturemodelmix2maxex.emthreshold"] = 255;
	this->parameters.map["mixturemodelmix2maxex.kmeansiterations"] = 256;
	this->parameters.map["mixturemodelmix2maxex.maxemiterations"] = 257;
	this->parameters.map["mixturemodelmix2maxex.maxerfloops"] = 258;
	this->parameters.map["mixturemodelmix2maxex.variancelimit"] = 259;
	this->parameters.map["mixturemodelmix2maxex.weightlimit"] = 260;
	this->parameters.map["modelhandler.backgroundmodeltype"] = 261;
	//this->parameters.map["modelann.desirederror"] = 262;
	//this->parameters.map["modelann.learningrate"] = 263;
	//this->parameters.map["modelann.maxepochs"] = 264;
	//this->parameters.map["modelann.networkstructure"] = 265;
	//this->parameters.map["modelann.patternmode"] = 266;
	//this->parameters.map["modelann.temperature"] = 267;
	//this->parameters.map["modelann.useannealing"] = 268;
	this->parameters.map["modelhandler.foregroundmodeltype"] = 269;
	this->parameters.map["modelhandler.maxnoisemodelorder"] = 270;
	this->parameters.map["modelhandler.maxspeakermodelorder"] = 271;
	this->parameters.map["modelhandler.mspergaussian"] = 272;
	this->parameters.map["modelhandler.orderguessemsteps"] = 273;
	this->parameters.map["modelhandler.orderguessmode"] = 274;
	this->parameters.map["modelhandler.snrthreshold"] = 275;
	this->parameters.map["modelhandler.onlythisspeaker"] = 276;
	this->parameters.map["modelhandler.speakermodelfeature"] = 277;
	this->parameters.map["modelpareto.usemarginaldistributions"] = 278;
	this->parameters.map["modelqgmm.deltabiclambda"] = 279;
	this->parameters.map["modelqgmm.maxmixtures"] = 280;
	this->parameters.map["modelqgmm.percentdifference"] = 281;
	this->parameters.map["resampling.fastconversion"] = 282;
	this->parameters.map["score.bbnmetriclambda"] = 283;
	this->parameters.map["segmentationaudiotypelzl.actionmodelfilename"] = 284;
	this->parameters.map["segmentationaudiotypelzl.classifierfilename"] = 285;
	this->parameters.map["segmentationaudiotypelzl.featurefilename"] = 286;
	this->parameters.map["segmentationaudiotypelzl.normalizationfilename"] = 287;
	this->parameters.map["segmentationaudiotypelzl.subcliplength"] = 288;
	this->parameters.map["segmentationaudiotypelzl.bbparameters.fftsize"] = 289;
	this->parameters.map["segmentationaudiotypelzl.bbparameters.framesize"] = 290;
	this->parameters.map["segmentationaudiotypelzl.bbparameters.framestep"] = 291;
	this->parameters.map["segmentationaudiotypelzl.bbparameters.preemphasizefactor"] = 292;
	this->parameters.map["segmentationaudiotypelzl.bbparameters.window"] = 293;
	this->parameters.map["segmentationaudiotypelzl.bbparameters.samplerate"] = 294;
	this->parameters.map["segmentationaudiotypelzl.bpparameters.fftsize"] = 295;
	this->parameters.map["segmentationaudiotypelzl.bpparameters.framesize"] = 296;
	this->parameters.map["segmentationaudiotypelzl.bpparameters.framestep"] = 297;
	this->parameters.map["segmentationaudiotypelzl.bpparameters.preemphasizefactor"] = 298;
	this->parameters.map["segmentationaudiotypelzl.bpparameters.window"] = 299;
	this->parameters.map["segmentationaudiotypelzl.bpparameters.samplerate"] = 300;
	this->parameters.map["segmentationaudiotypelzl.mfccparameters.fftsize"] = 301;
	this->parameters.map["segmentationaudiotypelzl.mfccparameters.framesize"] = 302;
	this->parameters.map["segmentationaudiotypelzl.mfccparameters.framestep"] = 303;
	this->parameters.map["segmentationaudiotypelzl.mfccparameters.mfccorder"] = 304; 
	this->parameters.map["segmentationaudiotypelzl.mfccparameters.adddeltadeltas"] = 305;
	this->parameters.map["segmentationaudiotypelzl.mfccparameters.adddeltas"] = 306;
	this->parameters.map["segmentationaudiotypelzl.mfccparameters.method"] = 307;
	this->parameters.map["segmentationaudiotypelzl.mfccparameters.denergy"] = 308;
	this->parameters.map["segmentationaudiotypelzl.mfccparameters.cmn"] = 309;
	this->parameters.map["segmentationaudiotypelzl.mfccparameters.filterbanksize"] = 310;
	this->parameters.map["segmentationaudiotypelzl.mfccparameters.preemphasizefactor"] = 311;
	this->parameters.map["segmentationaudiotypelzl.mfccparameters.sclib_frequencyscale"] = 312;
	this->parameters.map["segmentationaudiotypelzl.mfccparameters.sclib_smoothing"] = 313;
	this->parameters.map["segmentationaudiotypelzl.mfccparameters.sclib_maxfilterbankfrequency"] = 314;
	this->parameters.map["segmentationaudiotypelzl.mfccparameters.sclib_minfilterbankfrequency"] = 315;
	this->parameters.map["segmentationaudiotypelzl.mfccparameters.window"] = 316;
	this->parameters.map["segmentationaudiotypelzl.mfccparameters.coeffselection"] = 317;
	//this->parameters.map["segmentationaudiotypelzl.mfccparameters.coeffselection"] = 318;
	this->parameters.map["segmentationaudiotypelzl.mfccparameters.samplerate"] = 319;
	this->parameters.map["segmentationaudiotypelzl.nfrparameters.fftsize"] = 320;
	this->parameters.map["segmentationaudiotypelzl.nfrparameters.framesize"] = 321;
	this->parameters.map["segmentationaudiotypelzl.nfrparameters.framestep"] = 322;
	this->parameters.map["segmentationaudiotypelzl.nfrparameters.preemphasizefactor"] = 323;
	this->parameters.map["segmentationaudiotypelzl.nfrparameters.window"] = 324;
	this->parameters.map["segmentationaudiotypelzl.nfrparameters.nfrthreshold"] = 325;
	this->parameters.map["segmentationaudiotypelzl.nfrparameters.samplerate"] = 326;
	this->parameters.map["segmentationaudiotypelzl.sbpparameters.fftsize"] = 327;
	this->parameters.map["segmentationaudiotypelzl.sbpparameters.framesize"] = 328;
	this->parameters.map["segmentationaudiotypelzl.sbpparameters.framestep"] = 329;
	this->parameters.map["segmentationaudiotypelzl.sbpparameters.preemphasizefactor"] = 330;
	this->parameters.map["segmentationaudiotypelzl.sbpparameters.window"] = 331;
	this->parameters.map["segmentationaudiotypelzl.sbpparameters.samplerate"] = 332;
	this->parameters.map["segmentationaudiotypelzl.sfparameters.fftsize"] = 333;
	this->parameters.map["segmentationaudiotypelzl.sfparameters.framesize"] = 334;
	this->parameters.map["segmentationaudiotypelzl.sfparameters.framestep"] = 335;
	this->parameters.map["segmentationaudiotypelzl.sfparameters.preemphasizefactor"] = 336;
	this->parameters.map["segmentationaudiotypelzl.sfparameters.window"] = 337;
	this->parameters.map["segmentationaudiotypelzl.sfparameters.samplerate"] = 338;
	this->parameters.map["segmentationaudiotypelzl.zcrparameters.usechebyshev"] = 339;
	this->parameters.map["segmentationaudiotypelzl.zcrparameters.framesize"] = 340;
	this->parameters.map["segmentationaudiotypelzl.zcrparameters.framestep"] = 341;
	this->parameters.map["segmentationaudiotypelzl.zcrparameters.preemphasizefactor"] = 342;
	this->parameters.map["segmentationaudiotypelzl.zcrparameters.scaleresult"] = 343;
	this->parameters.map["segmentationaudiotypelzl.zcrparameters.samplerate"] = 344;
	this->parameters.map["segmentationchangeslz.adaptivethresholdalpha"] = 345;
	this->parameters.map["segmentationchangeslz.bayesianthreshold"] = 346;
	this->parameters.map["segmentationchangeslz.detectorwindowlength"] = 347;
	this->parameters.map["segmentationchangeslz.detectorwindowstep"] = 348;
	this->parameters.map["segmentationchangeslz.lastndistances"] = 349;
	this->parameters.map["segmentationchangeslz.priorsfilename"] = 350;
	this->parameters.map["segmentationchangeslz.time2changemodelfilename"] = 351;
	//this->parameters.map["segmentationchangestst.smallwindowlength"] = 352;
	//this->parameters.map["segmentationchangestst.mediumwindowlength"] = 353;
	//this->parameters.map["segmentationchangestst.largewindowlength"] = 354;
	//this->parameters.map["segmentationchangestst.windowstep"] = 355;
	//this->parameters.map["segmentationchangestst.smallclassifierfilename"] = 356;
	//this->parameters.map["segmentationchangestst.mediumclassifierfilename"] = 357;
	//this->parameters.map["segmentationchangestst.largeclassifierfilename"] = 358;
	//this->parameters.map["segmentationchangestst.whichsizetotrain"] = 359;
	//this->parameters.map["segmentationchangestst.positiveexamplesfilename"] = 360;
	//this->parameters.map["segmentationchangestst.negativeexamplesfilename"] = 361;
	//this->parameters.map["segmentationchangestst.msforadpativethresh"] = 362;
	//this->parameters.map["segmentationchangestst.posexcpaccuracydemand"] = 363;
	//this->parameters.map["segmentationchangestst.createadaptiveexamples"] = 364;
	//this->parameters.map["segmentationchangestst.fuserfilename"] = 365;
	//this->parameters.map["segmentationchangestst.msbetweenpeaks"] = 366;
	//this->parameters.map["segmentationchangestst.synthesizemoreexamples"] = 367;
	//this->parameters.map["segmentationchangestst.usepeakdetection"] = 368;
	//this->parameters.map["segmentationchangestst.featureconcatenationmode"] = 369;
	//this->parameters.map["segmentationchangestst.usegapforaverage"] = 370;
	//this->parameters.map["segmentationchangestst.pos2negcostratio"] = 371;
	//this->parameters.map["segmentationchangestst.smallwindowoverlap"] = 372;
	//this->parameters.map["segmentationchangestst.mediumwindowoverlap"] = 373;
	//this->parameters.map["segmentationchangestst.largewindowoverlap"] = 374;
	//this->parameters.map["segmentationchangestst.rectanglefeatures"] = 375;
	//this->parameters.map["segmentationchangestst.segdisfeatures"] = 376;
	this->parameters.map["segmentationhandler.audiotypemode"] = 377;
	this->parameters.map["segmentationhandler.changedetectormode"] = 378;
	this->parameters.map["segmentationhandler.silencedetectormode"] = 379;
	this->parameters.map["segmentationhandler.vuvdetectormode"] = 380;
	this->parameters.map["segmentationsilencelnk.energyquantizationlevel"] = 381;
	this->parameters.map["segmentationsilencelzl.energysilencethreshold"] = 382;
	this->parameters.map["segmentationsilencelzl.zcrsilencethreshold"] = 383;
	this->parameters.map["segmentationvuvlnk.energyquantizationlevel"] = 384;
	//this->parameters.map["segmentationvuvsam.framesize"] = 385;
	//this->parameters.map["segmentationvuvsam.period"] = 386;
	//this->parameters.map["segmentationvuvsam.firlength"] = 387;
	//this->parameters.map["segmentationvuvsam.fircount"] = 388;
	//this->parameters.map["segmentationvuvsam.fftlength"] = 389;
	//this->parameters.map["segmentationvuvsam.firfirstcenter"] = 390;
	//this->parameters.map["segmentationvuvsam.firbandwidth"] = 391;
	//this->parameters.map["segmentationvuvsisy.classifiertype"] = 392;
	//this->parameters.map["segmentationvuvsisy.classifiersel"] = 393;
	//this->parameters.map["segmentationvuvsisy.lpcclassifier"] = 394;
	//this->parameters.map["segmentationvuvsisy.mfccclassifier"] = 395;
	//this->parameters.map["segmentationvuvsisy.transmatrix"] = 396;
	//this->parameters.map["segmentationvuvas.intervalfortracking"] = 397;
	//this->parameters.map["segmentationvuvas.zcrparameters.framesize"] = 398;
	//this->parameters.map["segmentationvuvas.zcrparameters.framestep"] = 399;
	//this->parameters.map["segmentationvuvas.zcrparameters.preemphasizefactor"] = 400;
	//this->parameters.map["segmentationvuvas.zcrparameters.scaleresult"] = 401;
	//this->parameters.map["segmentationvuvas.zcrparameters.usechebyshev"] = 402;
	//this->parameters.map["segmentationvuvas.zcrparameters.lowcut"] = 403;
	//this->parameters.map["segmentationvuvas.zcrparameters.highcut"] = 404;
	//this->parameters.map["segmentationvuvas.steparameters.framesize"] = 405;
	//this->parameters.map["segmentationvuvas.steparameters.framestep"] = 406;
	//this->parameters.map["segmentationvuvas.steparameters.preemphasizefactor"] = 407;
	//this->parameters.map["segmentationvuvas.steparameters.scaleresult"] = 408;
	//this->parameters.map["segmentationvuvas.steparameters.usebutterworth"] = 409;
	//this->parameters.map["segmentationvuvas.steparameters.lowcut"] = 410;
	//this->parameters.map["segmentationvuvas.steparameters.highcut"] = 411;
	this->parameters.map["clusterer.maxiterations"] = 412;
	this->parameters.map["clusterer.numcluster"] = 413;
	this->parameters.map["clusterer.iterationstoinikmeanlist"] = 414;
	this->parameters.map["signalhandler.forcesamplerate"] = 415;
	this->parameters.map["signalmpeg.fastseeking"] = 416;
	this->parameters.map["signalmpeg.hqresampling"] = 417;
	this->parameters.map["signalmpeg.outputchannelcount"] = 418;
	this->parameters.map["signalmpeg.outputsamplerate"] = 419;
	this->parameters.map["speakerclusterer.doclustering"] = 420;
	this->parameters.map["speakerclusterer.speechseglengththreshold"] = 421;
	this->parameters.map["speakeridentification.dnormsamplecount"] = 422;
	this->parameters.map["speakeridentification.doidentification"] = 423;
	this->parameters.map["speakeridentification.normalizationmode"] = 424;
	this->parameters.map["speakeridentification.useubms"] = 425;
	this->parameters.map["transform.taperinglength"] = 426;
	this->parameters.map["segmentationaudiotypelzl.svmparameters.c"] = 427;
	this->parameters.map["segmentationaudiotypelzl.svmparameters.cache_size"] = 428;
	this->parameters.map["segmentationaudiotypelzl.svmparameters.coef0"] = 429;
	this->parameters.map["segmentationaudiotypelzl.svmparameters.cvfolds"] = 430;
	this->parameters.map["segmentationaudiotypelzl.svmparameters.cvmaxdatasetsize"] = 431;
	this->parameters.map["segmentationaudiotypelzl.svmparameters.cvcoarsecmin"] = 432;
	this->parameters.map["segmentationaudiotypelzl.svmparameters.cvcoarsecmax"] = 433; 
	this->parameters.map["segmentationaudiotypelzl.svmparameters.cvcoarsecstep"] = 434; 
	this->parameters.map["segmentationaudiotypelzl.svmparameters.cvcoarsegammamin"] = 435; 
	this->parameters.map["segmentationaudiotypelzl.svmparameters.cvcoarsegammamax"] = 436;
	this->parameters.map["segmentationaudiotypelzl.svmparameters.cvcoarsegammastep"] = 437;
	this->parameters.map["segmentationaudiotypelzl.svmparameters.cvfinecstep"] = 438;
	this->parameters.map["segmentationaudiotypelzl.svmparameters.cvfinecradius"] = 439;
	this->parameters.map["segmentationaudiotypelzl.svmparameters.cvfinegammastep"] = 440;
	this->parameters.map["segmentationaudiotypelzl.svmparameters.cvfinegammaradius"] = 441;
	this->parameters.map["segmentationaudiotypelzl.svmparameters.degree"] = 442;
	this->parameters.map["segmentationaudiotypelzl.svmparameters.docv"] = 443;
	this->parameters.map["segmentationaudiotypelzl.svmparameters.eps"] = 444; 
	this->parameters.map["segmentationaudiotypelzl.svmparameters.gamma"] = 445;
	this->parameters.map["segmentationaudiotypelzl.svmparameters.kernel_type"] = 446;
	this->parameters.map["segmentationaudiotypelzl.svmparameters.nu"] = 447;
	this->parameters.map["segmentationaudiotypelzl.svmparameters.p"] = 448;
	this->parameters.map["segmentationaudiotypelzl.svmparameters.probability"] = 449;
	this->parameters.map["segmentationaudiotypelzl.svmparameters.shrinking"] = 450;
	this->parameters.map["segmentationaudiotypelzl.svmparameters.svm_type"] = 451;
	this->parameters.map["segmentationaudiotypelzl.svmparameters.oneclassgammasearchmaxiterations"] = 452;
	this->parameters.map["segmentationaudiotypelzl.svmparameters.oneclassgammasearchrepeats"] = 453;
	this->parameters.map["segmentationsilencelzl.zcrparameters.framesize"] = 454;
	this->parameters.map["segmentationsilencelzl.zcrparameters.framestep"] = 455;
	this->parameters.map["segmentationsilencelzl.zcrparameters.samplerate"] = 456;
	this->parameters.map["segmentationsilencelzl.zcrparameters.preemphasizefactor"] = 457;
	this->parameters.map["segmentationsilencelzl.zcrparameters.usechebyshev"] = 458;
	this->parameters.map["segmentationsilencelzl.zcrparameters.scaleresult"] = 459;
	this->parameters.map["segmentationsilencelzl.steparameters.framesize"] = 460;
	this->parameters.map["segmentationsilencelzl.steparameters.framestep"] = 461;
	this->parameters.map["segmentationsilencelzl.steparameters.samplerate"] = 462;
	this->parameters.map["segmentationsilencelzl.steparameters.preemphasizefactor"] = 463;
	this->parameters.map["segmentationsilencelzl.steparameters.usebutterworth"] = 464;
	this->parameters.map["segmentationsilencelzl.steparameters.scaleresult"] = 465;
	this->parameters.map["segmentationvuvesps.pitchparameters.framesize"] = 466;
	this->parameters.map["segmentationvuvesps.pitchparameters.framestep"] = 467;
	this->parameters.map["segmentationvuvesps.pitchparameters.method"] = 468;
	this->parameters.map["segmentationvuvesps.pitchparameters.esps_cand_thresh"] = 469;
	this->parameters.map["segmentationvuvesps.pitchparameters.esps_lag_weight"] = 470;
	this->parameters.map["segmentationvuvesps.pitchparameters.esps_freq_weight"] = 471;
	this->parameters.map["segmentationvuvesps.pitchparameters.esps_trans_cost"] = 472;
	this->parameters.map["segmentationvuvesps.pitchparameters.esps_trans_amp"] = 473;
	this->parameters.map["segmentationvuvesps.pitchparameters.esps_trans_spec"] = 474;
	this->parameters.map["segmentationvuvesps.pitchparameters.esps_voice_bias"] = 475;
	this->parameters.map["segmentationvuvesps.pitchparameters.esps_double_cost"] = 476;
	this->parameters.map["segmentationvuvesps.pitchparameters.esps_min_f0"] = 477;
	this->parameters.map["segmentationvuvesps.pitchparameters.esps_max_f0"] = 478;
	this->parameters.map["segmentationvuvesps.pitchparameters.esps_n_cands"] = 479;
	this->parameters.map["segmentationvuvesps.pitchparameters.esps_wind_dur"] = 480;
	//this->parameters.map["transform.olaiterationcount"] = 481;
	this->parameters.map["featurespectrum.logarithmize"] = 482;
	this->parameters.map["featurespectrum.createphase"] = 483;
	this->parameters.map["modelhmm.statecount"] = 484;
	this->parameters.map["modelhmm.transitionstructure"] = 485;
	this->parameters.map["modelhmm.mixturesperstate"] = 486;
	this->parameters.map["modelhmm.useorthogonaltransform"] = 487;
	this->parameters.map["modelhmm.lefttoright"] = 488;
	this->parameters.map["modelhmm.maxiterations"] = 489;
	this->parameters.map["modelhmm.verbose"] = 490;
	this->parameters.map["mixturemodelbgmm.fullcovariance"] = 491;
	this->parameters.map["mixturemodelbgmm.randominitialization"] = 492;
	this->parameters.map["featurelpc.window"] = 493;
	this->parameters.map["featurelsp.window"] = 494;
	this->parameters.map["featureformant.framesize"] = 495;
	this->parameters.map["featureformant.framestep"] = 496;
	this->parameters.map["featureformant.lowcut"] = 497;
	this->parameters.map["featureformant.highcut"] = 498;
	this->parameters.map["featureformant.preemphasizefactor"] = 499;
	this->parameters.map["featureformant.samplerate"] = 500;
	this->parameters.map["featureformant.esps_lpc_ord"] = 501;
	this->parameters.map["featureformant.esps_lpc_type"] = 502;
	this->parameters.map["featureformant.esps_w_type"] = 503;
	this->parameters.map["featureformant.esps_ds_freq"] = 504;
	this->parameters.map["featureformant.esps_wdur"] = 505;
	this->parameters.map["featureformant.esps_nom_f1"] = 506;
	this->parameters.map["featureformant.esps_cor_wdur"] = 507;
	this->parameters.map["featureformant.esps_frame_int"] = 508;
	this->parameters.map["featureformant.esps_nform"] = 509;
	this->parameters.map["featurelpc.computegain"] = 510;
	this->parameters.map["featurewaveletenergydistribution.daubechiesorder"] = 511;
	//this->parameters.map["segmentationvuvjm.wedthreshold"] = 512;
	//this->parameters.map["segmentationvuvjm.zcrparameters.samplerate"] = 513;
	//this->parameters.map["segmentationvuvjm.zcrparameters.framesize"] = 514;
	//this->parameters.map["segmentationvuvjm.zcrparameters.framestep"] = 515; 
	//this->parameters.map["segmentationvuvjm.zcrparameters.preemphasizefactor"] = 516;
	//this->parameters.map["segmentationvuvjm.zcrparameters.scaleresult"] = 517;
	//this->parameters.map["segmentationvuvjm.zcrparameters.usechebyshev"] = 518;
	//this->parameters.map["segmentationvuvjm.zcrparameters.lowcut"] = 519;
	//this->parameters.map["segmentationvuvjm.zcrparameters.highcut"] = 520;
	//this->parameters.map["segmentationvuvjm.wedparameters.samplerate"] = 521;
	//this->parameters.map["segmentationvuvjm.wedparameters.framesize"] = 522;
	//this->parameters.map["segmentationvuvjm.wedparameters.framestep"] = 523;
	//this->parameters.map["segmentationvuvjm.wedparameters.preemphasizefactor"] = 524;
	//this->parameters.map["segmentationvuvjm.wedparameters.level"] = 525;
	//this->parameters.map["segmentationvuvjm.wedparameters.daubechiesorder"] = 526;
	//this->parameters.map["segmentationvuvjm.wedparameters.lowcut"] = 527;
	//this->parameters.map["segmentationvuvjm.wedparameters.highcut"] = 528;
	this->parameters.map["modelvq.codebooksize"] = 529;
	this->parameters.map["modelvq.splitmethod"] = 530;
	this->parameters.map["modelvq.maxiterations"] = 531;
	//this->parameters.map["modelgroup.submodeltype"] = 532;
	//this->parameters.map["modelgroup.groupidcolumn"] = 533;
	this->parameters.map["segmentationchangeskbk.r"] = 534;
	//this->parameters.map["segmentationchangeskbk.resolutionfactor"] = 535;
	this->parameters.map["segmentationchangeskbk.lambda"] = 536;
	this->parameters.map["segmentationchangeskbk.mfccparameters.adddeltadeltas"] = 537;
	this->parameters.map["segmentationchangeskbk.mfccparameters.adddeltas"] = 538;
	this->parameters.map["segmentationchangeskbk.mfccparameters.cmn"] = 539;
	this->parameters.map["segmentationchangeskbk.mfccparameters.denergy"] = 540;
	this->parameters.map["segmentationchangeskbk.mfccparameters.fftsize"] = 541;
	this->parameters.map["segmentationchangeskbk.mfccparameters.filterbanksize"] = 542;
	this->parameters.map["segmentationchangeskbk.mfccparameters.framesize"] = 543;
	this->parameters.map["segmentationchangeskbk.mfccparameters.framestep"] = 544;
	this->parameters.map["segmentationchangeskbk.mfccparameters.highcut"] = 545;
	this->parameters.map["segmentationchangeskbk.mfccparameters.lowcut"] = 546;
	this->parameters.map["segmentationchangeskbk.mfccparameters.coeffselection"] = 547;
	this->parameters.map["segmentationchangeskbk.mfccparameters.mfccorder"] = 548;
	this->parameters.map["segmentationchangeskbk.mfccparameters.preemphasizefactor"] = 549;
	this->parameters.map["segmentationchangeskbk.mfccparameters.window"] = 550;
	this->parameters.map["segmentationchangeskbk.mfccparameters.method"] = 551;
	this->parameters.map["segmentationchangeskbk.mfccparameters.sclib_frequencyscale"] = 552;
	this->parameters.map["segmentationchangeskbk.mfccparameters.sclib_minfilterbankfrequency"] = 553;
	this->parameters.map["segmentationchangeskbk.mfccparameters.sclib_maxfilterbankfrequency"] = 554;
	this->parameters.map["segmentationchangeskbk.mfccparameters.sclib_smoothing"] = 555;
	this->parameters.map["segmentationchangeskbk.tolerance"] = 556;
	this->parameters.map["speakerclusterer.constructnonoverlappingclusters"] = 557;
	this->parameters.map["distancemeasure.grounddistance"] = 558;
	this->parameters.map["distancemeasure.icrthreshold"] = 559;
	this->parameters.map["speakerclusterer.linkagemode"] = 560;
	this->parameters.map["modelhandler.outlierremovalmode"] = 561;
	this->parameters.map["featuresamples.framesize"] = 562;
	this->parameters.map["featuresamples.framestep"] = 563;
	this->parameters.map["featuresamples.highcut"] = 564;
	this->parameters.map["featuresamples.lowcut"] = 565;
	this->parameters.map["modelsvm.distancebasedtesting"] = 566;
	this->parameters.map["modelsvm.doparametersearch"] = 567;
	this->parameters.map["modeltime.syllablelength"] = 568;
	this->parameters.map["modeltime.trajectorystep"] = 569;
	this->parameters.map["modeltime.submodeltype"] = 570;
	this->parameters.map["modeltime.removetiming"] = 571;
	this->parameters.map["modeltime.templatecount"] = 572;
	this->parameters.map["modeltime.clusteringiterations"] = 573;
	this->parameters.map["modeltime.replacetrainingdata"] = 574;
	this->parameters.map["modeltime.checkfortrajectorization"] = 575;
	this->parameters.map["modeltime.worldmodelfile"] = 576;
	this->parameters.map["modeltime.normalizationfile"] = 577;
	this->parameters.map["featurepitch.sqrt"] = 578;
	this->parameters.map["segmentationsilencelzl.specificity"] = 579;
	this->parameters.map["segmentationaudiotypelzl.speechspecificity"] = 580;
	this->parameters.map["segmentationaudiotypelzl.musicspecificity"] = 581;
	this->parameters.map["speakerclusterer.specificity"] = 582;
	this->parameters.map["speakerclusterer.firstdistancematrixprefix"] = 583;

	//ATTENTION: lcase-names above!

  //set meaningful standard-parameters first:
	this->classifierAdaBoost.maxWeakClassifiers = 25;
	this->classifierAdaBoost.weakClassifierType = sclib::ctDecisionStump;
	this->classifierDecisionStump.splitMultiway = false;
  this->classifierMl.maxOrder = 1;
  this->classifierMl.minOrder = 128; 
  this->classifierMl.modelType = sclib::mtGMM_new;
  this->classifierSvm.C = 0; //i.e. use data-dependant default
  this->classifierSvm.cache_size = 100;
  this->classifierSvm.coef0 = 0;
  this->classifierSvm.cvFolds = 3;
  this->classifierSvm.cvMaxDatasetSize = 100000;
  this->classifierSvm.cvCoarseCMin = 0.0;
  this->classifierSvm.cvCoarseCMax = 10.0;
  this->classifierSvm.cvCoarseCStep = 2.0;
  this->classifierSvm.cvCoarseGammaMin = -15.0;
  this->classifierSvm.cvCoarseGammaMax = 3.0;
  this->classifierSvm.cvCoarseGammaStep = 2.0;
  this->classifierSvm.cvFineCStep = 0.25;
  this->classifierSvm.cvFineCRadius = 2.0;
  this->classifierSvm.cvFineGammaStep = 0.25;
  this->classifierSvm.cvFineGammaRadius = 2.0;
  this->classifierSvm.degree = 3;
  this->classifierSvm.doCV = true;
  this->classifierSvm.eps = 1e-3;
  this->classifierSvm.gamma = 0; //i.e. use data-dependant default
  this->classifierSvm.kernel_type = 2;
  this->classifierSvm.nu = 0.5;
  this->classifierSvm.p = 0.1;
  this->classifierSvm.probability = 1;
  this->classifierSvm.shrinking = 1;
  this->classifierSvm.svm_type = 0;
	this->classifierSvm.oneClassGammaSearchMaxIterations = 1000000;
	this->classifierSvm.oneClassGammaSearchRepeats = 10;
  this->cluster.mergeMode = sclib::mergeAddUp;
  this->debug.debugDir = NULL;
  this->originalDebugDir = NULL;
  setDebugDir("../data/test");
  this->debug.debugMode = sclib::dbClustering;
  this->debug.useDebugPrefix = true;
  this->distanceMeasure.BICpenaltyFactor = 1.0;
  this->distanceMeasure.mergeMode = sclib::mergeAddUp;
	this->distanceMeasure.WCDpenaltyFactor = 1.0;
	this->distanceMeasure.groundDistance = sclib::dmEuclid;
	this->distanceMeasure.ICRthreshold = 0.19547; //as in the paper; needs to be tuned to actual data: compute ICR measure for several homogeneous clusters, take mean+std of the measurements as threshold
	this->emd.debugLevel = 0;
	this->emd.maxIterations = 1000;
	this->emd.maxSigSize = 1000;
	this->enhancement.doEnhancement = false;
  this->enhancement.minNoiseDuration = 250;
  this->enhancement.noiseModelUpdateRate = 250;
  MArray_1D(this->enhancement.speechModelFile, 1, char, "SC_TweakableParameters: enhancement.speechModelFile");
  sprintf(this->enhancement.speechModelFile, "\0");
  this->featureBandPeriodicity.FFTsize = 512;
  this->featureBandPeriodicity.frameSize = 25;
  this->featureBandPeriodicity.frameStep = 25;
  this->featureBandPeriodicity.highCut = 7000.0;
  this->featureBandPeriodicity.lowCut = 50.0;
  this->featureBandPeriodicity.preEmphasizeFactor = 0.97;
  this->featureBandPeriodicity.window = sclib::wndHamming;
  this->featureBrightnessBandwidth.FFTsize = 512;
  this->featureBrightnessBandwidth.frameSize = 25;
  this->featureBrightnessBandwidth.frameStep = 25;
  this->featureBrightnessBandwidth.highCut = 7000.0;
  this->featureBrightnessBandwidth.lowCut = 50.0;
  this->featureBrightnessBandwidth.preEmphasizeFactor = 0.97;
  this->featureBrightnessBandwidth.window = sclib::wndHamming;
  this->featureFbe.addDeltaDeltas = false;
  this->featureFbe.addDeltas = false;
  this->featureFbe.CMN = false;
  this->featureFbe.dEnergy = true;
  this->featureFbe.FFTsize = 512;
  this->featureFbe.filterBankSize = 24;
  this->featureFbe.frameSize = 32;
  this->featureFbe.frameStep = 16;
  this->featureFbe.frequencyScale = sclib::scaleExpoLog;
  this->featureFbe.highCut = 7000.0;
  this->featureFbe.lowCut = 50.0;
	for (char i = 0; i < 8; i++) {
		this->featureFbe.MFCCcoeffSelection[i] = (char)(0xFF);
	}
  this->featureFbe.MFCCorder = 20;
  this->featureFbe.preEmphasizeFactor = 0.97;
  this->featureFbe.smoothing = sclib::smoothNone;
  this->featureFbe.window = sclib::wndHamming;
	this->featureFbe.minFilterBankFrequency = 0.0;
	this->featureFbe.maxFilterBankFrequency = 0.0;
	this->featureFbe.resultType = sclib::resultLog;
  this->featureMfcc.addDeltaDeltas = false;
  this->featureMfcc.addDeltas = false;
  this->featureMfcc.CMN = false;
  this->featureMfcc.dEnergy = true;
  this->featureMfcc.fftSize = 512;
  this->featureMfcc.filterBankSize = 40;
  this->featureMfcc.frameSize = 16;
  this->featureMfcc.frameStep = 8;
  this->featureMfcc.highCut = 7000.0;
  this->featureMfcc.lowCut = 50.0;
	for (char i = 0; i < 8; i++) {
		this->featureMfcc.coeffSelection[i] = (char)(0xFF);
	}
  this->featureMfcc.MFCCorder = 20;
  this->featureMfcc.preEmphasizeFactor = 0.97;
  this->featureMfcc.window = sclib::wndHamming;
	this->featureMfcc.method = sclib::modeSClib;
	this->featureMfcc.sclib_frequencyScale = sclib::scaleMel;
	this->featureMfcc.sclib_minFilterBankFrequency = 0.0;
	this->featureMfcc.sclib_maxFilterBankFrequency = 0.0;
	this->featureMfcc.sclib_smoothing = sclib::smoothNone;
  this->featureNfr.FFTsize = 512;
  this->featureNfr.frameSize = 25;
  this->featureNfr.frameStep = 25;
  this->featureNfr.highCut = 7000.0;
  this->featureNfr.lowCut = 50.0;
  this->featureNfr.NFRthreshold = 0.3;
  this->featureNfr.preEmphasizeFactor = 0.97;
  this->featureNfr.window = sclib::wndHamming;
  this->featureSpectrum.FFTsize = 512;
  this->featureSpectrum.frameSize = 32;
  this->featureSpectrum.frameStep = 16;
  this->featureSpectrum.highCut = 7000.0;
  this->featureSpectrum.lowCut = 50.0;
  this->featureSpectrum.preEmphasizeFactor = 0.97;
  this->featureSpectrum.window = sclib::wndHamming;
	this->featureSpectrum.logarithmize = true;
	this->featureSpectrum.createPhase = false;
  this->featureSpectrumFlux.FFTsize = 512;
  this->featureSpectrumFlux.frameSize = 25;
  this->featureSpectrumFlux.frameStep = 25;
  this->featureSpectrumFlux.highCut = 7000.0;
  this->featureSpectrumFlux.lowCut = 50.0;
  this->featureSpectrumFlux.preEmphasizeFactor = 0.97;
  this->featureSpectrumFlux.window = sclib::wndHamming;
  this->featureSubBandPower.FFTsize = 512;
  this->featureSubBandPower.frameSize = 25;
  this->featureSubBandPower.frameStep = 25;
  this->featureSubBandPower.highCut = 7000.0;
  this->featureSubBandPower.lowCut = 50.0;
  this->featureSubBandPower.preEmphasizeFactor = 0.97;
  this->featureSubBandPower.window = sclib::wndHamming;
  this->featureZcr.frameSize = 20; //  40 by Bing  20 by Nan
  this->featureZcr.frameStep = 10; //by Bing
  this->featureZcr.highCut = 3400.0; //by Bing
  this->featureZcr.lowCut = 200.0; //by Bing
  this->featureZcr.preEmphasizeFactor = 0.97; //by Bing
	this->featureZcr.useChebyshev = false;
	this->featureZcr.scaleResult = true;
  this->featureSte.frameSize = 40; //by Bing
  this->featureSte.frameStep = 10; //by Bing
  this->featureSte.highCut = 3400.0; //by Bing
  this->featureSte.lowCut = 200.0; //by Bing
  this->featureSte.preEmphasizeFactor = 0.97; //by Bing
	this->featureSte.useButterworth = false;
	this->featureSte.scaleResult = true;
	this->featureLpc.lowCut = 50.0;
	this->featureLpc.highCut = 7000.0;
	this->featureLpc.frameSize = 16; //by Jun
	this->featureLpc.frameStep = 8; //by Jun
	this->featureLpc.LPCorder = 4; //by Jun
	this->featureLpc.preEmphasizeFactor = 0.97; //by Jun
	this->featureLpc.window = sclib::wndHamming;
	this->featureLpcResidual.lowCut = 50.0;
	this->featureLpcResidual.highCut = 7000.0;
	this->featureLpcResidual.frameSize = 16; //by Jun
	this->featureLpcResidual.frameStep = 8; //by Jun
	this->featureLpcResidual.order = 4; //by Jun
	this->featureLpcResidual.preEmphasizeFactor = 0.97; //by Jun
	this->featureLpc.computeGain = false;
	this->featureLsp.lowCut = 50.0;
	this->featureLsp.highCut = 7000.0;
	this->featureLsp.frameSize = 32;
	this->featureLsp.frameStep = 16;
	this->featureLsp.LPCorder = 10;
	this->featureLsp.window = sclib::wndHamming;
	this->featureLsp.preEmphasizeFactor = 0.97;
	this->featureLsp.bisections = 4;
	this->featureLsp.delta = 0.00781250;
	this->featureLsp.maxLoops = 10;
	this->featureLsp.method = sclib::modeSpeex;
	this->featureLsp.minSeparation = 0.0;
	this->featurePitch.frameSize = 10;
	this->featurePitch.frameStep = 10;
	this->featurePitch.method = sclib::modeESPS;
	this->featurePitch.sqrt = true;
	this->featurePitch.esps_cand_thresh = (float)(0.3);
	this->featurePitch.esps_lag_weight = (float)(0.3);
	this->featurePitch.esps_freq_weight = (float)(0.02);
	this->featurePitch.esps_trans_cost = (float)(0.005);
	this->featurePitch.esps_trans_amp = (float)(0.5);
	this->featurePitch.esps_trans_spec = (float)(0.5);
	this->featurePitch.esps_voice_bias = (float)(0.0);
	this->featurePitch.esps_double_cost = (float)(0.35);
	this->featurePitch.esps_min_f0 = (float)(60);
	this->featurePitch.esps_max_f0 = (float)(400);
	this->featurePitch.esps_n_cands = 20;
	this->featurePitch.esps_wind_dur = (float)(0.0075);
	this->featureFormant.lowCut = 0.0;
	this->featureFormant.highCut = 0.0;
	this->featureFormant.frameSize = 49;
	this->featureFormant.frameStep = 10;
	this->featureFormant.sampleRate = 0.0;
	this->featureFormant.preEmphasizeFactor = .7;
	this->featureFormant.esps_cor_wdur = .01;
	this->featureFormant.esps_ds_freq = 10000.0;
	this->featureFormant.esps_frame_int = .01;
	this->featureFormant.esps_lpc_ord = 12;
	this->featureFormant.esps_lpc_type = 0;
	this->featureFormant.esps_nform = 4;
	this->featureFormant.esps_nom_f1 = -10.0;
	this->featureFormant.esps_w_type = 1;
	this->featureFormant.esps_wdur = .049;
	this->featureSdp.m = 6; // by Bing 
	this->featureSdp.lag = 1; // by Bing
	this->featureSdp.color = 24; // by Bing
	this->featureSdp.frameSize =50; // by Bing
	this->featureSdp.frameStep =50; // by Bing
	this->featureSdp.n = 0; // by Bing 
	this->featureSdp.pictureSize = 250; // by Bing
	this->featureSdp.tau = 50; // by Bing
	this->featureSdp.preEmphasizeFactor = 0.97; // by Bing 
  this->featureSamples.frameSize = 32;
  this->featureSamples.frameStep = 16;
  this->featureSamples.highCut = 0.0;
  this->featureSamples.lowCut = 0.0;
	this->general.firstScene = 1;
  this->general.lastScene = 99999;
  this->general.sceneSelection = 0; //(int)pow(2, 3)|(int)pow(2, 5)|(int)pow(2, 6)|(int)pow(2, 7);
  this->general.shortSpeechThreshold = 0;
  this->general.pauseSilenceThreshold = 999999;
  MArray_1D(this->general.preClusteringResultsPrefix, 1, char, "SC_TweakableParameters: general.preClusteringResultsPrefix");
  sprintf(this->general.preClusteringResultsPrefix, "\0");
  MArray_1D(this->general.featurePrefix, 1, char, "SC_TweakableParameters: general.featurePrefix");
  sprintf(this->general.featurePrefix, "\0");
  this->groundTruth.internalFrameSize = 10;
  this->groundTruth.pseudoSceneLength = 120000;
  this->groundTruth.videoFrameMachineOffset = 0;
	this->groundTruth.storeProbabilityInformation = true;
  this->mixtureModelBgmm.EMthreshold = 100;
  this->mixtureModelBgmm.maxEMiterations = 150;
  this->mixtureModelBgmm.varianceLimit = 0.01;
  this->mixtureModelBgmm.weightLimit = 0.01;
	this->mixtureModelBgmm.fullCovariance = false;
	this->mixtureModelBgmm.randomInitialization = true;
  this->mixtureModelGmm.EMthreshold = 100;
  this->mixtureModelGmm.kMeansIterations = 10;
  this->mixtureModelGmm.maxEMiterations = 150;
  this->mixtureModelGmm.varianceLimit = 0.01;
  this->mixtureModelGmm.weightLimit = 0.01;
  this->mixtureModelMixMax.bgModelCombination = true;
  this->mixtureModelMixMax.EMthreshold = 100;
  this->mixtureModelMixMax.kMeansIterations = 10;
  this->mixtureModelMixMax.maxEMiterations = 150;
  this->mixtureModelMixMax.maxERFloops = 100;
  this->mixtureModelMixMax.noiseCorruptionType = sclib::nctMax;
  this->mixtureModelMixMax.varianceLimit = 0.01;
  this->mixtureModelMixMax.weightLimit = 0.01;
  this->mixtureModelGmmubm.adaptMeans = true;
  this->mixtureModelGmmubm.adaptVariances = false;
  this->mixtureModelGmmubm.adaptWeights = false;
  this->mixtureModelGmmubm.relevanceFactor = 16.0;
  this->mixtureModelGmmubm.scoringMethod = sclib::scoringGMM_UBM;
  this->mixtureModelGmmubm.topCmixtures = 5;
  MArray_1D(this->mixtureModelGmmubm.ubmFileName, 1, char, "SC_TweakableParameters: mixtureModelGmmubm.ubmFileName");
  sprintf(this->mixtureModelGmmubm.ubmFileName, "\0");
  this->mixtureModelGmmubm.varianceLimit = 0.0;
  this->mixtureModelMix2Max.EMthreshold = 100;
  this->mixtureModelMix2Max.kMeansIterations = 10;
  this->mixtureModelMix2Max.maxEMiterations = 150;
  this->mixtureModelMix2Max.maxERFloops = 100;
  this->mixtureModelMix2Max.varianceLimit = 0.01;
  this->mixtureModelMix2Max.weightLimit = 0.01;
  this->mixtureModelMix2Max.bgModelCombination = false;
  this->mixtureModelMix2MaxEx.bgModelCombination = true;
  this->mixtureModelMix2MaxEx.EMthreshold = 100;
  this->mixtureModelMix2MaxEx.kMeansIterations = 10;
  this->mixtureModelMix2MaxEx.maxEMiterations = 150;
  this->mixtureModelMix2MaxEx.maxERFloops = 100;
  this->mixtureModelMix2MaxEx.varianceLimit = 0.01;
  this->mixtureModelMix2MaxEx.weightLimit = 0.01;
  this->modelHandler.backgroundModelType = sclib::mtGMM_new;
  this->modelHandler.foregroundModelType = sclib::mtGMM_new;
  this->modelHandler.maxNoiseModelOrder = 8;
  this->modelHandler.maxSpeakerModelOrder = 16;
  this->modelHandler.msPerGaussian = 1000;
  this->modelHandler.orderGuessEMsteps = 3;
  this->modelHandler.orderGuessMode = sclib::guessHeuristic;
  this->modelHandler.SNRthreshold = 15.0;
	this->modelHandler.outlierRemovalMode = sclib::outlierRemoveNone;
	this->modelHmm.stateCount = 4;
  MArray_1D(this->modelHandler.onlyThisSpeaker, 1, char, "SC_TweakableParameters: modelHandler.onlyThisSpeaker");
  sprintf(this->modelHandler.onlyThisSpeaker, "\0");
	this->modelHandler.speakerModelFeature = sclib::featureMFCC;
	MArray_1D(this->modelHmm.transitionStructure, sclib::bufferSize, char, "SC_TweakableParameters: modelHmm.transitionStructure");
	sprintf(this->modelHmm.transitionStructure, "1 1 1 1;1 1 1 1;1 1 1 1;1 1 1 1");
	this->modelHmm.mixturesPerState = 8;
	this->modelHmm.useOrthogonalTransform = false;
	this->modelHmm.leftToRight = false;
	this->modelHmm.maxIterations = 100;
	this->modelHmm.verbose = false;
	this->modelVq.codebookSize = 8;
	this->modelVq.splitMethod = sclib::modeLBG;
	this->modelVq.maxIterations = 100;
	this->modelSvm.distanceBasedTesting = true;
	this->modelSvm.doParameterSearch = false;
	this->modelTime.subModelType = sclib::mtSVM;
	this->modelTime.syllableLength = 100;
	this->modelTime.trajectoryStep = 1;
	this->modelTime.removeTiming = true;
	this->modelTime.templateCount = 0;
	this->modelTime.clusteringIterations = 100;
	this->modelTime.replaceTrainingData = true;
	this->modelTime.checkForTrajectorization = true;
  MArray_1D(this->modelTime.worldModelFile, 1, char, "SC_TweakableParameters: modelTime.worldModelFile");
  sprintf(this->modelTime.worldModelFile, "\0");
  MArray_1D(this->modelTime.normalizationFile, 1, char, "SC_TweakableParameters: modelTime.normalizationFile");
  sprintf(this->modelTime.normalizationFile, "\0");
	this->modelPareto.useMarginalDistributions = true;
  this->modelQgmm.deltaBIClambda = 1.0;
  this->modelQgmm.maxMixtures = 32;
  this->modelQgmm.percentDifference = 10.0;
  this->resampling.fastConversion = false;
  this->score.BBNmetricLambda = 0.5;
  MArray_1D(this->segmentationAudioTypeLzl.actionModelFileName, 1, char, "SC_TweakableParameters: segmentationAudioTypeLzl.actionModelFileName");
  sprintf(this->segmentationAudioTypeLzl.actionModelFileName, "\0");
  MArray_1D(this->segmentationAudioTypeLzl.classifierFileName, 1, char, "SC_TweakableParameters: segmentationAudioTypeLzl.classifierFileName");
  sprintf(this->segmentationAudioTypeLzl.classifierFileName, "\0");
  MArray_1D(this->segmentationAudioTypeLzl.featureFileName, 1, char, "SC_TweakableParameters: segmentationAudioTypeLzl.featureFileName");
  sprintf(this->segmentationAudioTypeLzl.featureFileName, "\0");
  MArray_1D(this->segmentationAudioTypeLzl.normalizationFileName, 1, char, "SC_TweakableParameters: segmentationAudioTypeLzl.normalizationFileName");
  sprintf(this->segmentationAudioTypeLzl.normalizationFileName, "\0");
	this->segmentationAudioTypeLzl.speechSpecificity = 0.5;
	this->segmentationAudioTypeLzl.musicSpecificity = 0.5;
	this->segmentationAudioTypeLzl.subClipLength = 1000;
	this->segmentationAudioTypeLzl.bbParameters.client = sclib::algorithm_ac_LZL;
	this->segmentationAudioTypeLzl.bbParameters.FFTsize = 512;
	this->segmentationAudioTypeLzl.bbParameters.frameSize = 25;
	this->segmentationAudioTypeLzl.bbParameters.frameStep = 25;
	this->segmentationAudioTypeLzl.bbParameters.preEmphasizeFactor = 0.97;
	this->segmentationAudioTypeLzl.bbParameters.window = sclib::wndHamming;
	this->segmentationAudioTypeLzl.bbParameters.sampleRate = 8000.0;
	this->segmentationAudioTypeLzl.bpParameters.client = sclib::algorithm_ac_LZL;
	this->segmentationAudioTypeLzl.bpParameters.FFTsize = 512;
	this->segmentationAudioTypeLzl.bpParameters.frameSize = 25;
	this->segmentationAudioTypeLzl.bpParameters.frameStep = 25;
	this->segmentationAudioTypeLzl.bpParameters.preEmphasizeFactor = 0.97;
	this->segmentationAudioTypeLzl.bpParameters.window = sclib::wndHamming;
	this->segmentationAudioTypeLzl.bpParameters.sampleRate = 8000.0;
	this->segmentationAudioTypeLzl.mfccParameters.client = sclib::algorithm_ac_LZL;
	this->segmentationAudioTypeLzl.mfccParameters.frameSize = 25;
	this->segmentationAudioTypeLzl.mfccParameters.frameStep = 25;
	this->segmentationAudioTypeLzl.mfccParameters.MFCCorder = 8;
	this->segmentationAudioTypeLzl.mfccParameters.addDeltaDeltas = false;
	this->segmentationAudioTypeLzl.mfccParameters.addDeltas = false;
	this->segmentationAudioTypeLzl.mfccParameters.method = sclib::modeSClib;
	this->segmentationAudioTypeLzl.mfccParameters.dEnergy = true;
	this->segmentationAudioTypeLzl.mfccParameters.CMN = false;
	this->segmentationAudioTypeLzl.mfccParameters.fftSize = 512;
	this->segmentationAudioTypeLzl.mfccParameters.filterBankSize = 24;
	this->segmentationAudioTypeLzl.mfccParameters.preEmphasizeFactor = 0.97;
	this->segmentationAudioTypeLzl.mfccParameters.sclib_frequencyScale = sclib::scaleExpoLog;
	this->segmentationAudioTypeLzl.mfccParameters.sclib_smoothing = sclib::smoothNone;
	this->segmentationAudioTypeLzl.mfccParameters.window = sclib::wndHamming;
	for (char i = 0; i < 8; i++) {
		this->segmentationAudioTypeLzl.mfccParameters.coeffSelection[i] = (char)(0xFF);
	}
	this->segmentationAudioTypeLzl.mfccParameters.sclib_maxFilterBankFrequency = 0.0;
	this->segmentationAudioTypeLzl.mfccParameters.sclib_minFilterBankFrequency = 0.0;
	this->segmentationAudioTypeLzl.mfccParameters.sampleRate = 8000.0;
	this->segmentationAudioTypeLzl.nfrParameters.client = sclib::algorithm_ac_LZL;
	this->segmentationAudioTypeLzl.nfrParameters.FFTsize = 512;
	this->segmentationAudioTypeLzl.nfrParameters.frameSize = 25;
	this->segmentationAudioTypeLzl.nfrParameters.frameStep = 25;
	this->segmentationAudioTypeLzl.nfrParameters.NFRthreshold = 0.3;
	this->segmentationAudioTypeLzl.nfrParameters.preEmphasizeFactor = 0.97;
	this->segmentationAudioTypeLzl.nfrParameters.window = sclib::wndHamming;
	this->segmentationAudioTypeLzl.nfrParameters.sampleRate = 8000.0;
	this->segmentationAudioTypeLzl.sbpParameters.client = sclib::algorithm_ac_LZL;
	this->segmentationAudioTypeLzl.sbpParameters.FFTsize = 512;
	this->segmentationAudioTypeLzl.sbpParameters.frameSize = 25;
	this->segmentationAudioTypeLzl.sbpParameters.frameStep = 25;
	this->segmentationAudioTypeLzl.sbpParameters.preEmphasizeFactor = 0.97;
	this->segmentationAudioTypeLzl.sbpParameters.window = sclib::wndHamming;
	this->segmentationAudioTypeLzl.sbpParameters.sampleRate = 8000.0;
	this->segmentationAudioTypeLzl.sfParameters.client = sclib::algorithm_ac_LZL;
	this->segmentationAudioTypeLzl.sfParameters.FFTsize = 512;
	this->segmentationAudioTypeLzl.sfParameters.frameSize = 25;
	this->segmentationAudioTypeLzl.sfParameters.frameStep = 25;
	this->segmentationAudioTypeLzl.sfParameters.preEmphasizeFactor = 0.97;
	this->segmentationAudioTypeLzl.sfParameters.window = sclib::wndHamming;
	this->segmentationAudioTypeLzl.sfParameters.sampleRate = 8000.0;
	this->segmentationAudioTypeLzl.zcrParameters.client = sclib::algorithm_ac_LZL;
	this->segmentationAudioTypeLzl.zcrParameters.frameSize = 25;
	this->segmentationAudioTypeLzl.zcrParameters.frameStep = 25;
	this->segmentationAudioTypeLzl.zcrParameters.preEmphasizeFactor = 0.97;
	this->segmentationAudioTypeLzl.zcrParameters.scaleResult = true;
	this->segmentationAudioTypeLzl.zcrParameters.useChebyshev = false;
	this->segmentationAudioTypeLzl.zcrParameters.sampleRate = 8000.0;
  this->segmentationAudioTypeLzl.svmParameters.C = 0; //i.e. use data-dependant default
  this->segmentationAudioTypeLzl.svmParameters.cache_size = 100;
  this->segmentationAudioTypeLzl.svmParameters.coef0 = 0;
  this->segmentationAudioTypeLzl.svmParameters.cvFolds = 5;
  this->segmentationAudioTypeLzl.svmParameters.cvMaxDatasetSize = 15000;
  this->segmentationAudioTypeLzl.svmParameters.cvCoarseCMin = -5.0;
  this->segmentationAudioTypeLzl.svmParameters.cvCoarseCMax = 15.0;
  this->segmentationAudioTypeLzl.svmParameters.cvCoarseCStep = 2.0;
  this->segmentationAudioTypeLzl.svmParameters.cvCoarseGammaMin = -15.0;
  this->segmentationAudioTypeLzl.svmParameters.cvCoarseGammaMax = 3.0;
  this->segmentationAudioTypeLzl.svmParameters.cvCoarseGammaStep = 2.0;
  this->segmentationAudioTypeLzl.svmParameters.cvFineCStep = 0.25;
  this->segmentationAudioTypeLzl.svmParameters.cvFineCRadius = 2.0;
  this->segmentationAudioTypeLzl.svmParameters.cvFineGammaStep = 0.25;
  this->segmentationAudioTypeLzl.svmParameters.cvFineGammaRadius = 2.0;
  this->segmentationAudioTypeLzl.svmParameters.degree = 3;
  this->segmentationAudioTypeLzl.svmParameters.doCV = true;
  this->segmentationAudioTypeLzl.svmParameters.eps = 1e-3;
  this->segmentationAudioTypeLzl.svmParameters.gamma = 0; //i.e. use data-dependant default
  this->segmentationAudioTypeLzl.svmParameters.kernel_type = 2;
  this->segmentationAudioTypeLzl.svmParameters.nu = 0.5;
  this->segmentationAudioTypeLzl.svmParameters.p = 0.1;
  this->segmentationAudioTypeLzl.svmParameters.probability = 1;
  this->segmentationAudioTypeLzl.svmParameters.shrinking = 1;
  this->segmentationAudioTypeLzl.svmParameters.svm_type = 0;
	this->segmentationAudioTypeLzl.svmParameters.oneClassGammaSearchMaxIterations = 1000000;
	this->segmentationAudioTypeLzl.svmParameters.oneClassGammaSearchRepeats = 10;
  this->segmentationChangesLz.adaptiveThresholdAlpha = 1.2;
	this->segmentationChangesLz.bayesianThreshold = 2.0;
  this->segmentationChangesLz.detectorWindowLength = 3000;
  this->segmentationChangesLz.detectorWindowStep = 500;
  this->segmentationChangesLz.lastNdistances = 4;
  MArray_1D(this->segmentationChangesLz.priorsFileName, 1, char, "SC_TweakableParameters: segmentationChangesLz.priorsFileName");
  sprintf(this->segmentationChangesLz.priorsFileName, "\0");
  MArray_1D(this->segmentationChangesLz.time2ChangeModelFileName, 2, char, "SC_TweakableParameters: segmentationChangesLz.time2ChangeModelFileName");
  sprintf(this->segmentationChangesLz.time2ChangeModelFileName, "\0");
	this->segmentationChangesKbk.r = 3.286;
	this->segmentationChangesKbk.lambda = 1.0;
	this->segmentationChangesKbk.tolerance = 1.0;
	this->segmentationChangesKbk.mfccParameters.client = sclib::algorithm_cd_KBK;
	this->segmentationChangesKbk.mfccParameters.frameSize = 16;
	this->segmentationChangesKbk.mfccParameters.frameStep = 10;
	this->segmentationChangesKbk.mfccParameters.MFCCorder = 36;
	this->segmentationChangesKbk.mfccParameters.addDeltaDeltas = true;
	this->segmentationChangesKbk.mfccParameters.addDeltas = true;
	this->segmentationChangesKbk.mfccParameters.method = sclib::modeSClib;
	this->segmentationChangesKbk.mfccParameters.dEnergy = false;
	this->segmentationChangesKbk.mfccParameters.CMN = false;
	this->segmentationChangesKbk.mfccParameters.fftSize = 512;
	this->segmentationChangesKbk.mfccParameters.filterBankSize = 40;
	this->segmentationChangesKbk.mfccParameters.preEmphasizeFactor = 0.97;
	this->segmentationChangesKbk.mfccParameters.sclib_frequencyScale = sclib::scaleMel;
	this->segmentationChangesKbk.mfccParameters.sclib_smoothing = sclib::smoothNone;
	this->segmentationChangesKbk.mfccParameters.window = sclib::wndHamming;
	this->segmentationChangesKbk.mfccParameters.coeffSelection[0] = 0xFD; //overall bit 1, 3, 4, 5, 6, 7, 8
	this->segmentationChangesKbk.mfccParameters.coeffSelection[1] = 0x97; //overall bit 9, 10, 11, 13, 16
	this->segmentationChangesKbk.mfccParameters.coeffSelection[2] = 0xE0; //overall bit 22, 23, 24
	this->segmentationChangesKbk.mfccParameters.coeffSelection[3] = 0x5F; //overall bit 25, 26, 27, 28, 29, 31
	this->segmentationChangesKbk.mfccParameters.coeffSelection[4] = 0x0D; //overall bit 33, 35, 36
	for (char i = 5; i < 8; i++) {
		this->segmentationChangesKbk.mfccParameters.coeffSelection[i] = (char)(0x00);
	}
	this->segmentationChangesKbk.mfccParameters.sclib_maxFilterBankFrequency = 0.0;
	this->segmentationChangesKbk.mfccParameters.sclib_minFilterBankFrequency = 0.0;
	this->segmentationHandler.audioTypeMode = sclib::algorithm_ac_LZL;
  this->segmentationHandler.changeDetectorMode = sclib::algorithm_cd_LZW;
  this->segmentationHandler.silenceDetectorMode = sclib::algorithm_sd_LNK;
  this->segmentationHandler.vUvDetectorMode = sclib::algorithm_vud_ESPS;
  this->segmentationSilenceLnk.energyQuantizationLevel = 10;
  this->segmentationSilenceLzl.zcrSilenceThreshold = 1.710257; //threshold was created by taking mean+2sd of zcr on /TRECVID2006/silence.wav
  this->segmentationSilenceLzl.energySilenceThreshold = 37.442281; //threshold was created by taking mean+2sd of ste on /TRECVID2006/silence.wav
	this->segmentationSilenceLzl.specificity = 0.5;
	this->segmentationSilenceLzl.zcrParameters.client = sclib::algorithm_sd_LZL;
  this->segmentationSilenceLzl.zcrParameters.frameSize = 25;
  this->segmentationSilenceLzl.zcrParameters.frameStep = 25;
	this->segmentationSilenceLzl.zcrParameters.sampleRate = 4000.0;
  this->segmentationSilenceLzl.zcrParameters.preEmphasizeFactor = 0.97;
	this->segmentationSilenceLzl.zcrParameters.useChebyshev = false;
	this->segmentationSilenceLzl.zcrParameters.scaleResult = true;
	this->segmentationSilenceLzl.steParameters.client = sclib::algorithm_sd_LZL;
  this->segmentationSilenceLzl.steParameters.frameSize = 25;
  this->segmentationSilenceLzl.steParameters.frameStep = 25;
  this->segmentationSilenceLzl.steParameters.sampleRate = 4000.0;
	this->segmentationSilenceLzl.steParameters.preEmphasizeFactor = 0.97;
	this->segmentationSilenceLzl.steParameters.useButterworth = false;
	this->segmentationSilenceLzl.steParameters.scaleResult = true;
	this->segmentationVUvEsps.pitchParameters.client = sclib::algorithm_vud_ESPS;
	this->segmentationVUvEsps.pitchParameters.frameSize = 10;
	this->segmentationVUvEsps.pitchParameters.frameStep = 10;
	this->segmentationVUvEsps.pitchParameters.method = sclib::modeESPS;
	this->segmentationVUvEsps.pitchParameters.esps_cand_thresh = (float)(0.3);
	this->segmentationVUvEsps.pitchParameters.esps_lag_weight = (float)(0.3);
	this->segmentationVUvEsps.pitchParameters.esps_freq_weight = (float)(0.02);
	this->segmentationVUvEsps.pitchParameters.esps_trans_cost = (float)(0.005);
	this->segmentationVUvEsps.pitchParameters.esps_trans_amp = (float)(0.5);
	this->segmentationVUvEsps.pitchParameters.esps_trans_spec = (float)(0.5);
	this->segmentationVUvEsps.pitchParameters.esps_voice_bias = (float)(0.0);
	this->segmentationVUvEsps.pitchParameters.esps_double_cost = (float)(0.35);
	this->segmentationVUvEsps.pitchParameters.esps_min_f0 = (float)(60);
	this->segmentationVUvEsps.pitchParameters.esps_max_f0 = (float)(400);
	this->segmentationVUvEsps.pitchParameters.esps_n_cands = 20;
	this->segmentationVUvEsps.pitchParameters.esps_wind_dur = (float)(0.0075);
  this->segmentationVUvLnk.energyQuantizationLevel = 10;
	this->signalHandler.forceSampleRate = 0;
	this->signalMpeg.fastSeeking = false;
	this->signalMpeg.hqResampling = true;
	this->signalMpeg.outputChannelCount = 1;
	this->signalMpeg.outputSampleRate = 16000;
	this->speakerClusterer.doClustering = true;
	this->speakerClusterer.constructNonOverlappingClusters = false;
  this->speakerClusterer.speechSegLengthThreshold = 0;
  this->speakerClusterer.distanceMeasure = sclib::dmGLR;
  this->speakerClusterer.globalCriterion = sclib::gcBIC;
  this->speakerClusterer.terminationCriterion = sclib::tcGc;
	this->speakerClusterer.linkageMode = sclib::linkageMerge;
	this->speakerClusterer.specificity = 0.5;
  MArray_1D(this->speakerClusterer.firstDistanceMatrixPrefix, 1, char, "SC_TweakableParameters: speakerClusterer.firstDistanceMatrixPrefix");
  sprintf(this->speakerClusterer.firstDistanceMatrixPrefix, "\0");
	this->speakerIdentification.dNormSampleCount = 2500; //approx. 30 sec.
	this->speakerIdentification.doIdentification = false;
  this->speakerIdentification.normalizationMode = sclib::normalizationDNorm;
	this->speakerIdentification.useUBMs = false;
  this->transform.taperingLength = 1.0;
	this->clusterer.maxIterations = 1; // by Bing
	this->clusterer.numCluster = 64; //by Bing
	this->clusterer.iterationsToIniKMeanList = 3; //by Bing

  //then, try to load user-defined values from the ini-file, if given:
  if (iniFileName != NULL && strncmp(iniFileName, "", 2) != 0) {
    pIni = new SC_Ini();

    if (pIni->openIni(iniFileName) == true) {
      while (pIni->readNextParameter(key, value) == true) {
        if (setByName(key, value) == false) {
          printf("\n Error in ini-file: Key-value-pair '%s=%s' is not valid", key, value);
        }
        MFree_1D(key);
        MFree_1D(value);
      }
      pIni->closeIni();
    } else {
      REPORT_ERROR(SVLIB_FileErr, "Given ini-file not found");
    }

    MFree_0D(pIni);
  }

  this->verbose = verbose;
}

//====================================================================================================================
//	The destructor
//====================================================================================================================
SC_TweakableParameters::~SC_TweakableParameters() {
	this->parameters.map.clear();
	MFree_1D(this->debug.debugDir);
	MFree_1D(this->enhancement.speechModelFile);
	MFree_1D(this->mixtureModelGmmubm.ubmFileName);
	MFree_1D(this->modelHandler.onlyThisSpeaker);
  MFree_1D(this->originalDebugDir);
	MFree_1D(this->general.preClusteringResultsPrefix);
	MFree_1D(this->general.featurePrefix);
  MFree_1D(this->segmentationAudioTypeLzl.actionModelFileName);
  MFree_1D(this->segmentationAudioTypeLzl.normalizationFileName);
  MFree_1D(this->segmentationAudioTypeLzl.featureFileName);
  MFree_1D(this->segmentationAudioTypeLzl.classifierFileName);
  MFree_1D(this->segmentationChangesLz.priorsFileName);
  MFree_1D(this->segmentationChangesLz.time2ChangeModelFileName);
	MFree_1D(this->modelHmm.transitionStructure);
	MFree_1D(this->modelTime.worldModelFile);
	MFree_1D(this->modelTime.normalizationFile);
	MFree_1D(this->speakerClusterer.firstDistanceMatrixPrefix);
}

//====================================================================================================================
//	An auxiliary function for setting the debugDir-parameter
//====================================================================================================================
void SC_TweakableParameters::setDebugDir(const char* dir) {
  MFree_1D(this->debug.debugDir);
  int len;
  
  if (((strrchr(dir, '/') - dir) / sizeof(char)) < strlen(dir)) {
    MArray_1D(this->debug.debugDir, strlen(dir)+2, char, "SC_TweakableParameters.setDebugDir: debugDir");
		sprintf(this->debug.debugDir, "%s/", dir);
	} else {
    MArray_1D(this->debug.debugDir, strlen(dir)+1, char, "SC_TweakableParameters.setDebugDir: debugDir");
		sprintf(this->debug.debugDir, "%s", dir);
	}

  //store a copy of the debug dir to make resetting the prefix by unsetDebugPrefix() possible
  MFree_1D(this->originalDebugDir);
  len = (int)(strlen(this->debug.debugDir));
  MArray_1D(this->originalDebugDir, len+1, char, "SC_TweakableParameters.setDebugDir: originalDebugDir");
  strncpy(this->originalDebugDir, this->debug.debugDir, len);
  this->originalDebugDir[len] = '\0';

	return;
}

//====================================================================================================================
//	Puts the prefix after the debugDir (that means: in front of each debug output)
//====================================================================================================================
void SC_TweakableParameters::setDebugPrefix(const char* prefix) {
  char *newDebugDir = NULL, *newPrefix = NULL;
  int len;

  //just make a copy of the given prefix and exchange ' ' with '_'
  len = (int)(strlen(prefix));
  MArray_1D(newPrefix, len+1, char, "SC_TweakableParameters.setDebugPrefix: newPrefix");
  newPrefix = strncpy(newPrefix, prefix, len);
  newPrefix[len] = '\0';
  sclib::strReplace(newPrefix, ' ', '_');

  if (this->debug.useDebugPrefix == true) {
    len = (int)(strlen(this->originalDebugDir) + strlen(newPrefix) + 2);
    //len = strlen(this->debugDir) + strlen(newPrefix) + 2;
    MArray_1D(newDebugDir, len, char, "SC_TweakableParameters.setDebugPrefix: newDebugDir");
    sprintf(newDebugDir, "%s%s_", this->originalDebugDir, newPrefix);
    //sprintf(newDebugDir, "%s%s_\0", this->debugDir, newPrefix);

    MFree_1D(this->debug.debugDir);
    this->debug.debugDir = newDebugDir;
  }

  MFree_1D(newPrefix);

  return;
}

//====================================================================================================================
//	Removes the previously setted prefix from debugDir by restoring the value in originalDebugDir
//====================================================================================================================
void SC_TweakableParameters::unsetDebugPrefix(void) {
  MFree_1D(this->debug.debugDir);
  MArray_1D(this->debug.debugDir, strlen(this->originalDebugDir)+1, char, "SC_TweakableParameters.unsetDebugPrefix: debugDir");
  sprintf(this->debug.debugDir, "%s", this->originalDebugDir);

  return;
}

//====================================================================================================================
//	A function for printing out the parameters
//====================================================================================================================
std::ostream& operator<< (std::ostream& OutS, SC_TweakableParameters& pTweak) {
  if (pTweak.verbose == true) {
		OutS << "classifierAdaBoost.maxWeakClassifiers=" << pTweak.classifierAdaBoost.maxWeakClassifiers << endl;
		OutS << "classifierAdaBoost.weakClassifierType=" << pTweak.classifierAdaBoost.weakClassifierType << endl;

		OutS << "classifierDecisionStump.splitMultiway=" << pTweak.classifierDecisionStump.splitMultiway << endl;

    OutS << "classifierMl.maxOrder=" << pTweak.classifierMl.maxOrder << endl;
    OutS << "classifierMl.minOrder=" << pTweak.classifierMl.minOrder << endl;
    OutS << "classifierMl.modelType=" << pTweak.classifierMl.modelType << endl;

    OutS << "classifierSvm.C=" << pTweak.classifierSvm.C << endl;
    OutS << "classifierSvm.cache_size=" << pTweak.classifierSvm.cache_size << endl;
    OutS << "classifierSvm.coef0=" << pTweak.classifierSvm.coef0 << endl;
    OutS << "classifierSvm.cvFolds=" << pTweak.classifierSvm.cvFolds << endl;
    OutS << "classifierSvm.cvMaxDatasetSize=" << pTweak.classifierSvm.cvMaxDatasetSize << endl;
    OutS << "classifierSvm.cvCoarseCMin=" << pTweak.classifierSvm.cvCoarseCMin << endl;
    OutS << "classifierSvm.cvCoarseCMax=" << pTweak.classifierSvm.cvCoarseCMax << endl;
    OutS << "classifierSvm.cvCoarseCStep=" << pTweak.classifierSvm.cvCoarseCStep << endl;
    OutS << "classifierSvm.cvCoarseGammaMin=" << pTweak.classifierSvm.cvCoarseGammaMin << endl;
    OutS << "classifierSvm.cvCoarseGammaMax=" << pTweak.classifierSvm.cvCoarseGammaMax << endl;
    OutS << "classifierSvm.cvCoarseGammaStep=" << pTweak.classifierSvm.cvCoarseGammaStep << endl;
    OutS << "classifierSvm.cvFineCStep=" << pTweak.classifierSvm.cvFineCStep << endl;
    OutS << "classifierSvm.cvFineCRadius=" << pTweak.classifierSvm.cvFineCRadius << endl;
    OutS << "classifierSvm.cvFineGammaStep=" << pTweak.classifierSvm.cvFineGammaStep << endl;
    OutS << "classifierSvm.cvFineGammaRadius=" << pTweak.classifierSvm.cvFineGammaRadius << endl;
    OutS << "classifierSvm.degree=" << pTweak.classifierSvm.degree << endl;
    OutS << "classifierSvm.doCV=" << pTweak.classifierSvm.doCV << endl;
    OutS << "classifierSvm.eps=" << pTweak.classifierSvm.eps << endl;
    OutS << "classifierSvm.gamma=" << pTweak.classifierSvm.gamma << endl;
    OutS << "classifierSvm.kernel_type=" << pTweak.classifierSvm.kernel_type << endl;
    OutS << "classifierSvm.nu=" << pTweak.classifierSvm.nu << endl;
    OutS << "classifierSvm.p=" << pTweak.classifierSvm.p << endl;
    OutS << "classifierSvm.probability=" << pTweak.classifierSvm.probability << endl;
    OutS << "classifierSvm.shrinking=" << pTweak.classifierSvm.shrinking << endl;
    OutS << "classifierSvm.svm_type=" << pTweak.classifierSvm.svm_type << endl;
		OutS << "classifierSvm.oneClassGammaSearchMaxIterations=" << pTweak.classifierSvm.oneClassGammaSearchMaxIterations << endl;
		OutS << "classifierSvm.oneClassGammaSearchRepeats=" << pTweak.classifierSvm.oneClassGammaSearchRepeats << endl;

    OutS << "cluster.mergeMode=" << pTweak.cluster.mergeMode << endl;

		OutS << "clusterer.maxIterations=" << pTweak.clusterer.maxIterations  << endl; // by bing
		OutS << "clusterer.numCluster=" <<pTweak.clusterer.numCluster<< endl; // by bing
		OutS << "clusterer.iterationsToIniKMeanList=" <<pTweak.clusterer.iterationsToIniKMeanList<< endl; // by bing
    
		OutS << "debug.debugDir=" << sclib::nullFilter(pTweak.debug.debugDir) << endl;
    OutS << "debug.debugMode=" << pTweak.debug.debugMode << endl;
    OutS << "debug.useDebugPrefix=" << pTweak.debug.useDebugPrefix << endl;
    
    OutS << "distanceMeasure.BICpenaltyFactor=" << pTweak.distanceMeasure.BICpenaltyFactor << endl;
    OutS << "distanceMeasure.mergeMode=" << pTweak.distanceMeasure.mergeMode << endl;
    OutS << "distanceMeasure.WCDpenaltyFactor=" << pTweak.distanceMeasure.WCDpenaltyFactor << endl;
    OutS << "distanceMeasure.groundDistance=" << pTweak.distanceMeasure.groundDistance << endl;
		OutS << "distanceMeasure.ICRthreshold=" << pTweak.distanceMeasure.ICRthreshold << endl;
		
		OutS << "emd.debugLevel=" << pTweak.emd.debugLevel << endl;
		OutS << "emd.maxIterations=" << pTweak.emd.maxIterations << endl;
		OutS << "emd.maxSigSize=" << pTweak.emd.maxSigSize << endl;
		
    OutS << "enhancement.doEnhancement=" << pTweak.enhancement.doEnhancement << endl;
    OutS << "enhancement.minNoiseDuration=" << pTweak.enhancement.minNoiseDuration << endl;
    OutS << "enhancement.noiseModelUpdateRate=" << pTweak.enhancement.noiseModelUpdateRate << endl;
    OutS << "enhancement.speechModelFile=" << sclib::nullFilter(pTweak.enhancement.speechModelFile) << endl;
    
    OutS << "featureBandPeriodicity.FFTsize=" << pTweak.featureBandPeriodicity.FFTsize << endl;
    OutS << "featureBandPeriodicity.frameSize=" << pTweak.featureBandPeriodicity.frameSize << endl;
    OutS << "featureBandPeriodicity.frameStep=" << pTweak.featureBandPeriodicity.frameStep << endl;
    OutS << "featureBandPeriodicity.highCut=" << pTweak.featureBandPeriodicity.highCut << endl;
    OutS << "featureBandPeriodicity.lowCut=" << pTweak.featureBandPeriodicity.lowCut << endl;
    OutS << "featureBandPeriodicity.preEmphasizeFactor=" << pTweak.featureBandPeriodicity.preEmphasizeFactor << endl;
    OutS << "featureBandPeriodicity.window=" << pTweak.featureBandPeriodicity.window << endl;
    
    OutS << "featureBrightnessBandwidth.FFTsize=" << pTweak.featureBrightnessBandwidth.FFTsize << endl;
    OutS << "featureBrightnessBandwidth.frameSize=" << pTweak.featureBrightnessBandwidth.frameSize << endl;
    OutS << "featureBrightnessBandwidth.frameStep=" << pTweak.featureBrightnessBandwidth.frameStep << endl;
    OutS << "featureBrightnessBandwidth.highCut=" << pTweak.featureBrightnessBandwidth.highCut << endl;
    OutS << "featureBrightnessBandwidth.lowCut=" << pTweak.featureBrightnessBandwidth.lowCut << endl;
    OutS << "featureBrightnessBandwidth.preEmphasizeFactor=" << pTweak.featureBrightnessBandwidth.preEmphasizeFactor << endl;
    OutS << "featureBrightnessBandwidth.window=" << pTweak.featureBrightnessBandwidth.window << endl;
   
    OutS << "featureFbe.addDeltaDeltas=" << pTweak.featureFbe.addDeltaDeltas << endl;
    OutS << "featureFbe.addDeltas=" << pTweak.featureFbe.addDeltas << endl;
    OutS << "featureFbe.CMN=" << pTweak.featureFbe.CMN << endl;
    OutS << "featureFbe.dEnergy=" << pTweak.featureFbe.dEnergy << endl;
    OutS << "featureFbe.FFTsize=" << pTweak.featureFbe.FFTsize << endl;
    OutS << "featureFbe.filterBankSize=" << pTweak.featureFbe.filterBankSize << endl;
    OutS << "featureFbe.frameSize=" << pTweak.featureFbe.frameSize << endl;
    OutS << "featureFbe.frameStep=" << pTweak.featureFbe.frameStep << endl;
    OutS << "featureFbe.frequencyScale=" << (int)(pTweak.featureFbe.frequencyScale) << endl;
    OutS << "featureFbe.highCut=" << pTweak.featureFbe.highCut << endl;
    OutS << "featureFbe.lowCut=" << pTweak.featureFbe.lowCut << endl;
    OutS << "featureFbe.MFCCcoeffSelection=" << "'" << pTweak.featureFbe.MFCCcoeffSelection << "'" << endl;
    OutS << "featureFbe.MFCCorder=" << pTweak.featureFbe.MFCCorder << endl;
    OutS << "featureFbe.preEmphasizeFactor=" << pTweak.featureFbe.preEmphasizeFactor << endl;
    OutS << "featureFbe.smoothing=" << (int)(pTweak.featureFbe.smoothing) << endl;
    OutS << "featureFbe.window=" << pTweak.featureFbe.window << endl;
    OutS << "featureFbe.minFilterBankFrequency=" << pTweak.featureFbe.minFilterBankFrequency << endl;
    OutS << "featureFbe.maxFilterBankFrequency=" << pTweak.featureFbe.maxFilterBankFrequency << endl;
    OutS << "featureFbe.resultType=" << (int)(pTweak.featureFbe.resultType) << endl;
    
    OutS << "featureMfcc.addDeltaDeltas=" << pTweak.featureMfcc.addDeltaDeltas << endl;
    OutS << "featureMfcc.addDeltas=" << pTweak.featureMfcc.addDeltas << endl;
    OutS << "featureMfcc.CMN=" << pTweak.featureMfcc.CMN << endl;
    OutS << "featureMfcc.dEnergy=" << pTweak.featureMfcc.dEnergy << endl;
    OutS << "featureMfcc.fftSize=" << pTweak.featureMfcc.fftSize << endl;
    OutS << "featureMfcc.filterBankSize=" << pTweak.featureMfcc.filterBankSize << endl;
    OutS << "featureMfcc.frameSize=" << pTweak.featureMfcc.frameSize << endl;
    OutS << "featureMfcc.frameStep=" << pTweak.featureMfcc.frameStep << endl;
    OutS << "featureMfcc.highCut=" << pTweak.featureMfcc.highCut << endl;
    OutS << "featureMfcc.lowCut=" << pTweak.featureMfcc.lowCut << endl;
    OutS << "featureMfcc.coeffSelection=" << "'" << pTweak.featureMfcc.coeffSelection << "'" << endl;
    OutS << "featureMfcc.MFCCorder=" << pTweak.featureMfcc.MFCCorder << endl;
    OutS << "featureMfcc.preEmphasizeFactor=" << pTweak.featureMfcc.preEmphasizeFactor << endl;
    OutS << "featureMfcc.window=" << pTweak.featureMfcc.window << endl;
    OutS << "featureMfcc.method=" << pTweak.featureMfcc.method << endl;
    OutS << "featureMfcc.sclib_frequencyScale=" << (int)(pTweak.featureMfcc.sclib_frequencyScale) << endl;
    OutS << "featureMfcc.sclib_minFilterBankFrequency=" << pTweak.featureMfcc.sclib_minFilterBankFrequency << endl;
    OutS << "featureMfcc.sclib_maxFilterBankFrequency=" << pTweak.featureMfcc.sclib_maxFilterBankFrequency << endl;
    OutS << "featureMfcc.sclib_smoothing=" << (int)(pTweak.featureMfcc.sclib_smoothing) << endl;
    
    OutS << "featureNfr.FFTsize=" << pTweak.featureNfr.FFTsize << endl;
    OutS << "featureNfr.frameSize=" << pTweak.featureNfr.frameSize << endl;
    OutS << "featureNfr.frameStep=" << pTweak.featureNfr.frameStep << endl;
    OutS << "featureNfr.highCut=" << pTweak.featureNfr.highCut << endl;
    OutS << "featureNfr.lowCut=" << pTweak.featureNfr.lowCut << endl;
    OutS << "featureNfr.NFRthreshold=" << pTweak.featureNfr.NFRthreshold << endl;
    OutS << "featureNfr.preEmphasizeFactor=" << pTweak.featureNfr.preEmphasizeFactor << endl;
    OutS << "featureNfr.window=" << pTweak.featureNfr.window << endl;
    
    OutS << "featureSpectrum.FFTsize=" << pTweak.featureSpectrum.FFTsize << endl;
    OutS << "featureSpectrum.frameSize=" << pTweak.featureSpectrum.frameSize << endl;
    OutS << "featureSpectrum.frameStep=" << pTweak.featureSpectrum.frameStep << endl;
    OutS << "featureSpectrum.highCut=" << pTweak.featureSpectrum.highCut << endl;
    OutS << "featureSpectrum.lowCut=" << pTweak.featureSpectrum.lowCut << endl;
    OutS << "featureSpectrum.preEmphasizeFactor=" << pTweak.featureSpectrum.preEmphasizeFactor << endl;
    OutS << "featureSpectrum.window=" << pTweak.featureSpectrum.window << endl;
    OutS << "featureSpectrum.logarithmize=" << pTweak.featureSpectrum.logarithmize << endl;
    OutS << "featureSpectrum.createPhase=" << pTweak.featureSpectrum.createPhase << endl;
    
		OutS << "featureSpectrumFlux.FFTsize=" << pTweak.featureSpectrumFlux.FFTsize << endl;
    OutS << "featureSpectrumFlux.frameSize=" << pTweak.featureSpectrumFlux.frameSize << endl;
    OutS << "featureSpectrumFlux.frameStep=" << pTweak.featureSpectrumFlux.frameStep << endl;
    OutS << "featureSpectrumFlux.highCut=" << pTweak.featureSpectrumFlux.highCut << endl;
    OutS << "featureSpectrumFlux.lowCut=" << pTweak.featureSpectrumFlux.lowCut << endl;
    OutS << "featureSpectrumFlux.preEmphasizeFactor=" << pTweak.featureSpectrumFlux.preEmphasizeFactor << endl;
    OutS << "featureSpectrumFlux.window=" << pTweak.featureSpectrumFlux.window << endl;
    
	  OutS << "featureSte.frameSize=" << pTweak.featureSte.frameSize << endl; //by Bing
    OutS << "featureSte.frameStep=" << pTweak.featureSte.frameStep << endl; //by Bing
    OutS << "featureSte.highCut=" << pTweak.featureSte.highCut << endl; //by Bing
    OutS << "featureSte.lowCut=" << pTweak.featureSte.lowCut << endl; //by Bing
    OutS << "featureSte.preEmphasizeFactor=" << pTweak.featureSte.preEmphasizeFactor << endl; //by Bing
		OutS << "featureSte.useButterworth=" << pTweak.featureSte.useButterworth << endl;
		OutS << "featureSte.scaleResult=" << pTweak.featureSte.scaleResult << endl;
    
    OutS << "featureSubBandPower.FFTsize=" << pTweak.featureSubBandPower.FFTsize << endl;
    OutS << "featureSubBandPower.frameSize=" << pTweak.featureSubBandPower.frameSize << endl;
    OutS << "featureSubBandPower.frameStep=" << pTweak.featureSubBandPower.frameStep << endl;
    OutS << "featureSubBandPower.highCut=" << pTweak.featureSubBandPower.highCut << endl;
    OutS << "featureSubBandPower.lowCut=" << pTweak.featureSubBandPower.lowCut << endl;
    OutS << "featureSubBandPower.preEmphasizeFactor=" << pTweak.featureSubBandPower.preEmphasizeFactor << endl;
    OutS << "featureSubBandPower.window=" << pTweak.featureSubBandPower.window << endl;

    OutS << "featureZcr.frameSize=" << pTweak.featureZcr.frameSize << endl;
    OutS << "featureZcr.frameStep=" << pTweak.featureZcr.frameStep << endl;
    OutS << "featureZcr.highCut=" << pTweak.featureZcr.highCut << endl;
    OutS << "featureZcr.lowCut=" << pTweak.featureZcr.lowCut << endl;
    OutS << "featureZcr.preEmphasizeFactor=" << pTweak.featureZcr.preEmphasizeFactor << endl;
		OutS << "featureZcr.useChebyshev=" << pTweak.featureZcr.useChebyshev << endl;
		OutS << "featureZcr.scaleResult=" << pTweak.featureZcr.scaleResult << endl;

		OutS << "featureLpc.highCut=" << pTweak.featureLpc.highCut << endl;
		OutS << "featureLpc.lowCut=" << pTweak.featureLpc.lowCut << endl;
		OutS << "featureLpc.frameSize=" << pTweak.featureLpc.frameSize << endl; //by Jun
		OutS << "featureLpc.frameStep=" << pTweak.featureLpc.frameStep << endl; //by Jun
		OutS << "featureLpc.preEmphasizeFactor=" << pTweak.featureLpc.preEmphasizeFactor << endl; //by Jun
		OutS << "featureLpc.LPCorder=" << pTweak.featureLpc.LPCorder << endl; //by Jun
		OutS << "featureLpc.window=" << pTweak.featureLpc.window << endl;
		OutS << "featureLpc.computeGain=" << pTweak.featureLpc.computeGain << endl;

		OutS << "featureLpcResidual.highCut=" << pTweak.featureLpcResidual.highCut << endl;
		OutS << "featureLpcResidual.lowCut=" << pTweak.featureLpcResidual.lowCut << endl;
		OutS << "featureLpcResidual.frameSize=" << pTweak.featureLpcResidual.frameSize << endl; //by Jun
		OutS << "featureLpcResidual.frameStep=" << pTweak.featureLpcResidual.frameStep << endl; //by Jun
		OutS << "featureLpcResidual.preEmphasizeFactor=" << pTweak.featureLpcResidual.preEmphasizeFactor << endl; //by Jun
		OutS << "featureLpcResidual.order=" << pTweak.featureLpcResidual.order << endl; //by Jun

		OutS << "featureLsp.highCut=" << pTweak.featureLsp.highCut << endl;
		OutS << "featureLsp.lowCut=" << pTweak.featureLsp.lowCut << endl;
		OutS << "featureLsp.frameSize=" << pTweak.featureLsp.frameSize << endl;
		OutS << "featureLsp.frameStep=" << pTweak.featureLsp.frameStep << endl;
		OutS << "featureLsp.preEmphasizeFactor=" << pTweak.featureLsp.preEmphasizeFactor << endl;
		OutS << "featureLsp.LPCorder=" << pTweak.featureLsp.LPCorder << endl;
		OutS << "featureLsp.window=" << pTweak.featureLsp.window << endl;
		OutS << "featureLsp.method=" << pTweak.featureLsp.method << endl;
		OutS << "featureLsp.delta=" << pTweak.featureLsp.delta << endl;
		OutS << "featureLsp.bisections=" << pTweak.featureLsp.bisections << endl;
		OutS << "featureLsp.minSeparation=" << pTweak.featureLsp.minSeparation << endl;
		OutS << "featureLsp.maxLoops=" << pTweak.featureLsp.maxLoops << endl;
		
		OutS << "featurePitch.frameSize=" << pTweak.featurePitch.frameSize << endl;
		OutS << "featurePitch.frameStep=" << pTweak.featurePitch.frameStep << endl;
		OutS << "featurePitch.method=" << pTweak.featurePitch.method << endl;
		OutS << "featurePitch.sqrt=" << pTweak.featurePitch.sqrt << endl;
		OutS << "featurePitch.esps_cand_thresh=" << pTweak.featurePitch.esps_cand_thresh << endl;
		OutS << "featurePitch.esps_lag_weight=" << pTweak.featurePitch.esps_lag_weight << endl;
		OutS << "featurePitch.esps_freq_weight=" << pTweak.featurePitch.esps_freq_weight << endl;
		OutS << "featurePitch.esps_trans_cost=" << pTweak.featurePitch.esps_trans_cost << endl;
		OutS << "featurePitch.esps_trans_amp=" << pTweak.featurePitch.esps_trans_amp << endl;
		OutS << "featurePitch.esps_trans_spec=" << pTweak.featurePitch.esps_trans_spec << endl;
		OutS << "featurePitch.esps_voice_bias=" << pTweak.featurePitch.esps_voice_bias << endl;
		OutS << "featurePitch.esps_double_cost=" << pTweak.featurePitch.esps_double_cost << endl;
		OutS << "featurePitch.esps_min_f0=" << pTweak.featurePitch.esps_min_f0 << endl;
		OutS << "featurePitch.esps_max_f0=" << pTweak.featurePitch.esps_max_f0 << endl;
		OutS << "featurePitch.esps_n_cands=" << pTweak.featurePitch.esps_n_cands << endl;
		OutS << "featurePitch.esps_wind_dur=" << pTweak.featurePitch.esps_wind_dur << endl;
		
		OutS << "featureFormant.lowCut=" << pTweak.featureFormant.lowCut << endl;
		OutS << "featureFormant.highCut=" << pTweak.featureFormant.highCut << endl;
		OutS << "featureFormant.frameSize=" << pTweak.featureFormant.frameSize << endl;
		OutS << "featureFormant.frameStep=" << pTweak.featureFormant.frameStep << endl;
		OutS << "featureFormant.sampleRate=" << pTweak.featureFormant.sampleRate << endl;
		OutS << "featureFormant.preEmphasizeFactor=" << pTweak.featureFormant.preEmphasizeFactor << endl;
		OutS << "featureFormant.esps_cor_wdur=" << pTweak.featureFormant.esps_cor_wdur << endl;
		OutS << "featureFormant.esps_ds_freq=" << pTweak.featureFormant.esps_ds_freq << endl;
		OutS << "featureFormant.esps_frame_int=" << pTweak.featureFormant.esps_frame_int << endl;
		OutS << "featureFormant.esps_lpc_ord=" << pTweak.featureFormant.esps_lpc_ord << endl;
		OutS << "featureFormant.esps_lpc_type=" << pTweak.featureFormant.esps_lpc_type << endl;
		OutS << "featureFormant.esps_nform=" << pTweak.featureFormant.esps_nform << endl;
		OutS << "featureFormant.esps_nom_f1=" << pTweak.featureFormant.esps_nom_f1 << endl;
		OutS << "featureFormant.esps_w_type=" << pTweak.featureFormant.esps_w_type << endl;
		OutS << "featureFormant.esps_wdur=" << pTweak.featureFormant.esps_wdur << endl;
		
		OutS << "featureSdp.m=" << pTweak.featureSdp.m << endl; // by bing
		OutS << "featureSdp.lag=" << pTweak.featureSdp.lag << endl; // by bing
		OutS << "featureSdp.color=" << pTweak.featureSdp.color << endl; // by bing
		OutS << "featureSdp.frameSize=" << pTweak.featureSdp.frameSize << endl; // by bing
		OutS << "featureSdp.frameStep=" << pTweak.featureSdp.frameStep << endl; // by bing
		OutS << "featureSdp.n=" << pTweak.featureSdp.n<< endl; // by bing
		OutS << "featureSdp.pictureSize=" << pTweak.featureSdp.pictureSize << endl; // by bing
		OutS << "featureSdp.tau=" << pTweak.featureSdp.tau<< endl; // by bing
		OutS << "featureSdp.preEmphasizeFactor=" << pTweak.featureSdp.preEmphasizeFactor << endl; // by bing
		
		OutS << "featureSamples.frameSize=" << pTweak.featureSamples.frameSize << endl;
		OutS << "featureSamples.frameStep=" << pTweak.featureSamples.frameStep << endl;
		OutS << "featureSamples.highCut=" << pTweak.featureSamples.highCut << endl;
		OutS << "featureSamples.lowCut=" << pTweak.featureSamples.lowCut << endl;
		
		OutS << "general.firstScene=" << pTweak.general.firstScene << endl;
    OutS << "general.lastScene=" << pTweak.general.lastScene << endl;
    OutS << "general.sceneSelection=" << pTweak.general.sceneSelection << endl;
		OutS << "general.shortSpeechThreshold=" << pTweak.general.shortSpeechThreshold << endl;
    OutS << "general.pauseSilenceThreshold=" << pTweak.general.pauseSilenceThreshold << endl;
    OutS << "general.preClusteringResultsPrefix=" << sclib::nullFilter(pTweak.general.preClusteringResultsPrefix) << endl;
    OutS << "general.featurePrefix=" << sclib::nullFilter(pTweak.general.featurePrefix) << endl;
    
    OutS << "groundTruth.internalFrameSize=" << pTweak.groundTruth.internalFrameSize << endl;
    OutS << "groundTruth.pseudoSceneLength=" << pTweak.groundTruth.pseudoSceneLength << endl;
    OutS << "groundTruth.videoFrameMachineOffset=" << pTweak.groundTruth.videoFrameMachineOffset << endl;
    OutS << "groundTruth.storeProbabilityInformation=" << pTweak.groundTruth.storeProbabilityInformation << endl;
    
    OutS << "mixtureModelBgmm.EMthreshold=" << pTweak.mixtureModelBgmm.EMthreshold << endl;
    OutS << "mixtureModelBgmm.maxEMiterations=" << pTweak.mixtureModelBgmm.maxEMiterations << endl;
    OutS << "mixtureModelBgmm.varianceLimit=" << pTweak.mixtureModelBgmm.varianceLimit << endl;
    OutS << "mixtureModelBgmm.weightLimit=" << pTweak.mixtureModelBgmm.weightLimit << endl;
    OutS << "mixtureModelBgmm.fullCovariance=" << pTweak.mixtureModelBgmm.fullCovariance << endl;
    OutS << "mixtureModelBgmm.randomInitialization=" << pTweak.mixtureModelBgmm.randomInitialization << endl;
    
    OutS << "mixtureModelGmm.EMthreshold=" << pTweak.mixtureModelGmm.EMthreshold << endl;
    OutS << "mixtureModelGmm.maxEMiterations=" << pTweak.mixtureModelGmm.maxEMiterations << endl;
    OutS << "mixtureModelGmm.kMeansIterations=" << pTweak.mixtureModelGmm.kMeansIterations << endl;
    OutS << "mixtureModelGmm.varianceLimit=" << pTweak.mixtureModelGmm.varianceLimit << endl;
    OutS << "mixtureModelGmm.weightLimit=" << pTweak.mixtureModelGmm.weightLimit << endl;
    
    OutS << "mixtureModelMixMax.bgModelCombination=" << pTweak.mixtureModelMixMax.bgModelCombination << endl;
    OutS << "mixtureModelMixMax.EMthreshold=" << pTweak.mixtureModelMixMax.EMthreshold << endl;
    OutS << "mixtureModelMixMax.kMeansIterations=" << pTweak.mixtureModelMixMax.kMeansIterations << endl;
    OutS << "mixtureModelMixMax.maxEMiterations=" << pTweak.mixtureModelMixMax.maxEMiterations << endl;
    OutS << "mixtureModelMixMax.maxERFloops=" << pTweak.mixtureModelMixMax.maxERFloops << endl;
    OutS << "mixtureModelMixMax.noiseCorruptionType=" << (int)(pTweak.mixtureModelMixMax.noiseCorruptionType) << endl;
    OutS << "mixtureModelMixMax.varianceLimit=" << pTweak.mixtureModelMixMax.varianceLimit << endl;
    OutS << "mixtureModelMixMax.weightLimit=" << pTweak.mixtureModelMixMax.weightLimit << endl;
    
    OutS << "mixtureModelGmmubm.adaptMeans=" << pTweak.mixtureModelGmmubm.adaptMeans << endl;
    OutS << "mixtureModelGmmubm.adaptVariances=" << pTweak.mixtureModelGmmubm.adaptVariances << endl;
    OutS << "mixtureModelGmmubm.adaptWeights=" << pTweak.mixtureModelGmmubm.adaptWeights << endl;
    OutS << "mixtureModelGmmubm.relevanceFactor=" << pTweak.mixtureModelGmmubm.relevanceFactor << endl;
    OutS << "mixtureModelGmmubm.scoringMethod=" << pTweak.mixtureModelGmmubm.scoringMethod << endl;
    OutS << "mixtureModelGmmubm.topCmixtures=" << pTweak.mixtureModelGmmubm.topCmixtures << endl;
    OutS << "mixtureModelGmmubm.ubmFileName=" << sclib::nullFilter(pTweak.mixtureModelGmmubm.ubmFileName) << endl;
    OutS << "mixtureModelGmmubm.varianceLimit=" << pTweak.mixtureModelGmmubm.varianceLimit << endl;
    
    OutS << "mixtureModelMix2Max.EMthreshold=" << pTweak.mixtureModelMix2Max.EMthreshold << endl;
    OutS << "mixtureModelMix2Max.kMeansIterations=" << pTweak.mixtureModelMix2Max.kMeansIterations << endl;
    OutS << "mixtureModelMix2Max.maxEMiterations=" << pTweak.mixtureModelMix2Max.maxEMiterations << endl;
    OutS << "mixtureModelMix2Max.maxERFloops=" << pTweak.mixtureModelMix2Max.maxERFloops << endl;
    OutS << "mixtureModelMix2Max.varianceLimit=" << pTweak.mixtureModelMix2Max.varianceLimit << endl;
    OutS << "mixtureModelMix2Max.weightLimit=" << pTweak.mixtureModelMix2Max.weightLimit << endl;
    OutS << "mixtureModelMix2Max.bgModelCombination=" << pTweak.mixtureModelMix2Max.bgModelCombination << endl;
    
    OutS << "mixtureModelMix2MaxEx.bgModelCombination=" << pTweak.mixtureModelMix2MaxEx.bgModelCombination << endl;
    OutS << "mixtureModelMix2MaxEx.EMthreshold=" << pTweak.mixtureModelMix2MaxEx.EMthreshold << endl;
    OutS << "mixtureModelMix2MaxEx.kMeansIterations=" << pTweak.mixtureModelMix2MaxEx.kMeansIterations << endl;
    OutS << "mixtureModelMix2MaxEx.maxEMiterations=" << pTweak.mixtureModelMix2MaxEx.maxEMiterations << endl;
    OutS << "mixtureModelMix2MaxEx.maxERFloops=" << pTweak.mixtureModelMix2MaxEx.maxERFloops << endl;
    OutS << "mixtureModelMix2MaxEx.varianceLimit=" << pTweak.mixtureModelMix2MaxEx.varianceLimit << endl;
    OutS << "mixtureModelMix2MaxEx.weightLimit=" << pTweak.mixtureModelMix2MaxEx.weightLimit << endl;
    
    OutS << "modelHandler.backgroundModelType=" << pTweak.modelHandler.backgroundModelType << endl;
    OutS << "modelHandler.foregroundModelType=" << pTweak.modelHandler.foregroundModelType << endl;
    OutS << "modelHandler.maxNoiseModelOrder=" << pTweak.modelHandler.maxNoiseModelOrder << endl;
    OutS << "modelHandler.maxSpeakerModelOrder=" << pTweak.modelHandler.maxSpeakerModelOrder << endl;
    OutS << "modelHandler.msPerGaussian=" << pTweak.modelHandler.msPerGaussian << endl;
    OutS << "modelHandler.orderGuessEMsteps=" << pTweak.modelHandler.orderGuessEMsteps << endl;
    OutS << "modelHandler.orderGuessMode=" << pTweak.modelHandler.orderGuessMode << endl;
    OutS << "modelHandler.SNRthreshold=" << pTweak.modelHandler.SNRthreshold << endl;
    OutS << "modelHandler.onlyThisSpeaker=" << sclib::nullFilter(pTweak.modelHandler.onlyThisSpeaker) << endl;
		OutS << "modelHandler.speakerModelFeature=" << pTweak.modelHandler.speakerModelFeature << endl;
		OutS << "modelHandler.outlierRemovalMode=" << (int)(pTweak.modelHandler.outlierRemovalMode) << endl;
		
		OutS << "modelHmm.stateCount=" << pTweak.modelHmm.stateCount << endl;
		OutS << "modelHmm.transitionStructure=" << pTweak.modelHmm.transitionStructure << endl;
		OutS << "modelHmm.mixturesPerState=" << pTweak.modelHmm.mixturesPerState << endl;
		OutS << "modelHmm.useOrthogonalTransform=" << pTweak.modelHmm.useOrthogonalTransform << endl;
		OutS << "modelHmm.leftToRight=" << pTweak.modelHmm.leftToRight << endl;
		OutS << "modelHmm.maxIterations=" << pTweak.modelHmm.maxIterations << endl;
		OutS << "modelHmm.verbose=" << pTweak.modelHmm.verbose << endl;
		
		OutS << "modelVq.codebookSize=" << pTweak.modelVq.codebookSize << endl;
		OutS << "modelVq.splitMethod=" << pTweak.modelVq.splitMethod << endl;
		OutS << "modelVq.maxIterations=" << pTweak.modelVq.maxIterations << endl;
		
		OutS << "modelSvm.distanceBasedTesting=" << pTweak.modelSvm.distanceBasedTesting << endl;
		OutS << "modelSvm.doParameterSearch=" << pTweak.modelSvm.doParameterSearch << endl;
		
		OutS << "modelTime.subModelType=" << (int)(pTweak.modelTime.subModelType) << endl;
		OutS << "modelTime.syllableLength=" << pTweak.modelTime.syllableLength << endl;
		OutS << "modelTime.trajectoryStep=" << pTweak.modelTime.trajectoryStep << endl;
		OutS << "modelTime.removeTiming=" << pTweak.modelTime.removeTiming << endl;
		OutS << "modelTime.templateCount=" << pTweak.modelTime.templateCount << endl;
		OutS << "modelTime.clusteringIterations=" << pTweak.modelTime.clusteringIterations << endl;
		OutS << "modelTime.replaceTrainingData=" << pTweak.modelTime.replaceTrainingData << endl;
		OutS << "modelTime.checkForTrajectorization=" << pTweak.modelTime.checkForTrajectorization << endl;
		OutS << "modelTime.worldModelFile=" << sclib::nullFilter(pTweak.modelTime.worldModelFile) << endl;
		OutS << "modelTime.normalizationFile=" << sclib::nullFilter(pTweak.modelTime.normalizationFile) << endl;
		
		OutS << "modelPareto.useMarginalDistributions=" << pTweak.modelPareto.useMarginalDistributions << endl;
		
    OutS << "modelQgmm.deltaBIClambda=" << pTweak.modelQgmm.deltaBIClambda << endl;
    OutS << "modelQgmm.maxMixtures=" << pTweak.modelQgmm.maxMixtures << endl;
    OutS << "modelQgmm.percentDifference=" << pTweak.modelQgmm.percentDifference << endl;
    
    OutS << "resampling.fastConversion=" << pTweak.resampling.fastConversion << endl;
    
    OutS << "score.BBNmetricLambda=" << pTweak.score.BBNmetricLambda << endl;
    
    OutS << "segmentationAudioTypeLzl.actionModelFileName=" << sclib::nullFilter(pTweak.segmentationAudioTypeLzl.actionModelFileName) << endl;
    OutS << "segmentationAudioTypeLzl.classifierFileName=" << sclib::nullFilter(pTweak.segmentationAudioTypeLzl.classifierFileName) << endl;
    OutS << "segmentationAudioTypeLzl.featureFileName=" << sclib::nullFilter(pTweak.segmentationAudioTypeLzl.featureFileName) << endl;
    OutS << "segmentationAudioTypeLzl.normalizationFileName=" << sclib::nullFilter(pTweak.segmentationAudioTypeLzl.normalizationFileName) << endl;
		OutS << "segmentationAudioTypeLzl.speechSpecificity=" << pTweak.segmentationAudioTypeLzl.speechSpecificity << endl;
		OutS << "segmentationAudioTypeLzl.musicSpecificity=" << pTweak.segmentationAudioTypeLzl.musicSpecificity << endl;
	  OutS << "segmentationAudioTypeLzl.subClipLength=" << pTweak.segmentationAudioTypeLzl.subClipLength << endl;
		OutS << "segmentationAudioTypeLzl.bbParameters.FFTsize=" << pTweak.segmentationAudioTypeLzl.bbParameters.FFTsize << endl;
		OutS << "segmentationAudioTypeLzl.bbParameters.frameSize=" << pTweak.segmentationAudioTypeLzl.bbParameters.frameSize << endl;
		OutS << "segmentationAudioTypeLzl.bbParameters.frameSize=" << pTweak.segmentationAudioTypeLzl.bbParameters.frameSize << endl;
		OutS << "segmentationAudioTypeLzl.bbParameters.preEmphasizeFactor=" << pTweak.segmentationAudioTypeLzl.bbParameters.preEmphasizeFactor << endl;
		OutS << "segmentationAudioTypeLzl.bbParameters.window=" << pTweak.segmentationAudioTypeLzl.bbParameters.window << endl;
		OutS << "segmentationAudioTypeLzl.bbParameters.sampleRate=" << pTweak.segmentationAudioTypeLzl.bbParameters.sampleRate << endl;
		OutS << "segmentationAudioTypeLzl.bpParameters.FFTsize=" << pTweak.segmentationAudioTypeLzl.bpParameters.FFTsize << endl;
		OutS << "segmentationAudioTypeLzl.bpParameters.frameSize=" << pTweak.segmentationAudioTypeLzl.bpParameters.frameSize << endl;
		OutS << "segmentationAudioTypeLzl.bpParameters.frameStep=" << pTweak.segmentationAudioTypeLzl.bpParameters.frameStep << endl;
		OutS << "segmentationAudioTypeLzl.bpParameters.preEmphasizeFactor=" << pTweak.segmentationAudioTypeLzl.bpParameters.preEmphasizeFactor << endl;
		OutS << "segmentationAudioTypeLzl.bpParameters.window=" << pTweak.segmentationAudioTypeLzl.bpParameters.window << endl;
		OutS << "segmentationAudioTypeLzl.bpParameters.sampleRate=" << pTweak.segmentationAudioTypeLzl.bpParameters.sampleRate << endl;
		OutS << "segmentationAudioTypeLzl.mfccParameters.fftSize=" << pTweak.segmentationAudioTypeLzl.mfccParameters.fftSize << endl;
		OutS << "segmentationAudioTypeLzl.mfccParameters.frameSize=" << pTweak.segmentationAudioTypeLzl.mfccParameters.frameSize << endl;
		OutS << "segmentationAudioTypeLzl.mfccParameters.frameStep=" << pTweak.segmentationAudioTypeLzl.mfccParameters.frameStep << endl;
		OutS << "segmentationAudioTypeLzl.mfccParameters.MFCCorder=" << pTweak.segmentationAudioTypeLzl.mfccParameters.MFCCorder << endl;
		OutS << "segmentationAudioTypeLzl.mfccParameters.addDeltaDeltas=" << pTweak.segmentationAudioTypeLzl.mfccParameters.addDeltaDeltas << endl;
		OutS << "segmentationAudioTypeLzl.mfccParameters.addDeltas=" << pTweak.segmentationAudioTypeLzl.mfccParameters.addDeltas << endl;
		OutS << "segmentationAudioTypeLzl.mfccParameters.method=" << pTweak.segmentationAudioTypeLzl.mfccParameters.method << endl;
		OutS << "segmentationAudioTypeLzl.mfccParameters.dEnergy=" << pTweak.segmentationAudioTypeLzl.mfccParameters.dEnergy << endl;
		OutS << "segmentationAudioTypeLzl.mfccParameters.CMN=" << pTweak.segmentationAudioTypeLzl.mfccParameters.CMN << endl;
		OutS << "segmentationAudioTypeLzl.mfccParameters.filterBankSize=" << pTweak.segmentationAudioTypeLzl.mfccParameters.filterBankSize << endl;
		OutS << "segmentationAudioTypeLzl.mfccParameters.preEmphasizeFactor=" << pTweak.segmentationAudioTypeLzl.mfccParameters.preEmphasizeFactor << endl;
		OutS << "segmentationAudioTypeLzl.mfccParameters.sclib_frequencyScale=" << (int)(pTweak.segmentationAudioTypeLzl.mfccParameters.sclib_frequencyScale) << endl;
		OutS << "segmentationAudioTypeLzl.mfccParameters.sclib_smoothing=" << (int)(pTweak.segmentationAudioTypeLzl.mfccParameters.sclib_smoothing) << endl;
		OutS << "segmentationAudioTypeLzl.mfccParameters.sclib_maxFilterBankFrequency=" << pTweak.segmentationAudioTypeLzl.mfccParameters.sclib_maxFilterBankFrequency << endl;
		OutS << "segmentationAudioTypeLzl.mfccParameters.sclib_minFilterBankFrequency=" << pTweak.segmentationAudioTypeLzl.mfccParameters.sclib_minFilterBankFrequency << endl;
		OutS << "segmentationAudioTypeLzl.mfccParameters.window=" << pTweak.segmentationAudioTypeLzl.mfccParameters.window << endl;
		OutS << "segmentationAudioTypeLzl.mfccParameters.coeffSelection=" << "'" << pTweak.segmentationAudioTypeLzl.mfccParameters.coeffSelection << "'" << endl;
		OutS << "segmentationAudioTypeLzl.mfccParameters.sampleRate=" << pTweak.segmentationAudioTypeLzl.mfccParameters.sampleRate << endl;
		OutS << "segmentationAudioTypeLzl.nfrParameters.FFTsize=" << pTweak.segmentationAudioTypeLzl.nfrParameters.FFTsize << endl;
		OutS << "segmentationAudioTypeLzl.nfrParameters.frameSize=" << pTweak.segmentationAudioTypeLzl.nfrParameters.frameSize << endl;
		OutS << "segmentationAudioTypeLzl.nfrParameters.frameStep=" << pTweak.segmentationAudioTypeLzl.nfrParameters.frameStep << endl;
		OutS << "segmentationAudioTypeLzl.nfrParameters.NFRthreshold=" << pTweak.segmentationAudioTypeLzl.nfrParameters.NFRthreshold << endl;
		OutS << "segmentationAudioTypeLzl.nfrParameters.preEmphasizeFactor=" << pTweak.segmentationAudioTypeLzl.nfrParameters.preEmphasizeFactor << endl;
		OutS << "segmentationAudioTypeLzl.nfrParameters.window=" << pTweak.segmentationAudioTypeLzl.nfrParameters.window << endl;
		OutS << "segmentationAudioTypeLzl.nfrParameters.sampleRate=" << pTweak.segmentationAudioTypeLzl.nfrParameters.sampleRate << endl;
		OutS << "segmentationAudioTypeLzl.sbpParameters.FFTsize=" << pTweak.segmentationAudioTypeLzl.sbpParameters.FFTsize << endl;
		OutS << "segmentationAudioTypeLzl.sbpParameters.frameSize=" << pTweak.segmentationAudioTypeLzl.sbpParameters.frameSize << endl;
		OutS << "segmentationAudioTypeLzl.sbpParameters.frameStep=" << pTweak.segmentationAudioTypeLzl.sbpParameters.frameStep << endl;
		OutS << "segmentationAudioTypeLzl.sbpParameters.preEmphasizeFactor=" << pTweak.segmentationAudioTypeLzl.sbpParameters.preEmphasizeFactor << endl;
		OutS << "segmentationAudioTypeLzl.sbpParameters.window=" << pTweak.segmentationAudioTypeLzl.sbpParameters.window << endl;
		OutS << "segmentationAudioTypeLzl.sbpParameters.sampleRate=" << pTweak.segmentationAudioTypeLzl.sbpParameters.sampleRate << endl;
		OutS << "segmentationAudioTypeLzl.sfParameters.FFTsize=" << pTweak.segmentationAudioTypeLzl.sfParameters.FFTsize << endl;
		OutS << "segmentationAudioTypeLzl.sfParameters.frameSize=" << pTweak.segmentationAudioTypeLzl.sfParameters.frameSize << endl;
		OutS << "segmentationAudioTypeLzl.sfParameters.frameStep=" << pTweak.segmentationAudioTypeLzl.sfParameters.frameStep << endl;
		OutS << "segmentationAudioTypeLzl.sfParameters.preEmphasizeFactor=" << pTweak.segmentationAudioTypeLzl.sfParameters.preEmphasizeFactor << endl;
		OutS << "segmentationAudioTypeLzl.sfParameters.window=" << pTweak.segmentationAudioTypeLzl.sfParameters.window << endl;
		OutS << "segmentationAudioTypeLzl.sfParameters.sampleRate=" << pTweak.segmentationAudioTypeLzl.sfParameters.sampleRate << endl;
		OutS << "segmentationAudioTypeLzl.zcrParameters.frameSize=" << pTweak.segmentationAudioTypeLzl.zcrParameters.frameSize << endl;
		OutS << "segmentationAudioTypeLzl.zcrParameters.frameStep=" << pTweak.segmentationAudioTypeLzl.zcrParameters.frameStep << endl;
		OutS << "segmentationAudioTypeLzl.zcrParameters.preEmphasizeFactor=" << pTweak.segmentationAudioTypeLzl.zcrParameters.preEmphasizeFactor << endl;
		OutS << "segmentationAudioTypeLzl.zcrParameters.scaleResult=" << pTweak.segmentationAudioTypeLzl.zcrParameters.scaleResult << endl;
		OutS << "segmentationAudioTypeLzl.zcrParameters.useChebyshev=" << pTweak.segmentationAudioTypeLzl.zcrParameters.useChebyshev << endl;
		OutS << "segmentationAudioTypeLzl.zcrParameters.sampleRate=" << pTweak.segmentationAudioTypeLzl.zcrParameters.sampleRate << endl;
		OutS << "segmentationAudioTypeLzl.svmParameters.C=" << pTweak.segmentationAudioTypeLzl.svmParameters.C << endl;
		OutS << "segmentationAudioTypeLzl.svmParameters.cache_size=" << pTweak.segmentationAudioTypeLzl.svmParameters.cache_size << endl;
		OutS << "segmentationAudioTypeLzl.svmParameters.coef0=" << pTweak.segmentationAudioTypeLzl.svmParameters.coef0 << endl;
		OutS << "segmentationAudioTypeLzl.svmParameters.cvFolds=" << pTweak.segmentationAudioTypeLzl.svmParameters.cvFolds << endl;
		OutS << "segmentationAudioTypeLzl.svmParameters.cvMaxDatasetSize=" << pTweak.segmentationAudioTypeLzl.svmParameters.cvMaxDatasetSize << endl;
		OutS << "segmentationAudioTypeLzl.svmParameters.cvCoarseCMin=" << pTweak.segmentationAudioTypeLzl.svmParameters.cvCoarseCMin << endl;
		OutS << "segmentationAudioTypeLzl.svmParameters.cvCoarseCMax=" << pTweak.segmentationAudioTypeLzl.svmParameters.cvCoarseCMax << endl;
		OutS << "segmentationAudioTypeLzl.svmParameters.cvCoarseCStep=" << pTweak.segmentationAudioTypeLzl.svmParameters.cvCoarseCStep << endl;
		OutS << "segmentationAudioTypeLzl.svmParameters.cvCoarseGammaMin=" << pTweak.segmentationAudioTypeLzl.svmParameters.cvCoarseGammaMin << endl;
		OutS << "segmentationAudioTypeLzl.svmParameters.cvCoarseGammaMax=" << pTweak.segmentationAudioTypeLzl.svmParameters.cvCoarseGammaMax << endl;
		OutS << "segmentationAudioTypeLzl.svmParameters.cvCoarseGammaStep=" << pTweak.segmentationAudioTypeLzl.svmParameters.cvCoarseGammaStep << endl;
		OutS << "segmentationAudioTypeLzl.svmParameters.cvFineCStep=" << pTweak.segmentationAudioTypeLzl.svmParameters.cvFineCStep << endl;
		OutS << "segmentationAudioTypeLzl.svmParameters.cvFineCRadius=" << pTweak.segmentationAudioTypeLzl.svmParameters.cvFineCRadius << endl;
		OutS << "segmentationAudioTypeLzl.svmParameters.cvFineGammaStep=" << pTweak.segmentationAudioTypeLzl.svmParameters.cvFineGammaStep << endl;
		OutS << "segmentationAudioTypeLzl.svmParameters.cvFineGammaRadius=" << pTweak.segmentationAudioTypeLzl.svmParameters.cvFineGammaRadius << endl;
		OutS << "segmentationAudioTypeLzl.svmParameters.degree=" << pTweak.segmentationAudioTypeLzl.svmParameters.degree << endl;
		OutS << "segmentationAudioTypeLzl.svmParameters.doCV=" << pTweak.segmentationAudioTypeLzl.svmParameters.doCV << endl;
		OutS << "segmentationAudioTypeLzl.svmParameters.eps=" << pTweak.segmentationAudioTypeLzl.svmParameters.eps << endl;
		OutS << "segmentationAudioTypeLzl.svmParameters.gamma=" << pTweak.segmentationAudioTypeLzl.svmParameters.gamma << endl;
		OutS << "segmentationAudioTypeLzl.svmParameters.kernel_type=" << pTweak.segmentationAudioTypeLzl.svmParameters.kernel_type << endl;
		OutS << "segmentationAudioTypeLzl.svmParameters.nu=" << pTweak.segmentationAudioTypeLzl.svmParameters.nu << endl;
		OutS << "segmentationAudioTypeLzl.svmParameters.p=" << pTweak.segmentationAudioTypeLzl.svmParameters.p << endl;
		OutS << "segmentationAudioTypeLzl.svmParameters.probability=" << pTweak.segmentationAudioTypeLzl.svmParameters.probability << endl;
		OutS << "segmentationAudioTypeLzl.svmParameters.shrinking=" << pTweak.segmentationAudioTypeLzl.svmParameters.shrinking << endl;
		OutS << "segmentationAudioTypeLzl.svmParameters.svm_type=" << pTweak.segmentationAudioTypeLzl.svmParameters.svm_type << endl;
		OutS << "segmentationAudioTypeLzl.svmParameters.oneClassGammaSearchMaxIterations=" << pTweak.segmentationAudioTypeLzl.svmParameters.oneClassGammaSearchMaxIterations << endl;
		OutS << "segmentationAudioTypeLzl.svmParameters.oneClassGammaSearchRepeats=" << pTweak.segmentationAudioTypeLzl.svmParameters.oneClassGammaSearchRepeats << endl;

    OutS << "segmentationChangesLz.detectorWindowLength=" << pTweak.segmentationChangesLz.detectorWindowLength << endl;
    OutS << "segmentationChangesLz.detectorWindowStep=" << pTweak.segmentationChangesLz.detectorWindowStep << endl;
    OutS << "segmentationChangesLz.adaptiveThresholdAlpha=" << pTweak.segmentationChangesLz.adaptiveThresholdAlpha << endl;
    OutS << "segmentationChangesLz.bayesianThreshold=" << pTweak.segmentationChangesLz.bayesianThreshold << endl;
		OutS << "segmentationChangesLz.lastNdistances=" << pTweak.segmentationChangesLz.lastNdistances << endl;
    OutS << "segmentationChangesLz.priorsFileName=" << sclib::nullFilter(pTweak.segmentationChangesLz.priorsFileName) << endl;
    OutS << "segmentationChangesLz.time2ChangeModelFileName=" << sclib::nullFilter(pTweak.segmentationChangesLz.time2ChangeModelFileName) << endl;
		
		OutS << "segmentationChangesKbk.r=" << pTweak.segmentationChangesKbk.r << endl;
		OutS << "segmentationChangesKbk.lambda=" << pTweak.segmentationChangesKbk.lambda << endl;
		OutS << "segmentationChangesKbk.tolerance=" << pTweak.segmentationChangesKbk.tolerance << endl;
		OutS << "segmentationChangesKbk.mfccParameters.frameSize=" << pTweak.segmentationChangesKbk.mfccParameters.frameSize << endl;
		OutS << "segmentationChangesKbk.mfccParameters.frameStep=" << pTweak.segmentationChangesKbk.mfccParameters.frameStep << endl;
		OutS << "segmentationChangesKbk.mfccParameters.MFCCorder=" << pTweak.segmentationChangesKbk.mfccParameters.MFCCorder << endl;
		OutS << "segmentationChangesKbk.mfccParameters.addDeltaDeltas=" << pTweak.segmentationChangesKbk.mfccParameters.addDeltaDeltas << endl;
		OutS << "segmentationChangesKbk.mfccParameters.addDeltas=" << pTweak.segmentationChangesKbk.mfccParameters.addDeltas << endl;
		OutS << "segmentationChangesKbk.mfccParameters.method=" << pTweak.segmentationChangesKbk.mfccParameters.method << endl;
		OutS << "segmentationChangesKbk.mfccParameters.dEnergy=" << pTweak.segmentationChangesKbk.mfccParameters.dEnergy << endl;
		OutS << "segmentationChangesKbk.mfccParameters.CMN=" << pTweak.segmentationChangesKbk.mfccParameters.CMN << endl;
		OutS << "segmentationChangesKbk.mfccParameters.fftSize=" << pTweak.segmentationChangesKbk.mfccParameters.fftSize << endl;
		OutS << "segmentationChangesKbk.mfccParameters.filterBankSize=" << pTweak.segmentationChangesKbk.mfccParameters.filterBankSize << endl;
		OutS << "segmentationChangesKbk.mfccParameters.preEmphasizeFactor=" << pTweak.segmentationChangesKbk.mfccParameters.preEmphasizeFactor << endl;
		OutS << "segmentationChangesKbk.mfccParameters.sclib_frequencyScale=" << (int)(pTweak.segmentationChangesKbk.mfccParameters.sclib_frequencyScale) << endl;
		OutS << "segmentationChangesKbk.mfccParameters.sclib_smoothing=" << (int)(pTweak.segmentationChangesKbk.mfccParameters.sclib_smoothing) << endl;
		OutS << "segmentationChangesKbk.mfccParameters.window=" << pTweak.segmentationChangesKbk.mfccParameters.window << endl;
		OutS << "segmentationChangesKbk.mfccParameters.coeffSelection=" << "'" << pTweak.segmentationChangesKbk.mfccParameters.coeffSelection << "'" << endl;
		OutS << "segmentationChangesKbk.mfccParameters.sclib_maxFilterBankFrequency=" << pTweak.segmentationChangesKbk.mfccParameters.sclib_maxFilterBankFrequency << endl;
		OutS << "segmentationChangesKbk.mfccParameters.sclib_minFilterBankFrequency=" << pTweak.segmentationChangesKbk.mfccParameters.sclib_minFilterBankFrequency << endl;
		
    OutS << "segmentationHandler.audioTypeMode=" << pTweak.segmentationHandler.audioTypeMode << endl;
    OutS << "segmentationHandler.changeDetectorMode=" << pTweak.segmentationHandler.changeDetectorMode << endl;
    OutS << "segmentationHandler.silenceDetectorMode=" << pTweak.segmentationHandler.silenceDetectorMode << endl;
    OutS << "segmentationHandler.vUvDetectorMode=" << pTweak.segmentationHandler.vUvDetectorMode << endl;

    OutS << "segmentationSilenceLnk.energyQuantizationLevel=" << pTweak.segmentationSilenceLnk.energyQuantizationLevel << endl;

    OutS << "segmentationSilenceLzl.energySilenceThreshold=" << pTweak.segmentationSilenceLzl.energySilenceThreshold << endl;
    OutS << "segmentationSilenceLzl.zcrSilenceThreshold=" << pTweak.segmentationSilenceLzl.zcrSilenceThreshold << endl;
		OutS << "segmentationSilenceLzl.specificity=" << pTweak.segmentationSilenceLzl.specificity << endl;
    OutS << "segmentationSilenceLzl.zcrParameters.frameSize=" << pTweak.segmentationSilenceLzl.zcrParameters.frameSize << endl;
    OutS << "segmentationSilenceLzl.zcrParameters.frameStep=" << pTweak.segmentationSilenceLzl.zcrParameters.frameStep << endl;
    OutS << "segmentationSilenceLzl.zcrParameters.sampleRate=" << pTweak.segmentationSilenceLzl.zcrParameters.sampleRate << endl;
    OutS << "segmentationSilenceLzl.zcrParameters.preEmphasizeFactor=" << pTweak.segmentationSilenceLzl.zcrParameters.preEmphasizeFactor << endl;
    OutS << "segmentationSilenceLzl.zcrParameters.useChebyshev=" << pTweak.segmentationSilenceLzl.zcrParameters.useChebyshev << endl;
    OutS << "segmentationSilenceLzl.zcrParameters.scaleResult=" << pTweak.segmentationSilenceLzl.zcrParameters.scaleResult << endl;
    OutS << "segmentationSilenceLzl.steParameters.frameSize=" << pTweak.segmentationSilenceLzl.steParameters.frameSize << endl;
    OutS << "segmentationSilenceLzl.steParameters.frameStep=" << pTweak.segmentationSilenceLzl.steParameters.frameStep << endl;
    OutS << "segmentationSilenceLzl.steParameters.sampleRate=" << pTweak.segmentationSilenceLzl.steParameters.sampleRate << endl;
    OutS << "segmentationSilenceLzl.steParameters.preEmphasizeFactor=" << pTweak.segmentationSilenceLzl.steParameters.preEmphasizeFactor << endl;
    OutS << "segmentationSilenceLzl.steParameters.useButterworth=" << pTweak.segmentationSilenceLzl.steParameters.useButterworth << endl;
    OutS << "segmentationSilenceLzl.steParameters.scaleResult=" << pTweak.segmentationSilenceLzl.steParameters.scaleResult << endl;

    OutS << "segmentationVUvLnk.energyQuantizationLevel=" << pTweak.segmentationVUvLnk.energyQuantizationLevel << endl;

		OutS << "segmentationVUvEsps.pitchParameters.frameSize=" << pTweak.segmentationVUvEsps.pitchParameters.frameSize << endl;
		OutS << "segmentationVUvEsps.pitchParameters.frameStep=" << pTweak.segmentationVUvEsps.pitchParameters.frameStep << endl;
		OutS << "segmentationVUvEsps.pitchParameters.method=" << pTweak.segmentationVUvEsps.pitchParameters.method << endl;
		OutS << "segmentationVUvEsps.pitchParameters.esps_cand_thresh=" << pTweak.segmentationVUvEsps.pitchParameters.esps_cand_thresh << endl;
		OutS << "segmentationVUvEsps.pitchParameters.esps_lag_weight=" << pTweak.segmentationVUvEsps.pitchParameters.esps_lag_weight << endl;
		OutS << "segmentationVUvEsps.pitchParameters.esps_freq_weight=" << pTweak.segmentationVUvEsps.pitchParameters.esps_freq_weight << endl;
		OutS << "segmentationVUvEsps.pitchParameters.esps_trans_cost=" << pTweak.segmentationVUvEsps.pitchParameters.esps_trans_cost << endl;
		OutS << "segmentationVUvEsps.pitchParameters.esps_trans_amp=" << pTweak.segmentationVUvEsps.pitchParameters.esps_trans_amp << endl;
		OutS << "segmentationVUvEsps.pitchParameters.esps_trans_spec=" << pTweak.segmentationVUvEsps.pitchParameters.esps_trans_spec << endl;
		OutS << "segmentationVUvEsps.pitchParameters.esps_voice_bias=" << pTweak.segmentationVUvEsps.pitchParameters.esps_voice_bias << endl;
		OutS << "segmentationVUvEsps.pitchParameters.esps_double_cost=" << pTweak.segmentationVUvEsps.pitchParameters.esps_double_cost << endl;
		OutS << "segmentationVUvEsps.pitchParameters.esps_min_f0=" << pTweak.segmentationVUvEsps.pitchParameters.esps_min_f0 << endl;
		OutS << "segmentationVUvEsps.pitchParameters.esps_max_f0=" << pTweak.segmentationVUvEsps.pitchParameters.esps_max_f0 << endl;
		OutS << "segmentationVUvEsps.pitchParameters.esps_n_cands=" << pTweak.segmentationVUvEsps.pitchParameters.esps_n_cands << endl;
		OutS << "segmentationVUvEsps.pitchParameters.esps_wind_dur=" << pTweak.segmentationVUvEsps.pitchParameters.esps_wind_dur << endl;

		OutS << "signalHandler.forceSampleRate=" << pTweak.signalHandler.forceSampleRate << endl;

		OutS << "signalMpeg.fastSeeking=" << pTweak.signalMpeg.fastSeeking << endl;
		OutS << "signalMpeg.hqResampling=" << pTweak.signalMpeg.hqResampling << endl;
		OutS << "signalMpeg.outputChannelCount=" << pTweak.signalMpeg.outputChannelCount << endl;
		OutS << "signalMpeg.outputSampleRate=" << pTweak.signalMpeg.outputSampleRate << endl;

    OutS << "speakerClusterer.doClustering=" << pTweak.speakerClusterer.doClustering << endl;
    OutS << "speakerClusterer.constructNonOverlappingClusters=" << pTweak.speakerClusterer.constructNonOverlappingClusters << endl;
    OutS << "speakerClusterer.speechSegLengthThreshold=" << pTweak.speakerClusterer.speechSegLengthThreshold << endl;
    OutS << "speakerClusterer.distanceMeasure=" << pTweak.speakerClusterer.distanceMeasure << endl;
    OutS << "speakerClusterer.globalCriterion=" << pTweak.speakerClusterer.globalCriterion << endl;
    OutS << "speakerClusterer.terminationCriterion=" << pTweak.speakerClusterer.terminationCriterion << endl;
		OutS << "speakerClusterer.linkageMode=" << (int)(pTweak.speakerClusterer.linkageMode) << endl;
		OutS << "speakerClusterer.specificity=" << pTweak.speakerClusterer.specificity << endl;
    OutS << "speakerClusterer.firstDistanceMatrixPrefix=" << sclib::nullFilter(pTweak.speakerClusterer.firstDistanceMatrixPrefix) << endl;

    OutS << "speakerIdentification.dNormSampleCount=" << pTweak.speakerIdentification.dNormSampleCount << endl;
    OutS << "speakerIdentification.doIdentification=" << pTweak.speakerIdentification.doIdentification << endl;
    OutS << "speakerIdentification.normalizationMode=" << pTweak.speakerIdentification.normalizationMode << endl;
    OutS << "speakerIdentification.useUBMs=" << pTweak.speakerIdentification.useUBMs << endl;

    OutS << "transform.taperingLength=" << pTweak.transform.taperingLength << endl;
	}

	return(OutS);
}

//====================================================================================================================
//	Set the parameter speciffied by it's name to the given value, which is converted implicitly to the needed type
//  The return-value indicates if the parameter existed or the value was valid (false, if one of this failed, else 
//  true)
//====================================================================================================================
bool SC_TweakableParameters::setByName(const char* parameterName, const char* value) {
  double dValue;
	float fValue;
	int iValue;
	char lCaseName[sclib::bufferSize];
	//static std::map<const char*, int, sclib::ltstr> this->parameters.map; //sclib::ltstr provides a less-then comparison on char*-strings, otherwise only adresses get compared, which fails always!

	//if (this->parameters.map.size() == 0) { //initialize the list once
	//
	//}

	sprintf(lCaseName, "%s", parameterName); //so that the const char remains const
	switch (this->parameters.map[sclib::lCase(lCaseName)]) {
		case 0:
			return false; //hits if the name wasn't found
		case 1:
			this->classifierAdaBoost.maxWeakClassifiers = atoi(value);
			return true;
		case 2:
			this->classifierAdaBoost.weakClassifierType = atoi(value);
			return true;
		case 3:
			this->classifierDecisionStump.splitMultiway = sclib::atob(value);
			return true;
		case 4:
			this->classifierMl.maxOrder = atoi(value);
			return true;
		case 5:
			this->classifierMl.minOrder = atoi(value);
			return true;
		case 6:
			this->classifierMl.modelType = atol(value);
			return true;
		case 7:
			this->classifierSvm.C = atof(value);
			return true;
		case 8:
			this->classifierSvm.cache_size = atof(value);
			return true;
		case 9:
			this->classifierSvm.coef0 = atof(value);
			return true;
		case 10:
			this->classifierSvm.cvFolds = atoi(value);
			return true;
		case 11:
			this->classifierSvm.cvMaxDatasetSize = atol(value);
			return true;
		case 12:
			this->classifierSvm.cvCoarseCMin = atof(value);
			return true;
		case 13:
			this->classifierSvm.cvCoarseCMax = atof(value);
			return true;
		case 14:
			this->classifierSvm.cvCoarseCStep = atof(value);
			return true;
		case 15:
			this->classifierSvm.cvCoarseGammaMin = atof(value);
			return true;
		case 16:
			this->classifierSvm.cvCoarseGammaMax = atof(value);
			return true;
		case 17:
			this->classifierSvm.cvCoarseGammaStep = atof(value);
			return true;
		case 18:
			this->classifierSvm.cvFineCStep = atof(value);
			return true;
		case 19:
			this->classifierSvm.cvFineCRadius = atof(value);
			return true;
		case 20:
			this->classifierSvm.cvFineGammaStep = atof(value);
			return true;
		case 21:
			this->classifierSvm.cvFineGammaRadius = atof(value);
			return true;
		case 22:
			this->classifierSvm.degree = atof(value);
			return true;
		case 23:
			this->classifierSvm.doCV = sclib::atob(value);
			return true;
		case 24:
			this->classifierSvm.eps = atof(value);
			return true;
		case 25:
			this->classifierSvm.gamma = atof(value);
			return true;
		case 26:
			this->classifierSvm.kernel_type = atoi(value);
			return true;
		case 27:
			this->classifierSvm.nu = atof(value);
			return true;
		case 28:
			this->classifierSvm.p = atof(value);
			return true;
		case 29:
			this->classifierSvm.probability = atoi(value);
			return true;
		case 30:
			this->classifierSvm.shrinking = atoi(value);
			return true;
		case 31:
			this->classifierSvm.svm_type = atoi(value);
			return true;
		case 32:
			this->classifierSvm.oneClassGammaSearchMaxIterations = atoi(value);
			return true;
		case 33:
			this->classifierSvm.oneClassGammaSearchRepeats = atoi(value);
			return true;
		case 34:
			if (sclib::pathExists(value) == true) {
				setDebugDir(value);
				return true;
			} else {
				return false; //value is not a valid parameter-value
			}
		case 35:
			this->cluster.mergeMode = atoi(value);
			return true;
		case 36:
			this->debug.useDebugPrefix = sclib::atob(value);
			return true;
		case 37:
			this->debug.debugMode = atol(value);
			return true;
		case 38:
			this->distanceMeasure.BICpenaltyFactor = atof(value);
			return true;
		//tata, one number missed...
		case 40:
			this->speakerClusterer.distanceMeasure = atoi(value);
			return true;
		case 41:
			this->speakerClusterer.globalCriterion = atoi(value);
			return true;
		case 42:
			this->distanceMeasure.mergeMode = atoi(value);
			return true;
		case 43:
			this->speakerClusterer.terminationCriterion = atoi(value);
			return true;
		case 44:
			this->distanceMeasure.WCDpenaltyFactor = atof(value);
			return true;
		case 45:
			this->emd.debugLevel = atoi(value);
			return true;
		case 46:
			this->emd.maxIterations = atol(value);
			return true;
		case 47:
			this->emd.maxSigSize = atol(value);
			return true;
		case 48:
			this->enhancement.doEnhancement = sclib::atob(value);
			return true;
		case 49:
			this->enhancement.minNoiseDuration = atoi(value);
			return true;
		case 50:
			this->enhancement.noiseModelUpdateRate = atoi(value);
			return true;
		case 51:
			if (value != NULL && strcmp(value, "") != 0) {
				MFree_1D(this->enhancement.speechModelFile);
				if (sclib::fileExists(value) == true) {
					MArray_1D(this->enhancement.speechModelFile, strlen(value)+1, char, "SC_TweakableParameters.setByName: enhancement.speechModelFile");
					sprintf(this->enhancement.speechModelFile, "%s", value);
					return true;
				} else {
					return false;
				}
			} else {
				return true;
			}
		case 52:
			this->featureBandPeriodicity.FFTsize = atoi(value);
			return true;
		case 53:
			this->featureBandPeriodicity.frameSize = atoi(value);
			return true;
		case 54:
			this->featureBandPeriodicity.frameStep = atoi(value);
			return true;
		case 55:
			this->featureBandPeriodicity.highCut = atof(value);
			return true;
		case 56:
			this->featureBandPeriodicity.lowCut = atof(value);
			return true;
		case 57:
			this->featureBandPeriodicity.preEmphasizeFactor = atof(value);
			return true;
		case 58:
			this->featureBandPeriodicity.window = atoi(value);
			return true;
		case 59:
			this->featureBrightnessBandwidth.FFTsize = atoi(value);
			return true;
		case 60:
			this->featureBrightnessBandwidth.frameSize = atoi(value);
			return true;
		case 61:
			this->featureBrightnessBandwidth.frameStep = atoi(value);
			return true;
		case 62:
			this->featureBrightnessBandwidth.highCut = atof(value);
			return true;
		case 63:
			this->featureBrightnessBandwidth.lowCut = atof(value);
			return true;
		case 64:
			this->featureBrightnessBandwidth.preEmphasizeFactor = atof(value);
			return true;
		case 65:
			this->featureBrightnessBandwidth.window = atoi(value);
			return true;
		case 66:
			this->featureFbe.addDeltaDeltas = sclib::atob(value);
			return true;
		case 67:
			this->featureFbe.addDeltas = sclib::atob(value);
			return true;
		case 68:
			this->featureFbe.CMN = sclib::atob(value);
			return true;
		case 69:
			this->featureFbe.dEnergy = sclib::atob(value);
			return true;
		case 70:
			this->featureFbe.FFTsize = atoi(value);
			return true;
		case 71:
			this->featureFbe.filterBankSize = atoi(value);
			return true;
		case 72:
			this->featureFbe.frameSize = atoi(value);
			return true;
		case 73:
			this->featureFbe.frameStep = atoi(value);
			return true;
		case 74:
			this->featureFbe.frequencyScale = atoi(value);
			return true;
		case 75:
			this->featureFbe.highCut = atof(value);
			return true;
		case 76:
			this->featureFbe.lowCut = atof(value);
			return true;
		case 77:
			if (value[0] != '\'') {
				return false; //we need this "string" encapsulated in single '-signs!!!
			} else {
				char i = 1;
				while (value[i] != '\'' && i < 9) {
					this->featureFbe.MFCCcoeffSelection[i-1] = value[i];
					i++;
				}
				for (char j = i; j < 8; i++) {
					this->featureFbe.MFCCcoeffSelection[i] = 0x00;
				}
			}
			return true;
		case 78:
			this->featureFbe.MFCCorder = atoi(value);
			return true;
		case 79:
			this->featureFbe.preEmphasizeFactor = atof(value);
			return true;
		case 80:
			this->featureFbe.smoothing = atoi(value);
			return true;
		case 81:
			this->featureFbe.window = atoi(value);
			return true;
		case 82:
			this->featureFbe.minFilterBankFrequency = atof(value);
			return true;
		case 83:
			this->featureFbe.maxFilterBankFrequency = atof(value);
			return true;
		case 84:
			this->featureFbe.resultType = atoi(value);
			return true;
		case 85:
			this->featureMfcc.addDeltaDeltas = sclib::atob(value);
			return true;
		case 86:
			this->featureMfcc.addDeltas = sclib::atob(value);
			return true;
		case 87:
			this->featureMfcc.CMN = sclib::atob(value);
			return true;
		case 88:
			this->featureMfcc.dEnergy = sclib::atob(value);
			return true;
		case 89:
			this->featureMfcc.fftSize = atoi(value);
			return true;
		case 90:
			this->featureMfcc.filterBankSize = atoi(value);
			return true;
		case 91:
			this->featureMfcc.frameSize = atoi(value);
			return true;
		case 92:
			this->featureMfcc.frameStep = atoi(value);
			return true;
		case 93:
			this->featureMfcc.highCut = atof(value);
			return true;
		case 94:
			this->featureMfcc.lowCut = atof(value);
			return true;
		case 95:
			if (value[0] != '\'') {
				return false; //we need this "string" encapsulated in single '-signs!!!
			} else {
				char i = 1;
				while (value[i] != '\'' && i < 9) {
					this->featureMfcc.coeffSelection[i-1] = value[i];
					i++;
				}
				for (char j = i; j < 8; i++) {
					this->featureMfcc.coeffSelection[i] = 0x00;
				}
			}
			return true;
		case 96:
			this->featureMfcc.MFCCorder = atoi(value);
			return true;
		case 97:
			this->featureMfcc.preEmphasizeFactor = atof(value);
			return true;
		case 98:
			this->featureMfcc.window = atoi(value);
			return true;
		case 99:
			this->featureMfcc.method = atoi(value);
			return true;
		case 100:
			this->featureMfcc.sclib_frequencyScale = atoi(value);
			return true;
		case 101:
			this->featureMfcc.sclib_minFilterBankFrequency = atoi(value);
			return true;
		case 102:
			this->featureMfcc.sclib_maxFilterBankFrequency = atoi(value);
			return true;
		case 103:
			this->featureMfcc.sclib_smoothing = atoi(value);
			return true;
		case 104:
			this->featureNfr.FFTsize = atoi(value);
			return true;
		case 105:
			this->featureNfr.frameSize = atoi(value);
			return true;
		case 106:
			this->featureNfr.frameStep = atoi(value);
			return true;
		case 107:
			this->featureNfr.highCut = atof(value);
			return true;
		case 108:
			this->featureNfr.lowCut = atof(value);
			return true;
		case 109:
			this->featureNfr.NFRthreshold = atof(value);
			return true;
		case 110:
			this->featureNfr.preEmphasizeFactor = atof(value);
			return true;
		case 111:
			this->featureNfr.window = atoi(value);
			return true;
		case 112:
			this->featureSpectrum.FFTsize = atoi(value);
			return true;
		case 113:
			this->featureSpectrum.frameSize = atoi(value);
			return true;
		case 114:
			this->featureSpectrum.frameStep = atoi(value);
			return true;
		case 115:
			this->featureSpectrum.highCut = atof(value);
			return true;
		case 116:
			this->featureSpectrum.lowCut = atof(value);
			return true;
		case 117:
			this->featureSpectrum.preEmphasizeFactor = atof(value);
			return true;
		case 118:
			this->featureSpectrum.window = atoi(value);
			return true;
		case 119:
			this->featureSpectrumFlux.FFTsize = atoi(value);
			return true;
		case 120:
			this->featureSpectrumFlux.frameSize = atoi(value);
			return true;
		case 121:
			this->featureSpectrumFlux.frameStep = atoi(value);
			return true;
		case 122:
			this->featureSpectrumFlux.highCut = atof(value);
			return true;
		case 123:
			this->featureSpectrumFlux.lowCut = atof(value);
			return true;
		case 124:
			this->featureSpectrumFlux.preEmphasizeFactor = atof(value);
			return true;
		case 125:
			this->featureSpectrumFlux.window = atoi(value);
			return true;
		case 126:
			this->featureSubBandPower.FFTsize = atoi(value);
			return true;
		case 127:
			this->featureSubBandPower.frameSize = atoi(value);
			return true;
		case 128:
			this->featureSubBandPower.frameStep = atoi(value);
			return true;
		case 129:
			this->featureSubBandPower.highCut = atof(value);
			return true;
		case 130:
			this->featureSubBandPower.lowCut = atof(value);
			return true;
		case 131:
			this->featureSubBandPower.preEmphasizeFactor = atof(value);
			return true;
		case 132:
			this->featureSubBandPower.window = atoi(value);
			return true;
		case 133:
			this->featureZcr.frameSize = atoi(value);
			return true;
		case 134:
			this->featureZcr.frameStep = atoi(value);
			return true;
		case 135:
			this->featureZcr.highCut = atof(value);
			return true;
		case 136:
			this->featureZcr.lowCut = atof(value);
			return true;
		case 137:
			this->featureZcr.preEmphasizeFactor = atof(value);
			return true;
		case 138:
			this->featureZcr.useChebyshev = sclib::atob(value);
			return true;
		case 139:
			this->featureZcr.scaleResult = sclib::atob(value);
			return true;
		case 140:
			this->featureSte.frameSize = atoi(value);
			return true;
		case 141:
			this->featureSte.frameStep = atoi(value);
			return true;
		case 142:
			this->featureSte.highCut = atof(value);
			return true;
		case 143:
			this->featureSte.lowCut = atof(value);
			return true;
		case 144:
			this->featureSte.preEmphasizeFactor = atof(value);
			return true;
		case 145:
			this->featureSte.useButterworth = sclib::atob(value);
			return true;
		case 146:
			this->featureSte.scaleResult = sclib::atob(value);
			return true;
		case 164:
			this->featureLpc.lowCut = atof(value);
			return true;
		case 165:
			this->featureLpc.highCut = atof(value);
			return true;
		case 166:
			this->featureLpc.frameSize = atoi(value);
			return true;
		case 167:
			this->featureLpc.frameStep = atoi(value);
			return true;
		case 168:
			this->featureLpc.preEmphasizeFactor = atof(value);
			return true;
		case 169:
			this->featureLpc.LPCorder = atoi(value);
			return true;
		case 170:
			this->featureLpcResidual.lowCut = atof(value);
			return true;
		case 171:
			this->featureLpcResidual.highCut = atof(value);
			return true;
		case 172:
			this->featureLpcResidual.frameSize = atoi(value);
			return true;
		case 173:
			this->featureLpcResidual.frameStep = atoi(value);
			return true;
		case 174:
			this->featureLpcResidual.preEmphasizeFactor = atof(value);
			return true;
		case 175:
			this->featureLpcResidual.order = atoi(value);
			return true;
		case 176:
			this->featureLsp.lowCut = atof(value);
			return true;
		case 177:
			this->featureLsp.highCut = atof(value);
			return true;
		case 178:
			this->featureLsp.frameSize = atoi(value);
			return true;
		case 179:
			this->featureLsp.frameStep = atoi(value);
			return true;
		case 180:
			this->featureLsp.preEmphasizeFactor = atof(value);
			return true;
		case 181:
			this->featureLsp.LPCorder = atoi(value);
			return true;
		case 182:
			this->featureLsp.minSeparation = atof(value);
			return true;
		case 183:
			this->featureLsp.bisections = atoi(value);
			return true;
		case 184:
			this->featureLsp.delta = atof(value);
			return true;
		case 185:
			this->featureLsp.method = atoi(value);
			return true;
		case 186:
			this->featureLsp.maxLoops = atoi(value);
			return true;
		case 187:
			this->featurePitch.frameSize = atoi(value);
			return true;
		case 188:
			this->featurePitch.frameStep = atoi(value);
			return true;
		case 189:
			this->featurePitch.method = atoi(value);
			return true;
		case 190:
			fValue = (float)(atof(value));
			if ((fValue < 0.01) || (fValue > 0.99)) {
				REPORT_ERROR(SVLIB_BadArg, "featurePitch.esps_cand_thresh parameter must be between [0.01, 0.99]");
				return false;
			} else {
				this->featurePitch.esps_cand_thresh = fValue;
				return true;
			}
		case 191:
			this->featurePitch.esps_lag_weight = (float)(atof(value));
			return true;
		case 192:
			this->featurePitch.esps_freq_weight = (float)(atof(value));
			return true;
		case 193:
			this->featurePitch.esps_trans_cost = (float)(atof(value));
			return true;
		case 194:
			this->featurePitch.esps_trans_amp = (float)(atof(value));
			return true;
		case 195:
			this->featurePitch.esps_trans_spec = (float)(atof(value));
			return true;
		case 196:
			this->featurePitch.esps_voice_bias = (float)(atof(value));
			return true;
		case 197:
			this->featurePitch.esps_double_cost = (float)(atof(value));
			return true;
		case 198:
			this->featurePitch.esps_min_f0 = (float)(atof(value));
			return true;
		case 199:
			this->featurePitch.esps_max_f0 = (float)(atof(value));
			return true;
		case 200:
			iValue = atoi(value);
			if ((iValue > 100) || (iValue < 3)) {
				REPORT_ERROR(SVLIB_BadArg, "featurePitch.esps_n_cands parameter must be between [3,100]");
				return false;
			} else {
				this->featurePitch.esps_n_cands = iValue;
				return true;
			}
		case 201:
			fValue = (float)(atof(value));
			if ((fValue > (float)(0.1)) || (fValue < (float)(0.0001))) {
				REPORT_ERROR(SVLIB_BadArg, "featurePitch.esps_wind_dur parameter must be between [0.0001, 0.1]");
				return false;
			} else {
				this->featurePitch.esps_wind_dur = fValue;
				return true;
			}
		case 202:
			this->featureSdp.color = atoi(value);
			return true;
		case 203:
			this->featureSdp.frameSize = atoi(value);
			return true;
		case 204:
			this->featureSdp.frameStep = atoi(value);
			return true;
		case 205:
			this->featureSdp.lag = atoi(value);
			return true;
		case 206:
			this->featureSdp.m = atoi(value);
			return true;
		case 207:
			this->featureSdp.n = atoi(value);
			return true;
		case 208:
			this->featureSdp.pictureSize = atoi(value);
			return true;
		case 209:
			this->featureSdp.tau = atoi(value);
			return true;
		case 210:
			this->featureSdp.preEmphasizeFactor = atof(value);
			return true;
		case 211:
			this->general.firstScene = atoi(value);
			return true;
		case 212:
			this->general.lastScene = atoi(value);
			return true;
		case 213:
			this->general.sceneSelection = atol(value);
			return true;
		case 214:
			this->general.shortSpeechThreshold = atol(value);
			return true;
		case 215:
			this->general.pauseSilenceThreshold = atoi(value);
			return true;
		case 216:
			if (value != NULL && strcmp(value, "") != 0) {
				MFree_1D(this->general.preClusteringResultsPrefix);
				MArray_1D(this->general.preClusteringResultsPrefix, strlen(value)+1, char, "SC_TweakableParameters.setByName: general.preClusteringResultsPrefix");
				sprintf(this->general.preClusteringResultsPrefix, "%s", value);
			}
			return true;
		case 217:
			if (value != NULL && strcmp(value, "") != 0) {
				MFree_1D(this->general.featurePrefix);
				MArray_1D(this->general.featurePrefix, strlen(value)+1, char, "SC_TweakableParameters.setByName: general.featurePrefix");
				sprintf(this->general.featurePrefix, "%s", value);
			}
			return true;
		case 218:
			this->groundTruth.internalFrameSize = atoi(value);
			return true;
		case 219:
			this->groundTruth.pseudoSceneLength = atol(value);
			return true;
		case 220:
			this->groundTruth.videoFrameMachineOffset = atoi(value);
			return true;
		case 221:
			this->groundTruth.storeProbabilityInformation = sclib::atob(value);
			return true;
		case 222:
			this->mixtureModelBgmm.EMthreshold = atof(value);
			return true;
		case 223:
			this->mixtureModelBgmm.maxEMiterations = atol(value);
			return true;
		case 224:
			this->mixtureModelBgmm.varianceLimit = atof(value);
			return true;
		case 225:
			this->mixtureModelBgmm.weightLimit = atof(value);
			return true;
		case 226:
			this->mixtureModelGmm.EMthreshold = atof(value);
			return true;
		case 227:
			this->mixtureModelGmm.maxEMiterations = atol(value);
			return true;
		case 228:
			this->mixtureModelGmm.kMeansIterations = atoi(value);
			return true;
		case 229:
			this->mixtureModelGmm.varianceLimit = atof(value);
			return true;
		case 230:
			this->mixtureModelGmm.weightLimit = atof(value);
			return true;
		case 231:
			this->mixtureModelMixMax.bgModelCombination = sclib::atob(value);
			return true;
		case 232:
			this->mixtureModelMixMax.EMthreshold = atof(value);
			return true;
		case 233:
			this->mixtureModelMixMax.kMeansIterations = atoi(value);
			return true;
		case 234:
			this->mixtureModelMixMax.maxEMiterations = atol(value);
			return true;
		case 235:
			this->mixtureModelMixMax.maxERFloops = atol(value);
			return true;
		case 236:
			this->mixtureModelMixMax.noiseCorruptionType = (unsigned char)atoi(value);
			return true;
		case 237:
			this->mixtureModelMixMax.varianceLimit = atof(value);
			return true;
		case 238:
			this->mixtureModelMixMax.weightLimit = atof(value);
			return true;
		case 239:
			this->mixtureModelGmmubm.adaptMeans = sclib::atob(value);
			return true;
		case 240:
			this->mixtureModelGmmubm.adaptVariances = sclib::atob(value);
			return true;
		case 241:
			this->mixtureModelGmmubm.adaptWeights = sclib::atob(value);
			return true;
		case 242:
			this->mixtureModelGmmubm.relevanceFactor = atof(value);
			return true;
		case 243:
			this->mixtureModelGmmubm.scoringMethod = atoi(value);
			return true;
		case 244:
			this->mixtureModelGmmubm.topCmixtures = atoi(value);
			return true;
		case 245:
			if (value != NULL && strcmp(value, "") != 0) {
				MFree_1D(this->mixtureModelGmmubm.ubmFileName);
				if (sclib::fileExists(value) == true) {
					MArray_1D(this->mixtureModelGmmubm.ubmFileName, strlen(value)+1, char, "SC_TweakableParameters.setByName: mixtureModelGmmubm.ubmFileName");
					sprintf(this->mixtureModelGmmubm.ubmFileName, "%s", value);
					return true;
				} else {
					return false;
				}
			} else {
				return true;
			}
		case 246:
			this->mixtureModelGmmubm.varianceLimit = atof(value);
			return true;
		case 247:
			this->mixtureModelMix2Max.EMthreshold = atof(value);
			return true;
		case 248:
			this->mixtureModelMix2Max.kMeansIterations = atoi(value);
			return true;
		case 249:
			this->mixtureModelMix2Max.maxEMiterations = atol(value);
			return true;
		case 250:
			this->mixtureModelMix2Max.maxERFloops = atol(value);
			return true;
		case 251:
			this->mixtureModelMix2Max.varianceLimit = atof(value);
			return true;
		case 252:
			this->mixtureModelMix2Max.weightLimit = atof(value);
			return true;
		case 253:
			this->mixtureModelMix2Max.bgModelCombination = sclib::atob(value);
			return true;
		case 254:
			this->mixtureModelMix2MaxEx.bgModelCombination = sclib::atob(value);
			return true;
		case 255:
			this->mixtureModelMix2MaxEx.EMthreshold = atof(value);
			return true;
		case 256:
			this->mixtureModelMix2MaxEx.kMeansIterations = atoi(value);
			return true;
		case 257:
			this->mixtureModelMix2MaxEx.maxEMiterations = atol(value);
			return true;
		case 258:
			this->mixtureModelMix2MaxEx.maxERFloops = atol(value);
			return true;
		case 259:
			this->mixtureModelMix2MaxEx.varianceLimit = atof(value);
			return true;
		case 260:
			this->mixtureModelMix2MaxEx.weightLimit = atof(value);
			return true;
		case 261:
			this->modelHandler.backgroundModelType = atol(value);
			return true;
		case 269:
			this->modelHandler.foregroundModelType = atol(value);
			return true;
		case 270:
			this->modelHandler.maxNoiseModelOrder = atoi(value);
			return true;
		case 271:
			this->modelHandler.maxSpeakerModelOrder = atoi(value);
			return true;
		case 272:
			this->modelHandler.msPerGaussian = atol(value);
			return true;
		case 273:
			this->modelHandler.orderGuessEMsteps = atol(value);
			return true;
		case 274:
			this->modelHandler.orderGuessMode = atoi(value);
			return true;
		case 275:
			this->modelHandler.SNRthreshold = atof(value);
			return true;
		case 276:
			if (value != NULL && strcmp(value, "") != 0) {
				MFree_1D(this->modelHandler.onlyThisSpeaker);
				MArray_1D(this->modelHandler.onlyThisSpeaker, strlen(value)+1, char, "SC_TweakableParameters.setByName: modelHandler.onlyThisSpeaker");
				sprintf(this->modelHandler.onlyThisSpeaker, "%s", value);
			}
			return true;
		case 277:
			this->modelHandler.speakerModelFeature = atoi(value);
			return true;
		case 278:
			this->modelPareto.useMarginalDistributions = sclib::atob(value);
			return true;
		case 279:
			this->modelQgmm.deltaBIClambda = atof(value);
			return true;
		case 280:
			this->modelQgmm.maxMixtures = atoi(value);
			return true;
		case 281:
			this->modelQgmm.percentDifference = atof(value);
			return true;
		case 282:
			this->resampling.fastConversion = sclib::atob(value);
			return true;
		case 283:
			this->score.BBNmetricLambda = atof(value);
			return true;
		case 284:
			if (value != NULL && strcmp(value, "") != 0) {
				MFree_1D(this->segmentationAudioTypeLzl.actionModelFileName);
				MArray_1D(this->segmentationAudioTypeLzl.actionModelFileName, strlen(value)+1, char, "SC_TweakableParameters.setByName: segmentationAudioTypeLzl.actionModelFileName");
				sprintf(this->segmentationAudioTypeLzl.actionModelFileName, "%s", value);
				return true;
			} else {
				return true;
			}
		case 285:
			if (value != NULL && strcmp(value, "") != 0) {
				MFree_1D(this->segmentationAudioTypeLzl.classifierFileName);
				//if (sclib::fileExists(value) == true) {
					MArray_1D(this->segmentationAudioTypeLzl.classifierFileName, strlen(value)+1, char, "SC_TweakableParameters.setByName: segmentationAudioTypeLzl.classifierFileName");
					sprintf(this->segmentationAudioTypeLzl.classifierFileName, "%s", value);
					return true;
				//} else {
				//  return false;
				//}
			} else {
				return true;
			}
		case 286:
			if (value != NULL && strcmp(value, "") != 0) {
				MFree_1D(this->segmentationAudioTypeLzl.featureFileName);
				//if (sclib::fileExists(value) == true) {
					MArray_1D(this->segmentationAudioTypeLzl.featureFileName, strlen(value)+1, char, "SC_TweakableParameters.setByName: segmentationAudioTypeLzl.featureFileName");
					sprintf(this->segmentationAudioTypeLzl.featureFileName, "%s", value);
					return true;
				//} else {
				//  return false;
				//}
			} else {
				return true;
			}
		case 287:
			if (value != NULL && strcmp(value, "") != 0) {
				MFree_1D(this->segmentationAudioTypeLzl.normalizationFileName);
				//if (sclib::fileExists(value) == true) {
					MArray_1D(this->segmentationAudioTypeLzl.normalizationFileName, strlen(value)+1, char, "SC_TweakableParameters.setByName: segmentationAudioTypeLzl.normalizationFileName");
					sprintf(this->segmentationAudioTypeLzl.normalizationFileName, "%s", value);
					return true;
				//} else {
				//  return false;
				//}
			} else {
				return true;
			}
		case 288:
			this->segmentationAudioTypeLzl.subClipLength = atol(value);
			return true;
		case 289:
			this->segmentationAudioTypeLzl.bbParameters.FFTsize = atoi(value);
			return true;
		case 290:
			this->segmentationAudioTypeLzl.bbParameters.frameSize = atoi(value);
			return true;
		case 291:
			this->segmentationAudioTypeLzl.bbParameters.frameStep = atoi(value);
			return true;
		case 292:
			this->segmentationAudioTypeLzl.bbParameters.preEmphasizeFactor = atof(value);
			return true;
		case 293:
			this->segmentationAudioTypeLzl.bbParameters.window = atoi(value);
			return true;
		case 294:
			this->segmentationAudioTypeLzl.bbParameters.sampleRate = atof(value);
			return true;
		case 295:
			this->segmentationAudioTypeLzl.bpParameters.FFTsize = atoi(value);
			return true;
		case 296:
			this->segmentationAudioTypeLzl.bpParameters.frameSize = atoi(value);
			return true;
		case 297:
			this->segmentationAudioTypeLzl.bpParameters.frameStep = atoi(value);
			return true;
		case 298:
			this->segmentationAudioTypeLzl.bpParameters.preEmphasizeFactor = atof(value);
			return true;
		case 299:
			this->segmentationAudioTypeLzl.bpParameters.window = atoi(value);
			return true;
		case 300:
			this->segmentationAudioTypeLzl.bpParameters.sampleRate = atof(value);
			return true;
		case 301:
			this->segmentationAudioTypeLzl.mfccParameters.fftSize = atoi(value);
			return true;
		case 302:
			this->segmentationAudioTypeLzl.mfccParameters.frameSize = atoi(value);
			return true;
		case 303:
			this->segmentationAudioTypeLzl.mfccParameters.frameStep = atoi(value);
			return true;
		case 304:
			this->segmentationAudioTypeLzl.mfccParameters.MFCCorder = atoi(value);
			return true;
		case 305:
			this->segmentationAudioTypeLzl.mfccParameters.addDeltaDeltas = sclib::atob(value);
			return true;
		case 306:
			this->segmentationAudioTypeLzl.mfccParameters.addDeltas = sclib::atob(value);
			return true;
		case 307:
			this->segmentationAudioTypeLzl.mfccParameters.method = atoi(value);
			return true;
		case 308:
			this->segmentationAudioTypeLzl.mfccParameters.dEnergy = sclib::atob(value);
			return true;
		case 309:
			this->segmentationAudioTypeLzl.mfccParameters.CMN = sclib::atob(value);
			return true;
		case 310:
			this->segmentationAudioTypeLzl.mfccParameters.filterBankSize = atoi(value);
			return true;
		case 311:
			this->segmentationAudioTypeLzl.mfccParameters.preEmphasizeFactor = atof(value);
			return true;
		case 312:
			this->segmentationAudioTypeLzl.mfccParameters.sclib_frequencyScale = atoi(value);
			return true;
		case 313:
			this->segmentationAudioTypeLzl.mfccParameters.sclib_smoothing = atoi(value);
			return true;
		case 314:
			this->segmentationAudioTypeLzl.mfccParameters.sclib_maxFilterBankFrequency = atof(value);
			return true;
		case 315:
			this->segmentationAudioTypeLzl.mfccParameters.sclib_minFilterBankFrequency = atof(value);
			return true;
		case 316:
			this->segmentationAudioTypeLzl.mfccParameters.window = atoi(value);
			return true;
		case 317:
			if (value[0] != '\'') {
				return false; //we need this "string" encapsulated in single '-signs!!!
			} else {
				char i = 1;
				while (value[i] != '\'' && i < 9) {
					this->segmentationAudioTypeLzl.mfccParameters.coeffSelection[i-1] = value[i];
					i++;
				}
				for (char j = i; j < 8; i++) {
					this->segmentationAudioTypeLzl.mfccParameters.coeffSelection[i] = 0x00;
				}
			}
			return true;
		//case 318:
			//this->segmentationAudioTypeLzl.mfccParameters.coeffSelection = atol(value);
			//return true;
		case 319:
			this->segmentationAudioTypeLzl.mfccParameters.sampleRate = atof(value);
			return true;
		case 320:
			this->segmentationAudioTypeLzl.nfrParameters.FFTsize = atoi(value);
			return true;
		case 321:
			this->segmentationAudioTypeLzl.nfrParameters.frameSize = atoi(value);
			return true;
		case 322:
			this->segmentationAudioTypeLzl.nfrParameters.frameStep = atoi(value);
			return true;
		case 323:
			this->segmentationAudioTypeLzl.nfrParameters.preEmphasizeFactor = atof(value);
			return true;
		case 324:
			this->segmentationAudioTypeLzl.nfrParameters.window = atoi(value);
			return true;
		case 325:
			this->segmentationAudioTypeLzl.nfrParameters.NFRthreshold = atof(value);
			return true;
		case 326:
			this->segmentationAudioTypeLzl.nfrParameters.sampleRate = atof(value);
			return true;
		case 327:
			this->segmentationAudioTypeLzl.sbpParameters.FFTsize = atoi(value);
			return true;
		case 328:
			this->segmentationAudioTypeLzl.sbpParameters.frameSize = atoi(value);
			return true;
		case 329:
			this->segmentationAudioTypeLzl.sbpParameters.frameStep = atoi(value);
			return true;
		case 330:
			this->segmentationAudioTypeLzl.sbpParameters.preEmphasizeFactor = atof(value);
			return true;
		case 331:
			this->segmentationAudioTypeLzl.sbpParameters.window = atoi(value);
			return true;
		case 332:
			this->segmentationAudioTypeLzl.sbpParameters.sampleRate = atof(value);
			return true;
		case 333:
			this->segmentationAudioTypeLzl.sfParameters.FFTsize = atoi(value);
			return true;
		case 334:
			this->segmentationAudioTypeLzl.sfParameters.frameSize = atoi(value);
			return true;
		case 335:
			this->segmentationAudioTypeLzl.sfParameters.frameStep = atoi(value);
			return true;
		case 336:
			this->segmentationAudioTypeLzl.sfParameters.preEmphasizeFactor = atof(value);
			return true;
		case 337:
			this->segmentationAudioTypeLzl.sfParameters.window = atoi(value);
			return true;
		case 338:
			this->segmentationAudioTypeLzl.sfParameters.sampleRate = atof(value);
			return true;
		case 339:
			this->segmentationAudioTypeLzl.zcrParameters.useChebyshev = sclib::atob(value);
			return true;
		case 340:
			this->segmentationAudioTypeLzl.zcrParameters.frameSize = atoi(value);
			return true;
		case 341:
			this->segmentationAudioTypeLzl.zcrParameters.frameStep = atoi(value);
			return true;
		case 342:
			this->segmentationAudioTypeLzl.zcrParameters.preEmphasizeFactor = atof(value);
			return true;
		case 343:
			this->segmentationAudioTypeLzl.zcrParameters.scaleResult = sclib::atob(value);
			return true;
		case 344:
			this->segmentationAudioTypeLzl.zcrParameters.sampleRate = atof(value);
			return true;
		case 345:
			this->segmentationChangesLz.adaptiveThresholdAlpha = atof(value);
			return true;
		case 346:
			this->segmentationChangesLz.bayesianThreshold = atof(value);
			return true;
		case 347:
			this->segmentationChangesLz.detectorWindowLength = atol(value);
			return true;
		case 348:
			this->segmentationChangesLz.detectorWindowStep = atol(value);
			return true;
		case 349:
			this->segmentationChangesLz.lastNdistances = atoi(value);
			return true;
		case 350:
			if (value != NULL && strcmp(value, "") != 0) {
				MFree_1D(this->segmentationChangesLz.priorsFileName);
				MArray_1D(this->segmentationChangesLz.priorsFileName, strlen(value)+1, char, "SC_TweakableParameters.setByName: segmentationChangesLz.priorsFileName");
				sprintf(this->segmentationChangesLz.priorsFileName, "%s", value);
			}
			return true;
		case 351:
			if (value != NULL && strcmp(value, "") != 0) {
				MFree_1D(this->segmentationChangesLz.time2ChangeModelFileName);
				MArray_1D(this->segmentationChangesLz.time2ChangeModelFileName, strlen(value)+1, char, "SC_TweakableParameters.setByName: segmentationChangesLz.time2ChangeModelFileName");
				sprintf(this->segmentationChangesLz.time2ChangeModelFileName, "%s", value);
			}
			return true;
		case 377:
			this->segmentationHandler.audioTypeMode = atoi(value);
			return true;
		case 378:
			this->segmentationHandler.changeDetectorMode = atoi(value);
			return true;
		case 379:
			this->segmentationHandler.silenceDetectorMode = atoi(value);
			return true;
		case 380:
			this->segmentationHandler.vUvDetectorMode = atoi(value);
			return true;
		case 381:
			this->segmentationSilenceLnk.energyQuantizationLevel = atoi(value);
			return true;
		case 382:
			this->segmentationSilenceLzl.energySilenceThreshold = atof(value);
			return true;
		case 383:
			this->segmentationSilenceLzl.zcrSilenceThreshold = atof(value);
			return true;
		case 384:
			this->segmentationVUvLnk.energyQuantizationLevel = atoi(value);
			return true;
		case 412:
			this->clusterer.maxIterations = atoi(value);
			return true;
		case 413:
			this->clusterer.numCluster = atoi(value);
			return true;
		case 414:
			this->clusterer.iterationsToIniKMeanList = atoi(value);
			return true;
		case 415:
			this->signalHandler.forceSampleRate = atoi(value);
			return true;
		case 416:
			this->signalMpeg.fastSeeking = sclib::atob(value);
			return true;
		case 417:
			this->signalMpeg.hqResampling = sclib::atob(value);
			return true;
		case 418:
			this->signalMpeg.outputChannelCount = atoi(value);
			return true;
		case 419:
			this->signalMpeg.outputSampleRate = atoi(value);
			return true;
		case 420:
			this->speakerClusterer.doClustering = sclib::atob(value);
			return true;
		case 421:
			this->speakerClusterer.speechSegLengthThreshold = atol(value);
			return true;
		case 422:
			this->speakerIdentification.dNormSampleCount = atol(value);
			return true;
		case 423:
			this->speakerIdentification.doIdentification = sclib::atob(value);
			return true;
		case 424:
			this->speakerIdentification.normalizationMode = atoi(value);
			return true;
		case 425:
			this->speakerIdentification.useUBMs = sclib::atob(value);
			return true;
		case 426:
			dValue = atof(value);

			if (dValue >= 0.0 && dValue <= 1.0) {
				this->transform.taperingLength = dValue;
				return true;
			} else {
				return false; //parameter value out of range
			}
		case 427:
			this->segmentationAudioTypeLzl.svmParameters.C = atof(value);
			return true;
		case 428:
			this->segmentationAudioTypeLzl.svmParameters.cache_size = atof(value);
			return true;
		case 429:
			this->segmentationAudioTypeLzl.svmParameters.coef0 = atof(value);
			return true;
		case 430:
			this->segmentationAudioTypeLzl.svmParameters.cvFolds = atoi(value);
			return true;
		case 431:
			this->segmentationAudioTypeLzl.svmParameters.cvMaxDatasetSize = atol(value);
			return true;
		case 432:
			this->segmentationAudioTypeLzl.svmParameters.cvCoarseCMin = atof(value);
			return true;
		case 433:
			this->segmentationAudioTypeLzl.svmParameters.cvCoarseCMax = atof(value);
			return true;
		case 434:
			this->segmentationAudioTypeLzl.svmParameters.cvCoarseCStep = atof(value);
			return true;
		case 435:
			this->segmentationAudioTypeLzl.svmParameters.cvCoarseGammaMin = atof(value);
			return true;
		case 436:
			this->segmentationAudioTypeLzl.svmParameters.cvCoarseGammaMax = atof(value);
			return true;
		case 437:
			this->segmentationAudioTypeLzl.svmParameters.cvCoarseGammaStep = atof(value);
			return true;
		case 438:
			this->segmentationAudioTypeLzl.svmParameters.cvFineCStep = atof(value);
			return true;
		case 439:
			this->segmentationAudioTypeLzl.svmParameters.cvFineCRadius = atof(value);
			return true;
		case 440:
			this->segmentationAudioTypeLzl.svmParameters.cvFineGammaStep = atof(value);
			return true;
		case 441:
			this->segmentationAudioTypeLzl.svmParameters.cvFineGammaRadius = atof(value);
			return true;
		case 442:
			this->segmentationAudioTypeLzl.svmParameters.degree = atof(value);
			return true;
		case 443:
			this->segmentationAudioTypeLzl.svmParameters.doCV = sclib::atob(value);
			return true;
		case 444:
			this->segmentationAudioTypeLzl.svmParameters.eps = atof(value);
			return true;
		case 445:
			this->segmentationAudioTypeLzl.svmParameters.gamma = atof(value);
			return true;
		case 446:
			this->segmentationAudioTypeLzl.svmParameters.kernel_type = atoi(value);
			return true;
		case 447:
			this->segmentationAudioTypeLzl.svmParameters.nu = atof(value);
			return true;
		case 448:
			this->segmentationAudioTypeLzl.svmParameters.p = atof(value);
			return true;
		case 449:
			this->segmentationAudioTypeLzl.svmParameters.probability = atoi(value);
			return true;
		case 450:
			this->segmentationAudioTypeLzl.svmParameters.shrinking = atoi(value);
			return true;
		case 451:
			this->segmentationAudioTypeLzl.svmParameters.svm_type = atoi(value);
			return true;
		case 452:
			this->segmentationAudioTypeLzl.svmParameters.oneClassGammaSearchMaxIterations = atoi(value);
			return true;
		case 453:
			this->segmentationAudioTypeLzl.svmParameters.oneClassGammaSearchRepeats = atoi(value);
			return true;
		case 454:
			this->segmentationSilenceLzl.zcrParameters.frameSize = atoi(value);
			return true;
		case 455:
			this->segmentationSilenceLzl.zcrParameters.frameStep = atoi(value);
			return true;
		case 456:
			this->segmentationSilenceLzl.zcrParameters.sampleRate = atof(value);
			return true;
		case 457:
			this->segmentationSilenceLzl.zcrParameters.preEmphasizeFactor = atof(value);
			return true;
		case 458:
			this->segmentationSilenceLzl.zcrParameters.useChebyshev = sclib::atob(value);
			return true;
		case 459:
			this->segmentationSilenceLzl.zcrParameters.scaleResult = sclib::atob(value);
			return true;
		case 460:
			this->segmentationSilenceLzl.steParameters.frameSize = atoi(value);
			return true;
		case 461:
			this->segmentationSilenceLzl.steParameters.frameStep = atoi(value);
			return true;
		case 462:
			this->segmentationSilenceLzl.steParameters.sampleRate = atof(value);
			return true;
		case 463:
			this->segmentationSilenceLzl.steParameters.preEmphasizeFactor = atof(value);
			return true;
		case 464:
			this->segmentationSilenceLzl.steParameters.useButterworth = sclib::atob(value);
			return true;
		case 465:
			this->segmentationSilenceLzl.steParameters.scaleResult = sclib::atob(value);
			return true;
		case 466:
			this->segmentationVUvEsps.pitchParameters.frameSize = atoi(value);
			return true;
		case 467:
			this->segmentationVUvEsps.pitchParameters.frameSize = atoi(value);
			return true;
		case 468:
			this->segmentationVUvEsps.pitchParameters.method = atoi(value);
			return true;
		case 469:
			this->segmentationVUvEsps.pitchParameters.esps_cand_thresh = (float)atof(value);
			return true;
		case 470:
			this->segmentationVUvEsps.pitchParameters.esps_lag_weight = (float)atof(value);
			return true;
		case 471:
			this->segmentationVUvEsps.pitchParameters.esps_freq_weight = (float)atof(value);
			return true;
		case 472:
			this->segmentationVUvEsps.pitchParameters.esps_trans_cost = (float)atof(value);
			return true;
		case 473:
			this->segmentationVUvEsps.pitchParameters.esps_trans_amp = (float)atof(value);
			return true;
		case 474:
			this->segmentationVUvEsps.pitchParameters.esps_trans_spec = (float)atof(value);
			return true;
		case 475:
			this->segmentationVUvEsps.pitchParameters.esps_voice_bias = (float)atof(value);
			return true;
		case 476:
			this->segmentationVUvEsps.pitchParameters.esps_double_cost = (float)atof(value);
			return true;
		case 477:
			this->segmentationVUvEsps.pitchParameters.esps_min_f0 = (float)atof(value);
			return true;
		case 478:
			this->segmentationVUvEsps.pitchParameters.esps_max_f0 = (float)atof(value);
			return true;
		case 479:
			this->segmentationVUvEsps.pitchParameters.esps_n_cands = atoi(value);
			return true;
		case 480:
			this->segmentationVUvEsps.pitchParameters.esps_wind_dur = (float)atof(value);
			return true;
		//case 481:
		//	
		//	return true;
		case 482:
			this->featureSpectrum.logarithmize = sclib::atob(value);
			return true;
		case 483:
			this->featureSpectrum.createPhase = sclib::atob(value);
			return true;
		case 484:
			this->modelHmm.stateCount = atoi(value);
			return true;
		case 485:
			if (value != NULL && strcmp(value, "") != 0) {
				MFree_1D(this->modelHmm.transitionStructure);
				MArray_1D(this->modelHmm.transitionStructure, strlen(value)+1, char, "SC_TweakableParameters.setByName: modelHmm.transitionStructure");
				sprintf(this->modelHmm.transitionStructure, "%s", value);
				return true;
			} else {
				MFree_1D(this->modelHmm.transitionStructure);
				return true;
			}
		case 486:
			this->modelHmm.mixturesPerState = atoi(value);
			return true;
		case 487:
			this->modelHmm.useOrthogonalTransform = sclib::atob(value);
			return true;
		case 488:
			this->modelHmm.leftToRight = sclib::atob(value);
			return true;
		case 489:
			this->modelHmm.maxIterations = atoi(value);
			return true;
		case 490:
			this->modelHmm.verbose = sclib::atob(value);
			return true;
		case 491:
			this->mixtureModelBgmm.fullCovariance = sclib::atob(value);
			return true;
		case 492:
			this->mixtureModelBgmm.randomInitialization = sclib::atob(value);
			return true;
		case 493:
			this->featureLpc.window = atoi(value);
			return true;
		case 494:
			this->featureLsp.window = atoi(value);
			return true;
		case 495:
			this->featureFormant.frameSize = atoi(value);
			return true;
		case 496:
			this->featureFormant.frameStep = atoi(value);
			return true;
		case 497:
			this->featureFormant.lowCut = atof(value);
			return true;
		case 498:
			this->featureFormant.highCut = atof(value);
			return true;
		case 499:
			this->featureFormant.preEmphasizeFactor = atof(value);
			return true;
		case 500:
			this->featureFormant.sampleRate = atof(value);
			return true;
		case 501:
			this->featureFormant.esps_lpc_ord = atoi(value);
			return true;
		case 502:
			this->featureFormant.esps_lpc_type = atoi(value);
			return true;
		case 503:
			this->featureFormant.esps_w_type = atoi(value);
			return true;
		case 504:
			this->featureFormant.esps_ds_freq = atof(value);
			return true;
		case 505:
			this->featureFormant.esps_wdur = atof(value);
			return true;
		case 506:
			this->featureFormant.esps_nom_f1 = atof(value);
			return true;
		case 507:
			this->featureFormant.esps_cor_wdur = atof(value);
			return true;
		case 508:
			this->featureFormant.esps_frame_int = atof(value);
			return true;
		case 509:
			this->featureFormant.esps_nform = atoi(value);
			return true;
		case 510:
			this->featureLpc.computeGain = sclib::atob(value);
			return true;
		case 529:
			this->modelVq.codebookSize = atoi(value);
			return true;
		case 530:
			this->modelVq.splitMethod = atoi(value);
			return true;
		case 531:
			this->modelVq.maxIterations = atoi(value);
			return true;
		case 534:
			this->segmentationChangesKbk.r = atof(value);
			return true;
		case 535:
			//this->segmentationChangesKbk.resolutionFactor = atof(value);
			return true;
		case 536:
			this->segmentationChangesKbk.lambda = atof(value);
			return true;
		case 537:
			this->segmentationChangesKbk.mfccParameters.addDeltaDeltas = sclib::atob(value);
			return true;
		case 538:
			this->segmentationChangesKbk.mfccParameters.addDeltas = sclib::atob(value);
			return true;
		case 539:
			this->segmentationChangesKbk.mfccParameters.CMN = sclib::atob(value);
			return true;
		case 540:
			this->segmentationChangesKbk.mfccParameters.dEnergy = sclib::atob(value);
			return true;
		case 541:
			this->segmentationChangesKbk.mfccParameters.fftSize = atoi(value);
			return true;
		case 542:
			this->segmentationChangesKbk.mfccParameters.filterBankSize = atoi(value);
			return true;
		case 543:
			this->segmentationChangesKbk.mfccParameters.frameSize = atoi(value);
			return true;
		case 544:
			this->segmentationChangesKbk.mfccParameters.frameStep = atoi(value);
			return true;
		case 545:
			this->segmentationChangesKbk.mfccParameters.highCut = atof(value);
			return true;
		case 546:
			this->segmentationChangesKbk.mfccParameters.lowCut = atof(value);
			return true;
		case 547:
			if (value[0] != '\'') {
				return false; //we need this "string" encapsulated in single '-signs!!!
			} else {
				char i = 1;
				while (value[i] != '\'' && i < 9) {
					this->segmentationChangesKbk.mfccParameters.coeffSelection[i-1] = value[i];
					i++;
				}
				for (char j = i; j < 8; i++) {
					this->segmentationChangesKbk.mfccParameters.coeffSelection[i] = 0x00;
				}
			}
			return true;
		case 548:
			this->segmentationChangesKbk.mfccParameters.MFCCorder = atoi(value);
			return true;
		case 549:
			this->segmentationChangesKbk.mfccParameters.preEmphasizeFactor = atof(value);
			return true;
		case 550:
			this->segmentationChangesKbk.mfccParameters.window = atoi(value);
			return true;
		case 551:
			this->segmentationChangesKbk.mfccParameters.method = atoi(value);
			return true;
		case 552:
			this->segmentationChangesKbk.mfccParameters.sclib_frequencyScale = atoi(value);
			return true;
		case 553:
			this->segmentationChangesKbk.mfccParameters.sclib_minFilterBankFrequency = atof(value);
			return true;
		case 554:
			this->segmentationChangesKbk.mfccParameters.sclib_maxFilterBankFrequency = atof(value);
			return true;
		case 555:
			this->segmentationChangesKbk.mfccParameters.sclib_smoothing = atoi(value);
			return true;
		case 556:
			this->segmentationChangesKbk.tolerance = atof(value);
			return true;
		case 557:
			this->speakerClusterer.constructNonOverlappingClusters = sclib::atob(value);
			return true;
		case 558:
			this->distanceMeasure.groundDistance = atoi(value);
			return true;
		case 559:
			this->distanceMeasure.ICRthreshold = atof(value);
			return true;
		case 560:
			this->speakerClusterer.linkageMode = atoi(value);
			return true;
		case 561:
			this->modelHandler.outlierRemovalMode = atoi(value);
			return true;
		case 562:
			this->featureSamples.frameSize = atoi(value);
			return true;
		case 563:
			this->featureSamples.frameStep = atoi(value);
			return true;
		case 564:
			this->featureSamples.highCut = atof(value);
			return true;
		case 565:
			this->featureSamples.lowCut = atof(value);
			return true;
		case 566:
			this->modelSvm.distanceBasedTesting = sclib::atob(value);
			return true;
		case 567:
			this->modelSvm.doParameterSearch = sclib::atob(value);
			return true;
		case 568:
			this->modelTime.syllableLength = atoi(value);
			return true;
		case 569:
			this->modelTime.trajectoryStep = atoi(value);
			return true;
		case 570:
			this->modelTime.subModelType = atoi(value);
			return true;
		case 571:
			this->modelTime.removeTiming = sclib::atob(value);
			return true;
		case 572:
			this->modelTime.templateCount = atoi(value);
			return true;
		case 573:
			this->modelTime.clusteringIterations = atoi(value);
			return true;
		case 574:
			this->modelTime.replaceTrainingData = sclib::atob(value);
			return true;
		case 575:
			this->modelTime.checkForTrajectorization = sclib::atob(value);
			return true;
		case 576:
			if (value != NULL && strcmp(value, "") != 0) {
				MFree_1D(this->modelTime.worldModelFile);
				MArray_1D(this->modelTime.worldModelFile, strlen(value)+1, char, "SC_TweakableParameters.setByName: modelTime.worldModelFile");
				sprintf(this->modelTime.worldModelFile, "%s", value);
				return true;
			} else {
				MFree_1D(this->modelTime.worldModelFile);
				return true;
			}
			return true;
		case 577:
			if (value != NULL && strcmp(value, "") != 0) {
				MFree_1D(this->modelTime.normalizationFile);
				MArray_1D(this->modelTime.normalizationFile, strlen(value)+1, char, "SC_TweakableParameters.setByName: modelTime.normalizationFile");
				sprintf(this->modelTime.normalizationFile, "%s", value);
				return true;
			} else {
				MFree_1D(this->modelTime.normalizationFile);
				return true;
			}
			return true;
		case 578:
			this->featurePitch.sqrt = sclib::atob(value);
			return true;
		case 579:
			this->segmentationSilenceLzl.specificity = sclib::getBetween(0.0, atof(value), 1.0);
			return true;
		case 580:
			this->segmentationAudioTypeLzl.speechSpecificity = sclib::getBetween(0.0, atof(value), 1.0);
			return true;
		case 581:
			this->segmentationAudioTypeLzl.musicSpecificity = sclib::getBetween(0.0, atof(value), 1.0);
			return true;
		case 582:
			if (atof(value) > 1.0) { //if the clustering-specificity is given >1, it is regarded as the wanted final number of speakers ;-)
				this->speakerClusterer.specificity = floor(atof(value));
			} else {
				this->speakerClusterer.specificity = sclib::getBetween(0.0, atof(value), 1.0);
			}
			return true;
		case 583:
			if (value != NULL && strcmp(value, "") != 0) {
				MFree_1D(this->speakerClusterer.firstDistanceMatrixPrefix);
				MArray_1D(this->speakerClusterer.firstDistanceMatrixPrefix, strlen(value)+1, char, "SC_TweakableParameters.setByName: speakerClusterer.firstDistanceMatrixPrefix");
				sprintf(this->speakerClusterer.firstDistanceMatrixPrefix, "%s", value);
			}
			return true;
		default:
			REPORT_ERROR(SVLIB_NoPara, "Erroneus parameter mapping in SC_TweakableParameters");
			return false;
	}
  
  return false; //parameter-name unknown
}
