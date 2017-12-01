/**************************************************************************/
/*    Responsibility:																											*/
/*      - Base class for voiced/unvoiced speech detectors                 */
/*																																				*/
/*    Author  : Thilo Stadelmann																					*/
/*    Date    : 03.03.2006																								*/
/**************************************************************************/

#include "SC_Segmentation_VUv.h"

//====================================================================================================================
//	constructor
//====================================================================================================================
SC_Segmentation_VUv::SC_Segmentation_VUv(SC_TweakableParameters* pTweak){
  this->pTweak = pTweak;
}

//====================================================================================================================
//	destructor
//====================================================================================================================
SC_Segmentation_VUv::~SC_Segmentation_VUv(){

}
