/**************************************************************************/
/*	This implements a stopwatch for testing which part of an algorithm  	*/
/*  is most time-consuming. Therefore, it provides an interface for as    */
/*  a specifyable amount of concurrent stopwatches, which can be started  */
/*  and stopped independantly within one class; outputting of the summed  */
/*  time is also encapsulated                                             */
/*																																				*/
/*    Author  : Thilo Stadelmann																					*/
/*    Date    : 03.11.2005																								*/
/**************************************************************************/

#include <stdio.h>
#include "SC_Timer.h"
#include "SC_Aux.h"
#include <SV_Error.h>

//====================================================================================================================
//	constructor; ptweak is only needed for storing the results in the specified debug-dir
//====================================================================================================================
SC_Timer::SC_Timer(unsigned long int timerCount, SC_TweakableParameters* pTweak) {
  this->pTweak = pTweak;
  this->timerCount = 1 + timerCount;
  this->alreadyLogged = false;

  MArray_1D(this->elapsedTime, this->timerCount, unsigned long int, "SC_Timer: elapsedTime");
  MArray_1D(this->startTime, this->timerCount, time_t, "SC_Timer: startTime");

  resetTimers();
  this->elapsedTime[0] = 0;
  this->startTime[0] = time(NULL);
}

//====================================================================================================================
//	destructor
//====================================================================================================================
SC_Timer::~SC_Timer() {
  stopTimers();
  logIt();

  MFree_1D(this->startTime);
  MFree_1D(this->elapsedTime);
}

//====================================================================================================================
//	start the specified stopwatch
//====================================================================================================================
void SC_Timer::startTimer(unsigned long int timerNr) {
  if (timerNr > 0 && timerNr < this->timerCount) {
    this->startTime[timerNr] = time(NULL);
  }

  return;
}

//====================================================================================================================
//	needed to start all stopwatches at the same time
//====================================================================================================================
void SC_Timer::startTimerEx(unsigned long int timerNr, time_t startTime) {
  if (timerNr > 0 && timerNr < this->timerCount) {
    this->startTime[timerNr] = startTime;
  }

  return;
}

//====================================================================================================================
//	start all stopwatches at the same time
//====================================================================================================================
void SC_Timer::startTimers(void) {
  time_t startTime = time(NULL);

  for (unsigned long int i = 1; i < this->timerCount; i++) {
    startTimerEx(i, startTime);
  }

  return;
}

//====================================================================================================================
//	stop the specified stopwatch, adding the elapsed time to the value already stored in elapseTime
//  the elapsed time since the last start (not the completete accumulated elapsd time) is returned
//  if this stopwatch wasn't started, 0 is returned
//  the startTime for this stopwatch is set to 0 to indicated that it is not running
//====================================================================================================================
unsigned long int SC_Timer::stopTimer(unsigned long int timerNr) {
  time_t stopTime = time(NULL);
  unsigned long int res = 0;

  if (timerNr > 0 && timerNr < this->timerCount) {
    if (this->startTime[timerNr] != 0) {
      res = (unsigned long int)(stopTime - this->startTime[timerNr]);
      this->elapsedTime[timerNr] += res;
      this->startTime[timerNr] = 0;
    }
  }

  return res;
}

//====================================================================================================================
//	needed to stop all stopwatches at the same time
//  the startTime for this stopwatch is set to 0 to indicated that it is not running
//====================================================================================================================
void SC_Timer::stopTimerEx(unsigned long int timerNr, time_t stopTime) {
  if (timerNr > 0 && timerNr < this->timerCount) {
    if (this->startTime[timerNr] != 0) {
      this->elapsedTime[timerNr] += (unsigned long int)(stopTime - this->startTime[timerNr]);
      this->startTime[timerNr] = 0;
    }
  }

  return;
}

//====================================================================================================================
//	stop all stopwatches at the same time
//====================================================================================================================
void SC_Timer::stopTimers(void) {
  time_t stopTime = time(NULL);

  for (unsigned long int i = 1; i < this->timerCount; i++) {
    stopTimerEx(i, stopTime);
  }

  return;
}

//====================================================================================================================
//	reset the specified stopwatch by zeroing the elapsedTime and the startTime
//====================================================================================================================
void SC_Timer::resetTimer(unsigned long int timerNr) {
  if (timerNr > 0 && timerNr < this->timerCount) {
    this->startTime[timerNr] = 0;
    this->elapsedTime[timerNr] = 0;
  }

  return;
}

//====================================================================================================================
//	reset all stopwatches
//====================================================================================================================
void SC_Timer::resetTimers(void) {
  for (unsigned long int i = 1; i < this->timerCount; i++) {
    resetTimer(i);
  }

  return;
}

//====================================================================================================================
//	get the accumulated some of all the elapsed time for the specified stopwatch in seconds [s]
//  if the stopwatch is currently running, return the time elapsed until the moment of this evaluation
//====================================================================================================================
unsigned long int SC_Timer::getElapsedTime(unsigned long int timerNr) {
  if (timerNr >= 0 && timerNr < this->timerCount) {
    return this->elapsedTime[timerNr] + ((this->startTime[timerNr] != 0) ? (unsigned long int)(time(NULL) - this->startTime[timerNr]) : 0);
  } else {
    return 0;
  }
}

//====================================================================================================================
//	log the elapsedTime of all stopWatches
//  if fileName=="", output is logged on screen (stdout)
//  TODO: not very sophisticated, needs even more improvement to log on screen/file and for single stopwatches...
//====================================================================================================================
void SC_Timer::logIt(const char *fileName) {
	char *buffer = new char[sclib::bufferSize*this->timerCount];
	
	if (alreadyLogged == false) { //make a newline before the first output
    sprintf(buffer, "\n");
    this->alreadyLogged = true;
	} else {
		sprintf(buffer, "");
	}

  sprintf(buffer, "%s  overall time: %ds\n", buffer, getElapsedTime(0));
  for (unsigned long int i = 1; i < this->timerCount; i++) {
    sprintf(buffer, "%s  elapsed time for timer %d: %ds (ratio to overall time: %.2f%%)\n", buffer, i, getElapsedTime(i), (double)getElapsedTime(i) / (double)getElapsedTime(0) * 100);
  }

	if (fileName==NULL  || strncmp(fileName, "", sclib::bufferSize)==0) {
		printf("%s", buffer);
	} else {
		sclib::stringOut(fileName, buffer, this->pTweak);
	}

	MFree_1D(buffer);

  return;
}
