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

#ifndef __SC_Timer_H__
#define __SC_Timer_H__

#include <time.h>
#include "SC_Api.h"
#include "SC_TweakableParameters.h"

class SCLIB_API SC_Timer {

	private :
    
    SC_TweakableParameters *pTweak;

    time_t *startTime;
    unsigned long int *elapsedTime;
    unsigned long int timerCount;
    bool alreadyLogged;

	protected:

    void startTimerEx(unsigned long int timerNr, time_t startTime);
    void stopTimerEx(unsigned long int timerNr, time_t stopTime);

	public :

	  SC_Timer(unsigned long int timerCount = 0, SC_TweakableParameters* pTweak = NULL);
		virtual ~SC_Timer();

    void startTimer(unsigned long int timerNr = 0);
    unsigned long int stopTimer(unsigned long int timerNr = 0);
    void resetTimer(unsigned long int timerNr = 0);
    unsigned long int getElapsedTime(unsigned long int timerNr = 0);

    void startTimers(void);
    void stopTimers(void);
    void resetTimers(void);

    void logIt(const char *fileName = ""); //if fileName=="", output is logged on screen (stdout)
};

#endif
