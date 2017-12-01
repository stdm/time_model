/**************************************************************************/
/*    Responsibility:																											*/
/*		  - Provides functionality to read an *.ini-file (like my           */
/*        mod_ini.bas in Visual Basic                                     */
/*																																				*/
/*    Author  : Thilo Stadelmann																					*/
/*    Date    : 14.09.2005																								*/
/**************************************************************************/

#ifndef __SC_Ini_H__
#define __SC_Ini_H__

#include <iostream>
#include "SC_Api.h"
#include "SC_Aux.h"

class SCLIB_API SC_Ini {
	private :

  protected :

    FILE *iniFile;
    long int lastPosition;
    char* lastSection;

  public :
		
    SC_Ini();
    virtual ~SC_Ini();
    
    //====================================================================================================================
    //	Open an ini-file for reading
    //====================================================================================================================
    bool openIni(const char *fileName);

    //====================================================================================================================
    //	Close the formerly opened ini-file
    //====================================================================================================================
    void closeIni(void);

    //====================================================================================================================
    //	return the value in the given section with the given key; return NULL if not found
    //
    //  TODO: - make it case-insensitive by converting all strings to upper case before comparing
    //        - allow whitespaces between the key and the '=' in key=value
    //====================================================================================================================
    char* readIni(const char* section, const char* key);

    //====================================================================================================================
    //	just return the next key/value pair out of the ini (last positions is rememberd by readNextParameter() and 
    //  readIni()); if section is given, the next pair out of this section is returned; if there is none, key and value 
    // are NULL and false is returned, otherwise they are filled and true is returned
    //====================================================================================================================
    bool readNextParameter(char* &key, char* &value, const char* section = NULL);
};

#endif

