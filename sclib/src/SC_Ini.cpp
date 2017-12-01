/**************************************************************************/
/*    Responsibility:																											*/
/*		  - Provides functionality to read an *.ini-file (like my           */
/*        mod_ini.bas in Visual Basic                                     */
/*																																				*/
/*    Author  : Thilo Stadelmann																					*/
/*    Date    : 14.09.2005																								*/
/**************************************************************************/

#include "SC_Ini.h"
#include <string.h>

//====================================================================================================================
//	The constructor
//====================================================================================================================
SC_Ini::SC_Ini() {
  this->iniFile = NULL;
  this->lastPosition = 0;
  this->lastSection = new char[sclib::bufferSize];
}

//====================================================================================================================
//	The destructor
//====================================================================================================================
SC_Ini::~SC_Ini() {
  if (this->iniFile != NULL) {
    closeIni();
  }
  MFree_1D(this->lastSection);
}

//====================================================================================================================
//	Open an ini-file for reading
//====================================================================================================================
bool SC_Ini::openIni(const char *fileName) {
  if (strncmp(fileName, "", sclib::bufferSize) != 0 && strncmp(fileName, "\0", sclib::bufferSize) != 0) {
	  this->iniFile = fopen(fileName, "r");
    this->lastPosition = 0;
    sprintf(this->lastSection, "\0");
  }
  
  if (this->iniFile == NULL) {
    return false;
  } else {
    return true;
  }
}

//====================================================================================================================
//	Close the formerly opened ini-file
//====================================================================================================================
void SC_Ini::closeIni(void) {
  if (this->iniFile != NULL) {
    fclose(this->iniFile);
    this->iniFile = NULL;
    this->lastPosition = 0;
    sprintf(this->lastSection, "\0");
  }

  return;
}

//====================================================================================================================
//	return the value in the given section with the given key; return NULL if not found
//
//  TODO: - make it case-insensitive by converting all strings to upper case before comparing
//        - allow whitespaces between the key and the '=' in key=value
//====================================================================================================================
char* SC_Ini::readIni(const char* section, const char* key) {
  char *buffer = new char[sclib::bufferSize];
  char *strippedLine;
  char *value = NULL;
  bool rightSection = false;
  char *pos;
  
  if ((this->iniFile == NULL) || (fseek(this->iniFile, 0, SEEK_SET) != 0)) {
    MFree_1D(buffer);
    return NULL;
  }
  
  while (!feof(this->iniFile)) {
    sclib::readline(iniFile, buffer, sclib::bufferSize);
    strippedLine = sclib::trim(buffer);
    
    if (sclib::like(strippedLine, "[\\[]%[\\]]") == true) { //handle section
      rightSection = (strncmp(strippedLine+1, section, strlen(strippedLine)-2) == 0) ? true : false;
      strncpy(this->lastSection, strippedLine+1, strlen(strippedLine)-2);
      this->lastSection[strlen(strippedLine)-2] = '\0';
    } else { //handle key
      if (rightSection == true) { //we are in the right section
        pos = strchr(strippedLine, '=');
        if (pos != NULL) { //strippedLine has the form "key=value"
          if (strncmp(strippedLine, key, pos-strippedLine) == 0) { //we have found the right key, so extract the value
            value = new char[strippedLine + strlen(strippedLine) - pos];
            sprintf(value, "%s\0", pos+1);
            lastPosition = ftell(this->iniFile);
          }
        }
      }
    }

    MFree_1D(strippedLine);
    if (value != NULL) {
      break;
    }
  }

  MFree_1D(buffer);
  return value;
}

//====================================================================================================================
//	just return the next key/value pair out of the ini (last positions is rememberd by readNextParameter() and 
//  readIni()); if section is given, the next pair out of this section is returned; if there is none, key and value 
//  are NULL and false is returned, otherwise they are filled and true is returned
//
//  lines beginning with '#' are regarded as comments and are ignored/overread!
//====================================================================================================================
bool SC_Ini::readNextParameter(char* &key, char* &value, const char* section) {
  char *buffer = new char[sclib::bufferSize];
  char *strippedLine;
  bool rightSection = (section == NULL) ? true : false;
  char *pos;

  MFree_1D(key);
  MFree_1D(value);

  if (section != NULL && strncmp(section, this->lastSection, sclib::bufferSize) != 0) {
    this->lastPosition = 0; //if we want the next entry from a section that differs from the last one, we want the first entry!
  }

  if ((this->iniFile == NULL) || (fseek(this->iniFile, this->lastPosition, SEEK_SET) != 0)) {
    MFree_1D(buffer);
    return false;
  }

  while (!feof(this->iniFile)) {
    sclib::readline(iniFile, buffer, sclib::bufferSize);
    strippedLine = sclib::trim(buffer);
    
    if (strippedLine[0] != '#') { //ignore comment lines
      if (sclib::like(strippedLine, "[\\[]%[\\]]") == true) { //handle section
        rightSection = (section == NULL || strncmp(strippedLine+1, section, strlen(strippedLine)-2) == 0) ? true : false;
        strncpy(this->lastSection, strippedLine+1, strlen(strippedLine)-2);
        this->lastSection[strlen(strippedLine)-2] = '\0';
      } else { //handle key
        if (rightSection == true) { //we are in the right section
          pos = strchr(strippedLine, '=');
          if (pos != NULL) { //strippedLine has the form "key=value"
            key = new char[pos - strippedLine + 1];
            strncpy(key, strippedLine, pos - strippedLine);
            key[pos - strippedLine] = '\0';

            value = new char[strippedLine + strlen(strippedLine) - pos];
            sprintf(value, "%s\0", pos+1);

            lastPosition = ftell(this->iniFile);
          } //form "key=value"
        } // right section
      } //was a key
    } //no comment

    MFree_1D(strippedLine);
    if (value != NULL && key != NULL) {
      break;
    }
  }

  MFree_1D(buffer);
  return (value != NULL && key != NULL);  
}

