//************************************************************************
//    Implementation of SV_DataIO class
//
//
//	  TODO : to read none native format (Big endian - little endian)
//
//    Author  : Jialong HE
//    Date    : April 27, 1999
//************************************************************************
#include <stdlib.h>
#include <stdio.h>
#include <fstream>
#include <string.h>
#include "SV_DataIO.h"
#include "SV_Data.h"
#include "SV_Error.h"

static char SV_LibID[] = "Copyright (c) by Jialong He";
//==========================================
// default constructor
// Fill File header with default values
//==========================================
SV_DataIO::SV_DataIO() {

}

//==========================================
// destructor
//==========================================
SV_DataIO::~SV_DataIO() {
//	DFile.close();	
}

//==========================================
// Close DFile, used for refopen 
//==========================================
void SV_DataIO::CloseFile(void) {
	this->DFile.clear(); //by thilo to kill errors on opening this file again if there where errors (e.g. reading after eof) in the last try
	DFile.close();
}

//===============================================
// Open stream DFile for READ_REC or WRITE_REC 
//===============================================
void SV_DataIO::OpenFile(const char *FName, int Mode) { //by thilo: added const

	if (Mode == WRITE_REC) {
		//DFile.open(FName, ios::app|ios::binary);
    DFile.open(FName, ios_base::out|ios_base::app|ios_base::binary); //by thilo: because of changes from vc6.0 to vc.net 2003
	}
	
	if (Mode == READ_REC) {
		//DFile.open(FName, ios::in|ios::binary);
    DFile.open(FName, ios_base::in|ios_base::binary); //by thilo: because of changes from vc6.0 to vc.net 2003
	}
	
	if (DFile.fail()) {
		REPORT_ERROR(SVLIB_FileErr, "OpenFile");
	}

}

//==========================================
// Put one data record into file
//==========================================
int	SV_DataIO::PutDataRec(SV_Data &Data) {

	//long DataBytes; //by thilo
	int bytes = 0; //by thilo

//	if (~Data.Valid()) {
//		REPORT_ERROR(SVLIB_BadArg, "PutDataRec");
//	} 

	//------------------------------------
	// Write RecordHdr into file
	//------------------------------------
	//DFile.write((char*)(&Data.Hdr), DHLen);
	//by thilo:
	SV_DataIO::SV_DatatypeSizes sizes;
	getCurrentDatatypeSizes(sizes);
	bytes += writeMachineHeader(&(this->DFile), sizes);
	bytes += writeArray(&(this->DFile), Data.Hdr.Signature, 8);
	bytes += writeScalar(&(this->DFile), Data.Hdr.ByteOrder);
	bytes += writeScalar(&(this->DFile), Data.Hdr.Version);	
	bytes += writeScalar(&(this->DFile), Data.Hdr.ID);
	bytes += writeArray(&(this->DFile), Data.Hdr.Name, 16);
	bytes += writeScalar(&(this->DFile), Data.Hdr.frameSize);
	bytes += writeScalar(&(this->DFile), Data.Hdr.frameStep);
	bytes += writeScalar(&(this->DFile), Data.Hdr.sampleRate);
	bytes += writeArray(&(this->DFile), Data.Hdr.NotUsed, NotUsedDHLen);
	//end by thilo

	//DFile.write((char*)(&Data.Row), sizeof(long));
	bytes += writeScalar(&(this->DFile), Data.Row); //by thilo

	//DFile.write((char*)(&Data.Col), sizeof(long));
	bytes += writeScalar(&(this->DFile), Data.Col); //by thilo
	
	if (DFile.good() != TRUE) {
		REPORT_ERROR(SVLIB_FileErr, "PutDataRec"); 
	}
	
	//------------------------------------
	// Write Data body into file
	//------------------------------------
	//DataBytes = Data.Row * Data.Col * sizeof(float);  //by thilo
	//DFile.write((char*)(Data.Mat[0]), DataBytes);
	bytes += writeMatrix(&(this->DFile), Data.Mat, Data.Row, Data.Col); //by thilo

	if (DFile.good() != TRUE) {
		REPORT_ERROR(SVLIB_FileErr, "PutDataRec"); 
	}

	//return(DHLen+DataBytes); //by thilo
	return bytes;
}

//==========================================
// Get one data record from stream
// if no record, return NULL
//==========================================
SV_Data *SV_DataIO::GetDataRec(void) {

	//long DataBytes; //by thilo
	SV_Data *Data;

	Data = new SV_Data;
	if (Data==NULL) {
		REPORT_ERROR(SVLIB_NoMem, "GetDataRec"); 
	}
	//------------------------------------
	// Read RecordHdr from file
	//------------------------------------
	//DFile.read((char*)&(Data->Hdr), DHLen);
	//by thilo
	SV_DataIO::SV_DatatypeSizes fileSizes, codeSizes;
	getCurrentDatatypeSizes(codeSizes);
	int bytes = readMachineHeader(&(this->DFile), fileSizes, true);
	if (bytes > 0) {
		consumeBytes(&(this->DFile), bytes);
	} else {
		bytes = 0;
	}
	bytes += readArray(&(this->DFile), Data->Hdr.Signature, 8, codeSizes, fileSizes);
	bytes += readScalar(&(this->DFile), Data->Hdr.ByteOrder, codeSizes, fileSizes);
	bytes += readScalar(&(this->DFile), Data->Hdr.Version, codeSizes, fileSizes);	
	bytes += readScalar(&(this->DFile), Data->Hdr.ID, codeSizes, fileSizes);
	bytes += readArray(&(this->DFile), Data->Hdr.Name, 16, codeSizes, fileSizes);
	bytes += readScalar(&(this->DFile), Data->Hdr.frameSize, codeSizes, fileSizes);
	bytes += readScalar(&(this->DFile), Data->Hdr.frameStep, codeSizes, fileSizes);
	bytes += readScalar(&(this->DFile), Data->Hdr.sampleRate, codeSizes, fileSizes);
	bytes += consumeArray(&(this->DFile), Data->Hdr.NotUsed, getNotUsedHeaderSize(fileSizes), codeSizes, fileSizes);
	//end by thilo

	//DFile.read((char*)&(Data->Row), sizeof(long));
	bytes += readScalar(&(this->DFile), Data->Row, codeSizes, fileSizes); //by thilo
	//DFile.read((char*)&(Data->Col), sizeof(long));
	bytes += readScalar(&(this->DFile), Data->Col, codeSizes, fileSizes); //by thilo
	if (DFile.good() != TRUE) {
    DFile.clear(); //by thilo: clear error from reading over EOF...
		delete Data;
		return(NULL);
	}

	Data->Alloc();  // allocate space for this record
	//------------------------------------
	// Read Data body into memory
	//------------------------------------
	//DataBytes = Data->Row * Data->Col * sizeof(float);  //by thilo
	//DFile.read((char*)(Data->Mat[0]), DataBytes);
	bytes += readMatrix(&(this->DFile), Data->Mat, Data->Row, Data->Col, codeSizes, fileSizes);

	if (DFile.good() != TRUE) {
		delete Data;
		Data = NULL;
	}

	return(Data);
};

//==========================================
// by thilo: Get one data record from stream
// if no record, return NULL (old method)
//==========================================
SV_Data *SV_DataIO::oldGetDataRec(void) {

	long DataBytes;
	SV_Data *Data;

	Data = new SV_Data;
	if (Data==NULL) {
		REPORT_ERROR(SVLIB_NoMem, "GetDataRec"); 
	}
	//------------------------------------
	// Read RecordHdr from file
	//------------------------------------
	DFile.read((char*)&(Data->Hdr), DHLen);
	DFile.read((char*)&(Data->Row), sizeof(long));
	DFile.read((char*)&(Data->Col), sizeof(long));
	if (DFile.good() != TRUE) {
    DFile.clear(); //by thilo: clear error from reading over EOF...
		delete Data;
		return(NULL);
	}

	Data->Alloc();  // allocate space for this record
	//------------------------------------
	// Read Data body into memory
	//------------------------------------
	DataBytes = Data->Row * Data->Col * sizeof(float); 
	DFile.read((char*)(Data->Mat[0]), DataBytes);
	
	if (DFile.good() != TRUE) {
		delete Data;
		Data = NULL;
	}

	return(Data);
};

//by thilo: returns the overall number of rows (and the number of cols, if equal in all recs, otherwise 0) of all the datarecs in the file
long int SV_DataIO::getAllRecRowCount(int &cols) {
	int bytes;
	long rowCount = 0, colCount = 0, count = 0;
	SV_Data *pData = new SV_Data();
	SV_DataIO::SV_DatatypeSizes fileSizes, codeSizes;

	//------------------------------------
	// Read RecordHdr into file
	//------------------------------------
  while (1) {
		// Read header
		getCurrentDatatypeSizes(codeSizes);
		bytes = readMachineHeader(&(this->DFile), fileSizes, false);
		bytes += readArray(&(this->DFile), pData->Hdr.Signature, 8, codeSizes, fileSizes);
		bytes += readScalar(&(this->DFile), pData->Hdr.ByteOrder, codeSizes, fileSizes);
		bytes += readScalar(&(this->DFile), pData->Hdr.Version, codeSizes, fileSizes);	
		bytes += readScalar(&(this->DFile), pData->Hdr.ID, codeSizes, fileSizes);
		bytes += readArray(&(this->DFile), pData->Hdr.Name, 16, codeSizes, fileSizes);
		bytes += readScalar(&(this->DFile), pData->Hdr.frameSize, codeSizes, fileSizes);
		bytes += readScalar(&(this->DFile), pData->Hdr.frameStep, codeSizes, fileSizes);
		bytes += readScalar(&(this->DFile), pData->Hdr.sampleRate, codeSizes, fileSizes);
		bytes += consumeArray(&(this->DFile), pData->Hdr.NotUsed, getNotUsedHeaderSize(fileSizes), codeSizes, fileSizes);
		bytes += readScalar(&(this->DFile), pData->Row, codeSizes, fileSizes);
		bytes += readScalar(&(this->DFile), pData->Col, codeSizes, fileSizes);

		if (DFile.good() != TRUE) {
			break;  // jump out while (1)
		}

    rowCount += pData->Row;
    colCount = (count == 0) ? pData->Col : ((colCount != pData->Col) ? 0 : colCount);

		//really read the data-matrix using the machine-independant functions instead of seeking... imposes less assumptions on file-architecture even if it is slower...
		//DataBytes = pData->Row * pData->Col * sizeof(float); 
    //DFile.seekg(DataBytes, ios::cur);  // ignore body
		bytes += consumeMatrix(&(this->DFile), pData->Mat, pData->Row, pData->Col, codeSizes, fileSizes);

    count++;
	}

  DFile.clear(); //clear error from reading over EOF...
  DFile.seekg(0, ios_base::beg); //jump back to start
  delete pData;
  cols = colCount;

  return rowCount;
}

//==========================================
// Get one data record with RecID
// if no record, return NULL
//==========================================
SV_Data *SV_DataIO::GetDataRec(int RecID) {

	//long DataBytes; //by thilo
	SV_Data *Data;
	int  Found = 0;

	Data = new SV_Data;
	if (Data==NULL) {
		REPORT_ERROR(SVLIB_NoMem, "GetDataRec(RecID)"); 
	}

	//------------------------------------
	// Read RecordHdr into file
	//------------------------------------
	while (1) {

		//--------------------------------
		// Read header,
		//--------------------------------
		//DFile.read((char*)&(Data->Hdr), DHLen);
		//by thilo:
		SV_DataIO::SV_DatatypeSizes fileSizes, codeSizes;
		getCurrentDatatypeSizes(codeSizes);
		int bytes = readMachineHeader(&(this->DFile), fileSizes, true);
		if (bytes > 0) {
			consumeBytes(&(this->DFile), bytes);
		} else {
			bytes = 0;
		}
		bytes += readArray(&(this->DFile), Data->Hdr.Signature, 8, codeSizes, fileSizes);
		bytes += readScalar(&(this->DFile), Data->Hdr.ByteOrder, codeSizes, fileSizes);
		bytes += readScalar(&(this->DFile), Data->Hdr.Version, codeSizes, fileSizes);	
		bytes += readScalar(&(this->DFile), Data->Hdr.ID, codeSizes, fileSizes);
		bytes += readArray(&(this->DFile), Data->Hdr.Name, 16, codeSizes, fileSizes);
		bytes += readScalar(&(this->DFile), Data->Hdr.frameSize, codeSizes, fileSizes);
		bytes += readScalar(&(this->DFile), Data->Hdr.frameStep, codeSizes, fileSizes);
		bytes += readScalar(&(this->DFile), Data->Hdr.sampleRate, codeSizes, fileSizes);
		bytes += consumeArray(&(this->DFile), Data->Hdr.NotUsed, getNotUsedHeaderSize(fileSizes), codeSizes, fileSizes);
		//end by thilo		

		//DFile.read((char*)&(Data->Row), sizeof(long));
		bytes += readScalar(&(this->DFile), Data->Row, codeSizes, fileSizes); //by thilo
		//DFile.read((char*)&(Data->Col), sizeof(long));
		bytes += readScalar(&(this->DFile), Data->Col, codeSizes, fileSizes); //by thilo

		if (DFile.good() != TRUE) {
			break;  // jump out while (1)
		}

		//DataBytes = Data->Row * Data->Col * sizeof(float);  //by thilo
		//----------------------------------
		// Found desired record, read body
		//----------------------------------
		if (Data->Hdr.ID == RecID) {
			Data->Alloc();  // allocate space for this record
			//DFile.read((char*)(Data->Mat[0]), DataBytes);
			bytes += readMatrix(&(this->DFile), Data->Mat, Data->Row, Data->Col, codeSizes, fileSizes); //by thilo
		  if (DFile.good() != TRUE) {
				break;
			} 
			else {
				Found = 1;
				break;
			}
		}
		else { // ignore body
			//DFile.seekg(DataBytes, ios::cur);
			bytes += consumeMatrix(&(this->DFile), Data->Mat, Data->Row, Data->Col, codeSizes, fileSizes); //by thilo
		}  
	}

	//-------------------------------
	//  
	//-------------------------------
	if (Found) {
		return(Data);
	}
	else {
		delete Data;
		return(NULL);
	}

};

//==========================================
// by thilo: Get one data record with RecID
// if no record, return NULL (old method)
//==========================================
SV_Data *SV_DataIO::oldGetDataRec(int RecID) {

	long DataBytes;
	SV_Data *Data;
	int  Found = 0;

	Data = new SV_Data;
	if (Data==NULL) {
		REPORT_ERROR(SVLIB_NoMem, "GetDataRec(RecID)"); 
	}

	//------------------------------------
	// Read RecordHdr into file
	//------------------------------------
	while (1) {

		//--------------------------------
		// Read header,
		//--------------------------------
		DFile.read((char*)&(Data->Hdr), DHLen);
		DFile.read((char*)&(Data->Row), sizeof(long));
		DFile.read((char*)&(Data->Col), sizeof(long));
		if (DFile.good() != TRUE) {
			break;  // jump out while (1)
		}

		DataBytes = Data->Row * Data->Col * sizeof(float);  //by thilo
		//----------------------------------
		// Found desired record, read body
		//----------------------------------
		if (Data->Hdr.ID == RecID) {
			Data->Alloc();  // allocate space for this record
		  if (DFile.good() != TRUE) {
				break;
			} 
			else {
				Found = 1;
				break;
			}
		}
		else { // ignore body
			DFile.seekg(DataBytes, ios::cur);
		}  
	}

	//-------------------------------
	//  
	//-------------------------------
	if (Found) {
		return(Data);
	}
	else {
		delete Data;
		return(NULL);
	}

};

//==========================================
// Get ALL data records from stream
// return linked list head
//==========================================
SV_Data *SV_DataIO::GetAllRec(void) {

	SV_Data *Head, *Curr;

	Head = GetDataRec();
	Curr = Head;
	while (Curr != NULL) {
		Curr->Next = GetDataRec();
		Curr = Curr->Next;
	}

	return(Head);
};

//==========================================
// by thilo: Get ALL data records from 
// stream return linked list head
//==========================================
SV_Data *SV_DataIO::oldGetAllRec(void) {

	SV_Data *Head, *Curr;

	Head = oldGetDataRec();
	Curr = Head;
	while (Curr != NULL) {
		Curr->Next = oldGetDataRec();
		Curr = Curr->Next;
	}

	return(Head);
};

//=============================================
// Get ALL data records from stream with RecID
// return linked list head
//=============================================
SV_Data *SV_DataIO::GetAllRec(int RecID) {

	SV_Data *Head, *Curr;

	Head = GetDataRec(RecID);
	Curr = Head;
	while (Curr != NULL) {
		Curr->Next = GetDataRec(RecID);
		Curr = Curr->Next;
	}

	return(Head);
};

//=============================================
// by thilo: Get ALL data records from stream 
// with RecID return linked list head
//=============================================
SV_Data *SV_DataIO::oldGetAllRec(int RecID) {

	SV_Data *Head, *Curr;

	Head = oldGetDataRec(RecID);
	Curr = Head;
	while (Curr != NULL) {
		Curr->Next = oldGetDataRec(RecID);
		Curr = Curr->Next;
	}

	return(Head);
};

//====================================================================================================================
//	By thilo: What follows is code that helps making *ALL* binary file i/o (not only of SV_Data objects, but also for
//            models, signals, ...) read-/writeable on systems with different endianness and/or 32bit/64bit
//            It has to be put here (in the SV_Lib) 'cause some binary i/o already happens here...
//====================================================================================================================

//====================================================================================================================
//  Each binary file should have a header containing information on the datatype-sizes (and endianness) of the machine 
//  that wrote it; this is extracted and returned (it is assumed that the file is at the correct position for the 
//  header, typically the very beginning). If the header is missing (to be compatibel with old binaries), a standard 
//  header is returned assuming that the file was written on a 32bit x86 machine. If retorePosition==true, the 
//  read-cursor of the file remains unaltered after the function call (also if no valid header was found)
//
//  The machine-header has the following form (22 bytes long):
//  |S|V|_|L|I|B|_|H|D|R|:|bool_size|char_size|short_size|int_size|long_size|long_long_size|float_size|double_size|long_double_size|endianness|\0|
//  all sizes are 1-byte binary values (not ASCII) telling the data-type's size in bytes
//====================================================================================================================
int SV_DataIO::readMachineHeader(istream *file, SV_DataIO::SV_DatatypeSizes &sizes, bool restorePosition) {
	char machineHeader[22];
	int initPos = file->tellg(); //save the position of the read-cursor
	bool foundHeader = false;

	file->read(machineHeader, 22);

	if (strncmp(machineHeader, "SV_LIB_HDR:", 11) == 0 && machineHeader[21] == 0) { //we found a valid header
		sizes.boolSize = machineHeader[11];
		sizes.charSize = machineHeader[12];
		sizes.shortSize = machineHeader[13];
		sizes.intSize = machineHeader[14];
		sizes.longSize = machineHeader[15];
		sizes.longLongSize = machineHeader[16];
		sizes.floatSize = machineHeader[17];
		sizes.doubleSize = machineHeader[18];
		sizes.longDoubleSize = machineHeader[19];
		sizes.endianness = machineHeader[20];
		foundHeader = true;
	} else { //no valid header found => assume that the file was written on as 32bit little endian (x86) machine
		sizes.boolSize = 1;
		sizes.charSize = 1;
		sizes.shortSize = 2;
		sizes.intSize = 4;
		sizes.longSize = 4;
		sizes.longLongSize = 8;
		sizes.floatSize = 4;
		sizes.doubleSize = 8;
		sizes.longDoubleSize = 8;
		sizes.endianness = SVLIB_LITTLE_ENDIAN;
	}

	if (restorePosition == true || foundHeader == false) {
		file->seekg(initPos);
	}

	return (foundHeader == true) ? 22 : SVLIB_Fail;	
}

//====================================================================================================================
//  The given machine-header is written into the given (write-opened) file at the current position
//====================================================================================================================
int SV_DataIO::writeMachineHeader(ostream *file, SV_DataIO::SV_DatatypeSizes sizes) {
	char machineHeader[22];

	machineHeader[0] = 'S';
	machineHeader[1] = 'V';
	machineHeader[2] = '_';
	machineHeader[3] = 'L';
	machineHeader[4] = 'I';
	machineHeader[5] = 'B';
	machineHeader[6] = '_';
	machineHeader[7] = 'H';
	machineHeader[8] = 'D';
	machineHeader[9] = 'R';
	machineHeader[10] = ':';
	machineHeader[11] = sizes.boolSize;
	machineHeader[12] = sizes.charSize;
	machineHeader[13] = sizes.shortSize;
	machineHeader[14] = sizes.intSize;
	machineHeader[15] = sizes.longSize;
	machineHeader[16] = sizes.longLongSize;
	machineHeader[17] = sizes.floatSize;
	machineHeader[18] = sizes.doubleSize;
	machineHeader[19] = sizes.longDoubleSize;
	machineHeader[20] = sizes.endianness;
	machineHeader[21] = '\0';

	file->write(machineHeader, 22);

	return 22;
}

//====================================================================================================================
//	Returns information about the endianness and datatype-sizes of the current machine/implementation
//====================================================================================================================
SV_DataIO::SV_DatatypeSizes& SV_DataIO::getCurrentDatatypeSizes(SV_DataIO::SV_DatatypeSizes& sizes) {
	sizes.boolSize = sizeof(bool);
	sizes.charSize = sizeof(char);
	sizes.doubleSize = sizeof(double);
	sizes.floatSize = sizeof(float);
	sizes.intSize = sizeof(int);
	sizes.longSize = sizeof(long int);
	sizes.shortSize = sizeof(short int);
	sizes.longLongSize = sizeof(long long int);
	sizes.longDoubleSize = sizeof(long double);
	sizes.endianness = getEndianness();

	return sizes;		
}

//====================================================================================================================
//  Returns endianness of current machine
//====================================================================================================================
unsigned char SV_DataIO::getEndianness(void) {
	unsigned char endianTest[2] = {1, 0}; //is the bit-pattern '0000000100000000' interpreted as 1 (=> little endian) or 256 (=> big endian) when casted to short?

	return (*(short*)endianTest == 1) ? SVLIB_LITTLE_ENDIAN : SVLIB_BIG_ENDIAN;
}

//====================================================================================================================
//  Inverts the byteorder of the given value (e.g. to correct for different endianness)
//====================================================================================================================
void SV_DataIO::invertByteOrder(char *value, unsigned char byteCount) {
	char tmp;

	for (int i = 0; i < byteCount/2; i++) {
		tmp = value[i];
		value[i] = value[byteCount-i];
		value[byteCount-i] = tmp;
	}

	return;
}

//====================================================================================================================
//  Regard swapLength bytes in the buffer as a unit and invert their order, for all those units in the buffer; 
//  return false if numberOfBytes/swapLength is not an integer
//====================================================================================================================
bool SV_DataIO::swapBytes(char* buffer, int swapLength, unsigned long int numberOfBytes) {
  unsigned long int numberOfUnits;
	bool res = false;

	if (numberOfBytes % swapLength == 0) {
		numberOfUnits = numberOfBytes / swapLength;

		for (unsigned long int x = 0; x < numberOfUnits; x += swapLength) {
			invertByteOrder(&buffer[x], swapLength);
		}
	}

	return res;
}

//====================================================================================================================
//  Read and discard 'bytes' bytes
//====================================================================================================================
int SV_DataIO::consumeBytes(istream *file, unsigned int bytes) {
	char *value;

	MArray_1D(value, bytes, char, "SV_DataIO.consumeBytes: value");
	file->read(value, bytes);
	MFree_1D(value);

	return (file->good() == TRUE) ? (int)(bytes) : SVLIB_Fail;
}

//====================================================================================================================
//  Returns the size of the NotUsed-member of the header-struct based on the knowledge of the number and type of 
//  members and their sizes as given in "sizes"
//====================================================================================================================
int SV_DataIO::getNotUsedHeaderSize(SV_DataIO::SV_DatatypeSizes sizes) {
	int bytes = DHLen - 24*sizes.charSize - 2*sizes.longSize - sizes.floatSize - 3*sizes.intSize; //attention: this uses explicit knowledge of SV_Data-internals!

	return bytes;
}

//====================================================================================================================
//	End by thilo
//====================================================================================================================
