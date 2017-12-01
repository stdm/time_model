/**************************************************************************/
/*	This class implements the cryptographic hash function	MD5 (message		*/
/*  digest 5); the implementation is due to Benjamin Neumann,             */
/*  http://www.ben-newman.de/com/MD5.php.                                 */
/*                                                                        */
/*  Usage: call Update() (maybe several times if data comes in not at     */
/*  once) to compute the checksum; call Final() to return the checksum    */
/*  when all data has been fed.                                           */
/*																																				*/
/*    Author  : Benjamin Neumann, Thilo Stadelman													*/
/*    Date    : 08.04.2009																								*/
/**************************************************************************/

#include "SC_MD5.h"
#include "SC_Aux.h" //by thilo
#include <string.h> //by thilo

/*****************************************************************************************
FUNKTION:              SC_MD5::RotateLeft
DETAILS:               privat
BESCHREIBUNG:          verschiebt die Bits eines 32 bit DWORDs nach Links um einen bestimmten WErt
RÜCKGABEWERT:          Der verschobene DWORD -Wert 
ARGUMENTE:             DWORD x : der zu rotierende Wert
                       int n   : die Anzahl der Bits, um die verschoben wird
*****************************************************************************************/
DWORD SC_MD5::RotateLeft(DWORD x, int n)
{
        //verschieben und x zurückgeben
        return (x << n) | (x >> (32-n));
}
 
 
/*****************************************************************************************
FUNKTION :             SC_MD5::FF
DETAILS:               protected
BESCHREIBUNG:          Implementierung des allgemeinen MD5 Transformationsalgorithmus
RÜCKGABEWERT:          keiner
ARGUMENTE:             DWORD &A, B, C, D : Aktuelle (Teil-) Checksum
                       DWORD X     : Eingabedaten
                       DWORD S     : MD5_SXX Transformationskonstante
                       DWORD T     : MD5_TXX Transformationskonstante
ANMERKUNGEN:           Keine
*****************************************************************************************/
void SC_MD5::FF( DWORD& A, DWORD B, DWORD C, DWORD D, DWORD X, DWORD S, DWORD T)
{
        DWORD F = (B & C) | (~B & D);
        A += F + X + T;
        A = RotateLeft(A, S);
        A += B;
}
 
 
/*****************************************************************************************
FUNKTION:              SC_MD5::GG
DETAILS:               protected
BESCHREIBUNG:          Implementierung des allgemeinen MD5 Transformationsalgorithmus
RÜCKGABEWERT:          keiner
ARGUMENTE:             DWORD &A, B, C, D : Aktuelle (Teil-) Checksum
                       DWORD X     : Eingabedaten
                       DWORD S     : MD5_SXX Transformationskonstante
                       DWORD T     : MD5_TXX Transformationskonstante
ANMERKUNGEN:           Keine
*****************************************************************************************/
void SC_MD5::GG( DWORD& A, DWORD B, DWORD C, DWORD D, DWORD X, DWORD S, DWORD T)
{
        DWORD G = (B & D) | (C & ~D);
        A += G + X + T;
        A = RotateLeft(A, S);
        A += B;
}
 
 
/*****************************************************************************************
FUNKTION:              SC_MD5::HH
DETAILS:               protected
BESCHREIBUNG:          Implementierung des allgemeinen MD5 Transformationsalgorithmus
RÜCKGABEWERT:          keiner
ARGUMENTE:             DWORD &A, B, C, D : Aktuelle (Teil-) Checksum
                       DWORD X     : Eingabedaten
                       DWORD S     : MD5_SXX Transformationskonstante
                       DWORD T     : MD5_TXX Transformationskonstante
ANMERKUNGEN:           Keine
*****************************************************************************************/
void SC_MD5::HH( DWORD& A, DWORD B, DWORD C, DWORD D, DWORD X, DWORD S, DWORD T)
{
        DWORD H = (B ^ C ^ D);
        A += H + X + T;
        A = RotateLeft(A, S);
        A += B;
}
 
 
/*****************************************************************************************
FUNKTION:              SC_MD5::II
DETAILS:               protected
BESCHREIBUNG:          Implementierung des allgemeinen MD5 Transformationsalgorithmus
RÜCKGABEWERT:          keiner
ARGUMENTE:             DWORD &A, B, C, D : Aktuelle (Teil-) Checksum
                       DWORD X     : Eingabedaten
                       DWORD S     : MD5_SXX Transformationskonstante
                       DWORD T     : MD5_TXX Transformationskonstante
ANMERKUNGEN:           Keine
*****************************************************************************************/
void SC_MD5::II( DWORD& A, DWORD B, DWORD C, DWORD D, DWORD X, DWORD S, DWORD T)
{
        DWORD I = (C ^ (B | ~D));
        A += I + X + T;
        A = RotateLeft(A, S);
        A += B;
}
 
 
/*****************************************************************************************
FUNKTION:              SC_MD5::ByteToDWord
DETAILS:               private
BESCHREIBUNG:          Transferiert die Daten eines 8 bit array in ein 32 bit array
RÜCKGABEWERTE:         void
ARGUMENTE:             DWORD* Output : 32 bit (unsigned long) Zielarray 
                       BYTE* Input   : 8 bit (unsigned char) Quellarray
                       UINT nLength  : Anzahl an 8 bit Dateneinträgen im Quellarray
BEMERKUNG:             Vier BYTES aus dem Quellarray werden in jeden DWORD Eintrag des Zielarrays
                       transferiert. Das erste BYTE wird in die bits (0-7) 
                       des Ziel DWORD kopiert, das zweite BYTE in die bits 8-15 usw. 
                       Der Algorithmus nimmt an, dass das Quellarray ein Vielfaches von 4 bytes lang ist,
                       so dass es perfekt in ein Array von 32 bit words passt.
*****************************************************************************************/
void SC_MD5::ByteToDWord(DWORD* Output, BYTE* Input, UINT nLength, int& error)
{
        //entry invariants
        if( nLength % 4 != 0 ) error = -23;
 
        //Initialisierung
        UINT i=0;       //index des Zeilarrays
        UINT j=0;       //index des Quellarrays
 
        //transferiere die Daten durch shifting und kopieren
        for ( ; j < nLength; i++, j += 4)
        {
                Output[i] = (DWORD)Input[j]               | 
                                        (DWORD)Input[j+1] << 8  | 
                                        (DWORD)Input[j+2] << 16 | 
                                        (DWORD)Input[j+3] << 24;
        }
}
 
/*****************************************************************************************
FUNKTION:              SC_MD5::Transform
DETAILS:               protected
BESCHREIBUNG:          MD5 Transformationsalgorithmus;  transformiert 'm_lMD5'
RÜCKGABEWERTE:         keiner
ARGUMENTE:             BYTE Block[64]
BEMERKUNG:             Eine MD5 checksum wird in vier 'Transformations-'Durchgängen berechnet.
                       Die MD5 checksum, aktuell in m_lMD5 gesichert, ist verflochten mit dem 
                       Transformationsprozess der Daten in 'Block'.  
*****************************************************************************************/
void SC_MD5::Transform(BYTE Block[64], int& error)
{
        //initialisiere die lokalen Daten mit der aktuellen checksum
        ULONG a = m_lMD5[0];
        ULONG b = m_lMD5[1];
        ULONG c = m_lMD5[2];
        ULONG d = m_lMD5[3];
 
        //Kopiert BYTES aus input 'Block' in ein Array von ULONGS 'X'
        ULONG X[16];
        ByteToDWord( X, Block, 64, error);
 
        //Runde 1 Transformation
        FF (a, b, c, d, X[ 0], MD5_S11, MD5_T01); 
        FF (d, a, b, c, X[ 1], MD5_S12, MD5_T02); 
        FF (c, d, a, b, X[ 2], MD5_S13, MD5_T03); 
        FF (b, c, d, a, X[ 3], MD5_S14, MD5_T04); 
        FF (a, b, c, d, X[ 4], MD5_S11, MD5_T05); 
        FF (d, a, b, c, X[ 5], MD5_S12, MD5_T06); 
        FF (c, d, a, b, X[ 6], MD5_S13, MD5_T07); 
        FF (b, c, d, a, X[ 7], MD5_S14, MD5_T08); 
        FF (a, b, c, d, X[ 8], MD5_S11, MD5_T09); 
        FF (d, a, b, c, X[ 9], MD5_S12, MD5_T10); 
        FF (c, d, a, b, X[10], MD5_S13, MD5_T11); 
        FF (b, c, d, a, X[11], MD5_S14, MD5_T12); 
        FF (a, b, c, d, X[12], MD5_S11, MD5_T13); 
        FF (d, a, b, c, X[13], MD5_S12, MD5_T14); 
        FF (c, d, a, b, X[14], MD5_S13, MD5_T15); 
        FF (b, c, d, a, X[15], MD5_S14, MD5_T16); 
 
        //Runde 2 Transformation
        GG (a, b, c, d, X[ 1], MD5_S21, MD5_T17); 
        GG (d, a, b, c, X[ 6], MD5_S22, MD5_T18); 
        GG (c, d, a, b, X[11], MD5_S23, MD5_T19); 
        GG (b, c, d, a, X[ 0], MD5_S24, MD5_T20); 
        GG (a, b, c, d, X[ 5], MD5_S21, MD5_T21); 
        GG (d, a, b, c, X[10], MD5_S22, MD5_T22); 
        GG (c, d, a, b, X[15], MD5_S23, MD5_T23); 
        GG (b, c, d, a, X[ 4], MD5_S24, MD5_T24); 
        GG (a, b, c, d, X[ 9], MD5_S21, MD5_T25); 
        GG (d, a, b, c, X[14], MD5_S22, MD5_T26); 
        GG (c, d, a, b, X[ 3], MD5_S23, MD5_T27); 
        GG (b, c, d, a, X[ 8], MD5_S24, MD5_T28); 
        GG (a, b, c, d, X[13], MD5_S21, MD5_T29); 
        GG (d, a, b, c, X[ 2], MD5_S22, MD5_T30); 
        GG (c, d, a, b, X[ 7], MD5_S23, MD5_T31); 
        GG (b, c, d, a, X[12], MD5_S24, MD5_T32); 
 
        //Runde 3 Transformation
        HH (a, b, c, d, X[ 5], MD5_S31, MD5_T33); 
        HH (d, a, b, c, X[ 8], MD5_S32, MD5_T34); 
        HH (c, d, a, b, X[11], MD5_S33, MD5_T35); 
        HH (b, c, d, a, X[14], MD5_S34, MD5_T36); 
        HH (a, b, c, d, X[ 1], MD5_S31, MD5_T37); 
        HH (d, a, b, c, X[ 4], MD5_S32, MD5_T38); 
        HH (c, d, a, b, X[ 7], MD5_S33, MD5_T39); 
        HH (b, c, d, a, X[10], MD5_S34, MD5_T40); 
        HH (a, b, c, d, X[13], MD5_S31, MD5_T41); 
        HH (d, a, b, c, X[ 0], MD5_S32, MD5_T42); 
        HH (c, d, a, b, X[ 3], MD5_S33, MD5_T43); 
        HH (b, c, d, a, X[ 6], MD5_S34, MD5_T44); 
        HH (a, b, c, d, X[ 9], MD5_S31, MD5_T45); 
        HH (d, a, b, c, X[12], MD5_S32, MD5_T46); 
        HH (c, d, a, b, X[15], MD5_S33, MD5_T47); 
        HH (b, c, d, a, X[ 2], MD5_S34, MD5_T48); 
 
        //Runde 4 Transformation
        II (a, b, c, d, X[ 0], MD5_S41, MD5_T49); 
        II (d, a, b, c, X[ 7], MD5_S42, MD5_T50); 
        II (c, d, a, b, X[14], MD5_S43, MD5_T51); 
        II (b, c, d, a, X[ 5], MD5_S44, MD5_T52); 
        II (a, b, c, d, X[12], MD5_S41, MD5_T53); 
        II (d, a, b, c, X[ 3], MD5_S42, MD5_T54); 
        II (c, d, a, b, X[10], MD5_S43, MD5_T55); 
        II (b, c, d, a, X[ 1], MD5_S44, MD5_T56); 
        II (a, b, c, d, X[ 8], MD5_S41, MD5_T57); 
        II (d, a, b, c, X[15], MD5_S42, MD5_T58); 
        II (c, d, a, b, X[ 6], MD5_S43, MD5_T59); 
        II (b, c, d, a, X[13], MD5_S44, MD5_T60); 
        II (a, b, c, d, X[ 4], MD5_S41, MD5_T61); 
        II (d, a, b, c, X[11], MD5_S42, MD5_T62); 
        II (c, d, a, b, X[ 2], MD5_S43, MD5_T63); 
        II (b, c, d, a, X[ 9], MD5_S44, MD5_T64); 
 
        //Füge die veränderten Werte zur aktuellen Checksum hinzu
        m_lMD5[0] += a;
        m_lMD5[1] += b;
        m_lMD5[2] += c;
        m_lMD5[3] += d;
}
 
 
/*****************************************************************************************
CONSTRUCTOR:     SC_MD5
BESCHREIBUNG:    Initialisiert die Member-Variablen
ARGUMENTE:       Keine
BEMERKUNG:       Keine
*****************************************************************************************/
SC_MD5::SC_MD5()
{
        // Initialisierung
        memset( m_lpszBuffer, 0, 64 );
        m_nCount[0] = m_nCount[1] = 0;
 
        // Lade magic state Initialisierungskonstanten
        m_lMD5[0] = MD5_INIT_STATE_0;
        m_lMD5[1] = MD5_INIT_STATE_1;
        m_lMD5[2] = MD5_INIT_STATE_2;
        m_lMD5[3] = MD5_INIT_STATE_3;
}
 
/*****************************************************************************************
FUNKTION:              SC_MD5::DWordToByte
DETAILS:               private
BESCHREIBUNG:          Transferiert die Daten eines 32 bit array in ein 8 bit array
RÜCKGABEWERTE:         void
ARGUMENTE:             BYTE* Output  : 8 bit Zielarray 
                       DWORD* Input  : 32 bit Quellarray
                       UINT nLength  : Anzahl an 8 bit Dateneinheiten im Quellarray
BEMERKUNG:             Ein DWORD aus dem Quellarray wird in vier Bytes im Zielarray transferiert. 
                        
                                
                       Der Algorithmus prüft, ob das Zielarray ein vielfaches von 4 bytes long ist,
                       sodass 8 bit BYTES perfekt in 32 bit DWORDs passen.
*****************************************************************************************/
void SC_MD5::DWordToByte(BYTE* Output, DWORD* Input, UINT nLength, int& error)
{
         //entry invariants
        if( nLength % 4 != 0 ) error = -22;
 
        //transferiere die Daten durch shifting und kopieren
        UINT i = 0;
        UINT j = 0;
        for ( ; j < nLength; i++, j += 4) 
        {
                Output[j] =   (BYTE)(Input[i] & 0xff);
                Output[j+1] = (BYTE)((Input[i] >> 8) & 0xff);
                Output[j+2] = (BYTE)((Input[i] >> 16) & 0xff);
                Output[j+3] = (BYTE)((Input[i] >> 24) & 0xff);
        }
}
 
 
/*****************************************************************************************
FUNKTION:              SC_MD5::Final
DETAILS:               protected
BESCHREIBUNG:          Implementierung des Haupt MD5 Algorithmus; Beendet die Checksum-Berechnung.
RÜCKGABEWERTE:         char* : Der endgültige hexadecimale MD5 checksum Wert (changed from CString to char* by thilo to avoid MFC)
ARGUMENTE:             None
BEMERKUNG:             Führt die finale MD5 checksum Berechnung durch('Update' erledigt die Hauptarbeit,
                       diese Funktion beendet nur die Berechnung.) 
*****************************************************************************************/
char* SC_MD5::Final(int& error)
{
        //Sichere Anzahl der Bits
        BYTE Bits[8];
        DWordToByte( Bits, m_nCount, 8, error);
 
        //Auffüllen bis 56 mod 64.
        UINT nIndex = (UINT)((m_nCount[0] >> 3) & 0x3f);
        UINT nPadLen = (nIndex < 56) ? (56 - nIndex) : (120 - nIndex);
        Update( PADDING, nPadLen, error );
 
        //Füge die Länge hinzu (vorm Auffüllen)
        Update( Bits, 8, error );
 
        //Sichere final state in 'lpszMD5'
        const int nMD5Size = 16;
        unsigned char lpszMD5[ nMD5Size ];
        DWordToByte( lpszMD5, m_lMD5, nMD5Size,error );
 
				/*
        //Konvertiere die hexadecimale Checksum in einen string
        CString strMD5;
        for ( int i=0; i < nMD5Size; i++) 
        {
                CString Str;
                if (lpszMD5[i] == 0) 
                {
                        Str = "00";
                }
                else if (lpszMD5[i] <= 15)      
                {
                        Str.Format("0%x",lpszMD5[i]);
                }
                else 
                {
                        Str.Format("%x",lpszMD5[i]);
                }
 
                if(Str.GetLength() != 2) 
                {
                        error = -20;
                        break;
                }
                strMD5 += Str;
        }
        if(strMD5.GetLength() != 32) error = -21;
        return strMD5;
				*/

				//by thilo:
				char *checksum = new char[sclib::bufferSize];
				char buffer[sclib::bufferSize];
				int len, count = 0;
				for (int i = 0; i < nMD5Size; i++) {
					if (lpszMD5[i] == 0) {
	          sprintf(buffer, "00");
          } else if (lpszMD5[i] <= 15) {
            sprintf(buffer, "0%x", lpszMD5[i]);
          } else {
            sprintf(buffer, "%x", lpszMD5[i]);
          } 
					len = (int)(strnlen(buffer, sclib::bufferSize));
					if(len != 2) {
            error = -20;
            break;
          }
					sprintf(checksum+count, "%s", buffer);
					count += len;
				}
				if(strnlen(checksum, sclib::bufferSize) != 32) error = -21;
				return checksum;
}
 
 
/*****************************************************************************************
FUNKTION:              SC_MD5::Update
DETAILS:               protected
BESCHREIBUNG:          Implementierung des Haupt MD5 Algorithmus
RÜCKGABEWERTE:         void
ARGUMENTE:             BYTE* Input    : Inputblock
                       UINT nInputLen : Länge des Inputblocks
BEMERKUNG:             Berechnet die Teil MD5 checksum von 'nInputLen' Bytes an Daten in 'Input'
*****************************************************************************************/
void SC_MD5::Update( BYTE* Input, ULONG nInputLen, int& error )
{
        //Anzahl der Bytes mod 64 berechnen
        UINT nIndex = (UINT)((m_nCount[0] >> 3) & 0x3F);
 
        //Anzahl der Bits aktualisieren
        if ((m_nCount[0] += nInputLen << 3) < (nInputLen << 3))
        {
                m_nCount[1]++;
        }
        m_nCount[1] += (nInputLen >> 29);
 
        //Transformiere so oft wie möglich
        UINT i=0;              
        UINT nPartLen = 64 - nIndex;
        if (nInputLen >= nPartLen)      
        {
                memmove( &m_lpszBuffer[nIndex], Input, nPartLen );
                
                Transform( m_lpszBuffer, error );
                for (i = nPartLen; i + 63 < nInputLen; i += 64) 
                {
                        Transform( &Input[i], error );
                }
                nIndex = 0;
        } 
        else 
        {
                i = 0;
        }
 
        // Übrig gebliebenen Input in Buffer
        memmove( &m_lpszBuffer[nIndex], &Input[i], nInputLen-i);
}
