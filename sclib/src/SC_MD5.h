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

#ifndef __SC_MD5_H__
#define __SC_MD5_H__

//Initialisierungskonstanten
#define MD5_INIT_STATE_0 0x67452301
#define MD5_INIT_STATE_1 0xefcdab89
#define MD5_INIT_STATE_2 0x98badcfe
#define MD5_INIT_STATE_3 0x10325476
 
//Konstanten für den Transformierungsprozess.
#define MD5_S11  7
#define MD5_S12 12
#define MD5_S13 17
#define MD5_S14 22
#define MD5_S21  5
#define MD5_S22  9
#define MD5_S23 14
#define MD5_S24 20
#define MD5_S31  4
#define MD5_S32 11
#define MD5_S33 16
#define MD5_S34 23
#define MD5_S41  6
#define MD5_S42 10
#define MD5_S43 15
#define MD5_S44 21
 
//Transformierungskonstanten - Runde 1
#define MD5_T01  0xd76aa478 //Transformationskonstante 1 
#define MD5_T02  0xe8c7b756 //Transformationskonstante 2
#define MD5_T03  0x242070db //Transformationskonstante 3
#define MD5_T04  0xc1bdceee //Transformationskonstante 4
#define MD5_T05  0xf57c0faf //Transformationskonstante 5
#define MD5_T06  0x4787c62a //Transformationskonstante 6
#define MD5_T07  0xa8304613 //Transformationskonstante 7
#define MD5_T08  0xfd469501 //Transformationskonstante 8
#define MD5_T09  0x698098d8 //Transformationskonstante 9
#define MD5_T10  0x8b44f7af //Transformationskonstante 10
#define MD5_T11  0xffff5bb1 //Transformationskonstante 11
#define MD5_T12  0x895cd7be //Transformationskonstante 12
#define MD5_T13  0x6b901122 //Transformationskonstante 13
#define MD5_T14  0xfd987193 //Transformationskonstante 14
#define MD5_T15  0xa679438e //Transformationskonstante 15
#define MD5_T16  0x49b40821 //Transformationskonstante 16
 
//Transformierungskonstanten - Runde 2
#define MD5_T17  0xf61e2562 //Transformationskonstante 17
#define MD5_T18  0xc040b340 //Transformationskonstante 18
#define MD5_T19  0x265e5a51 //Transformationskonstante 19
#define MD5_T20  0xe9b6c7aa //Transformationskonstante 20
#define MD5_T21  0xd62f105d //Transformationskonstante 21
#define MD5_T22  0x02441453 //Transformationskonstante 22
#define MD5_T23  0xd8a1e681 //Transformationskonstante 23
#define MD5_T24  0xe7d3fbc8 //Transformationskonstante 24
#define MD5_T25  0x21e1cde6 //Transformationskonstante 25
#define MD5_T26  0xc33707d6 //Transformationskonstante 26
#define MD5_T27  0xf4d50d87 //Transformationskonstante 27
#define MD5_T28  0x455a14ed //Transformationskonstante 28
#define MD5_T29  0xa9e3e905 //Transformationskonstante 29
#define MD5_T30  0xfcefa3f8 //Transformationskonstante 30
#define MD5_T31  0x676f02d9 //Transformationskonstante 31
#define MD5_T32  0x8d2a4c8a //Transformationskonstante 32
 
//Transformierungskonstanten - Runde 3
#define MD5_T33  0xfffa3942 //Transformationskonstante 33
#define MD5_T34  0x8771f681 //Transformationskonstante 34
#define MD5_T35  0x6d9d6122 //Transformationskonstante 35
#define MD5_T36  0xfde5380c //Transformationskonstante 36
#define MD5_T37  0xa4beea44 //Transformationskonstante 37
#define MD5_T38  0x4bdecfa9 //Transformationskonstante 38
#define MD5_T39  0xf6bb4b60 //Transformationskonstante 39
#define MD5_T40  0xbebfbc70 //Transformationskonstante 40
#define MD5_T41  0x289b7ec6 //Transformationskonstante 41
#define MD5_T42  0xeaa127fa //Transformationskonstante 42
#define MD5_T43  0xd4ef3085 //Transformationskonstante 43
#define MD5_T44  0x04881d05 //Transformationskonstante 44
#define MD5_T45  0xd9d4d039 //Transformationskonstante 45
#define MD5_T46  0xe6db99e5 //Transformationskonstante 46
#define MD5_T47  0x1fa27cf8 //Transformationskonstante 47
#define MD5_T48  0xc4ac5665 //Transformationskonstante 48
 
//Transformierungskonstanten - Runde 4
#define MD5_T49  0xf4292244 //Transformationskonstante 49
#define MD5_T50  0x432aff97 //Transformationskonstante 50
#define MD5_T51  0xab9423a7 //Transformationskonstante 51
#define MD5_T52  0xfc93a039 //Transformationskonstante 52
#define MD5_T53  0x655b59c3 //Transformationskonstante 53
#define MD5_T54  0x8f0ccc92 //Transformationskonstante 54
#define MD5_T55  0xffeff47d //Transformationskonstante 55
#define MD5_T56  0x85845dd1 //Transformationskonstante 56
#define MD5_T57  0x6fa87e4f //Transformationskonstante 57
#define MD5_T58  0xfe2ce6e0 //Transformationskonstante 58
#define MD5_T59  0xa3014314 //Transformationskonstante 59
#define MD5_T60  0x4e0811a1 //Transformationskonstante 60
#define MD5_T61  0xf7537e82 //Transformationskonstante 61
#define MD5_T62  0xbd3af235 //Transformationskonstante 62
#define MD5_T63  0x2ad7d2bb //Transformationskonstante 63
#define MD5_T64  0xeb86d391 //Transformationskonstante 64
 
 
//Null Daten (ausser des ersten BYTE) werden benutzt, um die Prüfsummenberechnung zu beenden
static unsigned char PADDING[64] = {
  0x80, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
};

//by thilo:
typedef unsigned char BYTE;
typedef unsigned long DWORD; 
typedef unsigned int UINT;
typedef unsigned long ULONG;
 
/*****************************************************************************************/
class SC_MD5 {
public:
        //KonstruKtor/Destruktor
        SC_MD5();
        virtual ~SC_MD5() {};
 
        //RSA MD5 Implementierung
        void Transform(BYTE Block[64], int& error);
        void Update(BYTE* Input, ULONG nInputLen, int& error);
        char* Final(int& iErrorCalculate); //by thilo: changed to char* to avoid MFC
 
protected:
        inline DWORD RotateLeft(DWORD x, int n);
        inline void FF( DWORD& A, DWORD B, DWORD C, DWORD D, DWORD X, DWORD S, DWORD T);
        inline void GG( DWORD& A, DWORD B, DWORD C, DWORD D, DWORD X, DWORD S, DWORD T);
        inline void HH( DWORD& A, DWORD B, DWORD C, DWORD D, DWORD X, DWORD S, DWORD T);
        inline void II( DWORD& A, DWORD B, DWORD C, DWORD D, DWORD X, DWORD S, DWORD T);
 
        //Nebenfunktionen
        inline void DWordToByte(BYTE* Output, DWORD* Input, UINT nLength, int& error);
        inline void ByteToDWord(DWORD* Output, BYTE* Input, UINT nLength, int& error);
        void MemoryMove(BYTE* from, BYTE* to, UINT size);
 
private:
        BYTE  m_lpszBuffer[64];   //Eingabepuffer
        ULONG m_nCount[2] ;        //Anzahl der bits, modulo 2^64 (lsb zuerst)
        ULONG m_lMD5[4] ;          //MD5 Prüfsumme
};

#endif
