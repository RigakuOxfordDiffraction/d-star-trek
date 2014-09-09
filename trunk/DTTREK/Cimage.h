#ifndef DT_CIMAGE_H
#define DT_CIMAGE_H
//
// Copyright (c) 1995 Molecular Structure Corporation
//
// RESTRICTED RIGHTS NOTICE SHORT FORM (JUNE 1987)
//
// Use, reproduction, or disclosure is subject to restrictions set
// forth in Contract No. W-31-109-ENG-38 and Contract No.
// 943072401 with the University of Chicago, Operator of
// Argonne National Laboratory.
//
// Cimage.h        Initial author: J.W. Pflugrath           03-Mar-1995
//    This file is the header file for class Cimage.
/*
 *
 * Copyright (C) 2014 Rigaku Americas Corporation
 *                    9009 New Trails Drive
 *                    The Woodlands, TX, USA  77381
 *
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *    * Redistributions of source code must retain the above copyright
 *      notice(s), this list of conditions and the following disclaimer.
 *    * Redistributions in binary form must reproduce the above copyright
 *      notice(s), this list of conditions and the following disclaimer in the
 *      documentation and/or other materials provided with the distribution.
 *    * Neither the name of the Rigaku Americas Corporation nor the 
 *      names of its contributors may be used to endorse or promote products
 *      derived from this software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED. IN NO EVENT SHALL RIGAKU AMERICAS CORPORATION BE LIABLE 
 * FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
 * DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
 * SERVICES; LOSS OF USE, DATA OR PROFITS; OR BUSINESS INTERUPTION) HOWEVER
 * CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
 * LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY
 * OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
 * DAMAGE.
 *
*/

//+Include files

#include "Dtrek.h"
#include "dtrekdefs.h"
#include "Cstring.h"
#include "Cimage_header.h"
#include "C3Ddata.h"
#include "Crefln.h"
#include "dskio.h"
#include "dtreksys.h"
#include "CGonioMask.h"

//+External prototypes

//+Definitions and constants
// For safety on the enum have the unknown state first.  This way if things
//  are initialized to 0, then the program will sort of know.

enum eImage_states {
  eImage_unknown_state,
  eImage_scrambled_state,
  eImage_unscrambled_state,
  eImage_available_state
};

enum eImage_orientations {
  eO_unknown,
  eO12,
  eO_12,
  eO1_2,
  eO_1_2,
  eO21,
  eO_21,
  eO2_1,
  eO_2_1
};

enum eImage_compression_types {
  eImage_other_compression,
  eImage_no_compression,
  eImage_MSC1_compression,
  eImage_BS1_compression,
  eImage_PCK1_compression,
  eImage_CBF_compression,
  eImage_BAS2000_compression,
  eImage_HIPIC_compression
};

typedef struct _tagPIXEL_INTENSITY_BIN
{
    int     m_nBegin;
    int     m_nEnd;
    int     m_nCount;
    
    _tagPIXEL_INTENSITY_BIN(int nBeg, int nEnd){m_nBegin=nBeg; m_nEnd=nEnd; m_nCount = 0;}
}PIXEL_INTENSITY_BIN;


//+Code begin

//class Crotation;  // Forward declaration of class Crotation
class C3Ddata;      // Forward declaration of class C3Ddata

class DTREK_EXPORT Cimage {

public:

  eImage_states            m_eThe_State;
  int                      m_nDim[2];
  eImage_data_types        m_eData_type;
  eImage_compression_types m_eCompression;
//
  int                m_nSequence_number;
  Cstring            m_sFilename;
  Cimage_header      m_oHeader;
//
  float              m_fSatVal;
  float              m_fMinRawPixOKValue;
  int                m_nData_size;    // Size needed for data in bytes
  int                m_nBytesAllocated;
  int                m_nDetNum;
  int                m_nPixValOffset;

  eCPU_byte_orders   m_eImgByteOrder;

  union {                          // Pointer to various data types
    short int           *pi;
    unsigned short int  *pui;
    int                 *pn;
    unsigned int        *pun;
    LONG                *pl;
    ULONG               *pul;
    char                *pc;
    unsigned char       *puc;
    float               *pf;
    void                *pv;
  } m_The_Data;

  union {                          // Pointer to the next pixel
    short int           *pi;
    unsigned short int  *pui;
    int                 *pn;
    unsigned int        *pun;
    LONG                *pl;
    ULONG               *pul;
    char                *pc;
    unsigned char       *puc;
    float               *pf;
    void                *pv;
  } m_Next_Pixel;

  int    m_nOverflowCount;                  // Num values in overflow table
  int   *m_pnOverflowOffset;                // Offset of overflow in data array
  float *m_pfOverflowValue;                 // Value of overflow

  unsigned short int *m_puiBitmapMask;
  int    m_nBitmapMaskByteCount;
  char*  m_pcHeader;

  CGonioMask*            m_poGonioMask;                     // Hopefully, a temporary member.  We would rather get this info from the embedded mask.
  bool                           m_bReadApplyGonioMask;         // Defaults to false.  We don't want to mess up others reading the images for display purposes.
  bool                           m_bReadApplyEmbeddedMask;  // Defaults to false.  We don't want to mess up others reading the images for display purposes.
  bool               m_bWriteRAXIS;             // Defaults to true.   Referenced when calling nWrite().

  Cstring            m_sDescription;
//
  Cstring            m_sScan_key;
  int                m_nSelection[10];

  float              m_fCompressInfo[5]; // Five floats for various uses

  eImage_orientations m_eOrientation;

  static int        ms_nVerbose;
  static int        ms_nPilatusPixelValueOffset;

// Special for speed: pointers to member functions.  The syntax of this
// is tricky, so watch out.

  float  (Cimage::*prfGetPixel)(const int nWhere0, const int nWhere1);
  float  (Cimage::*prfCvtUItoFloat)(const unsigned short int uiValue);
  int    (Cimage::*prnSetPixel)(const int nWhere0, const int nWhere1,
                                const float fValue);

// (class) Cscan              *poScan;
// (class) Crotation          *poRotation;
// (class) Cnonuniformity     *poNonunf;
// (class) Cspatial           *poSpatial;

////////////////////////////////////////////////////////////////////////
//+Constructors, destructors and assignments

Cimage ();                              // Construct an empty image

// Construct an empty image of this size

Cimage (const int nDim0, const int nDim1);

// Construct an empty image of this size  and datatype

Cimage (const int nDim0, const int nDim1, const eImage_data_types eD, const char* sName = "");

// Construct an image of this size and datatype with *Data

Cimage (const int nDim0, const int nDim1, const eImage_data_types eD,
        void *pvData);


Cimage (const Cstring& sFilename);      // Construct an image by reading from
                                        //    file sFilename

// Construct an image by reading from a file whose name is constructed
//    from a template and sequence num

Cimage (const Cstring& sTemplatename, const int nSeq);

// Construct an image from a C3Ddata object and a reflection

Cimage(C3Ddata *poShoebox, const int nDir=0, Crefln *poRefln=NULL);

// Cimage (Detector d);                 // Construct an image from Detector
                                        //    hardware
// Cimage (Detector d, Rotation r);     // Construct an image from Detector
                                        //    hardware while using goniometer
Cimage (Cimage& oImageIn, void *pData = NULL); // Construct an image from
                                               // another image and a pointer
                                               // to the data
// Construct image form header string and pointer to data
Cimage (const Cstring& sHeaderString, const void *pData);

// This routine is the assignment operator, which sets this Cimage to be
// equal to the passed Cimage
Cimage &operator=(const Cimage &oEqual);

~Cimage ();

////////////////////////////////////////////////////////////////////////
//+Member Function prototypes
 public:

inline int nGetDimension(const int nWhich) { return (m_nDim[nWhich]); }
inline int nGetDimension(void)  { return (m_nDim[0] * m_nDim[1]); }

       int nGetDimensions(int *pnDim0, int *pnDim1);

inline float fGetSatValue(void) { return (m_fSatVal); }
inline void  vSetSatValue(const float fSatVal) { m_fSatVal = fSatVal; }
inline float fGetMinRawPixOKValue(void) { return (m_fMinRawPixOKValue); }
inline void  vSetMinRawPixOKValue(const float fMinRawPixOKValue) { m_fMinRawPixOKValue = fMinRawPixOKValue; }
inline int nGetSize(void) {return (nDataSizeNeeded() + m_oHeader.nLength()); }

inline void vSetState(const eImage_states eState)
  { m_eThe_State = eState; }

int nSetHeader(Cimage_header *poHeader);

//

int nGetPixel(const int nWhere[2], short int *piPixel);
int nGetPixel(const int nWhere0, const int nWhere1, short int *piPixel);
int nGetPixel(const int nWhere0, const int nWhere1, unsigned short int *puiPixel);

int nGetPixel(const int nWhere[2], int *pnPixel);
int nGetPixel(const int nWhere0, const int nWhere1, int *pnPixel);
int nGetPixel(const int nWhere0, const int nWhere1, unsigned int *punPixel);

float fGetuiPixel(const int nWhere0, const int nWhere1);
float fGetfPixel(const int nWhere0, const int nWhere1);
float fGetPixel(const int nWhere0, const int nWhere1);
float fGetRAXISPixel(const int nWhere0, const int nWhere1);
float fGetCBFPixel(const int nWhere0, const int nWhere1);
float fGetDIP2030Pixel(const int nWhere0, const int nWhere1);
float fGetSiemensPixel(const int nWhere0, const int nWhere1);

float fCvtUItoFloat(const unsigned short int uiValue);
float fCvtRAXIStoFloat(const unsigned short int uiValue);
float fCvtCBFtoFloat(const unsigned short int uiValue);
float fCvtCBFtoFloat(const long int nValue);
float fCvtDIP2030toFloat(const unsigned short int uiValue);
float fCvtSiemenstoFloat(const unsigned short int uiValue);

int nSetNextPixel(const int nWhere1, const int Where2);
int nSetNextPixel(const int nWhere[2]);

inline Cstring sGetName(void) { return (m_sFilename); }

inline short int iGetNextPixel(void)
{   return (*m_Next_Pixel.pi++); }

inline unsigned short int uiGetNextPixel(void)
{   return (*m_Next_Pixel.pui++); }

inline unsigned short int uiGetNextPixelNoInc(void)
{   return (*m_Next_Pixel.pui); }

inline short int iGetNextPixelNoInc(void)
{   return (*m_Next_Pixel.pi); }

inline float fGetNextPixelNoInc(void)
{   return (*m_Next_Pixel.pf); }

inline float fGetNextPixel(void)
{   return (*m_Next_Pixel.pf++); }


// A few special inline functions for faster access in Ccalibrate

inline float fGetPixel(const int nOffset)
{
  return (*(m_The_Data.pf + nOffset)); }

#ifndef OSF1

int nGetPixel(const int nWhere0, const int nWhere1, LONG *plPixel);

inline void vSetPixel(const int nOffset, const LONG lPixel)
    { *(m_The_Data.pl + nOffset) = lPixel; }

inline void vSetNextPixel(const LONG lPixel)
{  *m_Next_Pixel.pl++ = lPixel; }

inline void vSetNextPixelNoInc(const LONG lPixel)
{  *m_Next_Pixel.pl = lPixel; }


#endif

inline LONG lGetNextPixel(void)
{   return (*m_Next_Pixel.pl++); }

inline LONG lGetPixel(const int nOffset)
{
  return (*(m_The_Data.pl + nOffset)); }

inline unsigned short int uiGetPixel(const int nOffset)
{
  return (*(m_The_Data.pui + nOffset)); }


inline void vSetPixel(const int nOffset, const unsigned short int uiPixel)
    { *(m_The_Data.pui + nOffset) = uiPixel; }


inline void vSetNextPixel(const short int iPixel)
{  *m_Next_Pixel.pi++ = iPixel; }

inline void vSetNextPixel(const unsigned short int uiPixel)
{  *m_Next_Pixel.pui++ = uiPixel; }

inline void vSetNextPixel(const unsigned short int *puiPixel)
{  *m_Next_Pixel.pui++ = *puiPixel; }

inline void vSetNextPixelNoInc(const short int iPixel)
{  *m_Next_Pixel.pi = iPixel; }

inline void vSetNextPixel(const int nPixel)
{  *m_Next_Pixel.pn++ = nPixel; }

inline void vSetNextPixelNoInc(const int nPixel)
{  *m_Next_Pixel.pn = nPixel; }

inline void vSetNextPixel(const float fPixel)
{  *m_Next_Pixel.pf++ = fPixel; }

inline void vSetNextPixelNoInc(const float fPixel)
{  *m_Next_Pixel.pf = fPixel; }

int nGetRect(const int nRect[4], char **paData);

int nSetPixel(const int nWhere[2], const short int iPixel);
int nSetPixel(const int nWhere0, const int nWhere1,
              const unsigned short int uiPixel);
int nSetPixel(const int nWhere0, const int nWhere1,
              const short int iPixel);
int nSetPixel(const int nWhere0, const int nWhere1,
              const float fPixel);
int nSetCBFPixel(const int nWhere0, const int nWhere1,
              const float fPixel);

int nSetPixel(const int nWhere0, const int nWhere1,
              const ULONG ulPixel);
int nSetPixel(const int nWhere0, const int nWhere1,
              const LONG lPixel);

void *pvGetData(void) {return m_The_Data.pv; }
Cimage_header oGetHeader(void) {return m_oHeader; }

inline static void vSetVerboseLevel(const int nLev) { ms_nVerbose = nLev; }
int nPutRect(int nArray);

int nCompress (void);
int nUncompress (void);
int nUnscramble (void);

int nConvertOrientation (void);
int nConvertIntToUIShort(const int nZeroConversion=1);

int nRead(void);
int nRead(const Cstring& sName);
int nRead(const Cstring& sTemplate, const int nNumber);

int nWrite(void);
int nWrite(const Cstring& sName, bool bMustWriteDTrekHeader=false);
//int nWrite(const Cstring& sTemplate, const int nNumber);

int nList(const int nRect[4]);
int nComputeApplyActiveMask();  

// Copy a Rect from another image
int nCopyRect(const Cimage& hFrom, const int nFrom[4], const int nTo[2]);

int nUnbin(void);

inline bool bIsAvailable(void)
{  return (eImage_available_state == m_eThe_State); }

inline eImage_states eGetState(void)
{  return (m_eThe_State); }

int nAvgSD(const int nPxArea[4],
           const float fExcludeMin,
           const float fExcludeMax,
           float *pfAvg, float *pfSD,
           float *pfMin, float *pfMax);

int    nSwapByteOrder(void);
int    nBytesPerPixel(void) const;
int    nDataSizeNeeded(void);
int    nGetDataType(void);
inline char *pcTheData(void) { return (m_The_Data.pc); }

Cstring sGetDataType(void);
inline int nGetDetNum(void) { return (m_nDetNum); }
inline int nGetBitmapMaskByteCount(void) { return (m_nBitmapMaskByteCount); }

private:
int     nInitValues(void);
int     m_nRef;
public:
int     nAddRef(){return ++m_nRef;}
int     nReleaseRef(){return --m_nRef;}
int     nGetRefCount()const{return m_nRef;};

#define DTREK_IMAGE_BIH_BITTEST             0x0001
#define DTREK_IMAGE_BIH_CUSTOM              0x0002

#define DTREK_IMAGE_BIH_RAXIS_II            0x0004
#define DTREK_IMAGE_BIH_CCD                 0x0008


#define PMT_OVERLAP_BEGIN                   8192
#define PMT_OVERLAP_END                     16383

void    vBuildIntensityHistogram(std::vector<_tagPIXEL_INTENSITY_BIN>& vecBins, 
                                 std::vector<int>& vecBitBins, 
                                 int& nNumberOfPixValuesInOverlapRange_PMT1,
                                 int& nNumberOfPixValuesInOverlapRange_PMT2,
                                 DTREK_WORD wCtrl);

eImage_format_types enGetOriginalFormat(){ return m_oHeader.enGetOriginalFormat(); }

};  // end of class Cimage

#endif   // DT_CIMAGE_H
