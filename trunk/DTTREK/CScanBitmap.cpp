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
// CScanBitmap.cc            Initial author: Thaddeus Niemeyer      1-17-2002
//  This file contains the member functions of class Cfind which implements
//    the spot finding encapsulation of d*TREK.
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

#include "CScanBitmap.h"
#include <string.h>
#include "dtsvd.h"



CScanBitmap::CScanBitmap() {
    vInit();
};

CScanBitmap::~CScanBitmap() {
    nCloseOutput();
    nCloseInput();
    if (m_pcLocalBitmap)
        delete[] m_pcLocalBitmap;
    m_pcLocalBitmap = NULL;
};


void  CScanBitmap::vInit() {
    m_a2nImageDim[0] = 0;
    m_a2nImageDim[1] = 0;
    m_nDataBits = 0;
    m_a2nScanRange[0] = -1;
    m_a2nScanRange[1] = -1;
    m_nDataBitMask = 0;   
    m_pFBitmapDataIn = NULL;
    m_pFBitmapDataOut = NULL;
    m_nInImagePointer = -1;
    m_nOutImagePointer = -1;
    m_pcLocalBitmap = NULL;
    m_nLocalBitmapSize = 0;
    strcpy(m_pcHeaderString,CSCAN_BITMAP_HEADERSTRING);
    -m_anBitmapStart;
};

int CScanBitmap::nUpdateHeader(Cimage_header& oHeader) {
    return oHeader.nReplaceValue((Cstring) "BITMAP_IMAGE_FILE", m_sFileName);
};


int CScanBitmap::nCloseOutput() {
    int nStat;
    if (m_pFBitmapDataOut) {
        nStat = nWriteTrailerInfo();
        fclose(m_pFBitmapDataOut);
    } else
        nStat = 0;
    m_pFBitmapDataOut = NULL;
    return nStat;
};

int CScanBitmap::nCloseInput() {
    int nStat;
    if (m_pFBitmapDataIn) {
        fclose(m_pFBitmapDataIn);
    } else
        nStat = 0;
    m_pFBitmapDataIn = NULL;
    return nStat;
};


int  CScanBitmap::nConsistencyCheck() {
    int nx;

    if ((m_a2nScanRange[0]<0) || (m_a2nScanRange[1] < 0) || (m_a2nScanRange[0] > m_a2nScanRange[1]))
        return 1;
    if ((m_a2nImageDim[0]<=0) || (m_a2nImageDim[1] <= 0))
        return 1;
    if ((m_nDataBits != 1) && (m_nDataBits != 2) && (m_nDataBits != 3) && (m_nDataBits != 4))
        return 1;
    if (strcmp(m_pcHeaderString,CSCAN_BITMAP_HEADERSTRING))
        return 1;
    if (m_nDataBitMask != (1 << m_nDataBits) - 1)
        return 1;

    // Make the real to compressed.
    -m_anRealToCompressed;
    for (nx=0;nx<= m_nDataBitMask + 1; nx++) {
        if (nx == m_nTypicalValue)
            m_anRealToCompressed + (-1);    // This value is not mapped.
        else if (nx>m_nTypicalValue)
            m_anRealToCompressed + (nx - 1);
        else
            m_anRealToCompressed + (nx);
    };
    for (nx=0;nx<= m_nDataBitMask + 1; nx++) {
        if (m_anRealToCompressed[nx]!=-1)
            m_anCompressedToReal[m_anRealToCompressed[nx]] = nx;
    };
    return 0;
};


int CScanBitmap::nOpenOutput(Cstring& sOutput,int nDim0,int nDim1,int nFirstInScan,int nLastInScan,int nDataBits,int nTypicalValue) {
    int nx;

    if (nCloseOutput())
        return 1;
    m_pFBitmapDataOut = fopen(sOutput.string(),"wb");
    if (!m_pFBitmapDataOut) {
        printf("ERROR:  Could not open '%s' for output.\n",sOutput.string());
        return 1;
    };
    m_sFileName = sOutput;

    // If nDim0 is >0 use settings, Otherwise, use the settings that we already have.
    if (nDim0 > 0) {
        m_a2nImageDim[0] = nDim0;
        m_a2nImageDim[1] = nDim1;
        m_nDataBits = nDataBits;
        m_a2nScanRange[0] = nFirstInScan;
        m_a2nScanRange[1] = nLastInScan;
        m_nDataBitMask = (1 << m_nDataBits) - 1;
        m_nTypicalValue = nTypicalValue;
    }    
    if (nConsistencyCheck())
        return 1;
    m_nOutImagePointer = m_a2nScanRange[0];

    nx = ((m_a2nImageDim[0]*m_a2nImageDim[1]*m_nDataBits) >> 3) + 1;
    if ((!m_pcLocalBitmap) || (nx != m_nLocalBitmapSize)) {
        delete[] m_pcLocalBitmap;
        m_pcLocalBitmap = new unsigned char[m_nLocalBitmapSize = nx];
    };
    vClearLocalBitmap();
    

    return 0;
};

void CScanBitmap::vClearLocalBitmap(unsigned char* pcBitmapOther) {

    int nx,ny;
    unsigned char* pcBitmap = pcBitmapOther?(pcBitmapOther):m_pcLocalBitmap;

    for (nx=0,ny=0; (ny==0) || (ny % 8);nx++,ny += m_nDataBits)
        vSetLocalBitmap(nx,0,m_nTypicalValue,pcBitmap);

    ny /= 8;
    for (nx=0;nx<m_nLocalBitmapSize;nx++)
        pcBitmap[nx] = pcBitmap[nx % ny];
};

int CScanBitmap::nGetLocalBitmap(int nPix0,int nPix1,unsigned char* pcBitmapOther) {
    int nx; 
    unsigned char* pcBitmap = pcBitmapOther?(pcBitmapOther):m_pcLocalBitmap;
    return ((((int) (pcBitmap[((nx = (m_nDataBits*(nPix1*m_a2nImageDim[0] + nPix0))) >> 3)])) + (((int) pcBitmap[(nx >> 3)+1]) << 8)) >> (nx % 8)) & m_nDataBitMask; 
}


void CScanBitmap::vSetLocalBitmap(int nPix0,int nPix1,int nValue,unsigned char* pcBitmapOther) {
    int nx,ny,nz;
    unsigned char* pcBitmap = pcBitmapOther?(pcBitmapOther):m_pcLocalBitmap;
    
    if (nPix1 == -1)
        nx = (m_nDataBits*nPix0);
    else
        nx = (m_nDataBits*(nPix1*m_a2nImageDim[0] + nPix0));
    nz = nx % 8;
    ny = pcBitmap[nx >> 3] + (((int) pcBitmap[(nx >> 3) + 1]) << 8);
    ny -= (((ny >> nz) & m_nDataBitMask) << nz);
    ny += ((nValue & m_nDataBitMask) << nz);
    pcBitmap[nx >> 3] = (ny & 255);
    pcBitmap[(nx >> 3) + 1] = ny >> 8;
    return;
};

int CScanBitmap::nWriteTrailerInfo() {
    int nStat;
    int nx;
    unsigned char cChar;

    // We should have as many entries as total files that need to be written.
    if (m_anBitmapStart.size() != m_a2nScanRange[1] - m_a2nScanRange[0] + 1) {
        printf("WARNING:  Expecting %d files to be added.  Only found %d.\n",
            m_a2nScanRange[1] - m_a2nScanRange[0] + 1,
            m_anBitmapStart.size());
        nStat = 1;
    } else
        nStat = 0;
    m_anBitmapStart + ((int) ftell(m_pFBitmapDataOut));

    for (nx=0;nx<m_anBitmapStart.size(); nx++) {
        cChar = m_anBitmapStart[nx] & 255;
        fwrite(&cChar,1,1,m_pFBitmapDataOut);
        cChar = (m_anBitmapStart[nx] >> 8) & 255;
        fwrite(&cChar,1,1,m_pFBitmapDataOut);
        cChar = (m_anBitmapStart[nx] >> 16) & 255;
        fwrite(&cChar,1,1,m_pFBitmapDataOut);
        cChar = (m_anBitmapStart[nx] >> 24) & 255;
        fwrite(&cChar,1,1,m_pFBitmapDataOut);
    };

    return nStat;
};

int CScanBitmap::nWriteHeaderInfo() {
    unsigned char  pcBuffer[200];
    int            nBufferPointer = 0;
    
    strcpy((char*) pcBuffer,CSCAN_BITMAP_HEADERSTRING);
    nBufferPointer = strlen(CSCAN_BITMAP_HEADERSTRING) + 1;
    pcBuffer[nBufferPointer++] = m_a2nImageDim[0] & 255;
    pcBuffer[nBufferPointer++] = (m_a2nImageDim[0] >> 8) & 255;
    pcBuffer[nBufferPointer++] = m_a2nImageDim[1] & 255;
    pcBuffer[nBufferPointer++] = (m_a2nImageDim[1] >> 8) & 255;
    pcBuffer[nBufferPointer++] = m_a2nScanRange[0] & 255;
    pcBuffer[nBufferPointer++] = (m_a2nScanRange[0] >> 8) & 255;
    pcBuffer[nBufferPointer++] = m_a2nScanRange[1] & 255;
    pcBuffer[nBufferPointer++] = (m_a2nScanRange[1] >> 8) & 255;
    pcBuffer[nBufferPointer++] = m_nDataBits;
    pcBuffer[nBufferPointer++] = m_nDataBitMask;
    pcBuffer[nBufferPointer++] = m_nTypicalValue;
    fwrite(pcBuffer,1,nBufferPointer,m_pFBitmapDataOut);
    -m_anBitmapStart;
    return 0;
};



int CScanBitmap::nWriteNextBitmap() {
    int nMaxAtypical;
    int nTypical;       // Flips between 0 and 1. 1==typical mode, 0 == atypical mode.

    unsigned char* pcDataPoint;     // Pointer to local bitmap.
    int            nBufCount;       // Pointer for m_acIO
    char*          pcBuffer;        // Pointer to m_acIO.

    int   nDataPointOffset;
    int   nDataPoint;
    int   nThisValue;
    int   nNextValue;
    int   nValueLength;
    int   nx;
    

    if ((m_nOutImagePointer > m_a2nScanRange[1]) || (m_nOutImagePointer < m_a2nScanRange[0]))
        return 1;

    // Write header if this is the first image.
    if (m_nOutImagePointer == m_a2nScanRange[0]) {
        if (nWriteHeaderInfo())
            return 1;       
    };
    m_anBitmapStart + (int) ( ftell(m_pFBitmapDataOut));

    // Compute maximum typical and atypical values.
    nMaxAtypical = ((1 << (8 - m_nDataBits))-1) - 4;


    // Write the data.
    pcDataPoint = m_pcLocalBitmap;
    nDataPointOffset = 0;
    nDataPoint = 0;
    nTypical = 0;
    nNextValue = (*pcDataPoint) & m_nDataBitMask;
    m_acIOBuffer[1000];
    nBufCount = 0;
    pcBuffer = &m_acIOBuffer[nBufCount];

    while (nDataPoint!=m_nLocalBitmapSize-1) {

        // This loop gets the next value, and finds out how often it is repeated.
        nThisValue = nNextValue;
        nValueLength = 0;
        do {
            nDataPointOffset += m_nDataBits;
            nValueLength++;
            if (nDataPointOffset >= 8) {
                nDataPoint++;
                pcDataPoint++;
                nDataPointOffset -= 8;
                if (nDataPoint==m_nLocalBitmapSize-1)
                    break;
            };
            nNextValue = ((((int) (*pcDataPoint)) + (((int) (*(pcDataPoint+1))) << 8)) >> nDataPointOffset) & m_nDataBitMask;
        } while (nThisValue == nNextValue);

        if (nThisValue == m_nTypicalValue) {
            // This is a typical value.
            if (!nTypical) {
                // we are not in typical mode.  Write out a single 0 byte
                *(pcBuffer++) = (char) 0;
                nBufCount++;
                // Flip to typical.
                nTypical = ! nTypical;
            };
            nx = nValueLength;
            do {
                nValueLength = nx;                
                nx = nValueLength >> 7;
                *(pcBuffer++) = (char) ((( nValueLength & 127) + 128*(nx!=0)));
                nBufCount++;
            } while (nx);

        } else {
            // This is an atypical value.
            if (nTypical) {
                // we are not in atypical mode.  Write out a single 0 byte.
                *(pcBuffer++) = (char) 0;
                nBufCount++;
                // Flip to atypical.
                nTypical = ! nTypical;
            };           

            nx = (m_anRealToCompressed[nThisValue]) << (8 - m_nDataBits);
            if (nValueLength <= nMaxAtypical) {
                // One byte will do.
                *(pcBuffer++) = (char) ( nx + nValueLength);
                nBufCount++;
            } else if (nValueLength < (1 << 8)) {
                // Two bytes.  
                *(pcBuffer++) = (char) ( nx + nMaxAtypical+1);
                nBufCount++;
                *(pcBuffer++) = (char) ( nValueLength);                
                nBufCount++;
            } else if (nValueLength < (1 << 16)) {
                // Three bytes.  
                *(pcBuffer++) = (char) ( nx + nMaxAtypical+2);
                nBufCount++;
                *(pcBuffer++) = (char) ( nValueLength & 255);
                nBufCount++;
                *(pcBuffer++) = (char) ( (nValueLength >> 8) & 255);
                nBufCount++;
            } else if (nValueLength < (1 << 24)) {
                // Four bytes.  
                *(pcBuffer++) = (char) ( nx + nMaxAtypical+3);
                nBufCount++;
                *(pcBuffer++) = (char) ( (nValueLength) & 255);
                nBufCount++;
                *(pcBuffer++) = (char) ( (nValueLength >> 8) & 255);
                nBufCount++;
                *(pcBuffer++) = (char) ( (nValueLength >> 16) & 255);
                nBufCount++;
            } else {
                // Five bytes.  
                *(pcBuffer++) = (char) ( nx + nMaxAtypical+4);
                nBufCount++;
                *(pcBuffer++) = (char) ( (nValueLength) & 255);
                nBufCount++;
                *(pcBuffer++) = (char) ( (nValueLength >> 8) & 255);
                nBufCount++;
                *(pcBuffer++) = (char) ( (nValueLength >> 16) & 255);
                nBufCount++;
                *(pcBuffer++) = (char) ( (nValueLength >> 24) & 255);
                nBufCount++;
            };
        };

        // Make sure there is enough buffer space.
        if (m_acIOBuffer.size() - nBufCount <=20) {
            m_acIOBuffer.setsize(m_acIOBuffer.size()*2);
            pcBuffer = &m_acIOBuffer[nBufCount];
        };

        // Flip typical mode.
        nTypical = !nTypical;
    };
    if (nBufCount != fwrite(&m_acIOBuffer[0],1,nBufCount,m_pFBitmapDataOut))
        return 1;

    m_nOutImagePointer++;
    return 0;
};

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

int CScanBitmap::nOpenInput(Cimage_header& oHeader) {
    Cstring sFileName;
    if (oHeader.nGetValue((Cstring) "BITMAP_IMAGE_FILE", &sFileName)) {
        sFileName = sTransSymbol(sDtrekGetPrefix());
        sFileName += CSCAN_BITMAP_FILE;
    };
    return nOpenInput(sFileName);
};

int CScanBitmap::nOpenInput(Cstring& sInput) {
    int nx;

    if (nCloseInput())
        return 1;
    m_pFBitmapDataIn = fopen(sInput.string(),"rb");
    if (!m_pFBitmapDataIn) {
        printf("ERROR:  Could not open '%s' for input.\n",sInput.string());
        return 1;
    };
    m_sFileName = sInput;

    if (nReadHeaderInfo()) {
        nCloseInput();
        return 1;
    };

    nx = ((m_a2nImageDim[0]*m_a2nImageDim[1]*m_nDataBits) >> 3) + 1;
    if ((!m_pcLocalBitmap) || (nx != m_nLocalBitmapSize)) {
        delete[] m_pcLocalBitmap;
        m_pcLocalBitmap = new unsigned char[m_nLocalBitmapSize = nx];
        vClearLocalBitmap();
    };    
    m_nInImagePointer = m_a2nScanRange[0];
    return 0;
};


int CScanBitmap::nReadHeaderInfo() {
    unsigned char pcBuffer[200];
    int nBufferPointer;

    memset(pcBuffer,0,200);
    fread(pcBuffer,1,strlen(CSCAN_BITMAP_HEADERSTRING) + 1,m_pFBitmapDataIn);
    if (strcmp(CSCAN_BITMAP_HEADERSTRING,(char*) pcBuffer))
        return 1;
    if (11 != fread(pcBuffer,1,11,m_pFBitmapDataIn))
        return 1;
    nBufferPointer = 0;
    m_a2nImageDim[0] = pcBuffer[nBufferPointer++];
    m_a2nImageDim[0] += (((int) pcBuffer[nBufferPointer++]) << 8);
    m_a2nImageDim[1] = pcBuffer[nBufferPointer++];
    m_a2nImageDim[1] += (((int) pcBuffer[nBufferPointer++]) << 8);
    m_a2nScanRange[0] = pcBuffer[nBufferPointer++];
    m_a2nScanRange[0] += (((int) pcBuffer[nBufferPointer++]) << 8);
    m_a2nScanRange[1] = pcBuffer[nBufferPointer++];
    m_a2nScanRange[1] += (((int) pcBuffer[nBufferPointer++]) << 8);
    m_nDataBits = (int) pcBuffer[nBufferPointer++];
    m_nDataBitMask = (int) pcBuffer[nBufferPointer++];
    m_nTypicalValue = (int) pcBuffer[nBufferPointer++];
    if (nConsistencyCheck())
        return 1;

    return nReadTrailerInfo();
};

int CScanBitmap::nReadTrailerInfo() {
    long int nCurrent;
    int nToRead;
    int nx,ny;
    itr<char> acBuffer;
    unsigned char* pcBuffer;

    nCurrent = ftell(m_pFBitmapDataIn);
    nToRead = (m_a2nScanRange[1] - m_a2nScanRange[0] + 1 + 1)*4;
    acBuffer[nToRead];
    pcBuffer = (unsigned char*) &acBuffer[0];
    if (0 != fseek(m_pFBitmapDataIn,-nToRead,SEEK_END))
        return 1;
    if (nToRead != fread(pcBuffer,1,nToRead,m_pFBitmapDataIn)) 
        return 1;
    -m_anBitmapStart;
    for (nx=0;nx<nToRead;) {
        m_anBitmapStart[nx/4] = 0;
        for (ny=0;ny<4;ny++)
            m_anBitmapStart[nx/4] += (((int) pcBuffer[nx++]) << (ny*8));
    };
    for (nx=1;nx<m_anBitmapStart.size();nx++) {
        if (m_anBitmapStart[nx]<=m_anBitmapStart[nx-1])
            return 1;
    };
    if (0 != fseek(m_pFBitmapDataIn,nCurrent,SEEK_SET))
        return 1;
    return 0;
};

int CScanBitmap::nReadNextBitmap() {
    int nRegisteredSize;
    int nMaxDataPoint;
    unsigned char* pcDataPoint;
    int nDataPoint;
    int nTypical;
    int nMaxAtypical;
    int nRepeat;
    int nValue;
    int nx,ny;
    unsigned char* pcBuffer;
  

    if ((m_nInImagePointer > m_a2nScanRange[1]) || (m_nInImagePointer < m_a2nScanRange[0]))
        return 1;

    // The assumption is made that the FILE* read position is in sync. with m_nInImagePointer.
    // 
    nRegisteredSize = m_anBitmapStart[m_nInImagePointer - m_a2nScanRange[0] + 1 ] - m_anBitmapStart[m_nInImagePointer - m_a2nScanRange[0]];
    m_acIOBuffer[nRegisteredSize];

    // Read in all the data at once.
    if (nRegisteredSize != fread(&m_acIOBuffer[0],1,nRegisteredSize,m_pFBitmapDataIn))
        return 1;

    vClearLocalBitmap();

    // Uncompress the data.


    // Write the data.
    nMaxDataPoint = m_a2nImageDim[0]*m_a2nImageDim[1];
    pcDataPoint = m_pcLocalBitmap;
    nDataPoint = 0;
    nTypical = 0;
    pcBuffer = (unsigned char*) &m_acIOBuffer[0];
    nMaxAtypical = ((1 << (8 - m_nDataBits))-1) - 4;

    while (nDataPoint != nMaxDataPoint) {

        // This loop gets the next value, and finds out how often it is repeated.
        if (nTypical) {
            nRepeat = 0;
            ny = 0;
            do {
                nx = (int) *(pcBuffer++);
                nRepeat += (nx & 127) << (7*ny);
                ny ++;
            } while (nx & 128);
            nValue = m_nTypicalValue;           
        } else {
            nValue = (int) *(pcBuffer) >> (8 - m_nDataBits);
            nx = (int) *(pcBuffer) - (nValue << (8 - m_nDataBits));
            pcBuffer++;
            if (nValue >= m_anCompressedToReal.size())
                return 1;
            else
                nValue = m_anCompressedToReal[nValue];
            nRepeat = 0;
            if (nx <= nMaxAtypical) {
                nRepeat = nx;
            } else if (nx == nMaxAtypical+1) {
                nRepeat = ((int) *(pcBuffer++)); 
            } else if (nx == nMaxAtypical+2) {
                nRepeat = ((int) *(pcBuffer++)); 
                nRepeat += ((int) *(pcBuffer++)) << 8; 
            } else if (nx == nMaxAtypical+3) {
                nRepeat = ((int) *(pcBuffer++)); 
                nRepeat += ((int) *(pcBuffer++)) << 8; 
                nRepeat += ((int) *(pcBuffer++)) << 16; 
            } else {
                nRepeat = ((int) *(pcBuffer++)); 
                nRepeat += ((int) *(pcBuffer++)) << 8; 
                nRepeat += ((int) *(pcBuffer++)) << 16; 
                nRepeat += ((int) *(pcBuffer++)) << 24; 
            };
        };
       
        if (nValue != m_nTypicalValue) {
            for (nx = 0; nx < nRepeat; nx++) {
                vSetLocalBitmap(nDataPoint + nx,-1,nValue);
            };
        };
        nDataPoint += nRepeat;
        if (nDataPoint > nMaxDataPoint)
            break;

        // Flip typical mode.
        nTypical = !nTypical;
    };
    if (nDataPoint != nMaxDataPoint)
        return 1;

    m_nInImagePointer++;

    return 0;
};


/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

int CScanBitmap::nReadBitmapForScan(int nScanNumber) {
    if ((nScanNumber > m_a2nScanRange[1]) || (nScanNumber < m_a2nScanRange[0]))
        return 1;
    if (!m_pFBitmapDataIn)
        return 1;
    m_nInImagePointer = nScanNumber;
    if (fseek(m_pFBitmapDataIn,m_anBitmapStart[m_nInImagePointer - m_a2nScanRange[0]],SEEK_SET))
        return 1;
    return nReadNextBitmap();
};
int CScanBitmap::nReadBitmapForScan(Cscan& oScan) {
    return nReadBitmapForScan(oScan.nGetSeqNum());
};

int CScanBitmap::nWriteMaskedImage(Cstring& sImageName,Cimage* poOtherImage) {
    int nx,ny;
    int nValue;
    int nLowValue;
    int nHighValue;
    int a2nDelta[2];
    bool bOnEdge;

    if (!m_pcLocalBitmap)
        return 0;
    

    Cimage oImageTemp(m_a2nImageDim[0],m_a2nImageDim[1],(poOtherImage)?(poOtherImage->m_eData_type):eImage_uI2);

    if (poOtherImage) {
        nLowValue = 0;
    } else {
        nLowValue = 1;
        nHighValue = 0;
    }
    
    for (nx=0;nx<m_a2nImageDim[0];nx++) {
        for (ny=0;ny<m_a2nImageDim[1];ny++) {
            if (nValue = nGetLocalBitmap(nx,ny)) {
                bOnEdge = false;

                // Check around the box and see if there is a zero pixel.                            
                for (a2nDelta[0] = -1; (a2nDelta[0] <= 1) && (!bOnEdge); a2nDelta[0]++) {
                    for (a2nDelta[1] = -1; (a2nDelta[1] <= 1) && (!bOnEdge); a2nDelta[1]++) {
                        if ((a2nDelta[0] + nx >= 0) && (a2nDelta[0] + nx < m_a2nImageDim[0]) &&
                            (a2nDelta[1] + ny >= 0) && (a2nDelta[1] + ny < m_a2nImageDim[1])) {
                            if (nGetLocalBitmap(nx + a2nDelta[0],ny + a2nDelta[1])!=nValue)
                                bOnEdge = true;
                        } else
                            bOnEdge = true;
                    };
                }; 

            } else
                bOnEdge = false;

                                        
            if (bOnEdge) {
                oImageTemp.nSetPixel(nx,ny,(unsigned short int) nValue);
            } else {
                if (poOtherImage)
                    (oImageTemp.*oImageTemp.prnSetPixel)(nx,ny,(poOtherImage->*poOtherImage->prfGetPixel)(nx,ny));
                else
                    oImageTemp.nSetPixel(nx,ny,(unsigned short int) nHighValue);
            };
        };
    };
    if (poOtherImage) {
        oImageTemp.m_oHeader.nCopyMask(poOtherImage->m_oHeader,"*");
    };
    oImageTemp.m_oHeader.nDelete(Cstring(D_K_BitmapSize));

    oImageTemp.nWrite(sImageName);
    return 0;    
};

int CScanBitmap::nMaskOut(Cimage& oImage) {
    int nx,ny;
    int a2nImageDim[2];


    if (!m_pcLocalBitmap)
        return 0;
    a2nImageDim[0] = oImage.nGetDimension(0);
    a2nImageDim[1] = oImage.nGetDimension(1);

    if ((a2nImageDim[0] != m_a2nImageDim[0]) || 
        (a2nImageDim[1] != m_a2nImageDim[1]))
        return 1;

    for (nx=0;nx<m_a2nImageDim[0];nx++) {
        for (ny=0;ny<m_a2nImageDim[1];ny++) {
            if (nGetLocalBitmap(nx,ny)) 
                (oImage.*oImage.prnSetPixel)(nx,ny,0.0);
        };
    };

    return 0;
};



void CScanBitmap::vSetLocalBitmap(C3Ddata& oShoebox,int nSlice)
{
    int n0,n1;
    int* pnMask;
    int nValue;
    int nPass;
    int nx;
    int anUsedValues[8];
    
    for (nx=0;nx<8;nx++)
        anUsedValues[nx] = 0;
    

    for (nPass = 0; nPass < 2; nPass++) {
        pnMask = oShoebox.pnGetMask() + nSlice*oShoebox.m_nExt[0]*oShoebox.m_nExt[1];

        for (nx=1;nx<6;nx++) {
            if (!anUsedValues[nx])
                break;
        };
        nValue = nx;
        
        for (n1 = 0; n1 < oShoebox.m_nExt[1];n1++) {
            for (n0 = 0; n0 < oShoebox.m_nExt[0]; n0++,pnMask++) {
                if ((n0 + oShoebox.m_nOrigOffset[0] < m_a2nImageDim[0]) &&
                    (n1 + oShoebox.m_nOrigOffset[1] < m_a2nImageDim[1])) {

                    nx = nGetLocalBitmap(n0 + oShoebox.m_nOrigOffset[0],n1 + oShoebox.m_nOrigOffset[1]);
                    anUsedValues[nx]=1;
                    
                    if ((nPass==1) && (*pnMask  & C3Ddata::ms_nEllipsoidMask)) {
                        vSetLocalBitmap(n0 + oShoebox.m_nOrigOffset[0],n1 + oShoebox.m_nOrigOffset[1],(nx==0)?nValue:7);
                    };
                };
            };
        };
    };
};

