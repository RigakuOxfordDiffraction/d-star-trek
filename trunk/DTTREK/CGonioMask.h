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
// CGonioMask.h           Initial author: Thaddeus Niemeyer Jan 2003
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


#ifndef CGONIOMASK_INCLUDE
#define CGONIOMASK_INCLUDE


#include "Dtrek.h"
#include "dtarray.h"
#include "Cstring.h"

class Cdetector;
class Cgoniometer;
class Cscan;
class Cimage;
class Cimage_header;

// Defines

#include "CGonioShadow.h"

class DTREK_EXPORT CGonioMask {
    itr<int> m_a4anBoxX[4];		//[Vertex][Box#]
    itr<int> m_a4anBoxY[4];
    itr<int> m_a4anBoxZ[4];
    itr<int> m_anStage;
    int      m_nStage;
    double   m_fScale;
    
    std::vector<CGonioShadow*>      m_vpShadows;

    enum     { m_nOmegaStage,m_nChiStage,m_nPhiStage };

    Cdetector* m_poDetector;
    Cgoniometer* m_poCrysGonio;
    Cscan* m_poScan;
    unsigned char* m_pcMask;
    double   m_dDetDistance;
    double   m_dDetSwing;

    int nFindPoints(double a3x3fRotTrianglePlate[3][3],int a3x2nPixel[3][2]);

public:
    int nAddBox(itr<double>& afPosX,itr<double>& afPosY,itr<double>& afPosZ);
    void vSetOmegaStage();
    void vSetChiStage();
    void vSetPhiStage();
    void vSetScale(double fScale);
    
	int nReadGonio(Cimage_header& oHeader);
    int nReadGonio(Cstring& sImageName);
    int nReadGonio(double fOmega,double fChi,double fPhi);
	int nApplyMask(Cimage& oImage);
    int nWriteMask(Cstring& sImageOut,Cstring& sImageIn,int nValue);
    int nWriteTestMask(Cstring& sInputImage,double fValue);
    int nBuildMask();
    int nInitShadows(Cimage *poImage);
    int nMaskShadows(Cimage *poImage);
    void vErasePixelsQuad(double a4x2dPoints[4][2], Cimage *poImage);
    void vErasePixelsCircle(double a2dCenter[2], double a2dRadius[2], Cimage *poImage);
    int  nPointInPolygon(int nPoints, 
			 double *xp, double *yp, double x, double y);
    CGonioMask();
    ~CGonioMask();
};

#endif




