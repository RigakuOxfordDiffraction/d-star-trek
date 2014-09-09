//
// Copyright (c) 2006 Rigaku/MSC, Inc.
//
// RESTRICTED RIGHTS NOTICE SHORT FORM (JUNE 1987)
//
// Use, reproduction, or disclosure is subject to restrictions set
// forth in Contract No. W-31-109-ENG-38 and Contract No.
// 943072401 with the University of Chicago, Operator of
// Argonne National Laboratory.
//
// CImageWait.cpp    Initial author: RB   07-Feb-2006
// This file contains the implementation of class CImageWait

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
#include <time.h>
#include "CImageWait.h"
#include "dtreksys.h"

const int cnWait_millisec = 15000;
//////////////////////////////////////////////////////////////////////////////////////
// A couple of global functions
void g_vSetImageWaitTime(const double fTime)
{
    CImageWait*     poImageWait = g_poGetImageWaitObjectPtr(DTREK_IMAGE_WAIT_GET);
    
    poImageWait->vSetWaitTime(fTime);
}

CImageWait* g_poGetImageWaitObjectPtr(DTREK_WORD wCtrl)
{
    static CImageWait*      s_poImageWait = NULL;

    if( wCtrl == DTREK_IMAGE_WAIT_GET && !s_poImageWait )
        s_poImageWait = new CImageWait();
    else if( wCtrl == DTREK_IMAGE_WAIT_DELETE && s_poImageWait )
    {
        delete s_poImageWait;
        s_poImageWait = NULL;
    }
    
    return s_poImageWait;
}
/////////////////////////////////////////////////////////////////////////////////////////
CImageWait::CImageWait()
{    
    m_nFileSize = 0L;
    m_fWaitTime = 0.0;

    // Determine required percent file size 
    m_nPercentFilesize = 90;
    if( "" != sGetEnv("DTREK_DTINT_PERCENTFILESIZE") )
    {
        m_nPercentFilesize = atoi(sGetEnv("DTREK_DTINT_PERCENTFILESIZE").string());
        printf("INFO: Image file size required to be at least %d%% of the first image size.\n", m_nPercentFilesize);
        fflush(stdout);
    }
}
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
CImageWait::CImageWait(const CImageWait& oWait)
{
    m_nFileSize = oWait.m_nFileSize;
    m_fWaitTime = oWait.m_fWaitTime;
    m_nPercentFilesize = oWait.m_nPercentFilesize;
}
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
CImageWait& CImageWait::operator=(const CImageWait& oWait)
{
    m_nFileSize = oWait.m_nFileSize;
    m_fWaitTime = oWait.m_fWaitTime;
    m_nPercentFilesize = oWait.m_nPercentFilesize;
    
    return *this;
}

CImageWait::~CImageWait()
{
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Wait until the first image appears and save its size in a class member m_nFileSize
bool CImageWait::bDoWaitForFirstImage(Cstring& sImagePath)
{
    if( 0.0 >= m_fWaitTime )
        return true;  // no wait required

    time_t       tStartTime = 0L;
    time(&tStartTime);
    
    time_t       tCurrentTime = 0L;
    double       fTimeElapsed = 0.0;

    // Wait until the image file appears
    while( m_fWaitTime * 2.0 > fTimeElapsed && lFileGetSize(sImagePath) <= 0 )
    {
        time(&tCurrentTime);
        fTimeElapsed = (double) tCurrentTime - (double) tStartTime;

        printf("... waiting for FIRST image, elapsed time: %f / %f\n", fTimeElapsed, m_fWaitTime);
        fflush(stdout);

        nWait(cnWait_millisec);
    }

    if( lFileGetSize(sImagePath) <= 0 )
    {
        printf("First image never arrived!\n");
        return false; // the image never showed up
    }
    // Wait until the image file size stops growing
    do
    {
        m_nFileSize = lFileGetSize(sImagePath);

        nWait(cnWait_millisec);
    } while ( lFileGetSize(sImagePath) > m_nFileSize);

    return true;
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//  Wait until the image (1) appears and (2) reaches a certain percentage of the previously-read image size saved. 
//  This helps prevent reading a partial image file that is still being written.                             
/// NOTE: This is not useful for compressed images
bool CImageWait::bDoWaitForNextImage(Cstring& sImagePath)
{
    if( 0.0 >= m_fWaitTime )
        return true; // no wait required

    time_t       tStartTime = 0L;
    (void) time(&tStartTime);
    
    time_t       tCurrentTime = 0L;
    double       fTimeElapsed = 0.0;

    LONG lCurrFileSize = 0;
    LONG lPrevFileSize = 0;
    lCurrFileSize = lFileGetSize(sImagePath);
    lPrevFileSize = lCurrFileSize;
    while(    ((lCurrFileSize / m_nPercentFilesize) < m_nFileSize/100)
           && (m_fWaitTime > fTimeElapsed) )
      {
        time(&tCurrentTime);
        fTimeElapsed = (double)tCurrentTime - (double) tStartTime;
        
        printf("... waiting for image, elapsed time: %f / %f\n", fTimeElapsed, m_fWaitTime);
        fflush(stdout);
	lPrevFileSize = lCurrFileSize;
        nWait(cnWait_millisec);
	lCurrFileSize = lFileGetSize(sImagePath);
//+JWP 2009-Apr-01
	if ( (0 < lPrevFileSize) && (lCurrFileSize == lPrevFileSize) )
	  break; // Leave this loop if the file size did not change from the previous positive value
//-JWP 2009-Apr-01
      }
    
    if( m_fWaitTime <= fTimeElapsed ) 
        return false; // image size hasn't reached the target within the allotted time

    return true;
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
bool CImageWait::bSaveImageSize(Cstring& sImagePath)
{
    m_nFileSize = lFileGetSize(sImagePath);

    if( m_nFileSize <= 0 )
    {
        printf("\nERROR: Cannot determine file size for %s\n",sImagePath.string());
        return false; // shouldn't happen
    }
    
    printf("\nINFO: Initial image file size is %d bytes.", m_nFileSize);
    printf("\nSubsequent image file must be at least %d%% of this size.", m_nPercentFilesize); 
    printf("\nbefore they will be read.\n");

    return true;
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
