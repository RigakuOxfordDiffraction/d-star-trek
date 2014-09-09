//
// Copyright (c) 2004 Rigaku/MSC, Inc.
//
// RESTRICTED RIGHTS NOTICE SHORT FORM (JUNE 1987)
//
// Use, reproduction, or disclosure is subject to restrictions set
// forth in Contract No. W-31-109-ENG-38 and Contract No.
// 943072401 with the University of Chicago, Operator of
// Argonne National Laboratory.
//
// CIceRingDetector.h     Initial author: RB        05-Feb-2003

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

#include "CIceRingLocator.h"
#include "CBeamData.h"
#include "Cscan.h"

CIceRingLocator::CIceRingLocator(Cimage_header* poHeader) 
{
	m_poHeader = poHeader;
}
/////////////////////////////////////////////////////////////////////////////////
void CIceRingLocator::vLocateRings(int nImageSeqNumber,
                                   double dMinTwoTheta_DEG,
                                   double dExpectedRingWidth_DEG,
                                   int nImageSizeCompressionFactor)
{
    if(NULL == m_poHeader)
    {
        printf("ERROR: No header for resolution ring detection.\n");
        return;
    }
    
    //////////////////////////////////////////////////////////////
    // Get image name according to header and sequence number
    Cscan           oScan(*m_poHeader);
    Cstring         strImageName = oScan.sGetImageName(nImageSeqNumber);
    
    printf("\nProcessing image '%s' for powder ring detection.\n", strImageName.string());
    printf("Search starts at 2Theta = %.1f deg. Expected ring width = %.1f deg\n\n", dMinTwoTheta_DEG, dExpectedRingWidth_DEG);

    Cimage		    oImage(strImageName);
    //+JWP_DBG
    if (500 >= oImage.nGetDimension(0))
	// If the image is already small in number of pixels, then do not compress
	nImageSizeCompressionFactor = 2;
    //-JWP_DBG
    Cbeamdata       oBeamData(oImage, nImageSizeCompressionFactor);
    
    static itr<double>     s_adKnownIceRingResolution_A;
    s_adKnownIceRingResolution_A  + 3.917 
                                + 3.684
                                + 3.458
                                + 2.683
                                + 2.261
                                + 2.081
                                + 1.9584
                                + 1.9272
                                + 1.8927
                                + 1.7292;

    if( 0 != oBeamData.nFindResoRings(dMinTwoTheta_DEG,
                                      dExpectedRingWidth_DEG,
                                      nImageSizeCompressionFactor,
                                      s_adKnownIceRingResolution_A,
                                      m_poHeader) )
    {
        printf("ERROR: Problem during resolution ring detection.\n");
    }
}
