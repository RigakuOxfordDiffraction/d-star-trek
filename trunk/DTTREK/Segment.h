//
// Copyright (c) 2005 Rigaku/MSC, Inc.
//
// RESTRICTED RIGHTS NOTICE SHORT FORM (JUNE 1987)
//
// Use, reproduction, or disclosure is subject to restrictions set
// forth in Contract No. W-31-109-ENG-38 and Contract No.
// 943072401 with the University of Chicago, Operator of
// Argonne National Laboratory.
//
// Segment.h       Initial author: RB               23-May-2005
//
// This files contains definitions of class CSegment

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

#ifndef DT_SEGMENT_H
#define DT_SEGMENT_H

#include "dtrekvec.h"

class DTREK_EXPORT CSegment
{
public:
    CSegment():m_dBeg(0.0), m_dEnd(0.0){}
    CSegment(double dBeg, double dEnd):m_dBeg(min(dBeg, dEnd)),m_dEnd(max(dBeg, dEnd)){} 
    CSegment(const CSegment& segment):m_dBeg(segment.m_dBeg), m_dEnd(segment.m_dEnd){}

    CSegment& operator=(const CSegment& segment){vSet(segment.m_dBeg, segment.m_dEnd); return *this;}

    // Return a segment that represents an intersection of *this* segment
    // with an input segment. If there is no intersection, return an empty
    // segment (0,0)
    CSegment operator& (const CSegment& segIn)const
    {
        CSegment    segIntersection;
        
        segIntersection.vSetBeg(max(m_dBeg, segIn.m_dBeg));
        segIntersection.vSetEnd(min(m_dEnd, segIn.m_dEnd));

        if( segIntersection.dGetBeg() >= segIntersection.dGetEnd() )
            segIntersection.vSet(0.0, 0.0);

        return segIntersection;
    }

    // Clamp an input value with this segment
    void  vClamp(double& dValue)
    {
        if( dValue < m_dBeg )
            dValue = m_dBeg;
        else if( dValue > m_dEnd )
            dValue = m_dEnd;
    }

    bool bIsIn(double dTest)const{return (m_dBeg <= dTest && m_dEnd >= dTest);}
    bool bIsIn(double dTest1, double dTest2)const{return (bIsIn(dTest1) && bIsIn(dTest2));}
    
    bool bIsInLeft(double dTest)const{return (m_dBeg <= dTest && m_dEnd > dTest);}

    void vSet(const double dBeg, const double dEnd){m_dBeg=min(dBeg, dEnd); m_dEnd=max(dBeg, dEnd);}
    void vSetBeg(const double dBeg){m_dBeg=dBeg;}
    void vSetEnd(const double dEnd){m_dEnd=dEnd;}

    double dGetBeg()const{return m_dBeg;}
    double dGetEnd()const{return m_dEnd;}

    bool bIsEmpty(){return (m_dBeg==m_dEnd);}
    double dGetWidth(){return (m_dEnd - m_dBeg);}

    bool bRoundOffTowardCenter(int nNumDecDigits)
    {
        double d1 = dRoundOff(m_dBeg, nNumDecDigits, DTREK_VEC_ROUNDOFF_CEILING);
        double d2 = dRoundOff(m_dEnd, nNumDecDigits, DTREK_VEC_ROUNDOFF_FLOOR);
        
        if( d1 <= d2 )
        {
            m_dBeg = d1;
            m_dEnd = d2;
            return true;
        }
        // Otherwise we cannot perform the requested rounding off
        return false;
    }
    
    // Does this segment overlap with the input one
//    bool bOverlap(segScanRotLimitsHeader) )

private:
    double  m_dBeg;
    double  m_dEnd;
};

#endif   // !DT_SEGMENT_H
////////////////////////////////////////////////////////////////////////

