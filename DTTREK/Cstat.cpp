//
// Copyright 1999 Molecular Structure Corporation
//                9009 New Trails Drive
//                The Woodlands, TX, USA  77381
//
// The contents are unpublished proprietary source
// code of Molecular Structure Corporation
//
// All rights reserved
//
// Cstat.cpp    Initial author: T.L.Hendrixson           Dec 1999
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

// Description:
//
//
// ToDo:
//
//


#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <minmax.h>

#include "Cstat.h"


CstatValue Cstat::ms_fErrRef=0.0;


/****************************************************************************
 ****************************************************************************
 **                                                                        **
 **                         Public Member Functions                        **
 **                                                                        **
 ****************************************************************************
 ****************************************************************************/

/****************************************************************************
 * This is the constructor for the class.                                   *
 ****************************************************************************/
Cstat::Cstat(int nInitialSize,int _nDefAlloc) :  
    m_nNumElements(0), m_nNumAllocated(0), m_fValues(NULL),m_fWeights(NULL), m_cMarkers(NULL), m_nDefAlloc(_nDefAlloc)
                      
    {
    if (nInitialSize)
        vGrowList(nInitialSize);
    m_fLastAverage=0;
    m_fLastDeviation=0;

return;
}

 /****************************************************************************
 ****************************************************************************/

Cstat::Cstat(const Cstat& oCopy,bool bCopyMarkedOnly,int _nDefAlloc) : 
    m_nNumElements(0), m_nNumAllocated(0), m_fValues(NULL),m_fWeights(NULL), m_cMarkers(NULL), m_nDefAlloc(_nDefAlloc)
    
    {
    roCopy(oCopy,bCopyMarkedOnly);
return;
};


/****************************************************************************
 ****************************************************************************/

Cstat& Cstat::roCopy(const Cstat& oCopy,bool bCopyMarkedOnly) {
    int nx;

    vClear();
    
    if (bCopyMarkedOnly) {
        for (nx=0;nx<oCopy.m_nNumElements;nx++) 
            if (oCopy.bIsMarked(nx)) 
                vAdd(oCopy.fGet(nx));
    } else {
        if (oCopy.nSize()) {
            vGrowList(oCopy.nSize());
            // For speedier computation, just copy memory right over.
            memcpy(&fGet(0),&oCopy.fGet(0),sizeof(CstatValue)*oCopy.nSize());
            memcpy(&m_cMarkers[0],&oCopy.m_cMarkers[0],sizeof(CstatMarker)*oCopy.nSize());
        };
    };
    m_fLastAverage=0;
    m_fLastDeviation=0;

    return (*this);
};




/****************************************************************************
 * This is the destructor for the class.                                    *
 ****************************************************************************/
Cstat::~Cstat()
{
    if (m_fValues)
      {
        delete [] m_fValues;
	m_fValues = NULL;
      }
    if (m_fWeights)
      {
        delete [] m_fWeights;
	m_fWeights = NULL;
      }
    if (m_cMarkers)
      {
        delete [] m_cMarkers;
	m_cMarkers = NULL;
      }
}

/****************************************************************************
 ****************************************************************************/
void
Cstat::vAdd(CstatValue d,CstatValue dWeight,bool bMarked)
{
   vGrowList(m_nNumElements+1);
   m_fValues[m_nNumElements] = d;
   m_fWeights[m_nNumElements] = dWeight;
   m_cMarkers[m_nNumElements] = bMarked ? 1 : 0;
   ++m_nNumElements;

   return;
}

bool Cstat::bMark(bool bValue, int nIndex)
{
    if( nIndex >= 0 && nIndex < m_nNumElements ) 
    {
        bool bReturn = m_cMarkers[nIndex] == 0 ? false : true;
        
        m_cMarkers[nIndex] = bValue;

        return bReturn;
    } 
    else if( nIndex == -1 )
    {
        for(int nx=0; nx < m_nNumElements; nx++)
            m_cMarkers[nx] = bValue ?  1 : 0;
        
        return true;
    }
    
    return false;
}

/****************************************************************************
 ****************************************************************************/
void
Cstat::vClear(void)
{
    m_nNumElements=0;
   return;
}

/****************************************************************************
 ****************************************************************************/
CstatValue
Cstat::fAverage(bool bRecompute)
{ 
    int nx;
    CstatValue fAverage;
    int nValCt;
    
    if( bRecompute )
    {      
        fAverage=0.0;
        m_fSumX=0.0;
        m_fSumXW = 0.0;
        m_fSumW = 0.0;
        nValCt=0;
        for (nx=0;nx<m_nNumElements;nx++) 
        {
            if( (bool)m_cMarkers[nx] && m_cMarkers[nx] < 2 )
            {
               m_fSumX+=m_fValues[nx];
               m_fSumXW+=m_fValues[nx]*m_fWeights[nx];
               m_fSumW+=m_fWeights[nx];
               nValCt++;
            }
        }
            
        if(0 == nValCt) 
           m_fLastAverage=0.0;
        else 
           m_fLastAverage=m_fSumXW/m_fSumW;
    }
    
    return m_fLastAverage;
};

/****************************************************************************
 ****************************************************************************/
CstatValue
Cstat::fStandardDeviation(bool bRecompute) {
    int nx;
    int nValCt;

    if (bRecompute) {
        fAverage();
        if (m_nNumElements<=1) m_fLastDeviation=0.0; 
        else {
            nValCt=0;
            m_fSumX2=0.0;
            for (nx=0; nx<m_nNumElements;nx++) if ((m_cMarkers[nx]) && (m_cMarkers[nx]<2)) {
                m_fSumX2+=m_fValues[nx]*m_fValues[nx];
                nValCt++;
            };
            if (nValCt<=1) 
                m_fLastDeviation=0.0;
            else
                m_fLastDeviation=sqrt((m_fSumX2*m_nNumElements-m_fSumX*m_fSumX)/m_nNumElements/(m_nNumElements-1));
        };
    };
return m_fLastDeviation;
}

/****************************************************************************
 ****************************************************************************/

CstatValue Cstat::fCovariance(Cstat& oOther) {
    CstatValue fCov;
    int nx;

    fCov = 0.0;
    fAverage();
    oOther.fAverage();
    for (nx=0;nx<m_nNumElements;nx++) {
        //if (m_cMarkers[nx] != oOther.m_cMarkers[nx]) {
        //    return 0.0;
        //};
        if ((m_cMarkers[nx]) && (m_cMarkers[nx]<2)) {
            fCov += (m_fValues[nx] - m_fLastAverage)*(oOther.m_fValues[nx] - oOther.m_fLastAverage);
        };
    };
    return fCov;
};


/****************************************************************************
 ****************************************************************************/

CstatValue Cstat::fMin() 
{
    if (m_nNumElements==0)
        return 0.0;
    
    CstatValue      fMin = 0.0;
    bool bMinDefined= false;

    for(int nx=0; nx<m_nNumElements;nx++) 
    {
        if ((m_cMarkers[nx]) && (m_cMarkers[nx]<2))
        {
            if (!bMinDefined)
                fMin=m_fValues[nx];
            else
                fMin=min(m_fValues[nx],fMin);
        
            bMinDefined=true;        
        }
    }

    return fMin;
}


CstatValue Cstat::fMax()
{
    if (m_nNumElements==0)
        return 0.0;
    
    CstatValue fMax = 0.0;
    bool bMaxDefined= false;

    for(int nx=0; nx<m_nNumElements;nx++)
    {
        if ((m_cMarkers[nx]) && (m_cMarkers[nx]<2)) 
        {
            if (!bMaxDefined)
                fMax=m_fValues[nx];
            else
                fMax=max(m_fValues[nx],fMax);
            
            bMaxDefined=true;        
        }
    }
    return fMax;
}

/****************************************************************************
 ****************************************************************************/

int Cstat::nMarkStat(bool bMark,eMarkTest eTest,CstatValue fVal) {
    int nx,ny;
    bool bMarkTest=!bMark;
    CstatValue fTest0;
    CstatValue fTest1;
    int nMark;
    int nChangeCount;

    nChangeCount=0;
    nMark=(bMark)?(1):(0);

    switch (eTest) {
    case S_ABOVE:
    case S_BELOW:
    case S_ABOVEBELOW:
        // Reject all outside sigma bounds.
        fTest0=fStandardDeviation()*fVal;
        fTest1=fAverage(0);
        for (nx=0;nx<m_nNumElements;nx++)
            if (m_cMarkers[nx]<2)
                if ((m_cMarkers[nx]!=0)==bMarkTest) 
                    if (eTest==S_ABOVE) {
                        if (m_fValues[nx]-fTest1>fTest0) {
                            nChangeCount+=(m_cMarkers[nx]!=nMark);
                            m_cMarkers[nx]=nMark;
                        };
                    } else if (eTest==S_BELOW) {
                        if (m_fValues[nx]-fTest1<-fTest0) {
                            nChangeCount+=(m_cMarkers[nx]!=nMark);
                            m_cMarkers[nx]=nMark;
                        };
                    } else if (eTest==S_ABOVEBELOW) {
                        if (fabs(m_fValues[nx]-fTest1)>fTest0) {
                            nChangeCount+=(m_cMarkers[nx]!=nMark);
                            m_cMarkers[nx]=nMark;
                        };                           
                    };              
        break;
    case S_ABOVEBELOW_ITER: {
        int nCount = 1;

        // First, all that are currently equal to bMark must be flagged with a 2.
        for (nx=0;nx<m_nNumElements;nx++)
            if (((bMark)?(1):(0))==m_cMarkers[nx])
                m_cMarkers[nx]=2;

        // Reject all outside sigma bounds iteratively. 

        while (nCount) {
            nCount = nMarkStat(bMark,S_ABOVEBELOW,fVal);
            nChangeCount+=nCount;
        };
        // Go back and check each one that was rejected.
        fTest0=fStandardDeviation(0)*fVal;
        fTest1=fAverage(0);
        for (nx=0;nx<m_nNumElements;nx++)
            if (((bMark)?(1):(0))==m_cMarkers[nx]) {
                if (fabs(m_fValues[nx]-fTest1)<=fTest0) {
                    m_cMarkers[nx]= !bMark ? 1 : 0;
                    nChangeCount--;
                };
            };                   
        };
        // Reset all flags of 2 to the original flag value.
        for (nx=0;nx<m_nNumElements;nx++)
            if (2 == m_cMarkers[nx])
                m_cMarkers[nx]=  bMark ? 1 : 0;
        break;
    case V_ABOVE:
        for (nx=0;nx<m_nNumElements;nx++)
            if (m_cMarkers[nx]<2)
                if ((m_cMarkers[nx]!=0)==bMarkTest)
                    if (m_fValues[nx]>fVal) {
                        nChangeCount+=(m_cMarkers[nx]!=nMark);
                        m_cMarkers[nx]= nMark;
                    };
        break;
    case V_BELOW:
        for (nx=0;nx<m_nNumElements;nx++)
            if (m_cMarkers[nx]<2)
                if ((m_cMarkers[nx]!=0)==bMarkTest)
                    if (m_fValues[nx]<fVal) {
                        nChangeCount+=(m_cMarkers[nx]!=nMark);
                        m_cMarkers[nx]= nMark;
                    };
        break;
    case ALL:
        for (nx=0;nx<m_nNumElements;nx++)
            if (m_cMarkers[nx]<2) {
                nChangeCount+=(m_cMarkers[nx]!=nMark);
                m_cMarkers[nx]= nMark;
            };
        break;
    case P_ABOVE:
    case P_BELOW:
        // Reject all but the top P_ABOVE or
        // the least P_BELOW reflections.

        int*        pnPercentReject;
        int         nReject;

        nReject = ((int) fVal)*m_nNumElements+1;

        pnPercentReject = new int[nReject];
        for (nx=0;nx<nReject;nx++)
            pnPercentReject[nx]=-1;
        for (nx=0;nx<m_nNumElements;nx++) 
            if (m_cMarkers[nx]<2)
                if ((m_cMarkers[nx]!=0)==bMarkTest) {
                    
                    int nMinMax;
                    
                    nMinMax=0;
                    for (ny=0;ny<nReject;ny++) {
                        if (pnPercentReject[ny]==-1)
                            break;
                        if (eTest==P_ABOVE) {
                            // Choose the least in the reject array.
                            if (m_fValues[pnPercentReject[nMinMax]]>m_fValues[pnPercentReject[ny]])
                                nMinMax=ny;
                        } else {
                            // Choose the greatest in the reject array.
                            if (m_fValues[pnPercentReject[nMinMax]]<m_fValues[pnPercentReject[ny]])
                                nMinMax=ny;
                        };
                    };
                    if (ny==nReject) {
                        if (eTest==P_ABOVE) {
                            if (m_fValues[nx]>m_fValues[pnPercentReject[nMinMax]])
                                pnPercentReject[nMinMax]=nx;
                        } else {
                            if (m_fValues[nx]<m_fValues[pnPercentReject[nMinMax]])
                                pnPercentReject[nMinMax]=nx;
                        };
                    } else
                        pnPercentReject[ny]=nx;
                };

        // Now,flag all the non-negative entries in pnPercentReject
        for (nx=0; nx<nReject;nx++)
            if (m_cMarkers[nx]<2)
                if (pnPercentReject[nx]!=-1) {
                    nChangeCount+=(m_cMarkers[nx]!=nMark);
                    m_cMarkers[pnPercentReject[nx]]= nMark;
            };
        delete[] pnPercentReject;
        break;

    };
    return nChangeCount;
};


  
/****************************************************************************
 ****************************************************************************
 **                                                                        **
 **                        Private Member Functions                        **
 **                                                                        **
 ****************************************************************************
 ****************************************************************************/

/****************************************************************************
 ****************************************************************************/
void
Cstat::vGrowList(int nSizeDesired)
{
    if (nSizeDesired>m_nNumAllocated) {
        if(0 == m_nNumAllocated)
            m_nNumAllocated = max(m_nDefAlloc,nSizeDesired);
        else
            m_nNumAllocated = max(m_nNumAllocated*2,nSizeDesired);
        
        CstatValue* p = new CstatValue[m_nNumAllocated];
        CstatValue* pd = new CstatValue[m_nNumAllocated];
        CstatMarker* pm = new CstatMarker[m_nNumAllocated];
        
        
        memcpy(p,m_fValues,m_nNumElements*sizeof(CstatValue));
        memcpy(pd,m_fWeights,m_nNumElements*sizeof(CstatValue));
        memcpy(pm,m_cMarkers,m_nNumElements*sizeof(CstatMarker));
        if (m_fValues)
            delete [] m_fValues;
        if (m_fWeights)
            delete [] m_fWeights;
        if (m_cMarkers)
            delete [] m_cMarkers;
        m_fValues = p;
        m_fWeights = pd;
        m_cMarkers = pm;

    };

   return;
}

/****************************************************************************
 ****************************************************************************
 **                                                                        **
 **                        Private Member Functions                        **
 **                                                                        **
 ****************************************************************************
 ****************************************************************************/

/****************************************************************************
 ****************************************************************************/



CstatValue Cstat::fNormalDeviate(CstatValue fAverage,CstatValue fSigma) {
    
    static int iset=0;
    static float gset;
    float fac,rsq,v1,v2;
    
    if (iset == 0) { // We don't have an extra deviate handy, so
        do 
        {
            v1 = (float)(2.0 * fUniformDeviate() - 1.0); // pick two uniform numbers in the square ex-tending from -1 to +1 in each direction, 
            v2 = (float)(2.0 * fUniformDeviate() - 1.0);
            rsq=v1*v1+v2*v2;       // see if they are in the unit circle,
        } while (rsq >= 1.0 || rsq == 0.0);     // and if they are not, try again.
        
        fac = (float)(sqrt(-2.0f * log(rsq)/rsq));

        // Now make the Box-Muller transformation to get two normal deviates. Return one and
        //    save the other for next time.
        gset=v1*fac;
        iset=1;             // Set flag
        return fAverage + fSigma*(v2*fac);
        
    } else {                // We have an extra deviate handy,
        iset=0;             // so unset the flag,
        return fAverage + fSigma*(gset);        // and return it.
    }
    
    return 0.0;
};

CstatValue Cstat::fPoissonDeviate(CstatValue fAverage) {
    
    static float sq,alxm,g,oldm=(-1.0); // oldm is a flag for whether fAverage has changed since last call. 
    static float  Gs_fPI                 = (3.14159265f);
    float em,t,y;
    
    if (fAverage < 12.0) { // Use direct method.
        if (fAverage != oldm)
        {
            oldm = (float)fAverage;
            
            g = (float)exp(-fAverage); // If fAverage is new, compute the exponential.
        }
        
        em = -1;
        
        t=1.0;
        
        do { // Instead of adding exponential deviates it is equivalent to multiply uniform deviates. 
            // We never actually have to take the log, merely compare to the pre-computed exponential.
            ++em;
            t *= (float)fUniformDeviate();
        } while (t > g);
        
    } else { // Use rejection method.
        if (fAverage != oldm)
        { // If fAverage has changed since the last call, then pre-
            // compute some functions that occur below. 
            oldm = (float)fAverage;
            sq = (float)sqrt(2.0*fAverage);
            alxm = (float)log(fAverage);
            g = (float)(fAverage*alxm-fGamma(fAverage+1.0));
            // The function gammln is the natural log of the gamma function, as given in x 6:1.
        }
        do {
            do {                    // y is a deviate from a Lorentzian comparison function. 
                y = (float)tan(Gs_fPI*fUniformDeviate());
                em = (float)( sq * y + fAverage);   // em is y, shifted and scaled.
            } while (em < 0.0);     // Reject if in regime of zero probability.
            em = (float)floor(em);           // The trick for integer-valued distributions.

            t = (float)(0.9*(1.0+y*y)*exp(em*alxm-fGamma(em+1.0)-g));
            // The ratio of the desired distribution to the comparison function; we accept or
            //    reject by comparing it to another uniform deviate. The factor 0.9 is chosen so
            //    that t never exceeds 1.
        } while (fUniformDeviate() > t);
    }
    return em;
};

CstatValue Cstat::fUniformDeviate(int nSeed ) {


    const int IA = 16807;
    const int IM = 2147483647;
    const double AM = (1.0/IM);
    const int IQ = 127773;
    const int IR = 2836;
    const int NTAB = 32;
    const int NDIV = (1+(IM-1)/NTAB);
    const double EPS = 1.2e-7;
    const double RNMX = (1.0-EPS);

    static int idum;

    if (nSeed<0)
        idum = nSeed;

    int j;
    long k;
    static long iy=0;
    static long iv[NTAB];
    float temp;

    if (idum <= 0 || !iy) {         //Initialize.
        if (-(idum) < 1) 
            idum=1;    // Be sure to prevent idum = 0.
        else 
            idum = -(idum);
        for (j=NTAB+7;j>=0;j--) {   // Load the shuffle table (after 8 warm-ups).
            k=(idum)/IQ;
            idum=IA*(idum-k*IQ)-IR*k;
            if (idum < 0) idum += IM;
            if (j < NTAB) iv[j] = idum;
        }
        iy=iv[0];
    }

    k=(idum)/IQ;                    // Start here when not initializing.
    idum=IA*(idum-k*IQ)-IR*k;   // Compute idum=(IA*idum) % IM without over-flows by Schrage's method. 
    if (idum < 0) idum += IM;
    j=iy/NDIV;                  // Will be in the range 0..NTAB-1.
    iy=iv[j];                   // Output previously stored value and recall the shufflee table. 
    iv[j] = idum;
    if( (temp = (float)(AM*iy)) > RNMX) 
        return RNMX;   // Because users don't expect endpoint values.
    else 
        return temp;
};





CstatValue Cstat::fGamma(CstatValue fX) {
    // Internal arithmetic will be done in double precision, 
    double x,y,tmp,ser;
    static double cof[6]={76.18009172947146,-86.50532032941677,
        24.01409824083091,-1.231739572450155,
        0.1208650973866179e-2,-0.5395239384953e-5};
    int j;
    y=x=fX;
    tmp=x+5.5;
    tmp -= (x+0.5)*log(tmp);
    ser=1.000000000190015;
    for (j=0;j<=5;j++) ser += cof[j]/++y;
    return -tmp+log(2.5066282746310005*ser/x);
}

void Cstat::vAddRaw(CstatValue fd,CstatValue fw) {
    m_fSumX+=fd;
    m_fSumXW+=fd*fw;
    m_fSumW+=fw;
    m_fSumX2+=fd*fd;
    m_nNumElements++;
};

void Cstat::vClearRaw() {
    m_fSumX = 0.0;
    m_fSumXW = 0.0;
    m_fSumW = 0.0;
    m_fSumX2 = 0.0;
    m_nNumElements = 0;
};

CstatValue Cstat::fAverageRaw() {
    if (m_fSumW)
        return m_fSumXW/m_fSumW;
    else
        return 0.0;
};

CstatValue Cstat::fStandardDeviationRaw() {
    if (m_nNumElements>1)
        return sqrt((m_fSumX2*m_nNumElements-m_fSumX*m_fSumX)/m_nNumElements/(m_nNumElements-1));
    else
        return 0.0;
};

