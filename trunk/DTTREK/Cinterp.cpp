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
// Cinterp.cpp      Initial author: tjn        August 1999
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

#include "Cimage.h"
#include "dtrekvec.h"

#ifdef ANL_TIMER
#include "anlTimer.h"
#endif



/* This function determines the points that will be used in interpolation, 
 * and returns the "code number" that indicates which pixels are used.
 * This "code number" is stored in the Flag field of the list, and 
 * will inhibit nInterpolateMethod() from being called in the future.
 */


#define CU16TOI8(x)  (((int) ((x) & 15)) - 7)
#define CI8TOU16(x)  ((unsigned int) ((x) + 7)) 
#define CU32TOI16(x) (((int) ((x) & 31)) - 15)
#define CI16TOU32(x) ((unsigned int) ((x) + 15)) 
#define CU64TOI32(x) (((int) ((x) & 63)) - 31)
#define CI32TOU64(x) ((unsigned int) ((x) + 31)) 

#if (0)
/* For 4 bit interpolation (24 bit code words) use the following */
#define MAXSEARCHLOWER 7
#define MAXSEARCHUPPER 7
#define MAKEINTERPCODEWORD(px00,px01,px10,px11,px20,px21) ((CI8TOU16(px00) << 0) + (CI8TOU16(px01) << 4) + (CI8TOU16(px10) << 8) + (CI8TOU16(px11) << 12) + (CI8TOU16(px20) << 16) + (CI8TOU16(px21) << 20))
#define INTERPPOINTX(p,x)  CU16TOI8(((unsigned int) (p)) >> ((x)*4))
#define BADFLAG MAKEINTERPCODEWORD(8,8,8,8,8,8)

#else
/* For 5 bit interpolation (30 (or 32) bit code words) use the following */
#define MAXSEARCHLOWER 15
#define MAXSEARCHUPPER 31

#define MAKEINTERPCODEWORD(px00,px01,px10,px11,px20,px21) ((CI16TOU32(px00) << 0) + (CI16TOU32(px01) << 5) + (CI16TOU32(px10) << 10) + (CI16TOU32(px11) << 15) + (CI32TOU64(px20) << 20) + (CI32TOU64(px21) << 26))
#define INTERPPOINTX(p,x)  (((x)>=4)?(CU64TOI32(((unsigned int) (p)) >>  ((x)*5+((x)==5)))):(CU32TOI16(((unsigned int) (p)) >>  ((x)*5))))
#define BADFLAG MAKEINTERPCODEWORD(16,16,16,16,32,32)
#endif

int ndir_cmp(const void* vpa,const void* vpb) {
   return (*((int*) vpa)-*((int*) vpb));
};

int nInterpolateMethod(int nPix0,int nPix1,int nOutDim0,int nOutDim1,int nBadValue,Cimage* poImgBadInterp,Cimage* poImageOut) {

   enum {MAXDIRS=8};
   enum {MAXTESTDIRS=28};

   enum eDirection { N,S,E,W, NE,NW,SE,SW};
   int nDirMins[MAXDIRS*3];
   int nDirDirs[MAXDIRS][2]={{1,0},{-1,0},{0,1},{0,-1}, {1,1},{1,-1},{-1,1},{-1,-1}};
   int nDirTests[MAXTESTDIRS][4]={ 
                                 {NE,NW,S,1},{N,NW,S,2},{N,NE,S,2},{N,S,E,1},
                                 {SE,SW,N,1},{S,SW,N,2},{S,SE,N,2},{N,S,W,1},
                                 {NE,SE,W,1},{E,NE,W,2},{E,SE,W,2},{E,W,S,1},
                                 {NW,SW,E,1},{W,NW,E,2},{W,SW,E,2},{E,W,N,1},
								 {SW,NE,N,2},{SW,NE,NW,1},{SW,NE,W,2},{SW,NE,S,2},{SW,NE,SE,1},{SW,NE,E,2},
								 {NW,SE,N,2},{NW,SE,NE,1},{NW,SE,E,2},{NW,SE,S,2},{NW,SE,SW,1},{NW,SE,W,2}
								 };
   int nDir;
   int nx,ny,nz;
   int n0,n1;
   int nOffset;
   unsigned short int uiBad = (unsigned short int)nBadValue;
   unsigned short int uiPix;
   float fBad;
   float fPix;

   fBad = (float) uiBad;

   // Search for the closest non-bad pixel in each of 8 directions of compass. N,S,E,W,NE,NW,SE,SW

   for (nDir=0;nDir<MAXDIRS;nDir++) {
       for (nx=0;nx<=MAXSEARCHUPPER;nx++) {
           n0=nPix0+nDirDirs[nDir][0]*nx;
           n1=nPix1+nDirDirs[nDir][1]*nx;
           if ((n0<0) || (n0>=nOutDim0)) { 
               nx=MAXSEARCHUPPER+1; 
               break; 
           };
           if ((n1<0) || (n1>=nOutDim1)) { 
               nx=MAXSEARCHUPPER+1; 
               break; 
           };
           nOffset = n0 + (nOutDim0 * n1);
	   fPix = (poImageOut->*poImageOut->prfGetPixel)(n0, n1);
           //uiPix=poImageOut->uiGetPixel(nOffset);
           //if (uiPix>uiBad) break; // Good pixel found.
           if (fPix > fBad) break; // Good pixel found.
       };
       nDirMins[nDir*3+0]=nx*nDirDirs[nDir][0];
       nDirMins[nDir*3+1]=nx*nDirDirs[nDir][1];
       // [0]=N/S [1]=E/W
       // sqrt(2)==1.5
       if (abs(nDirDirs[nDir][0])+abs(nDirDirs[nDir][1])==2)
       {
           nDirMins[nDir*3+2]=(int)(1.5*(abs(nx*nDirDirs[nDir][0])));
       } 
       else 
       {
           nDirMins[nDir*3+2]=(abs(nx*nDirDirs[nDir][0])+abs(nx*nDirDirs[nDir][1]));
       };
   };

   // Cycle through all MAXTESTDIRS directions to get the one with the minimum cumulative distance.
   int      nLowest = 0;
   int      nLowestVal = 1000;
   int      nLowestUpperTaken = 0;
   int      nUpperTaken = 0;

   
   for (nx=0;nx<MAXTESTDIRS;nx++) {
      nUpperTaken = -1;
      for (ny=0,nz=0;ny<3;ny++) {
         nz+=nDirMins[nDirTests[nx][ny]*3+2];
         if ((abs(nDirMins[nDirTests[nx][ny]*3+0])>MAXSEARCHLOWER) || (abs(nDirMins[nDirTests[nx][ny]*3+1])>MAXSEARCHLOWER)) {
             if ((nUpperTaken!=-1) || (abs(nDirMins[nDirTests[nx][ny]*3+0])>MAXSEARCHUPPER) || (abs(nDirMins[nDirTests[nx][ny]*3+1])>MAXSEARCHUPPER)) 
                 break;
             else
                 nUpperTaken = ny;
         };
      };
	  if (ny != 3)
		continue;
      if ((nz<nLowestVal) || ((nz==nLowestVal) && (nDirTests[nx][3]<nDirTests[nLowest][3])))
      { 
	    nLowest=nx; 
	    nLowestVal=nz; 
        nLowestUpperTaken = nUpperTaken;
	  };
   };
   // Test for invalid numbers.  This will indicate there were not three good pixels within MAXSEARCH pixels of the given bad pixel.
   if (nLowestVal==1000) 
       return BADFLAG;


   // Place the lowest direction numbers in the output number.
   int nCodeWordOut[6];
   for (nx=0;nx<3;nx++) {
        nCodeWordOut[nx*2+0]=nDirMins[nDirTests[nLowest][nx]*3+0];
        nCodeWordOut[nx*2+1]=nDirMins[nDirTests[nLowest][nx]*3+1];
      };
   // The longer code word needs to be represented by the last of the three tuples.
   if ((nLowestUpperTaken!=2) && (nLowestUpperTaken!=-1)) {
       std::swap(nCodeWordOut[nLowestUpperTaken*2 + 0],nCodeWordOut[2*2 + 0]);
       std::swap(nCodeWordOut[nLowestUpperTaken*2 + 1],nCodeWordOut[2*2 + 1]);
   };
   return MAKEINTERPCODEWORD(nCodeWordOut[0],nCodeWordOut[1],nCodeWordOut[2],nCodeWordOut[3],nCodeWordOut[4],nCodeWordOut[5]);
   
};

int nDoInterpolation(int nInterpCode,int nJ0,int nJ1,int nOutDim0,Cimage* poImageOut) {
   int nMat[3][3];
   int nCofMat[3][3];
   int nPoints[3];
   int nx,ny;
   int n0,n1;
   int nDet;
   int nOffset;
   float fPix;

   /* We are solving system of the form:
      [ p1 ]   [ x1 y1 1 ]    [ a ]
      [ p2 ] = [ x2 y2 1 ] *  [ b ]
      [ p3 ]   [ x3 y3 1 ]    [ c ]
        P    =     X       *    y

      The cofactor matrix of X looks like this:

      [ (y2 - y3)  (x3 - x2)  (x2*y3 - y2*x3) ] T 
      [ (y3 - y1)  (x1 - x3)  (x3*y1 - y3*x1) ] 
      [ (y1 - y2)  (x2 - x1)  (x1*y2 - y1*x2) ]

      The inverse of X is the cofactor matrix multiplied by 1/Det
      
      where Det= (Sum of last collumn of cofactor matrix ) = (x2*y3 -y2*x3) + (x3*y1 - y3*x1) + (x1*y2 - y1*x2)


  Once a,b,c are found, the value c gives the desired interpolated value.
  This model fits  f(x,y)= a*x + b*y + c

  The values of x1,y1,x2,y2,x3,y3 are all stored in the nInterpCode element, and can be retrieved by calling INTERPPOINTX
  */

   // Build the input matrix.
   for (nx=0;nx<3;nx++) nMat[2][nx]=1;
   for (nx=0;nx<3;nx++) {
      nMat[0][nx]=INTERPPOINTX(nInterpCode,nx*2+0);
      nMat[1][nx]=INTERPPOINTX(nInterpCode,nx*2+1);
   };
   // Find the values at the given points.
   for (nx=0;nx<3;nx++) {
         n0=INTERPPOINTX(nInterpCode,nx*2+0)+nJ0;
         n1=INTERPPOINTX(nInterpCode,nx*2+1)+nJ1;
         nOffset = n0 + (nOutDim0 * n1);
	 fPix = (poImageOut->*poImageOut->prfGetPixel)(n0, n1);
         //nPoints[nx]=poImageOut->uiGetPixel(nOffset);
	 nPoints[nx] = (int) fPix;
	 //cout << "fPix: " << fPix << endl;
   };

#define x1 nMat[0][0]
#define x2 nMat[0][1]
#define x3 nMat[0][2]
#define y1 nMat[1][0]
#define y2 nMat[1][1]
#define y3 nMat[1][2]
   /* Find the integral cofactor matrix */
   /* nCofMat[0][0]=(y2-y3); nCofMat[0][1]=(x3-x2); */ nCofMat[0][2]=x2*y3 -y2*x3;
   /* nCofMat[1][0]=(y3-y1); nCofMat[1][1]=(x1-x3); */ nCofMat[1][2]=x3*y1 -y3*x1;
   /* nCofMat[2][0]=(y1-y2); nCofMat[2][1]=(x2-x1); */ nCofMat[2][2]=x1*y2 -y1*x2; 
   /* We want to find "c", which is nPoints multiplied by the last row of nCofMat */
   
   for (nx=0,ny=0;nx<3;nx++) ny+=nCofMat[nx][2]*nPoints[nx];
   for (nx=0,nDet=0;nx<3;nx++) nDet+=nCofMat[nx][2];

   return (ny / nDet);
};
  

int nPrintBadFlag(int n0,int n1,int nInterpCode) {
   int nx;
   if (nInterpCode==BADFLAG) { printf("BADPIXEL [%d %d]\n",n0,n1); return 0; };
   printf("GOOD [%d %d]:  ",n0,n1);
   for (nx=0;nx<3;nx++) {
      printf("[[%d %d] ",INTERPPOINTX(nInterpCode,nx*2+0),INTERPPOINTX(nInterpCode,nx*2+1));
   };
   printf("]\n");
   return 0;
};


int nInterpolationFunction(int nOutDim0,int nOutDim1,int nBadValue,int nNumBadInterpOut,Cimage* poImgBadInterp,Cimage* poImageOut) {
   int ni;
   int nJ0,nJ1,nJV,nJF;


   #ifdef ANL_TIMER
   anl_reset_timer(0,"Matrix Interpolation With 3x3 Inversion");
   anl_start_timer(0);
   #endif
   poImgBadInterp->nSetNextPixel(0, 0);	  
   int nResult;
   for (ni = 0; ni < nNumBadInterpOut; ni++)  {
	      nJ0 = (int) poImgBadInterp->lGetNextPixel(); // px0
	      nJ1 = (int) poImgBadInterp->lGetNextPixel(); // px1
	      nJV = (int) poImgBadInterp->lGetNextPixel(); // value
	      nJF = (int) poImgBadInterp->lGetNextPixel(); // Interp method

         if (nJF==0) {
            nJF=nInterpolateMethod(nJ0,nJ1,nOutDim0,nOutDim1,nBadValue,poImgBadInterp,poImageOut);
	         poImgBadInterp->nSetPixel(3, ni, (LONG)nJF);
         };
         if (nJF!=BADFLAG) {
            // Use this line for printing.
            //nPrintBadFlag(nJ0,nJ1,nJF);
	   nResult = (LONG) nDoInterpolation(nJF,nJ0,nJ1,nOutDim0,poImageOut);
	   //cout << "nR: " << nResult << endl;
            poImgBadInterp->nSetPixel(2,ni,(LONG) nResult);
         };
   };
   #ifdef ANL_TIMER
   anl_stop_timer(0);
   anl_print_timer(0);
   #endif
   
   return 0;
};

