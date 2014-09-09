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
// dtsvd.h      Initial author: Thom Hendrixson        September 1999
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



/****************************************************************************
 * This file contains routines for calculating linear and non-linear        *
 * least-squares fits to data.  Most of the routines in this file were      *
 * taken from "Numerical Recipies in C", Press, Flannery, Teukolsky and     *
 * Vetterling, (Cambridge, 1991).  There are two versions of the linear     *
 * least-squares, svdfit() and least_squares().  svdfit is the routine      *
 * from Num. Rec. in C and least_squares() is a variant in which the matrix *
 * u[][] is precomputed and passed as x[][] (see below).                    *
 *                                                                          *
 * Caveat emptor.                                                           *
 *                                                                          *
 * T. Hendrixson, August 1993                                               *
 *                                                                          *
 * MODIFICATION HISTORY:                                                    *
 *   (NOTE: If you modify the programme, make a note both here and at the   *
 *    location of the modification.  the note here should be general and    *
 *    the note at the modification should be more specific.)                *
 *                                                                          *
 ****************************************************************************/

/****************************************************************************
 *                               Header Files                               *
 ****************************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <sys/types.h>

#ifdef CONIO
#include <conio.h>
#endif

#include "dtsvd.h"
#include <math.h>
#include "dtrekvec.h"
#include "string.h"
#include "dtreksys.h"
#include "dtarray.h"

#ifdef SUNOS
// This gets the SunOS swap routine
#include <algorithm>
#endif

#if !defined(VC6)
using std::cout;
using std::endl;
using std::flush;
#endif

/****************************************************************************
 *                            Macro Definitions                             *
 ****************************************************************************/
#define TRUE 1
#define FALSE 0
#define SINGULAR_LIMIT .00001

/****************************************************************************
 * This routine performs a linear equation solution by Gauss-Jordan elimi-  *
 * nation.  a[1..n][1..n] is an input matrix of n x n elements.             *
 * b[1..n][1..m] is an input matrix of n by m containing the m right-hand   *
 * side vectors.  On output, a is replaces by its matrix inverse and b is   *
 * replaced by the corresponding set of solution vectors.                   *
 ****************************************************************************/
int gaussj(double **a,
            int n,
            double **b,
            int m)
{
   int *indxc,*indxr,*ipiv;
   int i,icol,irow,j,k,l,ll;
   double big,dum,tmp,pivinv;

/*
 * The integer arrays indxr[1..n], indxc[1..n] and ipiv[1..n] are used for
 * bookkeeping on the pivoting.
 */
   indxr = ivector(1,n);
   indxc = ivector(1,n);
   ipiv = ivector(1,n);
   for(j = 1; j <= n; j++)
      ipiv[j] = 0;
/*
 * This is the main loop over the columns to be reduced.
 */
   for(i = 1; i <= n; i++){
      big = 0.;
/*
 * This is the outer loop of the search for a pivot element.
 */
      for(j = 1; j <= n; j++){
         if(ipiv[j] != 1){
            for(k = 1; k <= n; k++){
               if(ipiv[k] == 0){
                  if(fabs(a[j][k]) >= big){
                     big = fabs(a[j][k]);
                     irow = j;
                     icol = k;
                  }
               }
               else if(ipiv[k] > 1){
                  fprintf(stderr,"Singular matrix-1 in gaussj()\n");
#ifdef CONIO
                  getch();
#endif
                  return 1;
               }
            }
         }
      }
      (ipiv[icol])++;
/*
 * We now have the pivot element, so we interchange rows, if needed, to put
 * the pivot element on the diagonal.  The columns are not physically
 * interchanged, only relabelled: indec[i], the column of the ith pivot
 * element, is the ith column that is reduced, while indxr[i] is the row
 * in which that pivot element was originally located.  If indxr[i].ne.
 * indxc[i] there is an implied column interchange.  With this form of
 * bookkeeping, the solution b's will end up in the correct order, and the
 * inverse matrix will be scrambled by columns.
 */
      if(irow != icol){
         for(l = 1; l <= n; l++){
            tmp = a[irow][l]; a[irow][l] = a[icol][l]; a[icol][l] = tmp;
         }
         for(l = 1; l <= m; l++){
            tmp = b[irow][l]; b[irow][l] = b[icol][l]; b[icol][l] = tmp;
         }
      }
/*
 * We are ready to divide the pivot row by the pivot element located at 
 * irow and icol.
 */
      indxr[i] = irow;
      indxc[i] = icol;
      if(a[icol][icol] == 0.){
         fprintf(stderr,"Singular matrix-2 in gaussj(): icol %d\n",icol);
#ifdef CONIO
                  getch();
#endif
         return 1;
      }
      pivinv = 1./a[icol][icol];
      a[icol][icol] = 1.;
      for(l = 1; l < n; l++)
         a[icol][l] *= pivinv;
      for(l = 1; l <= m; l++)

		  b[icol][l] *= pivinv;
/*
 * Next we reduce the rows except for the pivot one.
 */
      for(ll = 1; ll <= n; ll++){
         if(ll != icol){
            dum = a[ll][icol];
            a[ll][icol] = 0.;
            for(l = 1; l <= n; l++)
               a[ll][l] -= a[icol][l]*dum;
            for(l = 1; l <= m; l++)
               b[ll][l] -= b[icol][l]*dum;
         }
      }
   }

/*
 * This is the end of the main loop over the columns of the reduction.  It 
 * only remains to unscramble the solution in view of the column inter-
 * changes.  We do this by interchanging pairs of columns in the reverse
 * order that the permutation was built up.
 */
   for(l = n; l >= 1; l--){
      if(indxr[l] != indxc[l]){
         for(k = 1; k <= n; k++){
            tmp = a[k][indxr[l]];
            a[k][indxr[l]] = a[k][indxc[l]];
            a[k][indxc[l]] = tmp;
         }
      }
   }

/*
 * All done.
 */
   free_ivector(ipiv,1,n);
   free_ivector(indxr,1,n);
   free_ivector(indxc,1,n);
   return 0;
}


/****************************************************************************
 * This routine allocates a double vector with the range [low_index,         *
 * high_index].                                                             *
 ****************************************************************************/
double *vectorDT(int nl,
              int nh)
{
   double *v = NULL;

   //   v = (double *)malloc((unsigned)(nh-nl+1)*sizeof(double));
   v = new double [nh-nl+1];
   if(!v){
      fprintf(stderr,"allocation failure in vectorDT()\n");
#ifdef CONIO
      getch();
#endif
      exit(1);
   }
   else
     memset((void*)v, 0, (unsigned)(nh-nl+1)*sizeof(double));
   return v-nl;
}

/****************************************************************************
 * This routine allocates an integer vector with the range [low_index,      *
 * high_index].                                                             *
 ****************************************************************************/
int *ivector(int nl,
             int nh)
{
   int *v = NULL;

   //   v = (int *)malloc((unsigned)(nh-nl+1)*sizeof(int));
   v = new int [nh-nl+1];
   if(!v){
      fprintf(stderr,"allocation failure in ivector()\n");
#ifdef CONIO
      getch();
#endif
      exit(1);
   }
   else
     memset((void*)v, 0, (unsigned)(nh-nl+1)*sizeof(int));
   return v-nl;
}

/****************************************************************************
 * This routine allocates a double matrix with the range [nrl,nrh][ncl,nch]. *
 ****************************************************************************/
double **matrix(int nrl,
               int nrh,
               int ncl,
               int nch)
{
   int i;
   double **m = NULL;

   //   m = (double **)malloc((unsigned)(nrh-nrl+1)*sizeof(double *));
   m = new double* [nrh-nrl+1];
   if (!m)
     {
       fprintf(stderr,"allocation failure 1 in matrix()\n");
#ifdef CONIO
       getch();
#endif
       exit(1);
     }
   else
     memset((void*)m, 0, (unsigned)(nrh-nrl+1)*sizeof(double*));
       
   m -= nrl;
   for(i = nrl; i <= nrh; i++)
     {
       // m[i] = (double *)malloc((unsigned)(nch-ncl+1)*sizeof(double));
       m[i] = new double [nch-ncl+1];
       if(!m[i])
	 {
	   fprintf(stderr,"allocation failure 2 in matrix()\n");
#ifdef CONIO
	   getch();
#endif
	   exit(1);
	 }
       else
	 memset((void*)m[i], 0, (unsigned)(nch-ncl+1)*sizeof(double));
      m[i] -= ncl;
     }
   return m;
}

/****************************************************************************
 * This routine frees a double vector allocated using vectorDT().              *
 ****************************************************************************/
void free_vector(double *v,
                 int nl,
                 int nh)
{
  //   free(v+nl);

  double *pdTemp;
  pdTemp = v + nl;
  delete [] pdTemp;
  pdTemp = NULL;
  return;
}

/****************************************************************************
 * This routine frees an integer vector allocated using ivector().          *
 ****************************************************************************/
void free_ivector(int *v,
                  int nl,
                  int nh)
{
   //   free(v+nl);

  int *pnTemp;
  pnTemp = v + nl;
  delete [] pnTemp;
  pnTemp = NULL;
  return;
}

/****************************************************************************
 * This routine frees a double matrix allocated using matrix().              *
 ****************************************************************************/
void free_matrix(double **m,
                 int nrl,
                 int nrh,
                 int ncl,
                 int nch)
{
   int i;
   double *pdTemp;
   double **ppdTemp;
   for(i = nrh; i >= nrl; i--)
     {
       pdTemp = (m[i]+ncl);
       delete [] pdTemp;
       pdTemp = NULL;
     }
     //      free(m[i]+ncl);
     //   free(m+nrl);
   ppdTemp = (m+nrl);
   delete [] ppdTemp;
   ppdTemp = NULL;
   return;
}

/****************************************************************************
 * This routine copies the sign of b to a, returning the new value of a.    *
 ****************************************************************************/

double copysignf(double a,
                double b)
{
   a = fabs(a);
//   a = (b < 0.) ? (-a) : (a);
   if (b < 0.) 
     a = -a;
   return (a);
}
      
/****************************************************************************
 * This routine computes the (a*a+b*b) without destructive overflow or  *
 * underflow.                                                               *
 ****************************************************************************/
double pythag(double a,
             double b)
{
   double c;

   a = fabs(a);
   b = fabs(b);
   if(a > b){
      c = b/a;
      c = sqrt(1.+c*c)*a;
   }
   else{
      if(b == 0.)
         c = 0.;
      else{
         c = a/b;
         c = sqrt(1+c*c)*b;
      }
   }

   return c;
}

/****************************************************************************
 * This routine solves AX=B for a vector X, where A is specified by the     *
 * arrays u[1..m][1..n], w[1..n] and v[1..n][1..n] as returned by svdcmp(). *
 * m and n are the dimensions of A and will be equal for square matrices.   *
 * b[1..m] is the input right-hand side.  x[1..n] is the output solution    *
 * vector.  No input quantities are destroyed, so the routine may be called *
 * sequentially with different b[]s.                                        *
 ****************************************************************************/
void svbksb(double **u,
            double  *w,
            double **v,
            int     m,
            int     n,
            double  *b,
            double  *x)
{
   int i,j,jj;
   double s,*tmp;

   tmp = vectorDT(1,n);

/*
 * Calculate U(transpose)*B
 */
   for(j = 1; j <= n; j++){
      s = 0.;
      if(w[j]){
         for(i = 1; i <= m; i++)
             s += u[i][j]*b[i];
         s /= w[j];
      }
      tmp[j] = s;
   }

/*
 * Multiply result by V to get answer
 */
   for(j = 1; j <= n; j++){
      s = 0.;
      for(jj = 1; jj <= n; jj++)
         s += v[j][jj]*tmp[jj];
      x[j] = s;
   }

   delete[] (tmp+1);
   return;
}

int svdcmp(double **a, int m, int n, double* w, double **v)
{
    double pythag(double a, double b);
    int flag,i,its,j,jj,k,l,nm;
    double anorm,c,f,g,h,s,scale,x,y,z,*rv1;
    rv1=vectorDT(1,n);
    g=scale=anorm=0.0; 
    for (i=1;i<=n;i++) {
        l=i+1;
        rv1[i]=scale*g;
        g=s=scale=0.0;
        if (i <= m) {
            for (k=i;k<=m;k++) scale += fabs(a[k][i]);
            if (scale) {
                for (k=i;k<=m;k++) {
                    a[k][i] /= scale;
                    s += a[k][i]*a[k][i];
                }
                f=a[i][i];
                g = -SIGND(sqrt(s),f);
                h=f*g-s;
                a[i][i]=f-g;
                for (j=l;j<=n;j++) {
                    for (s=0.0,k=i;k<=m;k++) s += a[k][i]*a[k][j];
                    f=s/h;
                    for (k=i;k<=m;k++) a[k][j] += f*a[k][i];
                }
                for (k=i;k<=m;k++) a[k][i] *= scale;
            }
        }
        w[i]=scale *g;
        g=s=scale=0.0;
        if(i<= m && i != n) {
            for (k=l;k<=n;k++) scale += fabs(a[i][k]);
            if (scale) {
                
                for (k=l;k<=n;k++) {
                    a[i][k] /= scale;
                    s += a[i][k]*a[i][k];
                }
                f=a[i][l];
                g = -SIGND(sqrt(s),f);
                h=f*g-s;
                a[i][l]=f-g;
                for (k=l;k<=n;k++) rv1[k]=a[i][k]/h;
                for (j=l;j<=m;j++) {
                    for (s=0.0,k=l;k<=n;k++) s += a[j][k]*a[i][k];
                    for (k=l;k<=n;k++) a[j][k] += s*rv1[k];
                }
                for (k=l;k<=n;k++) a[i][k] *= scale;
            }
        }
        anorm=max(anorm,(fabs(w[i])+fabs(rv1[i])));
    }
    for (i=n;i>=1;i--) { 
        if(i < n) {
            if (g) {
                for (j=l;j<=n;j++) 
                    v[j][i]=(a[i][j]/a[i][l])/g;
                for (j=l;j<=n;j++) {
                    for (s=0.0,k=l;k<=n;k++) s += a[i][k]*v[k][j];
                    for (k=l;k<=n;k++) v[k][j] += s*v[k][i];
                }
            }
            for (j=l;j<=n;j++) v[i][j]=v[j][i]=0.0;
        }
        v[i][i]=1.0;
        g=rv1[i];
        l=i;
    }
    for (i=min(m,n);i>=1;i--) { 
        l=i+1;
        g=w[i];
        for (j=l;j<=n;j++) a[i][j]=0.0;
        if (g) {
            g=1.0/g;
            for (j=l;j<=n;j++) {
                for (s=0.0,k=l;k<=m;k++) s += a[k][i]*a[k][j];
                f=(s/a[i][i])*g;
                for (k=i;k<=m;k++) a[k][j] += f*a[k][i];
            }
            for (j=i;j<=m;j++) a[j][i] *= g;
        } else for (j=i;j<=m;j++) a[j][i]=0.0;
        ++a[i][i];
    }
    for (k=n;k>=1;k--) { 
        for (its=1;its<=200;its++) {
            flag=1;
            for (l=k;l>=1;l--) { 
                nm=l-1; 
                if ((float)(fabs(rv1[l])+anorm) == (float) anorm) {
                    flag=0;
                    break;
                }
                if ((float)(fabs(w[nm])+anorm) == (float) anorm) 
                    break;
            }
            if (flag) {
                c=0.0; 
                s=1.0;
                for (i=l;i<=k;i++) {
                    f=s*rv1[i];
                    rv1[i]=c*rv1[i];
                    if ((float)(fabs(f)+anorm) == (float) anorm) 
                        break;
                    g=w[i];
                    h=pythag(f,g);
                    w[i]=h;
                    h=1.0/h;
                    c=g*h;
                    s = -f*h;
                    for (j=1;j<=m;j++) {
                        y=a[j][nm];
                        z=a[j][i];
                        a[j][nm]=y*c+z*s;
                        a[j][i]=z*c-y*s;
                    }
                }
            }
            z=w[k];
            if (l == k) { 
                if (z < 0.0) { 
                    w[k] = -z;
                    for (j=1;j<=n;j++) v[j][k] = -v[j][k];
                }
                break;
            }
            if (its == 200) {
                // nrerror("no convergence in 30 svdcmp iterations");
                return 1;
            };
            x=w[l]; 
            nm=k-1;
            y=w[nm];
            g=rv1[nm];
            h=rv1[k];
            f=((y-z)*(y+z)+(g-h)*(g+h))/(2.0*h*y);
            g=pythag(f,1.0);
            f=((x-z)*(x+z)+h*((y/(f+SIGND(g,f)))-h))/x;
            c=s=1.0; 
            for (j=l;j<=nm;j++) {
                i=j+1;
                g=rv1[i];
                y=w[i];
                h=s*g;
                g=c*g;
                z=pythag(f,h);
                rv1[j]=z;
                c=f/z;
                s=h/z;
                f=x*c+g*s;
                g = g*c-x*s;
                h=y*s;
                y *= c;
                for (jj=1;jj<=n;jj++) {
                    x=v[jj][j];
                    z=v[jj][i];
                    v[jj][j]=x*c+z*s;
                    v[jj][i]=z*c-x*s;
                }
                z=pythag(f,h);
                w[j]=z; 
                if (z) {
                    z=1.0/z;
                    c=f*z;
                    s=h*z;
                }
                f=c*g+s*y;
                x=c*y-s*g;
                for (jj=1;jj<=m;jj++) {
                    y=a[jj][j];
                    z=a[jj][i];
                    a[jj][j]=y*c+z*s;
                    a[jj][i]=z*c-y*s;
                }
            }
            rv1[l]=0.0;
            rv1[k]=f;
            w[k]=x;
        }
    }
    free_vector(rv1,1,n);
    return 0;
}

                



int nSolveAxB_svd(int nN, float* fANxN, float* fXN, float* fBN)
{
  // Return 0   on SUCCESS
  //      non-0 on FAILURE

  static double** fMatA=NULL;
  static double** fMatV=NULL;
  static double*  fVecW=NULL;
  static double*  fVecX=NULL;
  static double*  fVecB=NULL;
  static int nMatSize=0;
  int nx,ny;

  if (nN>nMatSize) 
    {
      // Allocate variables if they are too small.

      //if (bCheck_matrix(fMatA,1,nMatSize,1,nMatSize))
      if (fMatA)
	free_matrix(fMatA,1,nMatSize,1,nMatSize);
      //if (bCheck_matrix(fMatV,1,nMatSize,1,nMatSize))
      if (fMatV)
	free_matrix(fMatV,1,nMatSize,1,nMatSize);
      if (fVecW) free_vector(fVecW,1,nMatSize);
      if (fVecX) free_vector(fVecX,1,nMatSize);
      if (fVecB) free_vector(fVecB,1,nMatSize);
      fMatA = matrix(1,nN,1,nN);
      fMatV = matrix(1,nN,1,nN);
      fVecW = vectorDT(1,nN);
      fVecB = vectorDT(1,nN);
      fVecX = vectorDT(1,nN);
      nMatSize = nN;
    }

  float *pfTemp;
  double fMult;

  fMult = fScaleMatMN(nN,nN,fANxN,FALSE);

  pfTemp = fANxN;
  for (nx=1;nx<=nN;nx++)
    for (ny=1;ny<=nN;ny++,pfTemp++) 
      {
	fMatA[ny][nx]=*pfTemp*fMult;
      }

  if (svdcmp(fMatA,nN,nN,fVecW,fMatV))
    return (1);
  for (nx=0;nx<nN;nx++) 
    {
      fVecB[nx+1] = fBN[nx];
    }
  svbksb(fMatA,fVecW,fMatV,nN,nN,fVecB,fVecX);

  // Undo the affect of the multiplication.
  for (nx=0; nx < nN; nx++) 
    fXN[nx] = fVecX[nx+1] * fMult;
  return (0);
}

int nInvMatND_svd(int nN,float* fANxN,float* fANxNInv)
{
  // Return   0 on SUCCESS
  //      non-0 on FAILURE
    static double** fMatA=NULL;
	static double** fMatV=NULL;
	static double*  fVecW=NULL;
	static int nMatSize=0;
	int nx,ny,nz;
    double f0;

    if (0 >= nN)
      {
	cout << "WARNING in nInvMATND_svd, nN is 0!\n" << flush;
	return (1);
    } else if (nN == 1) {
        if (*fANxN==0.0)
            return (1);
        else
            *fANxNInv = 1.0/(*fANxN);

    } else if (nN == 2) {
        if (0.0==fInvMat2D(fANxN,fANxNInv))
            return (1);
        else
            return (0);

    } else if (nN == 3) {
        if (0.0==fInvMat3D(fANxN,fANxNInv))
            return (1);
        else
            return (0);
    } /* else if (nN == 4) {
        if (0.0==fInvMat4D(fANxN,fANxNInv))
            return (1);
        else
            return (0);
    }; */

	if (nN>nMatSize) {
	  //if (bCheck_matrix(fMatA,1,nMatSize,1,nMatSize))
	  if (fMatA)
	    free_matrix(fMatA,1,nMatSize,1,nMatSize);
	  //if (bCheck_matrix(fMatV,1,nMatSize,1,nMatSize))
	  if (fMatV)
	    free_matrix(fMatV,1,nMatSize,1,nMatSize);
	  if (fVecW) free_vector(fVecW,1,nMatSize);
	  fMatA=matrix(1,nN,1,nN);
	  fMatV=matrix(1,nN,1,nN);
	  fVecW=vectorDT(1,nN);
	  nMatSize=nN;
	};
    float *pfTemp;
    double fMult;

    fMult = fScaleMatMN(nN,nN,fANxN,FALSE);

    pfTemp = fANxN;
	for (nx=1;nx<=nN;nx++)
		for (ny=1;ny<=nN;ny++,pfTemp++) {
			fMatA[ny][nx]=*pfTemp*fMult;
		};

	if (svdcmp(fMatA,nN,nN,fVecW,fMatV))
        return (1);
    
    // Want to compute V [W]-1 [U]t
    
    // Multiply V by the diagonal matrix W
    for (nx=1;nx<=nN; nx++)
        for (ny=1; ny<=nN; ny++)
	  {
	    if (0.0f == fVecW[nx])
	      {
		cout << "WARNING in nInvMatND_svd, divide by 0 fVecW[nx]!\n" << flush;
		return (1);
	      }
            fMatV[ny][nx]/=fVecW[nx];
	  }
        
    // Multiply this by [U]t
    for (nx=1;nx<=nN; nx++) 
        for (ny=1; ny<=nN; ny++) {
            f0=0;
            for (nz=1; nz<=nN; nz++)
                f0+=fMatV[nx][nz]*fMatA[ny][nz];    // Note: A:=U after call to svdcmp
            fANxNInv[(nx-1)+(ny-1)*nN]=f0;
        };
    // Undo the affect of the multiplication.
    for (nx=0;nx<nN*nN;nx++)
        fANxNInv[nx]*=fMult;
	return (0);
};




int nSolveAxB_svd(int nN,double* fANxN,double* fXN,double* fBN)
{
  // Return 0   on SUCCESS
  //      non-0 on FAILURE

  static double** fMatA=NULL;
  static double** fMatV=NULL;
  static double*  fVecW=NULL;
  static double*  fVecX=NULL;
  static double*  fVecB=NULL;
  static int nMatSize=0;
  int nx,ny;

  if (nN>nMatSize) {
    //if (bCheck_matrix(fMatA,1,nMatSize,1,nMatSize))
    if (fMatA)
      free_matrix(fMatA,1,nMatSize,1,nMatSize);
    //if (bCheck_matrix(fMatV,1,nMatSize,1,nMatSize))
    if (fMatV)
      free_matrix(fMatV,1,nMatSize,1,nMatSize);
    if (fVecW) free_vector(fVecW,1,nMatSize);
    if (fVecX) free_vector(fVecX,1,nMatSize);
    if (fVecB) free_vector(fVecB,1,nMatSize);
    fMatA=matrix(1,nN,1,nN);
    fMatV=matrix(1,nN,1,nN);
    fVecW=vectorDT(1,nN);
    fVecB=vectorDT(1,nN);
    fVecX=vectorDT(1,nN);
    nMatSize=nN;
  };
  double *pfTemp;
  double fMult;

  fMult = fScaleMatMN(nN,nN,fANxN,FALSE);

  pfTemp = fANxN;
  for (nx=1;nx<=nN;nx++)
    for (ny=1;ny<=nN;ny++,pfTemp++) {
      fMatA[ny][nx]=*pfTemp*fMult;
    };

  if (svdcmp(fMatA,nN,nN,fVecW,fMatV))
    return (1);
  for (nx=0;nx<nN;nx++) {
    fVecB[nx+1] = fBN[nx]*fMult;
  };
  svbksb(fMatA,fVecW,fMatV,nN,nN,fVecB,fVecX);

  // Undo the affect of the multiplication.
  for (nx=0;nx<nN;nx++) 
    fXN[nx]=fVecX[nx+1];
    
  return (0);
};

int nInvMatND_svd(int nN,double* fANxN,double* fANxNInv,double* fANxNEigenVecs,double* fANEigenValues)
{
  // Return   0 on SUCCESS
  //      non-0 on FAILURE
    static double** fMatA=NULL;
    static double** fMatV=NULL;
    static double*  fVecW=NULL;
    static int nMatSize=0;
    int nx,ny,nz;
    double f0;

    if (0 >= nN)
    {
        cout << "WARNING in nInvMATND_svd, nN is 0!\n" << flush;
        return (1);
        
    } else if (fANxNEigenVecs || fANEigenValues) {
        // Don't simplify:  We want eigen vectors as well.
    } else if (nN == 1) {
        if (*fANxN==0.0)
            return (1);
        else
            *fANxNInv = 1.0/(*fANxN);

    } else if (nN == 2) {
        if (0.0==fInvMat2D(fANxN,fANxNInv))
            return (1);
        else
            return (0);

    } else if (nN == 3) {
        if (0.0==fInvMat3D(fANxN,fANxNInv))
            return (1);
        else
            return (0);
    };

    if (nN>nMatSize) {
      //if (bCheck_matrix(fMatA,1,nMatSize,1,nMatSize))
      if (fMatA)
	free_matrix(fMatA,1,nMatSize,1,nMatSize);
      //if (bCheck_matrix(fMatV,1,nMatSize,1,nMatSize))
      if (fMatV)
	free_matrix(fMatV,1,nMatSize,1,nMatSize);
      if (fVecW) free_vector(fVecW,1,nMatSize);
      fMatA=matrix(1,nN,1,nN);
      fMatV=matrix(1,nN,1,nN);
      fVecW=vectorDT(1,nN);
      nMatSize=nN;
    };
    double *pfTemp;
    double fMult;

    fMult = fScaleMatMN(nN,nN,fANxN,FALSE);

    pfTemp = fANxN;
    for (nx=1;nx<=nN;nx++)
      for (ny=1;ny<=nN;ny++,pfTemp++) {
	fMatA[ny][nx]=*pfTemp*fMult;
      };
    
    if (svdcmp(fMatA,nN,nN,fVecW,fMatV))
      return (1);
    
    // Want to compute V [W]-1 [U]t
    
    // Multiply V by the diagonal matrix W
    for (nx=1;nx<=nN; nx++) {
        if (0.0f == fVecW[nx])
        {
            
            cout << "WARNING in nInvMatND_svd, divide by 0 fVecW[nx]!\n" << flush;
            return (1);
        }
    
        
        for (ny=1; ny<=nN; ny++)
        {
            fMatV[ny][nx]/=max(1e-20,fVecW[nx]);
        }
    
    };

    if (fANEigenValues) {
        // Load Eigen values from the diagonal values of W
        for (nx=0;nx<nN;nx++)
            fANEigenValues[nx] = fVecW[nx+1]*fMult;
    };
    if (fANxNEigenVecs) {
        // Load Eigen vectors as collumn vectors in V
        for (nx=0;nx<nN;nx++) {
            for (ny=0;ny<nN;ny++) 
                fANxNEigenVecs[nx*nN+ny]= fMatV[ny+1][nx+1];
        };
    };

    
    // Multiply this by [U]t
    for (nx=1;nx<=nN; nx++) 
        for (ny=1; ny<=nN; ny++) {
            f0=0;
            for (nz=1; nz<=nN; nz++)
                f0+=fMatV[nx][nz]*fMatA[ny][nz];    // Note: A:=U after call to svdcmp
            fANxNInv[(nx-1)+(ny-1)*nN]=f0;
        };
    // Undo the affect of the multiplication.
    for (nx=0;nx<nN*nN;nx++)
        fANxNInv[nx]*=fMult;
	return (0);
};



int nSolveAxB_gauss(int nN,float* fANxN,float* fXN,float* fBN) 
{
  // Return 0   on SUCCESS
  //      non-0 on FAILURE

	static double** fMatA=NULL;
	static double** fMatB=NULL;
	static int nMatSize=0;
	int nx,ny;

	if (nN>nMatSize) {
	  //if (bCheck_matrix(fMatA,1,nMatSize,1,nMatSize))
	  if (fMatA)
	    free_matrix(fMatA,1,nMatSize,1,nMatSize);
	  //if (bCheck_matrix(fMatB,1,nMatSize,1,1))
	  if (fMatB)
	    free_matrix(fMatB,1,nMatSize,1,1);
	  fMatA=matrix(1,nN,1,nN);
	  fMatB=matrix(1,nN,1,1);
	  nMatSize=nN;
	};
	float *pfTempA, *pfTempB;
	pfTempA = fANxN;
	pfTempB = fBN;
	for (nx=1;nx<=nN;nx++) {
		for (ny=1;ny<=nN;ny++,pfTempA++) {
			fMatA[ny][nx]=*pfTempA;
		};
		fMatB[nx][1]=*pfTempB;
		pfTempB++;
	};
	if (gaussj(fMatA,nN,fMatB,1))
        return (1);
	for (nx=0;nx<nN;nx++) fXN[nx]=fMatB[nx+1][1];
	return (0);
};

static double ms_afDummyFuncs[ms_nNumDummyFuncs];

double fStandardFunction(double* pfFunction, int nPoint, int nNumPoints)
{
    int nx,ny;
    double f0,f1;

    if( pfFunction >= (double*)0x1 && pfFunction <= (double*)ms_nNumDummyFuncs ) 
    {
        f0 = (double)nPoint/nNumPoints; 
        if( pfFunction >= ms_pfXN )
        {
            nx = (int)(pfFunction - ms_pfXN);
            f1 = 1.0;
            for(ny = 0; ny < nx; ny++)
            {
                f1 *= f0;                   //RB (Thad) raising  nPoint/nNumPoints to the power 
            }
            
            return f1;
        }
        else if( pfFunction == ms_pfSqrtX )
        {
            return sqrt(f0);
        }
        else
            return 0.0;
    }
    else
        return pfFunction[nPoint];
}

int nSetLeastSquaresOutput(char* pcFile,int nMode)
{
    ms_pcLeastSquaresOutputFile = pcFile;
    ms_nLeastSquaresOutputMode = nMode;
    return 0;
}

int nSolveLeastSquares(int nNumFuncs,
                       int nNumPoints,
                       double* pfVals,
                       double* pfValSigmas,
                       double** ppfFuncs,
                       double* pfArgs,
                       double* pfSigmas,
                       int* pnUse,
                       double* pfTotalVariance,
                       itr<int>* panOutlierReject)
{
    int nx,ny,nz;
    int nNumFuncsUsed;
    int nCount;
    int nPass;
    double* pfTemp1;
    double fVariance;
    double fSumWeightsDiffSq;
    double f0;
    double fFunctionValue;      // Returned by the EvaluateFunction() macro below.
    double fMaxDeviation;       // Used to reject outliers.

#define EvaluateFunction(nPoint) { int nx; fFunctionValue = 0.0;\
                                   for (nx = 0; nx< nNumFuncsUsed;nx++)\
                                      fFunctionValue += fStandardFunction(ppfFuncs[pnRelIndex[nx]], nPoint, nNumPoints) * pfArgs[pnRelIndex[nx]]; }


    static double* pfMatA = NULL;
    static double* pfMatAInv = NULL;
    static double* pfMatB = NULL;
    static double* pfMatX = NULL;
    static int* pnRelIndex = NULL;
    static int nMatSize = 0;

    if (nNumFuncs>nMatSize) {
      if (pfMatA) delete[] pfMatA;
      if (pfMatAInv) delete[] pfMatAInv;
      if (pfMatB) delete[] pfMatB;
      if (pfMatX) delete[] pfMatX;
      if (pnRelIndex) delete[] pnRelIndex;
      pfMatA=new double[nNumFuncs*nNumFuncs];
      pfMatAInv = new double[nNumFuncs*nNumFuncs];
      pfMatB=new double[nNumFuncs];
      pfMatX=new double[nNumFuncs];
      pnRelIndex = new int[nNumFuncs];
      nMatSize=nNumFuncs;
    };
    
    // Figure out how many functions we will be using.
    if ( pnUse == NULL)
      {
        for (nx=0;nx<nNumFuncs;nx++)
	  pnRelIndex[nx] = nx;
        nNumFuncsUsed = nNumFuncs;
      } else {
      for (nx=0,nNumFuncsUsed=0;nx<nNumFuncs; nx++) 
	{
	  if (pnUse[nx]==1) 
	    pnRelIndex[nNumFuncsUsed++] = nx;
        };
    };

    for (nPass = 0; nPass < 2; nPass++)
      {
        // We now have the array pnRelIndex[0 .. (nNumFuncsUsed-1)] mapping the function.
        // number into the index of ppfFuincs;
        
        // Build the upper triangular portion of pfMatA

        for (nx=0;nx<nNumFuncsUsed;nx++) {
            for (ny=nx;ny<nNumFuncsUsed;ny++) {
                pfTemp1 = &pfMatA[nx+ny*nNumFuncsUsed];
                *pfTemp1 = 0.0;
            };
            pfTemp1 = &pfMatB[nx];
            *pfTemp1 = 0.0;
        };
        
        for (nz=0; nz<nNumPoints;nz++) {
            if (pfValSigmas == NULL)
                fVariance = 1.0;
            else
                fVariance = pfValSigmas[nz]*pfValSigmas[nz];
            
            // Skip if we are rejecting this point as an outlier.
            if (nPass>0) {
                EvaluateFunction(nz);
                if (fabs(fFunctionValue - pfVals[nz])>fMaxDeviation) {
					if (panOutlierReject)
						(*panOutlierReject)[nz] = 1;
                    continue;
				};
            } else {
	      if (panOutlierReject)
		(*panOutlierReject)[nz] = 0;
	    };
            
            
            for (nx=0;nx<nNumFuncsUsed;nx++) {
                for (ny=nx;ny<nNumFuncsUsed;ny++) {
                    pfTemp1 = &pfMatA[nx+ny*nNumFuncsUsed];
                   
                    *pfTemp1 += fStandardFunction(ppfFuncs[pnRelIndex[nx]],nz,nNumPoints) * 
                        fStandardFunction(ppfFuncs[pnRelIndex[ny]],nz,nNumPoints)/fVariance;
                };
                pfTemp1 = &pfMatB[nx];               
                f0 = pfVals[nz];
                // Add in contributions from functions that are not refining, but have non-zero values.
                if (pnUse) {
                    for (ny=0;ny<nNumFuncs;ny++) {
                        if (pnUse[ny]==2) {
                            f0 -= pfArgs[ny]*fStandardFunction(ppfFuncs[ny],nz,nNumPoints);
                        };
                    };
                };

                *pfTemp1 += f0*fStandardFunction(ppfFuncs[pnRelIndex[nx]],nz,nNumPoints)/fVariance;                        
            };
        };
        
        // Fill in other side of triangular matrix.
        
        for (nx=0;nx<nNumFuncsUsed;nx++) {
            for (ny=nx+1;ny<nNumFuncsUsed;ny++)
                pfMatA[ny + nx*nNumFuncsUsed] = pfMatA[nx+ny*nNumFuncsUsed];;
        };
        
        // Solve the system.
        if (nNumFuncsUsed) {
            nInvMatND_svd(nNumFuncsUsed,pfMatA,pfMatAInv);
            // Figure out the parameters.
            vMulMatNDVecND(nNumFuncsUsed,pfMatAInv,pfMatB,pfMatX);
        };
        
        // Set all parameters equal to zero.
        for (nx=0;nx<nNumFuncs;nx++) {
            if (pfArgs) {
                if ((!pnUse) || (pnUse[nx]!=2)) {
                    pfArgs[nx] = 0.0;
                };
            };
            if (pfSigmas)
                pfSigmas[nx] = 0.0;
        };
        
        // Stuff them back into their places.
        
        for (nx=0;nx<nNumFuncsUsed;nx++) {
            if (pfArgs)
                pfArgs[pnRelIndex[nx]] = pfMatX[nx];
            if (pfSigmas)
                pfSigmas[pnRelIndex[nx]] = 0.0;
        };

        
        // Produce maximum deviation for next pass.
        
        if (nPass==0) {
            double fSumX = 0.0;
            double fSumX2 = 0.0;
            for (nz =0; nz<nNumPoints;nz++) {
                EvaluateFunction(nz);
                fSumX2+=(fFunctionValue - pfVals[nz])*(fFunctionValue - pfVals[nz]);
                fSumX+= fFunctionValue - pfVals[nz];
            };
            fMaxDeviation = ms_fLeastSquaresOutlierRejectIoverSig*sqrt((fSumX2*nNumPoints-fSumX*fSumX)/nNumPoints/(nNumPoints-1));
        };
        if (nNumPoints <= nNumFuncs*5)
            break;
    };


    // Estimate sigma values.
    if (pfSigmas) {

        fSumWeightsDiffSq = 0.0;
        
        nCount = 0;
        for (nz =0; nz<nNumPoints;nz++) {
            if (pfValSigmas == NULL)
                fVariance = 1.0;
            else
                fVariance = pfValSigmas[nz]*pfValSigmas[nz];
            EvaluateFunction(nz);
            if (ABS(fFunctionValue - pfVals[nz])<=fMaxDeviation) {
                fSumWeightsDiffSq += (fFunctionValue - pfVals[nz])*(fFunctionValue - pfVals[nz])/fVariance;
                nCount++;
            };
        };
	if (0 < nCount)
	  fSumWeightsDiffSq/=nCount;
        for (nx=0;nx<nNumFuncsUsed;nx++) {
            pfSigmas[pnRelIndex[nx]] =  sqrt(fabs(pfMatAInv[nx*nNumFuncsUsed+nx])*fSumWeightsDiffSq);
        };

        if (pfTotalVariance)
            *pfTotalVariance = fSumWeightsDiffSq;
    };

    // print out results if there are files names.
    if (ms_pcLeastSquaresOutputFile) {
        char pcBuf[200];
        int* pnSort;
        float* pfValues;
        float* pfValuesInterp;
        FILE* pFOutput;
        
        pnSort = new int[nNumPoints];
        pfValues = new float[nNumPoints];
        pfValuesInterp = new float[nNumPoints];
            nCount = 0;
            for (nz =0; nz<nNumPoints;nz++) {
                if (pfValSigmas == NULL)
                    fVariance = 1.0;
                else
                    fVariance = pfValSigmas[nz]*pfValSigmas[nz];
                EvaluateFunction(nz);
                if (fabs(fFunctionValue - pfVals[nz])<=fMaxDeviation) {
                    pnSort[nCount] = nCount;
                    pfValuesInterp[nCount] = fFunctionValue;
                    pfValues[nCount] = pfVals[nz];
                    nCount++;
                };
            };
            if (ms_nLeastSquaresOutputMode == 0) {
                // Sort the points on the real data.
                g_pfCmpFloats = pfValues;
                qsort(&pnSort[0],nCount,sizeof(int),float_cmp_rel);
            } else if (ms_nLeastSquaresOutputMode == 1) {
                // Sort the points on the interpolated data.
                g_pfCmpFloats = pfValuesInterp;
                qsort(&pnSort[0],nCount,sizeof(int),float_cmp_rel);
            } else if ((ms_nLeastSquaresOutputMode == 2) || (ms_nLeastSquaresOutputMode==3)) {
                // Do not sort the data at all.
                // Use only one plot.
            } else
                printf("Illegal setting for ms_nLeastSquaresOutputMode.\n");


            if ((ms_nLeastSquaresOutputMode == 0) || (ms_nLeastSquaresOutputMode == 1) || (ms_nLeastSquaresOutputMode==2)) {
                // Set the output name.
                strcpy(pcBuf,ms_pcLeastSquaresOutputFile);
                strcat(pcBuf,".txt");
                pFOutput = fopen(pcBuf,"w+t");
                if (pFOutput) {
                    for (nz=0; nz<nCount;nz++) 
                        fprintf(pFOutput,"%d %f\n",nz,pfValues[pnSort[nz]]);
                    fclose(pFOutput);
                };
                
                strcpy(pcBuf,ms_pcLeastSquaresOutputFile);
                strcat(pcBuf,"_interp.txt");
                pFOutput = fopen(pcBuf,"w+t");
                if (pFOutput) {
                    for (nz=0; nz<nCount;nz++) 
                        fprintf(pFOutput,"%d %f\n",nz,pfValuesInterp[pnSort[nz]]);
                    fclose(pFOutput);
                };
            };

            if (ms_nLeastSquaresOutputMode == 3) {
                strcpy(pcBuf,ms_pcLeastSquaresOutputFile);
                strcat(pcBuf,".txt");
                pFOutput = fopen(pcBuf,"w+t");
                    for (nz=0; nz<nCount;nz++) 
                        fprintf(pFOutput,"%f %f\n",pfValues[pnSort[nz]],pfValuesInterp[pnSort[nz]]);
                    fclose(pFOutput);
            };

            delete[] pfValues;
            delete[] pfValuesInterp;
            delete[] pnSort;
    };

    return 0;      
}


// Routine only works for normal matrices (i.e., those that are unit. diagonalizable)
// Correct eigenvalues are returned for non-normal matrices, but the eigenvectors will be bogus.

void  vEigen(int nN,double* anxnfMat,double* anfEigenVals,double* anxnfEigenVecs) {
    static double** fMatA=NULL;
	static double** fMatV=NULL;
	static double*  fVecW=NULL;
	static int nMatSize=0;
	int nx,ny;

	if (nN>nMatSize) {
	  //if (bCheck_matrix(fMatA,1,nMatSize,1,nMatSize))
	  if (fMatA)
	    free_matrix(fMatA,1,nMatSize,1,nMatSize);
	  //if (bCheck_matrix(fMatV,1,nMatSize,1,nMatSize))
	  if (fMatV)
	    free_matrix(fMatV,1,nMatSize,1,nMatSize);
	  if (fVecW) free_vector(fVecW,1,nMatSize);
	  fMatA=matrix(1,nN,1,nN);
	  fMatV=matrix(1,nN,1,nN);
	  fVecW=vectorDT(1,nN);
	  nMatSize=nN;
	};
//+jwp 6-Jun-2001
	double *anxndMatTemp;
	anxndMatTemp = anxnfMat;
	for (nx=1;nx<=nN;nx++)
		for (ny=1;ny<=nN;ny++,anxndMatTemp++) {
			fMatA[ny][nx]=*anxndMatTemp;
		};
//-jwp 6-Jun-2001

	if (svdcmp(fMatA,nN,nN,fVecW,fMatV))
        return;

    if (anfEigenVals) {
        // Load Eigen values from the diagonal values of W
        for (nx=0;nx<nN;nx++)
            anfEigenVals[nx] = fVecW[nx+1];
    };
    if (anxnfEigenVecs) {
        // Load Eigen vectors as collumn vectors in V
        for (nx=0;nx<nN;nx++) {
            for (ny=0;ny<nN;ny++) 
                anxnfEigenVecs[nx*nN+ny]= fMatV[ny+1][nx+1];
        };
    };
    return;
};

int nSmartMatrixCompress(int& nDim,double* pfAMat,double* pfBVec,int* pnStatus, bool bTestOverlaps,int* pnOverlaps,double fPercentOverlap) {

    // nDim             = Dimension of original A matrix.
    // pfAMat, pfBVec   =  A matrix and b vector for the system.
    // pnStatus         = Status vector.  0 = Removed (singular ) 1 = Kept, 2 = Removed (non-singular) scale factor=1.0

    // We must remove the coefficients in the array that 
    // would lead to a singular matrix.  None of these had contributors,
    // so we can assume that they have scale factors of 1.0, and solve for a smaller set of unknowns.
    // We must also choose a set of variables to have values of 1.0.
    // If we sent in a zero B vector, this will cause the B vector to be non-zero.
    // Let F[x] be a heuristic function defining the "goodness" of a variable, 
    // Let {X} be a set of variables such that for each x in {X}, F[x] is maximal in it's connected component.    
    // Thus, we reduce the matrix by removing all columns y having ((F[y]==0) or (y in {X})).  
    // These collumns are all added together and negated to form the B vector.  For the resulting system,
    // all rows y having ((F[y]==0) or (y in {X})) are removed to form a square A matrix.  Similarly, the B
    // vector is reduced.

    // The heuristic function F[x] is defined as the total number of overlaps of a variable with other
    // variables in the group.
    // Groups are chosen as follows:  All members in the group intersect some other element of the group
    // with at least fPercentOverlaps overlaps.
    // If pnOverlaps is NULL, groups are chossen so that no block diagonals are formed.

    int nx,ny;
    double f0,f1;
    
    double* pfCounts;    // Sum of absolute values of entries in Row for each variable.
    int* pnConnect;     // Array used when finding the representitives.
    int  nNumReps;      // The number of representitives.
    int  nRep;

    // Flags
    int nFlagRemoved = 0;           
    int nFlagUnconsidered =1;
    int nFlagComponant = 2;
    int nFlagRemovedSingular = 3;
      

    pfCounts = new double[nDim];
    pnConnect = new int[nDim];

#define GETNAME(name,ni,nj) name[(nj)*nDim+(ni)] // Data stored in collumns
#define GETNAMED(name,d0,ni,nj) name[(nj)*(d0)+(ni)]

    // Compute pfCounts.
    for (nx=0;nx<nDim;nx++) {
        f0 = 0.0, f1=0.0;
        for (ny=0;ny<nDim;ny++)
            f0 += fabs(GETNAME(pfAMat,nx,ny));
        for (ny=0;ny<nDim;ny++)
            f1 += fabs(GETNAME(pfAMat,ny,nx));
        pfCounts[nx] = f0*f1;
    };

    // We must find a set of representitives from each connected componant.
    // The componants are selected sequentially.
    
    nNumReps = 0;

    // Mark each variable that has pfCounts[]==0.0 with an nFlagRemoved,
    for (nx=0;nx<nDim;nx++) {
        if ((pnOverlaps!=NULL) && (GETNAME(pnOverlaps,nx,nx)==0)) {
            pnConnect[nx] = nFlagRemovedSingular;
            pnStatus[nx] = 0;   // We know the status of this variable right now.                                
        } else
            if (bTestOverlaps)
                // If we are testing overlaps, we need yet determine if the variable
                // has a status of 1 or 2                
                pnConnect[nx] = nFlagUnconsidered;               
            else
                pnStatus[nx] = 1;
    };
    if (bTestOverlaps) while (1) {
        // Find a vertex that has not been considered yet.
        // Preferably, we choose a point with a maximal number of self refereneces,
        // or (if pnOverlaps==NULL) a variable with a maximal diagonal value.
        for (nx=0,nRep=-1;nx<nDim;nx++) {
            if (pnConnect[nx]==nFlagUnconsidered) {
                if (nRep==-1) 
                    nRep = nx;
                else if (pnOverlaps==NULL) {
                    if (GETNAME(pfAMat,nx,nx)>GETNAME(pfAMat,nRep,nRep))
                        nRep = nx; 
                } else {
                    if (GETNAME(pnOverlaps,nx,nx)>GETNAME(pnOverlaps,nRep,nRep))
                        nRep = nx;
                };
            };
        };

        // Did we find an unconsidered point?  If not, then we are done.
        if (nRep== -1)
            break;
        pnConnect[nRep] = nFlagComponant;

        // Begin finding neighbors for this new point from the remaining unconsidered points.
        // If any of these neighbors have been removed, then we are in trouble.
        bool bNewFound;
        do {
            bNewFound = FALSE;
            for (nx=0;nx<nDim;nx++) {
                if (pnConnect[nx] == nFlagComponant) {
                    for (ny=0;ny<nDim;ny++) {
                        if (pnConnect[ny]==nFlagUnconsidered) {
                            //We want the variable [ny] added to the componant only if
                            //it has an acceptable number of intersections with variable [nx]
                            if (pnOverlaps==NULL) {
                                if (((GETNAME(pfAMat,nx,ny)!=0.0) || (GETNAME(pfAMat,ny,nx)!=0.0))) {
                                    pnConnect[ny] = nFlagComponant;
                                    bNewFound = TRUE;
                                };
                            } else {
                                if (((double) (GETNAME(pnOverlaps,nx,ny)+GETNAME(pnOverlaps,ny,nx)))
                                    >=fPercentOverlap*((double) GETNAME(pnOverlaps,ny,ny))) {
                                    pnConnect[ny] = nFlagComponant;
                                    bNewFound = TRUE;
                                };
                            };
                        };
                    };
                };
            };
        } while (bNewFound);

        nNumReps++;

        // Finally, we set all variables in this componant to "removed"
        // We write out their status into the pnStatus array
        for (nx=0;nx<nDim;nx++) {
            if (pnConnect[nx] == nFlagComponant) {
                pnConnect[nx] = nFlagRemoved;
                if (nx==nRep)
                    pnStatus[nx] = 2;
                else
                    pnStatus[nx] = 1;
            };
        };
    };

    // Now we do the actual matrix compression.  This is producing a matrix that has more rows than collumns.
    // All variables that have pnStatus != 1 are compressed.
    // However, only singular variables are actually removed... All others will be compressed.

    int nParamIn0,nParamIn1;    // [0] = row, [1] = column
    int nParamOut0,nParamOut1;
    int nDimOut;
    double* pfAMatOut;


    for (nParamIn1= 0,nParamOut1= 0; nParamIn1<nDim; nParamIn1++)  {
        if (pnStatus[nParamIn1]==1) 
            nParamOut1++;
    };
    nDimOut = nParamOut1;
    pfAMatOut = new double[nDimOut*nDimOut];    

    for (nParamIn0 = 0, nParamOut0 = 0; nParamIn0 < nDim; nParamIn0++) {
        if (pnStatus[nParamIn0]==1) {
            pfBVec[nParamOut0] = pfBVec[nParamIn0];
            for (nParamIn1= 0,nParamOut1= 0; nParamIn1<nDim; nParamIn1++)  {
                if (pnStatus[nParamIn1]==1) {
                    GETNAMED(pfAMatOut,nDimOut,nParamOut0,nParamOut1)=GETNAME(pfAMat,nParamIn0,nParamIn1);
                    nParamOut1++;
                } else if (pnStatus[nParamIn1]!=0)
                    pfBVec[nParamOut0]-=GETNAME(pfAMat,nParamIn0,nParamIn1);
            };
            nParamOut0++;
        };
    };
    
    vCopyVecND(nParamOut1*nParamOut1,pfAMatOut,pfAMat);
    delete[] pfAMatOut;
    nDim = nParamOut1;
    
    return nNumReps;
};



int nInitQuadratic(double a9x9fCoeffA[9][9],double a9fCoeffb[9])
{
    // Initialize the coefficient matrices
    vZeroMat(9,9,(double*) a9x9fCoeffA);
    
    for(int nx=0; nx < 9; nx++) 
        a9fCoeffb[nx]=0.0;
    
    return 0;        
}

int nAddQuadratic(double fVal,
                  double fSigma,
                  double a3fVec[3],
                  double a9x9fCoeffA[9][9],
                  double a9fCoeffb[9])
{
    int nVar0,nVar1,nVarCt;
    int nLoop0,nLoop1,nLoopCt;
    double a3fNormVec[3];

    vCopyVec3D((double*) a3fVec,(double*) a3fNormVec);
    fNormVec3D(a3fNormVec);

    for (nVar0=0;nVar0<3;nVar0++)
        for (nVar1=0;nVar1<3;nVar1++) {
            // nVarCt maps (nVar0,nVar1) into 0..8
            nVarCt = nVar0 + nVar1*3;
            // constant coeff.
            a9fCoeffb[nVarCt] += a3fNormVec[nVar0]*a3fNormVec[nVar1]*fVal/(fSigma*fSigma);
            for (nLoop0=0;nLoop0<3;nLoop0++)
                for (nLoop1=0;nLoop1<3;nLoop1++) {
                    nLoopCt = nLoop0 + nLoop1*3;
                    a9x9fCoeffA[nLoopCt][nVarCt]+=a3fNormVec[nLoop0]*a3fNormVec[nLoop1]*a3fNormVec[nVar0]*a3fNormVec[nVar1]/(fSigma*fSigma);
                };
        };
    return 0;
};

int nSolveQuadratic(double a3x3fQuad[3][3], double a9x9fCoeffA[9][9], double a9fCoeffb[9])
{
    int nx,ny;

    float a6x6fCoeffA[6][6];
    float a6fCoeffb[6];
    float a6fSol[6];


    /* Convert the system to 6 variables instead of 9.

    0 = 0
    1 = 4
    2 = 8
    3 = 1 + 3
    4 = 2 + 6
    5 = 5 + 7
    */
   
    int nMap1[6]= { 0,  4,  8,  1,  2,  5};
    int nMap2[6]= { -1, -1, -1, 3, 6, 7};
    int nMapCt;

    nMapCt=0;
    for (nx=0;nx<6;nx++) {
        for (ny=0;ny<6;ny++) {
            a6x6fCoeffA[nx][ny]=a9x9fCoeffA[nMap1[nx]][nMap1[ny]];
            nMapCt++;
            if (nMap2[nx]!=-1) {
                a6x6fCoeffA[nx][ny]+=a9x9fCoeffA[nMap2[nx]][nMap1[ny]];
                nMapCt++;
            };
            if (nMap2[ny]!=-1) {
                a6x6fCoeffA[nx][ny]+=a9x9fCoeffA[nMap1[nx]][nMap2[ny]];
                nMapCt++;
            };
            if ((nMap2[nx]!=-1) && (nMap2[ny]!=-1)) {
                a6x6fCoeffA[nx][ny]+=a9x9fCoeffA[nMap2[nx]][nMap2[ny]];
                nMapCt++;
            };
        };
        a6fCoeffb[nx]=a9fCoeffb[nMap1[nx]];
        if (nMap2[nx]!=-1)
            a6fCoeffb[nx]+=a9fCoeffb[nMap2[nx]];
    };

  // solve the system using svd.
    if (nSolveAxB_svd(6,(float*) a6x6fCoeffA,(float*) a6fSol,(float*) a6fCoeffb))
        return 1;

    a3x3fQuad[0][0]=a6fSol[0];
    a3x3fQuad[1][1]=a6fSol[1];
    a3x3fQuad[2][2]=a6fSol[2];
    a3x3fQuad[1][0]=a3x3fQuad[0][1]=a6fSol[3];
    a3x3fQuad[2][0]=a3x3fQuad[0][2]=a6fSol[4];
    a3x3fQuad[1][2]=a3x3fQuad[2][1]=a6fSol[5];

    return 0;
};



#define QUADMATG(x,y) ((y)*nN+(x))          // Gives index for the G[x,y] element in the A matrix and in the b vector
#define QUADMATH(x) (nN*nN+(x))             // Gives index for the h[x] element in the A matrix and b vector
#define QUADMATC (nN*nN + nN)               // Gives index for the constant element in the A matrix and b vector.
#define QUADMATA(x,y) ((y)*nNumVars+(x))    // Gives index for A[x,y] element in A[[nN*(nN+1)+1][nN*(nN+1)+1]]

int nInitQuadratic(int nN,
                   bool bAllocateMemory,
                   double** pfCoeffA /* [nN*(nN+1)+1][nN*(nN+1)+1] */,
                   double** pfCoeffb /* [nN*(nN+1)+1] */)
{
    int nx;

    if (bAllocateMemory) {
        delete[] *pfCoeffA;
        delete[] *pfCoeffb;
        *pfCoeffA = new double[(nN*(nN+1)+1)*(nN*(nN+1)+1)];
        *pfCoeffb = new double[nN*(nN+1)+1];
    };
    // Initialize the coefficient matrices
    for (nx=0;nx<nN*(nN+1)+1;nx++) 
        (*pfCoeffb)[nx]=0.0;
    for (nx=0;nx<(nN*(nN+1)+1)*(nN*(nN+1)+1);nx++) 
        (*pfCoeffA)[nx] = 0.0;
    return 0;        
};

int nAddQuadratic(int nN,
                  double fVal,
                  double fSigma,
                  double* pfVec /* [nN] */,
                  double* pfCoeffA /* [nN*(nN+1)+1][nN*(nN+1)+1] */,
                  double* pfCoeffb  /* [nN*(nN+1)+1] */)
{
    int nVar0,nVar1;
    int nLoop0,nLoop1;
    int nNumVars = nN*(nN+1) + 1;

    /*  The variables are ordered so that the first nN*nN variables are for the [G] matrix coefficents in f(x) = xt*[G]*x + ht*x 
        The next nN variables are for the h vector.
        The last variable is for the constant.
    */

    // Derivatives w.r.t. G variables.
    for (nVar0=0;nVar0<nN;nVar0++) {
        for (nVar1=0;nVar1<nN;nVar1++) {
            // b vector
            pfCoeffb[QUADMATG(nVar0,nVar1)] += pfVec[nVar0]*pfVec[nVar1]*fVal/(fSigma*fSigma);
            // derivatives of G variables w.r.t. G[nVar0,nVar1] 
            for (nLoop0=0;nLoop0<nN;nLoop0++) {
                for (nLoop1=0;nLoop1<nN;nLoop1++) {
                    pfCoeffA[QUADMATA(QUADMATG(nVar0,nVar1),QUADMATG(nLoop0,nLoop1))]+=(pfVec[nLoop0]*pfVec[nLoop1]*pfVec[nVar0]*pfVec[nVar1])/(fSigma*fSigma);
                };
            };
            // derivatives of h variables w.r.t. G[nVar0,nVar1]
            for (nLoop0=0;nLoop0<nN;nLoop0++) {
                pfCoeffA[QUADMATA(QUADMATG(nVar0,nVar1),QUADMATH(nLoop0))]+=(pfVec[nLoop0]*pfVec[nVar0]*pfVec[nVar1])/(fSigma*fSigma);
            };
            // derivatives of const variable w.r.t. G[nVar0,nVar1]
            pfCoeffA[QUADMATA(QUADMATG(nVar0,nVar1),QUADMATC)]+=(pfVec[nVar0]*pfVec[nVar1])/(fSigma*fSigma);
        };
    };
    // Derivatives w.r.t. h variables.
    for (nVar0=0;nVar0<nN;nVar0++) {
        // b vector
        pfCoeffb[QUADMATH(nVar0)]+= pfVec[nVar0]*fVal/(fSigma*fSigma);
        // derivatives of G variables w.r.t. h[nVar0]
        for (nLoop0=0;nLoop0<nN;nLoop0++) {
            for (nLoop1=0;nLoop1<nN;nLoop1++) {
                pfCoeffA[QUADMATA(QUADMATH(nVar0),QUADMATG(nLoop0,nLoop1))]+=(pfVec[nLoop0]*pfVec[nLoop1]*pfVec[nVar0])/(fSigma*fSigma);
            };
        };
        // derivatives of h variables w.r.t. h[nVar0]
        for (nLoop0=0;nLoop0<nN;nLoop0++) {
            pfCoeffA[QUADMATA(QUADMATH(nVar0),QUADMATH(nLoop0))]+=(pfVec[nLoop0]*pfVec[nVar0])/(fSigma*fSigma);
        };
        // derivatives of constant variable w.r.t. h[nVar0]
        pfCoeffA[QUADMATA(QUADMATH(nVar0),QUADMATC)]+=(pfVec[nVar0])/(fSigma*fSigma);

    };
    // Derivatives w.r.t c variable.
    // b vector
    pfCoeffb[QUADMATC]+= fVal/(fSigma*fSigma);
    // derivatives of G variables w.r.t. constant
    for (nLoop0=0;nLoop0<nN;nLoop0++) {
        for (nLoop1=0;nLoop1<nN;nLoop1++) {
            pfCoeffA[QUADMATA(QUADMATC,QUADMATG(nLoop0,nLoop1))]+=(pfVec[nLoop0]*pfVec[nLoop1])/(fSigma*fSigma);
        };
    };
    // derivatives of h variables w.r.t. constant
    for (nLoop0=0;nLoop0<nN;nLoop0++) {
        pfCoeffA[QUADMATA(QUADMATC,QUADMATH(nLoop0))]+=(pfVec[nLoop0])/(fSigma*fSigma);
    };
    // derivatives of constant variable w.r.t. constant
    pfCoeffA[QUADMATA(QUADMATC,QUADMATC)]+=1.0/(fSigma*fSigma);

    return 0;
};

int nSolveQuadratic(int nN,
                    double* pfGMat,
                    double* pfhVec,
                    double* pfConst,
                    double* pfCoeffA,
                    double* pfCoeffb)
{
    int nx,ny;

    static double* pfCoeffSol = NULL;
    int nSolutionSize = 0;

    if (nN>nSolutionSize) {
        delete[] pfCoeffSol;
        pfCoeffSol = new double[nN*(nN+1)+1];
        nSolutionSize = nN;
    };

    // solve the system using svd.
    if (nSolveAxB_svd(nN*(nN+1)+1,pfCoeffA,pfCoeffSol,pfCoeffb))
        return 1;
    // Copy the answers back into the correct matrix elements.
    for (nx=0;nx<nN;nx++) {
        for (ny=0;ny<nN;ny++) {
            pfGMat[QUADMATG(nx,ny)] = 0.5*(pfCoeffSol[QUADMATG(nx,ny)]+pfCoeffSol[QUADMATG(ny,nx)]);
        };
        pfhVec[nx] = pfCoeffSol[QUADMATH(nx)];
    };
    *pfConst = pfCoeffSol[QUADMATC];
    return 0;
};


double fEvalQuadratic(int nN,double* pfCoeffA,double* pfCoeffb,double fCoeffC,double* pfVector) {
    double fSum;
    int nx,ny;

    fSum = 0;
    for (nx=0;nx<nN;nx++) {
        for (ny=0;ny<nN;ny++) {
            fSum += pfCoeffA[nx + nN*ny]*pfVector[nx]*pfVector[ny];
        };
        if (pfCoeffb)
            fSum += pfCoeffb[nx]*pfVector[nx];
    };
    fSum += fCoeffC;
    return fSum;
};



int test_svd(int nDim)
{
	int nx,ny;
	float* fMat;
	float* fB;
	float* fX;
	float* fBSolve;
	float f0;
	fMat = new float[nDim*nDim];
	fB = new float[nDim];
	fBSolve = new float[nDim];
	fX = new float[nDim];

    /*
	// This generates the Vandermond matrix 
	for (nx=0;nx<nDim;nx++)
		for (ny=0;ny<nDim;ny++) {
			f0=fMat[nx*nDim+ny]=(ny+1);		
			for (nz=0;nz<nx;nz++) fMat[nx*nDim+ny]*=f0;
		};
    */
	// This generates a near singular matrix 
	for (nx=0;nx<nDim;nx++) {
		for (ny=0;ny<nDim;ny++) {
			fMat[nx*nDim+ny]=1000+(rand() % 20);
		};
	};

	for (nx=0;nx<nDim;nx++) fB[nx]=(nx+3);
	printf("\n----Initial B Vector----\n");
	for (nx=0;nx<nDim;nx++) printf("%f ",fB[nx]);

	// Print the matrix 
	//
	//printf("\n----Matrix----\n");
	//for (ny=0;ny<nDim;ny++) {
	//	for (nx=0;nx<nDim;nx++) {
	//		printf("%8.8g ",fMat[nx*nDim+ny]);
	//	};
	//	printf("\n");
	//};
	//printf("\n");
	//
	
	nSolveAxB_svd(nDim,fMat,fX,fB);
	vMulMatNDVecND(nDim,fMat,fX,fBSolve);
	printf("\n----Solution B Vector (SVD)----\n");
	for (nx=0;nx<nDim;nx++) printf("%f ",fBSolve[nx]);
	for (nx=0,f0=0.0;nx<nDim;nx++) f0+=fabs(fBSolve[nx]-fB[nx]);
	printf("\nResid: %g\n",f0);


	printf("\n----Solution B Vector (GAUSS)----\n");
	nSolveAxB_gauss(nDim,fMat,fX,fB);
	vMulMatNDVecND(nDim,fMat,fX,fBSolve);
	for (nx=0;nx<nDim;nx++) printf("%f ",fBSolve[nx]);
	for (nx=0,f0=0.0;nx<nDim;nx++) f0+=fabs(fBSolve[nx]-fB[nx]);
	printf("\nResid: %g\n",f0);

	delete[] fMat;
	delete[] fB;
	delete[] fBSolve;
	delete[] fX;
    return 0;
};

  


float fScaleMatMN(int nN,int nM,float* fANxN,bool bApply) {
    float fAverage;
    int nx;
    fAverage = 0.0;
    for (nx=0;nx<nN*nM;nx++) {
            fAverage+=(float)fabs((double)fANxN[nx]);
    };
    if (fAverage==0.0)
        fAverage=1.0;
    else
        fAverage/=(nN*nM);
    if (bApply) {
        for (nx=0;nx<nN*nM;nx++)
            fANxN[nx]/=fAverage;
    };
    return 1.0/fAverage;
};


double fScaleMatMN(int nN,int nM,double* fANxN,bool bApply) {
    double fAverage;
    int nx;
    fAverage = 0.0;
    for (nx=0;nx<nN*nM;nx++) {
            fAverage+=fabs(fANxN[nx]);
    };
    if (fAverage==0.0)
        fAverage=1.0;
    else
        fAverage/=(nN*nM);
    if (bApply) {
        for (nx=0;nx<nN*nM;nx++)
            fANxN[nx]/=fAverage;
    };
    return 1.0/fAverage;
};

bool
bCheck_matrix(double **m, int nrl, int nrh, int ncl, int nch)
{
   int i;

   for(i = nrh; i >= nrl; i--)
     {
       if ( NULL == (m[i]+ncl) )
	 {
	   cout << "NULL pointer found in a scratch matrix." << flush;
	   return (TRUE);
	 }
     }
   return (FALSE);
}
/*
double ms_afDummyFuncs[ms_nNumDummyFuncs];
double* ms_pfConst = & ms_afDummyFuncs[0];
double* ms_pfX = & ms_afDummyFuncs[1];
double* ms_pfX2 = & ms_afDummyFuncs[2];
double* ms_pfSqrtX = & ms_afDummyFuncs[3];
char*   ms_pcLeastSquaresOutputFile = NULL;
int     ms_nLeastSquaresOutputMode = 0;
*/


int nBisection(int& nState,double fMin,double fMax,double& fTest,double fTestResult,double fResultTarget,double fTestTol) {
    static double afTestBounds[2];
    static double afTestResult[2];


    if (nState==0) {
        fTest = fMin;
        afTestBounds[0] = fMin;
        afTestBounds[1] = fMax;

    } else if (nState == 1) {
        fTest = fMax;
        afTestResult[0] = fTestResult;

    } else if (nState == 2) {
        afTestResult[1] = fTestResult;
        // We always want the lesser test result to be in [0], and the greater in [1]
        if (afTestResult[0] > afTestResult[1]) {
            std::swap(afTestBounds[0],afTestBounds[1]);
            std::swap(afTestResult[0],afTestResult[1]);
        };
        // If the user provided a target value, then make sure we have framed the value.
        if ((fResultTarget < afTestResult[0]) || (fResultTarget > afTestResult[1]))
            return 2;
        fTest = (afTestBounds[0] + afTestBounds[1])*0.5;

    } else {
        // In case it got changed.
        fTest = (afTestBounds[0] + afTestBounds[1])*0.5;
        
        if (ABS(fTest - afTestBounds[0])<fTestTol)
            return 1;

        // Examine last test value.
        if (fTestResult < fResultTarget) {
            afTestBounds[0] = fTest;
            afTestResult[0] = fTestResult;
        } else {
            afTestBounds[1] = fTest;
            afTestResult[1] = fTestResult;
        };
        
        fTest = (afTestBounds[0] + afTestBounds[1])*0.5;
    };
    nState++;
    return 0;
};


itr<double> gs_aafInterpData[9];

int nGridInterpClear(int nEntries) {
    int nx;
    for (nx=0;nx<9;nx++) {
        gs_aafInterpData[nx].setsize(nEntries);
        memset(&gs_aafInterpData[nx][0],0,sizeof(double)*nEntries);
    };
    return 0;
};

int nGridInterpAdd(double fPos0,double fPos1,int nEntries,double* pfData) {
    int nx;

    if (gs_aafInterpData[0].size() != nEntries) {
        return 1;
    };
    for (nx=0;nx<nEntries;nx++) {
        gs_aafInterpData[0][nx] += fPos0*fPos0;
        gs_aafInterpData[1][nx] += fPos0*fPos1;
        gs_aafInterpData[2][nx] += fPos1*fPos1;
        gs_aafInterpData[3][nx] += fPos0;
        gs_aafInterpData[4][nx] += fPos1;
        gs_aafInterpData[5][nx] += 1.0;
        gs_aafInterpData[6][nx] += pfData[nx]*fPos0;
        gs_aafInterpData[7][nx] += pfData[nx]*fPos1;
        gs_aafInterpData[8][nx] += pfData[nx];
    };
    return 0;
};

int nGridInterpSolve(double fPos0,double fPos1,int nEntries,double* pfData) {
    double a3x3fMatA[3][3];
    double a3x3fMatAInv[3][3];
    double a3fMatb[3];
    double a3fBackground[3];
    int nx;

    // Just in case we don't return.
    vZeroMat(nEntries,1,pfData);
    
    if ((!nEntries) || (gs_aafInterpData[0].size()!=nEntries) || (gs_aafInterpData[5][0]<=2.0))
        return 1;
    
    for (nx=0;nx<nEntries;nx++) {    
        vZeroMat(3,3,&a3x3fMatA[0][0]);
        vZeroMat(3,1,&a3fMatb[0]);
        
        a3x3fMatA[0][0] = gs_aafInterpData[0][nx];
        a3x3fMatA[1][0] = gs_aafInterpData[1][nx];
        a3x3fMatA[2][0] = gs_aafInterpData[3][nx];
        a3x3fMatA[1][1] = gs_aafInterpData[2][nx];
        a3x3fMatA[2][1] = gs_aafInterpData[4][nx];
        a3x3fMatA[2][2] = gs_aafInterpData[5][nx];
        a3x3fMatA[0][1] = a3x3fMatA[1][0];
        a3x3fMatA[0][2] = a3x3fMatA[2][0];
        a3x3fMatA[1][2] = a3x3fMatA[2][1];
        a3fMatb[0] = gs_aafInterpData[6][nx];
        a3fMatb[1] = gs_aafInterpData[7][nx];
        a3fMatb[2] = gs_aafInterpData[8][nx];
        
        if (0.0 == fInvMat3D(&a3x3fMatA[0][0],&a3x3fMatAInv[0][0])) {
            return 1;
        };   
        vMulMat3DVec3D(a3x3fMatAInv,a3fMatb,a3fBackground);
        pfData[nx] = a3fBackground[0]*fPos0 + a3fBackground[1]*fPos1 + a3fBackground[2];
    };
    return 0;
};


int nConvolve(itr<double>& afData,
              itr<double>& afResp,
              bool bConvolve,
              itr<double>& afResult)
{
    int nDataSize,nRespSize;
    itr<double> afTemp;
    int nSizeDiv2;
    double f0,f1;
    double fMin;
    int nx;
    
    nDataSize = afData.size();
    nRespSize = afResp.size();

    afTemp.setsize(nDataSize);
    afResult.setsize(nDataSize);
    memset(&afTemp[0],0,sizeof(double)*nDataSize);
    afTemp[0] = afResp[0];
    for (nx = 1; nx<(nRespSize+1)/2;nx++) {
        afTemp[nx] = afResp[nx];
        afTemp[nDataSize - nx] = afResp[nRespSize - nx];
    };
    -afResult + afData;

    nRealFFT1D(nDataSize,1,&afResult[0]);
    nRealFFT1D(nDataSize,1,&afTemp[0]);
    nSizeDiv2 = nDataSize >> 1;
    if (bConvolve) {
        for (nx=2;nx<nDataSize;nx+=2) {
            f0 = afResult[nx];
            afResult[nx] = (afResult[nx]*afTemp[nx] - afResult[nx+1]*afTemp[nx+1])/nSizeDiv2;
            afResult[nx+1] = (afResult[nx+1]*afTemp[nx] + f0*afTemp[nx+1])/nSizeDiv2;
        };
        afResult[0] = afResult[0]*afTemp[0]/nSizeDiv2;
        afResult[1] = afResult[1]*afTemp[1]/nSizeDiv2;
    } else {
        fMin = 1000;
        for (nx=2;nx<nDataSize;nx+=2) {
            f0 = afResult[nx];
            f1 = afTemp[nx]*afTemp[nx] + afTemp[nx+1]*afTemp[nx+1];
            if (f1 < fMin)
                fMin = f1;
            afResult[nx] = (afResult[nx]*afTemp[nx] + afResult[nx+1]*afTemp[nx+1])/f1/nSizeDiv2;
            afResult[nx+1] = (afResult[nx+1]*afTemp[nx] - f0*afTemp[nx+1])/f1/nSizeDiv2;
        };
        afResult[0] = afResult[0]/afTemp[0]/nSizeDiv2;
        afResult[1] = afResult[1]/afTemp[1]/nSizeDiv2;
    };

    nRealFFT1D(nDataSize,-1,&afResult[0]);
    return 0;
};
//////////////////////////////////////////////////////////////////////////////////////////////////////
//Given a set of data points x[1..ndata],y[1..ndata] with individual standard deviations
//sig[1..ndata], fit them to a straight line y = a + bx by minimizing chi-2. Returned are
//a,b and their respective probable uncertainties siga and sigb, the chi-square chi2, and the
//goodness-of-fit probability q (that the fit would have chi-2 this large or larger). If mwt=0 on
//input, then the standard deviations are assumed to be unavailable: q is returned as 1.0 and
//the normalization of chi2 is to unit standard deviation on all points.
bool bLinearFit(float x[], 
                float y[], 
                int ndata, 
                float sig[], 
                bool mwt, 
                float *a,
                float *b, 
                float *siga, 
                float *sigb, 
                float *chi2, 
                float *q)
{
    //FILE*   fp = fopen("C:\\Temp\\junk_bins.txt","w");
    //for(int i=1;i<ndata;i++)
      //fprintf(fp,"%f %f\n", x[i], y[i]);
    //fclose(fp);
    
    int     i=0;
	float   wt=0.0f, t=0.0f, sxoss=0.0f, 
            sx=0.0f, sy=0.0f, st2=0.0f, 
            ss=0.0f, sigdat=0.0f;

	*b=0.0;
	
    if( mwt )
    {
		ss=0.0;
		for (i=1;i<=ndata;i++) 

        {
			wt=1.0/pow(sig[i], 2.0f);
			ss += wt;
			sx += x[i]*wt;
			sy += y[i]*wt;
        }
	} 
    else
    {
		for (i=1;i<=ndata;i++)
        {
			sx += x[i];
			sy += y[i];
		}
		ss=ndata;
	}
	sxoss=sx/ss;
	if (mwt)
    {
		for (i=1;i<=ndata;i++)
        {
			t=(x[i]-sxoss)/sig[i];
			st2 += t*t;
			*b += t*y[i]/sig[i];
		}
	} 
    else
    {
		for (i=1;i<=ndata;i++)
        {
			t=x[i]-sxoss;
			st2 += t*t;
			*b += t*y[i];
		}
	}
	*b /= st2;
	*a=(sy-sx*(*b))/ss;
	*siga=sqrt((1.0f+sx*sx/(ss*st2))/ss);
	*sigb=sqrt(1.0f/st2);
	*chi2=0.0f;
	if (mwt == 0) {
		for (i=1;i<=ndata;i++)
			*chi2 += pow((y[i]-(*a)-(*b)*x[i]), 2.0f);
		*q=1.0;
		sigdat=sqrt((*chi2)/(ndata-2));
		*siga *= sigdat;
		*sigb *= sigdat;
	} 
    else
    {
		for (i=1;i<=ndata;i++)
			*chi2 += pow(((y[i]-(*a)-(*b)*x[i])/sig[i]), 2.0f);
		
        //*q=gammq(0.5*(ndata-2),0.5*(*chi2));
        if( !gammq(0.5f*(ndata-2), 0.5f*(*chi2), *q) )
            return false;
	}
    
    return true;
}
/////////////////////////////////////////////////////////////////////
bool gammq(float a, float x, float& fOut)
{
	if( x < 0.0 || a <= 0.0 ) 
    {
        printf("\nERROR: Invalid arguments in routine gammq\n");
        return false;
    }
    
    float gamser=0.0f,gammcf=0.0f,gln=0.0f;

    if (x < (a+1.0f))
    {
		if( !gser(&gamser,a,x,&gln) )
            return false;

		fOut = 1.0f - gamser;
	}
    else
    {
		if( !gcf(&gammcf,a,x,&gln) )
            return false;
		
        fOut = gammcf;
	}
    
    return true;
}
/////////////////////////////////////////////////////////////////////
#define ITMAX 100
#define EPS 3.0e-7

bool gser(float *gamser, float a, float x, float *gln)
{
	int n=0;
	float sum,del,ap;

	*gln=gammln(a);
	if( x <= 0.0f ) 
    {
        *gamser = 0.0f;
		
        if( x < 0.0 ) 
        {
            printf("ERROR: x less than 0 in routine gser");
            return false;
        }
        else
            return true;
	} 
    else
    {
		ap=a;
		del=sum=1.0f/a;
		for (n=1;n<=ITMAX;n++)
        {
			++ap;
			del *= x/ap;
			sum += del;
			if (fabs(del) < fabs(sum)*EPS)
            {
				*gamser=sum*exp(-x+a*log(x)-(*gln));
				
                return true;
			}
		}
		
        printf("ERROR: a too large, ITMAX too small in routine gser");
		
        return false;
	}
    
    return false;
}
/////////////////////////////////////////////////////////////////
float gammln(float xx)
{
	double x,y,tmp,ser;
	static double cof[6]={76.18009172947146,-86.50532032941677,
		24.01409824083091,-1.231739572450155,
		0.1208650973866179e-2,-0.5395239384953e-5};
	int j;

	y=x=xx;
	tmp=x+5.5;
	tmp -= (x+0.5)*log(tmp);
	ser=1.000000000190015;
	
    for (j=0;j<=5;j++)
        ser += cof[j]/++y;
	
    return (float)(-tmp+log(2.5066282746310005*ser/x));
}
//////////////////////////////////////////////////////////////////
#define FPMIN 1.0e-30
bool gcf(float *gammcf, float a, float x, float *gln)
{
	int i=0;
	float an=0.0f, b=0.0f, c=0.0f, d=0.0f, del=0.0f, h=0.0f;

	*gln=gammln(a);
	b=x+1.0-a;
	c=1.0f/FPMIN;
	d=1.0f/b;
	h=d;
	for (i=1;i<=ITMAX;i++)
    {
		an = -i*(i-a);
		b += 2.0f;
		d=an*d+b;
		
        if (fabs(d) < FPMIN)
            d=FPMIN;
		
        c=b+an/c;
		
        if (fabs(c) < FPMIN)
            c=FPMIN;
		
        d=1.0f/d;
		del=d*c;
		h *= del;
		
        if (fabs(del-1.0f) < EPS) 
            break;
	}
	if( i > ITMAX ) 
    {
        printf("ERROR: a too large, ITMAX too small in gcf");
        return false;
    }
	
    *gammcf=exp(-x+a*log(x)-(*gln))*h;
    
    return true;
}
///////////////////////////////////////////////////////////////////////
