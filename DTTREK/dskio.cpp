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
// Cimage.cc            Initial author: J.W. Pflugrath           03-Mar-1995
//    This file contains the member functions of class Cimage.
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
//+Description
//  dskio.c  has some file i/o operations in it.  They will be replaced as
//           we develop more code.
//+ToDo
//
//   Error messages need implementing
//
//+Include files

#include <stdio.h>
#include <string.h>
#include "dskio.h"

#define dskbor_ dskbor
#define dskbcr_ dskbcr
#define dskbr_  dskbr
#define dskbwr_ dskbwr
#define dskbow_ dskbow
#define dskbcw_ dskbcw
#define dskbw_  dskbw
#define dskbww_ dskbww

static FILE *fp[4] = {NULL, NULL, NULL, NULL};

/*
 * ***********************************
 * dskbor - Byte file Open for Reading
 * ***********************************
 */
extern "C"
#ifdef VAX11C
void dskbor(lun, filename, lfilename, istat)
     struct descr *filename;
     int *lun, *lfilename, *istat;
#else
void dskbor_(int *lun, char *filename, int *lfilename, int *istat)
#endif
{
//  FILE *fopen();
  char *temp;
  int  fx;

#ifdef VAX11C
   temp = filename->data;
#else
   temp = filename;
#endif
  fx = *lun - 1;
  temp[*lfilename] = '\0';
  if ( strncmp( temp, "stdin", 5 ) == 0 )
    fp[fx] = stdin;
  else
    {
      if ( (fp[fx] = fopen(temp, "rb")) != NULL)
	{
	  *istat = 0;
#ifdef DSKIO_DEBUG
	  fprintf(stderr,
		  "!***INFO in DSKBOR, opened file %s for read with id %d\n", 
		  temp, fx);
#endif
	}
      else
	{
	  fprintf(stderr,
		  "!***ERROR in DSKBOR, could not open file %s\n", temp);
	  *istat = -2;
	}
    }
  return;
}


/*
 * *******************************
 * dskbcr - CLOSE file for reading
 * ******************************
 */
extern "C"
#ifdef VAX11C
void dskbcr(lun, istat)
  int *lun, *istat;
#else
void dskbcr_(int *lun, int *istat)
#endif
{
  int fx;
  fx = *lun - 1;
  if (fp[fx] == NULL)
    {
      fprintf(stderr,
	      "!***ERROR in DSKBCR, file not open, cannot close!\n");
      *istat = -1;
    }
  else if ( fp[fx] != stdin )
    {
      fclose(fp[fx]);
#ifdef DSKIO_BEBUG
      fprintf(stderr,
	      "!***INFO in DSKBOR, closed file for read with id %d\n", 
	      fx);
#endif
      fp[fx] = NULL;
      *istat = 0;
    }
  return;
}


/*
 * ****************************************
 * dskcr - read from previously opened file
 * ****************************************
 */
extern "C"
#ifdef VAX11C
void dskbr(lun, data, ldata, istat)
     char *data;
     int  *lun, *ldata, *istat;
#else
void dskbr_(int *lun, char *data, int *ldata, int *istat)
#endif
{
  int fx;
  fx = *lun - 1;
  if (fp[fx] != NULL)
    {
      if ( (*istat = fread(data, sizeof(*data), *ldata, fp[fx])) == 0)
	{
	  if (0 != feof(fp[fx]))
	    fprintf(stderr, "!***INFO in DSKBR, end-of-file!\n");
	  else
	    fprintf(stderr, "!***ERROR in DSKBR, reading error.\n");
	  *istat = -1;
	}
      else if (*istat == *ldata)
	*istat = 0;
      else if (0 != feof(fp[fx]))
	{
	  fprintf(stderr, "!***INFO in DSKBR, end-of-file!\n");
	}
      else
	fprintf(stderr, "!***ERROR in DSKBR, file short!\n");
    }
  else
    {
      fprintf(stderr, "!***ERROR in DSKBR, no file open!\n");
      *istat = -1;
    }
  return;
}

/*
 * ****************************************
 * dskbwr - wait for read to complete (dummy routine)
 * ****************************************
 */
extern "C"
#ifdef VAX11C
void dskbwr(lun, lflag)
     int  *lun, *lflag;
#else
void dskbwr_(int *lun, int *lflag)
#endif
{
  return;
}

/*
 * ***********************************
 * dskbow - Byte file Open for Writing
 * ***********************************
 */
extern "C"
#ifdef VAX11C
void dskbow(lun, filename, lfilename, size, istat)
     struct descr *filename;
     char *filename;
     int *lun, *lfilename, *size, *istat;
#else
void dskbow_(int *lun, char *filename, int *lfilename, int *size, int *istat)
#endif

{
//  FILE *fopen();
  char *temp;
  int   fx;

  fx = *lun - 1;
#ifdef VAX11C
  temp = filename->data;
#else
  temp = filename;
#endif

  temp[*lfilename] = '\0';
  if ( strncmp( temp, "stdout", 6 ) == 0 )
    fp[fx] = stdout;
  else
    {
      if ( (fp[fx] = fopen(temp, "wb")) != NULL)
	{
#ifdef DSKIO_BEBUG
	  fprintf(stderr,
		  "!***INFO in DSKBOR, opened file %s for write with id %d\n", 
		  temp, fx);
#endif
	  *istat = 0;
	}
      else
	{
	  fprintf(stderr,
		  "!***ERROR in DSKBOW, could not open file %s\n", temp);
	  *istat = -2;
	}
    }
  return;
}


/*
 * ************************************
 * dskbcw - Byte Close file for Writing
 * ************************************
 */
extern "C"
#ifdef VAX11C
void dskbcw(lun, istat)
  int  *lun, *istat;
#else
void dskbcw_(int *lun, int *istat)
#endif
{
  int fx;
  fx = *lun - 1;
  if (fp[fx] == NULL)
    {
      fprintf(stderr,
	      "!***ERROR in DSKBCW, file not open, cannot close!\n");
      *istat = -1;
    }
  else
    {
      fclose(fp[fx]);
#ifdef DSKIO_BEBUG
      fprintf(stderr,
	      "!***INFO in DSKBOR, closed file for write with id %d\n", 
	      fx);
#endif
      fp[fx] = NULL;
      *istat = 0;
    }
  return;
}


/*
 * ********************************************
 * dskbw - Byte Write to previously opened file
 * ********************************************
 */
extern "C"
#ifdef VAX11C
void dskbw(lun, data, ldata, istat)
     char *data;
     int  *lun, *ldata, *istat;
#else
void dskbw_(int *lun, char *data, int *ldata, int *istat)
#endif
{
  int  fx;
  fx = *lun - 1;
  if (fp[fx] != NULL)
    {
      if ( (*istat = fwrite(data, sizeof(char), *ldata, fp[fx])) == 0)
	{
	  fprintf(stderr,
		  "!***ERROR in DSKBW, writing error\n");
	  *istat = -1;
	}
      else if (*istat == *ldata)
	  *istat = 0;
      else
	fprintf(stderr, "!***ERROR in DSKBW, file short!\n");
    }
  else
    {
      fprintf(stderr, "!***ERROR in DSKBW, no file open!\n");
      *istat = -1;
    }
  return;
}


/*
 * ****************************************
 * dskbww - wait for write to complete (dummy routine)
 * ****************************************
 */
extern "C"
#ifdef VAX11C
void dskbww(lun, lflag)
     int  *lun, *lflag;
#else
void dskbww_(int *lun, int *lflag)
#endif
{
  return;
}


/*
 * ****************************************
 * dskbww - get integer file descriptor
 * ****************************************
 */

#ifndef VAX11C
extern "C" FILE* pFdskfile(int *lun)
{
    return (fp[*lun - 1]);
};

#endif
