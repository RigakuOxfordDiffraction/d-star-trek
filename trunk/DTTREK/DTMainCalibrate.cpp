//
// Copyright (c) 1998-2006 Rigaku
//
//
// RESTRICTED RIGHTS NOTICE SHORT FORM (JUNE 1987)
//
// Use, reproduction, or disclosure is subject to restrictions set
// forth in Contract No. W-31-109-ENG-38 and Contract No.
// 943072401 with the University of Chicago, Operator of
// Argonne National Laboratory.
//
// DTMainCalibrate.cpp   Initial author: RB           08-Sep-2006
//
// Transferred from dtcalibrate.cc  Initial author: J.W. Pflugrath  02-June-1998
// This implements detector calibration.  Built from Marty Stanton's
// calibrate.

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
 
#include <stdio.h>        // for sprintf, until we figure out sstream
#include "Dtrek.h"
#include "dtrekvec.h"
#include "Cstring.h"
#include "Creflnlist.h"
#include "Cimage_header.h"

#include "DTMainCalibrate.h"

#include "Ccalibrate.h"

using std::cin;
using std::cout;
using std::endl;
using std::flush;


int CDTMainCalibrate::nExecute(unsigned int argc,      // Command line argument count
                               char        *argv[])    // Pointers to command line args
{
  int            i;
  int            nStat;
  Cstring        sTemp;
  Cstring        sCalibFiles;
  Cstring        sCalibOptions;
  Cstring        sProgram;
  Cimage_header *poHeader;
  Ccalibrate    *poCalib;
  bool           bNeedNewHeader = FALSE;

  // 
#ifndef  CORRECT_ONLY
#ifdef DISTOR_ONLY
  sProgram = "dtdistor";
#else
  sProgram = "dtcalibrate";
#endif
#else
  sProgram = "dtcorrect";
#endif
  vDtrekSetModuleName(sProgram);
  vPrintCopyrightInfo();

  //cout << D_K_DTREKVersion << endl;

  cout << "\nReference:\n"
       << "  Stanton, M., Phillips, W.C., Li, Y., and Kalata, K. (1992)\n"
       << "     'Correcting spatial distortions and nonuniform response in\n"
       << "     area detectors', J. Appl. Cryst. 25, 549-558.\n" << endl;
 
  // Copy command line files to sCalibFiles

  for (i = 1; (i < argc) && (i < 2); i++)
    {
      sCalibFiles = sCalibFiles + ' ' + (const char*) argv[i];
    }

  // Copy command line options to sCalibOptions

  for (i = 2; i < argc; i++)
    {
      sCalibOptions = sCalibOptions + ' ' + (const char*) argv[i];
    }

  bNeedNewHeader =    sCalibOptions.contains("-godistor") 
                   || sCalibOptions.contains("-gononunf")
                   || sCalibOptions.contains("-gomerge");

  // Get header file to start with
  
  argc--;
  argv++;

  nStat = 1;

  if (1 > argc)
    {
      // No header file argument and thats an error!
      DTREK_ERROR(1, "ERROR - no header filename!");
    }

  sTemp = (const char*) argv[0];
  if ( ("-help" == sTemp) || ("-h" == sTemp) ) 
    {
      DTREK_ERROR(0, "");
    }
  else if ("-noheader" == sTemp)
    {
      poHeader     = new Cimage_header();   // Create empty image header
    }
  else
    {
      // Read the name of a file and read the header only from that image file

      poHeader     = new Cimage_header(sTemp);   // Create (read) image header
    }

  if (!poHeader->bIsAvailable())
    {
      DTREK_ERROR(2, "Header not available!");
    }

  argc--;
  argv++;

  sTemp = (const char*) argv[0];

  // Create calibrate object

  poCalib   = new Ccalibrate(*poHeader);

  // Parse rest of command line arguments and do calibration

  nStat = poCalib->nDoCalibrate(sCalibOptions);

  if (sCalibOptions.contains("-prompt"))
    {
      cout << sProgram << " - PROMPT Write new output header? (y,n) [y]: ";
      getline(cin, sTemp);
      if ("" == sTemp) sTemp = 'y';
      if (!cin)
	nStat = 0;
      else if ( ('n' == sTemp.GetAt(0)) || ('N' == sTemp.GetAt(0)) )
	nStat = 1;
      else
	nStat = 0;
      cout << endl;
    }

  if (   sCalibOptions.contains("-h")
      || sCalibOptions.contains("-help"))
    {
      DTREK_ERROR(0, "");
    }

  if ( (0 == nStat) && (bNeedNewHeader) )
    {
      // Write results to output

      poCalib->nUpdateHeader(poHeader);
      poHeader->nReplaceValue(Ccalibrate::ms_sDtcalibFiles, sCalibFiles);
      poHeader->nReplaceValue(Ccalibrate::ms_sDtcalibOptions, 
			      sCalibOptions);
      nStat = poHeader->nWrite(poCalib->sGetOutHeader());
      if (0 == nStat)
	{
	  (void) poCalib->nNotifyDisplay("", sGetCWD() 
					  + poCalib->sGetOutHeader());
	  cout << sProgram << " - INFO wrote header file ";
	}
      else
	{
	  cout << sProgram << " - ERROR writing header file ";
	}
      cout << poCalib->sGetOutHeader() << endl;
    }
  else
    {
      cout << sProgram << " - WARNING did not write new header file.\n";
    }
  delete poHeader;
  delete poCalib;
  return (nStat);
}

void CDTMainCalibrate::vError(const int nErrorNum, const Cstring& sMessage)
{
  if (0 != nErrorNum)
    {
      cout << sMessage << endl;
    }
#ifndef CORRECT_ONLY
#ifdef DISTOR_ONLY
  cout << "\ndtdistor - Usage:\n"
       << "dtdistor header_file [ options ...]\n\n"
#else
  cout << "\ndtcalibrate - Usage:\n"
       << "dtcalibrate header_file [ options ...]\n\n"
#endif
#else
  cout << "\ndtcorrect - Usage:\n"
       << "dtcorrect -noheader [ options ...]\n\n"
#endif
       << "Command line options: Description:\n\n"
#ifndef CORRECT_ONLY
       << " header_file      Name of a file with a d*TREK header used to\n"
       << "                  initialize the calibration procedure.\n\n"
       << " -peakmin Int     Minimum integrated peak value for spot picking in the mask.\n"
       << " -peakcent fPx0 fPx1  First or center spot position in pixels the mask.\n"
       << " -dist  d         The sample to detector distance in mm.\n"
       << " -dump   file     Turn on dump mode and set output dump filename.\n"
       << " -mask   file     Sets input mask filename.\n"
#ifndef DISTOR_ONLY
       << " -flood  file     Sets input flood filename.\n"
       << " -badpix file     Sets input bad pixel filename.\n"
       << " -refer  file     Sets input reference filename.\n"
#endif
       << " -godistor        Run the 'Go DISTOR' command automatically.\n"
#ifndef DISTOR_ONLY
       << " -gononunf        Run the 'Go NONUNF' command automatically.\n"
       << " -dadd Dark_offset\n"
       << "                  When a dark image is subtracted, this offset or pedestal\n"
       << "                      keeps pixel values from becoming negative. Default: 0.\n"
       << " -dscale Dark_scale\n"
       << "                  When a dark image is subtracted, dark pixels are scaled\n"
       << "                      by Dark_scale before subtraction.  Default: 1.\n"
#endif
#endif
       << " -beamx fpx       The fast pixel coordinate of the beam position.\n"
       << " -beamy spx       The slow pixel coordinate of the beam position.\n"
       << " -distor file     Sets input/output spatial distortion files basename.\n"
#ifndef DISTOR_ONLY
       << " -transform file  Sets input/output transform filename. Default: TRANSFORM.\n"

       << " -nonunf file     Sets input/output nonuniformity filename.\n"
       << " -dark   file     Sets input dark filename. Default: DARK.\n"
#ifndef CORRECT_ONLY
       << " -tbadpix file    Sets output bad pixel filename.\n"
       << " -tlimit  nOrig0 nOrig1 nExt0 nExt1\n"
       << "                  Sets the limits of the transformed output image that appears\n"
       << "                  in the file written to disk.  To work correctly, it must be\n"
       << "                  used with -gotransform and -correct dt -gocorrect ... options.\n"
       << " -gotransform     Create a transform file for later use.\n"
#else
       << " -correct dns     Sets which correction to perform:\n"
       << "                          d: +dark current\n"
       << "                          n: +nonuniformity\n"
       << "                          s: +spatial distortion\n"
       << " -correct dt | t  Sets which transformation to perform:\n"
       << "                         dt: for dark current and transform\n"
       << "                          t: transform only\n"
#endif
       << " -gocorrect filein fileout\n"
       << "                  Input and output filenames to transform.  If\n"
       << "                      the filein contains any ? characters, then\n"
       << "                      it is treated as a scan template.\n\n"
       << " -seq Start End   Starting and ending sequence numbers if a scan\n"
       << "                      template is used with -gocorrect.  Default is\n"
       << "                      all available sequence numbers.\n"
       << " -verbose nVerlev Verbosity level of output.  Higher values give\n"
       << "                      more output.  Default: 1.\n"
#endif
       << " -h               Print some command line help text.\n\n";
#ifdef DISTOR_ONLY
  cout << "Example:\n"
       << " dtdistor dtcal01.head -mask grid.mask -godistor\n\n";
#else
#ifndef CORRECT_ONLY
  cout << "Examples:\n"
       << " dtcalibrate dtcalibrate.head -godistor -gononunf\n"
       << " dtcalibrate dtcalibrate.head -beamx 576 -beamy 582 -godistor -gononunf\n\n";
#else
  cout << "Example:\n"
       << " dtcorrect -noheader -dark dark.img -correct dt -gocorrect 'img.???' 'img.???c'\n\n";
#endif
#endif

}



