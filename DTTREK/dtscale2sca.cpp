//
// Copyright (c) 2011 Rigaku Americas, Corp.
//
// dtscale2sca.cc   Initial author: J.W. Pflugrath           
//    Convert a d*TREK reflnlist file output by dtscaleaverage -anome to a scalepack .sca file with anom 
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
//+ToDo
//
//   Error messages need implementing
//
//+Include files

#include <iostream>
#include <iomanip>
#include "Dtrek.h"
#include "Cscan.h"
#include "Ccrystal.h"
#include "Cimage_header.h"

#if !defined(VC6)
using std::cout;
using std::endl;
using std::flush;
using std::ifstream;
#endif

#define DTREK_ERROR(code,msg) {vError(code,msg);return (code);}

//+Code begin

//+Definitions, constants, and initialization of static member variables

//+Public functions

void vError(const int nErrorNum=0, const Cstring& sError=NULL);

//

int main(int   argc,      // Command line argument count
         char *argv[])    // Pointers to command line args
     
{
  int          nStat;
  Cstring      sInFile       = "output.sca";
  Cstring      sTemplate;
  Cstring      sOutFile      = "output.ref";
  Cstring      sBasename;

  int          nDim0, nDim1;

  int          i, j, k;  // Loop counters
  int          nTemp;
  float        fTemp;
  Cstring      sTemp, sTempTemp;
  Ccrystal     oCrystal;
  float        fPixel;
  float        a6fCellArg[6] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
  int          nSpacegroupArg = 0;
  int          nFI_IntPlus, nFI_IntMinus, nFI_SigmaPlus, nFI_SigmaMinus;

  Cimage_header *poHeader = NULL;
  Creflnlist *poReflnlist = NULL;
  Cstring        sMissingOption
                 = "ERROR: dtscale2sca - missing option argument!";
  Cstring        sInvalidArgument
                 = "ERROR: dtscale2sca - invalid argument!";

  vDtrekSetModuleName("dtscale2sca");
  vPrintCopyrightInfo();

  // 
  // Usage: dtscale2sca dtscale.ref 

  // Parse command line arguments
  
  argc--; argv++;

  nStat = 1;

  if (0 < argc)
    {
      sInFile =   (const char*) argv[0];
      if ( ("-h" == sInFile) || ("-help" == sInFile) )
	{
	  DTREK_ERROR(0, sTemp);
	}
      if (!bFileExists(sInFile))
	vError(6, "ERROR input file not found.");
      argc--; argv++;
    }
  else
    {
      vError(6, sMissingOption);
    }
  for (i = 0; i < argc; i++) 
    {
      sTemp = (const char*) argv[i];
      cout << "Command line string: >>" << sTemp << "<<" << endl;
      
      if ( ("-h" == sTemp) || ("-help" == sTemp) )
	{
	  DTREK_ERROR(0, sTemp);
	}
      else if ("-cell" == sTemp)
        {
          i++;
          for (j=0; (j < 6) && (i < argc); j++, i++)
            {
              nStat = sscanf(argv[i], "%f", &fTemp);
              if (1 == nStat)
                {
                  a6fCellArg[j] = fTemp;
                }
              else
                DTREK_ERROR(6, sInvalidArgument);
            }
          i--;
          if (6 != j)
            {
              DTREK_ERROR(5, sMissingOption);
            }
        }
      else if ("-spacegroup" == sTemp)
        {
          i++;
          if (i < argc)
            {
              nStat = sscanf(argv[i], "%d", &nTemp);
              if (1 == nStat)
		{
		  nSpacegroupArg = nTemp;
		}
              else
                DTREK_ERROR(6, sInvalidArgument);
            }
          else
            {
              DTREK_ERROR(5, sMissingOption);
            }
        }
      else if ("-help" == sTemp)
	{
	  vError(0, "");
	}
      else if ("-help" == sTemp)
	{
	  vError(0, "");
	}
      else
	{
	  vError(0,"" );
	}
    }

  ////////////////////////////////////////////////////////////////////////
  // Done with command line arguments so now go do the deed
  ////////////////////////////////////////////////////////////////////////

  // First list inputs

  sOutFile = sFileGetBasename(sInFile);
  // Remove extension and replace with .sca
  Cstring sExtension = sFileGetExtension(sOutFile);
  sBasename = sOutFile.before(sExtension); 
  // sBasename should end with a .
  sOutFile  = sBasename + "sca";  
  cout << "Input d*TREK dtscale file: " << sInFile
       << "\nOutput scalepack     file: " << sOutFile
       << "\n\n";

  // Some error checking

  nStat       = 0;

  // Try to open the file
  Cstring sMessage;

  if (!bFileExists(sTransSymbol(sInFile)))
    {
      nStat = -2;
      sMessage = "ERROR opening file: ";
    }
  else
    {
      // Read the d*TREK  reflnlist to hold the reflection information

      poReflnlist = new Creflnlist(sInFile);

      if (!poReflnlist->bIsAvailable())
	{
	  nStat = -2;
	  sMessage = "ERROR reading file: ";
	}
      else if (0 >= poReflnlist->nGetNumReflns())
	{
	  nStat = -2;
	  sMessage = "ERROR no reflections in the file: ";
	}
      else 
	{
	  // Check that the proper fields exist (I+, I-, etc) 
	  nFI_IntPlus    = -1;
	  nFI_IntMinus   = -1;
	  nFI_SigmaPlus  = -1;
	  nFI_SigmaMinus = -1;
	  nFI_IntPlus = poReflnlist->nGetFieldIndex(Creflnlist::ms_sfIntensityPlus);
	  nFI_SigmaPlus = poReflnlist->nGetFieldIndex(Creflnlist::ms_sfSigmaIPlus);
	  nFI_IntMinus = poReflnlist->nGetFieldIndex(Creflnlist::ms_sfIntensityMinus);
	  nFI_SigmaMinus = poReflnlist->nGetFieldIndex(Creflnlist::ms_sfSigmaIMinus);
	  if (    (0 < nFI_IntPlus)
	       && (0 < nFI_SigmaPlus)
               && (0 < nFI_IntMinus)
               && (0 < nFI_SigmaMinus)
             )
	    {
	      nStat = 0;
	    }
	  else
	    {
	      nStat = -3;
	      sMessage = "ERROR missing necessary columns in the reflection file: ";
	    }
	}
    }
  if (0 != nStat)
    {
      if (NULL != poReflnlist)
	{
	  delete poReflnlist;
	  poReflnlist = NULL;
	}
      cout << sMessage << sInFile << endl;
      exit (nStat);
    }

  // Reflnlist file exists, is read, has reflns and the necessary columns

  // Get the unit cell and spacegroup to write out
  // Only put the cell stuff in if we have a unit cell AND it was NOT overridden by the command line
  
  poHeader = poReflnlist->poGetHeader();
  if (0.0 >= (a6fCellArg[0] * a6fCellArg[1] * a6fCellArg[2] * a6fCellArg[3] * a6fCellArg[5] *a6fCellArg[5]))
    {
      // Get unit cell from .ref file
      nTemp = poHeader->nGetValue(Ccrystal::ms_sCrystalXUnitCell, 6, a6fCellArg);
      if (0 != nTemp)
	{
	  cout << "WARNING, no unit cell info in the input file: " << sInFile << endl;
	}
    }
  if (0 >= nSpacegroupArg)
    {
      nTemp = poHeader->nGetValue(Ccrystal::ms_sCrystalXSpacegroup, &nSpacegroupArg);
      if (0 != nTemp)
	{
	  cout << "WARNING, no spacegroup info in the input file: " << sInFile << endl;
	}
    }

  oCrystal.vSetCell(a6fCellArg);
  oCrystal.m_poSpacegroup->vSet(nSpacegroupArg);

  // Need to figure out if refln are centric or acentric for use later
  // May be easier to use the Cspacegroup::nReduceHKL(poRefln);
  
  // Convert spacegroup to a name in lowercase and without parentheses
  Cstring sCellName;
  sCellName = oCrystal.m_poSpacegroup->sGetName();
  sCellName.downcase();
  // Write out the unit cell and spacegroup

  int nCount = 0;
  int nErrorCount = 0;
  int nAbsentCount = 0;
  int nReduce;

  // Try to prevent overwrite of existing file.  WARNING: renamed version seems bogus
  nFileAppendVersion(sOutFile, FALSE);

  FILE   *pFTextOutFile = NULL;
  pFTextOutFile = fopen(sOutFile.string(), "w+t");
  if (!pFTextOutFile)
    {
      cout << "ERROR opening output file: " << sOutFile << endl;
      nStat = -4;
    }
  if (0 == nStat)
    {
      // Write out first line: 1 
      (void) fprintf(pFTextOutFile, "    1\n");
      // Write out the 2nd line: -985
      (void) fprintf(pFTextOutFile, " -985\n");
      // Write out the 3rd line:
      (void) fprintf(pFTextOutFile,"%10.3f%10.3f%10.3f%10.3f%10.3f%10.3f",
		     a6fCellArg[0], a6fCellArg[1], a6fCellArg[2],
		     a6fCellArg[3], a6fCellArg[4], a6fCellArg[5]);
      (void) fprintf(pFTextOutFile, " %s\n", sCellName.string());

      float dInt, dSig, dIntP, dIntM, dSigP, dSigM;

      Crefln oRefln(poReflnlist);   // Create a reflection object
      Crefln *poRefln;
      float fMax = -999999.99;
      float fMin = +999999.99;

      // Figure out how to scale the intensities and sigmas so there are no overflow of characters

      for (i = 0; i < poReflnlist->nGetNumReflns(); i++)
	{
	  poRefln = poReflnlist->poGetRefln(i);
	  dInt  = poRefln->fGetIntensity();
	  dSig  = poRefln->fGetSigmaI();
	  dIntP = poRefln->fGetField(poReflnlist->m_nFI_fIntensityPlus);
	  dSigP = poRefln->fGetField(poReflnlist->m_nFI_fSigmaIPlus);
	  dIntM = poRefln->fGetField(poReflnlist->m_nFI_fIntensityMinus);
	  dSigM = poRefln->fGetField(poReflnlist->m_nFI_fSigmaIMinus);
	  fMax = max(fMax, dInt);
	  fMax = max(fMax, dSig);
	  fMax = max(fMax, dIntP);
	  fMax = max(fMax, dSigP);
	  fMax = max(fMax, dIntM);
	  fMax = max(fMax, dSigM);
	  fMin = min(fMin, dInt);
	  fMin = min(fMin, dSig);
	  fMin = min(fMin, dIntP);
	  fMin = min(fMin, dSigP);
	  fMin = min(fMin, dIntM);
	  fMin = min(fMin, dSigM);
	}
      double dScale = 1.0;
      //         123456.1
      if (fMax > 999999.9)
	dScale = 999999.9 / fMax;

      if ( (dScale * fMin) < -99999.9)
	dScale = -99999.9 / (dScale * fMin);

      if (1.0 != dScale)
	cout << "INFO: refln intensities and sigmas scaled by " << dScale
             << "\n     in order to prevent overflow in the character field(s)." << endl;

      // Loop through the reflns and output in 

      nStat = 0;
      for (i = 0; (i < poReflnlist->nGetNumReflns()) && (0 == nStat); i++)
	{
	  nCount++;
	  poRefln = poReflnlist->poGetRefln(i);
	  dInt  = poRefln->fGetIntensity();
	  dSig  = poRefln->fGetSigmaI();
	  dIntP = poRefln->fGetField(poReflnlist->m_nFI_fIntensityPlus);
	  dSigP = poRefln->fGetField(poReflnlist->m_nFI_fSigmaIPlus);
	  dIntM = poRefln->fGetField(poReflnlist->m_nFI_fIntensityMinus);
	  dSigM = poRefln->fGetField(poReflnlist->m_nFI_fSigmaIMinus);

	  // Scale things before looking at negatives

	  dInt  *= dScale;
	  dIntM *= dScale;
	  dIntP *= dScale;
	  dSig  *= dScale;
	  dSigM *= dScale;
	  dSigP *= dScale;

	  // .sca files have an intensity of 0 if the corresponding sigma is negative
	  if (0.0 > dSig)  
	    {
	      dInt  = 0.0;
	      dSig  = -1.0;
	    }
	  if (0.0 > dSigP)
	    {
	      dIntP  = 0.0;
	      dSigP  = -1.0;
	    }
	  if (0.0 > dSigM) 
	    {
	      dIntM  = 0.0;
	      dSigM  = -1.0;
	    }

	  // Is the Miller index centric or acentric? 
	  nReduce = oCrystal.m_poSpacegroup->nReduceHKL(poRefln, NULL, NULL);
	  if (0 < nReduce)
	    {
	      // Centric
	      (void) fprintf(pFTextOutFile,"%4d%4d%4d%8.1f%8.1f\n",
			     poRefln->nGetH(), poRefln->nGetK(), poRefln->nGetL(),
			     dInt, dSig);
	    }
	  else if (0 == nReduce)
	    {
	      // Acentric
	      (void) fprintf(pFTextOutFile,"%4d%4d%4d%8.1f%8.1f%8.1f%8.1f\n",
			     poRefln->nGetH(), poRefln->nGetK(), poRefln->nGetL(),
			     dIntP, dSigP, dIntM, dSigM);
	    }
	  else if (-1 == nReduce)
	    {
	      cout << "WARNING Not output: systematically absent refln with HKL: " <<
		poRefln->nGetH() << ' ' << poRefln->nGetK() << ' ' << poRefln->nGetL() << endl;
	      nAbsentCount++;
	      nCount--; 
	    }
	  else
	    {
	      // Error
	      cout << "ERROR with refln number " << i << " with HKL: " <<
		poRefln->nGetH() << ' ' << poRefln->nGetK() << ' ' << poRefln->nGetL() << endl;
	      nErrorCount++;
	      nCount--; 
	    }
	}
    }
  // Close the output file
  if (pFTextOutFile) fclose(pFTextOutFile);

  //cout << "\nnStat = " << nStat << "\n";

  if (0 != nStat)
    {
      cout << "ERROR writing " << sOutFile << endl;
    }
  if (0 < nCount)
    {
      cout << "\nINFO, there are " << nCount << " reflns in the output file " << sOutFile << endl;
    }
  else if (0 >= nCount)
    {
      cout << "\nWARNING, no reflnlist written since no valid reflns!\n\n";
    }
  if (0 < nAbsentCount)
    {
      cout << "\nWARNING, there were " << nAbsentCount << " systematically absent reflns removed.\n\n";
    }
  if (0 < nErrorCount)
    {
      cout << "\nWARNING, there were " << nErrorCount << " lines that could not be parsed.\n\n";
    }
  if (NULL != poReflnlist)
    {
      delete poReflnlist;
      poReflnlist = NULL;
    }
  return (nStat);
}

void vError(const int nErrorNum, const Cstring& sMessage)
{
  if (0 != nErrorNum)
    {
      cout << sMessage << endl;
    }
  cout << "\ndtscale2sca - Convert a d*TREK dtscale.ref file to an HKL-2000 *.sca file.\n\n"
       << "Usage:\n"
       << "   dtscale2sca dtscale.ref [ options ...]\n\n"
       << "Command line options: Description:\n\n"

       << " dtscale.ref          Name of a d*TREK refln file to convert.  Please use the\n"
          "                      dtscale.ref file created by dtscaleaverage with the -anom option.\n\n" 
          " -spacegroup nSpace   Spacegroup to use.  Default from the dtscale.ref file or 0.\n\n"
          " -cell a b c alp bet gam\n"
          "                      Specify the unit cell.  Default from the dtscale.ref input file.\n\n"
          "Output file will be input file basename + .sca.\n\n"
          "Examples:\n\n"
          "   dtscale2sca dtscale.ref\n"
          "   dtscale2sca dtscale.ref -cell 78.8 78.8 37.9 90 90 90 -spacegroup 96\n"
       << endl;
  exit (nErrorNum);
}
