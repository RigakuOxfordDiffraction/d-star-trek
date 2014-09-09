//
// Copyright (c) 1996 Molecular Structure Corporation
//
// RESTRICTED RIGHTS NOTICE SHORT FORM (JUNE 1987)
//
// Use, reproduction, or disclosure is subject to restrictions set
// forth in Contract No. W-31-109-ENG-38 and Contract No.
// 943072401 with the University of Chicago, Operator of
// Argonne National Laboratory.
//
// dtscalemerge.cc     Initial author: J.W. Pflugrath           14-Jan-1996
//    This reads a D*TREK style reflection files, then scales and merges them.
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
//    dtscalemerge reads from the command line the filename of a
//    d*TREK image file to get crystal info from and the filename
//    of a reflection file to compute and apply scale factors for.
//    It reads in the input file,  and writes out the output file.
//    As a safety precaution, the output file will not be created or
//    overwritten if it already exists.
//
//    Example:  % dtscalemerge  image.head inputfile dtscalemerge.ref
//
//+ToDo
//
//   Error messages need implementing
//
//+Include files

#include <iostream.h>
#include <iomanip.h>
#include "Cstring.h"
#include "Crefln.h"
#include "Creflnlist.h"
#include "Cimage_header.h"
#include "Ccrystal.h"
#include "Cscalemerge.h"

//+Function prototypes

#ifdef SSI_PC
#define main   dtscalemerge_main
#define vError dtscalemerge_vError
#endif

void vError(const int nErrorNum=0, const Cstring& sError=NULL);

int main(int   argc,      // Command line argument count
         char *argv[])    // Pointers to command line args
{
  // Now parse and apply any output file options

  Cscalemerge *poScalemerge;
  Creflnlist  *poReflnlistIn;     // Pointers to input file
  Creflnlist   oReflnlistOut;
  Cstring      sInFile;
  Cstring      sOutScaledFile     = "?";  // dtscale.ref
  Cstring      sOutAverageFile    = "?";  // dtunavg.ref
  Cstring      sBatchFixed  = "";
  Cstring      sTemp;
  Cstring      sOut = "dtscalemerge.head";
  float        fScaleFixed  = 1.0;
  float        fBValueFixed = 0.0;
  float        fErrorAdd    = 0.0;
  float        fErrorMul    = 1.0;
  float        fRejectCrit  = (float)10.0E10;
  float        fSigma       = 0.0;
  int          nCycles      = 30;
  bool         bSigmaOverEstimate = FALSE;
  int          nOutputAnomFlag = 0;
  int          nScaleAnomFlag  = 0;
  int          nFixBFlag    = 0;
  int          nVerboseLevel= 3;
  int          nTexsan      = 0;
  int          nNoHeader    = 0;
  int          nNoUnavgHeader  = 0;
  int          nCountOverlap= 0;
  float        a2fResolution[2] = {999999.00, 0.00001f};

  Cimage_header *poHeader;
  poHeader  = NULL;

  int   j;
  int   nTemp;
  float fTemp;
  int         nStat = 0;

  vDtrekSetModuleName("dtscalemerge");
  cout << "\ndtscalemerge: Copyright (c) 1996 Molecular Structure"
       << " Corporation\n";
  cout << D_K_DTREKVersion << endl;

  // Copy command line to output log

  cout << "Command line:\n " << sGetCommandLine(argc, argv, 71) << endl << flush;
  // Parse command line arguments, if any (there had better be some!)

  argc--; argv++; // Skip program name

  int nArg = 0;
  nStat    = 0;

  if (nArg < argc)
    {
      // Read image header for some information

      poHeader = new Cimage_header(Cstring((const char*)argv[nArg]));
      if (!poHeader->bIsAvailable())
	{
	  delete poHeader;
	  vError(1, "Header not available!\n");
	}
    }

  nArg++;

  // By position ... this must be input reflnlist filename

  if (nArg < argc)
    {
      // Get input reflection file name, but do not read it in yet

      sInFile = argv[nArg];
    }
  else
    {
      vError(1,"Not enough command line arguments!\n");  // This calls exit
    }
  nArg++;

  while ( (nArg < argc) && ('-' == (char) *argv[nArg]) && (0 == nStat) )
    {
      // Remove leading '-'!
      sTemp = Cstring((argv[nArg]+1));
      if ("fix" == sTemp)
	{
	  // Batch name to fix is found, read next word to see what it is

	  nArg++;
	  if (nArg < argc)
	    {
	      sBatchFixed = (const char*)argv[nArg];
	    }
	  else
	    {
	      vError(1, "Missing batch name to fix!\n");  // This calls exit
	    }
	}
      else if ("add" == sTemp)
	{
	  nArg++;
	  if (nArg < argc)
	    {
	      nStat = -1 + sscanf((const char*)argv[nArg], "%f", &fErrorAdd);
	      if (0 != nStat)
		{
		  vError(1,"Bogus factor for -add!\n");
		}
	    }
	  else
	    {
	      vError(1,"Missing factor for -add!\n");
	    }
	}
      else if("out" == sTemp)
      {
         nArg++;
         if(nArg < argc)
            sOut = (const char *)argv[nArg];
      }
	  else if ("sigmaoverestimate" == sTemp) {
	  bSigmaOverEstimate = TRUE;

      } else if ("mul" == sTemp)
	{
	  nArg++;
	  if (nArg < argc)
	    {
	      nStat = -1 + sscanf((const char*)argv[nArg], "%f", &fErrorMul);
	      if (0 != nStat)
		{
		  vError(1,"Bogus factor for -mul!\n");
		}
	    }
	  else
	    {
	      vError(1,"Missing factor for -mul!\n");
	    }
	}
      else if ("reject" == sTemp)
	{
	  nArg++;
	  if (nArg < argc)
	    {
	      nStat = -1 + sscanf((const char*)argv[nArg], "%f", &fRejectCrit);
	      if (0 != nStat)
		{
		  vError(1,"Bogus factor for -reject!\n");
		}
	    }
	  else
	    {
	      vError(1,"Missing factor for -reject!\n");
	    }
	}
      else if ("cycles" == sTemp)
	{
	  nArg++;
	  if (nArg < argc)
	    {
	      nStat = -1 + sscanf((const char*)argv[nArg], "%d", &nCycles);
	      if (0 != nStat)
		{
		  vError(1,"Bogus factor for -cyc!\n");
		}
	    }
	  else
	    {
	      vError(1,"Missing factor for -cyc!\n");
	    }
	}
      else if ("sigma" == sTemp)
	{
	  nArg++;
	  if (nArg < argc)
	    {
	      nStat = -1 + sscanf((const char*)argv[nArg], "%f", &fSigma);
	      if (0 != nStat)
		{
		  vError(1,"Bogus factor for -sigma!\n");
		}
	    }
	  else
	    {
	      vError(1,"Missing factor for -sigma!\n");
	    }
	}
      else if ("scale" == sTemp)
	{
	  nArg++;
	  if (nArg < argc)
	    {
	      nStat = -1 + sscanf((const char*)argv[nArg], "%f", &fScaleFixed);
	      if (0 != nStat)
		{
		  vError(1,"Bogus scale factor for fixed batch!\n");
		}
	    }
	  else
	    {
	      vError(1,"Missing scale factor for fixed batch!\n");
	    }
	}
      else if ("bfac" == sTemp)
	{
	  nArg++;
	  if (nArg < argc)
	    {
	      nStat = -1 + sscanf((const char*)argv[nArg], "%f", &fBValueFixed);
	      if (0 != nStat)
		{
		  vError(1,"Bogus temperature factor for fixed batch\n");
		}
	    }
	  else
	    {
	      vError(1,"Missing temperature factor for fixed batch!\n");
	    }
	}
      else if ("reso" == sTemp)
	{
	  nArg++;
	  for (j=0; (j < 2) && (nArg < argc); j++, nArg++)
	    {
	      nTemp = sscanf(argv[nArg], "%f", &fTemp);
	      if (1 == nTemp)
		a2fResolution[j] = fTemp;
	      else
		vError(6, "Invalid resolution value for -reso option!\n");
	    }
	  nArg--;
	  if (2 != j)
	    {
	      vError(5, "Missing resolution value for -reso option!\n");
	    }
	}
      else if ( ("ref" == sTemp) || ("unavg" == sTemp) )
	{
	  nArg++;
	  if (nArg < argc)
	    {
	      sOutScaledFile = (const char*)argv[nArg];
	    }
	  else
	    {
	      vError(1,"Missing filename for -ref option!\n");
	    }
	}
      else if ("verbose" == sTemp)
	{
	  nArg++;
	  if (nArg < argc)
	    {
	      nStat = -1 + sscanf((const char*)argv[nArg], "%d", &nVerboseLevel);
	      if (0 != nStat)
		{
		  vError(1,"Bogus factor for -verbose!\n");
		}
	    }
	  else
	    {
	      vError(1,"Missing factor for -verbose!\n");
	    }
	}
      else if ("anom" == sTemp)
	{
	  nOutputAnomFlag = 1;
	}
      else if ("scaleanom" == sTemp)
	{
	  nScaleAnomFlag = 1;
	}
      else if ("fixB" == sTemp)
	{
	  nFixBFlag = 1;
	}
      else if ("fixb" == sTemp)
	{
	  nFixBFlag = 1;
	}
      else if ("countoverlap" == sTemp)
	{
	  nCountOverlap = 1;
	}
      else if ("texsan" == sTemp)
	{
	  nTexsan = 1;
	}
      else if ("texsan2" == sTemp)
	{
	  nTexsan = 2;
	}
      else if ("noheader" == sTemp)
	{
	  nNoHeader = 1;
	}
      else if ("nounavgheader" == sTemp)
	{
	  nNoUnavgHeader = 1;
	}
      else
	{
	  sTemp = "Unknown argument: " + sTemp + "\n";
//	  vError(1, (const char *) sTemp);
	  vError(1, sTemp.string());
	}
      nArg++;
    }

  if (nArg < argc)
    {
      // By position ... this must be output reflnlist filename

      sOutAverageFile = argv[nArg];
      nArg++;
    }
  if (nArg < argc)
    {
      vError(1,"Too many command line arguments!\n");  // This calls exit
    }

  cout << "Reading " << sInFile << " ..." << endl << flush;
  poReflnlistIn = new Creflnlist(sInFile);

  if (   !poReflnlistIn->bIsAvailable()
      || (1 > poReflnlistIn->nGetNumReflns()) )
    {
      sTemp = "Problem with input reflection list: " + sInFile;
      delete poHeader;
      delete poReflnlistIn;
      vError(2, sTemp);
    }

  // Have reflection list, now do something with it.

  poScalemerge = new Cscalemerge(poReflnlistIn, poHeader);

  if (!poScalemerge->bIsAvailable())
    {
      delete poHeader;
      delete poReflnlistIn;
      delete poScalemerge;
      vError(3, "Problem with input files!\n");
    }

  poScalemerge->vSetSigmaOverEstimate(bSigmaOverEstimate);
  poScalemerge->vSetMaxCycles(nCycles);
  poScalemerge->vSetCountOverlap(nCountOverlap);
  poScalemerge->vSetScale(fScaleFixed);
  poScalemerge->vSetTemp(fBValueFixed);
  poScalemerge->vSetError(fErrorMul, fErrorAdd);
  poScalemerge->vSetRejCrit(fRejectCrit);
  poScalemerge->vSetScaleAnomFlag(nScaleAnomFlag);
  if ( (0 != nScaleAnomFlag) && (1 == nOutputAnomFlag) )
    {
      cout << "\nWARNING!  Cannot have both -scaleanom and -anom options,"
	   << "\n          -anom option ignored!\n" << endl;
      nOutputAnomFlag = 0;
    }
  poScalemerge->vSetOutputAnomFlag(nOutputAnomFlag);
  poScalemerge->vSetFixBFlag(nFixBFlag);
  poScalemerge->vSetResolution(a2fResolution[0], a2fResolution[1]);
  poScalemerge->vSetVerbose(nVerboseLevel);
  if ("" != sBatchFixed) poScalemerge->vSetBatch(sBatchFixed);

  cout << "Inp ref file: " << sInFile << '\n' << flush;
  (void) poScalemerge->nList();

  nStat = poScalemerge->nScaleSetup();

  if (0 != nStat)
    {
      cerr << "ERROR in Cscalemerge::nScaleSetup!\n" << flush;
      delete poHeader;
      delete poReflnlistIn;
      delete poScalemerge;
      return (nStat);
    }

  // Deselect some reflns via a sigma cutoff
  // The check is before any modification of sigma via fErrorMul and fErrorAdd!
  // Do this after nScaleSetup because any new Crefln fields
  // have now all been added.

  char a255cTemp[255];
  sprintf(a255cTemp, "-fIntensity/fSigmaI<%.3f", fSigma);
  cout << "\nExcluding reflections from the scale factor calculation with"
       << "\nselection string: " << a255cTemp
       << endl << flush;

//  nStat = poReflnlistIn->nSelect((Cstring)a255cTemp);
  nStat = poScalemerge->m_poReflnlist->nSelect((Cstring)a255cTemp);

  if (0 > nStat)
    {
      cerr << "ERROR in selection mechanism!\n" << flush;
      delete poHeader;
      delete poReflnlistIn;
      delete poScalemerge;
      return (nStat);
    }
  cout << "Number of reflns which match above selection: " << nStat
       << "\nHowever, all reflections are used in the summary tables below.\n"
       << endl << flush;

  nStat = 0;

  if (0 == nStat)
    nStat = poScalemerge->nScale1(nCycles, &oReflnlistOut);// This does scaling!

  if ( (0 == nStat) && ("?" != sOutScaledFile) )
    {
      poScalemerge->vApplyWeights();
      if (0 != nTexsan)
	{
	  // Write out the unaveraged reflection list for texsan purposes:
	  // 1. Do not write out the header.
	  // 2. Do not write out reflns with sigmaI <= 0.
	  // -texsan:
	  // 3. Write out only h, k, l, I, sigmaI fields (the first 5 fields).
//+ 2-Mar-1999
	  // -texsan2:
	  // 4. Write out I/Lp, and, if available,  the absorption correction
	  //               factors 
//- 2-Mar-1999
	  int i;
	  int nFields;
	  char *pcSelect = NULL;
	  nFields  = poScalemerge->m_poReflnlist->nGetNumFields();
	  pcSelect = new char [nFields + 1];

	  // "Deselect" all fields

	  for (i = 0; i < nFields; i++)
	    pcSelect[i] = poScalemerge->m_poReflnlist->m_pcSelect[i] + 1;

	  // Keep the global selection the same

	  pcSelect[nFields] = poScalemerge->m_poReflnlist->m_pcSelect[nFields];

	  // Select only the 5 fields we want to keep

	  (void) poScalemerge->m_poReflnlist->nSelectField(
                                eReflnField_int_type,
			        poScalemerge->m_poReflnlist->m_nFI_nH,
			        char(0),
	                        pcSelect);
	  (void) poScalemerge->m_poReflnlist->nSelectField(
                                eReflnField_int_type,
				poScalemerge->m_poReflnlist->m_nFI_nK,
                                char(0),
				pcSelect);
	  (void) poScalemerge->m_poReflnlist->nSelectField(
                                eReflnField_int_type,
				poScalemerge->m_poReflnlist->m_nFI_nL,
                                char(0),
				pcSelect);
	  (void) poScalemerge->m_poReflnlist->nSelectField(
                                eReflnField_float_type,
				poScalemerge->m_poReflnlist->m_nFI_fIntensity,
                                char(0),
				pcSelect);
	  (void) poScalemerge->m_poReflnlist->nSelectField(
                                eReflnField_float_type,
				poScalemerge->m_poReflnlist->m_nFI_fSigmaI,
                                char(0),
				pcSelect);
//+ 2-Mar-1999
	  if (2 == nTexsan)
	    {
	      float fInt, fSigma, fTgeos, fScale;
	      float fLorentz, fPolarz, fLP;

	      int nFI_fTgeos  = poScalemerge->m_poReflnlist->nGetFieldIndex("fTgeos");
	      int nFI_fScale  = poScalemerge->m_poReflnlist->nGetFieldIndex("fScale");

	      int nFI_fLorentz = poScalemerge->m_poReflnlist->m_nFI_fLorentz;
	      int nFI_fPolarz  = poScalemerge->m_poReflnlist->m_nFI_fPolarz;
	      int nFI_fLP     = poScalemerge->m_poReflnlist->nGetFieldIndex("fLPfactor");

	      if (  (0 > nFI_fTgeos) || (0 > nFI_fScale) )
		{
		  cout << "WARNING: transmisson or scale factor not found.\n";
		}
	      if (  (0 > nFI_fLP)
		  && ( (0 > nFI_fLorentz) || (0 > nFI_fPolarz) ) )
		{
		  cout << "WARNING: Lorentz or Polarz factor not found.\n";
		}

	      // Use the 2nd field, whatever it is, to hold the pre-scaled
	      // Intensity, also with no LP correction

	      (void) poScalemerge->m_poReflnlist->nSelectField(
						       eReflnField_float_type,
							       2,
							       char(0),
							       pcSelect);

	      if (0 < nFI_fTgeos)
		(void) poScalemerge->m_poReflnlist->nSelectField(
                                eReflnField_float_type,
								 nFI_fTgeos,
								 char(0),
								 pcSelect);
	      if (0 < nFI_fScale)
		(void) poScalemerge->m_poReflnlist->nSelectField(
                                eReflnField_float_type,
								 nFI_fScale,
								 char(0),
								 pcSelect);

	      Crefln *poRefln;
	      for (i = 0; i < poScalemerge->m_poReflnlist->nGetNumReflns(); i++)
		{
		  poRefln = poScalemerge->m_poReflnlist->poGetRefln(i);
		  fTgeos  = poRefln->fGetField(nFI_fTgeos);
		  fScale  = poRefln->fGetField(nFI_fScale);
		  if (0 < nFI_fLP)
		    {
		      fLP = 1.0f / poRefln->fGetField(nFI_fLP);
		    }
		  else
		    {
		      fLorentz= poRefln->fGetField(nFI_fLorentz);
		      fPolarz = poRefln->fGetField(nFI_fPolarz);
		      fLP     = fLorentz * fPolarz;
		    }

		  fInt   = poRefln->fGetIntensity();
		  fSigma = poRefln->fGetSigmaI();
		  if (   (0.0 < fSigma)
		      && (0.0 < fLP)
		      && (0.0 < fScale)
		      && (0.0 < fTgeos) )
		    {
		      fInt   = fInt   / (fScale / fTgeos);
		      fSigma = fSigma / (fScale / fTgeos);
		      poRefln->vSetIntensity(fInt);
		      poRefln->vSetSigmaI(fSigma);
		      fInt   = fInt * fLP;
		      poRefln->vSetField(2, fInt);
		    }
		  else if (0.0 < fSigma)
		    {
		      cout << "WARNING -texsan2 problem: ";
		      poRefln->nList(2);
		    }
		}
	    }
//- 2-Mar-1999
	  poScalemerge->m_poReflnlist->vSetNoWriteHeader();
	  poScalemerge->m_poReflnlist->nSelect("-fSigmaI<=0.0");
	  cout << "Writing " << sOutScaledFile << " ..." << endl << flush;
	  nStat = poScalemerge->m_poReflnlist->nWrite(sOutScaledFile,
						      NULL,
						      pcSelect);
	  delete [] pcSelect;
	}
      else
	{
	  if (0 != nNoUnavgHeader)
	    poScalemerge->m_poReflnlist->vSetNoWriteHeader();
	  cout << "Writing " << sOutScaledFile << " ..." << endl << flush;
	  nStat = poScalemerge->m_poReflnlist->nWrite(sOutScaledFile);
	}
    }
  if ( (0 == nStat) && ("?" != sOutAverageFile) )
    {
      if (0 != nNoHeader)
	oReflnlistOut.vSetNoWriteHeader();
      cout << "Writing " << sOutAverageFile << " ..." << endl << flush;
      nStat = oReflnlistOut.nWrite(sOutAverageFile);
    }

  poHeader->nWrite(sOut);

  delete poReflnlistIn;
  delete poScalemerge;
  delete poHeader;

  cout << "\ndtscalemerge: Done.\n" << flush;
  return (nStat);
}

void vError(const int nErrorNum, const Cstring& sMessage)
{
  if (0 != nErrorNum)
    {
      cerr << sMessage << '\n';
    }
  cout << "dtscalemerge - Usage:\n"
       << "dtscalemerge input_header_file input_file\n"
       << "                               [ -fixB ]\n"
       << "                               [ -cycles nCycles ]\n"
       << "                               [ -reject fRejSigma ]\n"
       << "                               [ -reso fResoMin fResoMax ]\n"
       << "                               [ -sigma fExcludeSigma ]\n"
       << "                               [ -add   fErrorAdd ]\n"
       << "                               [ -mul   fErrorMul ]\n"
       << "                               [ -fix   sBatchName ]\n"
       << "                               [ -scale fFixedScale ]\n"
       << "                               [ -bfac  fFixedB ]\n"
       << "                               [ -ref   unavg_reflnlist_file ]\n"
       << "                               [ -out   headerfile ]\n"
       << "                               [ -countoverlap ]\n"
       << "                               [ -scaleanom | -anom ]\n"
       << "                               [ -texsan | -texsan2] [ -noheader ]\n"
       << "                               [ output_file ]\n\n"
       << "Command line options:\n\n"
       << "input_header_file\n"
       << "           The header of this file is used to create a crystal\n"
       << "           and spacegroup object.  The spacegroup is required in\n"
       << "           order to reduce reflection indices to the proper\n"
       << "           asymmetric unit.  The crystal unit cell dimensions are\n"
       << "           required to calculate the proper resolution and\n"
       << "           completeness.\n\n"
       << "input_file The name of an input reflection list file.\n"
       << "           If you wish to treat\n"
       << "           multiple reflection lists, first combine them with\n"
       << "           dtreflnmerge:  Default: There is no default, you must\n"
       << "           specify an input reflection list file.\n\n"
       << "-fixB      Fixes the so-called B-factors which are all set to 0.\n"
       << "           They are not refined.  Default: B-factors are refined.\n\n"
       << "-cycles nCycles\n"
       << "           Set the maximum number of non-linear least squares\n"
       << "           cycles to perform.  If the refinement converges, then\n"
       << "           the maximum number of cycles may not be reached.  No\n"
       << "           shifts are applied the last cycle.  Default: 30.\n\n";
  cout << "-reject fRejSigma\n"
       << "           Sets the rejection level for reflections.  Reflections\n"
       << "           with scaled intensities that differ by more than fRejSigma\n"
       << "           from the weighted average intensity calculate by other\n"
       << "           symmetry-related reflections are flagged as rejected by\n"
       << "           setting their observed standard deviations to be negative.\n"
       << "           Default: 1x10e11, so no reflections are rejected.\n\n"
       << "-reso fResoMin fResoMax\n"
       << "                      Resolution limits of input reflections to use.\n"
       << "                      Default: 999999 0.0\n\n"
       << "-sigma fExcludeSigma\n"
       << "           Input reflections with I/sigmaI less than fExcludeSigma\n"
       << "           are excluded from contributing to the scale factor refine-\n"
       << "           ment.  However, these reflections are included in the\n"
       << "           final statistics.  Default: 3.\n\n"
       << "-mul   fErrorMul\n"
       << "-add   fErrorAdd\n"
       << "           These factors alter the error model and weights applied\n"
       << "           to the observed reflections according to\n"
       << "           Whj = 1 / [ (Sighj * fErrorMul)^2 + (Ihj * fErrorAdd)^2]\n"
       << "           Defaults: fErrorMul = 1.  fErrorAdd = 0.\n\n"
       << "-fix   sBatchName\n"
       << "           sBatchName specifies the name of the batch whose scale\n"
       << "           factors will remain fixed.  The scale factors of all\n"
       << "           other batches will shift relative to this batch.\n"
       << "           Default: first batch in the input reflection list.\n\n"
       << "-scale fFixedScale\n"
       << "           fFixedScale specifies the scale factor kj of the fixed\n"
       << "           batch.  Default: 1.\n\n"
       << "-bfac  fFixedB\n"
       << "           fFixedB specifies the B-factor Bj of the fixed batch\n"
       << "           Default: 0.\n\n"
       << "-ref   unavg_reflnlist_file\n"
       << "           This options writes the scaled, but unaveraged, reflection\n"
       << "           list to the file unavg_reflnlist_file.  Rejected reflections\n"
       << "           will have standard deviations less than 0.\n"
       << "           Default: Do not write such a file.  See also -texsan.\n\n"
       << "-scaleanom Treat I+ and I- reflections separately during scaling and\n"
       << "           and rejection.  Output is hkl,Intensity,sigmaI with I- and\n"
       << "           I+ on separate lines (each only if present in the reflnlist).\n"
       << "           This is different from the -anom option which, if set, is\n"
       << "           turned off.  Default: I+ and I- are assumed to be equal.\n\n"
       << "-anom      Output anomalous information as Intensity,sigmaI,Intensity+,\n"
       << "           sigmaI+,Intensity-,sigmaI-.  Use -1 for missing information.\n"
       << "           Default: Output Intensity,sigmaI; do not output I+,sigI+,I-,sigI-\n"
       << "           on the same line.\n\n"
       << "-texsan    Output only h,k,l,Intensity,sigmaI fields with no header and\n"
       << "           no rejected measurements in the scaled, unaveraged file.\n\n"
       << "-texsan2   Output only h,k,l, unscaled Intensity, unscaled SigmaI, unscaled &\n"
       << "           non-LP-corrected Intensity, transmission factor, scale factor\n\n"
       << "-noheader  Do NOT write standard d*TREK reflnlist header on the scaled\n"
       << "           and averaged output file. Default: write a header.\n\n"
       << "-nounavgheader\n"
       << "           Do NOT write standard d*TREK reflnlist header on the unaveraged\n"
       << "           output file. Default: write a header.\n\n"
       << "-countoverlap\n"
       << "           Count and show overlaps among the scaling batches.\n"
       << "           Default: do not count nor show overlaps.\n\n"
       << "-out headerfile\n"
       << "           Write a header to the headerfile at the end of dtscalemerge\n\n"
       << "output_file\n"
       << "           The name of the output reflection list file. The name\n"
       << "           cannot begin with - character.  Default: Do not write\n"
       << "           an output scaled and averaged reflection list file.\n"
       << "           See also -anom and -noheader.\n\n"
       << endl << flush;
#ifndef SSI_PC
  exit (nErrorNum);
#endif
}
