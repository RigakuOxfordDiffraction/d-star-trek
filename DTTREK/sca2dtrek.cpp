//
// Copyright (c) 2010 Rigaku Americas, Corp.
//
// sca2dtrek.cc   Initial author: J.W. Pflugrath           
//    Convert a scalepack .sca file with anom to a d*TREK reflnlist file 
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
#include "Csource.h"

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

int main(int   argc,     // Command line argument count
         char *argv[])   // Pointers to command line args
     
{
  int          nStat;
  Cstring      sInFile       = "nm_output.sca";
  Cstring      sTemplate;
  Cstring      sOutFile      = "nm_output.ref";
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
  float        fWavelength = 0.0;
  float        fMosaicity = 0.0;
  bool         bDoScale = TRUE;

  enum eSca_filetype {
    eSca_unknown_type,
    eSca_noanom_type,
    eSca_anom_type,
    eSca_nomerge_type
  };

  eSca_filetype eSca_type = eSca_unknown_type;

  Cstring        sMissingOption
                 = "ERROR: sca2dtrek - missing option argument!";
  Cstring        sInvalidArgument
                 = "ERROR: sca2dtrek - invalid argument!";

  vDtrekSetModuleName("sca2dtrek");
  vPrintCopyrightInfo();

  // 
  // Usage: sca2dtrek 

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
  else if (0 == argc)
    {
      if (!bFileExists(sInFile))
	vError(6, "ERROR input file not found.");
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
      else if ("-noscale" == sTemp) 
	{
	  bDoScale = FALSE;
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
	  else
	    {
	      oCrystal.vSetCell(a6fCellArg);
	    }
        }
      else if ("-wavelength" == sTemp)
	{
          i++;
          if (i < argc)
            {
              nStat = sscanf(argv[i], "%f", &fTemp);
              if (1 == nStat)
		{
		  fWavelength = fTemp;
		}
              else
                DTREK_ERROR(6, sInvalidArgument);
            }
          else
            {
              DTREK_ERROR(5, sMissingOption);
            }
	}
      else if ("-mosaicity" == sTemp)
	{
          i++;
          if (i < argc)
            {
              nStat = sscanf(argv[i], "%f", &fTemp);
              if (1 == nStat)
		{
		  fMosaicity = fTemp;
		  oCrystal.vSetMosaicity((double)fMosaicity);
		}
              else
                DTREK_ERROR(6, sInvalidArgument);
            }
          else
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
		  oCrystal.m_poSpacegroup->vSet(nSpacegroupArg);
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
  // Remove extension and replace with .ref
  Cstring sExtension = sFileGetExtension(sOutFile);
  sBasename = sOutFile.before(sExtension); 
  // sBasename should end with a .
  sOutFile  = sBasename + "ref";  
  cout << "Input scalepack file: " << sInFile
       << "\nOutput d*TREK   file: " << sOutFile
       << "\n\n";

  // Some error checking

  nStat       = 0;

  // Try to open the file

  if (!bFileExists(sTransSymbol(sInFile)))
    {
      nStat = -2;
      cout << "ERROR failed to open file: " << sInFile << endl;
      return (nStat);
    }
  ifstream oIn( sTransSymbol(sInFile));
  if ((oIn.rdbuf()->is_open()))
    {
      if (!oIn)
	{
	  nStat = -2;
	  cout << "ERROR failed to open file: " << sInFile << endl;
	  return (nStat);
	}
    }
  int    nH, nK, nL;
  int    nScaReducedH, nScaReducedK, nScaReducedL;
  int    nHKLBatchNum, nHKLCentricFlag, nHKLSpindleFlag, nHKLAsymmUnit;

  float a6dCell[6];
  float dInt, dSig, dIntP, dIntM, dSigP, dSigM;
  char a30Spg[30];
  int   nNumSymOps = 0;
  strncpy(a30Spg, "None", 4);

  // Read a header ...

  getline(oIn, sTemp);
  //  Look at sTemp to try to discover which of the 3 kinds of .sca files 
  //  that it might be that sca2dtrek supports
  // no merge type has "num_sym_operators  spacegroup_name"
  // noanom and anom have the number 1
  Cstring asTemp[10];
  nStat = split(sTemp, asTemp, 3, ' ');
  if (3 <= nStat)
    {
      cout << "ERROR: Please check input file format.  It may be unsupported by sca2dtrek.\n";
      exit (1);
    }
  if (2 == nStat)
    {
      // Must be a nomerge
      // asTemp[0] has number of symmetry operator lines
      // asTemp[1] has the spacegroup name

      nNumSymOps = atoi(asTemp[0].string());
      strncpy(a30Spg, asTemp[1].string(), asTemp[1].length());
      eSca_type = eSca_nomerge_type;

      // Skip the symmetry operators
      for (i = 0; i < nNumSymOps; i++) 
	{
	  getline(oIn, sTemp);
	  getline(oIn, sTemp);
	}
      // Note: the unit cell is NOT in the file, but may be in a
      // nearby scalepack.sca file.
      // Next line will be the very first reflection
      nStat = 12;  // This makes sure we will enter the while loop below
    }
  else
    {
      // NOT nomerge format
      
      getline(oIn, sTemp);  // A negative number like -985
      getline(oIn, sTemp);  // The unit cell and spacegroup
      nStat = sscanf(sTemp.string(), "%10f%10f%10f%10f%10f%10f%s",
		     &a6dCell[0], &a6dCell[1], &a6dCell[2], 
		     &a6dCell[3], &a6dCell[4], &a6dCell[5], &a30Spg[0]);
      cout << "Unit cell in input file:  "
	   << a6dCell[0] << ", " << a6dCell[1] << ", " << a6dCell[2] << ", " 
	   << a6dCell[3] << ", " << a6dCell[4] << ", " << a6dCell[5] 
	   << "\nSpacegroup: " << a30Spg
	   << endl;
      if (7 != nStat)
	{
	  cout << "ERROR failed to read the header of file " << sInFile
	       << endl;
	  oIn.close();                 // Make sure we always close the stream!
	  return (-2);
	}
      // Next line is the first refln in the file
    }
  // Create a reflnlist to hold the reflection information

  Creflnlist *poReflnlist;

  poReflnlist = new Creflnlist(1);

  int nFI_nScaReducedH;
  int nFI_nScaReducedK;
  int nFI_nScaReducedL;
  int nFI_nHKLBatchNum;
  int nFI_nHKLCentricFlag;
  int nFI_nHKLSpindleFlag;
  int nFI_nHKLAsymmUnit;

  if (eSca_unknown_type == eSca_type)
    {
      poReflnlist->nExpandGetField(Creflnlist::ms_sfIntensityPlus);
      poReflnlist->nExpandGetField(Creflnlist::ms_sfSigmaIPlus);
      poReflnlist->nExpandGetField(Creflnlist::ms_sfIntensityMinus);
      poReflnlist->nExpandGetField(Creflnlist::ms_sfSigmaIMinus);
    }
  else if (eSca_nomerge_type == eSca_type)
    {
      nFI_nScaReducedH    = poReflnlist->nExpandGetField("nScaReducedH");
      nFI_nScaReducedK    = poReflnlist->nExpandGetField("nScaReducedK");
      nFI_nScaReducedL    = poReflnlist->nExpandGetField("nScaReducedL");
      nFI_nHKLBatchNum    = poReflnlist->nExpandGetField("nHKLBatchNum");
      nFI_nHKLCentricFlag = poReflnlist->nExpandGetField("nHKLCentricFlag");
      nFI_nHKLSpindleFlag = poReflnlist->nExpandGetField("nHKLSpindleFlag");
      nFI_nHKLAsymmUnit   = poReflnlist->nExpandGetField("nHKLAsymmUnit");

      // If neither unit cell nor spacegroup given on the command line,
      // see if a nearby scalepack.sca file has those parameters. 
      if (0.0 >= (a6fCellArg[0] * a6fCellArg[1] * a6fCellArg[2] * a6fCellArg[3] * a6fCellArg[5] *a6fCellArg[5]))
	{
	  // No cell specified, so look for one in the cwd.

	  Cstring sScalepackScaFile = "scalepack.sca";
	  ifstream oScalepackSca(sScalepackScaFile);
	  if ((oScalepackSca.rdbuf()->is_open()))
	    {
	      if (!oScalepackSca)
		{
		  nStat = -2;
		  cout << "WARNING failed to open file: " << sScalepackScaFile << endl;
		}
	      else
		{
		  // scalepack.sca file is like the NOT nomerge format
		  getline(oScalepackSca, sTemp);  // The number 1
		  getline(oScalepackSca, sTemp);  // A negative number like -985
		  if (-985 != atoi(sTemp.string()))
		    {
		      cout << "WARNING -985 not found in file: " << sScalepackScaFile << endl;
		    }
		  getline(oScalepackSca, sTemp);  // The unit cell and spacegroup
		  nStat = sscanf(sTemp.string(), "%10f%10f%10f%10f%10f%10f%s",
				 &a6dCell[0], &a6dCell[1], &a6dCell[2], 
				 &a6dCell[3], &a6dCell[4], &a6dCell[5], &a30Spg[0]);
		  cout << "Unit cell in " << sScalepackScaFile <<  " file:  "
		       << a6dCell[0] << ", " << a6dCell[1] << ", " << a6dCell[2] << ", " 
		       << a6dCell[3] << ", " << a6dCell[4] << ", " << a6dCell[5] 
		       << "\nSpacegroup: " << a30Spg
		       << endl;

		  oCrystal.vSetCell(a6dCell);
		}
	    }
	}
    }
  else
    {
      // Probably should give warning and exit
      cout << "ERROR: Please check input file format.  It may be unsupported by sca2dtrek.\n";
      exit (1);
    }

  Crefln oRefln(poReflnlist);   // Create a reflection object

  // Read in the scalepack refln information

  int nCount = 0;
  int nErrorCount = 0;
  int nDotCount = 50;
  bool bInsert;
  
  while ( (5 == nStat) || (7 == nStat) || (12 == nStat) )
    {
      bInsert = TRUE;
      
      getline(oIn, sTemp);
      if (oIn.eof())
	{
	  cout << "EOF\n";
	  break; // Get out of the while loop
	}

      // Parse the string

      if ( (eSca_nomerge_type == eSca_type) && (37 < sTemp.length()) ) 
	{
	  // There is a problem if the Intensity is LARGE and thus no space
	  // occurs in column 38 
	  // or if Sigma is LARGE and no space occurs in column 46
	  // so add spaces
	  // Test Sigma first, so column count for Intensity is not messed up
	  if (' ' != sTemp.GetAt(45))
	    {
	      sTempTemp = sTemp.before(45) + ' ' + sTemp.after(44);
	      sTemp = sTempTemp;
	    }
	  if (' ' != sTemp.GetAt(37))
	    {
	      sTempTemp = sTemp.before(37) + ' ' + sTemp.after(36);
	      sTemp = sTempTemp;
	    }
  //                                        1  2  3  4  5  6  7  8  9 10 11 12
  //123 123 123 123 123 123 12345 1 1 12 1234567 1234567
  //   0   0   6   0   0   6   512 0 1  1336831.1 51328.0
  //   0   0   7   0   0   7   174 0 0  1 66903.4 11114.7
	  nStat = sscanf(sTemp.string(), "%4d%4d%4d%4d%4d%4d%6d%2d%2d%3d%8f%8f",
			 &nH, &nK, &nL, &nScaReducedH, &nScaReducedK, &nScaReducedL,
			 &nHKLBatchNum, &nHKLCentricFlag, &nHKLSpindleFlag, &nHKLAsymmUnit,
			 &dInt, &dSig);
	  if (12 != nStat)
	    {
	      // Did not find enough values, probably EOF or error
	      bInsert = FALSE;
	      if (0 < sTemp.length())
		{
		  cout << "ERROR parsing this string:\n>>>" << sTemp << "<<<\n\n" << flush;
		  nErrorCount++;
		}
	    }
	  else
	    {
	      oRefln.vSetH(nH);
	      oRefln.vSetK(nK);
	      oRefln.vSetL(nL);
	      oRefln.vSetField(nFI_nScaReducedH, nScaReducedH);
	      oRefln.vSetField(nFI_nScaReducedK, nScaReducedK);
	      oRefln.vSetField(nFI_nScaReducedL, nScaReducedL);
	      oRefln.vSetField(nFI_nHKLBatchNum , nHKLBatchNum);
	      oRefln.vSetField(nFI_nHKLCentricFlag, nHKLCentricFlag);
	      oRefln.vSetField(nFI_nHKLSpindleFlag, nHKLSpindleFlag);
	      oRefln.vSetField(nFI_nHKLAsymmUnit, nHKLAsymmUnit);
	      oRefln.vSetIntensity(dInt);
	      oRefln.vSetSigmaI(dSig);
	    }
	}
      else // anom or noanom
	{
	  nStat = sscanf(sTemp.string(), "%4d%4d%4d%8f%8f%8f%8f",
			 &nH, &nK, &nL, &dIntP, &dSigP, &dIntM, &dSigM);
	  if (5 == nStat)
	    {
	      // Missing Int-, sigI-
	      dIntM = -1.0;
	      dSigM = -1.0;
	      dInt  = dIntP;
	      dSig  = dSigP; 
	    }
	  else if (7 == nStat)
	    {
	      // This means we had at least one (1) reflection with I+ and I- 
	      // so this is a  eSca_anom_type
	      eSca_type = eSca_anom_type;

	      if (0.0 < dSigP)
		{
		  if (0.0 < dSigM)
		    {
		      dInt = dIntP * dSigP  +  dIntM * dSigM;
		      dInt = dInt / (dSigP + dSigM);
		      dSig = (dSigP * dSigM) / (dSigP + dSigM);
		    }
		  else
		    {
		      dInt  = dIntP;
		      dSig  = dSigP; 
		    }
		}
	      else
		{
		  dInt  = dIntM;
		  dSig  = dSigM; 
		}
	    }
	  else
	    {
	      // Probably EOF or error
	      bInsert = FALSE;
	      if (0 < sTemp.length())
		{
		  cout << "ERROR parsing this string:\n>>>" << sTemp << "<<<\n\n" << flush;
		  nErrorCount++;
		}
	    }
	  if (bInsert)
	    {
	      oRefln.vSetH(nH);
	      oRefln.vSetK(nK);
	      oRefln.vSetL(nL);
	      oRefln.vSetIntensity(dInt);
	      oRefln.vSetSigmaI(dSig);
	      oRefln.vSetField(poReflnlist->m_nFI_fIntensityPlus, dIntP);
	      oRefln.vSetField(poReflnlist->m_nFI_fSigmaIPlus, dSigP);
	      oRefln.vSetField(poReflnlist->m_nFI_fIntensityMinus, dIntM);
	      oRefln.vSetField(poReflnlist->m_nFI_fSigmaIMinus, dSigM);
	    }
	}

      // For all refln types
      if (bInsert)
	{
	  poReflnlist->nInsert(&oRefln);  // Insert only if success
	  nCount++;
	  if (0 == (nCount % (nDotCount * 4)))
	    cout << ".";
	  if (0 == (nCount % (nDotCount * 280)))
	    {
	      cout << "\n";
	      //cout << " " << nCount << "\n";
	      nDotCount *= 10;
	    }
	}
    } // end of while loop

  //cout << "\nnStat = " << nStat << "\n";

  if (0 < nCount)
    {
      cout << "\nINF0: read in " << nCount << " reflections.\n\n";

      if (eSca_nomerge_type == eSca_type)
	{
	  a30Spg[0] = toupper(a30Spg[0]);  // Capitalize spacegroup name
	  if (0 >= nSpacegroupArg)
	    oCrystal.m_poSpacegroup->vSet(Cstring(a30Spg));

	  // Only put the cell stuff in if we have a unit cell AND it was NOT overridden by the command line

	  if (0.0 >= (a6fCellArg[0] * a6fCellArg[1] * a6fCellArg[2] * a6fCellArg[3] * a6fCellArg[5] *a6fCellArg[5]))
	    {
	      oCrystal.vSetCell(a6dCell);
	    }
	  else
	    {
	      // External unit cell is probably OK
	      poReflnlist->nPutCrystal(oCrystal);
	    }

	  poReflnlist->eSetWriteBinary(eReflnoutput_text);
	  nStat = poReflnlist->nWrite(sOutFile);
	}
      else if (eSca_anom_type != eSca_type)
	{
	  // It must be eSca_noanom_type since no I-'s were found in the file.
	  // Thus we can skip 4 fields when writing out the results
	  // Select only the 5 fields we want to keep

	  int nFields;
	  char *pcSelect = NULL;
	  nFields  = poReflnlist->nGetNumFields();
	  pcSelect = new char [nFields + 1];

	  // "Deselect" all fields

	  for (i = 0; i < nFields; i++)
	    pcSelect[i] = poReflnlist->m_pcSelect[i] + 1;

	  // Keep the global selection the same

	  pcSelect[nFields] = poReflnlist->m_pcSelect[nFields];

	  (void) poReflnlist->nSelectField(
                                eReflnField_int_type,
			        poReflnlist->m_nFI_nH,
			        char(0),
	                        pcSelect);
	  (void) poReflnlist->nSelectField(
                                eReflnField_int_type,
				poReflnlist->m_nFI_nK,
                                char(0),
				pcSelect);
	  (void) poReflnlist->nSelectField(
                                eReflnField_int_type,
				poReflnlist->m_nFI_nL,
                                char(0),
				pcSelect);
	  (void) poReflnlist->nSelectField(
                                eReflnField_float_type,
				poReflnlist->m_nFI_fIntensity,
                                char(0),
				pcSelect);
	  (void) poReflnlist->nSelectField(
                                eReflnField_float_type,
				poReflnlist->m_nFI_fSigmaI,
                                char(0),
				pcSelect);
				    
	  poReflnlist->eSetWriteBinary(eReflnoutput_text);
	  nStat = poReflnlist->nWrite(sOutFile, NULL, pcSelect);
	}
      else
	{
	  poReflnlist->eSetWriteBinary(eReflnoutput_text);
	  nStat = poReflnlist->nWrite(sOutFile);
	}
      if (0 != nStat)
	{
	  cout << "ERROR writing " << sOutFile << endl;
	}
    }
  else
    {
      cout << "\nWARNING, no reflnlist written since no reflns!\n\n";
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

  oIn.close();                 // Make sure we always close the stream!

  if ( (0 == nStat) && (bDoScale) )
    {
      // Now simply run dtscaleaverage that the user can run
      Cimage_header oHead;
      // We need a unit cell and a spacegroup in order to do scaling
      oCrystal.nUpdateHeader(&oHead);
      oCrystal.nList(2);

      // Update with a wavelength as well
      Csource oSource;
      if (0.0 < fWavelength)
	oSource.m_poWavelength->nSetWavelength(fWavelength);
      oSource.nUpdateHeader(&oHead);

      Cstring sHeadName;
      Cstring sLogFile = sBasename.substr(0, sBasename.length()-1) + "_dtscaleaverage.log";
      sHeadName = sBasename + "head";
      oHead.nWrite(sHeadName);
      sTemp = "dtscaleaverage " + sHeadName + " " + sOutFile 
	+ " -anom -errormodel 1 0.000001 -reject fraction .0000001 " 
	+ sBasename.left((int)sBasename.length()-1) + "_dtscale.ref" 
        + " 1> " + sLogFile + " 2>&1 &";
      cout << "\nINFO, starting dtscaleaverage script.  Log file is " << sLogFile << "\n\n";
      
      (void) nFileAppendVersion(sLogFile, TRUE);
      nStat = nDoSystemCommand(sTemp);
    }
  else
    {
    }
  return (nStat);
}

void vError(const int nErrorNum, const Cstring& sMessage)
{
  if (0 != nErrorNum)
    {
      cout << sMessage << endl;
    }
  cout << "\nsca2dtrek - Convert an HKL-2000 *.sca file to a d*TREK .ref file and\n"
       << "                      run dtscaleaverage to get nice statistics and a valid\n"
       << "                      Table 1 'Summary of data collection statistics'.\n\n"
       << "Usage:\n"
       << "   sca2dtrek [output_sca_file] [ options ...]\n\n"
       << "Command line options: Description:\n\n"

       << " output_sca           Optional name of a scalepack output file to convert.\n"
          "                      Please use the 'no merge original index' option/macro to\n"
          "                      scalepack if you want statistics.  Default input file is\n"
          "                      nm_output.sca.\n\n" 
          " -cell a b c alp bet gam\n"
          "                      Specify the unit cell.  Default from the .sca file if it\n"
          "                      has unit cell or from a scalepack.sca file in the current\n"
          "                      working directory.\n\n"
          " -spacegroup nSpace   Spacegroup to use.  Default from the input .sca file or\n"
          "                      scalepack.sca file or 0.\n\n"
          " -wavelength fWavelength  Wavelength to use when scaling.  Default: 1.\n\n"
          " -nostats             Do not run dtscaleaverage after conversion.  Default is\n"
          "                      to run dtscaleaverage with suitable arguments to calculate\n"
          "                      statistics.\n\n"
          "Output file will be input file basename + .ref.  A d*TREK .head file will be\n"
          "                      created when dtscaleaverage is run with the log file\n"
          "                      named with the same basename plus _dtscaleaverage.log\n"
          "                      appended.\n\n"
          "Examples:\n\n"
          "   sca2dtrek\n"
          "   sca2dtrek nm_output.sca\n"
          "   sca2dtrek nm_output.sca -nostats\n"
          "   sca2dtrek nm_output.sca -cell 78.8 78.8 37.9 90 90 90 -spacegroup 96\n\n"
          "Notes:\n"
          " HKL and scalepack include the systematically absent Bragg reflections which\n"
          "      they delete from the 'no merge original index' output file in the total\n"
          "      number of observations and unique reflections counts, so these numbers\n"
          "      may not match the results here.\n"
          " scalepack does not calculate <I/sigI>, but reports the somewhat different\n"
          "      <I>/<sigI>.  It is not obvious when scalepack reports the Rmerge or the\n"
          "      RmergeAnom (with I+ and I- kept separate). RmergeAnom is usually lower.\n"
          " If dtscaleaverage is run to get a Table 1, then a default mosaicity of 0.3\n"
          "      degrees and the wavelength (see above) is placed in the output\n"
          "      d*TREK file of the averaged unique reflections.\n\n"
       << endl;
  exit (nErrorNum);
}
