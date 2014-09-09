//
// Copyright (c) 2014 Rigaku Americas, Corp.
//
// xds2dtrek.cc   Initial author: J.W. Pflugrath           
//    Convert a XSCALE .hkl file with anom to a d*TREK reflnlist file 
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
  Cstring      sInFile       = "XDS_ASCII.HKL";
  Cstring      sTemplate;
  Cstring      sOutFile      = "XSCALE.ref";
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
  int          nSpacegroupHeader = 0;
  float        fWavelength = 0.0;
  float        fMosaicity = 0.0;
  bool         bDoScale = TRUE;

  enum eXDS_filetype {
    eXDS_unknown_type,
    eXDS_noanom_type,
    eXDS_anom_type,
    eXDS_nomerge_type
  };

  eXDS_filetype eXDS_type = eXDS_unknown_type;

  Cstring        sMissingOption
                 = "ERROR: xds2dtrek - missing option argument!";
  Cstring        sInvalidArgument
                 = "ERROR: xds2dtrek - invalid argument!";

  vDtrekSetModuleName("xds2dtrek");
  vPrintCopyrightInfo();

  // 
  // Usage: xds2dtrek 

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
  cout << "Input XDS file: " << sInFile
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

  int   nH, nK, nL;
  int   anFI_fArg[100];   // The first 5 (0,1,2,3,4) elements are NOT used
  int   nFieldOffset = 5;
  int   nFields = nFieldOffset;
  float a6dCell[6];
  float dInt, dSig, dIntP, dIntM, dSigP, dSigM;
  char  a30Spg[30];
  int   nNumSymOps = 0;
  strncpy(a30Spg, "None", 4);

  Cstring asTokens[100];
  
  // Read in the XDS refln information

  int nVerbose = 0;
  int nNumToks = 0;
  int nCount = 0;
  int nErrorCount = 0;
  int nDotCount = 50;
  bool bInsert;
  
  Creflnlist *poReflnlist = NULL;
  Crefln     *poRefln     = NULL;

  // Create an empty Creflnlist object

  poReflnlist = new Creflnlist(1);
  
  // Read the file line by line, hopefully it starts with a header ...

  getline(oIn, sTemp);
  while (oIn && !oIn.eof())
    {
      ;      bInsert = FALSE;

      // Successful getline
      if (3 < nVerbose)
	{
	  std::cout << "Read line: " << sTemp << std::endl;
	}
      nNumToks = split(sTemp, asTokens, 100, "= ,\t");
      if (0 == nNumToks)
	{
	}
      else if (0 == asTokens[0].find('!'))
	{
	  // It is a comment line.  Perform actions based on next token
	  if ("!UNIT_CELL_CONSTANTS" == asTokens[0])
	    {
	      // std::cout << "Unit cell: " << sTemp << std::endl;
	      a6dCell[0] = atof(asTokens[1].string());
	      a6dCell[1] = atof(asTokens[2].string());
	      a6dCell[2] = atof(asTokens[3].string());
	      a6dCell[3] = atof(asTokens[4].string());
	      a6dCell[4] = atof(asTokens[5].string());
	      a6dCell[5] = atof(asTokens[6].string());
	    }
	  else if ("!SPACE_GROUP_NUMBER" == asTokens[0])
	    {
	      // std::cout << "Spacegroup: " << sTemp << std::endl;
	      nSpacegroupHeader = atoi(asTokens[1].string());
	      oCrystal.m_poSpacegroup->vSet(nSpacegroupHeader);
	      
	    }
	  else if ("!END_OF_HEADER" == asTokens[0])
	    {
	      // std::cout << "End of header: " << sTemp << std::endl;	    
	      // Create an output reflection list with all the known fields
	      if (NULL == poReflnlist)
		{
		  // What happened?  It should have already been created, so ERRO
		  vError(0,"" );
		}
	    }
	  else if (0 == asTokens[0].find("!ITEM_"))
	    {
	      // std::cout << "Column title: " << sTemp << std::endl;
	      // std::cout << "asTokens[0].after() " << asTokens[0].after('_') << std::endl;
	      //  std::cout << "asTokens[1]: " << asTokens[1] << std::endl;
	      // Assume items 1-5 are h, k, l, I, sigI
	      if (5 < atoi(asTokens[1]))
		{
		  // Do something to save the fieldname to add
		  // And insert into poReflnlist

		  Cstring sField;
		  sField = Cstring('f') + asTokens[0].after('_');
		  if (   (NULL != poReflnlist) 
			 && (100 > nFields) )
		    {
		      anFI_fArg[nFields]   = poReflnlist->nExpandGetField(sField);
		      nFields++;
		    }

		}
	    }
	  else if ("!END_OF_DATA" == asTokens[0])
	    {
	      // std::cout << "End of data: " << sTemp << std::endl;
	      // So write out poReflnlist if it is available

	      if (NULL == poReflnlist)
		{
		  // Throw an error
		}
	      else if (0 >= poReflnlist->nGetNumReflns())
		{
		  // Throw an error
		  delete poReflnlist;
		  poReflnlist = NULL;
		}
	      else
		{
		  // Add the header

		  poReflnlist->nPutCrystal(oCrystal);

		  poReflnlist->eSetWriteBinary(eReflnoutput_text);
		  nStat = poReflnlist->nWrite(sOutFile);
		  if (0 != nStat)
		    {
		      cout << "ERROR writing " << sOutFile << endl;
		    }
		}
	    }
	  else
	    {
	      // It is a header line we do not parse just yet
	    }
	}
      else
	{
	  // It should be a reflection
	  if (NULL == poRefln)
	    {
	      // It is the 1st one in the list, so oR
	      poRefln = new Crefln(poReflnlist);   // Create a reflection object
	    }
	  poRefln->vSetH(atoi(asTokens[0].string()));
	  poRefln->vSetK(atoi(asTokens[1].string()));
	  poRefln->vSetL(atoi(asTokens[2].string()));
	  poRefln->vSetIntensity(atof(asTokens[3].string()));
	  poRefln->vSetSigmaI(atof(asTokens[4].string()));
	  for (int i=nFieldOffset; i < nFields; i++)
	    {
	      poRefln->vSetField(anFI_fArg[i], (float) atof(asTokens[i].string()));
	    }
	  poReflnlist->nInsert(poRefln);
	  nCount++;
	  if (0 == (nCount % (nDotCount * 4)))
	    std::cout << ".";
	  if (0 == (nCount % (nDotCount * 280)))
	    {
	      std::cout << "\n";
	      //cout << " " << nCount << "\n";
	      nDotCount *= 10;
	    }
	}
      getline(oIn, sTemp);
    }
  //cout << "\nnStat = " << nStat << "\n";

  if (0 < nCount)
    {
      cout << "\nINF0: read in " << nCount << " reflections.\n\n";

      // Only put in spacegroup from command line if non-zero
      if (0 < nSpacegroupArg)
	oCrystal.m_poSpacegroup->vSet(nSpacegroupArg);

      // Only put the cell stuff in if we have a unit cell AND it was NOT overridden by the command line

      if (0.0 >= (a6fCellArg[0] * a6fCellArg[1] * a6fCellArg[2] * a6fCellArg[3] * a6fCellArg[5] *a6fCellArg[5]))
	{
	  oCrystal.vSetCell(a6dCell);
	}
      // External unit cell is probably OK
      poReflnlist->nPutCrystal(oCrystal);

      poReflnlist->eSetWriteBinary(eReflnoutput_text);
      nStat = poReflnlist->nWrite(sOutFile);

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

  // Close the input file that we opened
  oIn.close();

  // Clean-up poReflnlist and poRefln
  if (NULL != poReflnlist)
    {
      delete poReflnlist;
      poReflnlist = NULL;
    }
  if (NULL != poRefln)
    {
      delete poRefln;
      poRefln = NULL;
    }

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
  cout << "\nxds2dtrek - Convert an XDS ASCII reflection file to a d*TREK .ref file and\n"
       << "                      run dtscaleaverage to get nice statistics and a valid\n"
       << "                      Table 1 'Summary of data collection statistics'.\n\n"
       << "Usage:\n"
       << "   xds2dtrek [XDS_ASCII_file] [ options ...]\n\n"
       << "Command line options: Description:\n\n"

       << " XDS_ASCII_file      Optional name of a XDS or XSCALE output file to convert.\n"
          "                      Please use the 'MERGE=FALSE' option to\n"
          "                      XDS or XSCALE if you want statistics.  Default input file is\n"
          "                      XDS_ASCII.HKL.\n\n" 
          " -cell a b c alp bet gam\n"
          "                      Specify the unit cell.  Default from the XDS file if it\n"
          "                      has unit cell.\n\n"
          " -spacegroup nSpace   Spacegroup to use.  Default from the input XDS file or 0.\n\n"
          " -wavelength fWavelength  Wavelength to use when scaling.  Default: 1.\n\n"
          " -nostats             Do not run dtscaleaverage after conversion.  Default is\n"
          "                      to run dtscaleaverage with suitable arguments to calculate\n"
          "                      statistics.\n\n"
          "Output file will be input file basename + .ref.  A d*TREK .head file will be\n"
          "                      created when dtscaleaverage is run with the log file\n"
          "                      named with the same basename plus _dtscaleaverage.log\n"
          "                      appended.\n\n"
          "Examples:\n\n"
          "   xds2dtrek\n"
          "   xds2dtrek XDS_ASCII.HKL\n"
          "   xds2dtrek XSCALE_ASCII.HKL -nostats\n"
          "   xds2dtrek XSCALE_ASCII.HKL -cell 78.8 78.8 37.9 90 90 90 -spacegroup 96\n\n"
          "Notes:\n"
          " If dtscaleaverage is run to get a Table 1, then a default mosaicity of 0.3\n"
          "      degrees and the wavelength (see above) is placed in the output\n"
          "      d*TREK file of the averaged unique reflections.\n\n"
       << endl;
  exit (nErrorNum);
}
