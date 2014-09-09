//
// Copyright (c) 2006 Rigaku
//
// RESTRICTED RIGHTS NOTICE SHORT FORM (JUNE 1987)
//
// Use, reproduction, or disclosure is subject to restrictions set
// forth in Contract No. W-31-109-ENG-38 and Contract No.
// 943072401 with the University of Chicago, Operator of
// Argonne National Laboratory.
//
// DTMainRefine.cpp  Initial author: RB           13-Jan-2005
//
//  This file contains the member functions of class CDTMainReflnMerge
//

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
//    dtreflnmerge reads from the command line the filenames of
//    D*TREK reflection files to merge into a single output file.
//    It reads in the input files,  and writes out the output file.
//+ToDo
//
//   Error messages need implementing
//
//+Include files

#include <stdlib.h>

#include <string.h>
#include "Cstring.h"
#include "Crefln.h"
#include "Creflnlist.h"
#include "Cimage_header.h"
#include "Ccrystal.h"
#include "Csource.h"

//+Function prototypes

#ifdef SSI_PC
#include "CrclHelper.h"
#include "CrclIncludes.h"
#endif

#include "DTMainReflnMerge.h"

#if !defined(VC6)
using std::cout;
using std::endl;
using std::flush;
#endif

/*  Any options that should be supported on both input and output need
    to be placed in this routine.  
*/


/*  Note that if sThisFile is of non-zero length, then we can assume that we loaded the
    reflection file without actually loading any data.  Options that do not have an output
    (such as -history) are therefore saved the effort of loading information that is not
    needed.  Options suce as -c= however must know the name of the file if they are expected
    to make modifications.
*/

int CDTMainReflnMerge::nAdditionalOptions(Cstring& sThisFile,Cstring& sOption, Creflnlist& oReflnlist) 
{
  int nx,ny;
  int nStat;
  float a9fMatrix[9];      // Matrix for potential reindexing
  static Cstring sTemp;

  if (0 == sOption.index("-reindex=")) 
    {
      // Reindexing requested!

      sOption = sOption.after("-reindex=");
        
      nStat = sscanf(sOption.string(),
                     "%f,%f,%f,%f,%f,%f,%f,%f,%f",
                     &a9fMatrix[0], &a9fMatrix[3],
                     &a9fMatrix[6], &a9fMatrix[1],
                     &a9fMatrix[4], &a9fMatrix[7],
                     &a9fMatrix[2], &a9fMatrix[5],
                     &a9fMatrix[8]);
      if (9 == nStat)
        {
          nStat = oReflnlist.nReindex(a9fMatrix);
        }
    } 
  else if (0 == sOption.index("-rebatch"))
    {
      Cstring sPatFrom[50];
      Cstring sPrintTo[50];
      Cstring sRemainder;
      Cstring sFrom;
      Cstring sTo;
      int nPatCount;
      int nFieldFrom;
      int nFieldTo;
      int nPat,nRef;
      Stringreg oReg;

      sTemp = sOption.after("-rebatch");
      sRemainder = sOption.after("=");
      sTemp = sTemp.before("=");
      if (!sTemp.length())
        sTemp = Creflnlist::ms_ssBatch;
      else 
        {
          sTemp = sTemp.after("-");
        }
      if (sTemp.after("-").length())
        {
          sFrom = sTemp.before("-");
          sTo = sTemp.after("-");
        } 
      else 
        {
          sFrom = sTemp;
          sTo = sTemp;
        }

      nFieldFrom = oReflnlist.nGetFieldIndex(sFrom);
      if (nFieldFrom<0)
        {
          cout << "Could not find field " << sFrom << "\n" << flush;
          return -1;
        }
      nFieldTo = oReflnlist.nExpandGetField(sTo);
      if (nFieldTo < 0)
        {
          cout << "Could not find/create field " << sTo << "\n" << flush;
          return -1;
        }

      for (nPat=0;sRemainder.length();nPat++) 
        {
          sPatFrom[nPat] = sRemainder.before(",");
          sRemainder = sRemainder.after(",");
          if (!sRemainder.length()) 
            {
              cout << "Failed to recieve output string for " << sPatFrom[nPat] << " in -rebatch command\n" << flush;
              return -1;
            }
          sPrintTo[nPat] = sRemainder.before(",");
          sRemainder = sRemainder.after(",");
        }
      nPatCount = nPat;

      for (nRef=0;nRef<oReflnlist.nGetNumReflns();nRef++) 
        {
          for (nPat=0;nPat<nPatCount;nPat++)
            {
              if (!oReg.nParse(oReflnlist[nRef].sGetField(nFieldFrom),sPatFrom[nPat]))
                {
                  if (oReg.nPrint(sPrintTo[nPat],sTemp)) 
                    {
                      cout << "Invalid printing string: " << sPrintTo[nPat] << "\n" << flush;
                      return -1;
                    }
                  oReflnlist[nRef].vSetField(nFieldTo,sTemp);
                  break;
                }
            }
        }
    } 
  else if (sOption=="-crystal")
    {
      Ccrystal* poCrystal;
      poCrystal=oReflnlist.poGetCrystal();
      if (!poCrystal)
        {
          cout << "WARNING: Could not find crystal information in file.\n" << flush;
        } 
      else 
        {
          cout << flush;
          poCrystal->nList(1);
          delete poCrystal;
        }
    } 
  else if (0 == sOption.index("-crystal="))
    {
      Cstring sHeaderFile = sOption.after("-crystal=");
      Ccrystal* poCrystal;
      
      if (sHeaderFile.length()) {
          Cimage_header oHead(sHeaderFile);
          if (!oHead.bIsAvailable())
          {
              cout << "ERROR: Could not load header '" << sOption.substr(9) << "'\n" << flush;
              return -1;
          }
          poCrystal = new Ccrystal(oHead);
          if (!poCrystal->bIsAvailable())
          {
              cout << "WARNING: Did not retrieve all relevant crystal information from '" << sOption.substr(9) << "'\n" << flush;
              delete poCrystal;
              return -1;
          }
      } else {
          poCrystal = new Ccrystal;
          poCrystal->vSetCell(0.0,0.0,0.0,90,90,90);
      };

      poCrystal->nList(1);
      
      if (sThisFile.length()) {
          // The output file does not exist.
          // Thus, the user wants us to modify this current file.
          if (oReflnlist.nUpdateCrystal(sThisFile,poCrystal)) {
              printf("ERROR:  Unable to update crystal in file '%s'.\n",sThisFile.string());
              return -1;
          };
      } else
          oReflnlist.nPutCrystal(*poCrystal);
  } 
  else if ((sOption=="-history") || 
           (sOption=="-history0") ||
           (sOption=="-history2") ||
               (sOption=="-history1"))
    {
      int nVerbose;
      char cc = sOption.string()[sOption.length()-1];
      
      if (cc=='0')
          nVerbose = 0;
      else if (cc=='1')
          nVerbose=1;
      else if (cc=='2')
          nVerbose=2;
      else
          nVerbose=1;
      oReflnlist.nDisplayLog(nVerbose);
    }
  else if ((!strncmp(sOption.string(),"-c=",3)) || (!strncmp(sOption,"-c+=",4)))
    {
      bool bAdd = FALSE;
      Cstring sData;
      Cstring sPrevData;
      if (!strncmp(sOption.string(),"-c+=",4)) {
          sData = sOption.substr(4);
          bAdd = TRUE;
      } 
          sData = sOption.substr(3);
      oReflnlist.nLogComment(sData,bAdd);          
  }
  else if (   ( 0 == sOption.find("-sort+")) || (0== sOption.find("-sort-")) 
           || (0 == sOption.find("-sort=")))
    {
      char* pcTok;
      int nIndex1;
      int nIndex2;
      bool bReverse;
      eReflnFieldType eType1;
      eReflnFieldType eType2;
        
      // Find the reverse orders.
      if (sOption.string()[5]=='+')
        bReverse=TRUE;
      else
        bReverse=FALSE;

      pcTok = strtok(&sOption.string()[1],"-+");
      pcTok = strtok(NULL,"-+");
      
      if (!pcTok)
        {
          cout << "Could not find sorting field for 'sort' option.\n" << flush;
          return -1;
        }
      switch (pcTok[0]) 
        {
        case 'n': eType1 = eReflnField_int_type; break;
        case 'f': eType1 = eReflnField_float_type; break;
        case 's': eType1 = eReflnField_Cstring_type; break;
        default:
          cout << "Invalid name '" << pcTok <<"' for sorting field.\n" << flush;
          return -1;            
        }
      nIndex1 = oReflnlist.nGetFieldIndex(sTemp = pcTok,eType1);
      if (nIndex1 < 0)
        {
          cout << "Field '" << pcTok << "' not found.\n" << flush;
          return -1;
        }
      if (NULL == (pcTok = strtok(NULL,"-")))
        {
          printf("Sorting ...\n");
          oReflnlist.vSort(eType1,nIndex1,NULL);            
        } 
      else
        {
          switch (pcTok[0])
            {
            case 'n': eType2 = eReflnField_int_type; break;
            case 'f': eType2 = eReflnField_float_type; break;
            case 's': eType2 = eReflnField_Cstring_type; break;
            default:
              cout << "Invalid name '" << pcTok <<"' for sorting field.\n" << flush;
              return -1;            
            }
          nIndex2 = oReflnlist.nGetFieldIndex(sTemp = pcTok,eType2);
          if (nIndex2 < 0)
            {
              cout << "Field '" << pcTok << "' not found.\n" << flush;
              return -1;
            }
          printf("Sorting ...\n");
          oReflnlist.vSort2(eType1,nIndex1,eType2,nIndex2,NULL);
        }
      if (bReverse) 
        {
          int* pnArray;
          // Invert the order of the sort.
          pnArray = oReflnlist.pnGetSortIndex();
          for (nx=0;nx<oReflnlist.nGetNumReflns()/2;nx++)
            {
              ny = pnArray[oReflnlist.nGetNumReflns()-1-nx];
              pnArray[oReflnlist.nGetNumReflns()-1-nx] = pnArray[nx];
              pnArray[nx] = ny;                    
            }
        }   
    } 
  else if (sOption == "-list") {
      eReflnoutput ePrev;
      cout << flush;
      //cerr << flush;
      
      ePrev = oReflnlist.eSetWriteBinary(eReflnoutput_text);
      oReflnlist.nWrite("",NULL,NULL,stdout);
      oReflnlist.eSetWriteBinary(ePrev);
/*****
  } else if (sOption == "-dering") {
        double afIceRingReso[10] = { 3.917, 3.684, 3.458, 2.683, 2.261, 2.081, 1.9584,
            1.9272, 1.8927, 1.7292 };
        Csource oSource(*poHeader);
        oReflnlistIn.nDeleteRing(afIceRingReso,7,1.0,oSource.m_poWavelength->fGetWavelength());
*****/
  } else
    {
      return 1;
    }
  return 0;
}


int CDTMainReflnMerge::nExecute(unsigned int argc,      // Command line argument count
                                char        *argv[])    // Pointers to command line args
{
  char**       ppcArgv;           // Crystal Clear does not want us to mess with argv[].  Instead, mess with ppcArgv !
  int          nNumInputs;        // Number of input lists
  Cstring     *psInFile;          // Names of input files
  Creflnlist  *poReflnlistIn;     // Pointer to input reflection list
  Creflnlist   oReflnlistOut;     // Output reflection list
  Cstring      sOutFile;

  int          nNumFieldsRemoved; // Number of fields to remove/ignore
  Cstring     *psFieldsRemoved;   // Names of fields to remove/ignore
  int          nNumOptions;
  Cstring     *psOptions;
  Cstring      sTemp;
  Cstring      sTemp2;
  int          nStat, i, j,nx,ny;
  int          nNumRefs;
  
  int          nMaxLinesToRead = -1;
  bool         bStructureFactor = FALSE;
  int          nIncludeTagField = 0;                // Do we include the sTag field?  This is induced by -diffs ,-tag, 
  bool         bDoDiffs = FALSE;                    // Are we using the -diffs  ?
  bool         bDummyOutputFile = FALSE;            // Are we assuming a dummy output file?
  Cstring      sDummy = "\t\t\t";
  Cstring      sNull;                               // Null file.
  char *pcRejectBatch = NULL;

  vDtrekSetModuleName("dtreflnmerge");
  vPrintCopyrightInfo();
  //cout << "\ndtreflnmerge: Copyright (c) 2006 Rigaku\n";
  //cout << D_K_DTREKVersion << endl << flush;

#ifdef SSI_PC
  Cstring sCCAppVersion = (const char*) argv[argc-1];
  if( sCCAppVersion.find("CrystalClear") != -1 || sCCAppVersion.find("SCXmini") != -1) {
	cout << sCCAppVersion;
	argc -= 2;
	cout << CCrclHelper::GetInstance()->GetTimestamp() << endl;
	argv[argc] = NULL;
  }
#endif

  cout << "Command line:\n " << sGetCommandLine(argc, argv, 71) << endl << endl << flush;

  // Parse command line arguments, into another buffer.  

  // Copy command line arguments into ppcArgv.
  ppcArgv = new char* [argc];
  for (i = 0; i < argc ; i++)
    {
      ppcArgv[i] = argv[i];
    }

  // Skip program name ... but make room for an extra argument at the end.
  // This is necc. since sometimes, we will be pasting a dummy, and we will require the 
  // extra field space at the end.

  j = 0;
  for (i = 0; i < argc - 1; i++)
    {
      ppcArgv[j++] = ppcArgv[i+1];

      //+jwp 15-Feb-2002, Also fix up argument string so that '-rejectbatch id blah'
      //                  becomes '-rejectbatchid=blah'

      if ( (j > 3) && (0 == strcmp(ppcArgv[j-2], "id"))
           && (0 == strcmp(ppcArgv[j-3], "-rejectbatch")) )
        {
          // Found a sequence of "-rejectbatch id template",
          // so modify it into a single token

          sTemp         = Cstring("-rejectbatchid=") + ppcArgv[i+1];
          pcRejectBatch = new char [sTemp.length()+1];
          strcpy(pcRejectBatch, sTemp.string());
          ppcArgv[j-3]  = pcRejectBatch;
          j             = j-2;
        }
    }
  argc--; 
  if (NULL != pcRejectBatch)
    argc = argc - 2;

  if (1 >= argc)
    {
      DTREK_ERROR(0,"Not enough arguments!\n");  // This calls exit
    }
  // Place a dummy output file at the end of the command line if there is only one 
  // file argument.
  // When we are parsing the data, the dummy will act like an output file, but will
  // not actually get written.
  for (i = 0, j= 0;i<argc;i++) {
      if (*ppcArgv[i]!='-')
          j++;
  };
  if (j<2) {
      argc++;
      ppcArgv[argc-1]=sDummy.string();
      nMaxLinesToRead = 0;
      bDummyOutputFile = TRUE;
  };

  // Allocate pointer arrays, this is excessive but it is simple and
  // we know we have enough space!

  psInFile          = new Cstring [argc];
  nNumInputs        = 0;
  psFieldsRemoved   = new Cstring [argc+10];
  nNumFieldsRemoved = 0;
  psOptions         = new Cstring [argc];
  nNumOptions       = 0;
  poReflnlistIn     = NULL;

  // Get fieldnames that should be removed from all input lists

  int nArg = 0;
  nStat    = 0;

  // Pre-scan to see if we are using the -diffs command.  If so, then we need to set some flags.
  // Also, scan to see if we can do an implicit -maxlines command.  This would happen for the
  // commands -history and -crystal
  for (nArg=argc-1;(nArg>=0) ;nArg--) {
      int nType = 0;

      if (!strncmp(ppcArgv[nArg],"-diffs",strlen("-diffs")))
          nType = 1;
      else if (!strncmp(ppcArgv[nArg],"-tag",strlen("-tag"))) 
          nType = 2;
      
      if (nType == 1) {
          nIncludeTagField = 1;
          bDoDiffs = TRUE;
          break;
      };
      if (nType == 2) {
          nIncludeTagField = 1;
      };
  };
  nArg = 0;
 

  // Remove leading '-'!
  while ( (nArg < argc) && ('-' == (char) *ppcArgv[nArg]) && (0 == nStat) )
    {
      sTemp = psFieldsRemoved[nNumFieldsRemoved] = Cstring((ppcArgv[nArg]+1));
      if (sTemp.before('=')=="maxlines") {
          if (1!=sscanf(sTemp.after('=').string(),"%d",&nMaxLinesToRead)) {
              printf("Could not read argument(s) for option '%s'\n",sTemp.string());
              nStat=1;
          };
      } else if ((sTemp.GetAt(0)!='n') && (sTemp.GetAt(0)!='f') 
                 && (sTemp.GetAt(0)!='s')) {
          printf("ERROR:  Field '%s' is not in correct format!\n",sTemp.string());
          nStat=1;
      } else {
          cout << "Removed field: " << psFieldsRemoved[nNumFieldsRemoved] << endl << flush;
          nNumFieldsRemoved++;
      };
      nArg++;
    }

  // Special processing.

  if (nIncludeTagField != 0) {
      oReflnlistOut.nExpandGetField(Creflnlist::ms_ssTag);
  };

  if (bDoDiffs) {
      oReflnlistOut.nExpandGetField(Creflnlist::ms_snPackedHKL);
  };


  // Now parse input filenames and possible field modifiers

  while ( (nArg < argc) && (0 == nStat) )
    {
      if ( ('-' != (char) *ppcArgv[nArg]) && ('+' != (char) *ppcArgv[nArg]) )
        {
      // See if this is the "dummy" name.
      // In that case, we must terminate after processing any input options. 
      // (i.e., there is no output file (as in dtreflnmerge thad.ref -history)

      if (ppcArgv[nArg]!=sDummy) {
            // This is the name of an input file (or the name of the outputfile)
            cout << "Filename read: " << ppcArgv[nArg] << '\n' << flush;
            psInFile[nNumInputs] = ppcArgv[nArg];
      };

          if (0 < nNumInputs)
            {
              // If there was a previous input file, now is the time to
              // process it.
              // Read in previous reflection file.
              // Treat options properly and place in output file

              cout << "Processing input file: " << psInFile[nNumInputs-1]
                   << endl << flush;

              // Read in the reflection list

              poReflnlistIn = new Creflnlist(psInFile[nNumInputs-1],nMaxLinesToRead);
              if (!poReflnlistIn->bIsAvailable())
                {
                  // Some kind of error with this reflnlist, so quit
                  nStat = 2;
                }
              else
                {
                  // OK, have reflection list
                  // Loop through options and set them in the
                  // reflection list (put them there if they do not exist)
                  // Also make sure they exist in output reflection list.

                  for (i = 0; (i < nNumOptions) && (0 <= nStat); i++)  
                    {
              nStat = nAdditionalOptions((bDummyOutputFile?psInFile[nNumInputs-1]:sNull),psOptions[i], *poReflnlistIn);
                      if (nStat > 0)
                        nStat = poReflnlistIn->nSelect(psOptions[i]);
                    }
                  if (0 <= nStat) nStat = 0;

                  nNumRefs = oReflnlistOut.nGetNumReflns();

                  // OK, have good poReflnlistIn, now insert into output list.
                  // First, make sure all fields exist in output list

                  if (0 == nStat)
                    nStat = oReflnlistOut.nExpandRefln(*poReflnlistIn);
                  if (0 == nStat)
                    {
                      nStat = oReflnlistOut.nInsertListFrom(*poReflnlistIn,
                                                            NULL,
                                                            poReflnlistIn->pnGetSortIndex());
                    }
          
                  // If we were tagging the source of each datum,
                  //   we place the tag here.

                  if (nIncludeTagField != 0)
                    {
                      sTemp = (nNumInputs);
                      for (nx=nNumRefs;nx<oReflnlistOut.nGetNumReflns();nx++)
                        {
                          oReflnlistOut[nx].vSetField(oReflnlistOut.m_nFI_sTag,
                                                      sTemp);
                        }
                    }
                }

              // Reset number of options back to 0,
              //  for next input reflection file

              nNumOptions = 0;
              delete poReflnlistIn;
              poReflnlistIn   = NULL;
            }

          if (ppcArgv[nArg] == sDummy)
            {
              delete [] psInFile;
              delete [] psFieldsRemoved;
              delete [] psOptions;
              delete []  ppcArgv;
              if (NULL != pcRejectBatch)
                delete [] pcRejectBatch;
              return (0);
            }
          nNumInputs++;
        }
      else
        {
          // This is a field to modify for this input file

          psOptions[nNumOptions] = Cstring((ppcArgv[nArg]));
          nNumOptions++;
        }
      nArg++;   // Onto the next argument in command line
    }


  // If we are doing differences, 
  // we need to change the sTag field so that it is correctly formated.

  if ((bDoDiffs) && (nNumInputs-1>0) && (nIncludeTagField==1))
    {
      sTemp = "";
      for (nx=0;nx<nNumInputs-1;nx++)
          sTemp += "0";
      sTemp += "_";

      for (nx=0;nx<oReflnlistOut.nGetNumReflns();nx++)
        {
          sTemp2 = oReflnlistOut[nx].sGetField(oReflnlistOut.m_nFI_sTag);
          sscanf(sTemp2.string(),"%d",&ny);
          sTemp2  = sTemp;
          sTemp2.string()[ny-1]='1';
          sTemp2 += ny;
          oReflnlistOut[nx].vSetField(oReflnlistOut.m_nFI_sTag,sTemp2);
        }
    }

  // Command line parsed

  if ( (0 == nStat) && (1 < nNumInputs) )
    {
      // Last input file was really output file, so write results

      sOutFile = psInFile[nNumInputs-1];

      nFileAppendVersion(sOutFile, TRUE);

      // See if output file already exists by trying to open it.  If it
      // does exist, then that is an error!
      if( bFileExists(sOutFile) )
      {
          sOutFile = "Output file " + sOutFile + " already exists!\n";
          
          DTREK_ERROR(2, sOutFile);
      }

      // Now parse and apply any output file options

      Cimage_header *poHeader;
      Ccrystal      *poCrystal;
      poHeader  = NULL;
      poCrystal = NULL;

      // Check for header information in the current reflection file.
      poHeader = oReflnlistOut.poGetHeader();
      poCrystal = oReflnlistOut.poGetCrystal();

      for (i = 0; i < nNumOptions; i++)
        {
          if ("-h" == psOptions[i].substr(0,2))
            {
              if (poHeader)
                delete poHeader;
              poHeader  = new Cimage_header(psOptions[i].from(2));
              if (poHeader->bIsAvailable())
                {
                  if (poCrystal)
                    delete poCrystal;
                  poCrystal = new Ccrystal (*poHeader);
                }
            }
          else if ("-oblique" == psOptions[i])
            {
              if (NULL == poHeader)
                {
                  DTREK_ERROR(3, "Cannot apply oblique incidence correction since no input header!\n");
                }
              else
                {
                  oReflnlistOut.vObliqueIncidenceCorrection(0, poHeader);
                }
            }
          else if ("-unoblique" == psOptions[i])
            {
              if (NULL == poHeader)
                {
                  DTREK_ERROR(3, "Cannot UNapply oblique incidence correction since no input header!\n");
                }
              else
                {
                  oReflnlistOut.vObliqueIncidenceCorrection(1, poHeader);
                }
            }
          else if (   ("-reduce" == psOptions[i])
                   || ("-reso"   == psOptions[i]) )
            {
              if (NULL == poCrystal)
                {
                  DTREK_ERROR(3, "Cannot reduce since no crystal information!\n");
                }
              else if (!poCrystal->bIsAvailable())
                {
                  DTREK_ERROR(3, "Cannot reduce since no crystal information!\n");
                }
              else if ("-reduce" == psOptions[i])
                {
                  // Reduce reflns to asymmetric unit

                  nStat = oReflnlistOut.nReduce(*poCrystal);
                }
              else if ("-reso" == psOptions[i])
                {
                  // Add fResolution field to all reflections

                  nStat = oReflnlistOut.nAddResol(*poCrystal);
                }
            }
          else if ("-detreso" == psOptions[i])
            {
              if (NULL == poHeader)
                {
                  DTREK_ERROR(3, "Cannot use -detreso since no header!\n");
                }
              else if (!poHeader->bIsAvailable())
                {
                  DTREK_ERROR(3, "Cannot use -detreso since no header!\n");
                }

              // Add fDetResolution field to all reflections

              nStat = oReflnlistOut.nAddDetResol(poHeader);
              if (0 != nStat)
                {
                  DTREK_ERROR(3, "Problem with -detreso!\n");
                }
            } 
          else if ("-binary" == psOptions[i])
            {
              oReflnlistOut.eSetWriteBinary(eReflnoutput_binary);
            }
          else if ("-ellipsoids" == psOptions[i])
          {
              if (oReflnlistOut.m_nFI_fEllipsoidA00>=0) {
                  (void) oReflnlistOut.nExpandGetField(oReflnlistOut.ms_sfEllipseMajorAxis);
                  (void) oReflnlistOut.nExpandGetField(oReflnlistOut.ms_sfEllipseMinorAxis);
                  (void) oReflnlistOut.nExpandGetField(oReflnlistOut.ms_sfEllipseAxisMajorOffset);
                  for (nx = 0; nx < oReflnlistOut.nGetNumReflns(); nx++) {
                      oReflnlistOut[nx].nConvertToAxialEllipse();
                  };
              };
          }
          else if ( ("-text" == psOptions[i]) || ("-ascii" == psOptions[i]) )
            {
              oReflnlistOut.eSetWriteBinary(eReflnoutput_text);
            }
          else if ("-nolog" == psOptions[i])
            {
              oReflnlistOut.nSetWriteLog(eLogoutput_nolog);
            } 
          else if ("-noheader" == psOptions[i])
            {
              oReflnlistOut.vSetNoWriteHeader();
            } 
          else if ("-f" == psOptions[i])
            {
              bStructureFactor = TRUE;
            } 
          else if (!strncmp("-diffs",psOptions[i].string(),strlen("-diffs")))
            {
              int               anConstantFields[50];
              eReflnFieldType   aeConstantFields[50];
              float             afTolerance[50];

              int nNumConstantFields = 0;
              float fDefaultTolerance = 0.10f;                   
              char* pc;

              // Remove cumbersome fields generated by the nReduce() command.

              if (oReflnlistOut.m_nFI_nPackedHKL)
                psFieldsRemoved[nNumFieldsRemoved++] = Creflnlist::ms_snPackedHKL;
              if (oReflnlistOut.m_nFI_nReducedH)
                psFieldsRemoved[nNumFieldsRemoved++] = Creflnlist::ms_snReducedH;
              if (oReflnlistOut.m_nFI_nReducedK)
                psFieldsRemoved[nNumFieldsRemoved++] = Creflnlist::ms_snReducedK;
              if (oReflnlistOut.m_nFI_nReducedL)
                psFieldsRemoved[nNumFieldsRemoved++] = Creflnlist::ms_snReducedL;
              if (oReflnlistOut.m_nFI_nFplusminusFlag)
                psFieldsRemoved[nNumFieldsRemoved++] = Creflnlist::ms_snFplusminusFlag;
              if (oReflnlistOut.m_nFI_nCentPhase)
                psFieldsRemoved[nNumFieldsRemoved++] = Creflnlist::ms_snCentPhase;

              // Discover fields that must remain constant.  
              // These also might have a tolerance associated 
              //    if they are floating point.

              pc=strtok(psOptions[i].string(),"-");
              for (pc=strtok(NULL,"-");pc;) 
                {
                  sTemp = pc;

                  if (pc[0]=='f')
                    aeConstantFields[nNumConstantFields] = eReflnField_float_type;
                  else if (pc[0] == 'n') 
                    aeConstantFields[nNumConstantFields] = eReflnField_int_type;
                  else if (pc[0] == 's') 
                    aeConstantFields[nNumConstantFields] = eReflnField_Cstring_type;
  
                  nx = oReflnlistOut.nGetFieldIndex(sTemp);
                  if (nx == -1)
                    {
                      printf("WARNING:  Could not find field '%s' in reflection list for -diffs.  Ignoring ...\n", sTemp.string());
                      continue;
                    }
                  anConstantFields[nNumConstantFields] = nx;
                  pc=strtok(NULL,"-");
                  if ((pc) && (strspn(pc,"0123456789.")==strlen(pc)))
                    {
                      sscanf(pc,"%f",&afTolerance[nNumConstantFields]);
                      pc = strtok(NULL,"-");
                    }
                  else
                    afTolerance[nNumConstantFields] = fDefaultTolerance;

                  if (nx != -1) 
                    nNumConstantFields++;              
                }

              // Call routine to do -diffs.
              if (oReflnlistOut.nDiffs(nNumConstantFields, aeConstantFields,
                                       anConstantFields, afTolerance))
                {
                  printf("ERROR:  Could not do -diffs on reflections.\n");
                  delete [] psInFile;
                  delete [] psFieldsRemoved;
                  delete [] psOptions;
                  delete []  ppcArgv;
                  if (NULL != pcRejectBatch)
                    delete [] pcRejectBatch;
                  return 1;
                }

              // Do a sort for convenience.
              printf("Sorting ...");
              oReflnlistOut.vSort2(eReflnField_int_type,
                                   oReflnlistOut.m_nFI_nPackedHKL,
                                   eReflnField_Cstring_type,
                                   oReflnlistOut.m_nFI_sTag, NULL);
              printf("done.\n");
            } 
          else 
            {
              nStat = nAdditionalOptions(sNull,psOptions[i],oReflnlistOut);
              if (nStat > 0) 
                nStat = oReflnlistOut.nSelect(psOptions[i]);
            }
        }

      // Do not write the fields that were flagged on command line

      char *pcSelect;
      pcSelect = new char [oReflnlistOut.m_nTotalFieldsPlus1];
      for (i = 0; i < oReflnlistOut.m_nTotalFieldsPlus1; i++)
        {
          pcSelect[i] = char(0);  // Everything selected
        }

      int j, k;

      for (i = 0; i < nNumFieldsRemoved; i++)
        {
          // Look for a matching field name
          // If found, deselect it

          k =  0;
          for (j = 0; j < oReflnlistOut.m_nIntReflnFields; j++)
            {
              sTemp = oReflnlistOut.sGetFieldName(j, eReflnField_int_type);
              if (sTemp == psFieldsRemoved[i])
                {
                  // Add 1 so different from m_pcSelect field

                  pcSelect[k] = oReflnlistOut.m_pcSelect[k] + 1;
                }
              k++;
            }
          for (j = 0; j < oReflnlistOut.m_nFloatReflnFields; j++)
            {
              sTemp = oReflnlistOut.sGetFieldName(j, eReflnField_float_type);
              if (sTemp == psFieldsRemoved[i])
                {
                  // Add 1 so different from m_pcSelect field

                  pcSelect[k] = oReflnlistOut.m_pcSelect[k] + 1;
                }
              k++;
            }
          for (j = 0; j < oReflnlistOut.m_nCstringReflnFields; j++)
            {
              sTemp = oReflnlistOut.sGetFieldName(j, eReflnField_Cstring_type);
              if (sTemp == psFieldsRemoved[i])
                {
                  // Add 1 so different from m_pcSelect field

                  pcSelect[k] = oReflnlistOut.m_pcSelect[k] + 1;
                }
              k++;
            }
        }
//mlw 03/27/98 add to output F and Sigma F
      if (bStructureFactor)
        {
          // Write F's and not Fsquared (intensities)

          Crefln *poRefln;
          float fSigmaI;
          float fIntensity;
          float fStructureFactor;
          float fSigmaF;
          int nSigmaFIndex, nStructureFactorIndex;
/*
          Cstring sNewNames[2];
          sNewNames[0] = "fStructureFactor";
          sNewNames[1] = "fSigmaF";
          oReflnlistOut.nExpandRefln(0, NULL, 2, sNewNames, 0, NULL);
          nStructureFactorIndex = oReflnlistOut.nGetFieldIndex(sNewNames[0])
          nSigmaFIndex          = oReflnlistOut.nGetFieldIndex(sNewNames[1]);
*/
          // Replace the fIntensity and fSigmaI fields with F's and SigmaF's

          nStructureFactorIndex = 0;
          nSigmaFIndex          = 1;

          for (i=0; i< oReflnlistOut.nGetNumReflns(); i++)
            {
              poRefln    = oReflnlistOut.poGetRefln(i);
              fIntensity = poRefln->fGetIntensity();
              fSigmaI    = poRefln->fGetSigmaI();
              if (0.0 < fIntensity)
                {
                  fStructureFactor = (float)sqrt((double)fIntensity);
                  fSigmaF          = fSigmaI / (2.0f * fStructureFactor);
                }
              else
                {
                  fStructureFactor = (float)0.0;
                  fSigmaF          = (float)0.0;
                }
              poRefln->vSetField(nStructureFactorIndex, fStructureFactor);
              poRefln->vSetField(nSigmaFIndex, fSigmaF);
            }
        }
//mlw end changes

      cout << "Writing file " << sOutFile << "..." << endl << flush;
     
      nStat = oReflnlistOut.nWrite(sOutFile, NULL, pcSelect);
      delete [] pcSelect;

      if (NULL != poHeader)
        delete poHeader;

      if (NULL != poCrystal)
        delete poCrystal;
    }

  // Delete new'd memory

  delete [] psInFile;
  delete [] psFieldsRemoved;
  delete [] psOptions;
  delete []  ppcArgv;
  if (NULL != pcRejectBatch)
    delete [] pcRejectBatch;
  if (NULL != poReflnlistIn)
    {
      delete poReflnlistIn;
    }

  return (nStat);
}


void CDTMainReflnMerge::vError(const int nErrorNum, const Cstring& sMessage)
    {
  if (0 != nErrorNum)
    {
      cout << sMessage << '\n';
  } 
  else 
  {
  cout << "dtreflnmerge - Usage:\n"
  "dtreflnmerge [output_field_options] input_file1 [input_file_options]\\\n"
  "                                 [ [input_file2 [input_file_options]\\\n"
  "                                   [...] ] \\\n"
  "                                    output_file [output_file_options]\n\n"
  "Command line options: Description:\n\n"
  "output_field_options\n"
  "     1)\n"
  "              Fields that will NOT appear in output_file.\n"
  "              The syntax is -field_name.  Example:\n"
  "              -fObs_pixel0 -fCalc_Xmm -nPackedHKL\n\n"
  "     2)\n"
  "              -maxlines=nLines"
  "              Specifiy the maximum lines to read in each file.  This is\n"
  "              useful in conjunction with -list.\n\n"
  "input_file1   The name of an input reflection list file.\n"
  "              The name cannot begin with a + or - character.\n\n"
  "input_file_options\n"
  "     1)\n"
  "              Options to specify which reflections in the\n"
  "              input_file are included or ignored.  By default\n"
  "              all reflections are included.  This takes the form:\n"
  "                -fieldresultOPERATORvalue  (no whitespace!)\n"
  "                +fieldresultOPERATORvalue \n"
  "              where\n"
  "                -      means exclude reflections which fit selection\n"
  "                +      means include reflections which fit selection\n"
  "                fieldresult is a single field name or two field names\n"
  "                       of the same type separated by a +, -, *, or /\n"
  "                       character.  With two field names, the math\n"
  "                       operation is performed and the result is used\n"
  "                       to compare to value.\n"
  "                OPERATOR  is one of ==, !=, <, >, >=, <=, &, |\n"
  "                value     is the value for comparison\n"
  "                For strings, the wildcards '*' and '?' can be used with\n"
  "                '==' or '!=' to select substrings.\n"
  "                .LE. .GE. .LT. .GT. .NE. .EQ. .AND. .OR. can be used\n"
  "                instead of <,>,<=,>=,!=, ==, &, or | \n"
  "                Example 1: to ignore reflections with intensities\n"
  "                         between 50 and 100 use:\n"
  "                         -fIntensity<=100.0 +fIntensity<=50.0\n"
  "                Example 2: to ignore reflections in an ice ring use:\n"
  "                         -fResolution<3.85 +fResolution<3.75\n"
  "                Example 3: to ignore reflections with I/sigmaI < 3 use:\n"
  "                         -fIntensity/fSigma<3\n"
  "                Example 4: to match all batches starting with 'X':\n"
  "                         -sBatch!=X*\n"
  "    2)\n"
  "                -fieldname=value\n"
  "              where fieldname is any fieldname, existing or new\n"
  "                    value is the value.\n"
  "              The field is set to value for all reflections in\n"
  "              the reflection list.  The type of field (int, float,\n"
  "              or string) is determined from the first letter of \n"
  "              fieldname (n, f or s) or from the value. Examples:\n"
  "                -sBatch=X001 -nDetNum=1\n\n"
  "    3)\n"
  "                -fieldnameOPERATORvalue\n"
  "              where fieldname is any existing integer or float field\n"
  "                    OPERATOR is one of +=, -=, *= or /=\n"
  "                    value is the value.\n"
  "              Value is added to, subtracted from, multiplied by or\n"
  "              divided into the current field value.  The field is\n"
  "              set to the result.  The type of field (int or float)\n"
  "              is determined from the first letter of \n"
  "              fieldname (n or f) or from the value. Examples:\n"
  "                -nH+=1 -nK-=1  -fSigmaI*=1.4142\n"
  "                -sBatch+=PRE -sBatch-=SUF   # Add prefix and suffix\n\n"
  "              Input_file_options are processed right to left, so the\n"
  "              order can be important.  By default, all reflections\n"
  "              are included.  As many input files as necessary may\n"
  "              be placed on the command line.  Also the input files\n"
  "              may be repeated with different input file options.\n\n"
  "    4)\n"
  "             -rejectbatch id sBatchTemplate\n"
  "             -rejectbatchid=sBatchTemplate\n"
  "              where reflns in specific sBatch id's are deleted.\n"
  "              sBatchTemplate can have wildcards ('?', '*') to reject\n"
  "              a number of batches.  Multiple batch id's must be separated\n"
  "              by a comma.  Ranges are valid as well (but not   \n"
  "              wildcards in ranges!). \n\n"

  "    5)\n"
  "             -reindex=MATRIX\n"
  "              where MATRIX is nine numbers separated by commas with\n"
  "              no intervening whitespace.  Reflection indices are\n"
  "              multiplied by this matrix and rounded-off to give new\n"
  "              reindexed indices.  Example to switch h and l:\n"
  "              -reindex=0,0,1,0,-1,0,1,0,0\n\n"
  "    6)\n"
  "             -rebatch=sPat1,sNew1,sPat2,sNew2...\n"
  "              OR\n"
  "             -rebatch-sFieldFrom=sPat1,sNew1,sPat2,sNew2...\n"
  "              OR\n"
  "             -rebatch-sFieldFrom-sFieldTo=sPat1,sNew1,sPat2,sNew2...\n"
  "              Rewrites batch names in field sFieldFrom into sFieldTo.\n"
  "              The sPat's  are pattern strings.  For each reflection, the\n"
  "              sPat's are scanned to see if any match the data in sFieldFrom.\n"
  "              If a match is found, the corresponding sNew value is written\n"
  "              to the sFieldTo field of the reflection. If sFieldFrom is omitted\n"
  "              it defaults to 'sBatch'. If sFieldTo is omitted it defaults to\n"
  "              the value of sFieldFrom\n"
  "              Each sNew can contain the macros %1 ... %9 which\n"
  "              map to substrings matching (respectively) each wildcard\n"
  "              or star found. %0 is the entire string.\n"
  "              For example:  ABCDEFGH matches A*D??G* with,\n"
  "                            %0 = ABCDEFGH  %1 = BC  %2 = E  %3 = F %4 = H\n\n"
  "    7)\n"
  "              -nonaxial\n"
  "              All non-axial reflections are excluded, leaving\n"
  "              only h00, 0k0, and 00l reflections.\n\n"
  "    8)\n"
  "              -history OR -history1 OR -history2 \n"
  "              Print reflection list history. -history1 and -history2\n"
  "              support level 1 and 2 verbosity in listing. \n"
  "    9)\n"
  "              -crystal OR -crystal= OR -crystal=sHeaderFile \n"
  "              Option -crystal lists crystal and spacegroup information.\n"
  "              Option -crystal=sHeaderFile extracts crystal information from\n"
  "              file sHeaderFile and places it in the reflection file.\n\n"
  "   10)\n"
  "              -c=sComment -OR c+=sComment\n"
  "              Option -c lists comments associated with a reflection file.\n"
  "              Options -c= (-c+=) set (add) comment information.\n"
  "              This command can be used alone, or with an operation \n"
  "              that is modifying the list.  \n"
  "              Example: L1.ref L2.ref -c=\"L1 merged with L2\" L3.ref\n\n"
  "   11)\n"
  "              -sort-sField OR -sort-sField1-sField2\n"
  "              Option -sort sorts on the specified field name \"sField\".\n"
  "              With two fields, -sort sorts primarily on field \"sField1\" and \n"
  "              secondarily on field \"sField2\"\n"
  "              Example: LIN.ref LOUT.ref -sort-fIntensity-fSigmaI\n"
  "              To reverse the order of sort, use a '+' instead of '-'.\n"
  "              Example: LIN.ref LOUT.ref -sort+fIntensity+fSigmaI\n"
  "              If sorting is specified with an input file, the file is\n"
  "              sorted separately before merging.\n\n"
  "   12)\n"
  "              -list\n"
  "              Dump reflection list to standard output.\n\n"
  "   NOTE:  with most shells the !, >, <, * characters must be escaped.\n"
  "output_file   The name of the output reflection list file.\n"
  "              The name cannot begin with a + or - character.  As a\n"
  "              safety precaution this file will not be written if it\n"
  "              already exists.\n\n"
  "output_file_options\n"
  "              As above for input_file_options, plus the following:\n"
  "              -hHEADER\n"
  "                 where the file HEADER contains a valid d*TREK header.\n"
  "                 Use this when crystal information is not available in\n"
  "                 the reflection file.\n"
  "              -reso\n"
  "                 compute resolution in Angstroms for each refln and\n"
  "                 add the field fResolution.  Resolution is computed\n"
  "                 from hkl and unit cell.  Prerequisite: -hHEADER.\n"
  "              -detreso\n"
  "                 compute resolution in Angstroms for each refln and\n"
  "                 add the field fDetResolution.  Resolution is computed\n"
  "                 from pixel position and detector position found in header.\n"
  "                 Prerequisite: -HEADER.\n"
  "              -reduce\n"
  "                 reduce the hkls to asymmetric unit and add the fields\n"
  "                 nPackedHKL, nReducedH, nReducedK, nReducedL, nAnomFlag,\n"
  "                 nCentPhase.  Sort the file on nPackedHKL.  These\n"
  "                 fields may be removed by naming them above in\n"
  "                 output_field_options. Prerequisite: -hHEADER.\n"
  "              -noheader\n"
  "                 do not write the standard d*TREK reflnlist file header\n"
  "                 at the beginning of the output file.  This is useful\n"
  "                 to convert to other formats, but makes the output file\n"
  "                 unusable to d*TREK programs.\n"
  "              -1to2\n"
  "                 convert a reflnlist with anomalous pairs in a single line\n"
  "                 with fields fIntensity+, fIntensity-, etc to a reflnlist\n"
  "                 with two lines per anomalous pair (one line for h+ and one\n"
  "                 line for h-.  Warning: centric reflns may get two lines.\n"
  "              -f\n"
  "                 convert fIntensity and fSigmaI fields to F =\n"
  "                 sqrt(fIntensity) and SigmaF = fSigmaI / (2 * F)\n"
  "                 Negative intensities are set to F=0, SigmaF=0.\n"
  "                 WARNING: The field labels remain fIntensity and fSigmaI.\n\n"
  "              -binary, -text\n"
  "                 Set binary or text output.  Default output set determined by\n"
  "                 the value of DTREK_REFLN_BINARY.\n\n"
  "              -nolog\n"
  "                 Do not log the operation of this call to dtreflnmerge\n\n"
  "              -tag\n"
  "                 Set sTag field in output list to a value ranging from 1 to\n"
  "                 the number of input files on the command line.\n\n"
  "              -diffs OR -diffs-sField1-... OR -diffs-sField1-fTol1-... \n"
  "                 Does a difference comparison of reflections in two or more\n"
  "                 files.  Differences between reflections are encoded in the\n"
  "                 sTag field.  Two reflections are regarded as unique whenever\n"
  "                 H,K, or L differ, OR  when one of the user supplied fields\n"
  "                 (sField1, sField2 ...) differ.  For floating point fields,\n"
  "                 an additional tolerances can be specified.\n\n"
  " Examples:                                                 \n"
  "    dtreflnmerge dtprofit.ref -sBatch+=1 dtprofit1.ref        \n"
  "    dtreflnmerge -rejectbatchid=0001,1001-1002,200? dtprofit1.ref\n"
  "    dtreflnmerge dtprofit.ref dtprofit.reft -text             \n"
  "    dtreflnmerge dtprofit.ref -rebatch=100,099 dtprofitr.ref  \n"
  "    dtreflnmerge file1.ref file2.ref -diffs                   \n"
  "    dtreflnmerge dtprofit.ref -crystal=sHeader                \n"
  "    dtreflnmerge dtprofit.ref -history                        \n"
  "    dtreflnmerge -maxlines=20 dtprofit.ref -list              \n"
  << endl << flush;
  }
}

