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
// DTMainAxial.cpp       Initial author: RB         14-Jan-2005
// This file contains the member functions of class CDTMainAxial

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

#include "Dtrek.h"

#include "DTMainAxial.h"
#include "CAxialPhoto.h"

#ifdef SSI_PC
#include "CrclHelper.h"
#include "CrclIncludes.h"
#endif

#if !defined(VC6)
using std::cout;
using std::endl;
using std::flush;
#endif

int CDTMainAxial::nExecute(unsigned int argc,     // Command line argument count
                           char        *argv[])   // Pointers to command line args
{
    
    
    
    Cimage_header *poHeader    = NULL;
    
    int          nStat;
    int          i, j;  // Loop counters
    int          nTemp;
    int          a2nTemp[2];
    float        fTemp;
    float        a2fTemp[2];
    Cstring      sTemp;
    
    Cstring      sMissingOption
        = "ERROR: dtaxial - missing option argument!";
    Cstring      sInvalidArgument
        = "ERROR: dtaxial - invalid argument!";
    
    
    
    vDtrekSetModuleName("dtaxial");
    cout << "\ndtaxial:  Copyright (c) 2006 Rigaku\n";
    cout << D_K_DTREKVersion << endl;

#ifdef SSI_PC
    Cstring sCCAppVersion = (const char*) argv[argc-1];
	if( sCCAppVersion.find("CrystalClear") != -1 || sCCAppVersion.find("SCXmini") != -1) {
		cout << sCCAppVersion;
		argc -= 2;
		cout << CCrclHelper::GetInstance()->GetTimestamp() << endl;
		argv[argc] = NULL;
	}
#endif
    
    // Copy command line to output log
    
    cout << "Command line:\n " << sGetCommandLine(argc, argv, 71) << endl << endl << flush;

    sTemp = (const char*) argv[1];
    if ("-help" == sTemp || "-h" == sTemp) {
        DTREK_ERROR(0, "");
    };
    
    // Parse command line arguments
    
    argc--; argv++;
    
    nStat = 1;
    
    // Create default scan and find objects;
    
    
    
    sTemp = (argc>=1)?argv[0]:"";
    if (sTemp != "")
        poHeader = new Cimage_header(sTemp);
    else
        DTREK_ERROR(0,"");
    
    if (!poHeader->bIsAvailable())
    {
        DTREK_ERROR(nStat, "ERROR: dtfind error reading files!");
    }
    CAxialPhoto oAxial(*poHeader);

    if (!oAxial.bIsAvailable()) {
        DTREK_ERROR(nStat,"ERROR:  Could not initialize axial photograph");
    };
    
    cout << "\n\n" << flush;
    
    // From here on we have a *poFind object and a poFind->m_poScan object!
    // Parse any additional command line options
    
    argv++;
    argc--;
    for (i = 0; i < argc; i++)
    {
		if( bIsSkipCommandLineArgument(argv[i]) )
			continue;

		sTemp = (const char*) argv[i];
        cout << "Command line string: >>" << sTemp << "<<" << endl;
        
        if ("-sigma" == sTemp)
        {
            i++;
            if (i < argc)
            {
                nStat = sscanf(argv[i], "%f", &fTemp);
                if (1 == nStat)
                    oAxial.vSetSigma(fTemp);
                else
                    DTREK_ERROR(6, sInvalidArgument);
            }
            else
            {
                DTREK_ERROR(5, sMissingOption);
            }
        }    
        else if ("-reso" == sTemp)
        {
            i++;

            Cstring       sResoParameters("");
            while( i < argc && '-' != argv[i][0] )
            {
                sResoParameters += argv[i];
                sResoParameters += " ";
                i++;
            }

            Cstring       sError("");
            if( !bParseCalculateResolution(sResoParameters, a2fTemp[0], a2fTemp[1], sError) )
            {
              DTREK_ERROR(6, sError);
            }

            i--;
            
            oAxial.vSetResolution(a2fTemp[0],a2fTemp[1]);
        }
        else if ("-axis" == sTemp)
        {
            i++;
            if ((i < argc) && 
                (1==sscanf(argv[i],"%d",&nTemp))) {
                oAxial.vSetAxisCode(nTemp);
            };
        }
        else if ("-rot" == sTemp)
        {
            i++;
            for (j=0; (j < 2) && (i < argc); j++, i++)
            {
                nStat = sscanf(argv[i], "%f", &a2fTemp[j]);
                if (1 != nStat)
                    DTREK_ERROR(6, sInvalidArgument);
            }
            
            i--;
            if (2 != j)
            {
                DTREK_ERROR(5,sMissingOption);
            };
            oAxial.vSetRotRange(a2fTemp[0],a2fTemp[1]);
        } else if ("-seq" == sTemp)
        {
            i++;
            for (j=0; (j < 2) && (i < argc); j++, i++)
            {
                nStat = sscanf(argv[i], "%d", &nTemp);
                if (1 == nStat)
                {
                    a2nTemp[j] = nTemp;
                }
                else
                    DTREK_ERROR(6, sInvalidArgument);
            }
            i--;
            if (2 != j)
            {
                DTREK_ERROR(5, sMissingOption);
            }
            oAxial.vAddSeq(a2nTemp[0],a2nTemp[1]);            
        }
        else if ("-zone" == sTemp) 
        {
            j = -1;
            if ((i+1<argc) && (argv[i+1]==(Cstring) "h"))
                j = 0;
            else if ((i+1<argc) && (argv[i+1]==(Cstring) "k"))
                j = 1;
            else if ((i+1<argc) && (argv[i+1]==(Cstring) "l"))
                j = 2;
            i+=2;
            if ((i<argc) && (1==sscanf(argv[i],"%f",&fTemp)) && (j!=-1)) {
                oAxial.vSetZone(j,(int) fTemp);
            } else
                DTREK_ERROR(6, sInvalidArgument);
        }
        else if ("-zonewidth" == sTemp)
        {
            i++;
            if ((i<argc) && (1==sscanf(argv[i],"%f",&fTemp))) {
                oAxial.vSetZoneWidth(fTemp);
            } else 
                DTREK_ERROR(6,sInvalidArgument);
        }
        else if ("-adjust" == sTemp) 
        {
            oAxial.vSetPixelAdjust(true);
        }        
        else if ("-help" == sTemp)
        {
            DTREK_ERROR(0, "");
        }
        else
        {
            DTREK_ERROR(0, "");
        }
    }

    if (oAxial.nRun()) 
        DTREK_ERROR(0, "");
    return 0;
}




void CDTMainAxial::vError(const int nErrorNum, const Cstring& sMessage)
{
  if (0 != nErrorNum)
    {
      cout << sMessage << endl;
    }
  cout << "\ndtaxial - Usage:\n"
       "dtaxial sHeaderfile [ options ...]\n\n"
       "Command line options: Description:\n\n"
       "sHeaderfile           Name of a file with a header that has the scan\n"
       "                      definition.\n"
       "                      unless overridden by the -2D option.\n\n"
       " -sigma  fSigma       Standard deviation (default 5.0) above back-\n"
       "                      ground for pixels to be mapped.\n\n"
       " -rot fRotMin fRotMax\n"
       "                      Restrict rotation ranges.  These are the goniometer\n"
       "                      axis rotations of the transformed pixels, NOT the\n"
       "                      input scan.\n\n"
       " -zone [h|k|l] fZone\n"
       "                      Restrict either h,k, or l to some specific value\n"
       "                      (+/- the zonewidth which can be adjusted using the\n"
       "                       -zonewidth option)\n"
       "                      Example:  To get hk0 zone use '-zone l 0'\n"
       "                      Default:  No restrictions on h,k, or l\n"
       " -zonewidth fZoneWidth\n"
       "                      Specify the width of the zone to investigate.\n"
       "                      Default:  If no -zone restrictions, has no effect\n\n"
       " -reso    fResoMin fResoMax\n"
       "                      Specifies the resolution range to map in the\n"
       "                      axial photograph.\n\n"
       " -reso edge           Set upper resolution to be equal to the highest\n"
       "                      *minimum* resolution found on any of the four\n"
       "                      edges (sides) of the image.\n\n"
       " -reso corner         Set upper resolution to be equal to the highest\n"
       "                      *maximum* resolution found on any of the four\n"
       "                      edges (sides) of the image.\n\n"
       " -seq    nStart nEnd  Specifies the starting and ending image sequence\n"
       "                      values in a scan to map in the axial photograph\n\n"
       " -adjust              For each spot found, adjust the pixels in the spot\n"
       "                      so that the spot maintains it's shape under the\n"
       "                      transformation.\n\n"
       " -axis nAxisCode\n"
       "                      Specify which (real) crystal axis to rotate about for\n"
       "                      the simulated axial photograph (see code interpretations\n"
       "                      below).\n"
       "                      Axial Codes (nAxisCode):\n"
       "                      100 = a, 010 = b, 001 = c, 110 = a + b\n"
       "                      Default:\n"
       "                      Rotate about the real a axis\n\n"
       "Information about scan angles and mask/non-uniformity of response is\n"
       "determined from the headers of the input image files.  If the\n"
       "headers are incorrect, they may be edited with dtheaderedit.\n\n"
       "Examples:\n\n" << flush;
}


