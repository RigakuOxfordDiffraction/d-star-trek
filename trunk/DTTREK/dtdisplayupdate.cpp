// Copyright (c) 2002 Rigaku/MSC, Inc.
//
// dtdisplayupdate.cc     Initial author: J.W. Pflugrath         07-Jun-2002
//   Send image filename and reflnlist filename updates to dtdisplay
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
// dtdisplayupdate [-d display] [-img image_filename] [-ref reflnlist_filename]
//+ToDo
//
//   Error messages need implementing
//
//+Include files

#include "Dtrek.h"
#include "dtreksys.h"
#include "CXprop.h"

using std::cin;
using std::cout;
using std::endl;
using std::flush;

//+Code begin

//+Definitions, constants, and initialization of static member variables

//+Public functions
void vError(const int nErrorNum=0, const Cstring& sError=NULL);

//

int main(int   argc,      // Command line argument count
         char *argv[])    // Pointers to command line args
{
  int            nStat;
  CXprop        *poXprop      = NULL;
  Cstring      sMissingOption
                 = "ERROR: dtdisplayupdate - missing option argument!";
  Cstring      sInvalidArgument
                 = "ERROR: dtdisplayupdate - invalid argument!";

  int          i, j, k;  // Loop counters
  int          nTemp;
  float        fTemp;
  Cstring      sTemp;
  Cstring      sImage     = "";
  Cstring      sReflnlist = "";
  Cstring      sDisplay   = "";

  float        fBeamCenter0 = 0;
  float        fBeamCenter1 = 0;
  bool         bSetBeamCenter = FALSE;
  Cstring      sBeamCenter;

  // Parse command line arguments

  argc--; argv++;

  nStat = 1;

  for (i = 0; i < argc; i++)
    {
      sTemp = (const char*) argv[i];
      cout << "Command line string: >>" << sTemp << "<<" << endl;

      if ( ("-help" == sTemp) || ("-h" == sTemp) )
        {
          vError(0, "");
	  exit (0);
        }
      else if ("-ref" == sTemp)
        {
          i++;
          if (i < argc)
            {
              sReflnlist = (const char *)argv[i];
            }
          else
            {
              vError(5, sMissingOption);
            }
        }
      else if ("-img" == sTemp)
        {
          i++;
          if (i < argc)
            {
              sImage = (const char *)argv[i];
            }
          else
            {
              vError(5, sMissingOption);
            }
        }
      else if ("-d" == sTemp)
        {
          i++;
          if (i < argc)
            {
              sDisplay = (const char *)argv[i];
            }
          else
            {
              vError(5, sMissingOption);
            }
        }
      else if ("-seq" == sTemp)
        {
          i++;
          if (i < argc)
            {
              nStat = sscanf(argv[i], "%d", &nTemp);
              if (1 == nStat)
		{
		}
              else
                vError(6, sInvalidArgument);
            }
          else
            {
              vError(5, sMissingOption);
            }
        }
      else if ("-beamcenter" == sTemp) 
	{
          if (   (i+2 >= argc)
	      || (1 !=sscanf(argv[i+1],"%f", &fBeamCenter0))
              || (1 !=sscanf(argv[i+2],"%f", &fBeamCenter1))) 
	    {
              fBeamCenter0 = -1;
              fBeamCenter1 = -1;
	    } 
	  else
	    {
	      sBeamCenter = "-beamcenter " + Cstring(fBeamCenter0, 0, 2) + Cstring(fBeamCenter1, 0, 2);
	      bSetBeamCenter = TRUE;
	      // cout << "Beam center>>>" << sBeamCenter << "<<<\n";
	      i+=2;
	    }
	}
    }

  //  cout << "imagefile:>>" << sImage << "<<\n"
  //       << "reffile>>" << sReflnlist << "<<\n" << flush;

  if ("" == sDisplay)
    poXprop = new CXprop("DTDISPLAY_WINDOW_ID", sDisplay.string());
  else
    poXprop = new CXprop("DTDISPLAY_WINDOW_ID");
  
  bool bChangedProp = FALSE;
  if ("" != sImage)
    {
      sTemp = sImage;
      if ("" != sReflnlist)
	{
	  sTemp = sTemp + " Reflnlist: " + sReflnlist;
	}
      poXprop->hSetProperty("DTDISPLAY_IMAGE_UPDATE",
			    sTemp);
      bChangedProp = TRUE;
    }
  else if ("" != sReflnlist)
    {
      poXprop->hSetProperty("DTDISPLAY_REFLN_UPDATE",
			    sReflnlist);
      //			  sGetCWD() + sTemp);
      bChangedProp = TRUE;
    }
  if (bSetBeamCenter) 
    {
      poXprop->hSetProperty("DTDISPLAY_BEAMCENTER_UPDATE",
			    sBeamCenter);
      bSetBeamCenter = FALSE;
      bChangedProp = TRUE;
    }

  // If there was an error changing a property, perhaps an unknown command line argument?
  if (!bChangedProp)
    {
      vError(0, "");
    }

  if (NULL != poXprop)
    delete poXprop;

  return (0);
}

void vError(const int nErrorNum, const Cstring& sMessage)
{
  if (0 != nErrorNum)
    {
      cout << sMessage << '\n';
    }
  cout << "\ndtdisplayupdate - Usage:\n"
       << "dtdisplayupdate [-d display] [-img image_filename] [-ref ref_filename] [Header: header_filename]"
    //       << " [-seq nSeqNum]"
       << "\n\n"
       << "Command line options: Description:\n\n"
       << " NOTE:  All filenames and templates should probably have the full path specification\n"
       << "        since the current working directory will likely not match the file directory.\n\n"
       << "-img image_filename  Name of an image file for dtdisplay to read and display\n\n"
       << "-img 'image_filename Template: /FULL/PATH/template???.img [New!] Seq: seq_num'\n"
       << "                     /FULL/PATH/template???.img is a d*TREK image filename \n"
       << "                     template with a full path specification and\n"
       << "                     seq_num is the integer used to fill in sequence number\n"
       << "                     in the template to get the image filename.  If the Template\n"
       << "                     is different from the one currently in use by dtdisplay, then\n"
       << "                     the file image_filename is read in.  The optional text 'New!'\n"
       << "                     indicates that the template should be treated as a new one and\n"
       << "                     and the image_filename is read in to get many dtdisplay\n"
       << "                     properties and the sequence number is ignored.\n\n"
//       << "-img next           Read in next image, same as pressing 'Next' in dtdisplay\n\n"
//       << "-img prev           Read in previous image, same as pressing 'Next' in dtdisplay\n\n"
       << "-ref   ref_filename  Name of an reflnlist file for dtdisplay to read and\n"
       << "                     display\n\n"
       << "Header: header_file.head\n"
       << "                     This optional argument specifies a valid d*TREK .head file to\n"
       << "                     replace information in the current header used by dtdisplay.  See\n"
       << "                     the example below.\n\n"
  //       << "-seq   nSeqNum       Image sequence number of an image file for dtdisplay\n"
    //<< "                     to read and display\n";
       << "Examples (note the use of the apostrophes):\n\n"
       << "   dtdisplayupdate -img ribo101026a_0001.img\n"
       << "   dtdisplayupdate -img '/usr/Rigaku/ribo101026a_0001.img Template: /usr/Rigaku/ribo101026a_????.img Seq: 15'\n"
       << "   dtdisplayupdate -img ribo101026a_0001.img -ref dtpredict0001.ref\n"
       << "   dtdisplayupdate -ref dtpredict0001.ref\n"
       << "   dtdisplayupdate -img 'ribo101026a_0001.img -ref dtpredict0001.ref Header: dtpredict.head'\n"
       << endl;
}

