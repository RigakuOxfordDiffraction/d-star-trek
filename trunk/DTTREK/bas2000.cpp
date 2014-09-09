//
// Copyright (c) 2007 Rigaku Corporation, Rigaku Americas Corporation
//
// RESTRICTED RIGHTS NOTICE SHORT FORM (JUNE 1987)
//
// Use, reproduction, or disclosure is subject to restrictions set
// forth in Contract No. W-31-109-ENG-38 and Contract No.
// 943072401 with the University of Chicago, Operator of
// Argonne National Laboratory.
//
// bas2000.cpp        Initial author: S.Yasukawa           2-Oct-2007
// This file contains BAS2000 style image file routines

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

#include <string.h>
#include "dtreksys.h"
#include "raxis.h"
#include "bas2000.h"
#include "Cimage.h"
#include "Cimage_header.h"
#include "dskio.h"
#include "dtrekdefs.h"
#include "dtrekvec.h"
#include "Cgoniometer.h"

int 
nReadBAS2000Header(const int* pnFile, const Cstring& rsFilename,
		 char *pcBuffer, Cimage_header *poHeader)
{
	int   nStat;       // Local status flag
	float a10fTemp[10];
	int   a10nTemp[10];

	static Cstring sBAS2000 = "BAS2000_";

	Cstring		sTemp;
	Cstring		infFilename = "";
	FILE*		fp;
	_BAS2000_header oHeader;

	char		dummy_buf[256];
	int			idummy;
	char		dummy[256];

	sTemp = rsFilename;
	if(-1 != sTemp.find(".img"))
	{
		infFilename = sTemp.before(".img");
		infFilename += ".inf";
	}
	else if(-1 != rsFilename.find(".try"))
	{
		infFilename = sTemp.before(".try");
		infFilename += ".trf";
	}

	if(!infFilename.length() || !bFileExists(infFilename)) 
	{
		//does not exist inf file
		return 1;
	}

	fp = fopen(infFilename.string(),"rt");
	if (!fp)
	{
		return 1;
	}

	oHeader.xpxl = 0;
	oHeader.zpxl = 0;
	oHeader.ratio = 202.73;
	oHeader.bswap = 1;
	oHeader.ccd = 0;

	fscanf(fp, "%s\n", dummy);
	if(strcmp(dummy, "CCD_IMAGE_FILE") == 0) {
		oHeader.ccd = 1;
		fscanf(fp, "%s\n", dummy);
   		fscanf(fp, "%d",   &idummy);  oHeader.xsiz = idummy*0.001;
   		fscanf(fp, "%d\n", &idummy);  oHeader.zsiz = idummy*0.001;
   		fscanf(fp, "%d\n", &oHeader.bit);
   		fscanf(fp, "%d\n", &oHeader.xpxl);
   		fscanf(fp, "%d\n", &oHeader.zpxl);
		oHeader.byte = (oHeader.bit + 7) / 8;
	} else {
		fscanf(fp, "%s\n", oHeader.date);
		fscanf(fp, "%s\n", dummy);

   		fscanf(fp, "%d",   &idummy);  oHeader.xsiz = idummy*0.001;
   		fscanf(fp, "%d\n", &idummy);  oHeader.zsiz = idummy*0.001;
   		fscanf(fp, "%d\n", &oHeader.bit);
   		fscanf(fp, "%d\n", &oHeader.xpxl);
   		fscanf(fp, "%d\n", &oHeader.zpxl);
		fscanf(fp, "%ld\n", &oHeader.sens);
		fscanf(fp, "%lf\n", &oHeader.lati);
		fgets(dummy_buf, 256, fp);
		fgets(dummy_buf, 256, fp);
		fgets(dummy_buf, 256, fp);
		fgets(dummy_buf, 256, fp);
		//fgets(dummy_buf, 256, fp);
		fscanf(fp, "%s\n", oHeader.cmnt);
		//sprintf(oHeader.cmnt, "%s", dummy_buf);
		oHeader.byte = (oHeader.bit + 7) / 8;
	}

	fclose(fp);

	if ( (0 >= oHeader.xpxl) || (0 >= oHeader.zpxl) )
	{
		return 1;
	}

	nStat = 0;

	a10nTemp[0] = oHeader.xpxl;
	a10nTemp[1] = oHeader.zpxl;

	poHeader->nReplaceValue(D_K_OriginalImageFormat, Cimage_header::sGetFormatTypeName(enImage_format_BAS2000));

  (void) poHeader->nReplaceValue(Cstring (D_K_Type), Cstring("mad"));
  (void) poHeader->nReplaceValue(Cstring (D_K_Dim), Cstring("2"));
  (void) poHeader->nReplaceValue(Cstring (D_K_Size1), a10nTemp[0]);
  (void) poHeader->nReplaceValue(Cstring (D_K_Size2), a10nTemp[1]);
  (void) poHeader->nReplaceValue(Cstring (D_K_DataType), 
				 Cstring(D_K_UnsignedShortInt));
  (void) poHeader->nReplaceValue(Cstring (D_K_Compression), 
				 Cstring("BAS2000"));

  (void) poHeader->nReplaceValue(Cstring (D_K_DtdisplayOrientation), 
				 Cstring("+X+Y"));

  poHeader->nReplaceValue(Cstring (D_K_ByteOrder), 
			  Cstring (D_K_BigEndian));

  (void) poHeader->nReplaceValue(Cstring (D_K_DetectorNumber), (int)1);
  (void) poHeader->nReplaceValue(Cstring (D_K_DetectorNames), sBAS2000);
  (void) poHeader->nReplaceValue(sBAS2000 + D_K_DetectorDimensions, 2, a10nTemp);

  a10fTemp[0] = oHeader.xsiz * oHeader.xpxl;
  a10fTemp[1] = oHeader.zsiz * oHeader.zpxl;
  (void) poHeader->nReplaceValue(sBAS2000 + D_K_DetectorSize, 2, a10fTemp, 5);
  (void) poHeader->nReplaceValue(sBAS2000 + D_K_DetectorDescription,
				 Cstring ("BAS2000 conversion"));

  // The following may not be correct:

  (void) poHeader->nReplaceValue(sBAS2000 + D_K_DetectorVectors,
				 Cstring ("1 0 0 0 1 0"));

  (void) poHeader->nReplaceValue(sBAS2000 + D_K_GonioNumValues,
				 (int) 6);
  (void) poHeader->nReplaceValue(sBAS2000 + D_K_GonioUnits,
				 Cstring ("deg deg deg mm mm mm"));
  (void) poHeader->nReplaceValue(sBAS2000 + D_K_GonioNames,
                 Cstring ("RotZ RotX/Swing RotY TransX TransY TransZ/Dist"));

  (void) poHeader->nReplaceValue(sBAS2000 + D_K_GonioVectors,
		 Cstring ("0 0 1 -1 0 0 0 1 0 1 0 0 0 1 0 0 0 -1"));

  a10fTemp[0] = 0.0;
  a10fTemp[1] = 0.0;   // Two theta
  a10fTemp[2] = 0.0;
  a10fTemp[3] = 0.0;
  a10fTemp[4] = 0.0;
  a10fTemp[5] = 100.0; // xtal_to_detector
      
  (void) poHeader->nReplaceValue(sBAS2000 + D_K_GonioValues, 6, a10fTemp, 3);

  (void) poHeader->nReplaceValue(sBAS2000 + D_K_NonunfType, 
				 Cstring(D_K_NonunfStateNone));

      
  (void) poHeader->nReplaceValue(sBAS2000 + D_K_SpatialDistortionType,
				 Cstring (D_K_SpatialTypeSimple));
  a10fTemp[0] = 0.5 * oHeader.xpxl;
  a10fTemp[1] = 0.5 * oHeader.zpxl;
  a10fTemp[2] = oHeader.xsiz;
  a10fTemp[3] = oHeader.zsiz;
  (void) poHeader->nReplaceValue(sBAS2000 + D_K_SpatialDistortionInfo,
				 4, a10fTemp, 5);
      
  (void) poHeader->nReplaceValue(sBAS2000 + D_K_SpatialDistortionVectors,
				 Cstring("1 0 0 -1"));

  (void) poHeader->nReplaceValue(Cstring (D_K_SaturatedValue), 500000);
  a10fTemp[0] = 1.0;
  a10fTemp[1] = 1.54178f;

  (void) poHeader->nReplaceValue(Cstring (D_K_SourcePrefix) + 
				 D_K_Wavelength, 2, a10fTemp);
  (void) poHeader->nReplaceValue(Cstring (D_K_SourcePolarz),
				 Cstring ("0.50 1 0 0"));
  (void) poHeader->nReplaceValue(Cstring (D_K_SourceSpectralDispersion),
				 Cstring ("0.0002 0.0002"));
  (void) poHeader->nReplaceValue(Cstring (D_K_SourceCrossfire),
				 Cstring ("0.0002 0.0002 0.0 0.0"));
  (void) poHeader->nReplaceValue(Cstring (D_K_SourceSize),
				 Cstring ("0.0 0.0 0.0 0.0"));
  (void) poHeader->nReplaceValue(Cstring (D_K_SourceVectors),
				 Cstring ("0 0 1 0 1 0 1 0 0"));
  (void) poHeader->nReplaceValue(Cstring (D_K_SourceValues),
				 Cstring ("0 0"));

  (void) poHeader->nReplaceValue(Cstring (D_K_CrystalPrefix) + 
				 D_K_GonioNumValues, (int) 3);
  (void) poHeader->nReplaceValue(Cstring (D_K_CrystalPrefix) + D_K_GonioNames,
				 Cstring ("Omega Kappa Phi"));
  (void) poHeader->nReplaceValue(Cstring (D_K_CrystalPrefix) + D_K_GonioUnits,
				 Cstring ("deg deg deg"));
  a10fTemp[0] =  1.0f;
  a10fTemp[1] =  0.0f;
  a10fTemp[2] =  0.0f;

  a10fTemp[3] = -cos(50.0 * Gs_dRADIANS_PER_DEGREE);
  a10fTemp[4] =  0.0f;
  a10fTemp[5] =  sin(50.0 * Gs_dRADIANS_PER_DEGREE);

  a10fTemp[6] =  1.0f;
  a10fTemp[7] =  0.0f;
  a10fTemp[8] =  0.0f;

  (void) poHeader->nReplaceValue(Cstring (D_K_CrystalPrefix) + D_K_GonioVectors,
				 9, a10fTemp, 4);

  a10fTemp[0] = 0.0;
  a10fTemp[1] = 0.0;
  a10fTemp[2] = 0.0;
  (void) poHeader->nReplaceValue(Cstring (D_K_CrystalPrefix) + D_K_GonioValues,
				 3, a10fTemp, 3);
  (void) poHeader->nReplaceValue(Cstring (D_K_CrystalPrefix) + 
				 D_K_GonioDescription,
				 Cstring ("Single w/2-arc goniometer head"));
  a10fTemp[0] = 0.0;   // Start_phi
  a10fTemp[1] = 1.0;   // End_phi
  a10fTemp[2] = a10fTemp[1] - a10fTemp[0];
  a10fTemp[3] = 300.0; // Exposure_time
  a10fTemp[4] = 1;
  a10fTemp[5] = 0.0;
  a10fTemp[6] = 0.0;
  a10fTemp[7] = 0.0;
  a10fTemp[8] = 0.0;
  a10fTemp[9] = 0.0;
  (void) poHeader->nReplaceValue(Cstring (D_K_Rotation), 10, a10fTemp, 3);
  (void) poHeader->nReplaceValue(Cstring (D_K_ScanPrefix) + D_K_Rotation,
				 10, a10fTemp, 3);
  (void) poHeader->nReplaceValue(Cstring (D_K_RotAxisName),
				 Cstring ("Phi"));
  (void) poHeader->nReplaceValue(Cstring (D_K_ScanPrefix) + 
				 D_K_RotAxisName,
				 Cstring ("Phi"));

  a10fTemp[0] = 1.0;
  a10fTemp[1] = 0.0;
  a10fTemp[2] = 0.0;

  (void) poHeader->nReplaceValue(Cstring (D_K_RotVector),
				 3, &a10fTemp[0], 4);
  (void) poHeader->nReplaceValue(Cstring (D_K_ScanPrefix) + D_K_RotVector,
				 3, &a10fTemp[0], 4);

  // Determine scan template and seq info from the filename

  // Change last 3 digits in rsFilename to ?s

  //Cstring sTemp;

  int a3nSeqInfo[3] = {0, 1, 0};

  sTemp = sBuildScanTemplate(rsFilename, 3, &a3nSeqInfo[0]);

  (void) poHeader->nReplaceValue(Cstring (D_K_ScanTemplate), sTemp);

  (void) poHeader->nReplaceValue(Cstring (D_K_ScanSeqInfo), 3, a3nSeqInfo);

  //Add

  (void) poHeader->nReplaceValue(Cstring(D_K_Comment), 
			     Cstring(oHeader.cmnt));

  (void) poHeader->nReplaceValue(sBAS2000 + "INF_BIT", oHeader.bit);//int

  (void) poHeader->nReplaceValue(sBAS2000 + "INF_SENS", oHeader.sens);//long

  (void) poHeader->nReplaceValue(sBAS2000 + "INF_LATI", oHeader.lati);//double

  (void) poHeader->nReplaceValue(sBAS2000 + "INF_RATIO", oHeader.ratio);//double

  (void) poHeader->nReplaceValue(sBAS2000 + "INF_BSWAP", oHeader.bswap);//int

  (void) poHeader->nReplaceValue(sBAS2000 + "INF_BYTE", oHeader.byte);//int

  (void) poHeader->nReplaceValue(sBAS2000 + "INF_CCD", oHeader.ccd);//int

   return (0);

}

int 
nReadBAS2000Data(const int *pnFile, unsigned short int *puiData,
	       Cimage_header *poHeader)
{
	int i;
	int nStat = 0;
	int nDim0, nDim1;    // Dimensions of the image.
	int nNumPixels;
	int nToRead;
	char *pcTemp;
	long *plTemp;
	short *psTemp;

	_BAS2000_header oHeader;

	double	dKval, dSens, dZeroOffset, dLinearVal;
	long	*pTable, lValue;
	int		nDepth;

	int		nValue;

	nStat = poHeader->nGetValue(D_K_Size1, &nDim0);
	nStat = nStat + poHeader->nGetValue(D_K_Size2, &nDim1);

	if (0 == nStat) 
		nStat = poHeader->nGetValue(Cstring("BAS2000_INF_BIT"), &oHeader.bit);
	if (0 == nStat) 
		nStat = poHeader->nGetValue(Cstring("BAS2000_INF_SENS"), &oHeader.sens);
	if (0 == nStat) 
		nStat = poHeader->nGetValue(Cstring("BAS2000_INF_LATI"), &oHeader.lati);

	if (0 == nStat) 
		nStat = poHeader->nGetValue(Cstring("BAS2000_INF_BYTE"), &oHeader.byte);
	if (0 == nStat) 
		nStat = poHeader->nGetValue(Cstring("BAS2000_INF_BSWAP"), &oHeader.bswap);
	if (0 == nStat) 
		nStat = poHeader->nGetValue(Cstring("BAS2000_INF_RATIO"), &oHeader.ratio);
	if (0 == nStat) 
		nStat = poHeader->nGetValue(Cstring("BAS2000_INF_CCD"), &oHeader.ccd);

	if (0 == nStat) 
    {
		//Number of pixels
		nNumPixels = nDim0 * nDim1;

		// Number of bytes to read
		nToRead = nNumPixels * oHeader.byte;

		(void) dskbr((int *)pnFile, (char *)puiData, &nToRead, &nStat);
	}
	
	if (0 == nStat) 
    {
		nDepth       = (int)pow(2.0, (double)oHeader.bit);
		dSens        = 4000.0 / oHeader.sens;
		dZeroOffset = dSens * pow(10.0, oHeader.lati*-0.5);

		pTable = new long [nDepth];
		//pTable = (long *)malloc(nDepth * sizeof(long));
		//InitArray(TYPE_LONG, pTable, nDepth, 0.0);

		for(lValue = 0 ; lValue < nDepth ; lValue++) {
			dKval = oHeader.lati*((double)lValue / nDepth - 0.5);
    		dLinearVal = (dSens * pow(10.0, dKval)) - dZeroOffset;
			pTable[lValue] = (long)(dLinearVal * oHeader.ratio);
		}

		pcTemp = new char [nToRead];

		if(oHeader.byte == 2 && oHeader.bswap == 1)
			swab((const void *)puiData, (void*) pcTemp,(size_t) nToRead);
		else
			memcpy(pcTemp, (char *)puiData, nToRead);

		
		plTemp = new long [nToRead];

		if(oHeader.ccd) {
			unsigned short  *pus;
			pus = (unsigned short *)pcTemp;

   			for(i = 0 ; i < nNumPixels ; i++){
				nValue = (unsigned short)*pus;
       			plTemp[i] = nValue;
       			pus++;
       		}
		} else {
			if(oHeader.bit == 8) { 
				unsigned char   *puc;
				puc = (unsigned char *)pcTemp;

    			for(i = 0 ; i < nNumPixels ; i++){
					nValue = *puc;
           			plTemp[i] = pTable[nValue];
            		puc++;
        		}
			} else  { //16bit
				unsigned short  *pus;
				pus = (unsigned short *)pcTemp;

    			for(i = 0 ; i < nNumPixels ; i++){
					nValue = (unsigned short)*pus;
	       			plTemp[i] = pTable[nValue];
        			pus++;
        		}
			}
		}

		psTemp = new short [nToRead];

    	for(i = 0 ; i < nNumPixels ; i++){
			if(plTemp[i] < 0x8000 || oHeader.ratio == 0)
				psTemp[i] = (unsigned short)plTemp[i];
			else
				psTemp[i] = ((short)((double)plTemp[i] / oHeader.ratio)) | 0x8000;
        }

		//swab(psTemp, (char *)puiData, nToRead);
		memcpy((char *)puiData, psTemp, nToRead);

		delete [] pTable;
		delete [] pcTemp;
		delete [] plTemp;
		delete [] psTemp;

	}


	return (nStat);
}

