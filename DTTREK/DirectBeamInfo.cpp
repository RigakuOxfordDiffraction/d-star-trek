//
// Copyright © 2000-2001 Molecular Structure Corporation
// Copyright © 2001-2002 Rigaku MSC
// Copyright © 2002-2005 Rigaku/MSC, Inc.
// Copyright © 2005      Rigaku Americas
//                       9009 New Trails Drive
//                       The Woodlands, TX, USA  77381
//
// The contents are unpublished proprietary source
// code of Rigaku Americas
//
// All rights reserved
//
// DirectBeamInfo.cpp   Initial author: T.L.Hendrixson               Sep 2000
//
// Description:
//
//
// ToDo:
//
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
 
/****************************************************************************
 *                              Include Files                               *
 ****************************************************************************/

#ifdef WIN32
#pragma warning(disable:4786) // so warning about debug term > 255 chars ignored
#endif

#include <sstream>
#include <iomanip>
#include <iostream>
#include <algorithm>

#include "DirectBeamInfo.h"

#include "Cstatistics.h"

#include "Cfind.h"
#include "C3Ddata.h"

#ifdef MSC_MSCSTRING_H
#include "StreamTemplate.h"
#endif

/****************************************************************************
 *                               Definitions                                *
 ****************************************************************************/

/****************************************************************************
 *                                 Typedefs                                 *
 ****************************************************************************/

/****************************************************************************
 *                          Structure definitions                           *
 ****************************************************************************/

/****************************************************************************
 *                               Enumerations                               *
 ****************************************************************************/

/****************************************************************************
 *                                Constants                                 *
 ****************************************************************************/

/****************************************************************************
 *                         Static Member Variables                          *
 ****************************************************************************/

int DirectBeamInfo::m_BoxSize = 51;
int DirectBeamInfo::m_BoxWidth = m_BoxSize/2;

/****************************************************************************
 *                             Global variables                             *
 ****************************************************************************/

/****************************************************************************
 *                            External variables                            *
 ****************************************************************************/

/****************************************************************************
 *                           Function prototypes                            *
 ****************************************************************************/

/****************************************************************************
 *                           Non-class functions                            *
 ****************************************************************************/

class comp_ptr_rank {
public:
   bool operator()(DirectBeamInfo::PeakInfo *a, DirectBeamInfo::PeakInfo *b)
   {
      return a->m_Rank < b->m_Rank;
   }
};

/****************************************************************************
 ****************************************************************************
 **                                                                        **
 **                         Public Member Functions                        **
 **                                                                        **
 ****************************************************************************
 ****************************************************************************/

/****************************************************************************
 * This is the constructor for the class.                                   *
 ****************************************************************************/
DirectBeamInfo::DirectBeamInfo(void) :
   m_WtI(0.25),
   m_WtR(0.75),
   m_PeakInfo(),
   m_pImage(NULL),
   m_Identifier(""),
   m_MaskFilename(""),
   m_OutputString(""),
   m_RankingSummary("")
{
   return;
}

/****************************************************************************
 * This is the destructor.                                                  *
 ****************************************************************************/
DirectBeamInfo::~DirectBeamInfo()
{
}

/****************************************************************************
 * This routine evaluates the direct beam image <VAR>pImage</VAR>.  The     *
 * parameters <VAR>StartLine</VAR> and <VAR>NumberOfLines</VAR> can be used *
 * to restrict the analysis of the direct beam to a smaller portion of the  *
 * image.  If these parameters are not specified, then the entire image     *
 * be examined.                                                             *
 *                                                                          *
 * During the evaluation, information about the direct beam, such as        *
 * position, size and intensity, is calculated and an output string         *
 * containing an ASCII representation of the direct beam and the            *
 * information determined is generated (but not written to any output       *
 * device).                                                                 *
 ****************************************************************************/
void DirectBeamInfo::Evaluate(Cimage *pImage,
                              const int StartLine,
                              const int NumberOfLines)
{
   m_PeakInfo = PeakInfo();
   m_OutputString = "";
   m_RankingSummary = "";

   if(NULL != pImage && pImage->bIsAvailable()){
      m_pImage = pImage;
   }
   else{
      m_pImage = NULL;
   }

   int i = (0 >= StartLine) ? 0 : StartLine;
   int j = (0 >= NumberOfLines && NULL != pImage) ? pImage->nGetDimension(1)-i : 
                                                    NumberOfLines;
   Analyze(i,j);

//   if(m_PeakInfo.m_IsUsed)
      GenerateOutputString();

   return;
}


/****************************************************************************
 * This routine places the integrated intensity, I, and the error in the    *
 * intensity, &sigma;(I), into the parameters <VAR>Intensity</VAR> and      *
 * <VAR>SigmaI</VAR>, respectively.                                         *
 ****************************************************************************/
void DirectBeamInfo::GetIntensity(double &Intensity, 
                                  double &SigmaI)
{
   Intensity = m_PeakInfo.m_Intensity;
   SigmaI = m_PeakInfo.m_SigmaI;
   return;
}

/****************************************************************************
 * This routine places the X (fastest index) and Y (slowest index)
 * positions of the most intense pixel in the direct beam in the 
 * parameters <VAR>XPosition</VAR> and <VAR>YPosition</VAR>,
 * respectively.  The values of the positions are in pixels.
 ****************************************************************************/
void DirectBeamInfo::GetMostIntensePixelPosition(long &XPosition, 
                                                 long &YPosition)
{
   XPosition = m_PeakInfo.m_MostIntensePixelX;
   YPosition = m_PeakInfo.m_MostIntensePixelY;
   return;
}

/****************************************************************************
 * This routine returns the intensity of the most intense pixel in
 * the direct beam.
 ****************************************************************************/
long DirectBeamInfo::GetMostIntensePixelIntensity(void)
{
   return m_PeakInfo.m_MostIntensePixelI;
}

/****************************************************************************
 ****************************************************************************/
long DirectBeamInfo::GetNetIntensity(void)
{
   return m_PeakInfo.m_NetIntensity;
}

/****************************************************************************
 * This routine returns the output string.                                  *
 ****************************************************************************/
std::string DirectBeamInfo::GetOutputString(void)
{
   return m_OutputString;
}

/****************************************************************************
 * This routine returns the intensity of the pixel closest to the position  *
 * of the direct beam.                                                      *
 ****************************************************************************/
long DirectBeamInfo::GetPixelIntensity(void)
{
   return m_PeakInfo.m_PixelIntensity;
}

/****************************************************************************
 ****************************************************************************/
DirectBeamInfo::PeakInfo DirectBeamInfo::GetPeakInfo(void)
{
   return m_PeakInfo;
}

/****************************************************************************
 * This routine places the X (fastest index) and Y (slowest index)          *
 * positions of the direct beam in the parameters <VAR>XPosition</VAR> and  *
 * <VAR>YPosition</VAR>, respectively.  The values of the positions are in  *
 * pixels.                                                                  *
 ****************************************************************************/
void DirectBeamInfo::GetPosition(double &XPosition, 
                                 double &YPosition)
{
   XPosition = m_PeakInfo.m_X;
   YPosition = m_PeakInfo.m_Y;
   return;
}

/****************************************************************************
 ****************************************************************************/
std::string DirectBeamInfo::GetRankingSummary(void)
{
   return m_RankingSummary;
}

/****************************************************************************
 * This routine places the width of the direct beam into the parameters     *
 * <VAR>XWidth</VAR> and <VAR>YWidth</VAR> for the X (fastest index) and Y  *
 * (slowest index), respectively.  The width is calculated using the full   *
 * width at half maximum pixel value for the peak.                          *
 ****************************************************************************/
void DirectBeamInfo::GetWidth_FWHM(double &XWidth, 
                                   double &YWidth)
{
   XWidth = m_PeakInfo.m_FWHM_X;
   YWidth = m_PeakInfo.m_FWHM_Y;
   return;
}

/****************************************************************************
 * This routine places the width of the direct beam into the parameters     *
 * <VAR>XWidth</VAR> and <VAR>YWidth</VAR> for the X (fastest index) and Y  *
 * (slowest index), respectively.  The width is calculated using the point  *
 * at which the intensity in the peak is 3&sigma; above the background.     *
 ****************************************************************************/
void DirectBeamInfo::GetWidth_3Sigma(double &XWidth, 
                                     double &YWidth)

{
   XWidth = m_PeakInfo.m_SigmaWidthX;
   YWidth = m_PeakInfo.m_SigmaWidthY;
   return;
}

/****************************************************************************
 * This routine sets the identifier used in the output string (in the first *
 * line).                                                                   *
 ****************************************************************************/
void DirectBeamInfo::SetIdentifier(const std::string &s)
{
   m_Identifier = s;
   return;
}

/****************************************************************************
 ****************************************************************************/
bool DirectBeamInfo::SetRankingWeights(const double &Intensity,
                                       const double &Resolution)
{
   bool Error = 0.0 > Intensity || 0.0 > Resolution;

   if(!Error){
      m_WtI = Intensity;
      m_WtR = Resolution;
   }

   return !Error;
}

/****************************************************************************
 * This routine sets the filename of the mask that will be used when
 * determining the direct beam position.  <VAR>Filename</VAR> is the
 * filename of the image file that will be used as the mask.  If
 * <VAR>Filename</VAR> is an empty string, then a mask will not be
 * used in the determination.
 ****************************************************************************/
void DirectBeamInfo::SetMaskFilename(const std::string &Filename)
{
   m_MaskFilename = Filename;
}

/****************************************************************************
 ****************************************************************************
 **                                                                        **
 **                       Protected Member Functions                       **
 **                                                                        **
 ****************************************************************************
 ****************************************************************************/
  
/****************************************************************************
 ****************************************************************************
 **                                                                        **
 **                        Private Member Functions                        **
 **                                                                        **
 ****************************************************************************
 ****************************************************************************/

/****************************************************************************
 * This routine analyzes the image to locate the position of the direct     *
 * beam and to generate information about the shape and intensity of the    *
 * direct beam.                                                             *
 ****************************************************************************/
void DirectBeamInfo::Analyze(const int StartLine,
                             const int ReadLines)
{

   // We need to change the image that we will be using so that we
   // make the saturated value bigger than what it actually is in the 
   // image, so if the direct beam is saturated, it will be found in
   // the 2D find of the image.

   float Saturation = 0.0;
   if(NULL != m_pImage){
      Saturation = m_pImage->fGetSatValue();
      m_pImage->vSetSatValue(2*Saturation);
   }

   // We need to change the direct beam position contained in the spatial
   // distortion info to be (0,0), so that the d*TREK spot finding routine
   // will find the direct beam if it is in the middle of the image.  Note
   // that we have to change both the SpatialDistortionKeyowrd and the 
   // SpatialBeamPosition keyword.

   std::string Prefix = GetDetectorPrefix();
   Cstring SpatialInfo("");
   Cstring SpatialBeam("");

   if(NULL != m_pImage){
      std::string Key = Prefix;
      Key += D_K_SpatialBeamPosition;
      m_pImage->m_oHeader.nGetValue(Key.c_str(),&SpatialBeam);
      if("" != SpatialBeam){
         m_pImage->m_oHeader.nReplaceValue(Key.c_str(),"0 0");
      }

      Key = Prefix;
      Key += D_K_SpatialDistortionInfo;
      m_pImage->m_oHeader.nGetValue(Key.c_str(),&SpatialInfo);
      if("" != SpatialInfo){
         double x,y,px,py;
         std::istringstream in(SpatialInfo.string());
         in >> x >> y >> px >> py;
         std::ostringstream out;
         out << "0 0 " << px << ' ' << py;
         m_pImage->m_oHeader.nReplaceValue(Key.c_str(),out.str().c_str());
      }
   }


   bool Found = LocateDirectBeam();
   if(NULL != m_pImage && Found && m_PeakInfo.m_IsUsed)
      CalculateDirectBeamSize();


   // Don't forget to change the values back to what they should be.

   if(NULL != m_pImage){

      m_pImage->vSetSatValue(Saturation);

      if("" != SpatialBeam){
         std::string Key = Prefix;
         Key += D_K_SpatialBeamPosition;
         m_pImage->m_oHeader.nReplaceValue(Key.c_str(),SpatialBeam);
      }

      if("" != SpatialInfo){
         std::string Key = Prefix;
         Key += D_K_SpatialDistortionInfo;
         m_pImage->m_oHeader.nReplaceValue(Key.c_str(),SpatialInfo);
      }
   }

   return;
}

/****************************************************************************
 ****************************************************************************/
void DirectBeamInfo::CalculateDirectBeamSize(void)
{

// Use a small portion of the image to do a full width at half maximum
// position calculation.  We will use a small region about the direct beam
// position for the calculation.

   C3Ddata *pBeamData = new C3Ddata(m_BoxSize,m_BoxSize,1,
                                    nint(m_PeakInfo.m_X)-m_BoxWidth,
                                    nint(m_PeakInfo.m_Y)-m_BoxWidth,0,
                                    e3Ddata_float);
   pBeamData->nFill2D(m_pImage,0);

// Get the 3sigma width of the spot.
// Also have to get the average background to use in the FWHM 
// determination.

   float BkgdAvg,dummy;
   float pos[3],size[3];

   pos[0] = pos[1] = m_BoxWidth;
   pos[2] = 0.5;  // since 2D spot, set center to be 1/2 a slice.

   pBeamData->nCalcGetPeakInfo(pos,&dummy,&dummy,&BkgdAvg,&dummy,size);

   m_PeakInfo.m_SigmaWidthX = size[0];
   m_PeakInfo.m_SigmaWidthY = size[1];

   CalculateFWHM(*pBeamData,BkgdAvg);

   delete pBeamData;

   return;
}

/****************************************************************************
 ****************************************************************************/
void DirectBeamInfo::CalculateFWHM(C3Ddata &BeamData,
                                   const float BkgdAvg)
{
   // Calculate the full width at half max in the x and y directions.  We
   // know that the most intense pixel is at (m_BoxWidth,m_BoxWidth), 
   // because we've defined our 3Ddata object that way.

   int n = m_BoxWidth+1;
   double *pProfile = new double [n];

   int i;
	// +ve x direction
   for(i = 0; i < n; i++)
      pProfile[i] = BeamData.fGetValue(m_BoxWidth+i,m_BoxWidth,0)-BkgdAvg;
   m_PeakInfo.m_FWHM_X = PositionOfHalfMaximum(pProfile,n);
   // -ve x direction
   for(i = 0; i < n; i++)
      pProfile[i] = BeamData.fGetValue(m_BoxWidth-i,m_BoxWidth,0)-BkgdAvg;
   m_PeakInfo.m_FWHM_X += PositionOfHalfMaximum(pProfile,n);

	// +ve y direction
   for(i = 0; i < n; i++)
      pProfile[i] = BeamData.fGetValue(m_BoxWidth,m_BoxWidth+i,0)-BkgdAvg;
   m_PeakInfo.m_FWHM_Y = PositionOfHalfMaximum(pProfile,n);
   // -ve y direction
   for(i = 0; i < n; i++)
      pProfile[i] = BeamData.fGetValue(m_BoxWidth,m_BoxWidth-i,0)-BkgdAvg;
   m_PeakInfo.m_FWHM_Y += PositionOfHalfMaximum(pProfile,n);

   delete [] pProfile;

   return;
}

/****************************************************************************
 ****************************************************************************/
void DirectBeamInfo::CalculateNetIntensity(C3Ddata &BeamData)
{
   // Do a box sum for the intensity.  Use the outer ring of pixels to
   // get the average pixel background, and then subtract the total pixel
   // background from the summed intensity.

   Cstatistics<float> bkgd;
   Cstatistics<float> box;

	for(int i = 0; i < m_BoxSize; i++){
		for(int j = 0; j < m_BoxSize; j++){
         float f = BeamData.fGetValue(i,j,0);
         box.Add(f);
			if(0 == j || m_BoxSize-1 == j || 0 == i || m_BoxSize-1 == i)
            bkgd.Add(f);
		}
	}

	m_PeakInfo.m_NetIntensity = nint(box.Sum()-bkgd.Average()*
                                    box.NumberOfEntries());

   return;
}

/****************************************************************************
 ****************************************************************************/
void DirectBeamInfo::CalculateProfileBox(int &Xmin,
                                         int &Ymin,
                                         int &Xsize,
                                         int &Ysize)
{
// Figure out the edges of the box for printing the profile

   int xsize = nint(m_PeakInfo.m_SigmaWidthX)/2+10;
   int ysize = nint(m_PeakInfo.m_SigmaWidthY)/2+10;
   int xmin = nint(m_PeakInfo.m_X-xsize);
   int ymin = nint(m_PeakInfo.m_Y-ysize);
   int xmax = nint(m_PeakInfo.m_X+xsize);
   int ymax = nint(m_PeakInfo.m_Y+ysize);


// Make sure that the entire box is in the image.  If not, adjust the
// edges so that it is.

   if(NULL != m_pImage){
      int x = m_pImage->nGetDimension(0);
      int y = m_pImage->nGetDimension(1);
      if(xmin < 0)
         xmin = 0;
      if(ymin < 0)
         ymin = 0;
      if(xmax >= x)
         xmax = x-1;
      if(ymax >= y)
         ymax = y-1;
   }

// Calculate the (possibly new) size of the box.

   xsize = xmax-xmin+1;
   ysize = ymax-ymin+1;
   if(16 > ysize)
      ysize = 16;

   Xsize = xsize;
   Ysize = ysize;
   Xmin = xmin;
   Ymin = ymin;

   return;
}

/****************************************************************************
 * This routine finds the most intense spot in the reflexion list,
 * <VAR>pRefnlist</VAR>, and stores the spot information as the
 * current direct beam information.
 ****************************************************************************/
bool DirectBeamInfo::FindMostIntenseSpot(Creflnlist &Reflnlist,
                                         Cfind &Find)
{
   int n = Reflnlist.nGetNumReflns();
   if(0 == n)
      return false;

   PeakInfo *peaks = new PeakInfo [n];
   PeakInfo **plist = new PeakInfo * [n];
   int i;
   for(i = 0; i < n; i++)
      plist[i] = peaks+i;

   double MaxIntensity = 0.0;
   GetPeakCentroids(Reflnlist,Find,peaks,MaxIntensity);

   for(i = 0; i < n; i++){
      peaks[i].m_RankI = GetIntensityRank(peaks[i],MaxIntensity);
      peaks[i].m_RankR = GetResolutionRank(peaks[i]);

      if(peaks[i].m_RankI > 0.05){
         peaks[i].m_IsUsed = true;
         peaks[i].m_Rank = m_WtI*peaks[i].m_RankI+m_WtR*peaks[i].m_RankR;
      }
      else{
         peaks[i].m_IsUsed = false;
         peaks[i].m_Rank = 0.0;
      }
   }

   SortByRank(plist,0,n);

   // Sorting places list from from low to high, so reverse it in place

   for(i = 0; i < n/2; i++){
      PeakInfo *p = plist[i];
      plist[i] = plist[n-i-1];
      plist[n-i-1] = p;
   }

   int k = 20;
   if(k > n)
      k = n;
   m_RankingSummary = RankingSummary(plist,k);

   m_PeakInfo = plist[0]->m_IsUsed ? *plist[0] : PeakInfo();

   delete [] plist;
   delete [] peaks;

   return m_PeakInfo.m_IsUsed;
}

/****************************************************************************
 ****************************************************************************/
DirectBeamInfo::PeakInfo DirectBeamInfo::GenerateBogusPeak(void)
{
   PeakInfo p;
   p.m_IsUsed = false;
   p.m_X = -1.0;//1502.3f+(rand()%100)/20.0f;
   p.m_Y = -1.0;//1499.6f+(rand()%100)/20.0f;
   p.m_PixelIntensity = 0;//32567+(rand()%1000)/200;
   p.m_Intensity = 0.0;//345678.5f+(rand()%1000)/100.0f;
   p.m_SigmaI = 0.0;//sqrt(m_PeakInfo.m_Intensity);
   p.m_MostIntensePixelX = 0;//nint(m_PeakInfo.m_X);
   p.m_MostIntensePixelY = 0;//nint(m_PeakInfo.m_Y);
   p.m_MostIntensePixelI = 0;//m_PeakInfo.m_PixelIntensity;
   p.m_NetIntensity = 0;//645678+(rand()%1000)/100;
	p.m_FWHM_X = 0;//5+(rand()%100)/100.0f;
	p.m_FWHM_Y = 0;//5+(rand()%100)/100.0f;
   p.m_SigmaWidthX = 0;//5;
   p.m_SigmaWidthY = 0;//5;
   return p;
}

/****************************************************************************
 * This routine generates the output string containing the direct beam      *
 * profile and information about the direct beam.  Note that if the image   *
 * contianing the direct beam is smaller than 7 pixels high, the profile    *
 * will probably look a little strange.                                     *
 ****************************************************************************/
void DirectBeamInfo::GenerateOutputString(void)
{

   int xsize = 16;
   int ysize = 16;
   int xmin = 0;
   int ymin = 0;
   int xorigin = 8;
   int yorigin = 8;

   if(NULL != m_pImage && m_PeakInfo.m_IsUsed){
      CalculateProfileBox(xmin,ymin,xsize,ysize);
      xorigin = nint(m_PeakInfo.m_X-xmin);
      yorigin = nint(m_PeakInfo.m_Y-ymin);
   }

// Create a data object containing just the region of interest.

   C3Ddata Data(xsize,ysize,1,xmin,ymin,0,e3Ddata_float);

   if(NULL == m_pImage || !m_PeakInfo.m_IsUsed){
   // If we don't have an image, then we must fill the data object with 
   // zeros and set the value of the origin peak.
      Data.vZero();
      Data.vSetValue(1.0f,xorigin,yorigin,0);
   }
   else
      Data.nFill2D(m_pImage,0);

// Build the output string.  The string looks like:
// 
// ...........
// ...11111...     Direct Beam Information
// ..128A821..     Beam Position(x,y), Intensity = (#,#),#
// ..18MZM81..     Most Intense Pixel Position(x,y), Intensity = (#,#),#
// ..1AZ+ZA1..     Integrated Intensity, sigma(I), I/sigma(I) = #,#,#
// ..81MZM81..     Net Intensity = #
// ..128A821..     3sigma width (x,y) = (#,#)
// ...11111...     Full width at half maximum (x,y) = (#,#)
// ...........
//
// Note that the string doesn't look like this anymore, but rather has a
// tabular form for the information.  However, the string does still have
// the ASCII representation of the profile.

   std::ostringstream oss;
   oss << std::fixed;

   for(int i = ysize-1; i >= 0; i--){
      oss << "   " << ProfileString(Data,i,xsize,xorigin,yorigin) <<
         "     " << TextInfo(i-yorigin) << std::endl;
   }

   m_OutputString = oss.str().c_str();

   return;
}

/****************************************************************************
 ****************************************************************************/
std::string DirectBeamInfo::GetDetectorPrefix(void)
{
   std::string Prefix("");

   if(NULL == m_pImage)
      return Prefix;

   // If there is more than one detector listed in the detector names
   // keyword, then we will use the first one.

   Cstring s;
   if(0 == m_pImage->m_oHeader.nGetValue(D_K_DetectorNames,&s)){
      std::string ss = s.string();
      Prefix = ss.substr(0,ss.find(' '));
   }

   return Prefix;
}

/****************************************************************************
 ****************************************************************************/
double DirectBeamInfo::GetIntensityRank(const DirectBeamInfo::PeakInfo &Peak,
                                        const double MaxIntensity)
{
   return Peak.m_Intensity/MaxIntensity;
}

/****************************************************************************
 ****************************************************************************/
void DirectBeamInfo::GetPeakCentroids(Creflnlist &Reflnlist,
                                      Cfind &Find,
                                      PeakInfo *pPeak,
                                      double &MaxIntensity)
{
   int n = Reflnlist.nGetNumReflns();
   if(0 == n)
      return;

   float (Cimage::*GetPixel)(const int nWhere0, const int nWhere1);
   GetPixel = m_pImage->prfGetPixel;
   
   MaxIntensity = 0.0;
   
   PeakInfo      oPeakInfoSav;
   const   int  cnTestBoxSize_beg = 20;
   const   int  cnTestBoxSize_incr = 20;
   const   int  cnTestBoxSize_end = 300;
   
   const double cfConvergenceRatioFactor = 1.05; // "100% plus" 
   const double cfConvergenceSigmaFactor = 3.0;  // number of sigmas

   Crefln*      pRefln = NULL;
    
   bool         bFoundAtLeastOnePositiveIntensity = false;
   
   // Process all found "peaks"
   for(int i = 0; i < n; i++)
   {
      pRefln = Reflnlist.poGetRefln(i);

      pPeak[i].m_MostIntensePixelX = nint(static_cast<double>(pRefln->fGetField(Reflnlist.m_nFI_fObsPx0)));
      pPeak[i].m_MostIntensePixelY = nint(static_cast<double>(pRefln->fGetField(Reflnlist.m_nFI_fObsPx1)));
      pPeak[i].m_MostIntensePixelI = nint(static_cast<double>(pRefln->fGetIntensity()));

      // Try a range of box sizes until (a) we get at least one positive intensity and 
      //                                (b) the intensity values converge or 
      //                                (c) the range of box sizes is exhausted
      bFoundAtLeastOnePositiveIntensity = false;
      for(int j = cnTestBoxSize_beg; j <= cnTestBoxSize_end; j+= cnTestBoxSize_incr)
      {
         Find.m_a2nSpotWindow[0] = Find.m_a2nSpotWindow[1] = j;
         Find.nFindCentroid2D(*m_pImage,pRefln);
         
         pPeak[i].m_nBoxSize = j;

         if( pRefln->fGetIntensity() > 0.0 )
         {           
            pPeak[i].m_Intensity = pRefln->fGetIntensity();
            pPeak[i].m_SigmaI = pRefln->fGetSigmaI();

            pPeak[i].m_X = pRefln->fGetField(Reflnlist.m_nFI_fObsPx0);
            pPeak[i].m_Y = pRefln->fGetField(Reflnlist.m_nFI_fObsPx1);
            pPeak[i].m_PixelIntensity = nint((m_pImage->*GetPixel)(nint(pPeak[i].m_X),nint(pPeak[i].m_Y)));

            if( bFoundAtLeastOnePositiveIntensity )
            {
                if( pPeak[i].m_Intensity < oPeakInfoSav.m_Intensity * cfConvergenceRatioFactor ||
                    pPeak[i].m_Intensity - oPeakInfoSav.m_Intensity < cfConvergenceSigmaFactor * 
                                                                      sqrt(pPeak[i].m_SigmaI*pPeak[i].m_SigmaI + 
                                                                      oPeakInfoSav.m_SigmaI*oPeakInfoSav.m_SigmaI) )
                {
                    // convergence achieved
                    pPeak[i] = oPeakInfoSav; // use the previous result
                    break;
                }
            }
            
            oPeakInfoSav = pPeak[i]; // save the current result
            
            bFoundAtLeastOnePositiveIntensity = true;
         }
      }

      if( pPeak[i].m_Intensity > MaxIntensity )
         MaxIntensity = pPeak[i].m_Intensity;
   }

   return;
}

/****************************************************************************
 ****************************************************************************/
double DirectBeamInfo::GetResolutionRank(const DirectBeamInfo::PeakInfo &Peak)
{
   int X = m_pImage->nGetDimension(0)/2;
   int Y = m_pImage->nGetDimension(1)/2;
   double R2 = X*X+Y*Y;
   double x = Peak.m_X-X;
   double y = Peak.m_Y-Y;
   return 1.0-sqrt((x*x+y*y)/R2);
}

/****************************************************************************
 ****************************************************************************/
bool DirectBeamInfo::LocateDirectBeam(void)
{
   bool Found = true;

   // If the image is bogus/non-existent, generate not-so-bogus values.

   if(NULL == m_pImage){
      m_PeakInfo = GenerateBogusPeak();
   }
   else{

      // Create a reflexion list and find all of the spots in the image

      Cnonunf *pMask = NULL;
      if("" != m_MaskFilename){
         Cnonunf *p = new Cnonunf(m_MaskFilename.c_str());
         if(p->bIsAvailable())
            pMask = p;
         else
            delete p;
      }

      Creflnlist        Reflnlist;
      Cfind     Find(m_pImage->m_oHeader,&Reflnlist,pMask);
      Find.m_nPeakFilter = 9;  // RB: use a larger than the default (6) filter size
      Find.nFind2D(*m_pImage,0.0,NULL,NULL,false);
      
      Reflnlist.nWrite("temp.ref");


      Found = (0 != Reflnlist.nGetNumReflns());
      if(Found)
         Found = FindMostIntenseSpot(Reflnlist,Find);

      if(Found){
         C3Ddata BeamData(m_BoxSize,m_BoxSize,1,
                          nint(m_PeakInfo.m_X)-m_BoxWidth,
                          nint(m_PeakInfo.m_Y)-m_BoxWidth,0,e3Ddata_float);
         BeamData.nFill2D(m_pImage,0);
         CalculateNetIntensity(BeamData);
      }

      delete pMask;
   }

   return Found;
}

/****************************************************************************
 * This routine returns a pair containing the lines that make up the
 * edges of the image.  Lines that onyl contain 0 are assumed to be 
 * paddng and not part of the actual image.  The first element of the 
 * pair returned is the top line in the actual image, and the second
 * element of the pair returned is the bottom line in the actual image.
 ****************************************************************************/
std::pair<int,int> DirectBeamInfo::LocateImageEdges(Cimage &Image)
{
   std::pair<int,int> a;
   a.first = a.second = 0;

   // Find edges of the image.

   int x = Image.nGetDimension(0);
   int y = Image.nGetDimension(1)-1;

   float (Cimage::*GetPixel)(const int, const int);
   GetPixel = Image.prfGetPixel;

   bool AllZeros = true;
   for(NULL; 0 <= y && AllZeros; y--){
      for(int i = 0; AllZeros && i < x; i++)
         AllZeros = (0.0 == (Image.*GetPixel)(i,y));
   }
   a.first = y;

   AllZeros = false;
   for(NULL; 0 <= y && !AllZeros; y--){
      AllZeros = true;
      for(int i = 0; AllZeros && i < x; i++)
         AllZeros = (0.0 == (Image.*GetPixel)(i,y));
   }
   a.second = y;

   return a;
}

/****************************************************************************
 * This routine pads a portion of the edges of the image, if they 
 * contain 0 values, with an average background value.  This helps
 * locate the direct beam when it is near the edge of the image.
 ****************************************************************************/
void DirectBeamInfo::PadImageEdges(void)
{
   if(NULL == m_pImage)
      return;

   // Find edges of the image.

   int maxx = m_pImage->nGetDimension(0);
   int maxy = m_pImage->nGetDimension(1);
   int upper = maxy/2;
   int lower = upper;
   int x = maxx/2;

   float (Cimage::*GetPixel)(const int nWhere0, const int nWhere1);
   GetPixel = m_pImage->prfGetPixel;

   for(NULL; upper < maxy; upper++){
      if(0.0 == (m_pImage->*GetPixel)(x,upper))
         break;
   }

   for(NULL; 0 <= lower; lower--){
      if(0.0 == (m_pImage->*GetPixel)(x,lower))
         break;
   }

   // Fill the a band of rows adjacent to the actual image with a 
   // constant, non-zero value so that the background calculation 
   // won't be hosed.  We will use the average value of the last row
   // the image.

   if(upper < maxx){

      int i;
      int j = upper-1;
      Cstatistics<float> stat;
      for(i = 0; i < maxx; i++)
         stat.Add((m_pImage->*GetPixel)(i,j));
      float value = stat.Average();

      int end = upper+100;
      for(i = upper; i < maxx && i < end; i++){
         for(j = 0; j < maxx; j++)
            m_pImage->nSetPixel(j,i,value);
      }
   }

   if(lower >= 0){
      int i;
      int j = lower+1;
      Cstatistics<float> stat;
      for(i = 0; i < maxx; i++)
         stat.Add((m_pImage->*GetPixel)(i,j));
      float value = stat.Average();

      int end = lower-100;
      for(i = lower; 0 <= i && end <= i; i--){
         for(j = 0; j < maxx; j++)
            m_pImage->nSetPixel(j,i,value);
      }
   }

   return;
}

/****************************************************************************
 * This routine determines the position of the half-height of one side of   *
 * a peak profile of length <VAR>Length</VAR>.  It is assumed that the      *
 * maximum of the peak is the first array element of <VAR>pProfile</VAR>.   *
 * The routine will return the postion of the half height relative to the   *
 * first array element, and will interpolate the position if required.      *
 ****************************************************************************/
double DirectBeamInfo::PositionOfHalfMaximum(const double *pProfile,
                                             const int Length)
{
   double HalfMax = pProfile[0]/2.0;

   // Find first point in profile below half maximum

   int i;
   for(i = 1; i < Length && pProfile[i] > HalfMax; i++)
      NULL;

   double Position;
   if(i == Length)
      Position = i-1;
   else if(pProfile[i] == HalfMax){
      Position = i;
   }
   else{

      // Find first point in profile above half maximum

      int j;
      for(j = i-1; j > 0 && pProfile[j] < HalfMax; j--)
         NULL;

      // Interpolate position of half maximum

      Position = j+(i-j)*(HalfMax-pProfile[j])/(pProfile[i]-pProfile[j]);
   }

   return Position;
}

/****************************************************************************
 ****************************************************************************/
std::string DirectBeamInfo::ProfileString(C3Ddata &Data,
                                          const int Line,
                                          const int Xsize,
                                          const int Xcenter,
                                          const int Ycenter)
{
   static char text[] = {'.','1','2','3','4','5','6','7','8','9','A','B','C',
                         'D','E','F','G','H','I','J','K','L','M','N','O','P',
                         'Q','R','S','T','U','V','W','X','Y','Z','+'};
   std::ostringstream oss;
   oss << std::fixed;

   for(int j = 0; j < Xsize; j++){
      if(Xcenter == j && Ycenter == Line)
         oss << text[sizeof(text)-1];
      else{
         double OutputFactor = Data.fGetValue(Xcenter,Ycenter,0)/36.0;
         if(0.0 == OutputFactor)
            OutputFactor = 1.0;
         int n = nint(Data.fGetValue(j,Line,0)/OutputFactor);
         if(35 < n)
            n = 35;
         else if(0 > n)
            n = 0;
         oss << text[n];
      }
   }

   return oss.str();
}

/****************************************************************************
 ****************************************************************************/
std::string DirectBeamInfo::RankingSummary(DirectBeamInfo::PeakInfo **a,
                                           const int N)
{
   std::ostringstream os;

// ....v....1....v....2....v....3....v....4....v....5....v....6....v....7....v....8
// peak   X    Y     I         X      Y        I       Rank = wt x Ri + wt x Rr
// #### #### #### ###### -> ####.# ####.# ##########.# #.## (#.## #.## #.## #.##)

   if(0.0 >= a[0]->m_Rank){
      os << std::endl << "No peaks meet minimum ranking criteria." << 
            std::endl << std::endl;
   }
   else{
      os << std::endl << std::fixed << std::setfill(' ') << std::right;
      os << "peak   X    Y     I         X      Y         I       Rank = wt x Ri + wt x Rr" << std::endl;
      for(int i = 0; i < N && 0.0 < a[i]->m_Rank; i++){
         os << std::setw(4) << i+1 << ' ' << 
               std::setw(4) << a[i]->m_MostIntensePixelX << ' ' <<
               std::setw(4) << a[i]->m_MostIntensePixelY << ' ' <<
               std::setw(6) << a[i]->m_MostIntensePixelI << ' ' <<
               " -> " <<
               std::setw(6) << std::setprecision(1) << a[i]->m_X << ' ' <<
               std::setw(6) << std::setprecision(1) << a[i]->m_Y << ' ' <<
               std::setw(12) << std::setprecision(1) << a[i]->m_Intensity << ' ' <<
               std::setw(4) << std::setprecision(2) << a[i]->m_Rank << " (" <<
               std::setw(4) << std::setprecision(2) << m_WtI << ' ' <<
               std::setw(4) << std::setprecision(2) << a[i]->m_RankI << ' ' <<
               std::setw(4) << std::setprecision(2) << m_WtR << ' ' <<
               std::setw(4) << std::setprecision(2) << a[i]->m_RankR << ')' <<
               std::endl;
      }
      os << "Ri = Rank(intensity)" << std::endl <<
            "Rr = Rank(radial distance from center of image)" << std::endl <<
            std::endl;
   }

   return os.str();
}

/****************************************************************************
 ****************************************************************************/
void DirectBeamInfo::SortByRank(DirectBeamInfo::PeakInfo **a,
                                const int lower,
                                const int upper)
{

   std::sort(a+lower,a+upper,comp_ptr_rank());
/*
   int h;

   for(h = 1; h <= (upper-lower)/9; h += 3*h+1)
      NULL;
   for(NULL; h > 0; h /= 3){
      for(int i = lower+h; i <= upper; i++){
         int j = i;
         PeakInfo *v = a[i];
         while(j >= lower+h && v->m_Rank < a[j-h]->m_Rank){
            a[j] = a[j-h];
            j -= h;
         }
         a[j] = v;
      }
   }
*/
   return;
}

/****************************************************************************
 ****************************************************************************/
std::string DirectBeamInfo::TextInfo(const int line)
{
   std::ostringstream oss;

   oss << std::fixed;
   oss << std::setfill(' ');
   oss << std::right;

   if(8 == line){
      oss << "Direct Beam Information";
      if("" != m_Identifier)
         oss << " for " << m_Identifier;
   }
   else if(6 == line)
      oss << "  Pixel Info  |    x        y        I";
   else if(5 == line)
      oss << "--------------+---------------------------";
   else if(4 == line){
      oss << "  Direct beam | ";
      if(m_PeakInfo.m_IsUsed){
         oss << std::setw(7) << std::setprecision(2) <<
                m_PeakInfo.m_X << "  " << 
                std::setw(7) << std::setprecision(2) <<
                m_PeakInfo.m_Y << "  " << 
                std::setw(7) <<
                m_PeakInfo.m_PixelIntensity;
      }
      else{
         oss << "  n/a      n/a      n/a";
      }
   }
   else if(3 == line){
      oss << "Highest pixel | ";
      if(m_PeakInfo.m_IsUsed){
         oss << std::setw(4) <<
                m_PeakInfo.m_MostIntensePixelX << "     " << 
                std::setw(4) <<
                m_PeakInfo.m_MostIntensePixelY << "     " << 
                std::setw(7) <<
                m_PeakInfo.m_MostIntensePixelI;
      }
      else{
         oss << "  n/a      n/a      n/a";
      }
   }
   else if(1 == line)
      oss << " Intensity |     I       sigmaI  I/sigmaI";
   else if(0 == line)
      oss << "-----------+------------------------------";
   else if(-1 == line){
      oss << "Integrated | ";
      if(m_PeakInfo.m_IsUsed){
         oss << std::setw(9) << 
                nint(m_PeakInfo.m_Intensity) << "  " <<
                std::setw(7) << std::setprecision(1) <<
                m_PeakInfo.m_SigmaI << "  " << 
                std::setw(7) << std::setprecision(1);
         double d = (0.0 != m_PeakInfo.m_SigmaI) ? 
                     m_PeakInfo.m_Intensity/m_PeakInfo.m_SigmaI : 0.0;
         oss << d;
      }
      else{
         oss << "   n/a       n/a       n/a";
      }
   }
   else if(-2 == line){
      double sigma = sqrt(static_cast<double>(m_PeakInfo.m_NetIntensity));
      oss << "       Net | ";
      if(m_PeakInfo.m_IsUsed){
         oss << std::setw(9) << 
                m_PeakInfo.m_NetIntensity << "  "; // <<
                //std::setw(7) << std::setprecision(1) <<
                //sigma << "  " << 
                //std::setw(7) << std::setprecision(1) <<
                //static_cast<double>(m_PeakInfo.m_NetIntensity)/sigma;
      }
      else{
         oss << "   n/a       n/a       n/a";
      }
   }
   else if(-4 == line)
      oss << "Spot width |   x       y";
   else if(-5 == line)
      oss << "-----------+----------------";
   else if(-6 == line){
      oss << "    3sigma | ";
      if(m_PeakInfo.m_IsUsed){
         oss << std::setw(6) << std::setprecision(2) <<
                m_PeakInfo.m_SigmaWidthX << "  " << 
                std::setw(6) << std::setprecision(2) <<
                m_PeakInfo.m_SigmaWidthY;
      }
      else{
         oss <<              " n/a     n/a";
      }
   }
   else if(-7 == line){
      oss << "      FWHM | ";
      if(m_PeakInfo.m_IsUsed){
         oss << std::setw(6) << std::setprecision(2) <<
                m_PeakInfo.m_FWHM_X << "  " << 
                std::setw(6) << std::setprecision(2) <<
                m_PeakInfo.m_FWHM_Y;
      }
      else{
         oss <<              " n/a     n/a";
      }
   }

   return oss.str();
}

/****************************************************************************
 ****************************************************************************
 **                                                                        **
 **                               Operators                                **
 **                                                                        **
 ****************************************************************************
 ****************************************************************************/

/****************************************************************************
 ****************************************************************************/
DirectBeamInfo &DirectBeamInfo::operator<<(std::ostream &o)
{
   o << m_OutputString;
   return *this;
}
