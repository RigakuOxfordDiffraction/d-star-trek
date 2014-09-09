#ifndef DT_MARIP_H
#define DT_MARIP_H
//
// Copyright (c) 1999 Molecular Structure Corporation
//
// marip.h        Initial author: J.W. Pflugrath           21-Feb-1999
//    This file is the header file for a MARIP style image file routines
//    I have included code given to me by Dr. Claudio Klein of MAR.
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

/***********************************************************************
 *
 * mar300_header.h
 *
 * Copyright by:        Dr. Claudio Klein
 *                      X-ray Research GmbH, Hamburg
 *
 * Version:     1.0
 * Date:        05/07/1996
 *
 **********************************************************************/

#include "dtrekdefs.h"

typedef struct {
        int	pixels_x;
	int	pixels_y;
        int	lrecl;
	int	max_rec;
        int	high_pixels;
	int	high_records;
        int	counts_start;
	int	counts_end;
        int	exptime_sec;
        int	exptime_units;

        float	prog_time;
	float	r_max;
        float	r_min;
	float	p_r;
        float	p_l;
	float	p_x;
        float	p_y;
	float	centre_x;
        float	centre_y;
	float	lambda;
        float	distance;
	float	phi_start;
        float	phi_end;
	float	omega;
        float	multiplier;

	char	date[24];

} MAR300_HEADER;

/***********************************************************************
 *
 * mar345_header.h
 *
 * Copyright by:        Dr. Claudio Klein
 *                      X-ray Research GmbH, Hamburg
 *
 * Version:     1.4
 * Date:        14/05/1998
 *
 * History:
 * Version    Date    	Changes
 * ______________________________________________________________________
 *
 * 1.4        14/05/98  Element gap introduced 
 *
 ***********************************************************************/

typedef struct {

	int	byteorder;		/* Always = 1234 */
        char    version[8];		/* Program version           */
        char    program[16];		/* Program name              */

	/* Scanner specific things */
	short	scanner;		/* Scanner serial no. */
	short	size;  			/* No. of pixels in 1 dimension */
	char	format;			/* Image format */
	char	mode;			/* Exposure mode */
	int	high;			/* No. high intensity pixels */
        int     pixels;			/* No. of pixels in image */
	int	adc_A;			/* Offset from channel A of ADC */
	int	adc_B;			/* Offset from channel B of ADC */
	int	add_A;			/* ADD to channel A of ADC */
	int	add_B;			/* ADD to channel B of ADC */
	int	gap;			/* GAP position seen by controller */

        float	pixel_length;		/* Length of 1 pixel */
        float	pixel_height;		/* Height of 1 pixel */
        float   multiplier;     	/* Multiplication factor */
        float   xcen;			/* Center x of transf. image */
        float   ycen;			/* Center y of transf. image */
        float   roff;			/* Radial offset             */
        float   toff;			/* Tangential offset         */
        float   gain;			/* Gain of detector          */

	/* Experimental conditions for this image */
        float   time;			/* Exposure time in secs */
        float   dosebeg;		/* Dose at start of expose */
        float   doseend;		/* Dose at end   of expose */
        float   dosemin;		/* Min. dose during expose */
        float   dosemax;		/* Max. dose during expose */
        float   doseavg;		/* Avg. dose during expose */
        float   dosesig;		/* Sig. dose during expose */
        float   wave;  			/* Wavelength [Ang.] */
        float   dist;			/* Distance [mm] */
        float   resol;			/* Max. resolution */
        float   phibeg;			/* Starting PHI */
        float   phiend;			/* Ending   PHI */
        float   omebeg;			/* Starting Omega */
        float   omeend;			/* Ending   Omega */
        float   theta;			/* Two theta */
        float   chi;			/* Chi */
	int	phiosc;			/* Phi oscillations */
	int	omeosc;			/* Omega oscillations */
        int     dosen;  		/* No. of X-ray readings   */

	/* Generator settings */
	char	source[32];		/* Type of source */
        float   kV;  			/* Generator: kV */
        float   mA;  			/* Generator: mA */
	
	/* Monochromator */
	char	filter[32];		/* Type of monochromator */
        float   polar; 			/* Beam polarization factor */
        float   slitx; 			/* Slit width               */
        float   slity; 			/* Slit height              */

	/* Image statistics  */
	int	valmin;			/* Min. pixel value */
	int	valmax;			/* Max. pixel value */
	float	valavg;			/* Avg. pixel value */
	float	valsig;			/* Sig. pixel value */
	int	histbeg;		/* Start of histogram */
	int	histend;		/* End   of histogram */
	int	histmax;		/* Max.  of histogram */

	/* Remark             */
	char	remark[56];		/* Remark */

	/* Time of production */
	char	date[24];		/* Creation date */

} MAR345_HEADER;

// Forward declaration of classes needed by the function prototypes below

class Cimage_header;

int nReadMARIPHeader(const int* pnFile, const Cstring& rsFilename,
			      char *pcBuffer, Cimage_header *poHeader);
extern "C" {

#define WORD unsigned short int

void get_pck		(FILE *, WORD *);
}

#endif   // DT_MARIP_H
