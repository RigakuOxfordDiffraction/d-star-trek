#ifndef DT_MARCCD_H
#define DT_MARCCD_H
//
// Copyright (c) 1997 Molecular Structure Corporation
//
// marccd.h        Initial author: J.W. Pflugrath           12-Feb-1998
//    This file is the header file for a MARCCD style image file routines
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

//+Include files

#include "Cstring.h"

//+Code begin

//  The following header definition comes from M. Blum via S. McSweeney

/*
   MarCCD Header Documentataion 

   from C code in frame.h  and types.h

   Documentation updated by M. Blum Fri Dec 12 10:56:49 CST 1997


   The full header, as written to the file, is a TIFF header.
   The initial 1024 bytes are a minimal TIFF header with a standard
   TIFF TAG pointing to the image data and a private TIFF TAG
   pointing to this header structure.  As written by mmx/marccd, the
   frame_header structure always begins at byte 1024 and is 3072 bytes long
   making the full header 4096 bytes.

   The meanings of the data types should be selfevident:
   (example:  UINT32 is an unsigned 32 bit integer)
   The exact C language definition is machine dependent

   Currently frames are always written:
         origin=UPPER_LEFT
         orientation=HFAST
         view_direction=FROM_SOURCE

*/

/* This number is  written into the byte_order fields in the
   native byte order of the machine writing the file */

#define MAR_LITTLE_ENDIAN   1234
#define MAR_BIG_ENDIAN      4321


/* possible orientations of frame data */
#define HFAST                   0
#define VFAST                   1

/* possible origins of frame data */
#define UPPER_LEFT              0
#define LOWER_LEFT              1
#define UPPER_RIGHT             2
#define LOWER_RIGHT             3

/* possible view directions of frame data for
   any given orientation and origin */

#define FROM_SOURCE             0
#define TOWARD_SOURCE           1

#define INT32  int
#define UINT16 unsigned short int
#define UINT32 unsigned int


typedef struct _MARCCD_header {
  char  acTIFF[1024];  // Space for the TIFF header

        /* File/header format parameters (256 bytes) */

  UINT32        header_type;      /* flag for header type  (can be used as magic number) */
  char header_name[16];           /* header name (MMX) */
  UINT32        header_major_version;     /* header_major_version (n.) */
  UINT32        header_minor_version;     /* header_minor_version (.n) */
  UINT32        header_byte_order;/* MAR_BIG_ENDIAN (Motorola,MIPS);
				     MAR_LITTLE_ENDIAN (DEC, Intel) */
  UINT32        data_byte_order;  /* MAR_BIG_ENDIAN (Motorola,MIPS); MAR_LITTLE_ENDIAN (DEC, Intel) */
  UINT32        header_size;      /* in bytes                    */
  UINT32        frame_type;       /* flag for frame type */
  UINT32        magic_number;     /* to be used as a flag -usually to indicate new file */
  UINT32        compression_type; /* type of image compression   */
  UINT32        compression1;     /* compression parameter 1 */
  UINT32        compression2;     /* compression parameter 2 */
  UINT32        compression3;     /* compression parameter 3 */
  UINT32        compression4;     /* compression parameter 4 */
  UINT32        compression5;     /* compression parameter 4 */
  UINT32        compression6;     /* compression parameter 4 */
  UINT32        nheaders;         /* total number of headers     */
  UINT32        nfast;            /* number of pixels in one line*/
  UINT32        nslow;            /* number of lines in image    */
  UINT32        depth;            /* number of bytes per pixel   */
  UINT32        record_length;    /* number of pixels between succesive rows */
  UINT32        signif_bits;      /* true depth of data, in bits */
  UINT32        data_type;        /* (signed,unsigned,float...) */
  UINT32        saturated_value;  /* value marks pixel as saturated */
  UINT32        sequence;         /* TRUE or FALSE */
  UINT32        nimages;          /* total number of images - size of each
				     is nfast*(nslow/nimages) */
  UINT32        origin;           /* corner of origin            */
  UINT32        orientation;      /* direction of fast axis   */
  UINT32        view_direction;   /* direction to view frame     */
  UINT32        overflow_location;/* FOLLOWING_HEADER,FOLLOWING_DATA */
  UINT32        over_8_bits;      /* # of pixels with counts > 255 */
  UINT32        over_16_bits;     /* # of pixels with count >65535 */
  UINT32        multiplexed;      /* multiplex flag */
  UINT32        nfastimages;      /* # of images in fast direction*/
  UINT32        nslowimages;      /* # of images in slow direction*/
  UINT32        background_applied; /* flags correction has been applied -
				       hold magic number ? */
  UINT32        bias_applied;       /* flags correction has been applied -
				       hold magic number ? */
  UINT32        flatfield_applied;  /* flags correction has been applied -
				       hold magic number ? */
  UINT32        distortion_applied; /* flags correction has been applied -
				       hold magic number ? */
  UINT32        original_header_type;     /* Header/frame type from file that frame is read from */
  UINT32        file_saved;         /* Flag that file has been saved, 
				       should be zeroed if modified */
  char reserve1[(64-40)*sizeof(INT32)-16];

  /* Data statistics (128) */
  UINT32        total_counts[2];  /* 64 bit integer range = 1.85E19*/
  UINT32        special_counts1[2];
  UINT32        special_counts2[2];
  UINT32        min;
  UINT32        max;
  UINT32        mean;
  UINT32        rms;
  UINT32        p10;
  UINT32        p90;
  UINT32        stats_uptodate;
#define MAXIMAGES 9
  UINT32        pixel_noise[MAXIMAGES];           /* 1000*base noise value (ADUs) */
  char reserve2[(32-13-MAXIMAGES)*sizeof(INT32)];

  /* More statistics (256) */

  UINT16 percentile[128];

  /* Goniostat parameters (128 bytes) */
  INT32 xtal_to_detector;         /* 1000*distance in millimeters*/
  INT32 beam_x;                   /* 1000*x beam position (pixels)*/
  INT32 beam_y;                   /* 1000*y beam position (pixels)*/
  INT32 integration_time;         /* integration time in milliseconds */
  INT32 exposure_time;            /* exposure time in milliseconds */
  INT32 readout_time;             /* readout time in milliseconds */
  INT32 nreads;                   /* number of readouts to get this image */
  INT32 start_twotheta;           /* 1000*two_theta angle */
  INT32 start_omega;              /* 1000*omega angle */
  INT32 start_chi;                /* 1000*chi angle */
  INT32 start_kapp;               /* 1000*kappa angle */
  INT32 start_phi;                /* 1000*phi angle */
  INT32 start_delta;              /* 1000*delta angle */
  INT32 start_gamma;              /* 1000*gamma angle */
  INT32 start_xtal_to_detector;   /* 1000*distance in mm (dist in um)*/
  INT32 end_twotheta;             /* 1000*two_theta angle */
  INT32 end_omega;                /* 1000*omega angle */
  INT32 end_chi;                  /* 1000*chi angle */
  INT32 end_kappa;                /* 1000*kappa angle */
  INT32 end_phi;                  /* 1000*phi angle */
  INT32 end_delta;                /* 1000*delta angle */
  INT32 end_gamma;                /* 1000*gamma angle */
  INT32 end_xtal_to_detector;     /* 1000*distance in mm (dist in um)*/
  INT32 rotation_axis;            /* active rotation axis */
  INT32 rotation_range;           /* 1000*rotation angle */
  INT32 detector_rotx;            /* 1000*rotation of detector around X */
  INT32 detector_roty;            /* 1000*rotation of detector around Y */
  INT32 detector_rotz;            /* 1000*rotation of detector around Z */

  char reserve3[(32-28)*sizeof(INT32)];

  /* Detector parameters (128 bytes) */

  INT32 detector_type;            /* detector type */
  INT32 pixelsize_x;              /* pixel size (nanometers) */
  INT32 pixelsize_y;              /* pixel size (nanometers) */
  INT32 mean_bias;                /* 1000*mean bias value */
  INT32 photons_per_100adu;       /* photons / 100 ADUs */
  INT32 measured_bias[MAXIMAGES]; /* 1000*mean bias value for each image*/
  INT32 measured_temperature[MAXIMAGES];  /* Temperature of each
					     detector in milliKelvins */
  INT32 measured_pressure[MAXIMAGES];     /* Pressure of each
					     chamber in icroTorr */
  /* Retired reserve4 when MAXIMAGES set to 9 from 16 and two
     fields removed, and temp and pressure added
     char reserve4[(32-(5+3*MAXIMAGES))*sizeof(INT32)];
  */

  /* X-ray source and optics parameters (128 bytes) */
  /* X-ray source parameters (8*4 bytes) */
  INT32 source_type;              /* (code) - target, synch. etc */
  INT32 source_dx;                        /* Optics param. - (size micron) */
  INT32 source_dy;                        /* Optics param. - (size micron) */
  INT32 source_wavelength;                /* wavelength (femtoMeters) */
  INT32 source_power;             /* (Watts) */
  INT32 source_voltage;           /* (Volts) */
  INT32 source_current;           /* (microAmps) */
  INT32 source_bias;              /* (Volts) */
  INT32 source_polarization_x;    /* () */
  INT32 source_polarization_y;    /* () */
  char reserve_source[4*sizeof(INT32)];

  /* X-ray optics_parameters (8*4 bytes) */
  INT32 optics_type;              /* Optics type (code)*/
  INT32 optics_dx;                /* Optics param. - (size micron) */
  INT32 optics_dy;                /* Optics param. - (size microns) */
  INT32 optics_wavelength;                /* Optics param. - (size microns) */
  INT32 optics_dispersion;                /* Optics param. -(*10E6) */
  INT32 optics_crossfire_x;       /* Optics param. - (microRadians) */
  INT32 optics_crossfire_y;       /* Optics param. - (microRadians) */
  INT32 optics_angle;             /* Optics param. - (monoch.
				     2theta - microradians) */
  INT32 optics_polarization_x;    /* () */
  INT32 optics_polarization_y;    /* () */
  char reserve_optics[4*sizeof(INT32)];

  char reserve5[((32-28)*sizeof(INT32))];

  /* File parameters (1024 bytes) */
  char filetitle[128];            /* Title                              */
  char filepath[128];             /* path name for data file            */
  char filename[64];              /* name of data file                  */
  char acquire_timestamp[32];     /* date and time of acquisition       */
  char header_timestamp[32];      /* date and time of header update     */
  char save_timestamp[32];        /* date and time file saved           */
  char file_comments[512];        /* comments  - can be used as desired */
  char reserve6[1024-(128+128+64+(3*32)+512)];

  /* Dataset parameters (512 bytes) */
  char dataset_comments[512];     /* comments  - can be used as
				     desired */

  char pad[3072-(256+128+256+(3*128)+1024+512)];     /* pad out to
							3072 bytes */

} tagMARCCD_header;

// Forward declaration of classes needed by the function prototypes below

class Cimage_header;

int nReadMARCCDHeader(const int* pnFile, const Cstring& rsFilename,
			      char *pcBuffer, Cimage_header *poHeader);
#endif   // DT_MARCCD_H
