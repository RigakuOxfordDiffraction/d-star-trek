#ifdef __cplusplus
extern "C" {
#endif

#ifndef WIN32
#include <unistd.h>
#endif
#include <stdlib.h>

/***********************************************************************
 *
 * mar345: pck.c
 *
 * Copyright by:        Dr. Claudio Klein
 *                      X-ray Research GmbH, Hamburg
 *
 * Version:     1.0
 * Date:        16/01/1997
 *
 ***********************************************************************/

#include <stdio.h>
#include <stddef.h>
#include <math.h>
#include <ctype.h>
#include <string.h>

#define BYTE char
#define WORD short int
#define LONG int

#define PACKIDENTIFIER "\nCCP4 packed image, X: %04d, Y: %04d\n"
#define PACKBUFSIZ BUFSIZ
#define DIFFBUFSIZ 16384L
#define max(x, y) (((x) > (y)) ? (x) : (y)) 
#define min(x, y) (((x) < (y)) ? (x) : (y)) 
#define abs(x) (((x) < 0) ? (-(x)) : (x))
const LONG setbits[33] = {0x00000000L, 0x00000001L, 0x00000003L, 0x00000007L,
			  0x0000000FL, 0x0000001FL, 0x0000003FL, 0x0000007FL,
			  0x000000FFL, 0x000001FFL, 0x000003FFL, 0x000007FFL,
			  0x00000FFFL, 0x00001FFFL, 0x00003FFFL, 0x00007FFFL,
			  0x0000FFFFL, 0x0001FFFFL, 0x0003FFFFL, 0x0007FFFFL,
			  0x000FFFFFL, 0x001FFFFFL, 0x003FFFFFL, 0x007FFFFFL,
			  0x00FFFFFFL, 0x01FFFFFFL, 0x03FFFFFFL, 0x07FFFFFFL,
			  0x0FFFFFFFL, 0x1FFFFFFFL, 0x3FFFFFFFL, 0x7FFFFFFFL,
                          0xFFFFFFFFL};
#define shift_left(x, n)  (((x) & setbits[32 - (n)]) << (n))
#define shift_right(x, n) (((x) >> (n)) & setbits[32 - (n)])

/***************************************************************************/

/*
 * Function prototypes
 */

LONG 		*diff_words	(WORD *, int, int, LONG *, LONG	);
void            get_pck         (FILE *,           WORD *	);
int 		put_pck		(WORD *, int, int, int		);
static void 	unpack_word	(FILE *, int, int, WORD *	);
  /*
static int  	pack_chunk	(				);
static void	pack_longs	(				);
static int	bits		(				);
  */
static int      bits(LONG *chunk, int n);
static void     pack_longs(LONG *lng, int n, BYTE **target, int *bit, int size);
static int      pack_chunk(LONG *lng, int nmbr, int bitsize, int fdesc);

/***************************************************************************
 * Function: Put_pck
 ***************************************************************************/
int  
put_pck(WORD *img, int x, int y, int fdesc) 
{ 
int 		chunksiz, packsiz, nbits, next_nbits, tot_nbits;
LONG 		buffer[DIFFBUFSIZ];
LONG 		*diffs = buffer;
LONG 		*end = diffs - 1;
LONG 		done = 0;

	while(done < (x * y)) {
	  
	    end = diff_words(img, x, y, buffer, done);
	    done += (end - buffer) + 1;

	    diffs = buffer;
	    while(diffs <= end) {
	      
	        packsiz = 0;
	        chunksiz = 1;
	        nbits = bits(diffs, 1);
	        while(packsiz == 0) {
		  
		    if(end <= (diffs + chunksiz * 2))
			packsiz = chunksiz;
		    else {
		
			  next_nbits = bits(diffs + chunksiz, chunksiz); 
			  tot_nbits = 2 * max(nbits, next_nbits);

			  if(tot_nbits >= (nbits + next_nbits + 6))
			      packsiz = chunksiz;
			  else {
			      
				nbits = tot_nbits;
				if(chunksiz == 64)
				    packsiz = 128;
				  else
				    chunksiz *= 2;
			  }

		    }
		}

		if ( pack_chunk(diffs, packsiz, nbits / packsiz, fdesc) == 0)
			return( 0 );
		diffs += packsiz;
	     }
	}
	if ( pack_chunk(NULL, 0, 0, fdesc) == 0 );
		return( 1 );

	return( 1 );
}

/***************************************************************************
 * Function: bits
 ***************************************************************************/
static int 
bits(LONG *chunk, int n)
{ 
  int size, maxsize, i;

  for (i = 1, maxsize = abs(chunk[0]); i < n; ++i)
    maxsize = max(maxsize, abs(chunk[i]));
  if (maxsize == 0)
    size = 0;
  else if (maxsize < 8)
    size = 4 * n;
  else if (maxsize < 16)
    size = 5 * n;
  else if (maxsize < 32)
    size = 6 * n;
  else if (maxsize < 64)
    size = 7 * n;
  else if (maxsize < 128)
    size = 8 * n;
  else if (maxsize < 65536)
    size = 16 * n;
  else
    size = 32 * n;
  return(size);
}

/***************************************************************************
 * Function: pack_chunk
 ***************************************************************************/
static int pack_chunk(LONG *lng, int nmbr, int bitsize, int fdesc)
{ 
static LONG 	bitsize_encode[33] = {0, 0, 0, 0, 1, 2, 3, 4, 5, 0, 0,
                                      0, 0, 0, 0, 0, 6, 0, 0, 0, 0, 0,
                                      0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 7};
LONG 		descriptor[2], i, j;
static BYTE 	*buffer = NULL;
static BYTE 	*buffree = NULL;
static int 	bitmark;

	if( buffer == NULL )
    {
	    buffree = buffer = (BYTE *) malloc(PACKBUFSIZ);
	    bitmark = 0;
	}

	if( lng != NULL ) 
    {
	    for (i = nmbr, j = 0; i > 1; i /= 2, ++j);
	    descriptor[0] = j;
	    descriptor[1] = bitsize_encode[bitsize];
	
        if((buffree - buffer) > (PACKBUFSIZ - (130 * 4)))
        {
		    if( write(fdesc, buffer, buffree - buffer) == -1 )
			    return 0;
		    
            buffer[0] = buffree[0];
	    	buffree = buffer;
	    }
	    
        pack_longs(descriptor, 2, &buffree, &bitmark, 3);
	    pack_longs(lng, nmbr, &buffree, &bitmark, bitsize);
    }
    else 
    {
        if( write(fdesc, buffer, (buffree - buffer) + 1) == -1 )
            return 0;
        
        free((void *) buffer);
        
        buffer = NULL;
    }

	return 1;
}

/***************************************************************************
 * Function: diff_words
 ***************************************************************************/
LONG 
*diff_words(WORD *word, int x, int y, LONG *diffs, LONG done)
{ 
LONG i = 0;
LONG tot = x * y;

	if(done == 0)
	  { 
	    *diffs = word[0];
	    ++diffs;
	    ++done;
	    ++i;
	  }
	while((done <= x) && (i < DIFFBUFSIZ))
	  {
	    *diffs = word[done] - word[done - 1];
	    ++diffs;
	    ++done;
	    ++i;
	  }
	while ((done < tot) && (i < DIFFBUFSIZ))
	  {
	    *diffs = word[done] - (word[done - 1] + word[done - x + 1] +
                     word[done - x] + word[done - x - 1] + 2) / 4;
	    ++diffs;
	    ++done;
	    ++i;
	  }
	return(--diffs);
}

/***************************************************************************
 * Function: pack_longs
 ***************************************************************************/
static void
pack_longs(LONG *lng, int n, BYTE **target, int *bit, int size)
  { 
	LONG mask, window;
	int valids, i, temp;
	int temp_bit = *bit;
	BYTE *temp_target = *target;

	if (size > 0)
	  {
	    mask = setbits[size];
	    for(i = 0; i < n; ++i)
	      {
		window = lng[i] & mask;
		valids = size;
		if(temp_bit == 0)
			*temp_target = (BYTE) window;
		  else
		    {
		      temp = shift_left(window, temp_bit);
        	      *temp_target |= temp;
		    }
		 window = shift_right(window, 8 - temp_bit);
		valids = valids - (8 - temp_bit);
		if(valids < 0)
		    temp_bit += size;
		  else
		    {
		      while (valids > 0)
			{ 
			  *++temp_target = (BYTE) window;
          		  window = shift_right(window, 8);
          		  valids -= 8;
			}
        	      temp_bit = 8 + valids;
		    }
      		if(valids == 0)
      		  { 
		    temp_bit = 0;
        	    ++temp_target;
		  }
	      }
  	    *target = temp_target;
  	    *bit = (*bit + (size * n)) % 8;
	  }
}

/***************************************************************************
 * Function: get_pck
 ***************************************************************************/
void get_pck(FILE *fp, WORD *img)
{ 
    int x = 0, y = 0, i = 0;
    
    char    c = '0';
    
    char header[BUFSIZ];

	if( fp == NULL ) 
        return;
	
    rewind (fp);
	
    header[0] = '\n';
    header[1] = 0;

    while( c != EOF && ( x == 0 || y == 0) ) 
    {
        i = x = y = 0;

        c = (char)i;

        while((++i < BUFSIZ) && (c != EOF) && (c != '\n') && (x==0) && (y==0))
        {
            if((header[i] = c = (char)getc(fp)) == '\n')
          	    sscanf(header, PACKIDENTIFIER, &x, &y);
        }
    }

    unpack_word(fp, x, y, img);
}

/***************************************************************************
 * Function: unpack_word
 ***************************************************************************/
static void unpack_word(FILE *packfile, int x, int y, WORD *img)
{
int 		valids = 0, spillbits = 0, usedbits, total = x * y;
LONG 		window = 0L, spill = 0L, pixel = 0L, nextint = 0L, bitnum = 0L, pixnum = 0L;
static int 	bitdecode[8] = {0, 4, 5, 6, 7, 8, 16, 32};

    while( pixel < total )
    {
      	if (valids < 6)
        {
          if (spillbits > 0)
          {
      		window |= shift_left(spill, valids);
        	valids += spillbits;
        	spillbits = 0;
	      }
          else 
          {
      		spill = (LONG) getc(packfile);
        	spillbits = 8;
	      }
    }
    else {
    	pixnum = 1 << (window & setbits[3]);
      	window = shift_right(window, 3);
      	bitnum = bitdecode[window & setbits[3]];
      	window = shift_right(window, 3);
      	valids -= 6;
      	while ((pixnum > 0) && (pixel < total)) {
      		if (valids < bitnum) {
        		if (spillbits > 0) {
          			window |= shift_left(spill, valids);
            			if ((32 - valids) > spillbits) {
            				valids += spillbits;
              				spillbits = 0;
				}
            			else {
            				usedbits = 32 - valids;
              				spill = shift_right(spill, usedbits);
              				spillbits -= usedbits;
              				valids = 32;
				}
			}
          		else {
          			spill = (LONG) getc(packfile);
            			spillbits = 8;
			}
		}
        	else {
			--pixnum;
          		if (bitnum == 0) 
				nextint = 0;
          		else {
          			nextint = window & setbits[bitnum];
            			valids -= bitnum;
            			window = shift_right(window, bitnum);
            			if ((nextint & (1 << (bitnum - 1))) != 0)
              				nextint |= ~setbits[bitnum];}
          			if (pixel > x) {
         				img[pixel] = (WORD) (nextint +
                                      		(img[pixel-1] + img[pixel-x+1] +
                                       		img[pixel-x] + img[pixel-x-1] + 2) / 4);
            				++pixel;
				}
          			else if (pixel != 0) {
          				img[pixel] = (WORD) (img[pixel - 1] + nextint);
            				++pixel;
				}
          			else
            				img[pixel++] = (WORD) nextint;
			}
		}
	}
   }
}


#ifdef __cplusplus
} /* end extern "C" */
#endif

