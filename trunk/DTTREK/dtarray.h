//
// Copyright 1999 Molecular Structure Corporation
//                9009 New Trails Drive
//                The Woodlands, TX, USA  77381
//
// The contents are unpublished proprietary source
// code of Molecular Structure Corporation
//
// All rights reserved
//
// Cstat.cpp    Initial author: Thaddeus J Niemeyer           Dec 1999
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
// Description:
//
//
// ToDo:
//
//

#ifndef DTARRAY_DEFINED
#define DTARRAY_DEFINED

#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <limits.h>
#include <stdarg.h>
#include <ctype.h>
#include "minmax.h"

#include "dtrekdefs.h"

#include "Dtrek.h"


//////////////////////////////////
//////////////////////////////////
/***********
Apparently not used by d*TREK
template <class T_DEL> class deleter {
	T_DEL** tdp;
	public:
	int clear() { tdp=NULL;return 0; };
	deleter(T_DEL*& _tdp) {tdp=&_tdp;};
  ~deleter() { if (tdp) delete[] (*tdp); }
	};

template <class T_TYPE> T_TYPE boundR(T_TYPE& x,T_TYPE low,T_TYPE high) { if (x<low) x=low; else if (x>high) x=high; return x; };
template <class T_TYPE> T_TYPE boundNR(T_TYPE x,T_TYPE low,T_TYPE high) { if (x<low) x=low; else if (x>high) x=high; return x; };
**************/

//////////////////////////////////
//////////////////////////////////

#define TEMPLATEERROR { exit(0); }

template <class T> class staticarray {
	T* poType;
	int nAllocated;
	
public:

	T* poGet() { return poType; };
	T* poGet(int nSize) { nSetSize(nSize); return poType; };
	int nSetSize(int nSize) { if (nAllocated<nSize) { delete[] poType; poType = new T[nSize]; nAllocated = nSize; }; return 0; };

	staticarray() { nAllocated = 0; poType = NULL; };
    ~staticarray() { delete[] poType; };
};

//////////////////////////////////
//////////////////////////////////

DTREK_EXPORT void* itr_malloc(int nSize);
DTREK_EXPORT void itr_free(void* pvPointer);



template <class T_ITER> class itr
{
protected:

	T_ITER* d;
	int             maxd;                   // The number of elements allocated in array.
	int             dpt;                    // Pointer for iterating first,next loops.
	int             maxwi;                  // Maximum index written

	int increase(int);

public:
	T_ITER& operator [](int index)
    {
      //Old line of code which gives core dumps with newer compiler:
      //return d[ maxd < index + 1 ? increase(index+1) : 0, maxwi = maxwi > index ? maxwi : index, index];
      // was replace by the following which seems to work:
      // However, it would be nice to re-write the "nTemp = ...." line to be more clear, would it not?
      int nTemp;
      nTemp = maxd < index + 1 ? increase(index+1) : 0, maxwi = maxwi > index ? maxwi : index, index;
      // In the above line, nTemp is always 0.
      // debugging: if (index > maxd) printf("index, maxd are::  %d %d\n", index, maxd);      
      return d[index];
    }

	int            push(const T_ITER& x) {  (*this)[dpt+1]=(T_ITER&)x;dpt++; return 0;};
	T_ITER& pop()         { if (dpt==-1L) TEMPLATEERROR; return /* (*this) */d[dpt--];};
	T_ITER&  top() { return (*this)[dpt]; };

	int     index() const { return dpt;};
	T_ITER&  now() { return (*this)[dpt]; };
	int      isempty() const { return ((dpt==-1) || (dpt>maxwi));};
	T_ITER& 	last() { return (*this)[maxwi]; };

	int		size() const { return maxwi+1; };
	int           setsize(int newsize)	{ (*this)[newsize-1]; maxwi=newsize-1; return 0;};
	int            clear() { dpt=-1; maxwi=-1; return 0;};
	int				reset() { maxd=0;d=NULL;clear(); return 0;};
    int				dupc(char val,int sizein = -1) { (*this)[max(0,((sizein==-1)?size():sizein)-1)]; memset(d,val,sizeof(T_ITER)*((sizein==-1)?size():sizein)); return 0;};

	int            add(const T_ITER& x) { (*this)[maxwi+1]=(T_ITER&) x; return 0;};
	int                              copy(itr<T_ITER>& x);     // Beware of this function!!!!!

	int             insert(const itr<T_ITER>& x,int insertspot,int replace,int ignoreins=0);
    int             insert(const T_ITER& x,int insertspot);
	int                              remove(int removespot,int deleteamount) { return insert(*this,removespot,deleteamount,1); };
	int             append(const itr<T_ITER>& x);
	int             find(T_ITER& x,int oldindex=-1);
	int             findNR(T_ITER x,int oldindex=-1);

    int          operator==(const itr<T_ITER>& x) const { return 0; };
	itr<T_ITER>&       operator+=(const itr<T_ITER>& x) { append(x); return (*this);};
	itr<T_ITER>&       operator+(const itr<T_ITER>& x) { append(x); return (*this);};
	itr<T_ITER>&        operator+(const T_ITER& x) {  add(x); return (*this);};
    itr<T_ITER>&        operator+(T_ITER* x) { int p = max(0,maxwi); do { (*this)[p++] = (*x); } while (!((*(x++))==0)); return (*this); };

	itr<T_ITER>&		operator -() { clear(); return (*this);};
	itr<T_ITER>&       operator=(const itr<T_ITER>& t);
    itr<T_ITER>&       operator=(T_ITER* t) { clear(); return ((*this) + t); };          // Used for example with strings.

    int operator=(int t) { if (t!=0) TEMPLATEERROR; itr_free((void*) d); d=NULL; return 0;};

	itr();
  virtual ~itr();
	};

inline int operator == (itr<int>& oOp1,int nInt) { return 0; };
inline int operator == (itr<double>& oOp1,int nInt) { return 0; };


template <class T_ITER> itr<T_ITER>::itr()
{
  maxwi = -1;  // Maximum index written (i.e. actually set to a value)
  dpt   = -1;  // Pointer for iterating first, next loops
  maxd  =  0;  // Current number of elements in this iterator
  d     =  NULL; // A pointer to the first element in the array, if NULL there is no array yet

  //+2012-05-22 jwp add one line
  //increase(1);  // Always create at least 1 entry
  //-2012-05-22 jwp 
};

template <class T_ITER> itr<T_ITER>::~itr() {
    itr_free((void*) d);
	d=NULL;
	};

template <class T_ITER> int itr<T_ITER>::append(const itr<T_ITER>& t)
{
	int x,y,z;
	
    if (&t==this) 
        TEMPLATEERROR;
	
    for (x = 0, y = size(), z = t.size(); x < z; x++)
        (*this)[y+x]=((itr<T_ITER>&)t)[x];

	return 0;
}

template <class T_ITER> int itr<T_ITER>::increase(int _maxlocal) {
  T_ITER*  _d;
  int _maxd;
  int x = 0;
  _maxd = _maxlocal;
  // JWP: Increase to size _maxd (that is: do not increaseyby _maxd
  //printf("T_ITER::increase called with _maxd: %d\n\n", _maxd);
  //cout << "T_ITER::increase called with _maxd: " << _maxd << endl << flush; 

  if (maxd > _maxd) return 0;  // Already have allocated enough memory for this size array 

  _maxd = max(10L, max(maxd*2, _maxd) );  // Always have at least an array of size 10 
    
  _d = (T_ITER*) itr_malloc(sizeof(T_ITER)*_maxd); // Allocate memory for an array of size _maxd
  memset(_d, 0, sizeof(T_ITER)*(_maxd));  // Set all values of the newly allocated memory to zero 

  if ( (sizeof(T_ITER) <= 4) && (0 < maxd) && (NULL != d) )
    {        
      // Only copy old stuff if there was memory allocated for it
	memcpy(_d, d, sizeof(T_ITER) * maxd);
	x = maxd;
    } 
  else 
    {
      for (x=0; x < maxd; x++)  // Bug maybe? was here if maxd was 0 because x was not initialized? ...
	_d[x] = d[x];           // Copy old values in old location to new location 
      // When x leaves this loop it should have the proper value needed
    }
  //printf("x is:  %d\n\n", x);
  T_ITER dtemp;             // Create a zero value of this element
  dtemp = 0;
  for (; x < _maxd; x++)     // Note that the starting value of x is set previously
    {
      //printf("in loop, x is: %d\n\n", x);
      _d[x] = dtemp;           // Initialize all the remaining allocated elements to the zero value
    }


  if (!d)                 
    { 
      maxwi = -1;  // I am not sure if these should be reset to -1
      dpt   = -1;  // And this set to -1 here as well.
    };		// In case this was a "cold increase" in which case all data members are zero
  itr_free((void*) d);  // if d == NULL, itr_free should do nothing

  d    = _d;
  maxd = _maxd;
  return 1;  // 1 means size of array was increased
};


template <class T_ITER> itr<T_ITER>& itr<T_ITER>::operator=(const itr<T_ITER>& t) {
	int     ct,x;
	if (&t==this) return *this;
	dpt=-1;
	increase(t.maxd);
	for (ct=0,x=t.size();ct<x;ct++) d[ct]=t.d[ct];
	maxwi=t.maxwi;
	dpt=t.dpt;
	return *this;
	};

// This copier copies memory, not objets.
template <class T_ITER> int itr<T_ITER>::copy(itr<T_ITER>& x) {
	T_ITER* t;
	if (maxd<x.maxd) {
        itr_free((void*) d);
        t = (T_ITER*) itr_malloc(sizeof(T_ITER)*x.maxd);
		d=t;
		maxd=x.maxd;
		};
	maxwi=x.maxwi;
	dpt=x.dpt;
	memcpy(d,x.d,sizeof(T_ITER)*x.maxd);
	return 0;
	};

template <class T_ITER> int itr<T_ITER>::insert(const itr<T_ITER>& ins,int insertspot,int remove,int ignoreins) {
	 int x,y;
	 int this_size;
	 int ins_size;

	 this_size=size();
	 if (ignoreins) ins_size=0; else ins_size=ins.size();

	 increase(this_size+ins_size-remove);
	 maxwi=this_size-1+ins_size-remove;
	 if (maxwi<-1) TEMPLATEERROR;
	 if (remove>ins_size) {
		for (x=insertspot+remove;x<this_size;x++) d[x+ins_size-remove]=d[x];
		} else {
		for (x=this_size-1;x>=insertspot+remove;x--) d[x+ins_size-remove]=d[x];
		};
	 for (x=insertspot,y=0;y<ins_size;x++,y++) d[x]=((itr<T_ITER>&) ins)[y];
	 return (0);
	};

template <class T_ITER> int itr<T_ITER>::insert(const T_ITER& ins,int insertspot) {
    int x;
    int this_size;

    this_size = size();
    if (insertspot>this_size)
        return 1;
    increase(size()+1);
    for (x=this_size;x-1>=insertspot;x--)
        d[x] = d[x-1];
    maxwi = this_size;
    d[insertspot] = ins;
    return 0;
};

template <class T_ITER> int itr<T_ITER>::find(T_ITER& x, int oldindex)
{
    for(oldindex++; oldindex < size(); oldindex++) 
    {    
        if( d[oldindex] == x ) 
            return oldindex; 
    }
    
    return -1;
}

template <class T_ITER> int itr<T_ITER>::findNR(T_ITER x,int oldindex) {
 for(oldindex++;oldindex<size();oldindex++) if (d[oldindex]==x) return oldindex; return -1;
 };

typedef itr<double> doublearray;
typedef itr<int> HITR;
typedef itr<float> floatarray;


#endif
