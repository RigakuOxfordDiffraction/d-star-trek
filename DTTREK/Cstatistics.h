#ifndef CSTATISTICS_H
#define CSTATISTICS_H

//# 
//#Copyright 1999 Molecular Structure Corporation  
//#               9009 New Trails Drive  
//#               The Woodlands, TX, USA  77381  
//# 
//#The contents are unpublished proprietary source  
//#code of Molecular Structure Corporation  
//# 
//#All rights reserved   
//# 
//#Cstatistics.h      Initial author: T.L.Hendrixson           Dec 1999
//#
//# This was previously named Cstat, but had to be changed because there
//# is a Cstat class defined in the dTREK library
//#

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

#include <math.h>

#include <vector>

/****************************************************************************
 *                               Definitions                                *
 ****************************************************************************/

/****************************************************************************
 *                          Structure definitions                           *
 ****************************************************************************/

/****************************************************************************
 *                               Enumerations                               *
 ****************************************************************************/

/****************************************************************************
 *                                 Typedefs                                 *
 ****************************************************************************/

/****************************************************************************
 *                                Constants                                 *
 ****************************************************************************/

/****************************************************************************
 *                            External variables                            *
 ****************************************************************************/

/****************************************************************************
 *                           Forward declarations                           *
 ****************************************************************************/

/****************************************************************************
 *                            Cstatistics class                             *
 ****************************************************************************/

// Description of class goes here

template <class T>
class Cstatistics
{

public:

   typedef typename std::vector<T>::iterator iter;

//# Constructor

      // Creates an instance of the class.
   Cstatistics(void) : m_SumX(0.0), m_SumX2(0.0)
   {}  

//# Destructor  

      // The destructor for this class.  This will free any memory that  
      // was allocated by the class. 
   virtual ~Cstatistics()
   {}

      // Adds a new entry to the list
   void Add(T d)
   { m_Values.push_back(d); m_SumX += (double)d; m_SumX2 += (double)d*(double)d; }


      // Returns the average of the elements in the list.
   double Average(void) const
   { return (0 == m_Values.size()) ? 0.0 : m_SumX/m_Values.size(); }

      // Clear the list of entries
   void Clear(void)
   { m_Values.clear(); m_SumX = m_SumX2 = 0.0; }

      // Returns a vector containing all of the images in the list
   std::vector<T> GetValues(void) const
   { return m_Values; }

      // Returns the value of the minimum element in the list
   T Min(void) const
   { T a = m_Values.front(); 
     for(iter i = m_Values.begin(); m_Values.end() != i; ++i){
        if(a > *i)
           a = *i;
     }
     return a;
   }

      // Returns the value of the maximum element in the list
   T Max(void) const
   { T a = m_Values.front(); 
     for(iter i = m_Values.begin(); m_Values.end() != i; ++i){
        if(a < *i)
           a = *i;
     }
     return a;
   }

      // Returns the number of entries in the list.
   int NumberOfEntries(void) const
   { return m_Values.size(); }

      // Returns the standard deviation of the elements in the list.
   double StandardDeviation(void) const
   { return (1 >= m_Values.size()) ? 0.0 :  sqrt((m_SumX2*m_Values.size()-
     m_SumX*m_SumX)/m_Values.size()/(m_Values.size()-1)); }

      // Returns the sum of the elements in the list.
   double Sum(void) const
   { return (0 == m_Values.size()) ? 0.0 : m_SumX; }

      //
	   // <HR>
   T operator[](unsigned long i) const
   { return (m_Values.size() <= i) ? T(0) : m_Values[i]; }

protected:

//# These variables are defined here because the derived classes need to  
//# be able to access them. 

private:

//# There is no copy constructor  
      // A non-existant copy constructor.  This is defined as a private
      // member function so that user code cannot contain the following:
      // <PRE><CODE>
      // Cstatistics instance1;
      // Cstatistics instance2(instance1);
      // </CODE></PRE>
      // in an attempt to create a copy of an instance of this class.
      // There is no reason why this must be so, it just isn't implemented yet.
   Cstatistics(const Cstatistics &Copy);  

//# There is no assignment operator  
      // A non-existant assignment operator.  This is defined as a private
      // member function so that user code cannot contain the following:
      // <PRE><CODE>
      // Cstatistics instance1,instance2;
      // instance2 = instance1;
      // </CODE></PRE>
      // in an attempt to create a copy of an instance of this class.
      // There is no reason why this must be so, it just isn't implemented yet.
   Cstatistics &operator=(const Cstatistics &AssignEqual);  


//# Member variables

   std::vector<T> m_Values;     // the list

   double m_SumX;           // sum [Xi]
   double m_SumX2;          // sum [(Xi)2]

};

#endif /* CSTATISTICS_H */
