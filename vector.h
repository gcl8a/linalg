/*
 * 2011-06-21: operator >> now indexes using operator [], and does not read the vector size
 * [this last part is a hack that needs to be undone]
 * TBD: Average() and StDev() functions
 * 2011-05-23: added statistical calculations Max(), Min(), and Range() [nb, Max and Min don't use fabs]
 * 2011-05-01: added Zero() function that I thought I already had
 * 2011-04-28: added Last() to vector
10/10/04: added Norm funcitons. CHANGED CalcMax to return infinity norm
09/15/04: added Average function
09/15/04: added Correl function, which emulates Excel CORREL (Covariance/<PI>Variances)
09/09/04: added CalcMax function [uses fabs]
06/16/03: removed XVectorError in favor of XError*/
//01/09/03: changed NO_VAILDATEs to VALIDATE
//11/11/02: fixed some ANSI incompatibilities
//1/22/02: made -1 default upper for null vector
//to avoid error when copying null vector
//12/04/01: added transpose

/*My very own vector template. Incorporates vector addition, etc.
Allows for vectors that start at something other than zero
Takes care of the memory management

Note that at the end, I define a few frequently used classes:

vector --> vector of doubles
fvector --> vector of floats
ivector --> vector of ints
*/

#ifndef __LINALG_VECTOR_H
#define __LINALG_VECTOR_H

#include <iostream>

#include <template/xerror.h>

using namespace std;

template < class T >
T MIN( T t1, T t2 ) { return t1 < t2 ? t1 : t2; }

template < class T >
T MAX( T t1, T t2 )
   { return t1 > t2 ? t1 : t2; }


//see matrix.h


template < class T >
class TMatrix;

/*
Define a vector that runs from indices _lower->_upper, inclusive.
_tData is a pointer to the block allocated in memory, but
_tIndex is offset so that it can be accessed properly
*/


template < class T >
class TVector
   {
private:
   int _lower;
   int _upper;
   T * _tData;
   T * _tIndex;

public:
   TVector( void ) : _lower( 0 ), _upper( -1 ) { _tData = _tIndex = 0; }
   TVector( int length, int lower = 0 );
   TVector( T * data, int length, int lower = 0 );

   //copy constructor and assignment routines
   TVector( const TVector < T > & );

   TVector < T > & operator = ( const TVector & );
   TVector < T > & operator = ( const T& );

   //destructor
   ~TVector( void )
      {
      if ( _tData ) delete[] _tData;
      _tData = 0;
      }

   //access to an entry is by reference to allow assignment as an lvalue
   T & operator[] ( int i ) const
      {
#ifdef __VALIDATE__
      if ( i < _lower || i > _upper )
         throw XError( "Array out-of-bounds" );
#endif
      return _tIndex[i];
      }

   int Length( void )const { return _upper - _lower + 1; }
   T Last(void) const {return (*this)[Length()-1];}
   void Zero(void)
   {
	   if(_tData)
	   {
		   for ( int i = _lower; i <= _upper; i++ ) _tIndex[i] = 0;
	   }
   }

   //checks to see if two vectors are compatible for comparison, addition, whatever
   int IsCompatible( const TVector < T > & t ) const
      { return ( _lower == t._lower && _upper == t._upper ); }

   //comparison operators
   int operator == ( const TVector < T > & ) const;
   //overload an operator for "nearly" equals
   int operator %= ( const TVector < T > & ) const;

   TVector < T > operator + ( const TVector < T > & ) const;
   TVector < T > operator - ( const TVector < T > & ) const;
   TVector < T > & operator += ( const TVector < T > & );
   TVector < T > & operator -= ( const TVector < T > & );

   TVector < T > operator * ( const T & ) const;
   TVector < T > operator / ( const T & ) const;
   TVector < T > & operator *= ( const T & );
   TVector < T > & operator /= ( const T & );

   TVector < T > DotMultiply( const TVector < T > & ) const;
   TVector < T > DotDivide( const TVector < T > & ) const;

   //TVector < T > operator / ( const TMatrix < T > & tM ) const;
   //{TVector x(Length()); tM.Solve(x, *this); return x;}

   TMatrix < T > MakeTranspose( void );
   T Dot( const TVector < T > & ) const;
   T Average( void )const;
   T Correl( const TVector < T > & ) const;

   T Range (void) const {return Max()-Min();}
   T Min( void )const;
   int IndexOfMax( void )const;
   T Max( void )const;

   T AbsMax( void )const;
   T CalcNorm( int )const;

   T CalcL1Norm( void )const { return CalcNorm( 1 ); } //temp -- make faster later
   T CalcL2Norm( void )const { return CalcNorm( 2 ); }

   ostream& Write (ostream&) const;
   istream& Read (istream&);
 
   friend class TMatrix < T >;
   };

#define dvector TVector<double>
#define fvector TVector<float>
#define ivector TVector<int>

template < class T >
TVector < T >::TVector( int length, int lower )
: _lower( lower ), _upper( length + lower - 1 )
   {
   _tData = new T[length];
   _tIndex = _tData - lower;
   for ( int i = _lower; i <= _upper; i++ ) _tIndex[i] = 0;
   }

//create a vector from an existing array, data.
//memory management of data is taken over by this object

template < class T >
TVector < T >::TVector( T * data, int length, int lower )
: _lower( lower ), _upper( length + lower - 1 )
   {
   _tData = 0;
   _tIndex = data - lower;
   }

//creates a copy

template < class T >
TVector < T >::TVector( const TVector < T > & t )
: _lower( t._lower ), _upper( t._upper )
   {
   _tData = new T[_upper - _lower + 1];
   _tIndex = _tData - _lower;

   for ( int i = _lower; i <= _upper; i++ ) _tIndex[i] = t._tIndex[i];
   }

//assignment operator

template < class T >
TVector < T > & TVector < T >::operator = ( const TVector < T > & tV )
   {
   _lower = tV._lower;
   _upper = tV._upper;

   //if we have memory already allocated, delete it
   if ( _tData ) delete[] _tData;

   //now start anew
   _tData = new T[_upper - _lower + 1];
   _tIndex = _tData - _lower;

   for ( int i = _lower; i <= _upper; i++ ) _tIndex[i] = tV._tIndex[i];
   return * this;
   }

template < class T >
TVector < T > & TVector < T >::operator = ( const T & t )
{
	for ( int i = _lower; i <= _upper; i++ )
			_tData[i] = t;

	return * this;
}


template < class T >
int TVector < T >::operator == ( const TVector < T > & tV ) const
   {
#ifdef __VALIDATE__
   if ( !IsCompatible( tV ) ) return 0;
#endif

   for ( int i = _lower; i <= _upper; i++ )
      if ( ( * this ) [i] != tV[i] ) return 0;
   return 1;
   }

#define ALLOWABLE_ERROR 1E-4

//checks for "nearly" equal, within factor of ALLOWABLE_ERROR
//scaled by the infinity norm of this vector

template < class T >
int TVector < T >::operator %= ( const TVector < T > & t ) const
//checks to see if two vectors nearly equal each other
   {
   if ( !IsCompatible( t ) ) throw XError("Incompatible vector comparison!");

   //find largest element
   T tBiggest = 0;
   for ( int i = _lower; i <= _upper; i++ )
      if ( fabs( ( * this ) [i] ) >= tBiggest ) tBiggest = fabs( ( * this ) [i] );

   //don't divide by zero, just check if they're both zero
   if ( tBiggest == 0 ) return operator == ( t );

   for ( int j = _lower; j <= _upper; j++ )
      if ( fabs( ( ( * this ) [j] - t[j] ) / tBiggest ) > ALLOWABLE_ERROR )
         return 0;

   return 1;
   }

template < class T >
TVector < T > TVector < T >::operator + ( const TVector < T > & tV ) const
   {
   TVector < T > tSum( * this );
   return tSum += tV;
   }

template < class T >
TVector < T > TVector < T >::operator - ( const TVector < T > & tV ) const
   {
#ifdef __VALIDATE__
   if ( !IsCompatible( tV ) )
       throw XError( "Invalid vector addition!" );
#endif

   TVector < T > tDiff( Length(), _lower );
   for ( int i = _lower; i <= _upper; i++ ) tDiff[i] = _tIndex[i] - tV[i];
   return tDiff;
   }

template < class T >
TVector < T > & TVector < T >::operator += ( const TVector < T > & tV )
   {
#ifdef __VALIDATE__
   if ( !IsCompatible( tV ) )
       throw XError( "Incompatible vector addition!" );
#endif

   for ( int i = _lower; i <= _upper; i++ ) _tIndex[i] += tV[i];
   return * this;
   }

template < class T >
TVector < T > & TVector < T >::operator -= ( const TVector < T > & tV )
   {
#ifdef __VALIDATE__
   if ( !IsCompatible( tV ) )
       throw XError( "Invalid vector subtraction!" );
#endif

   for ( int i = _lower; i <= _upper; i++ ) _tIndex[i] -= tV[i];
   return * this;
   }

template < class T >
TVector < T > TVector < T >::operator * ( const T & t ) const
   {
   TVector < T > tProduct( * this );
   return tProduct *= t;
   }

template < class T >
TVector < T > TVector < T >::operator / ( const T & t ) const
   {
   TVector < T > tProduct( * this );
   return tProduct /= t;
   }

template < class T >
TVector < T > & TVector < T >::operator *= ( const T & t )
   {
   for ( int i = _lower; i <= _upper; i++ ) _tIndex[i] *= t;
   return * this;
   }

template < class T >
TVector < T > & TVector < T >::operator /= ( const T & t )
   {
   for ( int i = _lower; i <= _upper; i++ ) _tIndex[i] /= t;
   return * this;
   }

template < class T >
TVector < T > TVector < T >::DotMultiply( const TVector < T > & tV ) const
   {
   TVector < T > tProduct( * this );
   for ( int i = _lower; i <= _upper; i++ ) tProduct[i] *= tV[i];
   return tProduct;
   }

template < class T >
TVector < T > TVector < T >::DotDivide( const TVector < T > & tV ) const
   {
   TVector < T > tProduct( * this );
   for ( int i = _lower; i <= _upper; i++ ) tProduct[i] /= tV[i];
   return tProduct;
   }

template < class T >
TMatrix < T > TVector < T >::MakeTranspose( void )
   {
   TMatrix < T > t( 1, 0, Length(), _lower );
   for ( int i = _lower; i <= _upper; i++ ) t[0] [i] = _tIndex[i];
   return t;
   }

template < class T >
T TVector < T >::Dot( const TVector < T > & tV ) const
   {
#ifdef __VALIDATE__
   if ( !IsCompatible( tV ) )
       throw XError( "Invalid vector dot product!" );
#endif

   T tDot = 0;
   for ( int i = _lower; i <= _upper; i++ ) tDot += _tIndex[i] * tV[i];
   return tDot;
   }

template < class T >
T TVector < T >::Average( void )const
   {
   T avg = 0;
   for ( int i = _lower; i <= _upper; i++ )
      { avg += _tIndex[i]; }

   return avg / Length();
   }

template < class T >
T TVector < T >::Correl( const TVector < T > & tV ) const
   {
   if ( !IsCompatible( tV ) )
       throw XError( "Invalid vector correlation!" );

   T cov = 0;
   T var1 = 0;
   T var2 = 0;

   T avg1 = Average();
   T avg2 = tV.Average();


   for ( int i = _lower; i <= _upper; i++ )
      {
      cov += ( _tIndex[i] - avg1 ) * ( tV[i] - avg2 );
      var1 += ( _tIndex[i] - avg1 ) * ( _tIndex[i] - avg1 );
      var2 += ( tV[i] - avg2 ) * ( tV[i] - avg2 );
      }

   return cov / sqrt( var1 * var2 );
   }

template < class T >
T TVector < T >::Min( void )const
   {
   T tMin = _tIndex[_lower];
   for ( int i = _lower; i <= _upper; i++ )
      {
      tMin = _tIndex[i] < tMin ? _tIndex[i] : tMin;
      }
   return tMin;
   }

template < class T >
int TVector < T >::IndexOfMax( void )const
   {
   int max = _lower;
   for ( int i = _lower+1; i <= _upper; i++ )
      {
      max = _tIndex[i] > _tIndex[max] ? i : max;
      }
   return max;
   }

template < class T >
T TVector < T >::Max( void )const
   {
   T tMax = 0;
   for ( int i = _lower; i <= _upper; i++ )
      {
      tMax = _tIndex[i] > tMax ? _tIndex[i] : tMax;
      }
   return tMax;
   }

template < class T >
T TVector < T >::AbsMax( void )const
   {
   T tMax = 0;
   for ( int i = _lower; i <= _upper; i++ )
      {
      T fabsi = fabs( _tIndex[i] );
      tMax = fabsi > tMax ? fabsi : tMax;
      }
   return tMax;
   }

template < class T >
T TVector < T >::CalcNorm( int order ) const
   {
   T tNorm = 0;
   if ( !order ) //zero is surrogate for infinity norm
      return AbsMax();

   for ( int i = _lower; i <= _upper; i++ ) tNorm += pow( fabs( _tIndex[i] ), order );
   return pow( tNorm, 1.0 / ( double )order );
   }

template < class T >
ostream & operator << ( ostream & stream, const TVector < T > & t )
   {
   for ( int i = 0; i < t.Length(); i++ )
      {
      if ( i ) stream << '\t';
      stream << t[i];
      }

   return stream;
   }

template < class T >
ostream & TVector < T >::Write( ostream & stream ) const
   {
   for ( int i = _lower; i <= _upper; i++ )
      {
      if ( i - _lower ) stream << '\t';
      stream << ( * this ) [i];
      }

   return stream;
   }

template < class T >
istream & operator >> ( istream & stream, TVector < T > & t )
   {
   //int len;
   //stream >> len;
   //t = TVector < T > ( len );
   for ( int i = 0; i < t.Length(); i++ )
      stream >> t[i];

   return stream;
   }

template < class T >
istream & TVector < T >::Read( istream & stream )
   {
   int len;
   stream >> len;
   * this = TVector < T > ( len );
   for ( int i = _lower; i <= _upper; i++ )
      stream >> _tIndex[i];

   return stream;
   }


#endif
