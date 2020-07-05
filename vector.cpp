//create a vector of length length, starting at lower
//automatically zeros all elements

#include <linalg/vector.h>

class Test {protected: int test;};
class DTest : public Test {public: DTest(void) {tst=1;};

using namespace std;

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
int TVector < T >::operator == ( const TVector < T > & tV ) const
   {
   if ( !IsCompatible( tV ) ) return 0;
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
   if ( !IsCompatible( tV ) )
       throw XError( "Invalid vector addition!" );

   TVector < T > tDiff( Length(), _lower );
   for ( int i = _lower; i <= _upper; i++ ) tDiff[i] = _tIndex[i] - tV[i];
   return tDiff;
   }

template < class T >
TVector < T > & TVector < T >::operator += ( const TVector < T > & tV )
   {
   if ( !IsCompatible( tV ) )
       throw XError( "Incompatible vector addition!" );

   for ( int i = _lower; i <= _upper; i++ ) _tIndex[i] += tV[i];
   return * this;
   }

template < class T >
TVector < T > & TVector < T >::operator -= ( const TVector < T > & tV )
   {
   if ( !IsCompatible( tV ) )
       throw XError( "Invalid vector subtraction!" );

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
   if ( !IsCompatible( tV ) )
       throw XError( "Invalid vector dot product!" );

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
T TVector < T >::CalcMax( void )const
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
      return CalcMax();

   for ( int i = _lower; i <= _upper; i++ ) tNorm += pow( fabs( _tIndex[i] ), order );
   return pow( tNorm, 1.0 / ( double )order );
   }

/*template < class T >
ostream & operator << ( ostream & stream, const TVector < T > & t )
   {
   for ( int i = t.GetLower(); i <= t.GetUpper(); i++ )
      {
      if ( i - t.GetLower() ) stream << '\t';
      stream << t[i];
      }

   return stream;
   }*/

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
   int len;
   stream >> len;
   t = TVector < T > ( len );
   for ( int i = t._lower; i <= t._upper; i++ )
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
