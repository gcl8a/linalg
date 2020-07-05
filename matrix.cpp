#include <linalg/matrix.h>

using namespace std;

template < class T >
TMatrix < T >::TMatrix( int m, int n, int ml, int nl )
: _mLower( ml ), _nLower( nl ), _mUpper( ml + m - 1 ), _nUpper( nl + n - 1 )
/*makes a matrix and initializes to zero
*/
{
	_pptData = new T * [m];
	_pptIndeces = new T * [m];

	_tMatrix = _pptIndeces - ml;
	for ( int i = 0; i < m; i++ )
	{
		_pptData[i] = new T[n];
		_pptIndeces[i] = & _pptData[i] [-nl];
	}

	Zero();
}

template < class T >
TMatrix < T >::TMatrix( const TMatrix & tMatrix )
: _mLower( tMatrix._mLower ), _nLower( tMatrix._nLower ),
_mUpper( tMatrix._mUpper ), _nUpper( tMatrix._nUpper )
{
	int m = _mUpper - _mLower + 1;
	int n = _nUpper - _nLower + 1;

	_pptData = new T * [m];
	_pptIndeces = new T * [m];

	_tMatrix = _pptIndeces - _mLower;
	for ( int i = 0; i < m; i++ )
	{
		_pptData[i] = new T[n];
		_pptIndeces[i] = & _pptData[i] [-_nLower];
		for ( int j = 0; j < n; j++ ) _pptData[i] [j] = tMatrix._pptData[i] [j];
	}
}

template < class T >
TMatrix < T >::TMatrix( const TVector < T > & tV )
: _mLower( tV._lower ), _nLower( 0 ),
_mUpper( tV._upper ), _nUpper( 0 )
{
	int m = _mUpper - _mLower + 1;
	int n = _nUpper - _nLower + 1;

	_pptData = new T * [m];
	_pptIndeces = new T * [m];

	_tMatrix = _pptIndeces - _mLower;
	for ( int i = 0; i < m; i++ )
	{
		_pptData[i] = new T[n];
		_pptIndeces[i] = & _pptData[i] [-_nLower];
		_pptData[i] [0] = tV._tData[i];
	}
}

template < class T >
TMatrix < T > & TMatrix < T >::operator = ( const TMatrix < T > & tMatrix )
{
	if ( _pptData )
	{
		for ( int i = 0; i <= _mUpper - _mLower; i++ )
		{ delete[] _pptData[i]; }

		delete[] _pptData;
		delete[] _pptIndeces;
	}

	_mLower = tMatrix._mLower;
	_nLower = tMatrix._nLower;
	_mUpper = tMatrix._mUpper;
	_nUpper = tMatrix._nUpper;

	int m = _mUpper - _mLower + 1;
	int n = _nUpper - _nLower + 1;

	_pptData = new T * [m];
	_pptIndeces = new T * [m];

	_tMatrix = _pptIndeces - _mLower;
	for ( int i = 0; i < m; i++ )
	{
		_pptData[i] = new T[n];
		_pptIndeces[i] = & _pptData[i] [-_nLower];
		for ( int j = 0; j < n; j++ ) _pptData[i] [j] = tMatrix._pptData[i] [j];
	}

	return * this;
}

//added 4/13/01

template < class T >
TMatrix < T > & TMatrix < T >::operator = ( const T & t )
{
	for ( int i = _mLower; i <= _mUpper; i++ )
		for ( int j = _nLower; j <= _nUpper; j++ )
			_tMatrix[i] [j] = t;

	return * this;
}

template < class T >
TMatrix < T >::~TMatrix( void )
{
	if ( _pptData )
	{
		for ( int i = 0; i <= _mUpper - _mLower; i++ )
		{ delete[] _pptData[i]; }

		delete[] _pptData;
		delete[] _pptIndeces;
	}
}

template < class T >
int TMatrix < T >::IsSymmetric( void )const
{
	if ( !IsSquare() ) throw XError("Non-square matrix!");
	if ( !BaseEqual() ) throw XError("Non-symmetric matrix!");
	for ( int i = _mLower; i <= _mUpper; i++ )
		for ( int j = _nLower; j <= _nUpper; j++ )
			if ( _tMatrix[i] [j] != _tMatrix[j] [i] ) return 0;
	return 1;
}

template < class T >
int TMatrix < T >::operator == ( const TMatrix < T > & tM ) const
{
	if ( !IsCompatibleEqual( tM ) ) throw XError("Incompatible matrix equation!");
	for ( int i = _mLower; i <= _mUpper; i++ )
		for ( int j = _nLower; j <= _nUpper; j++ )
			if ( _tMatrix[i] [j] != tM[i] [j] ) return 0;

	return 1;
}

template < class T >
int TMatrix < T >::Zero( void )
{
	for ( int i = _mLower; i <= _mUpper; i++ )
	{
		for ( int j = _nLower; j <= _nUpper; j++ )
		{
#ifdef __VALIDATE__
			if ( i < _mLower || i > _mUpper || j < _nLower || j > _nUpper )
				throw XError( "Illegal index" );
#endif
			_tMatrix[i] [j] = 0;
		}
	}
	return 1;
}

template < class T >
int TMatrix < T >::ExchangeRows( int i1, int i2 )
{
	for ( int j = _nLower; j <= _nUpper; j++ )
	{
		T temp = _tMatrix[i1] [j];
		_tMatrix[i1] [j] = _tMatrix[i2] [j];
		_tMatrix[i2] [j] = temp;
	}
	return 1;
}

template < class T >
int TMatrix < T >::ExchangeColumns( int j1, int j2 )
{
	for ( int i = _mLower; i <= _mUpper; i++ )
	{
		T temp = _tMatrix[i] [j1];
		_tMatrix[i] [j1] = _tMatrix[i] [j2];
		_tMatrix[i] [j2] = temp;
	}
	return 1;
}

template < class T >
T TMatrix < T >::Max( void )const
{
	T d = 0;
	for ( int i = _mLower; i <= _mUpper; i++ )
		for ( int j = _nLower; j <= _nUpper; j++ )
			if ( fabs( _tMatrix[i] [j] ) > fabs( d ) ) d = _tMatrix[i] [j];
	return d;
}

template < class T >
TMatrix < T > TMatrix < T >::operator * ( const TMatrix < T > & tM )
const
{
	if ( !IsCompatibleMult( tM ) )
	      throw XError( "Error multiplying matrices!" );

	TMatrix < T > tSolution( _mUpper - _mLower + 1, tM._nUpper - tM._nLower + 1,
		_mLower, tM._nLower );
	for ( int i = tSolution._mLower; i <= tSolution._mUpper; i++ )
	{
		for ( int j = tSolution._nLower; j <= tSolution._nUpper; j++ )
		{
			tSolution[i] [j] = 0;
			for ( int k = _nLower; k <= _nUpper; k++ )
			{ tSolution[i] [j] += _tMatrix[i] [k] * tM[k] [j]; }
		}
	}

	return tSolution;
}

template < class T >
TVector < T > TMatrix < T >::operator * ( const TVector < T > & tV )
const
{
	if ( !IsCompatible( tV ) )
		throw XError( "Error multiplying matrix!" );

	TVector < T > tSol( _mUpper - _mLower + 1, _mLower );
	for ( int i = _mLower; i <= _mUpper; i++ )
	{
		tSol[i] = 0;
		for ( int k = _nLower; k <= _nUpper; k++ )
		{ tSol[i] += _tMatrix[i] [k] * tV[k]; }
	}

	return tSol;
}

template < class T >
TMatrix < T > TMatrix < T >::operator * ( const T & t ) const
{
	TMatrix < T > sol( * this );
	return sol *= t;
}

template < class T >
TMatrix < T > & TMatrix < T >::operator *= ( const T & t )
{
	for ( int i = _mLower; i <= _mUpper; i++ )
		for ( int j = _nLower; j <= _nUpper; j++ )
			_tMatrix[i] [j] *= t;
	return * this;
}

template < class T >
TMatrix < T > TMatrix < T >::operator / ( const T & t ) const
{
	TMatrix < T > sol( * this );
	return sol /= t;
}

template < class T >
TMatrix < T > & TMatrix < T >::operator /= ( const T & t )
{
	for ( int i = _mLower; i <= _mUpper; i++ )
		for ( int j = _nLower; j <= _nUpper; j++ )
			_tMatrix[i] [j] /= t;
	return * this;
}

template < class T >
TMatrix < T > TMatrix < T >::operator + ( const TMatrix < T > & tM ) const
{
	TMatrix < T > tSolution( * this );
	return tSolution += tM;
}

template < class T >
TMatrix < T > & TMatrix < T >::operator += ( const TMatrix < T > & tM )
{
	if ( !IsCompatibleEqual( tM ) )
		throw XError( "Error adding matricies!" );

	for ( int i = _mLower; i <= _mUpper; i++ )
		for ( int j = _nLower; j <= _nUpper; j++ )
			_tMatrix[i] [j] += tM[i] [j];
	return * this;
}

template < class T >
TMatrix < T > TMatrix < T >::operator - ( const TMatrix < T > & tM ) const
{
	TMatrix < T > tSolution( * this );
	return tSolution -= tM;
}

template < class T >
TMatrix < T > & TMatrix < T >::operator -= ( const TMatrix < T > & tM )
{
	if ( !IsCompatibleEqual( tM ) )
		throw XError( "Error subtracting matricies!" );

	for ( int i = _mLower; i <= _mUpper; i++ )
		for ( int j = _nLower; j <= _nUpper; j++ )
			_tMatrix[i] [j] -= tM[i] [j];
	return * this;
}

template < class T >
TMatrix < T > TMatrix < T >::DotMult( const TMatrix < T > & tM ) const
{
	if ( !IsCompatibleEqual( tM ) )
		throw XError( "Error multiplying matrices!" );

	TMatrix < T > tSolution( * this );
	for ( int i = _mLower; i <= _mUpper; i++ )
		for ( int j = _nLower; j <= _nUpper; j++ )
			tSolution[i] [j] *= tM[i] [j];

	return tSolution;
}

template < class T >
TMatrix < T > TMatrix < T >::DotDivide( const TMatrix < T > & tM ) const
{
	if ( !IsCompatibleEqual( tM ) )
		throw XError( "Error dividing matrices!" );

	TMatrix < T > tSolution( * this );
	for ( int i = _mLower; i <= _mUpper; i++ )
		for ( int j = _nLower; j <= _nUpper; j++ )
			tSolution[i] [j] /= tM[i] [j];

	return tSolution;
}

template < class T >
int TMatrix < T >::Solve( TVector < T > & tUnknowns,
						 const TVector < T > & tRHS ) const
{
	TMatrix < T > a( * this );
	a.Invert();
	tUnknowns = a * tRHS;
#ifdef __VALIDATE__
	if ( !( ( * this ) * tUnknowns %= tRHS ) )
		throw XError( "Error solving matrix!" );
#endif
	return 1;
}
/*
template < class T >
TVector < T > & TMatrix < T >::SolveBandDiag( TVector < T > & tUnknowns,
{
	if ( !IsCompatibleSolve( rhs ) )
		throw XError( "Error solving equation!" );

	//n.b. requires square and ml=nl, mu=nu
	int lower = _mLower;
	int upper = _mUpper;

	for ( int i = lower; i < upper; i++ )
	{
		int kmax = i + m1 > upper ? upper : i + m1;
		for ( int k = i + 1; k <= kmax; k++ )
		{
			T factor = _tMatrix[k] [i] / _tMatrix[i] [i];
			int jmax = i + m2 > upper ? upper : i + m2;
			for ( int j = i; j <= jmax; j++ )
			{ _tMatrix[k] [j] -= _tMatrix[i] [j] * factor; }
			rhs[k] -= rhs[i] * factor;
		}
	}

	//now backsubstitute
	for ( i = upper; i > lower; i-- )
	{
		int kmin = i - m2 < lower ? lower : i - m2;
		for ( int k = i - 1; k >= kmin; k-- )
		{
			T factor = _tMatrix[k] [i] / _tMatrix[i] [i];
			_tMatrix[k] [i] = 0;
			rhs[k] -= rhs[i] * factor;
		}
	}

	//now solve for unknowns
	tUnknowns = TVector < T > ( _nUpper - _nLower + 1, _nLower );
	for ( i = lower; i <= upper; i++ )
	{ tUnknowns[i] = rhs[i] / _tMatrix[i] [i]; }

	return tUnknowns;
}
 */
template < class T >
TMatrix < T > TMatrix < T >::FindInverse( void )const
{
	TMatrix < T > t = * this;
	return t.Invert();
}

template < class T >
TMatrix < T > & TMatrix < T >::Invert( void )
{
	if ( !IsSquare() )  throw XError("Non-square matrix!");

	//for speed and convenience, work in the 0-based matrix frame
	int _M = _mUpper - _mLower;
	int _N = _nUpper - _nLower;

	TVector < T > tPivot( _M + 1 ); //tracks pivots
	int iRow = 0;
	int iCol = 0;

	TVector < int > tRowIndex( _M + 1 ); //for keeping track of row operations
	TVector < int > tColIndex( _M + 1 );

	for ( int i = 0; i <= _M; i++ )
	{
		T tBiggest = 0;
		for ( int j = 0; j <= _M; j++ )
		{
			if ( tPivot[j] != 1 )
			{
				for ( int k = 0; k <= _M; k++ )
				{
					if ( tPivot[k] == 0 )
					{
						if ( fabs( _pptData[j] [k] ) >= tBiggest )
						{
							tBiggest = fabs( _pptData[j] [k] );
							iRow = j;
							iCol = k;
						}
					}
					else if ( tPivot[k] > 1 )
						throw XError("Singular matrix!"); ;
				}
			}
		}
		++( tPivot[iCol] );
		if ( iRow != iCol )
		{ ExchangeRows( iRow + _mLower, iCol + _mLower ); }

		tRowIndex[i] = iRow;
		tColIndex[i] = iCol;

		if ( _pptData[iCol] [iCol] == 0 )
						throw XError("Singular matrix!"); ;

		T tPivotInverse = 1.0 / _pptData[iCol] [iCol];
		_pptData[iCol] [iCol] = 1.0;

		for ( int jj = 0; jj <= _N; jj++ )
			_pptData[iCol] [jj] *= tPivotInverse;

		for ( int l = 0; l <= _M; l++ )
		{
			if ( l != iCol )
			{
				T tDummy = _pptData[l] [iCol];
				_pptData[l] [iCol] = 0.0;
				for ( int jj = 0; jj <= _N; jj++ )
					_pptData[l] [jj] -= _pptData[iCol] [jj] * tDummy;
			}
		}
	}
	for ( int m = _M; m >= 0; m-- )
	{
		if ( tRowIndex[m] != tColIndex[m] )
		{ ExchangeColumns( tRowIndex[m] + _nLower, tColIndex[m] + _nLower ); }
	}

	return * this;
}

template < class T >
TMatrix < T > TMatrix < T >::MakeTranspose( void )const
{
	TMatrix < T > tTranspose( _nUpper - _nLower + 1, _mUpper - _mLower + 1, _nLower, _mLower );
	for ( int i = _mLower; i <= _mUpper; i++ )
	{
		for ( int j = _nLower; j <= _nUpper; j++ )
		{ tTranspose[j] [i] = _tMatrix[i] [j]; }
	}
	return tTranspose;
}

template < class T >
int TMatrix < T >::CalcEigens( TMatrix < T > & eigVals, TMatrix < T > & eigVects )
{
	if ( !IsSquare() ) throw XError( "Non-square matrix!" );
	if ( !BaseZero() ) throw XError( "Eig analysis only allowed for base 0 matrices" );

	eigVects = * this;
	eigVects.Zero();

	eigVals = * this;
	eigVals.Zero();

	int N = _nUpper;
	TVector < T > b( N + 1 );
	TVector < T > z( N + 1 );

	for ( int i = 0; i <= N; i++ )
	{
		eigVects[i] [i] = 1.;
		eigVals[i] [i] = _tMatrix[i] [i];

		b[i] = _tMatrix[i] [i];
		z[i] = 0;
	}

	for ( int k = 0; k <= 50; k++ )
	{
		T sm = 0;
		for ( int i = 0; i <= N; i++ )
		{
			for ( int j = i + 1; j <= N; j++ )
			{ sm += fabs( _tMatrix[i] [j] ); }
		}

		if ( sm == 0 )
		{
			//we're done--let's sort things first
			for ( int i = 0; i < N; i++ )
			{
				int k = i;
				T p = eigVals[i] [i];
				for ( int j = i + 1; j <= N; j++ )
				{
					if ( eigVals[j] [j] >= p )
					{
						k = j;
						p = eigVals[j] [j];
					}
				}
				if ( k != i )
				{
					eigVals[k] [k] = eigVals[i] [i];
					eigVals[i] [i] = p;
					for ( int j = 0; j <= N; j++ )
					{
						T v = eigVects[j] [k];
						eigVects[j] [k] = eigVects[j] [i];
						eigVects[j] [i] = v;
					}
				}
			}

			return 1;
		}

		T thresh = ( k < 4 ) ? 0.2 * sm / ( N * N ) : 0.0;

		for ( int i = 0; i <= N; i++ )
		{
			for ( int j = i + 1; j <= N; j++ )
			{
				T g = 100.0 * fabs( _tMatrix[i] [j] );
				if ( k > 4 && ( float )( fabs( eigVals[i] [i] ) + g ) == ( float )fabs( eigVals[i] [i] )
					&& ( float )( fabs( eigVals[j] [j] ) + g ) == ( float )fabs( eigVals[j] [j] ) )
					_tMatrix[i] [j] = 0.0;

				else if ( fabs( _tMatrix[i] [j] ) > thresh )
				{
					T h = eigVals[j] [j] - eigVals[i] [i];
					T t;
					if ( ( float )( fabs( h ) + g ) == ( float )fabs( h ) ) t = _tMatrix[i] [j] / h;
					else
					{
						T theta = 0.5 * h / _tMatrix[i] [j];
						t = 1.0 / ( fabs( theta ) + sqrt( 1.0 + theta * theta ) );
						if ( theta < 0 ) t = -t;
					}

					T c = 1.0 / sqrt( 1 + t * t );
					T s = t * c;
					T tau = s / ( 1.0 + c );
					h = t * _tMatrix[i] [j];
					z[i] -= h;
					z[j] += h;

					eigVals[i] [i] -= h;
					eigVals[j] [j] += h;
					_tMatrix[i] [j] = 0.0;

					for ( int l = 0; l < i; l++ )
					{
						T g = _tMatrix[l] [i];
						T h = _tMatrix[l] [j];
						_tMatrix[l] [i] = g - s * ( h + g * tau );
						_tMatrix[l] [j] = h + s * ( g - h * tau );
					}

					for ( int l = i + 1; l < j; l++ )
					{
						T g = _tMatrix[i] [l];
						T h = _tMatrix[l] [j];
						_tMatrix[i] [l] = g - s * ( h + g * tau );
						_tMatrix[l] [j] = h + s * ( g - h * tau );
					}

					for ( int l = j + 1; l <= N; l++ )
					{
						T g = _tMatrix[i] [l];
						T h = _tMatrix[j] [l];
						_tMatrix[i] [l] = g - s * ( h + g * tau );
						_tMatrix[j] [l] = h + s * ( g - h * tau );
					}

					for ( int l = 0; l <= N; l++ )
					{
						T g = eigVects[l] [i];
						T h = eigVects[l] [j];
						eigVects[l] [i] = g - s * ( h + g * tau );
						eigVects[l] [j] = h + s * ( g - h * tau );
					}
				}
			}
		}

		for ( int i = 0; i <= N; i++ )
		{
			b[i] += z[i];
			eigVals[i] [i] = b[i];
			z[i] = 0.0;
		}
	}

	return 0;
}

template < class T >
TMatrix < T > TMatrix < T >::absolute( void ) const
{
 TMatrix<T> result = *this;
 for ( int i = _mLower; i <= _mUpper; i++ )
    for ( int j = _nLower; j <= _nUpper; j++ )
       if(result[i] [j]<0) result[i][j]=- result[i] [j];

 return result;
}


template < class T >
TVector < T > TMatrix < T >::GetColumn( int i ) const
{
	if ( i < _nLower || i > _nUpper ) throw XError( "Illegal address in matrix" );
	else
	{
		TVector < T > tV( _mUpper + 1, _mLower );
		for ( int j = _mLower; j <= _mUpper; j++ ) tV[j] = _tMatrix[j] [i];
		return tV;
	}
}

template < class T >
int TMatrix < T >::SetColumn( int i, const TVector < T > & tV )
{
	if ( i < _nLower || i > _nUpper ) throw XError( "Illegal address in matrix" );
	if ( !IsEqualRows( tV ) ) throw XError( "Illegal address in matrix" );
	else
		for ( int j = _mLower; j <= _mUpper; j++ ) _tMatrix[j] [i] = tV[j];

	return 1;
}

template < class T >
TVector < T > TMatrix < T >::GetRow( int i ) const
{
	if ( i < _mLower || i > _mUpper ) throw XError( "Illegal address in matrix" );
	else
	{
		TVector < T > tV( _nUpper + 1, _nLower );
		for ( int j = _nLower; j <= _nUpper; j++ ) tV[j] = _tMatrix[i] [j];
		return tV;
	}
}

template < class T >
int TMatrix < T >::SetRow( int i, const TVector < T > & tV )
{
	if ( i < _mLower || i > _mUpper ) throw XError( "Illegal address in matrix" );
	if ( !IsEqualColumns( tV ) ) throw XError( "Illegal address in matrix" );
	else
		for ( int j = _nLower; j <= _nUpper; j++ ) _tMatrix[i] [j] = tV[j];

	return 1;
}
/*
template <class T> int TMatrix<T>::CalcMultiStats(TVector<T>& mean, TMatrix<T>& covariance) const
//calculates the mean and covariance of [columns] number of samples of [rows] dimensions
{
int N=CountColumns();
int M=CountRows();

//find the mean
mean=TVector<T>(M);
for(int i=0; i<N; i++)
{
mean+=GetColumn(i);
}
mean/=N;

TMatrix<T> cov(M, N);
//find the covariance
for(i=0; i<N; i++)
{
//subtract the mean
cov.SetColumn(i, GetColumn(i)-mean);
}

covariance=cov*cov.MakeTranspose()*(1.0/(double)(N-1));

return N;
}      */

template < class T >
TMatrix < T > TMatrix < T >::Eye( int size )
{
	TMatrix < T > tM( size, size );
	tM.Zero();
	for ( int i = 0; i < size; i++ ) tM[i] [i] = 1;
	return tM;
}

template < class T >
ostream & TMatrix < T >::Write( ostream & stream ) const
{
	for ( int i = _mLower; i <= _mUpper; i++ )
	{
		for ( int j = _nLower; j <= _nUpper; j++ )
		{
			if ( j != _nLower ) stream << ", ";
			stream.precision( 8 );
			stream << _tMatrix[i] [j];
		}
		stream << endl;
	}
	return stream;
}

/*
template < class T >
ostream & operator << ( ostream & stream, const TMatrix < T > & tm )
{
for ( int i = tm._mLower; i <= tm._mUpper; i++ )
{
for ( int j = tm._nLower; j <= tm._nUpper; j++ )
{
if ( j != tm._nLower ) stream << ", ";
stream.precision( 8 );
stream << tm[i] [j];
}
stream << endl;
}
return stream;
}*/


class SampleMatrix : public dmatrix
	/*SampleMatrix assumes that the samples vectors are in the columns. i.e., there
	are <rows> number of samples
	*/
{
protected:
	dvector mean;
	dmatrix covariance;

public:
	SampleMatrix( void ) { }

	SampleMatrix( int vec_len, int samples )
		: dmatrix( vec_len, samples ), mean( vec_len ), covariance( vec_len, vec_len ) { }

	SampleMatrix( const SampleMatrix & sm ) : dmatrix( sm ), mean( sm.mean ),
		covariance( sm.covariance )
	{ }

	SampleMatrix & operator = ( const SampleMatrix & sm )
	{
		dmatrix::operator = ( sm );
		mean = sm.mean;
		covariance = sm.covariance;
		return * this;
	}

	~SampleMatrix( void ) { }

	/*		int CalcMultiStats(dvector&, dmatrix&);

	};  */

	double CalcMahalanobis( const dvector & t ) const
	{
		dvector mean;
		dmatrix cov;
		CalcMultiStats( mean, cov );
		double mahal = ( mean - t ).Dot( cov.FindInverse() * ( mean - t ) );

		return mahal;
	}

	int SampleMatrix::CalcMultiStats( dvector & mean, dmatrix & covariance ) const
		//calculates the mean and covariance of [columns] number of samples of [rows] dimensions
	{
		int N = CountColumns();
		int M = CountRows();

		//find the mean
		mean = dvector( M );
		for ( int i = 0; i < N; i++ )
		{ mean += GetColumn( i ); }
		mean /= N;

		dmatrix cov( M, N );
		//find the covariance
		for ( int i = 0; i < N; i++ )
		{
			//subtract the mean
			cov.SetColumn( i, GetColumn( i ) - mean );
		}

		covariance = cov * cov.MakeTranspose() * ( 1.0 / ( double )( N - 1 ) );

		return N;
	}
};

