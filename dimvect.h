//01/09/03: changed NO_VAILDATEs to VALIDATE
//last modified 2/28/01

#ifndef __LINALG_DIMVECT_H
#define __LINALG_DIMVECT_H

#ifndef __TEMPLATE_XERROR_H
#include <template/xerror.h>
#endif

#include <iostream.h>

template <class T, int NDIM> class TDimVector
/*A vector is just an array that provides services for linear algebra
*/
{
 private:
  T* _tData;
 public:
  TDimVector(void);
  TDimVector(const T& t);
  TDimVector(const TDimVector<T, NDIM>& tV);
  TDimVector<T, NDIM>& operator = (const TDimVector<T, NDIM>& tV);
  ~TDimVector(void) {if(_tData) delete[] _tData;}
  
  T& operator [] (const int i) const
    {
#ifdef __VALIDATE__
      if(i<0 || i>=NDIM)
	throw XError("Array out-of-bounds");
#endif
      return _tData[i];
    }
  
  int GetLength(void) const {return NDIM;}
  
  TDimVector<T, NDIM> operator * (const T&) const;
  TDimVector<T, NDIM> operator / (const T&) const;
  
  TDimVector<T, NDIM> operator + (const TDimVector<T, NDIM>&) const;
  TDimVector<T, NDIM> operator - (const TDimVector<T, NDIM>&) const;
  
  T Dot(const TDimVector<T, NDIM>&) const;
  TDimVector<T, 3> Cross(const TDimVector<T, 3>&) const;
  T Cross(const TDimVector<T, 2>&) const;
  
  T Magnitude(void) const;
  
  //special constructors
  TDimVector<T, NDIM>(const T&, const T&);
};

template <class T, int NDIM> TDimVector<T, NDIM>::TDimVector(void)
  : _tData(new T[NDIM])
{
  for(int i=0; i<NDIM; i++) _tData[i]=0;
}

template <class T, int NDIM>
TDimVector<T, NDIM>::TDimVector(const T& t)
  : _tData(new T[NDIM])
{
  for(int i=0; i<NDIM; i++) _tData[i]=t;
}

template <class T, int NDIM>
TDimVector<T, NDIM>::TDimVector(const TDimVector<T, NDIM>& tV)
  : _tData(new T[NDIM])
{
  for(int i=0; i<NDIM; i++) _tData[i]=tV[i];
}

template <class T, int NDIM> TDimVector<T, NDIM>&
TDimVector<T, NDIM>::operator = (const TDimVector<T, NDIM>& tV)
{
  if(_tData) delete[] _tData;
  _tData=new T[NDIM];
  for(int i=0; i<NDIM; i++) _tData[i]=tV[i];
  return *this;
}

template <class T, int NDIM> TDimVector<T, NDIM>::TDimVector(const T& t0, const T& t1)
{
  (*this)[0]=t0;
  (*this)[1]=t1;
}

/*template <class T> TDimVector<T> operator - (const LTPoint<T>& p1, const LTPoint<T>& p2)
  {
  TDimVector<T> tV(2);
  tV[0]=p2._tX-p1._tX;
  tV[1]=p2._tY-p1._tY;
  return tV;
  }         */

template <class T, int NDIM> TDimVector<T, NDIM> TDimVector<T, NDIM>::operator * (const T& t) const
{
  TDimVector<T, NDIM> tProd;
  for(int i=0; i<NDIM; i++)
    {
      tProd[i]=(*this)[i]*t;
    }
  return tProd;
}

template <class T, int NDIM> TDimVector<T, NDIM> TDimVector<T, NDIM>::operator / (const T& t) const
{
  TDimVector<T, NDIM> tQuot;
  for(int i=0; i<NDIM; i++)
    {
      tQuot[i]=(*this)[i]/t;
    }
  return tQuot;
}

template <class T, int NDIM> TDimVector<T, NDIM> TDimVector<T, NDIM>::operator +
(const TDimVector<T, NDIM>& tV) const
{
  TDimVector<T, NDIM> tSum;
  for(int i=0; i<NDIM; i++)
    {
      tSum[i]=(*this)[i]+tV[i];
    }
  return tSum;
}

template <class T, int NDIM> TDimVector<T, NDIM> TDimVector<T, NDIM>::operator -
(const TDimVector<T, NDIM>& tV) const
{
  TDimVector<T, NDIM> tDiff;
  for(int i=0; i<NDIM; i++)
    {
      tDiff[i]=(*this)[i]+tV[i];
    }
  return tDiff;
}

template <class T, int NDIM> T TDimVector<T, NDIM>::Dot(const TDimVector<T, NDIM>& tV) const
{
  T tDot=0;
  for(int i=0; i<NDIM; i++)
    {
      tDot+=(*this)[i]*tV[i];
    }
  return tDot;
}

template <class T, int NDIM> TDimVector<T, 3> TDimVector<T, NDIM>::Cross(const TDimVector<T, 3>& tV) const
{
  TDimVector<T, 3> tCross;
  
  tCross[0]=(*this)[1]*tV[2]-(*this)[2]*tV[1];
  tCross[1]=(*this)[2]*tV[0]-(*this)[0]*tV[2];;
  tCross[2]=(*this)[0]*tV[1]-(*this)[1]*tV[0];
  
  return tCross;
}

template <class T, int NDIM> T TDimVector<T, NDIM>::Cross(const TDimVector<T, 2>& tV) const
{
  return (*this)[0]*tV[1]-(*this)[1]*tV[0];
}

template <class T, int NDIM> T TDimVector<T, NDIM>::Magnitude(void) const
{
  T tMag=0;
  for(int i=0; i<NDIM; i++)
    {
      tMag+=(*this)[i]*(*this)[i];
    }
  return sqrt(tMag);
}

/*template <class T, int NDIM> TDimVector<T, NDIM> TDimVector<T, NDIM>::Subtract
  (TDimPoint<T, NDIM> t1)
  {
  TDimVector<T, NDIM> tV;
  for(int i=0; i<NDIM; i++)
  {
  tV[i]+=tP1[i]-tP2[i];
  }
  return tV;
  }   */

#endif
