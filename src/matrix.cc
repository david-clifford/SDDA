#include "matrix.h"

template <class Element> 
TMatrix<Element> :: TMatrix(const int nrow, const int ncol)
{
  _nrows = nrow;
  _ncols = ncol;
  _length = _nrows*_ncols;
  data = new Element[_length];
  _mypointer = 1;
}

template <class Element> 
TMatrix<Element> :: TMatrix(Element *xref, const int nrow, const int ncol)
{
  _nrows = nrow;
  _ncols = ncol;
  _length = _nrows*_ncols;
  data = xref;
  _mypointer = 0;
}

template <class Element> 
TMatrix<Element> :: ~TMatrix()
{
  if (_mypointer) delete [] data;
}

template <class Element> 
void TMatrix<Element> :: operator=(const Element &val)
{
  int i;
  for(i=0; i< _length; i++) data[i] = val;
}

template class TMatrix<double>;
template class TMatrix<int>;
template class TVector<double>;
template class TVector<int>;

