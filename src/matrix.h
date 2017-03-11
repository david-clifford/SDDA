#include <math.h>

template <class Element>
class TMatrix {
 public: 
  TMatrix(const int, const int);
  virtual ~TMatrix();

  TMatrix(Element *, const int, const int);

  void operator=(const Element &);
  inline Element &operator()(const int, const int);
  inline const Element &operator()(const int, const int) const;

  inline int nrows() const {return(_nrows);}
  inline int ncols() const {return(_ncols);}
  inline int length() const {return(_length);}


 private:
  TMatrix();
  TMatrix(const TMatrix &);
  void operator=(const TMatrix &);

  Element *data;
  int _nrows, _ncols;
  int _mypointer;

  int index(const int i,const int j) const {return((j-1)*_nrows+(i-1)) ;}
  int _length;
};

template <class Element> 
Element & TMatrix<Element> ::operator()(const int i, const int j)
{
  return( data[index(i,j)] );
}

template <class Element> 
const Element & TMatrix<Element> ::operator()(const int i, const int j) const
{
  return( data[index(i,j)] );
}

typedef TMatrix<double> Matrix;
typedef TMatrix<int> IMatrix;

template <class Element>
class TVector : public TMatrix<Element> {
 public:
  explicit TVector(const int nrow) : TMatrix<Element>(nrow, 1) {}
  TVector(Element *x, const int nrow) : TMatrix<Element>(x,nrow, 1) {}

  void operator=(const Element &val) { TMatrix<Element>::operator=(val); }
  Element &operator()(const int i) { return(TMatrix<Element>::operator()(i,1));}
  const Element &operator()(const int i) const { return(TMatrix<Element>::operator()(i,1));}

 
private:
  TVector();
  TVector(const TVector &);
  void operator=(const TVector &);

};

typedef TVector<double> Vector;
typedef TVector<int> IVector;




