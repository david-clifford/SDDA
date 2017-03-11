#include "matrix.h"

void dlda(const Matrix &X, const IVector &y, 
	  Matrix &means, Vector &vars, IVector &counts);

void dqda(const Matrix &X, const IVector &y, 
	  Matrix &means, Matrix &vars, IVector &counts);

 
void pred_dlda(const Matrix &means, const Vector &vars, 
	       const IVector &counts, const Vector &priors,
	       const Matrix &X, Matrix &R);


void pred_dqda(const Matrix &means, const Matrix &vars, 
	       const IVector &counts, const Vector &priors,
	       const Matrix &X, Matrix &R);

void sdda(const Matrix &X, const IVector &y,
	  const Matrix &means, const Matrix &vars,
	  const IVector &counts, const Vector &priors,
	  IVector &S, const IVector &Start, const IVector &Never,
	  const bool prob, const bool cache, const bool xval,
	  const int jmax,
	  Vector &estore, Vector &pstore);


