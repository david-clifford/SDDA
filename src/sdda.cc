#include "matrix.h"
#include "extra.h"

static void ConvToProb(Matrix &P) 
{
  int i,j;
  double rmax, rsum;
  //
  // Takes a matrix of log likelihoods (+- const) and converts each
  // row to a vector of probabilities, being careful to avoid
  // overflow.
  //
  for(i=1; i<= P.nrows(); i++) {
    rmax = P(i,1);
    for(j=2; j<= P.ncols(); j++) if(P(i,j) > rmax) rmax = P(i,j);
    rsum =0.0;
    for(j=1; j<= P.ncols(); j++) {
      P(i,j) = exp(P(i,j)-rmax);
      rsum += P(i,j);
    }
    for(j=1; j<= P.ncols(); j++) P(i,j) /= rsum;
  }
}

void dlda(const Matrix &X, const IVector &y, 
	  Matrix &means, Vector &vars, IVector &counts) 
{
  int n = X.nrows(), p = X.ncols(), G = counts.length();
  int i,j,k;
  //
  // Basic diagonal linear discriminant analysis calculates means
  // and counts for each variable by class and a pooled variance for
  // each variable
  //
  

  means = 0.0; vars = 0.0;
  counts = 0;
    
  for(i=1; i<=n; i++) counts(y(i))++;

  for(j=1; j<=p; j++) {
    for(i=1; i<=n; i++) means(y(i),j) += X(i,j);
    for(k=1; k<=G; k++) means(k,j) /= counts(k);
      
    for(i=1; i<=n; i++) vars(j) += sqr(X(i,j)-means(y(i),j));
    vars(j) /= (n-G);
  }
}

void dqda(const Matrix &X, const IVector &y, 
	  Matrix &means, Matrix &vars, IVector &counts) 
{
  int n = X.nrows(), p = X.ncols(), G = counts.length();

  int i,j,k;
  //
  // Basic diagonal linear discriminant analysis calculates means,
  // variances and counts for each variable by class
  //

  means = 0.0; vars = 0.0;
  counts = 0;
    
  for(i=1; i<=n; i++) counts(y(i))++;

  for(j=1; j<=p; j++) {
    for(i=1; i<=n; i++) means(y(i),j) += X(i,j);
    for(k=1; k<=G; k++) means(k,j) /= counts(k);
      
    for(i=1; i<=n; i++) vars(y(i),j) += sqr(X(i,j)-means(y(i),j));
    for(k=1; k<=G; k++) vars(k,j) /= (counts(k)-1);
  }
}

void pred_dlda(const Matrix &means, const Vector &vars, 
	       const IVector &counts, const Vector &priors,
	       const Matrix &X, Matrix &R)
{ 
  int n = X.nrows(), p = X.ncols(), G = counts.length();

  int i,j,k;
  double ptmp;
  //
  // compute predicted class probabilities for dlda
  //
    
  for(k=1; k<=G; k++) {
    ptmp = -2*log(priors(k));
    for(i=1; i<=n; i++) {
      R(i,k) = ptmp;
      for(j=1; j<=p; j++) R(i,k) += sqr(X(i,j)-means(k,j))/vars(j);
      R(i,k) = -0.5*R(i,k);
    }
  }
    
  ConvToProb(R);
}	       

void pred_dqda(const Matrix &means, const Matrix &vars, 
	       const IVector &counts, const Vector &priors,
	       const Matrix &X, Matrix &R)
{ 
  int n = X.nrows(), p = X.ncols(), G = counts.length();

  int i,j,k;
  double ptmp;
  //
  // compute predicted class probabilities for dqda
  //

    
  for(k=1; k<=G; k++) {
    ptmp = -2*log(priors(k));
    for(i=1; i<=n; i++) {
      R(i,k) = ptmp;
      for(j=1; j<=p; j++) R(i,k) += log(vars(k,j)) + sqr(X(i,j)-means(k,j))/vars(k,j);
      R(i,k) = -0.5*R(i,k);
    }
  }
    
  ConvToProb(R);
}



class Robject {
  // 
  // A helper object the computes (and caches?) the i-th obs, j-th
  // variable, k-th class contribution to the x-validated log
  // likelihood
  //
public :
  Robject(const Matrix &X0, const IVector &y0, 
	  const Matrix &means0, const Matrix &vars0, const IVector &counts0,
	  const bool cache, const bool xval);
  virtual ~Robject() { if(cached) delete RddaMatrix;}
  double value(const int i, const int j, const int k);

private : 
  Robject();
  Robject(const Robject &);
  void operator=(const Robject &);

  double RddaCompute(const int i, const int j, const int k,
		     const int h, const double fac);

  const Matrix &X; 
  const IVector &y; 
  const Matrix &means;
  const Matrix &vars;
  const IVector &counts;
  double *RddaMatrix;
  int n, p, G;
  int nG, nG1;
  bool LDA, cached, xvalid;
};

Robject::Robject(const Matrix &X0, const IVector &y0, 
		 const Matrix &means0, const Matrix &vars0, const IVector &counts0,
		 const bool cache, const bool xval) 
  : X(X0), y(y0), means(means0), vars(vars0), counts(counts0),
    n(X0.nrows()), p(X0.ncols()), G(counts0.length())
{
  // Constructor

  nG = n-G;
  nG1 = n-G-1;

  if(vars.nrows()==1) LDA = true;
  else LDA = false;

  if(xval) xvalid = true;
  else xvalid = false;

  if(cache) {
    int i,j,k,h, pntr=0;
    double fac;
    RddaMatrix = new double[n*p*G];

    for(i = 1; i<=n; i++) {
      h = y(i);
      fac = ((double) counts(h))/(counts(h)-1);
      for(j=1; j<=p; j++) for(k=1; k<=G; k++) 
	RddaMatrix[pntr++] = RddaCompute(i,j,k,h,fac);
    }
    cached = true;
  } else {
    cached = false;
  }
}

double Robject::value(const int i, const int j, const int k) 
{
  // Interrogator

  if(cached) {
    return(RddaMatrix[((i-1)*p+(j-1))*G+k-1]);
  } else {
    int h = y(i);
    double fac = ((double) counts(h))/(counts(h)-1);
    return(RddaCompute(i,j,k,h,fac));
  }
  
}

double Robject::RddaCompute(const int i, const int j, const int k,
			    const int h, const double fac)
{
  double Rtmp;
  if(xvalid) {

    // Primary compute function.  Computes the individual terms (i-th
    // obs, j-th var, k-th class) to the x-validated log likelihood

    if(LDA) {
      // LDA
      Rtmp = (nG1)*sqr(X(i,j)-means(k,j))/((nG)*vars(1,j)-fac*sqr(X(i,j)-means(h,j)));
      if(h==k) Rtmp = fac*fac*Rtmp;

    } else {
      if ( k!=h) {
	Rtmp = log(vars(k,j)) + sqr(X(i,j)-means(k,j))/vars(k,j);
      } else {
	double vtmp;
	vtmp = ((counts(h)-1)*vars(h,j) - fac*sqr(X(i,j)-means(h,j)))/(counts(h)-2);
	Rtmp = log(vtmp) + fac*fac*sqr(X(i,j)-means(h,j))/vtmp;    
      }
    }
  } else {
    // Non x-validated code
    if(LDA) {
      Rtmp = sqr(X(i,j)-means(k,j))/vars(1,j);
    } else {
      Rtmp = log(vars(k,j)) + sqr(X(i,j)-means(k,j))/vars(k,j);
    }
  
  }
  return(Rtmp);
    
}
   
static Robject *Rdda;

  

static void BestAdd(const IVector &y, const IVector &S, const Matrix &T, 
		    int &jres, int &err)
{
  int i,j,k;
  int etmp, kmin;
  double pmin, ptmp;

  // Find the best variable to add using an x-validated error rate
  // criterion

  err = T.nrows();
  jres = 0;
  for(j=1; j<= S.length(); j++) if(!S(j)) {
    etmp = 0;
    for(i=1; i<= T.nrows(); i++) {
      pmin = T(i,1) + Rdda->value(i,j,1);
      kmin = 1;
      for(k=2; k<= T.ncols(); k++) {
	ptmp = T(i,k) + Rdda->value(i,j,k);
	if(ptmp<pmin) {
	  pmin=ptmp; kmin=k;
	}
      }
      if(kmin != y(i)) etmp++;
    }
    if (etmp < err) {
      err = etmp; jres=j;
    }
  }
}

static void BestAddProb(const IVector &y,
			const IVector &S, const Matrix &T, 
			int &jres, int &err, double &minsum)
{
  int i,j,k;
  double cfac, sump;
  int kmin, etmp;
  double pmin, ptmp;
  Vector Prob(T.ncols());

  // Find the best variable to add using an x-validated likelihood
  // criterion. Also compute the x-validated error rate.

  jres = 0;
  for(j=1; j<= S.length(); j++) if(!S(j)) {
    sump = 0.0;
    etmp = 0;
    for(i=1; i<= T.nrows(); i++) {
      pmin = T(i,1) + Rdda->value(i,j,1);
      kmin = 1;
      Prob(1) = pmin;
      for(k=2; k<= T.ncols(); k++) {
	ptmp = T(i,k) + Rdda->value(i,j,k);
	if(ptmp<pmin) {
	  pmin=ptmp; kmin=k;
	}
	Prob(k) = ptmp;
      }
      if(kmin != y(i)) etmp++;
      cfac = 0.0;
      for(k=1; k<= T.ncols(); k++) cfac += exp(-0.5*(Prob(k)-pmin));
      sump += Prob(y(i))-pmin + 2*log(cfac);
    }

    if (jres==0 || sump<minsum) {
	jres = j;
	minsum = sump;
	err = etmp;
    }
  }
}

static void AddIn(IVector &S, Matrix &T, 
		  const int j)
{
  int i,k;

  // Add in variable j  

  for(i=1; i<=T.nrows(); i++) for(k=1; k<=T.ncols(); k++) 
    T(i,k) += Rdda->value(i,j,k);
  S(j) = 1;
}

void sdda(const Matrix &X, const IVector &y,
	  const Matrix &means, const Matrix &vars,
	  const IVector &counts, const Vector &priors,
	  IVector &S, const IVector &Start, const IVector &Never,
	  const bool prob, const bool cache, const bool xval,
	  const int jmax,
	  Vector &estore, Vector &pstore)
{
  int n = X.nrows(), p = X.ncols(), G = counts.length();
  int i,j,k;
  double ptmp;
  int err, etmp, kflg, jin;

  Matrix T(n,G);


  jin = 0;
  estore = 0.0;
  pstore = 0.0;

  // Instantiate a x-validated log likelihood helper object.
  // If cache is true precomputes just about everything.

  Rdda = new Robject(X,y,means,vars,counts,cache, xval);

  // Initialise the log-likelihood matrix to the priors

  for(k=1; k<=G; k++) {
    ptmp = -2*log(priors(k));
    for(i=1; i<=n; i++) T(i,k) = ptmp;
  }

  // Initialise the variable flags
  // S(j) = 0 means variable is not in.
  // S(j) = 1 means variable is in.

  // Start is a vector of flags for variables that should be put in
  // at the start. Never is a vector of variables that should never
  // go in.

  S = 0;
  for(j=1; j<=p; j++) if(Start(j)) {
    AddIn(S,T,j);
    jin++;
  }
  for(j=1; j<=p; j++) if(Never(j)) S(j) = 1;
    
  kflg = 1;
  for(k=2; k<=G; k++) if(priors(k)>priors(kflg)) kflg = k;

  err = 0;
  for(i=1; i<=n; i++) if(y(i) != kflg) err++;

  ptmp = 0.0;
  while(1) {
    if(prob) BestAddProb(y,S,T,j,etmp, ptmp);
    else BestAdd(y,S,T,j,etmp);
    jin++;
    if(jin<=SDDA_abs(jmax)) {
      estore(jin) = etmp;
      pstore(jin) = ptmp;
    }
    if (jmax<=0) if (etmp >= err) break;
    AddIn(S,T,j);
    if( jmax>0 && jin==jmax) break;
    err = etmp;
  }
      
  for(j=1; j<=p; j++) if(Never(j)) S(j) = 0;
  delete Rdda;
}
	  

	    
