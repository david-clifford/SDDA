#include "sdda.h"
#include "extra.h"

extern "C" {

  void dldawrap(int *nref, int *pref, int *Gref,
		double *Xref, int *yref,
		double *mref, double *vref, int *cref)
  {
    int n = *nref, p = *pref, G = *Gref;
    
    Matrix X(Xref,n,p);
    IVector y(yref,n);
    Matrix means(mref,G,p);
    Vector vars(vref,p);
    IVector counts(cref,G);

    dlda(X,y,means,vars,counts);
  }


  void dqdawrap(int *nref, int *pref, int *Gref,
		double *Xref, int *yref,
		double *mref, double *vref, int *cref)
  {
    int n = *nref, p = *pref, G = *Gref;
    
    Matrix X(Xref,n,p);
    IVector y(yref,n);
    Matrix means(mref,G,p);
    Matrix vars(vref,G,p);
    IVector counts(cref,G); 
    dqda(X,y,means,vars,counts);
  }

  void preddldawrap(int *pref, int *Gref,
		    double *mref, double *vref, int *cref, double *qref,
		    int *nref, double *Xref, double *Rref)
  {
    int n = *nref, p = *pref, G = *Gref;
    
    Matrix means(mref,G,p);
    Vector vars(vref,p);
    IVector counts(cref,G);
    Vector priors(qref,G);

    Matrix X(Xref,n,p);
    Matrix R(Rref,n,G);

    pred_dlda(means,vars,counts,priors,X,R);
  }

  void preddqdawrap(int *pref, int *Gref,
		    double *mref, double *vref, int *cref, double *qref,
		    int *nref, double *Xref, double *Rref)
  {
    int n = *nref, p = *pref, G = *Gref;
    
    Matrix means(mref,G,p);
    Matrix vars(vref,G,p);
    IVector counts(cref,G);
    Vector priors(qref,G);

    Matrix X(Xref,n,p);
    Matrix R(Rref,n,G);

    pred_dqda(means,vars,counts,priors,X,R);
  }


  void sddawrap(int *nref, int *pref, int *Gref,
		double *Xref, int *yref,
		double *mref, double *vref, int *cref, double *qref,
		int *startref, int *neverref, int *sref, 
		int *lda, int *iprob, int *icache, int *ixval, int *jmref,
		double *estref, double *pstref )
  {
    int n = *nref, p = *pref, G = *Gref;
    int vdim;

    // 
    // Main Stepwise Diagonal Discriminant Analysis Function
    //

    if (*lda) vdim = 1;
    else vdim = G;
    
    Matrix X(Xref,n,p);
    IVector y(yref,n);
    Matrix means(mref,G,p);
    Matrix vars(vref,vdim,p);
    IVector counts(cref,G);
    Vector priors(qref,G);

    IVector S(sref, p);
    IVector Start(startref, p);
    IVector Never(neverref, p);


    int jmax;
    bool prob, cache, xval;

    jmax = *jmref;

    Vector estore(estref, SDDA_max(1,SDDA_abs(jmax)));
    Vector pstore(pstref, SDDA_max(1,SDDA_abs(jmax)));


    if(*iprob==0) prob=false;
    else prob=true;

    if(*icache==0) cache = false;
    else cache=true;

    if(*ixval==0) xval = false;
    else xval = true;

    sdda(X,y,means,vars,counts,priors,
	 S,Start,Never,
	 prob, cache, xval, jmax,
	 estore, pstore);

  }

}
