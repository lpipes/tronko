#include "opt.h"

#define NMAX 10000000
#define NR_END 1
#define FREE
#define FREE_ARG char*
#define ALF 1.0e-4
#define EPS 3.0e-8
#define FMAX(a,b) (maxarg1=(a),maxarg2=(b),(maxarg1) > (maxarg2) ?\
        (maxarg1) : (maxarg2))
#define SQR(a) ((sqrarg=(a)) == 0.0 ? 0.0 : sqrarg*sqrarg)

double ITMAX;
double TOLX;
double TOLX2;
double STPMX;
double *dg,*g,*hdg,*pnew,*xi,**hessin; 
int npar, CENTRALMODE, PRECISIONLEVEL;
FILE *tempfile;
static double maxarg1,maxarg2;
static double sqrarg;
double oldf0;

void nrerror(char error_text[])
/* Numerical Recipes standard error handler */
{
	fprintf(stderr,"Numerical Recipes run-time error...\n");
	fprintf(stderr,"%s\n",error_text);
	fprintf(stderr,"...now exiting to system...\n");
	exit(1);
}

double *dvector(long nl, long nh)
/* allocate a double vector with subscript range v[nl..nh] */
{
	double *v;

	v=(double *)malloc((size_t) ((nh-nl+1+NR_END)*sizeof(double)));
	if (!v) nrerror("allocation failure in dvector()");
	return v-nl+NR_END;
}

void free_dvector(double *v, long nl)
/* free a double vector allocated with dvector() */
{
    free((FREE_ARG) (v+nl-NR_END));
    
}


double **dmatrix(long nrl, long nrh, long ncl, long nch)
/* allocate a double matrix with subscript range m[nrl..nrh][ncl..nch] */
{
	long i, nrow=nrh-nrl+1,ncol=nch-ncl+1;
	double **m;

	/* allocate pointers to rows */
	m=(double **) malloc((size_t)((nrow+NR_END)*sizeof(double*)));
	if (!m) nrerror("allocation failure 1 in matrix()");
	m += NR_END;
	m -= nrl;

	/* allocate rows and set pointers to them */
	m[nrl]=(double *) malloc((size_t)((nrow*ncol+NR_END)*sizeof(double)));
	if (!m[nrl]) nrerror("allocation failure 2 in matrix()");
	m[nrl] += NR_END;
	m[nrl] -= ncl;

	for(i=nrl+1;i<=nrh;i++) m[i]=m[i-1]+ncol;

	/* return pointer to array of pointers to rows */
	return m;
}


void free_dmatrix(double **m, long nrl, long ncl)
/* free a double matrix allocated by dmatrix() */
{
    free((FREE_ARG) (m[nrl]+ncl-NR_END));
    free((FREE_ARG) (m+nrl-NR_END));
}

/*int lnsrch(int n, double xold[], double fold, double g[], double p[], double x[],
	double *f, double stpmax, int *check, double (*func)(double []), double lowbound[], double upbound[])
{
	int i;
	double a,alam,alam2,alamin,b,disc,f2,fold2,rhs1,rhs2,slope,sum,temp,
		test,tmplam;
	double new1;

	//printf ("B: "); for (i=1;i<=n;i++) printf("%f ",p[i]);
	//printf("\n");
	*check=0;
	for (sum=0.0,i=1;i<=n;i++) sum += p[i]*p[i];
	sum=sqrt(sum);
	if (sum > stpmax)
		for (i=1;i<=n;i++) p[i] *= stpmax/sum;
	for (slope=0.0,i=1;i<=n;i++)
		slope += g[i]*p[i];
	test=0.0;
	for (i=1;i<=n;i++) {
		temp=fabs(p[i])/FMAX(fabs(xold[i]),1.0);
		if (temp > test) test=temp;
	}
        if (test==0.0)
                return 1;
	alamin=TOLX2/test;
	alam=1.0;
	for (;;) {
		for (i=1;i<=n;i++)
			{
			new1=xold[i]+alam*p[i];
			if (new1 < lowbound[i])	//my bloody code
				new1 = lowbound[i];
			if  (new1 > upbound[i])
				new1 = upbound[i];
			x[i] = new1;
			//printf("x[i]: %f,newin[i]: %f,xold[i]: %f,alam: %f,p[i]: %f\n",x[i],newin[i],xold[i],alam,p[i]);
			}
		*f=(*func)(x);
		if (alam < alamin) {
			for (i=1;i<=n;i++) x[i]=xold[i];
			*check=1;
			return 10;
		} else if (*f <= fold+ALF*alam*slope) return 10;
		else {
			if (alam == 1.0)
				tmplam = -slope/(2.0*(*f-fold-slope));
			else {
				rhs1 = *f-fold-alam*slope;
				rhs2=f2-fold2-alam2*slope;
				a=(rhs1/(alam*alam)-rhs2/(alam2*alam2))/(alam-alam2);
				b=(-alam2*rhs1/(alam*alam)+alam*rhs2/(alam2*alam2))/(alam-alam2);
				if (a == 0.0) tmplam = -slope/(2.0*b);
				else {
					disc=b*b-3.0*a*slope;
					if (disc<0.0) {*//*nrerror("Roundoff problem in lnsrch.") return -1;*//*disc=0.0;}
					*//*else*//*tmplam=(-b+sqrt(disc))/(3.0*a);
				}
				if (tmplam>0.5*alam)
					tmplam=0.5*alam;
			}
		}
		alam2=alam;
		f2 = *f;
		fold2=fold;
		alam=FMAX(tmplam,0.1*alam);
	}
}*/
int lnsrch_Arr(int n, double xold[], double fold, double g[], double p[], double x[],
	double *f, double stpmax, int *check, double (*func)(double [] , int), double lowbound[], double upbound[], int whichRoot)
{
	int i;
	double a,alam,alam2,alamin,b,disc,f2,fold2,rhs1,rhs2,slope,sum,temp,
		test,tmplam;
	double new1;

	/*printf ("B: "); for (i=1;i<=n;i++) printf("%f ",p[i]);
	printf("\n");*/
	*check=0;
	for (sum=0.0,i=1;i<=n;i++) sum += p[i]*p[i];
	sum=sqrt(sum);
	if (sum > stpmax)
		for (i=1;i<=n;i++) p[i] *= stpmax/sum;
	for (slope=0.0,i=1;i<=n;i++)
		slope += g[i]*p[i];
	test=0.0;
	for (i=1;i<=n;i++) {
		temp=fabs(p[i])/FMAX(fabs(xold[i]),1.0);
		if (temp > test) test=temp;
	}
        if (test==0.0)
                return 1;
	alamin=TOLX2/test;
	alam=1.0;
	for (;;) {
		for (i=1;i<=n;i++)
			{
			new1=xold[i]+alam*p[i];
			if (new1 < lowbound[i])	/*my bloody code*/
				new1 = lowbound[i];
			if  (new1 > upbound[i])
				new1 = upbound[i];
			x[i] = new1;
			/*printf("x[i]: %f,newin[i]: %f,xold[i]: %f,alam: %f,p[i]: %f\n",x[i],newin[i],xold[i],alam,p[i]);*/
			}
		*f=(*func)(x,whichRoot);
		if (alam < alamin) {
			for (i=1;i<=n;i++) x[i]=xold[i];
			*check=1;
			return 10;
		} else if (*f <= fold+ALF*alam*slope) return 10;
		else {
			if (alam == 1.0)
				tmplam = -slope/(2.0*(*f-fold-slope));
			else {
				rhs1 = *f-fold-alam*slope;
				rhs2=f2-fold2-alam2*slope;
				a=(rhs1/(alam*alam)-rhs2/(alam2*alam2))/(alam-alam2);
				b=(-alam2*rhs1/(alam*alam)+alam*rhs2/(alam2*alam2))/(alam-alam2);
				if (a == 0.0) tmplam = -slope/(2.0*b);
				else {
					disc=b*b-3.0*a*slope;
					if (disc<0.0) {/*nrerror("Roundoff problem in lnsrch.") return -1;*/disc=0.0;}
					/*else*/tmplam=(-b+sqrt(disc))/(3.0*a);
				}
				if (tmplam>0.5*alam)
					tmplam=0.5*alam;
			}
		}
		alam2=alam;
		f2 = *f;
		fold2=fold;
		alam=FMAX(tmplam,0.1*alam);
	}
}

void doNRinits(int n)

	{
	dg=dvector(1,n);
	g=dvector(1,n);
	hdg=dvector(1,n);
	hessin=dmatrix(1,n,1,n);
	pnew=dvector(1,n);
	xi=dvector(1,n);
	}

void freeNRinits(int n)

{
    free_dvector(dg, 1);
    free_dvector(g,1);
    free_dvector(hdg,1);
    free_dmatrix(hessin,1,1);
    free_dvector(pnew,1);
    free_dvector(xi,1);
}

	
/*void dfpmin(double p[], int n, double gtol, int *iter, double *fret,
	double(*func)(double []), void (*dfunc)(double [], double [],double [], double [], double(*fu)(double [])), double lowbound[], double upbound[])
{
	int lnsrch(int n, double xold[], double fold, double g[], double p[], double x[],
		 double *f, double stpmax, int *check, double (*func)(double []), double lowbound[], double upbound[]);
	int check,i,its,j;
	double den,fac,fad,fae,fp,stpmax,sum=0.0,sumdg,sumxi,temp,test;

	fp=(*func)(p);
	(*dfunc)(p,g, lowbound, upbound, func);
	for (i=1;i<=n;i++) {
		for (j=1;j<=n;j++) hessin[i][j]=0.0;
		hessin[i][i]=1.0;
		xi[i] = -g[i];
		sum += p[i]*p[i];
	}
	stpmax=STPMX*FMAX(sqrt(sum),(double)n);
	for (its=1;its<=ITMAX;its++) {
		*iter=its;
		printf ("A: "); for (i=1;i<=n;i++) printf("%f ",g[i]); printf("\n");
		if (lnsrch(n,p,fp,g,xi,pnew,fret,stpmax,&check,func,lowbound,upbound) == -1) //MY CODE
			return;		fp = *fret;
		for (i=1;i<=n;i++) {
			xi[i]=pnew[i]-p[i];
			p[i]=pnew[i];
		}
		test=0.0;
		for (i=1;i<=n;i++) {
			temp=fabs(xi[i])/FMAX(fabs(p[i]),1.0);
			if (temp > test) test=temp;
		}
		if (test < TOLX) {
			//FREEALL
			return;
		}
		for (i=1;i<=n;i++) dg[i]=g[i];
		(*dfunc)(p,g,lowbound,upbound, func);
	//	printf ("C: "); for (i=1;i<=n;i++) printf("%f ",g[i]); printf("\n");
		test=0.0;
		den=FMAX(*fret,1.0);
		for (i=1;i<=n;i++) {
			temp=fabs(g[i])*FMAX(fabs(p[i]),1.0)/den;
			if (temp > test) test=temp;
		}
		if (test < gtol) {
			//FREEALL
			return;
		}
		for (i=1;i<=n;i++) dg[i]=g[i]-dg[i];
		for (i=1;i<=n;i++) {
			hdg[i]=0.0;
			for (j=1;j<=n;j++) hdg[i] += hessin[i][j]*dg[j];
		}
		fac=fae=sumdg=sumxi=0.0;
		for (i=1;i<=n;i++) {
			fac += dg[i]*xi[i];
			fae += dg[i]*hdg[i];
			sumdg += SQR(dg[i]);
			sumxi += SQR(xi[i]);
		}
		if (fac*fac > EPS*sumdg*sumxi) {
			fac=1.0/fac;
			fad=1.0/fae;
			for (i=1;i<=n;i++) dg[i]=fac*xi[i]-fad*hdg[i];
			for (i=1;i<=n;i++) {
				for (j=1;j<=n;j++) {
					hessin[i][j] += (double)(fac*xi[i]*xi[j]
					-fad*hdg[i]*hdg[j]+fae*dg[i]*dg[j]);
				}
			}
		}
		for (i=1;i<=n;i++) {
			xi[i]=0.0;
			for (j=1;j<=n;j++) xi[i] -= hessin[i][j]*g[j];
		}
	}
	printf("too many iterations in dfpmin");
	//FREEALL
}*/
void dfpmin_Arr(double p[], int n, double gtol, int *iter, double *fret,
	double(*func)(double [],int), void (*dfunc)(double [], double [],double [], double [], double(*fu)(double [], int), int), double lowbound[], double upbound[], int whichRoot)
{
	int lnsrch_Arr(int n, double xold[], double fold, double g[], double p[], double x[],
		 double *f, double stpmax, int *check, double (*func)(double [],int), double lowbound[], double upbound[],int whichRoot);
	int check,i,its,j;
	double den,fac,fad,fae,fp,stpmax,sum=0.0,sumdg,sumxi,temp,test;

	fp=(*func)(p,whichRoot);
	(*dfunc)(p,g, lowbound, upbound, func,whichRoot);
	for (i=1;i<=n;i++) {
		for (j=1;j<=n;j++) hessin[i][j]=0.0;
		hessin[i][i]=1.0;
		xi[i] = -g[i];
		sum += p[i]*p[i];
	}
	stpmax=STPMX*FMAX(sqrt(sum),(double)n);
	for (its=1;its<=ITMAX;its++) {
		*iter=its;
		/*printf ("A: "); for (i=1;i<=n;i++) printf("%f ",g[i]); printf("\n");*/
		if (lnsrch_Arr(n,p,fp,g,xi,pnew,fret,stpmax,&check,func,lowbound,upbound,whichRoot) == -1) /*MY CODE*/
			return;		fp = *fret;
		for (i=1;i<=n;i++) {
			xi[i]=pnew[i]-p[i];
			p[i]=pnew[i];
		}
		test=0.0;
		for (i=1;i<=n;i++) {
			temp=fabs(xi[i])/FMAX(fabs(p[i]),1.0);
			if (temp > test) test=temp;
		}
		if (test < TOLX) {
			/*FREEALL*/
			return;
		}
		for (i=1;i<=n;i++) dg[i]=g[i];
		(*dfunc)(p,g,lowbound,upbound, func, whichRoot);
	/*	printf ("C: "); for (i=1;i<=n;i++) printf("%f ",g[i]); printf("\n");*/
		test=0.0;
		den=FMAX(*fret,1.0);
		for (i=1;i<=n;i++) {
			temp=fabs(g[i])*FMAX(fabs(p[i]),1.0)/den;
			if (temp > test) test=temp;
		}
		if (test < gtol) {
			/*FREEALL*/
			return;
		}
		for (i=1;i<=n;i++) dg[i]=g[i]-dg[i];
		for (i=1;i<=n;i++) {
			hdg[i]=0.0;
			for (j=1;j<=n;j++) hdg[i] += hessin[i][j]*dg[j];
		}
		fac=fae=sumdg=sumxi=0.0;
		for (i=1;i<=n;i++) {
			fac += dg[i]*xi[i];
			fae += dg[i]*hdg[i];
			sumdg += SQR(dg[i]);
			sumxi += SQR(xi[i]);
		}
		if (fac*fac > EPS*sumdg*sumxi) {
			fac=1.0/fac;
			fad=1.0/fae;
			for (i=1;i<=n;i++) dg[i]=fac*xi[i]-fad*hdg[i];
			for (i=1;i<=n;i++) {
				for (j=1;j<=n;j++) {
					hessin[i][j] += (double)(fac*xi[i]*xi[j]
					-fad*hdg[i]*hdg[j]+fae*dg[i]*dg[j]);
				}
			}
		}
		for (i=1;i<=n;i++) {
			xi[i]=0.0;
			for (j=1;j<=n;j++) xi[i] -= hessin[i][j]*g[j];
		}
	}
	printf("too many iterations in dfpmin");
	/*FREEALL*/
}

/*void Yanggradient (int n, double x[], double f0, double g[],
    double (*fun)(double x[]), double space[], int central, double lowbound[], double upbound[])
{

*/	/*f0=fun(x) is given for Central=0*/

  /* int i,j;
   double *x0=space, *x1=space+n, eh0=1e-8, eh01=1e-8, eh;

   if (central) {
      for (i=1;i<=n;i++)  {
	 for (j=1;j<=n;j++)  x0[j]=x1[j]=x[j];
	 eh=pow(eh01*(fabs(x[i])+1), 0.67);
	 x0[i]-=eh; x1[i]+=eh;
	if (x0[i]<lowbound[i])
		{x1[i]+=eh; g[i] = ((*fun)(x1)-f0)/(eh*2.0);}
	else if (x1[i]>upbound[i])
		{x0[i]-=eh; g[i] = (f0-(*fun)(x0))/(eh*2.0);}
	else
	 	g[i] = ((*fun)(x1) - (*fun)(x0))/(eh*2.0);
	if (x[i] <= lowbound[i] && g[i] > 0.0)
		g[i] = 0.0;
	else if (x[i] >= upbound[i] && g[i] < 0.0)
		g[i] = 0.0; 
     }
   }
   else {
      for (i=1;i<=n;i++)  {
	 for (j=1;j<=n;j++)  x1[j]=x[j];
	*/ /*eh=eh0*(fabs(x[i])+1);*/
	/* eh=2.0*pow(eh0*(fabs(x[i])+1), 0.67);
	if (x1[i]+eh>upbound[i])
		{
	 	x1[i]-=eh;
	 	g[i] = (f0-(*fun)(x1))/eh;
	 	}
	else
		{
	 	x1[i]+=eh;
	 	g[i] = ((*fun)(x1)-f0)/eh;
	 	}
	if (x[i] <= lowbound[i] && g[i] > 0.0)
		g[i] = 0.0;
	else if (x[i] >= upbound[i] && g[i] < 0.0)
		g[i] = 0.0;
      }
   }
}*/
void Yanggradient_Arr (int n, double x[], double f0, double g[],
    double (*fun)(double x[] ,int), double space[], int central, double lowbound[], double upbound[],int whichRoot)
{

	/*f0=fun(x) is given for Central=0*/

   int i,j;
   double *x0=space, *x1=space+n, eh0=1e-8, eh01=1e-8, eh;

   if (central) {
      for (i=1;i<=n;i++)  {
	 for (j=1;j<=n;j++)  x0[j]=x1[j]=x[j];
	 eh=pow(eh01*(fabs(x[i])+1), 0.67);
	 x0[i]-=eh; x1[i]+=eh;
	if (x0[i]<lowbound[i])
		{x1[i]+=eh; g[i] = ((*fun)(x1,whichRoot)-f0)/(eh*2.0);}
	else if (x1[i]>upbound[i])
		{x0[i]-=eh; g[i] = (f0-(*fun)(x0,whichRoot))/(eh*2.0);}
	else
	 	g[i] = ((*fun)(x1,whichRoot) - (*fun)(x0,whichRoot))/(eh*2.0);
	if (x[i] <= lowbound[i] && g[i] > 0.0)
		g[i] = 0.0;
	else if (x[i] >= upbound[i] && g[i] < 0.0)
		g[i] = 0.0; 
     }
   }
   else {
      for (i=1;i<=n;i++)  {
	 for (j=1;j<=n;j++)  x1[j]=x[j];
	 /*eh=eh0*(fabs(x[i])+1);*/
	 eh=2.0*pow(eh0*(fabs(x[i])+1), 0.67);
	if (x1[i]+eh>upbound[i])
		{
	 	x1[i]-=eh;
	 	g[i] = (f0-(*fun)(x1,whichRoot))/eh;
	 	}
	else
		{
	 	x1[i]+=eh;
	 	g[i] = ((*fun)(x1,whichRoot)-f0)/eh;
	 	}
	if (x[i] <= lowbound[i] && g[i] > 0.0)
		g[i] = 0.0;
	else if (x[i] >= upbound[i] && g[i] < 0.0)
		g[i] = 0.0;
      }
   }
}

/*void getgradient(double invec[], double outvec[], double lowbound[], double upbound[], double(*func)(double []))

	{
	int i;
	static double space[200];
	double f0, arbitrarysum = 0.0;

	f0 = func(invec);
    if (CENTRALMODE > 0 && PRECISIONLEVEL >1)
		Yanggradient(npar, invec, f0, outvec, func, space, 1, lowbound, upbound);
	else 
		Yanggradient(npar, invec, f0, outvec, func, space, 0, lowbound, upbound);
//	printf("GRADIENT:\n");
	 *//*if (NULL==(tempfile=fopen("tempfile","w"))){
  		puts ("Kan ikke aabne tempfile!");
   		exit(-1);}*/
	/*for(i=1; i<=npar; i++)
		{
	*//*	fprintf(tempfile,"%f ",invec[i]);  */
//		printf(" %5f",outvec[i]);
	/*	if (invec[i] > lowbound[i] && invec[i] < upbound[i])
			arbitrarysum = arbitrarysum + fabs(outvec[i]);
	*/	/*else printf("*");*/
	//	}
  /* 	fclose(tempfile);*/
	/*if (CENTRALMODE == 0 && fabs(oldf0-f0) < 0.1)
		{
		CENTRALMODE = 1;
	*/	/*printf("Changing to central method\n");*/
	//	}
	/*else if (CENTRALMODE > 0 && fabs(oldf0-f0) > 1.0)
		CENTRALMODE = 0;
	else if (arbitrarysum < 0.01*npar)
		CENTRALMODE = 2;
	oldf0 = f0;
//	printf(" Like (grad): %f\n",-oldf0);
	}*/
void getgradient_Arr(double invec[], double outvec[], double lowbound[], double upbound[], double(*func)(double [], int), int whichRoot)

	{
	int i;
	static double space[200];
	double f0, arbitrarysum = 0.0;

	f0 = func(invec,whichRoot);
    if (CENTRALMODE > 0 && PRECISIONLEVEL >1)
		Yanggradient_Arr(npar, invec, f0, outvec, func, space, 1, lowbound, upbound, whichRoot);
	else 
		Yanggradient_Arr(npar, invec, f0, outvec, func, space, 0, lowbound, upbound, whichRoot);
//	printf("GRADIENT:\n");
	 /*if (NULL==(tempfile=fopen("tempfile","w"))){
  		puts ("Kan ikke aabne tempfile!");
   		exit(-1);}*/
	for(i=1; i<=npar; i++)
		{
	/*	fprintf(tempfile,"%f ",invec[i]);  */
//		printf(" %5f",outvec[i]);
		if (invec[i] > lowbound[i] && invec[i] < upbound[i])
			arbitrarysum = arbitrarysum + fabs(outvec[i]);
		/*else printf("*");*/
		}
  /* 	fclose(tempfile);*/
	if (CENTRALMODE == 0 && fabs(oldf0-f0) < 0.1)
		{
		CENTRALMODE = 1;
		/*printf("Changing to central method\n");*/
		}
	else if (CENTRALMODE > 0 && fabs(oldf0-f0) > 1.0)
		CENTRALMODE = 0;
	else if (arbitrarysum < 0.01*npar)
		CENTRALMODE = 2;
	oldf0 = f0;
//	printf(" Like (grad): %f\n",-oldf0);
	}
/*Call nrinits before calling this the first time.  Precision level can be 0, 1 or 2*/
/*double findmax(double newinvecter[], double lowbound[], double upbound[], int n, double (*fun)(double x[]), int precisionlevel)

	{

	double gtol, fret;
	int i, iter;
	
	oldf0 = 0.0;
	npar = n;
    PRECISIONLEVEL = precisionlevel;
        
        if (PRECISIONLEVEL == 2)
            {
            gtol = 0.0000000001;
            ITMAX = 500;
            TOLX = 4*EPS;
            TOLX2  = 1.0e-7;
            STPMX = 200.0;
            }
        else if (PRECISIONLEVEL == 1){
            gtol = 0.0001;
            ITMAX = 100;
            TOLX = 10*EPS;
            TOLX2  = 1.0e-5;
            STPMX = 50.0;
            }
        else {
            gtol = 0.01;
            ITMAX = 50;
            TOLX = 100*EPS;
            TOLX2  = 1.0e-3;
            STPMX = 10.0;
            }
		for (i=1; i<=n; i++){
            if (newinvecter[i]==lowbound[i])
                newinvecter[i]= newinvecter[i]+fabs(1.05*newinvecter[i]);
            else if (newinvecter[i]==upbound[i])
                newinvecter[i]= newinvecter[i]-fabs(1.05*newinvecter[i]);
        }
        
	//do    {
	dfpmin(newinvecter, npar, gtol, &iter, &fret, fun, getgradient, lowbound, upbound);
	//}while (CENTRALMODE < 2 && PRECISIONLEVEL == 2);
	return -fret;
	}*/
double findmax_Arr(double newinvecter[], double lowbound[], double upbound[], int n, double (*fun)(double x[], int), int precisionlevel, int whichRoot)

	{

	double gtol, fret;
	int i, iter;
	
	oldf0 = 0.0;
	npar = n;
    PRECISIONLEVEL = precisionlevel;
        
        if (PRECISIONLEVEL == 2)
            {
            gtol = 0.0000000001;
            ITMAX = 500;
            TOLX = 4*EPS;
            TOLX2  = 1.0e-7;
            STPMX = 200.0;
            }
        else if (PRECISIONLEVEL == 1){
            gtol = 0.0001;
            ITMAX = 100;
            TOLX = 10*EPS;
            TOLX2  = 1.0e-5;
            STPMX = 50.0;
            }
        else {
            gtol = 0.01;
            ITMAX = 50;
            TOLX = 100*EPS;
            TOLX2  = 1.0e-3;
            STPMX = 10.0;
            }
		for (i=1; i<=n; i++){
            if (newinvecter[i]==lowbound[i])
                newinvecter[i]= newinvecter[i]+fabs(1.05*newinvecter[i]);
            else if (newinvecter[i]==upbound[i])
                newinvecter[i]= newinvecter[i]-fabs(1.05*newinvecter[i]);
        }
        
	//do    {
	dfpmin_Arr(newinvecter, npar, gtol, &iter, &fret, fun, getgradient_Arr, lowbound, upbound, whichRoot);
	//}while (CENTRALMODE < 2 && PRECISIONLEVEL == 2);
	return (double) -fret;
	}
