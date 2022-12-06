#include <assert.h>
#include "likelihood.h"


/*void makecon(int node)

{
  int i, j;
  double L1, L2, max, m;

  if (tree[node].up[0]==-1 && tree[node].mrca < 0) {printf("Error in identication of MRCAs\n"); exit(-1);}*//*this can eventually be removed*/
  /*else if (tree[node].up[0]==-1){
    for (i=0; i<STATESPACE; i++)
      tree[node].like[i] = 1.0;
  }
  else if (tree[node].mrca > -1){
    for (i=0; i<STATESPACE; i++)
      tree[node].like[i] = pdata(statevector[i], tree[node].s, tree[node].numsites, numbase);*//*here we would also need to change 'numbases to match the number of bases in the alignemnt*/
 /* }
  else {
    makecon(tree[node].up[0]);
    makecon(tree[node].up[1]);
    maketransitionmatrix(0, tree[tree[node].up[0]].bl);
    maketransitionmatrix(1, tree[tree[node].up[1]].bl);
    //printf("node %i\n",node);
    max=0.0;
    for (i=0; i<STATESPACE; i++){
      L1=L2=0;
      for (j=0; j<STATESPACE; j++){
        L1 += PMAT1[i][j]*tree[tree[node].up[0]].like[j];
        L2 += PMAT2[i][j]*tree[tree[node].up[1]].like[j];
      }
      //printf("(%lf, %lf) ",L1, L2);
      if ((tree[node].like[i]=L1*L2)>max)
        max = tree[node].like[i];
    } //printf("\n");
    //underflow control
    for (i=0; i<STATESPACE; i++)
      tree[node].like[i]=tree[node].like[i]/max;
    m = log(max);
    UFC = UFC + m;
  }

}*/
void makecon_Arr(int node, int whichRoot){
  int i, j;
  double L1, L2, max, m;

  if (treeArr[whichRoot][node].up[0]==-1 && treeArr[whichRoot][node].mrca < 0) {printf("Error in identication of MRCAs\n"); exit(-1);}/*this can eventually be removed*/
  else if (treeArr[whichRoot][node].up[0]==-1 && node != 200){
    for (i=0; i<STATESPACE; i++)
      treeArr[whichRoot][node].like[i] = 1.0;
  }
  else if (treeArr[whichRoot][node].mrca > -1 || node == 200){
    for (i=0; i<STATESPACE; i++)
      treeArr[whichRoot][node].like[i] = pdata(statevector[i], treeArr[whichRoot][node].s, treeArr[whichRoot][node].numsites, numbaseArr[whichRoot]);/*here we would also need to change 'numbases to match the number of bases in the alignemnt*/
  }
  else {
    makecon_Arr(treeArr[whichRoot][node].up[0],whichRoot);
    makecon_Arr(treeArr[whichRoot][node].up[1],whichRoot);
    maketransitionmatrix(0, treeArr[whichRoot][treeArr[whichRoot][node].up[0]].bl);
    maketransitionmatrix(1, treeArr[whichRoot][treeArr[whichRoot][node].up[1]].bl);
    printf("node %i\n",node);
    max=0.0;
    for (i=0; i<STATESPACE; i++){
      L1=L2=0;
      for (j=0; j<STATESPACE; j++){
        L1 += PMAT1[i][j]*treeArr[whichRoot][treeArr[whichRoot][node].up[0]].like[j];
        L2 += PMAT2[i][j]*treeArr[whichRoot][treeArr[whichRoot][node].up[1]].like[j];
      }
      //printf("(%lf, %lf) ",L1, L2);
      if ((treeArr[whichRoot][node].like[i]=L1*L2)>max)
        max = treeArr[whichRoot][node].like[i];
    } //printf("\n");
    //underflow control
    for (i=0; i<STATESPACE; i++)
      treeArr[whichRoot][node].like[i]=treeArr[whichRoot][node].like[i]/max;
    m = log(max);
    UFC = UFC + m;
  }

}

/*double tlike()

{ *//*assumes uniform prior distribution on the state space*/
/*
  int i, j, s;
  double L1, L2, Lsum=0;

  if (tree[root].mrca>-1) {printf("Data contain only one species\n"); exit(-1);}
  UFC=0.0;
  makecon(tree[root].up[0]);
  makecon(tree[root].up[1]);
  maketransitionmatrix(0, tree[tree[root].up[0]].bl);
  maketransitionmatrix(1, tree[tree[root].up[1]].bl);
  for (i=0; i<STATESPACE; i++){
    L1=L2=0;
    for (j=0; j<STATESPACE; j++){
      L1 += PMAT1[i][j]*tree[tree[root].up[0]].like[j];
      L2 += PMAT2[i][j]*tree[tree[root].up[1]].like[j];
    }
    Lsum += L1*L2;
  }
  return log(Lsum)+UFC; *//*the stationary distribution is uniform after scaling so we do not need to multiply with the stationary probabilities*/
/*}*/
double tlike_Arr(int whichRoot)

{ /*assumes uniform prior distribution on the state space*/

  int i, j, s;
  double L1, L2, Lsum=0;

  if (treeArr[whichRoot][rootArr[whichRoot]].mrca>-1) {printf("Data contain only one species\n"); exit(-1);}
  UFC=0.0;
  makecon_Arr(treeArr[whichRoot][rootArr[whichRoot]].up[0],whichRoot);
  makecon_Arr(treeArr[whichRoot][rootArr[whichRoot]].up[1],whichRoot);
  maketransitionmatrix(0, treeArr[whichRoot][treeArr[whichRoot][rootArr[whichRoot]].up[0]].bl);
  maketransitionmatrix(1, treeArr[whichRoot][treeArr[whichRoot][rootArr[whichRoot]].up[1]].bl);
  for (i=0; i<STATESPACE; i++){
    L1=L2=0;
    for (j=0; j<STATESPACE; j++){
      L1 += PMAT1[i][j]*treeArr[whichRoot][treeArr[whichRoot][rootArr[whichRoot]].up[0]].like[j];
      L2 += PMAT2[i][j]*treeArr[whichRoot][treeArr[whichRoot][rootArr[whichRoot]].up[1]].like[j];
    }
    Lsum += L1*L2;
  }
  return log(Lsum)+UFC; /*the stationary distribution is uniform after scaling so we do not need to multiply with the stationary probabilities*/
}
//optimization function starts counting at 1 so arrays have dimensionality n+1
/*double getlikelihood(double par[4])
{

  double like, locpar[2];

  locpar[0]  = par[1];
  locpar[1] = par[2];

  definegammaquantiles(STATESPACE, locpar);
  inittransitionmatrix(par[3]);
  like = tlike();
  return -like+log((double)STATESPACE);
}
*/
double getlikelihood_Arr(double par[4], int whichRoot){

  double like, locpar[2];

  locpar[0]  = par[1];
  locpar[1] = par[2];

  definegammaquantiles(STATESPACE, locpar);
  inittransitionmatrix(par[3]);
  like = tlike_Arr(whichRoot);
  return -like+log((double)STATESPACE);
}

/*double maximizelikelihood(double parameters[4], int precision)
{

  //optimization function starts counting at 1 so arrays have dimensionality n+1
  double L, lowbound[4],  upbound[4];

  doNRinits(4);
  parameters[1]=parameters[2]=parameters[3]=1.0;
  lowbound[1]=lowbound[2]=lowbound[3]=0.000001;
  upbound[1]=upbound[2]=100.0; upbound[3]=1000.0;
  L=findmax(parameters, lowbound, upbound, 3, getlikelihood, precision);
  freeNRinits(4);
  return L;
}*/
double maximizelikelihood_Arr(double parameters[4], int precision, int whichRoot)
{

  //optimization function starts counting at 1 so arrays have dimensionality n+1
  double L, lowbound[4],  upbound[4];

  doNRinits(4);
  parameters[1]=parameters[2]=parameters[3]=1.0;
  lowbound[1]=lowbound[2]=lowbound[3]=0.000001;
  upbound[1]=upbound[2]=100.0; upbound[3]=1000.0;
  //L=findmax(parameters, lowbound, upbound, 3, getlikelihood, precision);
  L=findmax_Arr(parameters, lowbound, upbound, 3, getlikelihood_Arr, precision, whichRoot);
  freeNRinits(4);
  return L;
}

//this algorithm runs from the root to calculate the conditional likelihood given all data towards the parental node for each node
/*void makeposterior(int node)
{
  int i,j, parent, otherb;
  double bl, sum, templike[STATESPACE];

  parent = tree[node].down;
  bl = tree[node].bl;
  maketransitionmatrix(0, bl);
  if ((otherb = tree[parent].up[0])==node)
    otherb = tree[parent].up[1];
  maketransitionmatrix(1, tree[otherb].bl);
  for (i=0; i<STATESPACE; i++){
    templike[i]=0;
    for (j=0; j<STATESPACE; j++)
      templike[i]=templike[i]+tree[otherb].like[j]*PMAT2[i][j];
    templike[i]=templike[i]*tree[parent].posterior[i];
  }
  for (i=0; i<STATESPACE; i++){
    tree[node].posterior[i]=0.0;
    for (j=0; j<STATESPACE; j++)
      tree[node].posterior[i] = tree[node].posterior[i] + PMAT1[i][j]*templike[j];
  }
  if (tree[node].up[0]>-1 && tree[node].mrca == -1)
  {
    makeposterior(tree[node].up[0]);
    makeposterior(tree[node].up[1]);
  }
}
*/
//finds the posterior for all nodes that are mrca for a species
//for deeper nodes it gives the conditional likelihood given all data towards the parental node

/*void getposterior(double par[4])
{
  int i, j;
  double sum;


  par[1]=5.0; par[2]=1.0; par[3]=500.0;


  getlikelihood(par); *//*need to call likelihood again*/
 /* for (i=0; i<STATESPACE; i++)
    tree[root].posterior[i]=1.0;
  makeposterior(tree[root].up[0]);
  makeposterior(tree[root].up[1]);
  for (j=0; j<2*numspec-1; j++){
    if (tree[j].mrca > -1){
      sum = 0.0;
      for (i=0; i<STATESPACE; i++)
      {
        if (tree[j].up[0]>-1)
          tree[j].posterior[i] = tree[j].posterior[i]*pdata(statevector[i], tree[j].s, tree[j].numsites, numbase); *//*here we would also need to change 'numbases'*/
        /*sum = sum + tree[j].posterior[i];
      }
      for (i=0; i<STATESPACE; i++)
        tree[j].posterior[i] = tree[j].posterior[i]/sum;
    }
  }
}*/
/*
void maketransitionmatrixnc(int n, double t){
	int i, j, k;
	double EXPOS[4];
	for (k=0; k<4; k++){
		//printf("RRVALnc[0]=%lf,RRVALnc[1]=%lf,RRVALnc[2]=%lf,RRVALnc[3]=%lf\n",RRVALnc[0],RRVALnc[1],RRVALnc[2],RRVALnc[3]);
		EXPOS[k] = exp(t*RRVALnc[k]);
		//assert(EXPOS[k]>=0.0 && EXPOS[k]<=1.0);
		if ( EXPOS[k] < 0.000000 ){
			printf("EXPOS[k]=%e\n",EXPOS[k]);
		}
		//if (EXPOS[k] > 1.000000 ){
		//	printf("1-EXPOS[k]=%e\n",1-EXPOS[k]);
		//}
		assert(EXPOS[k]>=0.000000);
		//assert(EXPOS[k]<=1.000000);
	}

	for (i=0; i<4; i++)
  {
    for (j=0; j<4; j++)
    {
      PMATnc[n][i][j] = 0.0;
      for (k=0; k<4; k++){
		  //printf("PMATnc[%d][%d][%d]=%lf\tRRVECnc[%d][%d]=%lf\tLRVECnc[%d][%d]=%lf\tEXPOS[%d]=%lf\n",n,i,j,PMATnc[n][i][j],k,j,RRVECnc[k][j],i,k,LRVECnc[i][k],k,EXPOS[k]);
		  PMATnc[n][i][j] =  PMATnc[n][i][j] + RRVECnc[k][j]*LRVECnc[i][k]*EXPOS[k];
		  //printf("RRVECnc[%d][%d]*LRVECnc[%d][%d]*EXPOS[%d]=%lf\n",k,j,i,k,k,RRVECnc[k][j]*LRVECnc[i][k]*EXPOS[k]);
		  //printf("PMATnc[%d][%d][%d]=%lf\n",n,i,j,PMATnc[n][i][j]);
		  //assert(PMATnc[n][i][j]>=0);
		  //assert(PMATnc[n][i][j]<=1);
	  }
	}
  }

  /*  if (t<0.0001){
      printf("PMAT[%i]nc (t=%lf)\n",n,t);
      for (i=0; i<4; i++){
      for (j=0; j<4; j++)
      printf("%.10lf ",PMATnc[n][i][j]);
      printf("\n");printf("\n");
      }}*/

//}*/
void maketransitionmatrixnc_Arr(int n, double t, int whichRoot){
	int i, j, k;
	double EXPOS[4];
	for (k=0; k<4; k++){
		EXPOS[k] = exp(t*RRVALnc[k]);
		if ( EXPOS[k] < 0.000000 ){
			printf("EXPOS[k]=%e\n",EXPOS[k]);
		}
		//assert(EXPOS[k]>=0.000000);
	}
	for (i=0; i<4; i++){
    	for (j=0; j<4; j++){
      		PMATnc[n][i][j] = 0.0;
      		for (k=0; k<4; k++){
		  		PMATnc[n][i][j] =  PMATnc[n][i][j] + RRVECnc[k][j]*LRVECnc[i][k]*EXPOS[k];
	  		}
		}
  	}
}
/*void makeconnc(int node, double lambda)

{
  int i, j, child, seqn, site;
  double L, max;

  child = tree[node].up[0];

    //if (node ==root) printf("node: %i, child: %i, UP0: %i\n",node,child,tree[child].up[0]);

  if (tree[child].up[0]==-1){
    maketransitionmatrixnc(0, lambda*tree[child].bl);
    seqn=child-numspec+1;
    for (site=0; site<numbase; site++){
      for (i=0; i<4; i++)
        tree[node].likenc[site][i] = PMATnc[0][i][seq[seqn][site]];
       if (site==0 && node ==root) printf("node 10 (b=%i): %lf %lf %lf %lf\n",seq[seqn][site],PMATnc[0][0][seq[seqn][site]],PMATnc[0][1][seq[seqn][site]],PMATnc[0][2][seq[seqn][site]],PMATnc[0][3][seq[seqn][site]]);
    }
  }
  else {
    makeconnc(child, lambda);
    maketransitionmatrixnc(0, lambda*tree[child].bl);
    for (site=0; site<numbase; site++)
    {
      for (i=0; i<4; i++){
        tree[node].likenc[site][i]=0.0;
        for (j=0; j<4; j++)
          tree[node].likenc[site][i] += PMATnc[0][i][j]*tree[child].likenc[site][j];
      }
    }
  }
  child = tree[node].up[1];
  if (tree[child].up[1]==-1){
    maketransitionmatrixnc(0,lambda*tree[child].bl);
    seqn=child-numspec+1;
    for (site=0; site<numbase; site++)
      for (i=0; i<4; i++)
        tree[node].likenc[site][i] = tree[node].likenc[site][i]*PMATnc[0][i][seq[seqn][site]];
  }
  else {
    makeconnc(child, lambda);
    maketransitionmatrixnc(0,lambda*tree[child].bl);
    for (site=0; site<numbase; site++)
    {
      max=0.0;
      for (i=0; i<4; i++){
        L=0.0;
        for (j=0; j<4; j++)
          L += PMATnc[0][i][j]*tree[child].likenc[site][j];
        //underflow control
        if ((tree[node].likenc[site][i] = tree[node].likenc[site][i]*L)>max)
          max = tree[node].likenc[site][i];
      }
      //might be worth with some code here dealing with the case of max=0+epsilon

      if (max<0.00000000001) printf("Warning, max = %lf\n",max);

      for (i=0; i<4; i++)
        tree[node].likenc[site][i]=tree[node].likenc[site][i]/max;
      //  if (node ==root && site==0) printf("Node %i, Conditional likes: %lf %lf %lf %lf, ",node,tree[node].likenc[site][0],tree[node].likenc[site][1],tree[node].likenc[site][2],tree[node].likenc[site][3]);
      UFCnc[site] = UFCnc[site] + log(max);
      //  printf("max: %lf, UFC: %lf\n", max,UFCnc[site]);
    }
  }
}*/
void makeconnc_Arr(int node, double lambda, int whichRoot){
  int i, j, child, seqn, site;
  double L, max;
  child = treeArr[whichRoot][node].up[0];
  if (treeArr[whichRoot][child].up[0]==-1){
    maketransitionmatrixnc_Arr(0, lambda*treeArr[whichRoot][child].bl,whichRoot);
    seqn=child-numspecArr[whichRoot]+1;
    for (site=0; site<numbaseArr[whichRoot]; site++){
      for (i=0; i<4; i++)
        treeArr[whichRoot][node].likenc[site][i] = PMATnc[0][i][seqArr[whichRoot][seqn][site]];
       if (site==0 && node ==rootArr[whichRoot]) printf("node 10 (b=%i): %lf %lf %lf %lf\n",seqArr[whichRoot][seqn][site],PMATnc[0][0][seqArr[whichRoot][seqn][site]],PMATnc[0][1][seqArr[whichRoot][seqn][site]],PMATnc[0][2][seqArr[whichRoot][seqn][site]],PMATnc[0][3][seqArr[whichRoot][seqn][site]]);
    }
  }
  else {
    makeconnc_Arr(child, lambda, whichRoot);
    maketransitionmatrixnc_Arr(0, lambda*treeArr[whichRoot][child].bl,whichRoot);
    for (site=0; site<numbaseArr[whichRoot]; site++){
      for (i=0; i<4; i++){
        treeArr[whichRoot][node].likenc[site][i]=0.0;
        for (j=0; j<4; j++)
          treeArr[whichRoot][node].likenc[site][i] += PMATnc[0][i][j]*treeArr[whichRoot][child].likenc[site][j];
      }
    }
  }
  child = treeArr[whichRoot][node].up[1];
  if (treeArr[whichRoot][child].up[1]==-1){
    maketransitionmatrixnc_Arr(0,lambda*treeArr[whichRoot][child].bl,whichRoot);
    seqn=child-numspecArr[whichRoot]+1;
    for (site=0; site<numbaseArr[whichRoot]; site++)
      for (i=0; i<4; i++)
        treeArr[whichRoot][node].likenc[site][i] = treeArr[whichRoot][node].likenc[site][i]*PMATnc[0][i][seqArr[whichRoot][seqn][site]];
  }
  else {
    makeconnc_Arr(child, lambda,whichRoot);
    maketransitionmatrixnc_Arr(0,lambda*treeArr[whichRoot][child].bl,whichRoot);
    for (site=0; site<numbaseArr[whichRoot]; site++)
    {
      max=0.0;
      for (i=0; i<4; i++){
        L=0.0;
        for (j=0; j<4; j++)
          L += PMATnc[0][i][j]*treeArr[whichRoot][child].likenc[site][j];
        //underflow control
        if ((treeArr[whichRoot][node].likenc[site][i] = treeArr[whichRoot][node].likenc[site][i]*L)>max)
          max = treeArr[whichRoot][node].likenc[site][i];
      }
      //might be worth with some code here dealing with the case of max=0+epsilon

      if (max<0.00000000001) printf("Warning, max = %lf\n",max);

      for (i=0; i<4; i++)
        treeArr[whichRoot][node].likenc[site][i]=treeArr[whichRoot][node].likenc[site][i]/max;
      //  if (node ==root && site==0) printf("Node %i, Conditional likes: %lf %lf %lf %lf, ",node,tree[node].likenc[site][0],tree[node].likenc[site][1],tree[node].likenc[site][2],tree[node].likenc[site][3]);
      UFCnc[site] = UFCnc[site] + log(max);
      //  printf("max: %lf, UFC: %lf\n", max,UFCnc[site]);
    }
  }
}

void inittransitionmatrixnc(double pi[4], double par[10])

{
  int i, j;
  double sum, piT, RIVAL[4], RIVEC[4][4],  A[4][4], workspace[8];
  for (i=0; i<8; i++){
	  workspace[i]=0;
  }
  A[0][1]=pi[1]*par[4];
  A[0][2]=pi[2]*par[5];
  A[0][3]=pi[3]*par[6];
  A[1][0]=pi[0]*par[4];
  A[1][2]=pi[2]*par[7];
  A[1][3]=pi[3]*par[8];
  A[2][0]=pi[0]*par[5];
  A[2][1]=pi[1]*par[7];
  A[2][3]=pi[3]; //unscaled rate of GT = 1.0
  A[3][0]=pi[0]*par[6];
  A[3][1]=pi[1]*par[8];
  A[3][2]=pi[2]; //unscaled rate of GT = 1.0

  for (i=0; i<4; i++)
  {
    A[i][i]=0.0;
    sum=0.0;
    for (j=0; j<4; j++)
      sum = sum + A[i][j];
    A[i][i] = -sum;
  }

     /*for (i=0; i<4; i++){
       for (j=0; j<4; j++)
       printf("%lf ",A[i][j]);
       printf("\n\n");
       }*/
  if (eigen(1, A[0], 4, RRVALnc, RIVAL, RRVECnc[0], RIVEC[0], workspace) != 0)
  {
    printf("Transitions matrix did not converge or contained non-real values!\n");
    exit(-1);
  }

  for (i=0; i<4; i++)
    for (j=0; j<4; j++){
	//	printf("LRVECnc[%d][%d]=%lf\t",i,j,LRVECnc[i][j]);
		LRVECnc[i][j] = RRVECnc[i][j];
	//	printf("RRVECnc[%d][%d]=%lf\n",i,j,RRVECnc[i][j]);
		//assert(LRVECnc[i][j]>=0);
		//assert(LRVECnc[i][j]<=1);
	}
  if (matinv(RRVECnc[0],4, 4, workspace) != 0)
    printf("Could not invert matrix!\nResults may not be reliable!\n");
  //printf("after matinv back to inittransitionmatrixnc\n");
  /*for (i=0; i<4; i++){
	  for (j=0; j<4; j++){
		printf("LRVECnc[%d][%d]=%lf\t",LRVECnc[i][j],i,j);
		printf("RRVECnc[%d][%d]=%lf\n",RRVECnc[i][j],i,j);
		}
	}*/
}
/*
double getlike_gamma(double par[10]){*/
  /*GTR + gamma model
    par[1]: inverse function for piA;
    par[2]: inverse function for piC;
    par[3]: inverse function for piG;
    par[4]: AC;
    par[5]: AG;
    par[6]: AT;
    par[7]: CG;
    par[8]: CT;
    par[9]: alpha;
    unscaled rate of G<->T defined to be 1.0;
    */
/*
double stand, L, loclike, **locloglike, max, pi[4], gampar[2], d, like = 0.0;
  int i, j, k;

  COUNT2++;

  stand = 1.0+par[1]+par[2]+par[3];
  pi[0]=par[1]/stand;
  pi[1]=par[2]/stand;
  pi[2]=par[3]/stand;
  pi[3]=1.0-pi[0]-pi[1]-pi[2];
  gampar[0]=gampar[1]=par[9]; //We are setting alpha=beta to keep a constant mean to avoid identifiability issues.  This is not the same as a standard gammma.
	//printf("piA: %.2lf, piC: %.2lf, piG: %.2lf, piT: %.2lf, pAC: %.2lf, ",pi[0],pi[1],pi[2],pi[3],par[4]);
	//printf("numbase is %d\n",numbase);
	//printf("pAG: %.2lf, pAT: %.2lf, pCG: %.2lf, pCT: %.2lf, pGT: 1.00, a: %.2lf ",par[5],par[6],par[7],par[8],par[9]);
  
	//printf("numbase is %d\n",numbase);
  UFCnc = malloc((numbase)*(sizeof(double)));
  statevector = malloc(NUMCAT*(sizeof(double)));
  locloglike = malloc(numbase*(sizeof(double *)));
  for (i=0; i<numbase; i++)
    locloglike[i] = malloc(NUMCAT*(sizeof(double)));
  definegammaquantiles(NUMCAT, gampar);

  //FIX THIS IF USING MULTIPLE CATEGORIES!!!!

  //  printf("Gamma_0 set to 1.0 ");
  statevector[0]=1.0;



  inittransitionmatrixnc(pi, par);
  for (j=0; j<NUMCAT; j++){
    for (i=0; i<numbase; i++)
      UFCnc[i]=0.0;
    makeconnc(root, statevector[j]);
    for (i=0; i<numbase; i++){
      L=0.0;
      for (k=0;k<4;k++)
        L += tree[root].likenc[i][k]*pi[k];
      if (L>0.0) locloglike[i][j] = log(L) + UFCnc[i];
    }
  }
  for (i=0; i<numbase; i++){
    loclike=0.0;
    max = -100000000000.0;
    for (j=0; j<NUMCAT; j++){
      if (locloglike[i][j]>max) //underflow protection
        max=locloglike[i][j];
    }
    for (j=0; j<NUMCAT; j++){
      d=locloglike[i][j]-max;
      if (d>-100)
        loclike += exp(d);
    }
    like = like + log(loclike) + max;
  }
  free(statevector);
  for (i=0; i<numbase; i++)
    free(locloglike[i]);
  free(locloglike);
  //printf("LIKE: %lf\n",like - (double)numbase*log((double)NUMCAT));
  //printf("\n");
  free(UFCnc);
  return -like + (double)numbase*log((double)NUMCAT);
}*/
double getlike_gamma_Arr(double par[10], int whichRoot){
  /*GTR + gamma model
    par[1]: inverse function for piA;
    par[2]: inverse function for piC;
    par[3]: inverse function for piG;
    par[4]: AC;
    par[5]: AG;
    par[6]: AT;
    par[7]: CG;
    par[8]: CT;
    par[9]: alpha;
    unscaled rate of G<->T defined to be 1.0;
    */

double stand, L, loclike, **locloglike, max, pi[4], gampar[2], d, like = 0.0;
  int i, j, k;
  COUNT2++;
  stand = 1.0+par[1]+par[2]+par[3];
  pi[0]=par[1]/stand;
  pi[1]=par[2]/stand;
  pi[2]=par[3]/stand;
  pi[3]=1.0-pi[0]-pi[1]-pi[2];
  gampar[0]=gampar[1]=par[9]; //We are setting alpha=beta to keep a constant mean to avoid identifiability issues.  This is not the same as a standard gammma.
	//printf("piA: %.2lf, piC: %.2lf, piG: %.2lf, piT: %.2lf, pAC: %.2lf, ",pi[0],pi[1],pi[2],pi[3],par[4]);
	//printf("numbase is %d\n",numbase);
	//printf("pAG: %.2lf, pAT: %.2lf, pCG: %.2lf, pCT: %.2lf, pGT: 1.00, a: %.2lf ",par[5],par[6],par[7],par[8],par[9]);
  
	//printf("numbase is %d\n",numbase);
  UFCnc = malloc((numbaseArr[whichRoot])*(sizeof(double)));
  statevector = malloc(NUMCAT*(sizeof(double)));
  locloglike = malloc(numbaseArr[whichRoot]*(sizeof(double *)));
  for (i=0; i<numbaseArr[whichRoot]; i++)
    locloglike[i] = malloc(NUMCAT*(sizeof(double)));
  definegammaquantiles(NUMCAT, gampar);

  //FIX THIS IF USING MULTIPLE CATEGORIES!!!!

  //  printf("Gamma_0 set to 1.0 ");
  statevector[0]=1.0;



  inittransitionmatrixnc(pi, par);
  for (j=0; j<NUMCAT; j++){
    for (i=0; i<numbaseArr[whichRoot]; i++)
      UFCnc[i]=0.0;
    makeconnc_Arr(rootArr[whichRoot], statevector[j],whichRoot);
    for (i=0; i<numbaseArr[whichRoot]; i++){
      L=0.0;
      for (k=0;k<4;k++)
        L += treeArr[whichRoot][rootArr[whichRoot]].likenc[i][k]*pi[k];
      if (L>0.0) locloglike[i][j] = log(L) + UFCnc[i];
    }
  }
  for (i=0; i<numbaseArr[whichRoot]; i++){
    loclike=0.0;
    max = -100000000000.0;
    for (j=0; j<NUMCAT; j++){
      if (locloglike[i][j]>max) //underflow protection
        max=locloglike[i][j];
    }
    for (j=0; j<NUMCAT; j++){
      d=locloglike[i][j]-max;
      if (d>-100)
        loclike += exp(d);
    }
    like = like + log(loclike) + max;
  }
  free(statevector);
  for (i=0; i<numbaseArr[whichRoot]; i++)
    free(locloglike[i]);
  free(locloglike);
  //printf("LIKE: %lf\n",like - (double)numbase*log((double)NUMCAT));
  //printf("\n");
  free(UFCnc);
  return -like + (double)numbaseArr[whichRoot]*log((double)NUMCAT);
}
/*
double maximizelikelihoodnc_globals(double parameters[10], int precision)
{


  //optimization function starts counting at 1 so arrays have dimensionality n+1
  int i;
  double L, lowbound[10], upbound[10];

  doNRinits(10);
  for (i=1; i<10; i++){
    lowbound[i]=0.05;
    upbound[i]=20.0;
  }
  lowbound[9]=0.3;
  //L=findmax(parameters, lowbound, upbound, 10, getlike_gamma, precision);
  L=findmax(parameters, lowbound, upbound, 9, getlike_gamma, precision);
  freeNRinits(10);
  return L;
} */
double maximizelikelihoodnc_globals_Arr(double parameters[10], int precision, int whichRoot){
  //optimization function starts counting at 1 so arrays have dimensionality n+1
  int i;
  double L, lowbound[10], upbound[10];
  doNRinits(10);
  for (i=1; i<10; i++){
    lowbound[i]=0.05;
    upbound[i]=20.0;
  }
  lowbound[9]=0.3;
  //L=findmax(parameters, lowbound, upbound, 10, getlike_gamma, precision);
  L=(double) findmax_Arr(parameters, lowbound, upbound, 9, getlike_gamma_Arr, precision,whichRoot);
  freeNRinits(10);
  return L;
} 
/*double like_bl(double par[2]){

  int i, j, s, base;
  double b, p, L=0.0;

  maketransitionmatrixnc(0, par[1]);


  for (s=0; s<numbase; s++){
    p=0.0;
    if (tree[localnode].up[0]==-1)
    {
      base = seq[localnode-numspec+1][s];
	  assert(base >= 0 && base <= 4);
      if (base<4){
        for (j=0; j<4; j++){
          p += localpi[base]*PMATnc[0][base][j]*templike_nc[s][j];
           //printf("%i: %lf %lf %lf: %lf\n",j, p,localpi[j],PMATnc[0][base][j],templike_nc[s][j]);
        }
      }
      else p=1.0;
      //printf("site %i base: %i (%lf)\n",s,base,p);
    }
    else {
      for (i=0; i<4; i++){
        b=0;
        for (j=0; j<4; j++)
          b += PMATnc[0][i][j]*templike_nc[s][j];
        p += b*localpi[i]*tree[localnode].likenc[s][i];
      }
    }
    L += log(p);
  	//printf("L is %lf\n",L);
  } COUNT++;
    
   // printf("%lf\n",L); exit(-1);
  return -L;
}*/
double like_bl_Arr(double par[2], int whichRoot){

  int i, j, s, base;
  double b, p, L=0.0;

  maketransitionmatrixnc_Arr(0, par[1],whichRoot);


  for (s=0; s<numbaseArr[whichRoot]; s++){
    p=0.0;
    if (treeArr[whichRoot][localnode].up[0]==-1)
    {
      base = seqArr[whichRoot][localnode-numspecArr[whichRoot]+1][s];
	  assert(base >= 0 && base <= 4);
      if (base<4){
        for (j=0; j<4; j++){
          p += localpi[base]*PMATnc[0][base][j]*templike_nc[s][j];
           //printf("%i: %lf %lf %lf: %lf\n",j, p,localpi[j],PMATnc[0][base][j],templike_nc[s][j]);
        }
      }
      else p=1.0;
      //printf("site %i base: %i (%lf)\n",s,base,p);
    }
    else {
      for (i=0; i<4; i++){
        b=0;
        for (j=0; j<4; j++)
          b += PMATnc[0][i][j]*templike_nc[s][j];
        p += b*localpi[i]*treeArr[whichRoot][localnode].likenc[s][i];
      }
    }
    L += log(p);
  	//printf("L is %lf\n",L);
  } COUNT++;
    
   // printf("%lf\n",L); exit(-1);
  return -L;
}
/*
void maxbl_nc(int node, int parent, double pi[4], int precision)
{
  double par[2], minpar[2], maxpar[2], L;

  par[1]=tree[node].bl; //This stuff should probably be cleaned up
  minpar[1]=MINBL;
  maxpar[1]=MAXBL;
  localpi=pi;
  localnode=node;
  //printf("Node %i likelihood before: %lf (%lf)\n",node, getlike_gamma(currentestimate),tree[node].bl);
  L = findmax(par, minpar, maxpar, 1, like_bl, precision);
  //printf("Node %i likelihood after: %lf (%lf) ",node, getlike_gamma(currentestimate),tree[node].bl);
 //printf("Node %i LIKE BL: %lf\n",node, L);
  tree[node].bl = par[1];

}*/
void maxbl_nc_Arr(int node, int parent, double pi[4], int precision, int whichRoot){
  double par[2], minpar[2], maxpar[2], L;
  par[1]=treeArr[whichRoot][node].bl; /*This stuff should probably be cleaned up*/
  minpar[1]=MINBL;
  maxpar[1]=MAXBL;
  localpi=pi;
  localnode=node;
  L = findmax_Arr(par, minpar, maxpar, 1, like_bl_Arr, precision, whichRoot);
  treeArr[whichRoot][node].bl = par[1];
}

/*void recurse_estimatebracnhlengths(int node, double pi[4], int precision)
{
  int i,j, s, parent, otherb, child1, child2;
  double max, bl;

  child1 = tree[node].up[0];
  parent = tree[node].down;
  bl = tree[node].bl;

    
  if ((otherb = tree[parent].up[0])==node)
    otherb = tree[parent].up[1];
  maketransitionmatrixnc(1, tree[otherb].bl);
  for (s=0; s<numbase; s++){
    for (i=0; i<4; i++){
      max=0.0;
      if (tree[otherb].up[0]==-1){
        templike_nc[s][i] = PMATnc[1][i][seq[otherb-numspec+1][s]]; //if(s==0) printf("temp[s][%i]: %lf (b=%i, %lf), ",i,templike_nc[s][i],seq[otherb-numspec+1][s], PMATnc[1][i][seq[otherb-numspec+1][s]]);
		}
      else {
        templike_nc[s][i]=0.0;
        for (j=0; j<4; j++)
          templike_nc[s][i]=templike_nc[s][i]+tree[otherb].likenc[s][j]*PMATnc[1][i][j];
      }
      if ((templike_nc[s][i]=templike_nc[s][i]*tree[parent].posteriornc[s][i])>max)
        max=templike_nc[s][i];
      //      if(s==0 && node ==0) {printf("node %i (%i): temp[s][%i]: %lf, other/parent: %i/%i like: %lf, parentlike: %lf  ",node,s,i,templike_nc[s][i], otherb, parent,tree[otherb].likenc[s][i],tree[parent].posteriornc[s][i]);}
    }
    //for (i=0; i<4; i++)//underflow protection
      //templike_nc[s][i]=templike_nc[s][i]/max;
  }
  maxbl_nc(node,  parent, pi, precision);
  maketransitionmatrixnc(0, tree[node].bl);
  for (s=0; s<numbase; s++){
    max=0.0;
    for (i=0; i<4; i++){
      tree[node].posteriornc[s][i]=0.0;
      for (j=0; j<4; j++)
        tree[node].posteriornc[s][i] = tree[node].posteriornc[s][i] + PMATnc[0][i][j]*templike_nc[s][j];
      if (tree[node].posteriornc[s][i]>max)//more hysterical underflow protection
        max=tree[node].posteriornc[s][i];
    }
    for (i=0; i<4; i++)
      tree[node].posteriornc[s][i]=tree[node].posteriornc[s][i]/max;
  }
  if (child1>-1)
  {
    child2 = tree[node].up[1];
    recurse_estimatebracnhlengths(child1, pi, precision);
    recurse_estimatebracnhlengths(child2, pi, precision);
  }
}*/
void recurse_estimatebranchlengths_Arr(int node, double pi[4], int precision, int whichRoot){
  int i,j, s, parent, otherb, child1, child2;
  double max, bl;
  child1 = treeArr[whichRoot][node].up[0];
  parent = treeArr[whichRoot][node].down;
  bl = treeArr[whichRoot][node].bl;
  if ((otherb = treeArr[whichRoot][parent].up[0])==node)
    otherb = treeArr[whichRoot][parent].up[1];
  maketransitionmatrixnc_Arr(1, treeArr[whichRoot][otherb].bl,whichRoot);
  for (s=0; s<numbaseArr[whichRoot]; s++){
    for (i=0; i<4; i++){
      max=0.0;
      if (treeArr[whichRoot][otherb].up[0]==-1){
        templike_nc[s][i] = PMATnc[1][i][seqArr[whichRoot][otherb-numspecArr[whichRoot]+1][s]]; /*if(s==0) printf("temp[s][%i]: %lf (b=%i, %lf), ",i,templike_nc[s][i],seq[otherb-numspec+1][s], PMATnc[1][i][seq[otherb-numspec+1][s]]);*/}
      else {
        templike_nc[s][i]=0.0;
        for (j=0; j<4; j++)
          templike_nc[s][i]=templike_nc[s][i]+treeArr[whichRoot][otherb].likenc[s][j]*PMATnc[1][i][j];
      }
      if ((templike_nc[s][i]=templike_nc[s][i]*treeArr[whichRoot][parent].posteriornc[s][i])>max)
        max=templike_nc[s][i];
      //      if(s==0 && node ==0) {printf("node %i (%i): temp[s][%i]: %lf, other/parent: %i/%i like: %lf, parentlike: %lf  ",node,s,i,templike_nc[s][i], otherb, parent,tree[otherb].likenc[s][i],tree[parent].posteriornc[s][i]);}
    }
    //for (i=0; i<4; i++)//underflow protection
      //templike_nc[s][i]=templike_nc[s][i]/max;
  }
  maxbl_nc_Arr(node,  parent, pi, precision,whichRoot);
  maketransitionmatrixnc_Arr(0, treeArr[whichRoot][node].bl,whichRoot);
  for (s=0; s<numbaseArr[whichRoot]; s++){
    max=0.0;
    for (i=0; i<4; i++){
      treeArr[whichRoot][node].posteriornc[s][i]=0.0;
      for (j=0; j<4; j++)
        treeArr[whichRoot][node].posteriornc[s][i] = treeArr[whichRoot][node].posteriornc[s][i] + PMATnc[0][i][j]*templike_nc[s][j];
      if (treeArr[whichRoot][node].posteriornc[s][i]>max)//more hysterical underflow protection
        max=treeArr[whichRoot][node].posteriornc[s][i];
    }
    for (i=0; i<4; i++)
      treeArr[whichRoot][node].posteriornc[s][i]=treeArr[whichRoot][node].posteriornc[s][i]/max;
  }
  if (child1>-1)
  {
    child2 = treeArr[whichRoot][node].up[1];
    recurse_estimatebranchlengths_Arr(child1, pi, precision, whichRoot);
    recurse_estimatebranchlengths_Arr(child2, pi, precision, whichRoot);
  }
}

/*this does one pass at estimating branch lengths */
/*FORGOT ABOUT MULTIPLE CATEGORIES!!!*/
/*void estimatebracnhlengths(double par[10], int precision)
{
  int i, j, s, child1, child2;
  double stand, pi[4];

  for (i=0; i<10; i++)
    currentestimate[i]=par[i];

  doNRinits(1);
  stand = 1.0+par[1]+par[2]+par[3];
  pi[0]=par[1]/stand;
  pi[1]=par[2]/stand;
  pi[2]=par[3]/stand;
  pi[3]=1.0-pi[0]-pi[1]-pi[2];
  child1 = tree[root].up[0];
  child2 = tree[root].up[1];
  templike_nc = malloc(numbase*(sizeof(double *)));
  for (i=0; i<numbase; i++)
    templike_nc[i]=malloc(4*(sizeof(double)));

  //Root is different because only one parameter can be estimated for the two child nodes. We arbitrarily fix one of the branchlengths to equal MINBL and only estimate the other one

  if ((tree[child2].bl=tree[child2].bl+tree[child1].bl-MINBL)<MINBL)
    tree[child2].bl=MINBL;
  tree[child1].bl=MINBL;
  getlike_gamma(par); //this is not necessary if the likelihood fucntion has already been called



  //  printf("Root likelihoods: %lf %lf %lf %lf\n",tree[root].likenc[s][0],tree[root].likenc[s][1],tree[root].likenc[s][2],tree[root].likenc[s][3]);
  //  printf("Child1 likelihoods: %lf %lf %lf %lf\n",tree[child1].likenc[s][0],tree[child1].likenc[s][1],tree[child1].likenc[s][2],tree[child1].likenc[s][3]);
  //  printf("Child2 likelihoods: %lf %lf %lf %lf\n",tree[child2].likenc[s][0],tree[child2].likenc[s][1],tree[child2].likenc[s][2],tree[child2].likenc[s][3]);

  maketransitionmatrixnc(0, tree[child2].bl+tree[child1].bl);
  for (s=0; s<numbase; s++){
    for (i=0; i<4; i++){
      tree[root].posteriornc[s][i] = 1.0;//tree[child1].likenc[s][i];
        if (tree[child2].up[0]>-1){
            tree[child1].posteriornc[s][i]=0.0;
            for (j=0; j<4; j++)
                tree[child1].posteriornc[s][i] += tree[child2].likenc[s][j]*PMATnc[0][i][j];
            }
        else
                tree[child1].posteriornc[s][i] = PMATnc[0][i][seq[child2-numspec+1][s]];
            }
    //  if(s==0){
    //    printf("node %i initialization: %lf %lf %lf %lf\n",root,tree[root].posteriornc[s][0],tree[root].posteriornc[s][1],tree[root].posteriornc[s][2],tree[root].posteriornc[s][3]);
    //   printf("Initial child 1 like: %lf %lf %lf %lf\n",tree[child1].likenc[s][0],tree[child1].likenc[s][1],tree[child1].likenc[s][2],tree[child1].likenc[s][3]);
    // }
  }
  if (tree[child1].up[0]>-1)
    recurse_estimatebracnhlengths(tree[child1].up[0], pi, precision);
  if (tree[child1].up[1]>-1)
    recurse_estimatebracnhlengths(tree[child1].up[1], pi, precision);
  recurse_estimatebracnhlengths(child2, pi, precision);
  for (i=0; i<numbase; i++)
    free(templike_nc[i]);
  free(templike_nc);
  freeNRinits(1);
}*/
void estimatebranchlengths_Arr(double par[10], int precision, int whichRoot){
  int i, j, s, child1, child2;
  double stand, pi[4];
  //for (i=0; i<10; i++){
    //currentestimate[i]=par[i];
//	}
  doNRinits(1);
  stand = 1.0+par[1]+par[2]+par[3];
  pi[0]=par[1]/stand;
  pi[1]=par[2]/stand;
  pi[2]=par[3]/stand;
  pi[3]=1.0-pi[0]-pi[1]-pi[2];
  child1 = treeArr[whichRoot][rootArr[whichRoot]].up[0];
  child2 = treeArr[whichRoot][rootArr[whichRoot]].up[1];
  templike_nc = malloc(numbaseArr[whichRoot]*(sizeof(double *)));
  for (i=0; i<numbaseArr[whichRoot]; i++)
    templike_nc[i]=malloc(4*(sizeof(double)));

  /*Root is different because only one parameter can be estimated for the two child nodes. We arbitrarily fix one of the branchlengths to equal MINBL and only estimate the other one*/

  if ((treeArr[whichRoot][child2].bl=treeArr[whichRoot][child2].bl+treeArr[whichRoot][child1].bl-MINBL)<MINBL)
    treeArr[whichRoot][child2].bl=MINBL;
  treeArr[whichRoot][child1].bl=MINBL;
  getlike_gamma_Arr(par,whichRoot); /*this is not necessary if the likelihood fucntion has already been called*/
  maketransitionmatrixnc_Arr(0, treeArr[whichRoot][child2].bl+treeArr[whichRoot][child1].bl,whichRoot);
  for (s=0; s<numbaseArr[whichRoot]; s++){
    for (i=0; i<4; i++){
      treeArr[whichRoot][rootArr[whichRoot]].posteriornc[s][i] = 1.0;/*tree[child1].likenc[s][i];*/
        if (treeArr[whichRoot][child2].up[0]>-1){
            treeArr[whichRoot][child1].posteriornc[s][i]=0.0;
            for (j=0; j<4; j++)
                treeArr[whichRoot][child1].posteriornc[s][i] += treeArr[whichRoot][child2].likenc[s][j]*PMATnc[0][i][j];
            }
        else
                treeArr[whichRoot][child1].posteriornc[s][i] = PMATnc[0][i][seqArr[whichRoot][child2-numspecArr[whichRoot]+1][s]];
            }
    //  if(s==0){
    //    printf("node %i initialization: %lf %lf %lf %lf\n",root,tree[root].posteriornc[s][0],tree[root].posteriornc[s][1],tree[root].posteriornc[s][2],tree[root].posteriornc[s][3]);
    //   printf("Initial child 1 like: %lf %lf %lf %lf\n",tree[child1].likenc[s][0],tree[child1].likenc[s][1],tree[child1].likenc[s][2],tree[child1].likenc[s][3]);
    // }
  }
  if (treeArr[whichRoot][child1].up[0]>-1)
    recurse_estimatebranchlengths_Arr(treeArr[whichRoot][child1].up[0], pi, precision,whichRoot);
  if (treeArr[whichRoot][child1].up[1]>-1)
    recurse_estimatebranchlengths_Arr(treeArr[whichRoot][child1].up[1], pi, precision,whichRoot);
  recurse_estimatebranchlengths_Arr(child2, pi, precision,whichRoot);
  for (i=0; i<numbaseArr[whichRoot]; i++)
    free(templike_nc[i]);
  free(templike_nc);
  freeNRinits(1);
}

//Other schemes could be used for iterative optimizsiton.  This one seems to provide a reasobale compromise between time and accuracy
/*void estimatenucparameters(double parameters[10]){
  double L;
  COUNT=COUNT2=0;
	printf("gets in estimatenucparameters()\n");
  printf("Initial likelihoodL value = %lf\n",-getlike_gamma(parameters));
  estimatebracnhlengths(parameters, 0);
  L=maximizelikelihoodnc_globals(parameters, 0);
  printf("Current ML value = %lf\n",L);
  estimatebracnhlengths(parameters, 0);
  L=maximizelikelihoodnc_globals(parameters, 0);
  printf("Current ML value = %lf\n",L);
  estimatebracnhlengths(parameters,1);
  L=maximizelikelihoodnc_globals(parameters, 0);
  printf("Current ML value = %lf\n",L);
  estimatebracnhlengths(parameters,2);
  L=maximizelikelihoodnc_globals(parameters, 2);
  printf("Current ML value = %lf\n",L);
  estimatebracnhlengths(parameters,2);
  estimatebracnhlengths(parameters,2);
  printf("Current ML value = %lf\n",-getlike_gamma(parameters));
  printf("Number of calls to likelihood functions: %i %i\n",COUNT,COUNT2);

}*/
void estimatenucparameters_Arr(double parameters[10], int whichRoot){
  double L;
  COUNT=COUNT2=0;
  clearGlobals();
  printf("Initial likelihoodL value = %lf\n",-getlike_gamma_Arr(parameters,whichRoot));
  //estimatebracnhlengths(parameters, 0);
  estimatebranchlengths_Arr(parameters, 0, whichRoot);
  L=maximizelikelihoodnc_globals_Arr(parameters, 0,whichRoot);
  printf("Current ML value = %lf\n",L);
  estimatebranchlengths_Arr(parameters, 0, whichRoot);
  L=maximizelikelihoodnc_globals_Arr(parameters, 0, whichRoot);
  printf("Current ML value = %lf\n",L);
  estimatebranchlengths_Arr(parameters,1,whichRoot);
  L=maximizelikelihoodnc_globals_Arr(parameters, 0,whichRoot);
  printf("Current ML value = %lf\n",L);
  estimatebranchlengths_Arr(parameters,2,whichRoot);
  L=maximizelikelihoodnc_globals_Arr(parameters, 2,whichRoot);
  printf("Current ML value = %lf\n",L);
  estimatebranchlengths_Arr(parameters,2,whichRoot);
  estimatebranchlengths_Arr(parameters,2,whichRoot);
  printf("Current ML value = %lf\n",-getlike_gamma_Arr(parameters, whichRoot));
}
void clearGlobals(){
	int i,j, k;
	for(i=0;i<STATESPACE;i++){
		RRVAL[i]=0;
		for(j=0; j<STATESPACE;j++){
			LRVEC[i][j]=0;
			RRVEC[i][j]=0;
			PMAT1[i][j]=0;
			PMAT2[i][j]=0;
		}
	}
	for(i=4;i<4;i++){
		RRVALnc[i]=0;
		for(i=4;i<4;i++){
			LRVECnc[i][j]=0;
			RRVECnc[i][j]=0;
			}
	}
  	for (i=0;i<4;i++){
    	PMATnc[0][i][4]=1.0;//missing data
    	PMATnc[1][i][4]=1.0;//missing data
  	}
	/*for(i=0;i<2;i++){
		for(j=0;j<4;j++){
			for(k=0;k<5;k++){
				PMATnc[i][j][k]=0;
			}
		}
	}*/
	parameters[0]=0.0;
	for(i=1;i<10;i++){
		parameters[i]=1.0;
	}
}
			
/*
void makeposterior_nc(int node)
{
  int i,j, s, parent, otherb, child1, child2, b;
  double bl, max;


  child1 = tree[node].up[0];
  child2 = tree[node].up[1];
  parent = tree[node].down;
  bl = tree[node].bl;
  maketransitionmatrixnc(0, bl);
  if ((otherb = tree[parent].up[0])==node)
    otherb = tree[parent].up[1];
  maketransitionmatrixnc(1, tree[otherb].bl);
  for (s=0; s<numbase; s++){
    if (tree[otherb].up[0]>-1){
      for (i=0; i<4; i++){
        templike_nc[s][i]=0;
        for (j=0; j<4; j++)
          templike_nc[s][i] += tree[otherb].likenc[s][j]*PMATnc[1][i][j];
        templike_nc[s][i]=templike_nc[s][i]*tree[parent].posteriornc[s][i];
      }
    }
    else{
      b=seq[otherb-numspec+1][s];
      for (i=0; i<4; i++)
        templike_nc[s][i] = PMATnc[1][i][b]*tree[parent].posteriornc[s][i];
    }
    for (i=0; i<4; i++){
      tree[node].posteriornc[s][i]=0.0;
      max=0.0;
      for (j=0; j<4; j++){
        if ((tree[node].posteriornc[s][i] = tree[node].posteriornc[s][i] + PMATnc[0][i][j]*templike_nc[s][j])>max)
          max=tree[node].posteriornc[s][i];//more underflow protection
      }
    }
    for (i=0; i<4; i++)
      tree[node].posteriornc[s][i]=tree[node].posteriornc[s][i]/max;
  }
  //printf("In make posterior node %i: (pa=%lf,pt=%lf) (templikea&t: %lf %lf\n",node,tree[node].posteriornc[0][0],tree[node].posteriornc[0][3],templike_nc[0][0],templike_nc[0][3]);
  if (tree[child1].up[0]>-1)
    makeposterior_nc(child1);
  if (tree[child2].up[0]>-1)
    makeposterior_nc(child2);
}*/
void makeposterior_nc_Arr(int node, int whichRoot)
{
  int i,j, s, parent, otherb, child1, child2, b;
  double bl, max;


  child1 = treeArr[whichRoot][node].up[0];
  child2 = treeArr[whichRoot][node].up[1];
  parent = treeArr[whichRoot][node].down;
  bl = treeArr[whichRoot][node].bl;
  maketransitionmatrixnc_Arr(0, bl,whichRoot);
  if ((otherb = treeArr[whichRoot][parent].up[0])==node)
    otherb = treeArr[whichRoot][parent].up[1];
  maketransitionmatrixnc_Arr(1, treeArr[whichRoot][otherb].bl,whichRoot);
  for (s=0; s<numbaseArr[whichRoot]; s++){
    if (treeArr[whichRoot][otherb].up[0]>-1 || strcmp(treeArr[whichRoot][otherb].name,"OU061397_1")==0){
      for (i=0; i<4; i++){
        templike_nc[s][i]=0;
        for (j=0; j<4; j++)
          templike_nc[s][i] += treeArr[whichRoot][otherb].likenc[s][j]*PMATnc[1][i][j];
        templike_nc[s][i]=templike_nc[s][i]*treeArr[whichRoot][parent].posteriornc[s][i];
      }
    }
    else{
      b=seqArr[whichRoot][otherb-numspecArr[whichRoot]+1][s];
      for (i=0; i<4; i++)
        templike_nc[s][i] = PMATnc[1][i][b]*treeArr[whichRoot][parent].posteriornc[s][i];
    }
    for (i=0; i<4; i++){
      treeArr[whichRoot][node].posteriornc[s][i]=0.0;
      max=0.0;
      for (j=0; j<4; j++){
        if ((treeArr[whichRoot][node].posteriornc[s][i] = treeArr[whichRoot][node].posteriornc[s][i] + PMATnc[0][i][j]*templike_nc[s][j])>max)
          max=treeArr[whichRoot][node].posteriornc[s][i];//more underflow protection
      }
    }
    for (i=0; i<4; i++)
      treeArr[whichRoot][node].posteriornc[s][i]=treeArr[whichRoot][node].posteriornc[s][i]/max;
  }
  //printf("In make posterior node %i: (pa=%lf,pt=%lf) (templikea&t: %lf %lf\n",node,tree[node].posteriornc[0][0],tree[node].posteriornc[0][3],templike_nc[0][0],templike_nc[0][3]);
  if (child1==-1){ return;}
  if (treeArr[whichRoot][child1].up[0]>-1)
    makeposterior_nc_Arr(child1,whichRoot);
  if (treeArr[whichRoot][child2].up[0]>-1)
    makeposterior_nc_Arr(child2,whichRoot);
}
/*
//THIS FUNCTION IS IMPLEMENTED WITHOUT GAMMA CATEGORIES
void getposterior_nc(double par[10])
{
  int i, j, s, k, parent, b, notdonebefore;
  double p, sum, pi[4], stand, **templike;

   // printf("AA\n");

  getlike_gamma(par); //need to call likelihood again
  templike_nc = malloc(numbase*(sizeof(double *)));
  for (i=0; i<numbase; i++)
    templike_nc[i]=malloc(4*(sizeof(double)));
  stand = 1.0+par[1]+par[2]+par[3];
  pi[0]=par[1]/stand;
  pi[1]=par[2]/stand;
  pi[2]=par[3]/stand;
  pi[3]=1.0-pi[0]-pi[1]-pi[2];
  for (s=0; s<numbase; s++)
    for (i=0; i<4; i++)
      tree[root].posteriornc[s][i] = 1.0;
  if (tree[tree[root].up[0]].up[0]>-1) makeposterior_nc(tree[root].up[0]);
  if (tree[tree[root].up[1]].up[0]>-1) makeposterior_nc(tree[root].up[1]);
  for (j=0; j<2*numspec-1; j++){
    if (tree[j].up[0]>-1){
      //printf("Node %i: L1=(%lf, %lf, %lf, %lf), alike2=%lf, ",j,tree[j].likenc[0][0],tree[j].likenc[0][1],tree[j].likenc[0][2],tree[j].likenc[0][3],tree[j].posteriornc[0][0]);
      for (s=0; s<numbase; s++) {
        sum = 0.0;
        for (i=0; i<4; i++)
          sum = sum + (tree[j].posteriornc[s][i]=tree[j].likenc[s][i]*tree[j].posteriornc[s][i]*pi[i]);
        for (i=0; i<4; i++)
          tree[j].posteriornc[s][i] = tree[j].posteriornc[s][i]/sum;
      }
     // printf("posterior(a)=%lf\n",tree[j].posteriornc[0][0]);
    }
    else { //printf("B\n");
      for (s=0; s<numbase; s++) {
        notdonebefore=1;
        b=seq[j-numspec+1][s];
        if (b==4){
          if (notdonebefore==1) {maketransitionmatrixnc(0, tree[j].bl); notdonebefore=0;}
          parent=tree[j].down;
          sum = 0.0;
          for (i=0; i<4; i++){
            tree[j].posteriornc[s][i]=0.0;
            for (k=0; k<4; k++)
              tree[j].posteriornc[s][i] += tree[parent].posteriornc[s][k]*PMATnc[0][i][k];
            sum = sum + (tree[j].posteriornc[s][i]=tree[j].posteriornc[s][i]*pi[i]);
          }
          for (i=0; i<4; i++)
            tree[j].posteriornc[s][i] = tree[j].posteriornc[s][i]/sum;
        }
        else {
          for (i=0; i<4; i++){
            if (i==b) tree[j].posteriornc[s][i]=1.0;
            else tree[j].posteriornc[s][i]=0.0;
          }
        }
      }
  //  printf("Node %i: aposterior=%lf, base=%i\n",j,tree[j].posteriornc[0][0],seq[j-numspec+1][0]);
    }
  }
  for (i=0; i<numbase; i++)
    free(templike_nc[i]);
  free(templike_nc);

}*/
void getposterior_nc_Arr(double par[10], int whichRoot){
  int i, j, s, k, parent, b, notdonebefore;
  double p, sum, pi[4], stand, **templike;
   // printf("AA\n");
  getlike_gamma_Arr(par,whichRoot); /*need to call likelihood again*/
  templike_nc = malloc(numbaseArr[whichRoot]*(sizeof(double *)));
  for (i=0; i<numbaseArr[whichRoot]; i++)
    templike_nc[i]=malloc(4*(sizeof(double)));
  stand = 1.0+par[1]+par[2]+par[3];
  pi[0]=par[1]/stand;
  pi[1]=par[2]/stand;
  pi[2]=par[3]/stand;
  pi[3]=1.0-pi[0]-pi[1]-pi[2];
  for (s=0; s<numbaseArr[whichRoot]; s++)
    for (i=0; i<4; i++)
      treeArr[whichRoot][rootArr[whichRoot]].posteriornc[s][i] = 1.0;
  if (treeArr[whichRoot][treeArr[whichRoot][rootArr[whichRoot]].up[0]].up[0]>-1) makeposterior_nc_Arr(treeArr[whichRoot][rootArr[whichRoot]].up[0],whichRoot);
  if (treeArr[whichRoot][treeArr[whichRoot][rootArr[whichRoot]].up[1]].up[0]>-1) makeposterior_nc_Arr(treeArr[whichRoot][rootArr[whichRoot]].up[1],whichRoot);
  for (j=0; j<2*numspecArr[whichRoot]-1; j++){
    if (treeArr[whichRoot][j].up[0]>-1){
      //printf("Node %i: L1=(%lf, %lf, %lf, %lf), alike2=%lf, ",j,tree[j].likenc[0][0],tree[j].likenc[0][1],tree[j].likenc[0][2],tree[j].likenc[0][3],tree[j].posteriornc[0][0]);
      for (s=0; s<numbaseArr[whichRoot]; s++) {
        sum = 0.0;
        for (i=0; i<4; i++)
          sum = sum + (treeArr[whichRoot][j].posteriornc[s][i]=treeArr[whichRoot][j].likenc[s][i]*treeArr[whichRoot][j].posteriornc[s][i]*pi[i]);
        for (i=0; i<4; i++)
          treeArr[whichRoot][j].posteriornc[s][i] = treeArr[whichRoot][j].posteriornc[s][i]/sum;
      }
     // printf("posterior(a)=%lf\n",tree[j].posteriornc[0][0]);
    }
    else { //printf("B\n");
      for (s=0; s<numbaseArr[whichRoot]; s++) {
        notdonebefore=1;
        b=seqArr[whichRoot][j-numspecArr[whichRoot]+1][s];
        if (b==4 || j==200){
          if (notdonebefore==1) {maketransitionmatrixnc_Arr(0, treeArr[whichRoot][j].bl,whichRoot); notdonebefore=0;}
          parent=treeArr[whichRoot][j].down;
          sum = 0.0;
          for (i=0; i<4; i++){
            treeArr[whichRoot][j].posteriornc[s][i]=0.0;
            for (k=0; k<4; k++)
              treeArr[whichRoot][j].posteriornc[s][i] += treeArr[whichRoot][parent].posteriornc[s][k]*PMATnc[0][i][k];
            sum = sum + (treeArr[whichRoot][j].posteriornc[s][i]=treeArr[whichRoot][j].posteriornc[s][i]*pi[i]);
          }
          for (i=0; i<4; i++)
            treeArr[whichRoot][j].posteriornc[s][i] = treeArr[whichRoot][j].posteriornc[s][i]/sum;
        }
        else {
          for (i=0; i<4; i++){
            if (i==b) treeArr[whichRoot][j].posteriornc[s][i]=1.0;
            else treeArr[whichRoot][j].posteriornc[s][i]=0.0;
          }
        }
      }
  //  printf("Node %i: aposterior=%lf, base=%i\n",j,tree[j].posteriornc[0][0],seq[j-numspec+1][0]);
    }
  }
  for (i=0; i<numbaseArr[whichRoot]; i++)
    free(templike_nc[i]);
  free(templike_nc);

}
