/* file: rrgen.c
     Entry number  : 201
     Authors       : P. E. McSharry, G. D. Clifford
     Organization  : Dept Maths & Dept Engineering, University of Oxford, UK
     email address : mcsharry at maths.ox.ac.uk, gari at robots.ox.ac.uk

   This program should compile without errors or warnings using:
	gcc -Wall rrgen.c -lm

   See http://www.physionet.org/challenge/2002/ for further information on
   the CinC Challenge 2002.

   This program was used to generate series rr19 and rr31 of the challenge
   dataset.
*/

/* NOTE - some algorithms from numerical recipes in C are used. 
   If you have not done so, please purchase a set of disks from 
   the appropriate source before copying/running this code */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define	SPS	128	/* sampling frequency (determines quantization of
			   RR intervals;  see main()) */
#define SWAP(a,b) tempr=(a);(a)=(b);(b)=tempr
#define MIN(a,b) (a < b ? a : b)
#define MAX(a,b) (a > b ? a : b)
#define PI (2.0*asin(1.0))
#define NR_END 1
#define FREE_ARG char*

#define IA 16807
#define IM 2147483647
#define AM (1.0/IM)
#define IQ 127773
#define IR 2836
#define NTAB 32
#define NDIV (1+(IM-1)/NTAB)
#define EPS 1.2e-7
#define RNMX (1.0-EPS)

long rseed;

#define AWAKE  1.0
#define ASLEEP 0.0
#define YES    1
#define NO     0

int first_pass  = 0;
float last_rr = 0.0; /* intialise the first rr interval to be zero */
float last_Rt = 0.0; /* intialise the first R-peak time stamp to be zero */

float *trpeaks_new;
long N; /* number of seconds */
long count = 1; /* global counter */

/*---------------------------------------------------------------------------*/
/* modify beat for ectopy                                                    */
/*---------------------------------------------------------------------------*/

float modify_ectopic(float rr2, float rr1)
{ /* shorten beat by 80% +/- 10% randomly chosen */
  float new_rr;
  float delta;
  float rrint;

  delta = ((float)rand()/RAND_MAX);
  delta = 0.2*(delta-0.5);
  /* length of last rr interval */
  rrint = rr2 - rr1;
  /* shorten it by 0.7 +/-  a bit */
  rrint = (0.7 - delta)*rrint; 
  /* work out new timing */
  new_rr = rr1 + rrint;

  return(new_rr); 
}

/*---------------------------------------------------------------------------*/
/* generate an extra beat (artefact)                                         */
/*---------------------------------------------------------------------------*/

float make_noise(float rr2, float rr1)
{ /* shorten beat by 80% +/- 10% randomly chosen */
  float new_rr;
  float delta;
  float rrint;

  delta = ((float)rand()/RAND_MAX);
  delta = ((delta*(SPS-(2.0/SPS))) + (1.0/SPS))/SPS; 

  /* length of last rr interval */
  rrint = rr2 - rr1;
  /* shorten it by 0.5 +/-  a bit */
  rrint = delta*rrint; 
  /* work out new timing */
  new_rr = rr1 + rrint;

  return(new_rr); 
}

/*--------------------------------------------------------------------------*/
/* CHECK IF ECTOPICS SHOULD BE ADDED                                        */ 
/*--------------------------------------------------------------------------*/

int check_for_ectopic(float Pe)
{
  /* float Pe=0.0003;*/ /* probability of ectopic (~1 per hour)*/

  /* no state or hr dependency */
  int rtn=NO; 
  float rand_no;

  /* check to see if we really want ectopic beats */
  if(Pe==0) return(NO);

  /* if we do ... */

  /* pick a random number between 0 and 1*/
  rand_no =  ((float)rand()/RAND_MAX);

  if(rand_no<Pe)    rtn = YES;
  else              rtn = NO;

  return(rtn);
}

/*---------------------------------------------------------------------------*/
/*      CHECK IF ARTEFACT SHOULD BE ADDED                                    */
/*---------------------------------------------------------------------------*/

int check_for_noise(int state, float hr, float Pn)
{
  float rand_no;
  int   rtn =NO; 
  float Pa0 = 0.0004; /* probability of artefact state 0 - ASLEEP*/
  float Pa1 = 0.0048; /* probability of artefact state 1 - AWAKE */

  Pa1 = Pn; 
  Pa0 = Pn/12;  /* Probability of noise in sleep is 12 times less than when awake */

  /* check to see if we really want to add noise beats */
  if(Pn==0) return(NO);

  /* if we do ... */

  /* pick a random number */
  rand_no =  ((float)rand()/RAND_MAX);

  /* if asleep, flag an artefact randomly (scaling using heart rate) */
  if(state==ASLEEP){
    if((rand_no*60/hr)<Pa0)    rtn = YES;
    else                       rtn = NO;
  }

  /* if awake, flag an artefact randomly with higher prob*/
  if(state!=ASLEEP){
    if((rand_no*60/hr)<Pa1)    rtn = YES;
    else                       rtn = NO;
  }

  return(rtn);
}

/*---------------------------------------------------------------------------*/
/*      ALLOCATE MEMORY FOR VECTOR                                           */
/*---------------------------------------------------------------------------*/
 
float *vector(long nl, long nh)
/* allocate a float vector with subscript range v[nl..nh] */
{
        float *v;
 
        v=(float *)malloc((size_t) ((nh-nl+1+NR_END)*sizeof(float)));
        if (!v) printf("allocation failure in vector");
        return v-nl+NR_END;
}

/*---------------------------------------------------------------------------*/
/*      FREE MEMORY FOR VECTOR                                               */
/*---------------------------------------------------------------------------*/
 
void free_vector(float *v, long nl, long nh)
/* free a float vector allocated with vector() */
{
        free((FREE_ARG) (v+nl-NR_END));
}

/*---------------------------------------------------------------------------*/
/*      MEAN CALCULATOR                                                      */
/*---------------------------------------------------------------------------*/
 
float mean(float *x, int n)
/* n-by-1 vector, calculate mean */
{
        int j;
        float add;
 
        add = 0.0;
        for(j=1;j<=n;j++)  add += x[j];
 
        return (add/n);
}

/*---------------------------------------------------------------------------*/
/*      STANDARD DEVIATION CALCULATOR                                        */
/*---------------------------------------------------------------------------*/
 
float std(float *x, int n)
/* n-by-1 vector, calculate standard deviation */
{
        int j;
        float add,mean,diff,total;
 
        add = 0.0;
        for(j=1;j<=n;j++)  add += x[j];
        mean = add/n;
 
        total = 0.0;
        for(j=1;j<=n;j++)
        {
           diff = x[j] - mean;
           total += diff*diff;
        }
 
        return (sqrt(total/(n-1)));
}

/*--------------------------------------------------------------------------*/
/*    UNIFORM DEVIATES                                                      */
/*--------------------------------------------------------------------------*/

float ran1(long *idum)
{
	int j;
	long k;
	static long iy=0;
	static long iv[NTAB];
	float temp;

	if (*idum <= 0 || !iy) {
		if (-(*idum) < 1) *idum=1;
		else *idum = -(*idum);
		for (j=NTAB+7;j>=0;j--) {
			k=(*idum)/IQ;
			*idum=IA*(*idum-k*IQ)-IR*k;
			if (*idum < 0) *idum += IM;
			if (j < NTAB) iv[j] = *idum;
		}
		iy=iv[0];
	}
	k=(*idum)/IQ;
	*idum=IA*(*idum-k*IQ)-IR*k;
	if (*idum < 0) *idum += IM;
	j=iy/NDIV;
	iy=iv[j];
	iv[j] = *idum;
	if ((temp=AM*iy) > RNMX) return RNMX;
	else return temp;
}

/*--------------------------------------------------------------------------*/
/*    GAUSSIAN DEVIATES                                                     */
/*--------------------------------------------------------------------------*/

float gasdev(long *idum)
{
	float ran1(long *idum);
	static int iset=0;
	static float gset;
	float fac,rsq,v1,v2;

	if  (iset == 0) {
		do {
			v1=2.0*ran1(idum)-1.0;
			v2=2.0*ran1(idum)-1.0;
			rsq=v1*v1+v2*v2;
		} while (rsq >= 1.0 || rsq == 0.0);
		fac=sqrt(-2.0*log(rsq)/rsq);
		gset=v1*fac;
		iset=1;
		return v2*fac;
	} else {
		iset=0;
		return gset;
	}
}


/*--------------------------------------------------------------------------*/
/*    FFT                                                                   */
/*--------------------------------------------------------------------------*/

void four1(float data[], unsigned long nn, int isign)
{
	unsigned long n,mmax,m,j,istep,i;
	double wtemp,wr,wpr,wpi,wi,theta;
	float tempr,tempi;

	n=nn << 1;
	j=1;
	for (i=1;i<n;i+=2) {
		if (j > i) {
			SWAP(data[j],data[i]);
			SWAP(data[j+1],data[i+1]);
		}
		m=n >> 1;
		while (m >= 2 && j > m) {
			j -= m;
			m >>= 1;
		}
		j += m;
	}
	mmax=2;
	while (n > mmax) {
		istep=mmax << 1;
		theta=isign*(6.28318530717959/mmax);
		wtemp=sin(0.5*theta);
		wpr = -2.0*wtemp*wtemp;
		wpi=sin(theta);
		wr=1.0;
		wi=0.0;
		for (m=1;m<mmax;m+=2) {
			for (i=m;i<=n;i+=istep) {
				j=i+mmax;
				tempr=wr*data[j]-wi*data[j+1];
				tempi=wr*data[j+1]+wi*data[j];
				data[j]=data[i]-tempr;
				data[j+1]=data[i+1]-tempi;
				data[i] += tempr;
				data[i+1] += tempi;
			}
			wr=(wtemp=wr)*wpr-wi*wpi+wr;
			wi=wi*wpr+wtemp*wpi+wi;
		}
		mmax=istep;
	}
}

/*--------------------------------------------------------------------------*/
/*    BLOCK LENGTH                                                          */
/*--------------------------------------------------------------------------*/

int blocklength()
{
   long L;
   double alpha,beta,u;

   alpha = 5466.8;
   beta = - 2.2;

   u = (double)rand()/RAND_MAX;
   while(u < 1e-3)  u = (double)rand()/RAND_MAX;
   
   L = (int)pow( u/alpha, 1.0/beta);
 
   return L;
}

/*--------------------------------------------------------------------------*/
/*    TRANS LENGTH                                                          */
/*--------------------------------------------------------------------------*/

int translength()
{
   int L;
   float u;

   u = (float)rand()/RAND_MAX;
   L = 5 + 25*u;

   return L;
}

/*--------------------------------------------------------------------------*/
/*    INTERP                                                                */
/*--------------------------------------------------------------------------*/

void interp(float *y, float *x, int n, int r)
{
   int i,j;
   float a;

   for(i=1;i<=n-1;i++)
   {
      for(j=1;j<=r;j++) 
      {
         a = (j-1)*1.0/r;
         y[(i-1)*r+j] = (1.0-a)*x[i] + a*x[i+1];
      }
   }
   

}

/*--------------------------------------------------------------------------*/
/*    GENERATE TR BLOCK                                                     */
/*--------------------------------------------------------------------------*/

void trblock(float *tr, float flo, float fhi, 
float flostd, float fhistd, float lfhfratio,  
float rrmean, float rrtrend, float rrstd, int nbeats)
{
   int i,j,n;
   float c1,c2,w1,w2,sig1,sig2,xstd,ratio;
   float sfrr,dtrr;
   float df,dw1,dw2,*w,*Hw,*Sw,*ph0,*ph,*SwC,*x,*rr;

   sfrr = 1.0;
   dtrr = 1.0;
   n = (int)pow(2.0, ceil(log(1.2*nbeats*rrmean*sfrr)/log(2.0))); 
   if(n < 256) n = 256;

   w = vector(1,n);
   Hw = vector(1,n);
   Sw = vector(1,n);
   ph0 = vector(1,n/2-1);
   ph = vector(1,n);
   SwC = vector(1,2*n);
   x = vector(1,n);
   rr = vector(1,n*10);

   w1 = 2.0*PI*flo;
   w2 = 2.0*PI*fhi;

   c1 = 2.0*PI*flostd;
   c2 = 2.0*PI*fhistd;
   sig2 = 1.0;
   sig1 = lfhfratio;


   df = sfrr/n;
   for(i=1;i<=n;i++) w[i] = (i-1)*2.0*PI*df;
   for(i=1;i<=n;i++) 
   {
      dw1 = w[i]-w1;
      dw2 = w[i]-w2;
      Hw[i] = sig1*exp(-dw1*dw1/(2.0*c1*c1))/sqrt(2*PI*c1*c1) 
            + sig2*exp(-dw2*dw2/(2.0*c2*c2))/sqrt(2*PI*c2*c2); 
   }

   for(i=1;i<=n/2;i++) Sw[i] = (sfrr/2.0)*Hw[i];
   for(i=n/2+1;i<=n;i++) Sw[i] = (sfrr/2.0)*Hw[n-i+1];

   /* randomise the phases */
   for(i=1;i<=n/2-1;i++) ph0[i] = 2.0*PI*(float)rand()/RAND_MAX;
   ph[1] = 0.0;
   for(i=1;i<=n/2-1;i++) ph[i+1] = ph0[i];
   ph[n/2+1] = 0.0;
   for(i=1;i<=n/2-1;i++) ph[n-i+1] = - ph0[i]; 



   /* make complex spectrum */
   for(i=1;i<=n;i++) SwC[2*i-1] = Sw[i]*cos(ph[i]);
   for(i=1;i<=n;i++) SwC[2*i] = Sw[i]*sin(ph[i]);


   /* calculate inverse fft */
   four1(SwC,n,-1);

   /* extract real part */
   for(i=1;i<=n;i++) x[i] = (1.0/n)*SwC[2*i-1];


   xstd = std(x,n);
   ratio = rrstd/xstd; 

   for(i=1;i<=n;i++) x[i] *= ratio;
   for(i=1;i<=n;i++) x[i] += rrmean;

   /* upsample */
   interp(rr, x, n, 10);

   /* incorporate trend */
   for(i=1;i<=n*10;i++) rr[i] = rr[i] + 0.1*rrtrend*(-0.5 + 1.0*i/(n*10));

   /* add noise */
   for(i=1;i<=n*10;i++) rr[i] += (rr[i]/10.0)*(float)rand()/RAND_MAX;

   /* determine r-peak times */
   tr[1] = rr[1];
   for(i=2;i<=nbeats;i++)
   {
      j = (int)ceil( (tr[i-1] + rr[i])/0.1); 
      tr[i] = tr[i-1] + rr[j];
   }


   free_vector(w,1,n);
   free_vector(Hw,1,n);
   free_vector(Sw,1,n);
   free_vector(ph0,1,n/2-1);
   free_vector(ph,1,n);
   free_vector(SwC,1,2*n);
   free_vector(x,1,n);
   free_vector(rr,1,n*10);
}


/*--------------------------------------------------------------------------*/
/*    GENERATE TRANSIENT                                                    */
/*--------------------------------------------------------------------------*/

void trtrans(float *tr, int n, float rr1, float rr2, float *trb)
{
   int i,iturn;
   float rr12,rrthrs,rrdrop,rrmin,a,b,*rr,*rrb,*rrdev,rrbmean;

   rr = vector(1,n);
   rrb = vector(1,n);
   rrdev = vector(1,n);

   rr12 = MIN(rr1,rr2);


   rrthrs = 0.4 + 0.2*(float)rand()/RAND_MAX;
   rrdrop = 0.1 + 0.1*(float)rand()/RAND_MAX;
   rrmin = MAX(rrthrs, rr12-rrdrop); 

   /* choose turning point at U(0.3,0.6) x (n+1) */
   iturn = (int)ceil((0.3 + 0.3*(float)rand()/RAND_MAX)*(n+1)); 
   
   /* create triangle using slopes a and b */
   a = (rrmin - rr1)/(iturn - 0.0);
   b = (rr2 - rrmin)/(n+1 - iturn);

   for(i=1;i<=iturn;i++) rr[i] = rr1 + a*i;
   for(i=iturn+1;i<=n;i++) rr[i] = rrmin + b*(i-iturn);

   /* make rr dynamics from block mean zero */
   rrb[1] = trb[1];
   for(i=2;i<=n;i++) rrb[i] = trb[i] - trb[i-1];
   rrbmean = mean(rrb, n);
   for(i=1;i<=n;i++) rrdev[i] = rrb[i] - rrbmean;

   tr[1] = rr[1] + rrdev[1];
   for(i=2;i<=n;i++) tr[i] = tr[i-1] + rr[i] + rrdev[i];

   free_vector(rr,1,n);
   free_vector(rrb,1,n);
   free_vector(rrdev,1,n);
}


/* This function is called once, before generate() is called.  See main()
   for a description of its arguments.
*/

void initialize(long seed, long tmax, float Pe, float Pn)
{
  /*  srand((unsigned int)seed); 
      return; */
   
   long i,j,time,Nb,Nt,Ndays,Nbeats;
   float RRmean,RRamp,dRRamp,RRstdmin,RRstdmax,RRtrendamp,LFHFmin,LFHFmax;
   float Tc,dRR,rrmean,rrstd,rrtrend,lfhfratio,rra,rrb;
   long Tsleepstart,Tsleepdur,Tsleepfinish;
   float flo,fhi,flostd,fhistd,sfrr,dtrr;
   float *RR0,*tr0,*tr1,*trb,*trpeaks,u,RRsleep,acir;
   float T1, T2, T3, T4;

   float hr, *state;
   int   test;

   /*   printf("seed = %ld\n",seed); */
   /*   printf("tmax = %ld\n",tmax); */

   /* set size in seconds and intialise global vector of r peaks */
   N = tmax;
   trpeaks_new = vector(1,2*N);


   /* this memory allocation works for total averaged HR <= 120 bpm */
   RR0 = vector(1,2*N);
   trpeaks = vector(1,2*N);
   tr0 = vector(1,2*N);
   tr1 = vector(1,2*N);
   trb = vector(1,2*N);
   state = vector(1,2*N);

   rseed = -seed;

   /* define RR long term characteristics */
   RRmean = 0.7 + 0.3*(float)rand()/RAND_MAX;
   RRamp = 0.075 + 0.2*(float)rand()/RAND_MAX;
   RRsleep = 0.1 + 0.1*(float)rand()/RAND_MAX;
   dRRamp = 0.03 + 0.1*(float)rand()/RAND_MAX;
   RRstdmin = 0.01;
   RRstdmax = 0.02;
   RRtrendamp = 1.0 + 0.25*(float)rand()/RAND_MAX;
   LFHFmin = 0.5;
   LFHFmax = 8.0;

    /* define characteristics of rr power spectrum */
   flo = 0.1;
   fhi = 0.25;
   flostd = 0.01;
   fhistd = 0.01;
   sfrr = 1;
   dtrr = 1/sfrr;


   /* define mean RR-interval using circadium rhythm period */
   Tc = (24.0 + gasdev(&rseed))*3600.0;
   acir = (float)rand()/RAND_MAX;
   T1 = 10.0 + 3.0*(float)rand()/RAND_MAX;
   T2 = 7.0 + 2.0*(float)rand()/RAND_MAX;
   T3 = 4.0 + 2.0*(float)rand()/RAND_MAX;
   T4 = 2.0 + 2.0*(float)rand()/RAND_MAX;
   for(i=1;i<=N;i++){

     RR0[i] = RRmean + RRamp*sin(PI + (2.0*PI/Tc)*i) 
              + 0.25*acir*RRamp*sin((2.0*PI/(Tc/T1))*i) 
              + 0.25*RRamp*(1.0 - acir)*cos((2.0*PI/(Tc/T2))*i)
              + 0.25*acir*RRamp*sin((2.0*PI/(Tc/T3))*i) 
              + 0.25*RRamp*(1.0 - acir)*cos((2.0*PI/(Tc/T4))*i);

     state[i]=AWAKE;
   }


   /* find number of 24 hour periods and insert sleep stages */
   Ndays = (long)ceil(N/(24*3600));
   /* printf("Ndays = %ld\n",Ndays); */
   for(j=1;j<=Ndays;j++)
   {

      /* define sleep state Tsleepstart ~ U(14,16) Tsleepdur ~ U(6,8) */
      Tsleepstart = (j-1)*24*3600 + (long)((14.0 + 2.0*(float)rand()/RAND_MAX)*3600.0);
      Tsleepdur = (long)((6.0 + 2.0*(float)rand()/RAND_MAX)*3600.0);    
      Tsleepfinish = Tsleepstart + Tsleepdur;
      Tsleepfinish = MIN(Tsleepfinish, N);
      
      /*     printf("Tsleepstart = %ld\n",Tsleepstart);  */
      /*      printf("Tsleepfinish = %ld\n",Tsleepfinish); */
      for(i=Tsleepstart;i<=Tsleepfinish;i++){
	RR0[i] = RRmean + 0.5*RRsleep*(1.0+sin(2*PI/(20*60) )) ;
	/* set sleep state */
	state[i] = ASLEEP;
      }
   }


   /* add noise to mean RR value */
   for(i=1;i<=N;i++) RR0[i] += 0.2*RRamp*gasdev(&rseed);


   /* initialise count in time and beats */
   time = 1;
   /* i indexes beats rather than seconds now ! */
   i = 1;
 
   /* generate first block */
   Nb = blocklength();


   /* choose rr generation states: rrmean, rrstd, lfhfratio */
   /* deviation from RR0 is given by dRR ~ U(0, dRRamp); */
   /* dRR = dRRamp*(float)(2.0*rand()/RAND_MAX - 1.0); */
   
   /* mimic segment switching in PRL 87, p168105 */
   u = (float)rand()/RAND_MAX; 
   dRR = 0.75*dRRamp*(1.0 + exp(gasdev(&rseed))/10.0)*u/fabs(u);
   while(fabs(dRR) > dRRamp)
   {
      u = (float)rand()/RAND_MAX;
      dRR = 0.5*dRRamp*(1.0 + (float)exp(gasdev(&rseed))/10.0)*u/fabs(u);
   }
      
   rrmean = RR0[time] + dRR;
   rrtrend = RRtrendamp*(2*(float)rand()/RAND_MAX - 1); 
   rrstd = RRstdmin + (RRstdmax - RRstdmin)*(float)rand()/RAND_MAX;
   lfhfratio = LFHFmin + (LFHFmax-LFHFmin)*(float)rand()/RAND_MAX;


   /* generate RR over this block */
   trblock(tr0, flo, fhi, flostd, fhistd, lfhfratio, rrmean, rrtrend, rrstd, Nb);

   for(j=1;j<=Nb;j++) trpeaks[i+j] = trpeaks[i] + tr0[j];
   i += Nb; 

   /* convert beat time stamp back to nearest second */
   time = (long)(trpeaks[i]); 
   rra = tr0[Nb] - tr0[Nb-1];

   /* printf("time = %ld\n",time); */

   /*
   for(j=2;j<=Nb;j++) printf("%d %e\n",j,tr0[j]- tr0[j-1]);
   exit(1); 
   */

   /* generate record */
   while(time < N)
   {
    
      fhi = 0.2 + 0.1*(float)rand()/RAND_MAX;
      
      /* choose a new block time in beats */
      Nb = blocklength();
      
      /* choose a transition duration in beats */
      Nt = translength();

      /* choose rr generation states: rrmean, rrstd, lfhfratio */
      dRR = dRRamp*(float)rand()/RAND_MAX;
      rrmean = RR0[time] + dRR;
      rrtrend = RRtrendamp*(2.0*(float)rand()/RAND_MAX - 1.0);
      rrstd = RRstdmin + (RRstdmax - RRstdmin)*(float)rand()/RAND_MAX;
      lfhfratio = LFHFmin + (LFHFmax - LFHFmin)*(float)rand()/RAND_MAX;

      /* generate RR over this block */
    trblock(tr0, flo, fhi, flostd, fhistd, lfhfratio, rrmean, rrtrend, rrstd, Nb);
      rrb = tr0[2] - tr0[1];

     /* generate RR during transition between blcoks with rr-intervals rra, rrb */
      rrtrend = 0.0;
    trblock(trb, flo, fhi, flostd, fhistd, lfhfratio, rrmean, rrtrend, rrstd, Nt);
    trtrans(tr1, Nt, rra, rrb, trb);


     /* fill in transition and new block */
     for(j=1;j<=Nt;j++) trpeaks[i+j] = trpeaks[i] + tr1[j];
     i += Nt;
     for(j=1;j<=Nb;j++) trpeaks[i+j] = trpeaks[i] + tr0[j];
     i += Nb;

     /* remember last rr-interval */
     rra = rrb;

     time = (long)(trpeaks[i]);

   }
   Nbeats = i;


   /* check for adding ectopics or artefacts */
   /* loop over entire data set */
   j=1;
   trpeaks_new[1] = trpeaks[1];
   for(i=2;i<=Nbeats;i++){
     /* find out what heart rate is at a particular rpeak time */
     hr=60.0/RR0[(long)(trpeaks[i])];     

     /* check so see if we (randomly) want to make this ectopic */
     test = check_for_ectopic(Pe); /* no state dependency */
     if (test!=NO){
       trpeaks_new[j] = modify_ectopic(trpeaks[i],trpeaks[i-1]);
     }    
     if (test == NO){ /* just keep the current beat */
       test = check_for_noise(state[i], hr, Pn); /* hr and state dependent */
       if (test!=NO){ /* we have noise, add an extra beat */
	 trpeaks_new[j] = make_noise(trpeaks[i],trpeaks[i-1]);
	 j++;
       }
       /* regardless, keep the current beat */
       trpeaks_new[j] = trpeaks[i];
     }
     j++;
   }

   free_vector(RR0,1,2*N);
   free_vector(trpeaks,1,2*N);
   free_vector(tr0,1,2*N);
   free_vector(tr1,1,2*N);
   free_vector(trb,1,2*N);
}


/* This function is called once per RR interval.  It should return the length
   of the next (simulated) RR interval in seconds.
*/

float generate(void)
{
    float rr;
    
    count++;
    rr= trpeaks_new[count+1]-trpeaks_new[count];

    /* trap rounding errors */
    if(rr<(2.0/SPS)) rr = 2.0/SPS;

    return(rr);
}
    
int main(int argc, char **argv)
{
    float t = 0.0;	    /* sum of intervals since the beginning of the
			       simulation, in seconds */
    long ts = 0;	    /* t, in sample intervals */
    long tsp = 0;	    /* previous value of ts */
    long tmax = 24*60*60;   /* 24 hours, in seconds */
    long seed;		    /* a 32-bit random variable that can be used to
			       initialize the generator */
    float Pe   = 0.0003;    /* Probability of ectopy ~ 1 per hr */
    float Pn   = 0.0048;    /* Probability of noise ~ 16 per hr */
    long atol();

    if (argc < 2) {
	fprintf(stderr, "usage: %s seed [tmax] P(ectopy) P(noise)\n", argv[0]);
	exit(1);
    }
    seed = atol(argv[1]);

    if (argc > 2)
	tmax = atol(argv[2]);
    else 
      tmax=115;

    if(tmax<115) 
        tmax=115;

    if(argc>3) Pe = atof(argv[3]);

    if(argc>4) Pn = atof(argv[4]);

    /* Initialise and runn RR interval generator, pass in length and
       probability of artefacts and nosie */
    initialize(seed, tmax, Pe, Pn);

    /* generate RR interval */
    while ((t += generate()) < tmax) {	/* add RR interval to running time */
	/* calculate and output a quantized RR interval */
	ts = (long)(SPS*t + 0.5);
	printf("%5.3f\n", (ts - tsp)/((float)SPS));
	tsp = ts;
    }

    exit(0);
}
