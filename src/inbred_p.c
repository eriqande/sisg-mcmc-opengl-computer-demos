/* 
	Just a little blurb about the general sequence in which the chain goes about things:
	
	1. Propose a change to p  (draw a colored triangle on p-relevant screens.  Use PROPOSE_COLOR) 
	2. Accept or reject the change to p (make the triangle ACCEPT_COLOR for accept and REJECT_COLOR for reject.)
	3. Propose a change to f (draw a PROPOSE_COLOR triangle on f-relevant screens)
	4. Accept or reject the change to f (draw the accept/reject decision as above)
	5. cycling over the z's	
		-Propose a change to z_i
		-Accept or reject the change to z_i
	6. accumulate in the average.
	
I will use a system of bitmasking to control all this because sometimes I will want to toggle certain update types.  

*/

#define UN_EXTERN



#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "ranlib.h"

#include "MathStatRand.h"
#include "MCTypesEtc.h"
#include "ECA_Opt3.h"

#include "GFMCMC.h"
#include "GFMCMC_DrawingUtilities.h"
#include "ECA_MemAlloc.h"

/* bitmasks to toggle certain operations */
#define PROP_F  0x1
#define PROP_P  0x2
#define PROP_Z  0x4
#define ACC_REJ_F 0x8
#define ACC_REJ_P 0x10
#define ACC_REJ_Z 0x12
#define UPDATE_AVES 0x14



#define MAX_TRACK_LEN 1000000


/* an enum to tell us what type of simulation to do */
#define ip_NUM_SIM_TYPES 4
typedef enum {
	P_AND_F_JOINTLY,
	P_AND_F_COMPONENT_WISE,
	P_F_Z_SILLY,
	P_F_Z_GIBBS
} TypeOfSim;

/* for roughly keeping track of where we are in a particular update */
typedef enum {
	PROPOSED,
	REJECTED,
	ACCEPTED,
	STAY_THERE,
	MOVE_TO_NEW,
	ACCUMULATING,
	GIBBSING
} StepStatus;

/* a data structure to hold everything */
typedef struct {

	/* the data */
	int N;  /* total number of individuals observed */
	int n[3]; /* number of inds have 0, 1, or 2 copies of the a allele. So prob of an [0] is p^2 */
	int *y; /* vector of genotypic states of the N individuals */


	/* the parameters and latent variables */
	dval *p; /* frequency of allele A */
	dval *f; /* inbreeding parameter */
	ival **z; /* array of indicators for whether indivs are inbred or not.  1 = inbred.  0 = not. */
	
	/* parameters of the prior distributions for p and f */
	double p_pri[2];  /* parameters of a beta dsn */
	double f_pri[2]; 
	
	/* standard deviations of the proposal dsns for p and f */
	double p_sd;
	double f_sd;
	
	/* the joint data probability (prior times full-data likelihood)  */
	dval *joint;
	
	/* the values of the proposals */
	double p_p;
	double f_p;
	int *z_p;
	
	/* whether or not the proposals get accepted.  0=reject and 1=accept the proposal */
	ival *p_a;
	ival *f_a;
	int *z_a;
	
	/* to count the number of 0's and 1's amongst the z's */
	int nz[2];
	
	/* count the number of the two different allelic types */
	int az[2]; 	/* az[0] is number of A alleles az[0] is number of A alleles. i think
				 looking at my code it appears that az[0] corresponds to p and az[1] to 1-p */
	
	
	/* arrays to store values for plotting scatter plots */
	float *p_track;
	float *f_track;
	int f_p_scatter_idx;
	
	/* a trace struct for the joint prob trace */
	sliding_trace_struct *joint_trace;
	
	/* trace structs for p and f */
	sliding_trace_struct *p_trace;
	sliding_trace_struct *f_trace;
	
	
	/* some variables that will control whether or not certain variables will get updated in a sweep */
	int Update_p;
	int Update_f;
	int Update_z;
	
	/* a variable that tells us what kind of simulation we are doing */
	TypeOfSim SimType;
	
	/* a variable for keeping track of which step we are on in any given sweep */
	int Step;
	
	/* a string that says what the current step is */
	char StepString[1000];
	
	StepStatus StStat;
	
	unsigned int OpsMask;  /* bit mask telling us which operations to do */
	unsigned int DoneMask; /* bit mask telling us what we have done so far in this sweep */
	
	/* this is for cycling over the z proposal/accept cycles one at a time */
	int z_p_current; 
	int z_a_current;
	
	int NumAccum;
	
	long Seed1;
	long Seed2;
	
	int SampleOnlyF;  /* freezes all values in the model except f */
	int SampleOnlyP;
	
	int ShowNormalF;
	int ShowBetaOnF;
	
	int ShowNormalP;
	int ShowBetaOnP;
	
	
} inbred_vars_struct;


/* global variables */
inbred_vars_struct *gIBV;

void ResetAverages(inbred_vars_struct *V)
{
	InitDvalSummaryToZero(V->p);
	InitDvalSummaryToZero(V->f);
	
	InitIvalSummaryToZero(V->p_a);
	InitIvalSummaryToZero(V->f_a);
	
	gfduInitSlidingTraceToZero(V->p_trace);
	gfduInitSlidingTraceToZero(V->f_trace);
	
	V->f_p_scatter_idx=0;
	
	V->NumAccum=0;
	
}

void InitChain(inbred_vars_struct *V) 
{

	ResetAverages(V);
	
	/* this needs to be much more filled out */
	V->p->v = ranf();
	V->f->v = ranf();

}



/* returns the joint distribution of all parameters.
	The value it returns depends on the type of sim.  If there
	are no latent variables in the sim, they get ignored.
	If Proposed==1 it computes the JointDensity for the proposed
	values.  Otherwise not.  Note that I drop some inconsequential
	normalizing constant terms from the result */
double LogJointDensity(inbred_vars_struct *V, int Proposed)
{
	int i,*z;
	double p,f,q;
	double sum;
	
	if(Proposed==0) {
		p = V->p->v;
		f = V->f->v;
		z = V->z_p;  /* note, we only call this function for f and p updates, and we treat the z's as given */
	}
	else if(Proposed==1) {
		p = V->p_p;
		f = V->f_p;
		z = V->z_p;
	}	
	else {
		fprintf(stderr,"Error! Argument \"Proposed\" to JointDensity must be 0 or 1, not %d.  Exiting.\n",Proposed);
		exit(1);
	}
	
	/* here is our shorthand for 1-p */
	q=1-p;

	/* start with the priors for p and f */
	sum = LogBetaPDF(V->p_pri[0], V->p_pri[1], p) + LogBetaPDF(V->f_pri[0], V->f_pri[1], f);
	 
	/* then, depending on the type of sim we are doing, multiply in the other factors */
	if(V->SimType==P_AND_F_JOINTLY || V->SimType==P_AND_F_COMPONENT_WISE) {
		sum += V->n[0]*log( (1-f)*p*p + f*p)   +    V->n[1]*log( (1-f)*(2*p*q))   +    V->n[2]*log( (1-f)*q*q + f*q);
	}
	else if(V->SimType==P_F_Z_SILLY) {
		/* here is our shorthand for 1-p */
		q=1-p;
		
		/* start with the priors for p and f */
		sum = LogBetaPDF(V->p_pri[0], V->p_pri[1], p) + LogBetaPDF(V->f_pri[0], V->f_pri[1], f);

		/* then cycle over the sampled individuals and keep putting them on the product (sum of logs).
		This will be mightily inefficient */
		for(i=0;i<V->N;i++)  { double temp;
			switch(V->y[i]) {
				case(0):
					temp = p;
					if(V->z[i]->v==0) {temp *= p*(1-f);}
					else {temp *= f;}
					break;
				case(1):
					temp = (1-f)*2.0*p*q;
					break;
				case(2):
					temp = q;
					if(V->z[i]->v==0) {temp *= q*(1-f);}
					else {temp *= f;}
					break;
				default:
					fprintf(stderr,"Error! unknown genotype in StepAheadGibbs(). Exiting\n");
					exit(1);
			}
			sum += log(temp);
		}
	}
	return(sum);
}





/*
	With the joint updates, here is what happens for each value of V->Step
	0. Propose new values of p and f
	1. Accept/reject proposed values
	2. Sit there for a step, record values in averages
	
	I increment the values of V->Step in here, so, we have to use the 
	incremented versions to figure out what to draw on the screens.
*/
void StepAheadJoint_P_F_Jointly(inbred_vars_struct *V)
{
	if(V->Step==0) {
		V->p_p = gennor((float)V->p->v, (float)V->p_sd);
		V->f_p = gennor((float)V->f->v, (float)V->f_sd);
		V->Step++;
		sprintf(V->StepString,"Propose New Value for (p,f)  p*=%.4f  f*=%.4f",V->p_p,V->f_p);
		V->StStat=PROPOSED;
		return;
	}
	else if(V->Step==1) {
		if(V->p_p<=0 || V->p_p>=1 || V->f_p<=0 || V->f_p>=1) {
			V->p_a->v=0;
			V->f_a->v=0;
			V->StStat=REJECTED;
			sprintf(V->StepString,"REJECT New Value for (p,f)");
			V->Step++;
			return;
		}
		else {double R;
			R = LogJointDensity(V,1) - LogJointDensity(V,0);  /* the proposals are symmetrical, and so drop out */
			if(log(ranf()) < R) {
				V->p_a->v=1;
				V->f_a->v=1;
				V->StStat=ACCEPTED;
				sprintf(V->StepString,"ACCEPT New Value for (p,f)");
			}
			else {
				V->p_a->v=0;
				V->f_a->v=0;
				V->StStat=REJECTED;
				sprintf(V->StepString,"REJECT New Value for (p,f)");
			}
			V->Step++;
			return;
		}
	}
	else if(V->Step==2) {
		if(V->StStat==REJECTED) {
			sprintf(V->StepString,"Stay Where you Are!");
			V->StStat=STAY_THERE;
		}
		else if(V->StStat==ACCEPTED) {
			sprintf(V->StepString,"Move to the new value of (p,f)!");
			V->StStat=MOVE_TO_NEW;
		}
		if(V->p_a->v==1) {V->p->v = V->p_p;  }
		if(V->f_a->v==1) {V->f->v = V->f_p;  }
		V->Step++;
	}
	else if(V->Step==3) {
		V->StStat=ACCUMULATING;
		sprintf(V->StepString,"Adding Current Values to Running Averages");
		V->Step=0;
		return;
	}
}	



/* this one updates P and F separately.  Averages are still accumulated for everything at the very end */
void StepAheadJoint_P_F_Separately(inbred_vars_struct *V)
{
	if(V->Step==0) {  /* propose update of p */
		if(V->SampleOnlyF==1) {
			V->Step+=3;
			sprintf(V->StepString,"Skip p update.  Only sampling f",V->p_p);
			return;
		}
		V->p_p = gennor((float)V->p->v, (float)V->p_sd);
		V->f_p = V->f->v;
		V->Step++;
		sprintf(V->StepString,"Propose New Value for p  p*=%.4f",V->p_p);
		V->StStat=PROPOSED;
		return;
	}
	else if(V->Step==1) {
		if(V->p_p<=0 || V->p_p>=1) {
			V->p_a->v=0;
			V->StStat=REJECTED;
			sprintf(V->StepString,"REJECT New Value for p");
			V->Step++;
			return;
		}
		else {double R;
			R = LogJointDensity(V,1) - LogJointDensity(V,0);  /* the proposals are symmetrical, and so drop out */
			if(log(ranf()) < R) {
				V->p_a->v=1;
				V->StStat=ACCEPTED;
				sprintf(V->StepString,"ACCEPT New Value for p");
			}
			else {
				V->p_a->v=0;
				V->StStat=REJECTED;
				sprintf(V->StepString,"REJECT New Value for p");
			}
			V->Step++;
			return;
		}
	}
	else if(V->Step==2) {
		if(V->StStat==REJECTED) {
			sprintf(V->StepString,"Stay Where you Are (Don't Change p)!");
			V->StStat=STAY_THERE;
		}
		else if(V->StStat==ACCEPTED) {
			sprintf(V->StepString,"Move to the new value of p!");
			V->StStat=MOVE_TO_NEW;
		}
		if(V->p_a->v==1) {V->p->v = V->p_p;  }
		V->Step++;
		return;
	}
	else if(V->Step==3) {  /* propose a new value for f */
		if(V->SampleOnlyP==1) {
			V->Step+=3;
			sprintf(V->StepString,"Skip f update.  Only sampling p",V->p_p);
			return;
		}
		V->f_p = gennor((float)V->f->v, (float)V->f_sd);
		V->p_p = V->p->v;
		V->Step++;
		sprintf(V->StepString,"Propose New Value for f  f*=%.4f",V->f_p);
		V->StStat=PROPOSED;
		return;
	}
	else if(V->Step==4) {
		if(V->f_p<=0 || V->f_p>=1) {
			V->f_a->v=0;
			V->StStat=REJECTED;
			sprintf(V->StepString,"REJECT New Value for f");
			V->Step++;
			return;
		}
		else {double R;
			R = LogJointDensity(V,1) - LogJointDensity(V,0);  /* the proposals are symmetrical, and so drop out */
			if(log(ranf()) < R) {
				V->f_a->v=1;
				V->StStat=ACCEPTED;
				sprintf(V->StepString,"ACCEPT New Value for f");
			}
			else {
				V->f_a->v=0;
				V->StStat=REJECTED;
				sprintf(V->StepString,"REJECT New Value for f");
			}
			V->Step++;
			return;
		}
	}
	else if(V->Step==5) {
		if(V->StStat==REJECTED) {
			sprintf(V->StepString,"Stay Where you Are (Don't Change f)!");
			V->StStat=STAY_THERE;
		}
		else if(V->StStat==ACCEPTED) {
			sprintf(V->StepString,"Move to the new value of f!");
			V->StStat=MOVE_TO_NEW;
		}
		if(V->f_a->v==1) {V->f->v = V->f_p;  }
		V->Step++;
		return;
	}
	else if(V->Step==6) {
		V->StStat=ACCUMULATING;
		sprintf(V->StepString,"Adding Current Values to Running Averages");
		V->Step=0;
		return;
	}
}	


/* given genotype g (number of a alleles) and frequency p
of A alleles, and inbreeding parameter f, this returns the
full conditional probability that an individual is inbred
*/
double FullCondProbInbred(int g, double p, double f) 
{
	double q=1.0-p;
	
	if(g==1) {
		return(0.0);
	}
	else if(g==2) {
		return( q*f / (  q*f + (1-f)*q*q ) );
	}
	else if(g==0) {
		return( p*f / (  p*f + (1-f)*p*p ) );
	}
}


/* this one samples the z's via a gibbs sampler (which, in this case also turns out to 
be what you would get if you proposed 0 or 1 with probability 1/2 and then accepted it
based on the hastings ratio.  */ 
void StepAheadJoint_P_F_Silly(inbred_vars_struct *V)
{
	int i;
	
	/* Step 0 is update the z variables */
	if(V->Step==0) {
		if(V->SampleOnlyF==1 || V->SampleOnlyP==1) {
			V->Step++;
			sprintf(V->StepString,"Skip z update.  Only sampling f or p",V->p_p);
			return;
		}
		V->nz[0]=V->nz[1]=0; /* sum these up to use later */
		for(i=0;i<V->N;i++) {
			if(V->y[i]==1) {
				V->z[i]->v = 0;
			}
			else if( ranf() < FullCondProbInbred(V->y[i], V->p->v, V->f->v) ) {
				V->z[i]->v = 1;
			}
			else {
				V->z[i]->v = 0;
			}
			V->nz[ V->z[i]->v ]++;
		}
		V->Step++;
		sprintf(V->StepString,"Simulate New Values for the z's");
		V->StStat=GIBBSING;
		
		/* now, count the number of alleles "sampled" while we are at it */
		V->az[0]=V->az[1]=0;
		for(i=0;i<V->N;i++)  {
			switch(V->y[i]) {
				case(0):
					if(V->z[i]->v==0) {V->az[0]+=2;}
					else {V->az[0]+=1;}
					break;
				case(1):
					V->az[0]+=1;
					V->az[1]+=1;
					break;
				case(2):
					if(V->z[i]->v==0) {V->az[1]+=2;}
					else {V->az[1]+=1;}
					break;
				default:
					fprintf(stderr,"Error! unknown genotype in StepAheadGibbs(). Exiting\n");
					exit(1);
			}
		}

		
		
		return;
	}

	if(V->Step==1) {  /* propose update of p */
		if(V->SampleOnlyF==1) {
			V->Step+=3;
			sprintf(V->StepString,"Skip p update.  Only sampling f",V->p_p);
			return;
		}
		V->p_p = gennor((float)V->p->v, (float)V->p_sd);
		V->f_p = V->f->v;
		V->Step++;
		sprintf(V->StepString,"Propose New Value for p  p*=%.4f",V->p_p);
		V->StStat=PROPOSED;
		return;
	}
	else if(V->Step==2) {
		if(V->p_p<=0 || V->p_p>=1) {
			V->p_a->v=0;
			V->StStat=REJECTED;
			sprintf(V->StepString,"REJECT New Value for p");
			V->Step++;
			return;
		}
		else {double R;
			R = LogJointDensity(V,1) - LogJointDensity(V,0);  /* the proposals are symmetrical, and so drop out */
			if(log(ranf()) < R) {
				V->p_a->v=1;
				V->StStat=ACCEPTED;
				sprintf(V->StepString,"ACCEPT New Value for p");
			}
			else {
				V->p_a->v=0;
				V->StStat=REJECTED;
				sprintf(V->StepString,"REJECT New Value for p");
			}
			V->Step++;
			return;
		}
	}
	else if(V->Step==3) {
		if(V->StStat==REJECTED) {
			sprintf(V->StepString,"Stay Where you Are (Don't Change p)!");
			V->StStat=STAY_THERE;
		}
		else if(V->StStat==ACCEPTED) {
			sprintf(V->StepString,"Move to the new value of p!");
			V->StStat=MOVE_TO_NEW;
		}
		if(V->p_a->v==1) {V->p->v = V->p_p;  }
		V->Step++;
		return;
	}
	else if(V->Step==4) {  /* propose a new value for f */
		if(V->SampleOnlyP==1) {
			V->Step+=3;
			sprintf(V->StepString,"Skip f update.  Only sampling p",V->p_p);
			return;
		}
		V->f_p = gennor((float)V->f->v, (float)V->f_sd);
		V->p_p = V->p->v;
		V->Step++;
		sprintf(V->StepString,"Propose New Value for f  f*=%.4f",V->f_p);
		V->StStat=PROPOSED;
		return;
	}
	else if(V->Step==5) {
		if(V->f_p<=0 || V->f_p>=1) {
			V->f_a->v=0;
			V->StStat=REJECTED;
			sprintf(V->StepString,"REJECT New Value for f");
			V->Step++;
			return;
		}
		else {double R;
			R = LogJointDensity(V,1) - LogJointDensity(V,0);  /* the proposals are symmetrical, and so drop out */
			if(log(ranf()) < R) {
				V->f_a->v=1;
				V->StStat=ACCEPTED;
				sprintf(V->StepString,"ACCEPT New Value for f");
			}
			else {
				V->f_a->v=0;
				V->StStat=REJECTED;
				sprintf(V->StepString,"REJECT New Value for f");
			}
			V->Step++;
			return;
		}
	}
	else if(V->Step==6) {
		if(V->StStat==REJECTED) {
			sprintf(V->StepString,"Stay Where you Are (Don't Change f)!");
			V->StStat=STAY_THERE;
		}
		else if(V->StStat==ACCEPTED) {
			sprintf(V->StepString,"Move to the new value of f!");
			V->StStat=MOVE_TO_NEW;
		}
		if(V->f_a->v==1) {V->f->v = V->f_p;  }
		V->Step++;
		return;
	}
	else if(V->Step==7) {
		V->StStat=ACCUMULATING;
		sprintf(V->StepString,"Adding Current Values to Running Averages");
		V->Step=0;
		return;
	}
}	






void StepAheadGibbs(inbred_vars_struct *V)
{
	int i;
	
	/* Step 0 is update the z variables */
	if(V->Step==0) {
		V->nz[0]=V->nz[1]=0; /* sum these up to use later */
		for(i=0;i<V->N;i++) {
			if(V->y[i]==1) {
				V->z[i]->v = 0;
			}
			else if( ranf() < FullCondProbInbred(V->y[i], V->p->v, V->f->v) ) {
				V->z[i]->v = 1;
			}
			else {
				V->z[i]->v = 0;
			}
			V->nz[ V->z[i]->v ]++;
		}
		V->Step++;
		sprintf(V->StepString,"Simulate New Values for the z's");
		V->StStat=GIBBSING;
		return;
	}
	
	/* step 1 is update p */
	if(V->Step==1) {
	
		/* first, count the number of alleles "sampled" */
		V->az[0]=V->az[1]=0;
		for(i=0;i<V->N;i++)  {
			switch(V->y[i]) {
				case(0):
					if(V->z[i]->v==0) {V->az[0]+=2;}
					else {V->az[0]+=1;}
					break;
				case(1):
					V->az[0]+=1;
					V->az[1]+=1;
					break;
				case(2):
					if(V->z[i]->v==0) {V->az[1]+=2;}
					else {V->az[1]+=1;}
					break;
				default:
					fprintf(stderr,"Error! unknown genotype in StepAheadGibbs(). Exiting\n");
					exit(1);
			}
		}
		
		/* then, sample a new value of p */
		V->p->v = genbet(V->p_pri[0]+V->az[0], V->p_pri[1]+V->az[1]);
		sprintf(V->StepString,"Simulate a New Value for p");
		V->StStat=GIBBSING;
		V->p_a->v = 1; /* this is always the case with Gibbs sampling */
		V->Step++;
	}
	
	else if(V->Step==2)  { /* step 2 is update f */
		V->f->v = genbet(V->f_pri[1]+V->nz[1],V->f_pri[0]+V->nz[0]); 
		sprintf(V->StepString,"Simulate a New Value for f");
		V->StStat=GIBBSING;
		V->f_a->v = 1; /* always the case for Gibbs sampling */
		V->Step++;
		return;
	}
	else if(V->Step==3) {
		sprintf(V->StepString,"Rest and accumulate averages");
		V->StStat=ACCUMULATING;
		V->Step=0;
		return;
	}
	
}



void AccumulateAverages(inbred_vars_struct *V)
{
	/* under all circumstances, we accumulate the average for p and f */
	IncrementDval(V->p);
	IncrementDval(V->f);
	IncrementIval(V->p_a);
	IncrementIval(V->f_a);
	V->p_track[V->f_p_scatter_idx]=V->p->v;
	V->f_track[V->f_p_scatter_idx++]=V->f->v;
	
	gfduIncrementSlidingTraceNoGlobalSweepCounter(V->p_trace, V->p->v,V->p->NumAved);
	gfduIncrementSlidingTraceNoGlobalSweepCounter(V->f_trace, V->f->v,V->p->NumAved);
}


/*  allocate memory to as many fields as you can before knowing what the data look like
	and set some default values  */
inbred_vars_struct *AllocInbredVars(void)
{
	int i;
	
	inbred_vars_struct *ret = (inbred_vars_struct *)malloc(sizeof(inbred_vars_struct));
	
	ret->p = AllocDval(0,1,.01);
	ret->f = AllocDval(0,1,.01);
	ret->joint = AllocDval(0,1,.01);
	
	for(i=0;i<2;i++) {  /* here the default is to set them to uniform */
		ret->p_pri[i]=1;
		ret->f_pri[i]=1;
		ret->az[i]=0;
		ret->nz[i]=0;
	}
	
	ret->p_p=0;
	ret->f_p=0;
	
	ret->p_a = AllocIval(0,1,-1);
	ret->f_a = AllocIval(0,1,-1);

	ret->Update_p = 1;
	ret->Update_f = 1;
	ret->Update_z = 1;
	
	
	ret->p_track = (float *)calloc(MAX_TRACK_LEN,sizeof(float));
	ret->f_track = (float *)calloc(MAX_TRACK_LEN,sizeof(float));

	ret->NumAccum = 0;
	
	ret->p_trace = gfduAllocSlidingTrace(500,1);
	ret->f_trace = gfduAllocSlidingTrace(500,1);
	
	ret->SampleOnlyF = 0;
	ret->SampleOnlyP = 0;
	
	ret->ShowNormalF = 0;
	ret->ShowBetaOnF = 0;

	ret->ShowNormalP = 0;
	ret->ShowBetaOnP = 0;

	return(ret);
}


inbred_vars_struct *getInbredP_Options(int argc, char *argv[])
{
	inbred_vars_struct *ret;
	int TraceLength;
	int ns_f = 0,
		simtype_f=0;
	
	DECLARE_ECA_OPT_VARS  
 
	/* allocate memory and set some defaults */
	ret = AllocInbredVars(); 
	TraceLength = 2000000;  /* some day might give an option to change this */
	ret->p_sd = .05;
	ret->f_sd = .05;
	ret->p_pri[0] = ret->p_pri[1] = 1.0;
	ret->f_pri[0] = ret->f_pri[1] = 1.0;
	
	
	/* some information about the program, author, etc */
	SET_PROGRAM_NAME("inbred_p");  /* Note:  QUOTED string */
	SET_PROGRAM_SHORT_DESCRIPTION("illustrate MCMC by estimating p and f from a simple inbreeding model"); /* QUOTED string */
	SET_PROGRAM_LONG_DESCRIPTION(This program implements a variety of different MCMC schemes to estimate the parameters
								p and f in a simple model of inbreeding.  p is the frequency of the \042A\042 allele at 
								a diallelic locus having alleles \042A\042 and \042a\042.  f is the inbreeding parameter---i.e.\054 
								the probability that the two gene copies at the locus are identical by descent.  The user supplies the
								data which is the number of aa individuals\054 the number of heterozygotes (aA or Aa)\054 and the number 
								of AA individuals.  This program was written for the MCMC in Statistical Genetics module of the 
								Summer Insitute in Statistical Genetics at the University of Washington.);  /* UN-QUOTED string! */
	SET_PROGRAM_AUTHOR_STRING("Eric C. Anderson"); /* QUOTED string */
	SET_VERSION("Version XX");  /* QUOTED string */
	SET_VERSION_HISTORY("Version XX written Nov. 8, 2007\nGrew out of simpler programs.\nblah, blah..."); /*QUOTED string */
	 
	 
	BEGIN_OPT_LOOP  
	 
	 
	/* use this to start a block of related options.  Close the block with the CLOSE_SUBSET macro (below)  */
	OPEN_SUBSET(Data Input Options,  /* Name of Subset for output in --help, --help-full, and --help-nroff */
	Data Input, /* Name of subset for lumping them together in GuiLiner */
	Data input /* Description of subset.  Not really used currently */
	)   /* NOTE.  These are all UNQUOTED strings.  Beware of commas. */
	
	if(REQUIRED_OPTION(
		Genotype Data,
		ns_f,
		n,
		geno-counts,
		J0 J1 J2,
		Number of indivs with 0\054 1\054 or 2 copies of the \042a\042 allele,
		Number of individuals with 0\054 1\054 or 2 copies of the \042a\042 allele.  Note that because
		  p is the frequency of the \042A\042 allle\054 the probability (w/o inbreeding) that an 
		  individual as 0 copies of the \042a\042 allele is p^2. )) {
		if(ARGS_EQ(3)) {int i; int j; int k;
			ret->n[0] = GET_INT;
			ret->n[1] = GET_INT;
			ret->n[2] = GET_INT;
			ret->N = ret->n[0] + ret->n[1] + ret->n[2];
			
			/* allocate memory to the z's, etc */
			ret->z = IvalVector(0, ret->N-1, 0, 1, 2);
			ret->z_p = (int *)calloc(ret->N, sizeof(int));
			ret->z_a = (int *)calloc(ret->N, sizeof(int));
			ret->y = (int *)calloc(ret->N, sizeof(int));
			
			/* then fill in the y's */
			for(k=0,i=0;i<3;i++) {
				for(j=0;j<ret->n[i];j++)  {
					ret->y[k++] = i;
				}
			}
		}
	}
	
	if(OPTION(
		Type Of Simulation,
		simtype_f,
		t,
		sim-type,
		S,
		set the type of MCMC simulation to start with,
		Set the type of MCMC simulation to start with.  Available choices are: \042joint\042\054 \042component\042\054
		 \042latent-mh\042\054 or \042latent-gibbs\042.)) {
			if(ARGS_EQ(1)) { char tempstr[5000];
				GET_STR(tempstr);
				if(strcmp(tempstr,"joint")==0) {
					ret->SimType = P_AND_F_JOINTLY;
				}
				else if(strcmp(tempstr,"component")==0) {
					ret->SimType = P_AND_F_COMPONENT_WISE;
				}
				else if(strcmp(tempstr,"latent-mh")==0) {
					ret->SimType = P_F_Z_SILLY;
				}
				else if(strcmp(tempstr,"latent-gibbs")==0) {
					ret->SimType = P_F_Z_GIBBS;
				}
				else {
					fprintf(stderr, "Error! Unrrecognized string argument \"%s\" to option -t/--sim-type. Exiting.\n",tempstr);
					exit(1);
				}
				
			}
		
	}
	 
	CLOSE_SUBSET;  /* done with the greeting-related options */
	 
	END_OPT_LOOP   /* Macro that loses the Main ECA_Opt Loop and does some error checking 
				  behind the scenes */
				  
	/* down here we set some more values and do some more allocation of memory */
	ret->joint_trace = gfduAllocSlidingTrace(TraceLength,1);
	
	
	return(ret);
		
}




/* handle the keyboard input to the console window, allowing 't'
   to toggle between Jeffreys and Uniform priors on Theta and
   'p' to do the same for Pi.  This overrides whatever 'p' and 't'
   might be in the Normal Key Processing Function.  */
void ipConsoleKeys(unsigned char key, int x, int y)
{
	switch(key) {
		case('g'):
			gIBV->SimType=P_F_Z_GIBBS;
			gIBV->Step=0;
			return;
		case('c'):
			gIBV->SimType=P_AND_F_COMPONENT_WISE;
			gIBV->Step=0;
			return;
		case('j'):
			gIBV->SimType=P_AND_F_JOINTLY;
			gIBV->Step=0;
			return;
		case('S'):
			gIBV->SimType=P_F_Z_SILLY;
			gIBV->Step=0;
			return;
		case('f'):
			gIBV->SampleOnlyF = (gIBV->SampleOnlyF + 1) % 2;
			return;
		case('p'):
			gIBV->SampleOnlyP = (gIBV->SampleOnlyP + 1) % 2;
			return;
		
		default:
			gfmProcessNormalKeys(key, x, y);
	}
}


void ipFHistKeyFunc(unsigned char key, int x, int y)
{
	switch(key) {
		case('b'):
			gIBV->ShowBetaOnF = (gIBV->ShowBetaOnF + 1) % 2;
			return;
		case('n'):
			gIBV->ShowNormalF = (gIBV->ShowNormalF + 1) % 2;
			return;
		default:
			gfmProcessNormalKeys(key, x, y);
	}
}

void ipPHistKeyFunc(unsigned char key, int x, int y)
{
	switch(key) {
		case('b'):
			gIBV->ShowBetaOnP = (gIBV->ShowBetaOnP + 1) % 2;
			return;
		case('n'):
			gIBV->ShowNormalP = (gIBV->ShowNormalP + 1) % 2;
			return;
		default:
			gfmProcessNormalKeys(key, x, y);
	}
}



void ipDrawConsoleWindow(void)
{
	char temp[250];
	
	glColor3fv(gWindowsSettings[glutGetWindow()]->ColorScheme->Text);
	gfmRenderBitmapString(-.1f,1.0f,GLUT_BITMAP_HELVETICA_12, "Program: inbred_p");
	gfmRenderBitmapString(-.1f,.9f,GLUT_BITMAP_HELVETICA_12, "Author: E.C. Anderson");
	
	sprintf(temp,"Data: \"%s\"", "boing");
	gfmRenderBitmapString(-.1f,.8f,GLUT_BITMAP_HELVETICA_12, temp);
	
	switch(gIBV->SimType) {
		case(P_AND_F_JOINTLY):
			sprintf(temp,"Type of Sim: P and F Jointly");
			break;
		case(P_AND_F_COMPONENT_WISE):
			sprintf(temp,"Type of Sim: P and F Component-wise");
			break;
		case(P_F_Z_SILLY):
			sprintf(temp,"Type of Sim: P, F, Z Silly");
			break;
		case(P_F_Z_GIBBS):
			sprintf(temp,"Type of Sim: P, F, Z Gibbs");
			break;
		default:
			sprintf(temp,"Type of Sim: UNKNOWN!");
			break;
	}
	gfmRenderBitmapString(-.1f,.7f,GLUT_BITMAP_HELVETICA_12, temp);

	if(gIBV->SampleOnlyF==1) {
		sprintf(temp,"Fix all but F (f)");
		gfmRenderBitmapString(-.1f,.6f,GLUT_BITMAP_HELVETICA_12, temp);
	}
	
	if(gIBV->SampleOnlyP==1) {
		sprintf(temp,"Fix all but P (p)");
		gfmRenderBitmapString(-.1f,.5f,GLUT_BITMAP_HELVETICA_12, temp);
	}
	
	sprintf(temp,"Options: (j)oint, (c)omponent, (g)ibbs, (S)illy, fix(f), fix(p)");
		gfmRenderBitmapString(-.1f,.4f,GLUT_BITMAP_HELVETICA_12, temp);
	
	sprintf(temp,"Num of Iterations: %d",gIBV->NumAccum);
	gfmRenderBitmapString(-.1f,.2f,GLUT_BITMAP_HELVETICA_12, temp);
	
	if(gGo) 
		sprintf(temp,"RUNNING!");
	else 
		sprintf(temp,"STOPPED!");
		
	gfmRenderBitmapString(-.1f,.3f,GLUT_BITMAP_HELVETICA_12, temp);
	
}

void ipDrawStatusWindow(void)
{
	glColor3fv(gWindowsSettings[glutGetWindow()]->ColorScheme->Text);
	gfmRenderBitmapString(-.1f,.4f,GLUT_BITMAP_TIMES_ROMAN_24, gIBV->StepString);

	
}

void ipDrawPF_Scatter(void)
{
	int i;
	gsCS;
	
	/* first, draw the axes */
	gfmDrawXYAxes();

		
	glColor3fv(CS->Axes);
	glPointSize(1.0);
	glBegin(GL_POINTS);
		
		for(i=0;i<gIBV->f_p_scatter_idx;i++)  {

			glVertex2f((GLfloat)(gIBV->f_track[i]),(GLfloat)(gIBV->p_track[i])  ) ;
			
		}
	glEnd();
	
	/* then draw the current point */
	glColor3fv(CS->Series[1]);
	glPointSize(7.0);
	glBegin(GL_POINTS);
		glVertex2f(gIBV->f->v,gIBV->p->v);
	glEnd();

	
	/* then, if we just proposed something, draw that in---as a point with a line to it */
	if(gIBV->StStat==PROPOSED) {
		glColor3fv(CS->Series[2]);
		glPointSize(7.0);
		glBegin(GL_POINTS);
			glVertex2f(gIBV->f_p,gIBV->p_p);  /* note! when doing component-wise, we have to set the proposed value of the unchanged variable to
												the current values of those variables! */
		glEnd();
		
		glColor3fv(CS->Series[2]);
		glBegin(GL_LINE_STRIP);
			glVertex2f(gIBV->f->v,gIBV->p->v);
			glVertex2f(gIBV->f_p,gIBV->p_p);
		glEnd();
	}
	
	
	/* now, do it in a different color if it is accepted etc */
	if(gIBV->StStat==ACCEPTED || gIBV->StStat==REJECTED) {
		if(gIBV->StStat==ACCEPTED) {
			glColor3fv(CS->Series[3]);
		}
		else {
			glColor3fv(CS->Series[4]);
		}
		glPointSize(7.0);
		glBegin(GL_POINTS);
			glVertex2f(gIBV->f_p,gIBV->p_p);  /* note! when doing component-wise, we have to set the proposed value of the unchanged variable to
													   the current values of those variables! */
		glEnd();
		
		glBegin(GL_LINE_STRIP);
			glVertex2f(gIBV->f->v,gIBV->p->v);
			glVertex2f(gIBV->f_p,gIBV->p_p);
		glEnd();
	}
	
		
	
}



void ipDrawPHistogram(void)
{
	double d;
	gsCS;
	
	glColor3fv(CS->Series[4]);
	gfduDrawDoubleHistAsRec(gIBV->p->Hist, gIBV->p->Hist->dTot, gIBV->p->v);
	
	
	/* draw the normal proposal dsn */
	if(gIBV->ShowNormalP==1) {
		glColor3fv(CS->Series[6]);
		glLineWidth(3.0);
		glBegin(GL_LINE_STRIP);
		for(d=.001;d<1.000;d+=.005) {
				glVertex2f(d,NormalPDF(gIBV->p->v, gIBV->p_sd*gIBV->p_sd,d));
		}
		glEnd();
	}
	
	/* draw the full conditional for p */
	if( (gIBV->SimType==P_F_Z_GIBBS || gIBV->SimType==P_F_Z_SILLY)  && gIBV->ShowBetaOnP==1) {
		glColor3fv(CS->Series[7]);
		glLineWidth(3.0);
		glBegin(GL_LINE_STRIP);
		for(d=.001;d<1.000;d+=.005) {
				glVertex2f(d,exp(LogBetaPDF(gIBV->p_pri[0]+gIBV->az[0], gIBV->p_pri[1]+gIBV->az[1],d)));
		}
		glEnd();
	}
	
	glLineWidth(1.0); /* go back to the previous line width */


	
	gfmDrawXYAxes();	
}


void ipDrawPTrace(void)
{	
	gsCS;
	
	glColor3fv(CS->Series[4]);

	gfduDrawSlidingTrace(gIBV->p_trace);
	
	gfmDrawXYAxes();

}


void ipDrawFHistogram(void)
{
	double d;
	gsCS;
	
	glColor3fv(CS->Series[5]);
	gfduDrawDoubleHistAsRec(gIBV->f->Hist, gIBV->f->Hist->dTot, gIBV->f->v);
	
	
	/* draw the normal proposal dsn */
	if(gIBV->ShowNormalF==1) {
		glColor3fv(CS->Series[6]);
		glLineWidth(3.0);
		glBegin(GL_LINE_STRIP);
		for(d=.001;d<1.000;d+=.005) {
				glVertex2f(d,NormalPDF(gIBV->f->v, gIBV->f_sd*gIBV->f_sd,d));
		}
		glEnd();
	}
	
	/* draw the full conditional for f */
	if( (gIBV->SimType==P_F_Z_GIBBS || gIBV->SimType==P_F_Z_SILLY)  && gIBV->ShowBetaOnF==1) {
		glColor3fv(CS->Series[7]);
		glLineWidth(3.0);
		glBegin(GL_LINE_STRIP);
		for(d=.001;d<1.000;d+=.005) {
				glVertex2f(d,exp(LogBetaPDF(gIBV->f_pri[1]+gIBV->nz[1], gIBV->f_pri[0]+gIBV->nz[0],d)));
		}
		glEnd();
	}
	
	glLineWidth(1.0); /* go back to the previous line width */
	
	gfmDrawXYAxes();
	
	
}


void ipDrawFTrace(void)
{	
	gsCS;
	
	glColor3fv(CS->Series[5]);

	gfduDrawSlidingTrace(gIBV->f_trace);
	
	gfmDrawXYAxes();

}


void ipDrawIndivs(void)
{
	int i,l,c;
	GLfloat xs = .8f, ys = .8f;
	GLfloat x,y, yy;
	ColorScheme3f *CS = gWindowsSettings[glutGetWindow()]->ColorScheme;
	XY_axis_struct *XA = gWindowsSettings[glutGetWindow()]->Xaxis;
	clipping_struct *C = gWindowsSettings[glutGetWindow()]->Clips;
	GLfloat WoH = gfmWoH();
	char temp[100];
	GLfloat texth;
		
	for(i=0;i<gIBV->N;i++) {
		y = 2.0f * (GLfloat)gIBV->N - (2.0f * (GLfloat)i);
		
		/*  draw text telling the individual's number */
		glColor3fv(CS->Text);
		switch(gIBV->y[i]) {
			case(0):
				if(gIBV->z[i]->v == 0)  { sprintf(temp,"%d    AA   Not Inbred",i+1); }
				else { sprintf(temp,"%d    AA   Inbred",i+1); }
				break;
			case(1):
				if(gIBV->z[i]->v == 0)  { sprintf(temp,"%d    Aa   Not Inbred",i+1); }
				else { sprintf(temp,"%d    Aa   Inbred",i+1); }
				break;
			case(2):
				if(gIBV->z[i]->v == 0)  { sprintf(temp,"%d    aa   Not Inbred",i+1); }
				else { sprintf(temp,"%d    aa   Inbred",i+1); }				
				break;
		}
		gfmStrokeString(temp, 1.9f ,-.1f, y-1.8f, 0, 0.0f,WoH);

#ifdef NOTDEFD		
		for(l=0;l<3;l++)  {
			x = (GLfloat)l;
			yy = y;
			for(c=0;c<2;c++)  {
				switch(c) {
					case(0):  /* in this case we are doing the observed data */
						if(gIBV->y[i]==0 || gIBV->y[i]==1) { /* homozygous for A */
							glColor3fv(CS->Series[0]);  /*  set the color for the A allele */
						}
						else {
							glColor3fv(CS->Series[1]);  /*  set the color for the a allele */
						}
						glRectf(x,yy,x+xs,yy-ys);
						
						/*  move the second allele over half of a unit: */
						x += 1.0f;
						
						if(gIBV->y[i]==0) { /* homozygous for A */
							glColor3fv(CS->Series[0]);  /*  set the color for the A allele */
						}
						else {
							glColor3fv(CS->Series[1]);  /*  set the color for the a allele */
						}
						glRectf(x,yy,x+xs,yy-ys);
						break;
				}
			}
		}
#endif
	}
		/*  draw a separating line  */
/*		if(XA->DrawIt==1)  {
			glColor3fv(CS->Axes);
			glBegin(GL_LINES);
				glVertex2f(0.0f,y-2.5f);
				glVertex2f((GLfloat)gC->Dat->L-(1.0f-xs),y-2.5f);
			glEnd();
		}
*/			
	
/*	gsHANDLE_EXTREMA(0.0f,(GLfloat)gC->Dat->L,0.0f, 3.0f * (GLfloat)gC->Dat->M); */
	
	
}

void ipDrawF_accept(void)
{
	gsCS;
	glColor3fv(CS->Series[5]);
	glRectf(.25,0.0,.75,gIBV->f_a->Ave);
	gfmDrawXYAxes();
}


void ipDrawP_accept(void)
{
	gsCS;
	glColor3fv(CS->Series[4]);
	glRectf(.25,0.0,.75,gIBV->p_a->Ave);
	gfmDrawXYAxes();
}


void ipDrawZ_A_Counts(void)
{
	char temp[250];
	gsCS;
	
	glColor3fv(CS->Text);
	
	sprintf(temp,"IBD: %d",gIBV->nz[1]);
	gfmRenderBitmapString(-.1f,.9f,GLUT_BITMAP_TIMES_ROMAN_24, temp);
	
	
	sprintf(temp,"~IBD: %d",gIBV->nz[0]);
	gfmRenderBitmapString(-.1f,.7f,GLUT_BITMAP_TIMES_ROMAN_24, temp);
	
	
	
	sprintf(temp,"A: %d",gIBV->az[1]);
	gfmRenderBitmapString(-.1f,.4f,GLUT_BITMAP_TIMES_ROMAN_24, temp);
	
	
	sprintf(temp,"a: %d",gIBV->az[0]);
	gfmRenderBitmapString(-.1f,.2f,GLUT_BITMAP_TIMES_ROMAN_24, temp);
}

void gfmUserDefd_ResetAllAverages(void)
{
	ResetAverages(gIBV);
}

void gfmUserDefd_InitializeChain(void)
{
	long temp1, temp2;
	
	if(gUseSameSeeds == 1) {
		setall((long)gIBV->Seed1, (long)gIBV->Seed2);
	}
	else {
		getsd(&temp1, &temp2);
		gIBV->Seed1 = temp1;
		gIBV->Seed2 = temp2;
	}
	InitChain(gIBV);
	
	
}


void gfmUserDefd_OneStepByGlobals(void)
{
	
	if(gIBV->SimType==P_AND_F_JOINTLY) {
		StepAheadJoint_P_F_Jointly(gIBV);
	}
	else if(gIBV->SimType==P_AND_F_COMPONENT_WISE) {
		StepAheadJoint_P_F_Separately(gIBV);
	}
	else if(gIBV->SimType==P_F_Z_SILLY) {
		StepAheadJoint_P_F_Silly(gIBV);
	}
	else if(gIBV->SimType==P_F_Z_GIBBS) {
		StepAheadGibbs(gIBV);
	}
	
	
	if(gIBV->StStat==ACCUMULATING) {
		AccumulateAverages(gIBV);
		gIBV->NumAccum++;
	}
	
	
	/*
	SingleSweep(gC);
	IncrementValues(gC,1);
	/* update the sliding Trace */
	/*gfduIncrementSlidingTrace(gCompleteDataLogLike,gC->Lat->TheCompleteDataLogLike); */
	

}




/*  this is the function that defines what all the windows are 
#define gfnhCAT_PROBS_SETUP gsDRAW_FUNC(gfnhDrawCategoryProbs); gsNUM_COLOR_KEYS( &(gC->Dat->Gn->v)); gsCOLOR_KEYS(gC->Dat->CategoryNames); gsMIDDLE_CLICK_MENU(CATEGORY_RIGHT_CLICK); gsCOLOR_SCHEME(DEEP_BLUE);
#define gfnhCAT_UNIPRI_PROBS_SETUP gsDRAW_FUNC(gfnhDrawCategoryUniPriProbs); gsNUM_COLOR_KEYS( &(gC->Dat->Gn->v)); gsCOLOR_KEYS(gC->Dat->CategoryNames); gsMIDDLE_CLICK_MENU(CATEGORY_RIGHT_CLICK); gsCOLOR_SCHEME(DEEP_BLUE);

#define gfnhALLELE_FREQ_SETUP gsDRAW_FUNC(gfnhDrawAlleleFreqs); gsNUM_COLOR_KEYS( &(gMaxAlleNum)); gsCOLOR_KEYS(gAlleStrings); gsMIDDLE_CLICK_MENU(ALLELE_RIGHT_CLICK); gsCOLOR_SCHEME(FISHER_PRICE);
*/

void gfmUserDefd_DefineWindows(void)
{


	gsOUTPUT_FILE_PREFIX("InbredP");
	
		
	/*  go ahead and do the definitions: */
	gsCONSOLE("\"Info\" Window");
		gsDRAW_FUNC(ipDrawConsoleWindow);
		gsKEYBOARD_FUNC(ipConsoleKeys);
	
	
	gsNEW_WINDOW(1,"P and F Scatter-Trace");
		gsDRAW_FUNC(ipDrawPF_Scatter);
		gsAXES_POS_UNIT_X;
		
	gsNEW_WINDOW(2,"Current Status");
		gsDRAW_FUNC(ipDrawStatusWindow);
		
	gsNEW_WINDOW(3,"P Histogram");
		gsDRAW_FUNC(ipDrawPHistogram);
		gsKEYBOARD_FUNC(ipPHistKeyFunc);
		gsXLO(-0.1);
		gsXHI(1.1);
		gsYHI(15.0);
		gsYLO(-2.0);
		gsAXES_POS_UNIT_X;
		gsCOLOR_SCHEME(FISHER_PRICE);
		
	gsNEW_WINDOW(4,"F Histogram");
		gsDRAW_FUNC(ipDrawFHistogram);
		gsKEYBOARD_FUNC(ipFHistKeyFunc);
		gsXLO(-0.1);
		gsXHI(1.1);
		gsYHI(15.0);
		gsYLO(-2.0);
		gsAXES_POS_UNIT_X;
		gsCOLOR_SCHEME(FISHER_PRICE);
	
	gsNEW_WINDOW(5,"F Trace");
		gsY_AXIS_CROSSES(0.0);
		gsYLO(-.1)
		gsYHI(1.1)
		gsDRAW_FUNC(ipDrawFTrace);
		gsCOLOR_SCHEME(FISHER_PRICE);
		
		
	gsNEW_WINDOW(6,"P Trace");
		gsY_AXIS_CROSSES(0.0);
		gsYLO(-.1)
		gsYHI(1.1)
		gsDRAW_FUNC(ipDrawPTrace);
		gsCOLOR_SCHEME(FISHER_PRICE);
		
		
	gsNEW_WINDOW(7,"F Prop SD");
		gsDOUBLE_SLIDER("f* SD",gIBV->f_sd,0.0001,.9999);
		
	gsNEW_WINDOW(8,"P Prop SD");
		gsDOUBLE_SLIDER("p* SD",gIBV->p_sd,0.0001,.9999);
		
	gsNEW_WINDOW(9,"Indivs");
		gsDRAW_FUNC(ipDrawIndivs);
		gsCOLOR_SCHEME(FISHER_PRICE);
		gsXLO(-.15f * 1.5);
		gsXHI(1.15f * 1.5);
		gsYLO(-.01f * 2.0f * gIBV->N);
		gsYHI(1.00f * 2.0f * gIBV->N);
		
		gsPADDING(-.17,1.17,-.17,1.17);
		
		
	gsNEW_WINDOW(10,"Ppn Acc F");
		gsDRAW_FUNC(ipDrawF_accept);
		gsAXES_POS_UNIT_X;
		gsNO_XAXIS_NUMS;
		gsXLO(-.5);

	gsNEW_WINDOW(11,"Ppn Acc P");
		gsDRAW_FUNC(ipDrawP_accept);
		gsAXES_POS_UNIT_X;
		gsNO_XAXIS_NUMS;
		gsXLO(-.5);
		
	gsNEW_WINDOW(12,"Z and A counts");
		gsDRAW_FUNC(ipDrawZ_A_Counts);
		
		
	/*
	gsNEW_WINDOW(2,"Average Category Probabilities");
		gfnhCAT_PROBS_SETUP
		gsSET(GFNH_CURRENT_OR_AVES, GFNH_AVES);
		
	
	gsNEW_WINDOW(3,"Current Allele Freqs");
		gfnhALLELE_FREQ_SETUP
		gsSET(GFNH_CURRENT_OR_AVES, GFNH_CURRENT);
		
	gsNEW_WINDOW(4,"Average Allele Freqs");
		gfnhALLELE_FREQ_SETUP
		gsSET(GFNH_CURRENT_OR_AVES, GFNH_AVES);
		
	gsNEW_WINDOW(5,"Allele Frequency Histogram");
		gsDRAW_FUNC(gfnhDrawAlleleFreqsHist);
		gsCOLOR_KEYS(gAlleStrings);
		gsCOLOR_SCHEME(FISHER_PRICE);
		gsAXES_POS_UNIT_X;
*/
		/*  note that NumColorKeys gets set depending on the Displayed Locus from */
		/*  within the function gfnhDrawAlleleFreqsHist() */

/*
	gsNEW_WINDOW(6,"Category Probs Histogram");
		gsDRAW_FUNC(gfnhDrawCatProbsHist);
		gsCOLOR_KEYS(gC->Dat->CategoryNames);
		gsNUM_COLOR_KEYS(&(gC->Dat->Gn->v));
		gsCOLOR_SCHEME(DEEP_BLUE);
		gsAXES_POS_UNIT_X;
		
	gsNEW_WINDOW(7,"Observed Data");
		gsDRAW_FUNC(gfnhDrawObservedData);
		gsCOLOR_SCHEME(FISHER_PRICE);
		gsXLO(-.15f * gC->Dat->L);
		gsXHI(1.15f * gC->Dat->L);
		gsYLO(-.15f * 3.0f * gC->Dat->M);
		gsYHI(1.15f * 3.0f * gC->Dat->M);
		
		gsPADDING(-.17,1.17,-.17,1.17);
		
	
	gsNEW_WINDOW(11,"Current Origins Of Alleles");
		gsDRAW_FUNC(gfnhDrawWs);
		gsCOLOR_SCHEME(FISHER_PRICE);
		gsXLO(-.15f * gC->Dat->L);
		gsXHI(1.15f * gC->Dat->L);
		gsYLO(-.15f * 3.0f * gC->Dat->M);
		gsYHI(1.15f * 3.0f * gC->Dat->M);
		
		gsPADDING(-.17,1.17,-.17,1.17);
		
		
	gsNEW_WINDOW(12,"Average Origins Of Alleles");
		gsDRAW_FUNC(gfnhDrawW_Aves);
		gsCOLOR_SCHEME(FISHER_PRICE);
		gsXLO(-.15f * gC->Dat->L);
		gsXHI(1.15f * gC->Dat->L);
		gsYLO(-.15f * 3.0f * gC->Dat->M);
		gsYHI(1.15f * 3.0f * gC->Dat->M);
		gsPADDING(-.17,1.17,-.17,1.17);
		
		
		
	gsNEW_WINDOW(8,"Complete Data LogL Trace");
		gsDRAW_FUNC(gfnhDrawCompleteDataLogLTrace);
		gsCOLOR_SCHEME(FISHER_PRICE);
		
	gsNEW_WINDOW(10,"Current Kullback Leibler Div By Locus")
		gsDRAW_FUNC(gfnhDrawKullbLeib);
		gsXLO(-.15 * 12.0f);
		gsXHI(1.15 * 12.0f);
		gsYLO(-.15f * gC->Dat->L);
		gsYHI(1.15f * gC->Dat->L);
		gsCOLOR_SCHEME(FISHER_PRICE);
		gsSET(GFNH_CURRENT_OR_AVES, GFNH_CURRENT);
		gsAXES_POS_QUADRANT
		gsNO_Y_AXIS
		
	gsNEW_WINDOW(9,"Average Kullback Leibler Div By Locus")
		gsDRAW_FUNC(gfnhDrawKullbLeib);
		gsXLO(-.15 * 12.0f);
		gsXHI(1.15 * 12.0f);
		gsYLO(-.15f * (GLfloat)gC->Dat->L);
		gsYHI(1.15f * (GLfloat)gC->Dat->L);
		gsCOLOR_SCHEME(FISHER_PRICE);
		gsSET(GFNH_CURRENT_OR_AVES, GFNH_AVES);
		gsAXES_POS_QUADRANT
		gsNO_Y_AXIS
		
	gsNEW_WINDOW(13, "Pi Values")
		gsDRAW_FUNC(gfnhDrawPiValues);
		gsCOLOR_SCHEME(DEEP_BLUE);
		gsCOLOR_KEYS(gC->Dat->CategoryNames);
		gsNUM_COLOR_KEYS(&(gC->Dat->Gn->v));
		gsYLO(-2)
		gsYHI(2)
		gsXLO(-.6)
		gsXHI(1.15)
	
	gsNEW_WINDOW(14, "Pi Histogram")
		gsDRAW_FUNC(gfnhDrawPiHists);
		gsAXES_POS_QUADRANT;
		gsCOLOR_SCHEME(DEEP_BLUE);
		gsCOLOR_KEYS(gC->Dat->CategoryNames);
		gsNUM_COLOR_KEYS(&(gC->Dat->Gn->v));
		
	gsNEW_WINDOW(15,"Current Scaled Likelihood");
		gfnhCAT_UNIPRI_PROBS_SETUP
		gsSET(GFNH_CURRENT_OR_AVES, GFNH_CURRENT);
			
	
	gsNEW_WINDOW(16,"Average Scaled Likelihoods");
		gfnhCAT_UNIPRI_PROBS_SETUP
		gsSET(GFNH_CURRENT_OR_AVES, GFNH_AVES);
		
	gsNEW_WINDOW(20, "Allele Freq Clines");
		gsDRAW_FUNC(gfnhDrawAlleFreqClines);
		gsXLO( gClines->ClineX[0] - .1 * (gClines->ClineX[gClines->NumClineXs-1] - gClines->ClineX[0])  ); 
		gsXHI( gClines->ClineX[gClines->NumClineXs-1] );
		gsYLO(-.20);
		gsYHI(1.20);
		gsCOLOR_SCHEME(FISHER_PRICE);
		gsY_AXIS_CROSSES(0.0f);
		gsX_AXIS_CROSSES(0.0f);
*/		

}


void gfmUserDefd_DefineMenus(void)
{
	/*gUserDefdMenus[ALLELE_RIGHT_CLICK] = glutCreateMenu(gfnhOpenAlleleRightClickWindow);
	glutAddMenuEntry("Open Histogram View For Selected Locus", 5);


	gUserDefdMenus[CATEGORY_RIGHT_CLICK] = glutCreateMenu(gfnhOpenCategoryRightClickWindow);
	glutAddMenuEntry("Open Histogram View For Selected Individual", 0); */
	return;
}





/*  this is a place that a user can do some last minute things like  */
/*  make some new menus and the like.  Here I don't have anything to do. */
void gfmUserDefd_LastWords(void)
{
	return;
}




int main(int argc, char *argv[])
{
	inbred_vars_struct *IBV;
	
	IBV = getInbredP_Options(argc, argv);
	
	InitChain(IBV);
	
	gIBV = IBV;
	
	glutInit(&argc, argv);
	gfmInitGFM();

	
		return(0);
}