// this is a little file to hold all the functions related to the Pritchard et al
// method for dealing with admixed populations.

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>

#include "ECA_MemAlloc.h"
#include "MCTypesEtc.h"
#include "ECA_utilities.h"
#include "MathStatRand.h"
#include "ranlib.h"
#include "GFMCMC.h"
#include "GFMCMC_DrawingUtilities.h"

#include "NewHybrids.h"
#include "MAZ4NewHybs.h"
#include "GLUT_for_PritchEtc.h"


// here is a little function to update the W's at all loci in all individuals
// under the Pritchard model.  Does not store anything in PofW.  HOWEVER, it does
// store the number of genes from each of the different components for
// each individual.  That info is then used for updating Q.
void PritchUpdateW(hyb_chain *C)
{
	int c,s,i,l;
	double normo;  // to accumulate a sum
	double Numer[MAX_NUM_SOURCES_PRITCH]; 
	
	
	
	// cycle over all the individuals
	CYCLE_i(C->Dat)
	
		// cycle over the sources and initialize NumW to 0 to accumulate a sum
		CYCLE_s(C->Lat->PritLat)
			C->Lat->Ind[i]->PritInd->NumW[s]->v = 0;
		END1CYCLE
	
		// cycle over all loci
		CYCLE_l(C->Dat)
		
			// CODOM ROUTINE
			if(C->Dat->LocTypes[l] == CODOM )  {
			
				// cycle over the first and second alleles
				for(c=0;c<2;c++)  {
					normo = 0.0;  // initialize to accumulate a sum
					// then cycle over the possible sources
					CYCLE_s(C->Lat->PritLat)
						if(YnotMissing(C,i,l))  {  // if there is data at the locus, include the allele freq
							// numerator is the product of allele freq and Q
							Numer[s] = C->Lat->Theta[ s ][ l ][ C->Lat->Ind[i]->Y[l][c]->v ]->v  *
										C->Lat->Ind[i]->PritInd->Q[s]->v;
						}
						else {  // otherwise, if the locus is missing, we just have the prior probability
								// that a gene in the indiv is from one or another of the sources
							Numer[s] = C->Lat->Ind[i]->PritInd->Q[s]->v;
						}			
							
						normo += Numer[s];
					END1CYCLE
					
					// then cycle again to normalize those
					CYCLE_s(C->Lat->PritLat)
						Numer[s] /= normo;
					END1CYCLE
					
					// then draw a new w for it
					C->Lat->Ind[i]->W[l][c]->v = IntFromProbsRV(Numer,0,C->Lat->PritLat->NumSources->v);
				
					// and increment the appropriate NumW
					C->Lat->Ind[i]->PritInd->NumW[ C->Lat->Ind[i]->W[l][c]->v  ]->v++;
					
				}   // end cycle over the two gene copies "c"
			
			}   // closes the CODOM routine conditional
			else {
				printf("\nSorry my friend.  You are using a non-codominant marker");
				printf("\nwith the Pritchard method.  No can do!  Exiting...\n\n");
				exit(1);
			}
			
		END1CYCLE	// loci
				
	END1CYCLE  // individuals
	
}


// little function to update all the Q's of the individuals.  Pretty straightforward---
// Q has a Dirichlet dsn given the "prior" (alpha) and the "data" w.
// You have to call PritchUpdateW before this, so that the NumW's are stored properly
void PritchUpdateQ(hyb_chain *C)
{
	int i,s;
	double TempRV[MAX_NUM_SOURCES_PRITCH], TempPars[MAX_NUM_SOURCES_PRITCH];  // for interfacing with the DirichletRV function
	
	// cycle over individuals
	CYCLE_i(C->Dat)
		
		// cycle over s to set the parameters (NumW + alpha's)
		CYCLE_s(C->Lat->PritLat)
			TempPars[s] = (double)C->Lat->Ind[i]->PritInd->NumW[s]->v + C->Lat->PritLat->Alpha[s]->v;	
		END1CYCLE
		
		// draw the new random vector which is the Q
		DirichletRV(TempPars,C->Lat->PritLat->NumSources->v, TempRV);
		
		// then cycle over s again to copy from TempRV back into Q
		CYCLE_s(C->Lat->PritLat)
			C->Lat->Ind[i]->PritInd->Q[s]->v = TempRV[s];	
		END1CYCLE
		
	END1CYCLE
	
	// that's it!  Gee, that was easy!
	
}


Pritch_Indiv_Latent_struct *AllocPritchInd(int NumSources)
{
	Pritch_Indiv_Latent_struct *temp;
	
	temp = (Pritch_Indiv_Latent_struct *)ECA_MALLOC(sizeof(Pritch_Indiv_Latent_struct));
	
	temp->Q = DvalVector(0,NumSources-1,0.0,1.0,.02);
	temp->NumW = IvalVector(0,NumSources-1,0,0,0);
	
	temp->GtypProb = DvalVector(0,1,0.0,0.0,0.0);
	
	temp->V = AllocIval(0,0,0);
	temp->PofV = DvalVector(0,1,0.0,1.0,.02);
	
	return(temp);
}

// incomplete function still
Pritch_Overall_Latent_struct *AllocPritchOverall(int NumSources)
{
	int i;
	Pritch_Overall_Latent_struct *temp;
	double Min[] = {0.0,0.0,0.0,0.0};
	double Max[] = {3.0, 3.0, 3.0, 3.0};
	double SD[] = {.1, .1, .1, .1};
	
	temp = (Pritch_Overall_Latent_struct *)ECA_MALLOC(sizeof(Pritch_Overall_Latent_struct));
	
	temp->NumSources = AllocIval(0,0,0);
	// then set to initial value
	temp->NumSources->v = NumSources;
	
	// allocate to and set values of the mins and maxes for alpha.
	temp->MinAlpha = DvalVector(0,NumSources-1,0.0,0.0,0.0);
	temp->MaxAlpha = DvalVector(0,NumSources-1,0.0,0.0,0.0);
	temp->Sigma_SD = DvalVector(0,NumSources-1,0.0,0.0,0.0);
	for(i=0;i<NumSources;i++)  {
		temp->MinAlpha[i]->v = Min[i];
		temp->MaxAlpha[i]->v = Max[i];
		temp->Sigma_SD[i]->v = SD[i];
	}
	
	
	temp->Alpha = DvalVector(0,NumSources-1,Min[0],Max[0],(Max[0]-Min[0])/50.0);
	
	// then just set them all to one for now
	for(i=0;i<NumSources;i++)  {
		temp->Alpha[i]->v = 1.0;
	}
	
	// allocated to and Xi
	temp->Xi = DvalVector(0,1,0.0,1.0,1.0/50.0);
	// for now, set them all to .5
	temp->Xi[0]->v = .5;
	temp->Xi[1]->v = .5;
	
	temp->Vcounts = IvalVector(0,1,0,0,0);
	
	
	
	return(temp);
}

// do all the allocation to the pritchard stuff
void AllocatePritchardParts(hyb_chain *C)
{
	int i;
	
	// get the overall stuff allocated to
	C->Lat->PritLat = AllocPritchOverall(2);
	
	// get the individuals allocated to
	CYCLE_i(C->Dat)
		C->Lat->Ind[i]->PritInd = AllocPritchInd(2);
	END1CYCLE
	
	
}



// initialize the variables that are specific to the pritchard method
void InitializePritchVars(hyb_chain *C)
{
	int i,s,l,k,jv;
	double InArray[MAX_NUM_SOURCES_PRITCH], tempPar1[2],tempPar2[2];
	double OutArray[MAX_NUM_SOURCES_PRITCH], tempRV1[2], tempRV2[2];
	double tempTheta[MAX_ALLELES], tempPrior[MAX_ALLELES];
	
	
	// draw thetas from their priors
	// cycle over species and loci, initializing Theta
	CYCLE_s(C->Lat->PritLat)
		CYCLE_l(C->Dat) 
			// store the priors in the temporary array
			CYCLE_k(C->Dat)
				tempPrior[k] = C->Pri->Lambda[s][l][k]->v;
			END1CYCLE
			// simulate the new value
			DirichletRV(tempPrior,C->Dat->Kl[l],tempTheta);
			
			// copy it into the theta struct
			CYCLE_k(C->Dat)
				C->Lat->Theta[s][l][k]->v = tempTheta[k];
			END1CYCLE
		END1CYCLE
	END1CYCLE
	
	
	// draw new alphas from their uniform prior
	CYCLE_s(C->Lat->PritLat)
		C->Lat->PritLat->Alpha[s]->v = C->Lat->PritLat->MinAlpha[s]->v + (double)ranf() *
			C->Lat->PritLat->MaxAlpha[s]->v;
		if(gAlphasConstrainedEqual==1)  {  // if they are constrained equal, then set
										   // them all equal to the first
										   
			C->Lat->PritLat->Alpha[s]->v = C->Lat->PritLat->Alpha[0]->v;
		}
	END1CYCLE
	
	
	
	
	// draw values for the Q's
	CYCLE_i(C->Dat)
		// fill the parameter array for the dirichlet dsn
		CYCLE_s(C->Lat->PritLat)
			InArray[s] = C->Lat->PritLat->Alpha[s]->v;
		END1CYCLE
		
		// draw the new dirichlet rv
		DirichletRV(InArray, C->Lat->PritLat->NumSources->v, OutArray);
		
		// then copy over the result
		CYCLE_s(C->Lat->PritLat)
			C->Lat->Ind[i]->PritInd->Q[s]->v = OutArray[s];
		END1CYCLE
	END1CYCLE
	
	// draw new Pi's and Xi's from their priors, and also 
	// new V's and Z's while we are at it.
	// these will always only have two components
	// cycle to store in the InArray
	CYCLE_jv	
		tempPar1[jv] = C->Pri->XiPrior[jv]->v;
		tempPar2[jv] = C->Pri->Zeta[jv]->v;
	END1CYCLE
	// then draw the new RV's
	DirichletRV(tempPar1, 2, tempRV1);
	DirichletRV(tempPar2, 2, tempRV2);
	// New Z's and V's:
	// then simulate Z for each individual by drawing directly from Pi
	CYCLE_i(C->Dat)
		C->Lat->Ind[i]->Z->v = IntFromProbsRV(tempRV2,0,2);
		C->Lat->Ind[i]->PritInd->V->v = IntFromProbsRV(tempRV1,0,2);
	END1CYCLE
	
	// then cycle and copy Xi and Pi from the temporary arrays
	CYCLE_jv	
		C->Lat->PritLat->Xi[jv]->v = tempRV1[jv];
		C->Lat->Pi[jv]->v = tempRV2[jv];
	END1CYCLE
	
	
}

// here is a little function to update the alpha variable in the Pritchard et al.
// method.  Metropolis-update of all the components by making independent
// proposals to each.  Pretty 
// straightforward stuff
/*void PritchUpdateAlpha(hyb_chain *C, int comp)
{
	int s;
	double a[MAX_NUM_SOURCES_PRITCH];  // the proposed values
	
	// cycle over and propose the values
	CYCLE_s(C->Lat->PritLat)
		
	END1CYCLE
}
*/



/*
	updates the r-th component of alpha by drawing a folded normal for a proposal then
	accepting or not based on the M-H ratio.  
*/
void PritchUpdateAlphaRComp(hyb_chain *C, int r)
 {
 	int i,s;
 	double NewAlpha,OldAlpha; 
 	double rando;
 	double numerator, denominator, ratio;
 	double std_dev;
 	double 	UPPER,  
 			LOWER;   // the upper and lower bounds we want Alpha to range in. 
 			// Note that there will be problems iF std_dev is close to UPPER-LOWER.
 			// I am going to not worry about that now and just keep std_dev small and 
 			// UPPER rather large.
 	double forward_density, backward_density;  // the normal densities for the forward proposal and the backward proposal
 	double CurrentPars[MAX_NUM_SOURCES_PRITCH];
 	double ProposedPars[MAX_NUM_SOURCES_PRITCH];
 	double tempQ[MAX_NUM_SOURCES_PRITCH];
 	
 		
 	UPPER = C->Lat->PritLat->MaxAlpha[r]->v;
 	LOWER = C->Lat->PritLat->MinAlpha[r]->v;
 	std_dev = C->Lat->PritLat->Sigma_SD[r]->v;
 	OldAlpha = C->Lat->PritLat->Alpha[r]->v;
 	
 	// propose a new Alpha between upper and lower:
 	NewAlpha = gennor((float)OldAlpha,(float)std_dev);
 	
 	// Reflect (fold) it if necessary
 	if(NewAlpha < LOWER) {
 		NewAlpha = LOWER + (LOWER - NewAlpha);
 	}
 	if(NewAlpha > UPPER)  {
 		NewAlpha = UPPER - (NewAlpha-UPPER);
 	}
 	
 	// copy the proposed and current values around
 	CYCLE_s(C->Lat->PritLat)
 		CurrentPars[s] = C->Lat->PritLat->Alpha[s]->v;
 		ProposedPars[s] = C->Lat->PritLat->Alpha[s]->v;
 	END1CYCLE
 	// here we make the change to ProposedPars
 	ProposedPars[r] = NewAlpha;
 	
 	if(gAlphasConstrainedEqual==1)  {   // we are constraining them to all be the same...
 		CYCLE_s(C->Lat->PritLat)  
 			ProposedPars[s] = NewAlpha;
 		END1CYCLE
 	}
 	
 	// now we compute the numerator and denominator of the ratio of limit dsns part of
 	// the MH ratio.  Each will be the sum of logs
 	numerator = 0.0;
 	denominator = 0.0;
 	
  	CYCLE_i(C->Dat)  // cycle over all the fish in the mixture
 		 		
 		 	// store the fish's Q's in the temp array of doubles
 		 	CYCLE_s(C->Lat->PritLat)
 		 		tempQ[s] = C->Lat->Ind[i]->PritInd->Q[s]->v;
 		 	END1CYCLE
 		 	
	 		denominator += LogDirichletPDF(CurrentPars, tempQ, C->Lat->PritLat->NumSources->v);
	 		
	 		
	 		// then add the proper part to the numerator:
	 		numerator += LogDirichletPDF(ProposedPars, tempQ, C->Lat->PritLat->NumSources->v);
	 		
	 		
 	END1CYCLE
 
 	
 	// Now we have to include the ratio of the proposals
 	forward_density = NormalPDF(OldAlpha, std_dev*std_dev,   NewAlpha) +   // the straightforward part
 					  NormalPDF(OldAlpha, std_dev*std_dev,   2.0*LOWER - NewAlpha) + // the part reflected around LOWER
 					  NormalPDF(OldAlpha, std_dev*std_dev,   2.0*UPPER - NewAlpha);  // the part reflected around UPPER
 	   
 	backward_density = NormalPDF(NewAlpha, std_dev*std_dev,   OldAlpha) +   // the straightforward part
 					   NormalPDF(NewAlpha, std_dev*std_dev,   2.0*LOWER - OldAlpha) + // the part reflected around LOWER
 					   NormalPDF(NewAlpha, std_dev*std_dev,   2.0*UPPER - OldAlpha);  // the part reflected around UPPER
 	
 	// include these in the MH ratio.  The forward_density part goes on bottom, of course.
 	numerator += log(backward_density);
 	denominator += log(forward_density);
 	
 	// now, down at the end of this, we decide if we accept the proposal or not.
 	rando = log(ranf());
 	ratio = numerator - denominator;
 	
 	if(rando < ratio) {  // then we accept the change
 		CYCLE_s(C->Lat->PritLat)
	 		C->Lat->PritLat->Alpha[s]->v = ProposedPars[s];
 		END1CYCLE
 	}
 	
 	
 }
 
 
 
 /*
	updates all components of alpha, assuming that they are all the same.  Does so with a proposal distribution
	that rejects automatically if the proposal is outside the bounds.  
*/
void PritchUpdateAlphaSymmet(hyb_chain *C)
 {
 	int i,s;
 	double NewAlpha,OldAlpha; 
 	double rando;
 	double numerator, denominator, ratio;
 	double std_dev;
 	double 	UPPER,  
 			LOWER;   // the upper and lower bounds we want Alpha to range in. 
 			// Note that there will be problems iF std_dev is close to UPPER-LOWER.
 			// I am going to not worry about that now and just keep std_dev small and 
 			// UPPER rather large.
 	double forward_density, backward_density;  // the normal densities for the forward proposal and the backward proposal
 	double CurrentPars[MAX_NUM_SOURCES_PRITCH];
 	double ProposedPars[MAX_NUM_SOURCES_PRITCH];
 	double tempQ[MAX_NUM_SOURCES_PRITCH];
 	
 		
 	UPPER = 2;
 	LOWER = 0;
 	std_dev = .05;
 	OldAlpha = C->Lat->PritLat->Alpha[0]->v;
 	
 	// propose a new Alpha between upper and lower:
 	NewAlpha = gennor((float)OldAlpha,(float)std_dev);
 	
 	// Reflect (fold) it if necessary
 	if(NewAlpha < LOWER || NewAlpha > UPPER) {
 		return;
 	}
 	
 	// copy the proposed and current values around
 	CYCLE_s(C->Lat->PritLat)
 		CurrentPars[s] = C->Lat->PritLat->Alpha[s]->v;
 		ProposedPars[s] = NewAlpha;
 	END1CYCLE
 	
 	
 	// now we compute the numerator and denominator of the ratio of limit dsns part of
 	// the MH ratio.  Each will be the sum of logs
 	numerator = 0.0;
 	denominator = 0.0;
 	
  	CYCLE_i(C->Dat)  // cycle over all the fish in the mixture
 		 		
 		 	// store the fish's Q's in the temp array of doubles
 		 	CYCLE_s(C->Lat->PritLat)
 		 		tempQ[s] = C->Lat->Ind[i]->PritInd->Q[s]->v;
 		 	END1CYCLE
 		 	
	 		denominator += LogDirichletPDF(CurrentPars, tempQ, C->Lat->PritLat->NumSources->v);
	 		
	 		
	 		// then add the proper part to the numerator:
	 		numerator += LogDirichletPDF(ProposedPars, tempQ, C->Lat->PritLat->NumSources->v);
	 		
	 		
 	END1CYCLE
 
 	
 	// Now we have to include the ratio of the proposals
 	forward_density = NormalPDF(OldAlpha, std_dev*std_dev,   NewAlpha);
 	   
 	backward_density = NormalPDF(NewAlpha, std_dev*std_dev,   OldAlpha);
 	
 	// include these in the MH ratio.  The forward_density part goes on bottom, of course.
 	numerator += log(backward_density);
 	denominator += log(forward_density);
 	
 	// now, down at the end of this, we decide if we accept the proposal or not.
 	rando = log(ranf());
 	ratio = numerator - denominator;
 	
 	if(rando < ratio) {  // then we accept the change
 		CYCLE_s(C->Lat->PritLat)
	 		C->Lat->PritLat->Alpha[s]->v = ProposedPars[s];
 		END1CYCLE
 	}
 	
 	
 }


void PritchUpdateRandomCompOfAlpha(hyb_chain *C)
{
	int r;
	
	// first choose r at random:
	r = UniformRV(0,C->Lat->PritLat->NumSources->v-1);
	
	// then update it
	PritchUpdateAlphaRComp(C,r);
}



void IncrementPritchVars(hyb_chain *C, int DoThetas)
{
	int i,s,a,l,k,jv;
	
	
	// increment the thetas
	if(DoThetas)  {
		CYCLE_alk(C->Dat)
			IncrementDval(C->Lat->Theta[a][l][k]);
		END3CYCLE
	}
	
	CYCLE_i(C->Dat)
	
		// only increment the Q's of those individuals that are JA_MIXED or if the simulation is
		// includes no pure individuals
		if(C->Lat->Ind[i]->PritInd->V->v == JA_ADMIXED || TypeOfPritchEtcSim == PRITCHARD
					 || TypeOfPritchEtcSim == JABES_NO_PURES)  {
			CYCLE_s(C->Lat->PritLat)  
				IncrementDval(C->Lat->Ind[i]->PritInd->Q[s]);
			END1CYCLE
		}
		
		// and increment the PofV's
		CYCLE_jv
			IncrementDval(C->Lat->Ind[i]->PritInd->PofV[jv]);
		END1CYCLE
	END1CYCLE
	
	// increment the alphas
	CYCLE_s(C->Lat->PritLat)
		IncrementDval(C->Lat->PritLat->Alpha[s]);
	END1CYCLE
	
	// increment xi and pi
	CYCLE_jv
		IncrementDval(C->Lat->PritLat->Xi[jv]);
		IncrementDval(C->Lat->Pi[jv]);
	END1CYCLE
	
	// increment the kullback-leibler distances on the loci, too
	CYCLE_l(gC->Dat)
		IncrementDval(C->Lat->Locus_KB[l]);
	END1CYCLE
	

	
}



// Reset all the averages and things for all of the variables in a chain
// Currently it only deals with the few variables I have been keeping track of
void PritchResetAllAveragesEtc(hyb_chain *C)
{
	int a,l,k,i,s,jv;
	// first do all the Thetas
	CYCLE_alk(C->Dat)
		InitDvalSummaryToZero(C->Lat->Theta[a][l][k]);
	END3CYCLE
	
	
	// reset all the Q's...
	CYCLE_i(C->Dat)
		CYCLE_s(C->Lat->PritLat)  
			InitDvalSummaryToZero(C->Lat->Ind[i]->PritInd->Q[s]);
		END1CYCLE
		
		// and reset all the PofV's
		CYCLE_jv
			InitDvalSummaryToZero(C->Lat->Ind[i]->PritInd->PofV[jv]);
		END1CYCLE
	END1CYCLE
	
	
	// reset all the Alpha's
	CYCLE_s(C->Lat->PritLat)  
		InitDvalSummaryToZero(C->Lat->PritLat->Alpha[s]);
	END1CYCLE
	
	// reset the Xi's and the Pi's
	CYCLE_jv
		InitDvalSummaryToZero(C->Lat->PritLat->Xi[jv]);
		InitDvalSummaryToZero(C->Lat->Pi[jv]);
	END1CYCLE

/*	
	// then initialize the sliding trace structs, too
	gfduInitSlidingTraceToZero(C->Lat->CompleteDataLogLike);
*/
	
	// reset the averages of the kullback-leibler stuff too
	CYCLE_l(C->Dat)
		InitDvalSummaryToZero(C->Lat->Locus_KB[l]);
	END1CYCLE

}


// updates the W's by the Baum method
// AND IT FILLS OUT NumWs as well, so that the next call to PritchUpdateQ will work.   
void UpdateWbyBaum(hyb_chain *C)
{
	int i,t,l,c,s;
	int NumZs = C->Dat->L * 2;
	MAZ_var_struct MV;
	
	if(C->Lat->PritLat->NumSources->v != 2)  {
		printf("\n\nCalling UpdateWbyBaum when NumSources = %d, Not 2",C->Lat->PritLat->NumSources);
		printf("\nExiting to System!...\n\n");
		exit(1);
	}
	
	// set the C field in MV
	MV.C = C;
	
	
	// cycle over all the individuals
	CYCLE_i(C->Dat)
	
		// set the Indiv field in MV for the current individual i
		MV.Indiv = i;
	
		// for the i-th individual, call MargAllocZs
		C->Lat->Ind[i]->PritInd->GtypProb[JA_ADMIXED]->v =   MargAllocZ(
														NumZs,
														gBaumIdx,
														0,
														1,
														C->Lat->PritLat->Alpha[0]->v,
														C->Lat->PritLat->Alpha[1]->v,
														MargAllocDataProb, 
														DoubleRanf,
														&MV,
														gBaumZOutputParm
													);
													
		// now we have to copy those new values to the individual, and count NumWs while we are at it
		// cycle over the sources and initialize NumW to 0 to accumulate a sum
		CYCLE_s(C->Lat->PritLat)
			C->Lat->Ind[i]->PritInd->NumW[s]->v = 0;
		END1CYCLE
		for(t=0;t<NumZs;t++) {
			l = t/2;
			c = t%2;
			
			// only do this if the indiv is currently JA_ADMIXED or if the type of simulation
			// is JABES_NO_PURE
			if(C->Lat->Ind[i]->PritInd->V->v == JA_ADMIXED  || TypeOfPritchEtcSim == JABES_NO_PURES)  {
				// increment NumWs
				C->Lat->Ind[i]->PritInd->NumW[ gBaumZOutputParm[t] ]->v += 1;
				
				// copy the z into the appropriate w
				C->Lat->Ind[i]->W[l][c]->v = gBaumZOutputParm[t];
			}
			
		}
		
	END1CYCLE
}


// updates all the variables in a single chain
// I do this in the order of W and Z first, then theta and
// pi, because it is easier to initialize Pi and Theta (and the Z's)
// than initializing W's and Z's
void PritchSingleSweep(hyb_chain *C)
{
	int i;
	FILE *out;
	
	
	TypeOfPritchEtcSim = PRITCHARD;
	gAlphasConstrainedEqual = 1;
	// Do all the updates specific for each type of simulation
	
	
	if(TypeOfPritchEtcSim == PRITCHARD)  {
	
		
	
		PritchUpdateW(C);
		PritchUpdateQ(C);
		PritchUpdateAlphaSymmet(C);
		
		
		
		C->Lat->PritLat->Alpha[0]->v= .1;
		C->Lat->PritLat->Alpha[1]->v = .1;
		
		if(gNumSweepsAfterBurnIn % 100 == 0) {
			out = fopen("PritchQsOut.txt","w");
			fprintf(out,"ID     Q");
			CYCLE_i(C->Dat)
				fprintf(out,"\n%d    %.8f",i+1,C->Lat->Ind[i]->PritInd->Q[1]->Ave);
			END1CYCLE
			fclose(out);
		}
	}
//	
	else if(TypeOfPritchEtcSim == JABES_NO_PURES)  {
		UpdateWbyBaum(C);
		PritchUpdateQ(C);
		PritchUpdateRandomCompOfAlpha(C);
	}
	
	else if(TypeOfPritchEtcSim == JABES)  {
		
		
		UpdateWbyBaum(C);
		JABESComputePofZ(C);
		JABES_UpdateZs(C);
		JABES_UpdateWForPure(C);
		PritchUpdateQ(C);
		JABESUpdateRandomCompOfAlpha(C);
		JABES_UpdatePi(C);
		JABESUpdateV(C);
		JABES_UpdateXi(C);
		
		
	}
	
	
	
	

	// these get done the same regardless
	UpdateTheta(C);
	KullbLeib(C);
	
}


double DoubleRanf(void)
{
	return((double)ranf());
}

void JABES_UpdatePi(hyb_chain *C)
{
	int i,g;
	double temp[MAX_GN];
	
	
	// have to count the Z's first, but only on the ones that are currently
	// allocated to the Pure Pile.
	// cycle over 0 and 1 first and set the s's to zero
	CYCLE_g(C->Dat)
		C->Lat->s[g] = 0.0;
	END1CYCLE
	
	CYCLE_i(C->Dat)
		if(C->Lat->Ind[i]->PritInd->V->v == JA_PURE) {
			C->Lat->s[  C->Lat->Ind[i]->Z->v ] += 1.0;
		}
	END1CYCLE
	
	// then  add the priors
	// then add the priors to those
	CYCLE_g(C->Lat)
		C->Lat->s[g] += C->Pri->Zeta[g]->v;
	END1CYCLE
	
	// then simulate a new value
	DirichletRV(C->Lat->s,C->Lat->Gn->v,temp);
	
	// and copy it into the Pi struct
	CYCLE_g(C->Lat)
		C->Lat->Pi[g]->v = temp[g];
	END1CYCLE
}




void JABES_UpdateXi(hyb_chain *C)
{
	double tempPars[2], tempRVs[2];
	int jv;
	
	// first count up the Vcounts
	JABEScount_Vcounts(C);
	
	// then make the tempPars by adding the XiPriors to that
	CYCLE_jv
		tempPars[jv] = (double)C->Lat->PritLat->Vcounts[jv]->v + C->Pri->XiPrior[jv]->v;
	END1CYCLE
	
	// then draw the Dirichlet random vector
	DirichletRV(tempPars,2,tempRVs);
	
	// then cycle one last time to copy those values over:
	CYCLE_jv
		C->Lat->PritLat->Xi[jv]->v = tempRVs[jv];
	END1CYCLE
	
}


// count up the number of individuals in the pure vs admixed categories
void JABEScount_Vcounts(hyb_chain *C)
{
	int i,jv;
	
	// initialize to accumulate a sum
	CYCLE_jv
		C->Lat->PritLat->Vcounts[jv]->v = 0;
	END1CYCLE
	
	// now cycle over the individuals and count them up
	CYCLE_i(C->Dat)
		C->Lat->PritLat->Vcounts [ C->Lat->Ind[i]->PritInd->V->v]->v++;
	END1CYCLE
}

// for each individual, this computes the probability of its genotype given 
// that it is pure.  This quantity gets stored in the PritInd field GtypProb[0].
// at the same time, the gtyp probs that it is pure 0 or pure 1 get recorded, AND
// the probability, given Pi that it is pure 0 or pure 1 are recorded in 
// PofZ
void JABESComputePofZ(hyb_chain *C) 
{
	int i,g;
	double normo;
	
	CYCLE_i(C->Dat)
		// compute the genotype probs (and the full conditionals for the w's)
		// for each individual
		
		// the following puts the values into LogPofY[g] for the individual
		for(g=0;g<2;g++)  {
			JABES_StoreLogPofY_Pure(C,i,g);
		}
		
		
		// then we must apply Bayes law to use PofY and Pi to compute 
		// PofZ.  In the future I may want to do some fancy stuff to prevent
		// underflow here (i.e. add in and later pull out a constant of log something before
		// something before exponentiating LogPofY to become a probability). 
		
		normo = 0.0;  // initialize to accumulate a sum of probabilities
		
		// cycle over the two populations
		for(g=0;g<2;g++)  {
			// first store them as unnormalized probabilities
			C->Lat->Ind[i]->PofZ[g]->v = exp(C->Lat->Ind[i]->LogPofY[g]->v) * 
										 C->Lat->Pi[g]->v;
										 
			// and add that to the normalization factor
			normo += C->Lat->Ind[i]->PofZ[g]->v;
		}
		
		// down here, normo holds the prob of the indiv's gtyp given he is PURE
		C->Lat->Ind[i]->PritInd->GtypProb[JA_PURE]->v = normo;
		
		// now cycle again over g and normalize PofZ
		for(g=0;g<2;g++)  {
			// divide each one by normo
			C->Lat->Ind[i]->PofZ[g]->v /= normo;
		}
		
	END1CYCLE
}


// updates the W's in the individuals with V == JA_PURE.
// this is very simple.  It just uses the value of the individual's
// Z and fills that in for the W's
void JABES_UpdateWForPure(hyb_chain *C)
{
	int l,i,c;
	
	CYCLE_i(C->Dat)
		if(C->Lat->Ind[i]->PritInd->V->v == JA_PURE) {  // only do this for individuals that are currently pure
			CYCLE_l(C->Dat)
				for(c=0;c<2;c++)  {
					C->Lat->Ind[i]->W[l][c]->v = C->Lat->Ind[i]->Z->v;
				}
			
			END1CYCLE
		}  // closes the "if pure" conditional
	END1CYCLE
}





// update the Z in each individual, that is currently JA_PURE.
// Requires that PofZ has been set already (usually done in JABES_ComputePofZ.
void JABES_UpdateZs(hyb_chain *C)
{
	int i;
	CYCLE_i(C->Dat)
		if(C->Lat->Ind[i]->PritInd->V->v == JA_PURE)  {
			C->Lat->Ind[i]->Z->v = IntegerFromProbsRV(C->Lat->Gn->v,C->Lat->Ind[i]->PofZ);
		}
	
	END1CYCLE
	
}

// this just normalizes the GtypProbs by weighting by the Xi
// assumes that the GtypProbs have already been computed
void JABESComputePofV(hyb_chain *C)
{
	int jv,i;  
	double normo;
	
	CYCLE_i(C->Dat)
		normo = 0.0;  // initialize to accumulate a sum
		// cycle over once and fill with the unnormalized values
		CYCLE_jv
			C->Lat->Ind[i]->PritInd->PofV[jv]->v = C->Lat->PritLat->Xi[jv]->v *
												C->Lat->Ind[i]->PritInd->GtypProb[jv]->v;
												
			normo += C->Lat->Ind[i]->PritInd->PofV[jv]->v;
		END1CYCLE
		
		// then cycle again and normalize
		CYCLE_jv
			C->Lat->Ind[i]->PritInd->PofV[jv]->v /= normo;
		END1CYCLE
		
	END1CYCLE
}


// this samples new values of the V's.  It does this by first calling
// JABESComputePofV, then drawing new V's.  Remember, this assumes that
// there has already been a call to JABESComputePofZ
void JABESUpdateV(hyb_chain *C)
{
	int i;
	
	
	// first compute PofV for everyone
	JABESComputePofV(C);
	
	// then cycle over the individuals and draw new V's
	CYCLE_i(C->Dat)
		C->Lat->Ind[i]->PritInd->V->v = IntegerFromProbsRV(2,C->Lat->Ind[i]->PritInd->PofV);
	END1CYCLE
	
}	


/*
	updates the r-th component of alpha by drawing a folded normal for a proposal then
	accepting or not based on the M-H ratio.  This uses the Q carried only by the indivs
	with V == JA_ADMIXED.  This is a silly function.  Will soon replace it with 
	one that does the MH ratio more better (as pointed out in JABES paper)
*/
void JABESUpdateAlphaUsingQ(hyb_chain *C, int r)
 {
 	int i,s;
 	double NewAlpha,OldAlpha; 
 	double rando;
 	double numerator, denominator, ratio;
 	double std_dev;
 	double 	UPPER,  
 			LOWER;   // the upper and lower bounds we want Alpha to range in. 
 			// Note that there will be problems iF std_dev is close to UPPER-LOWER.
 			// I am going to not worry about that now and just keep std_dev small and 
 			// UPPER rather large.
 	double forward_density, backward_density;  // the normal densities for the forward proposal and the backward proposal
 	double CurrentPars[MAX_NUM_SOURCES_PRITCH];
 	double ProposedPars[MAX_NUM_SOURCES_PRITCH];
 	double tempQ[MAX_NUM_SOURCES_PRITCH];
 	
 		
 	UPPER = C->Lat->PritLat->MaxAlpha[r]->v;
 	LOWER = C->Lat->PritLat->MinAlpha[r]->v;
 	std_dev = C->Lat->PritLat->Sigma_SD[r]->v;
 	OldAlpha = C->Lat->PritLat->Alpha[r]->v;
 	
 	// propose a new Alpha between upper and lower:
 	NewAlpha = gennor((float)OldAlpha,(float)std_dev);
 	
 	// Reflect (fold) it if necessary
 	if(NewAlpha < LOWER) {
 		NewAlpha = LOWER + (LOWER - NewAlpha);
 	}
 	if(NewAlpha > UPPER)  {
 		NewAlpha = UPPER - (NewAlpha-UPPER);
 	}
 	
 	// copy the proposed and current values around
 	CYCLE_s(C->Lat->PritLat)
 		CurrentPars[s] = C->Lat->PritLat->Alpha[s]->v;
 		ProposedPars[s] = C->Lat->PritLat->Alpha[s]->v;
 	END1CYCLE
 	// here we make the change to ProposedPars
 	ProposedPars[r] = NewAlpha;
 	
 	if(gAlphasConstrainedEqual==1)  {   // we are constraining them to all be the same...
 		CYCLE_s(C->Lat->PritLat)  
 			ProposedPars[s] = NewAlpha;
 		END1CYCLE
 	}
 	
 	
 	
 	// now we compute the numerator and denominator of the ratio of limit dsns part of
 	// the MH ratio.  Each will be the sum of logs
 	numerator = 0.0;
 	denominator = 0.0;
 	
  	CYCLE_i(C->Dat)  // cycle over all the fish in the mixture
 		 		
 		 if(C->Lat->Ind[i]->PritInd->V->v == JA_ADMIXED)  {  // only do this for admixed ones.
 		 	// store the fish's Q's in the temp array of doubles
 		 	CYCLE_s(C->Lat->PritLat)
 		 		tempQ[s] = C->Lat->Ind[i]->PritInd->Q[s]->v;
 		 	END1CYCLE
 		 	
	 		denominator += LogDirichletPDF(CurrentPars, tempQ, C->Lat->PritLat->NumSources->v);
	 		
	 		
	 		// then add the proper part to the numerator:
	 		numerator += LogDirichletPDF(ProposedPars, tempQ, C->Lat->PritLat->NumSources->v);
	 	}
	 		
	 		
 	END1CYCLE
 	
 	// Now we have to put include the ratio of the proposals
 	forward_density = NormalPDF(OldAlpha, std_dev*std_dev,NewAlpha) +   // the straightforward part
 					  NormalPDF(OldAlpha, std_dev*std_dev,2.0*LOWER - NewAlpha) + // the part reflected around LOWER
 					  NormalPDF(OldAlpha, std_dev*std_dev,2.0*UPPER - NewAlpha);  // the part reflected around UPPER
 	
 	backward_density = NormalPDF(NewAlpha, std_dev*std_dev,OldAlpha) +   // the straightforward part
 					   NormalPDF(NewAlpha, std_dev*std_dev,2.0*LOWER - OldAlpha) + // the part reflected around LOWER
 					   NormalPDF(NewAlpha, std_dev*std_dev,2.0*UPPER - OldAlpha);  // the part reflected around UPPER
 	
 	// include these in the MH ratio.  The forward_density part goes on bottom, of course.
 	numerator += log(backward_density);
 	denominator += log(forward_density);
 	
 	// now, down at the end of this, we decide if we accept the proposal or not.
 	rando = log(ranf());
 	ratio = numerator - denominator;
 	
 	if(rando < ratio) {  // then we accept the change
 		CYCLE_s(C->Lat->PritLat)
	 		C->Lat->PritLat->Alpha[s]->v = ProposedPars[s];
 		END1CYCLE
 	}
 	
 	
 }

// currently uses JABESUpdateAlphaUsingQ.  Must update when I have 
// the better function to use.
void JABESUpdateRandomCompOfAlpha(hyb_chain *C)
{
	int r;
	
	// first choose r at random:
	r = UniformRV(0,C->Lat->PritLat->NumSources->v-1);
	
	// then update it
	JABESUpdateAlphaUsingQ(C,r);
}


// returns the probabability of the i-th individual's gtyp given it is PURE and from
// source a.
void JABES_StoreLogPofY_Pure(hyb_chain *C, int i, int a)
{
	int c,l;
	double LogProb = 0.0;  // initialize to accumulate a sum
	
	
	CYCLE_l(C->Dat)
		
		if(YnotMissing(C,i,l))  {  // only count non-missing loci
			for(c=0;c<2;c++)  {
				LogProb += log(C->Lat->Theta[a][l][ C->Lat->Ind[i]->Y[l][c]->v ]->v);
			}
		}
	
	END1CYCLE
	
	// at the end of that, we store the value in LogPofY[a]
	C->Lat->Ind[i]->LogPofY[a]->v = LogProb;
}



