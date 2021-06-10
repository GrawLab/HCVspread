// ###################################################################
// PROJECT: SPATIAL SIMULATION OF VIRAL SPREAD WITHIN A 2D LAYER OF CELLS
//
// Main file for the simulations -> called from R
// @Authors: Peter Kumberger & Frederik Graw
// ###################################################################

// ============================================
// Standard header-files

#include <stdio.h>
#include <iostream>
#include <fstream>
#include <list>
#include <stdlib.h>
#include <math.h>
#include <R.h>
#include <Rmath.h>
#include <algorithm>    // std::random_shuffle
#include <vector>       // std::vector
#include <ctime>        // std::time
#include <cstdlib>      // std::rand, std::srand
#include <iterator>     // std::distance

using namespace std;

// Include the corefile with all the main functions + parameters previously defined
// (paths may have to be adapted)
#include HCVSpread_corefile.h

// ============================================


/* all parameters given below should be taken from an overview-file allowing
 different simulation runs:
 sim = number of different simulations, timesteps = maximal number of timesteps one simulation will run,
 prob_budDIFF = probability of an infected hepatocyte to release a viron to the diffusing viral load,
 prob_infectCC = probability to infect by cell-to-cell transmission,
 prob_infectDIFF = probability of a free diffusing viral particle to infect a hepatocyte
 cellcell = fraction of cell-to-cell transmission of total mode of transmission (cellcell=1 - all transmission via cell-to-cell)
// replTime = fixed time a hepatocyte needs to create a virion (after budding/infection this time is taken for reduction)
 prob_dieUI,I,J = probability of uninfected,infected and infectious hepatocytes to die
 prob_infinf = probability of an infected cell to become infectious
 max_init_nbr = maximal nbr of infectious hepatocytes after initiation time
 init_exp_rate = Rate for initiation of infectious hepatocytes
 trunc_time = Time at which initiation of infection is over
 coupling = rate of coupling between adjacent sites of the viral grid (see Funk et al., JTB 2005)
 rate_dieV = death rate of free diffusing viral particles
 runsV = loops through the viral grid to allow for more widespread diffusion
*/


// =============================================
// MAIN FUNCTION FOR R - called by the R-file
// =============================================
extern "C"{
    void Simulate_infection(int *sim,int *timesteps,double *dummy,double *dummyDIST,double *prob_infectCC,double *prob_infectDIFF,double *frac_HCV,double *prob_dieUI,double *prob_dieJ,double *prob_dieI,
                            double *prob_infinf,int *max_init_nbr,double *init_exp_rate,int *trunc_time,int *media_change,int *prolif_time,double *perc_empty,double *coupling,double *rate_dieV,
                            double *runsV, double *R0, double *prodRNA,double *Rcap,double *rate_exportRNA,double *rate_degradeRNA, double *E2_0, double *degr_e2, int *initINF
							){
		/* write the savefile for the hepatocytes and virions and the corresponding headers */
        Rprintf("initINF: %d \n",*initINF);
        //	outHEPstat<<"sim"<<","<<"time"<<","<<"uninfected"<<","<<"infected"<<","<<"infectious"<<","<<"dead"<<","<<"empty"<<","<<"viralLoad"<<","<<"IFN"<<endl;
		
		/* write the savefile for the hepatocytes */
//		outHEPINFO<<"sim"<<","<<"time"<<","<<"x"<<","<<"y"<<","<<"id"<<","<<"status"<<","<<"IT"<<","<<"IIT"<<","<<"hcv"<<endl;
        int no_clust=0;
        
        /* initialize hexagonal grid */
        vector<int> grid_id;
        grid_id=initGrid();
        //initGrid();
        Rprintf("Hexagonal Grid initialized");
        /* initialize the 2D hepatic layer */
		initHepLayer(grid_id);
		Rprintf(", Hepatic Layer initialized");
		
		/* initialize the corresponding viral grid */
		initViralGrid(grid_id);
		Vector Loads,infType;
		Loads=setZerovector(Loads);
        infType=setZerovector(infType);
		Rprintf(", Viral Grid initialized \n");
		
		/* introduce list of infectious cells (I) (at the start the length is 0) */
		list<int> I;
		int length_I=0;
		/* introduce list of newly infected cells (J) (all infected cells at the start are newly infected) */
		list<int> J;
        /* introduce list of infectious cells in a state of eclipse after successful cell-to-cell transmission (at the start the length is 0) */
        list<int> K;
        /* introduce list of empty spaces/cells (E) (here, this has length 0) */
        list<int> E;
		/* introduce list of death cells (D)
		list<int> D;*/
	
		/* when should be saved - every 60min*/
		int saveT=60;
		double testS=0;
		int msim=*sim-1;
		// for each simulation
		for (int m=msim;m<*sim;m++){
			
			/* define the vector for the counts */
			Vector stat_count;
			
			/* reset the cell layer anatomy to all hepatocytes uninfected/ and viral concentration to 0 and add empty spaces */
			E=resetHepLayer(grid_id,E,*perc_empty,dummy,*prolif_time);
			resetViralGrid(grid_id);
						
            int t=*initINF;                         // set the clock to initINF (initiation time of infection) and start simulation run
            saveList(I,m,t,1);
            
            int length_infectious=0;
            double E2=*E2_0;     // initialize E2 concentration
            
            while (t<*timesteps){
                /* initialize infection (only one infected per timestep) and save them into a list */
                if (t <= *trunc_time) {
                    double prob_init;
                    prob_init = *max_init_nbr*probInit(t,*init_exp_rate,*trunc_time);
                    randbin(dummy,prob_init);
                    //outTEST<<t<<","<<*dummy<<endl;
                    if(*dummy==1){
                        J=initInfected(1,J,grid_id,&no_clust,t);
                    }
                }
                
                //* --------------------------------------------------------------------
                // PROLIFERATION OF UNINFECTED HEPATOCYTES
                E=prolif_cells(E,dummy,*prolif_time,t);
				
				//* --------------------------------------------------------------------
				// INFECTIONS BY CELL-TO-CELL TRANSMISSION AND DIFFUSING VIRUS
				
				/* possible infection of hepatocytes */
				length_infectious=I.size();
				
				/* (a) by cell-to-cell transmission */
				if (length_infectious>0){
                    I=updateIntraVirusConc(I,t,*prodRNA,*Rcap,*rate_exportRNA,*rate_degradeRNA,*frac_HCV);
                    K=updateIntraVirusConc(K,t,*prodRNA,*Rcap,*rate_exportRNA,*rate_degradeRNA,*frac_HCV);
                    /* pointer to list J */
                    list<int> *pJ;
                    pJ = &J;
                    /* pointer to list K */
                    list<int> *pK;
                    pK = &K;
					I=spreadInfection_CC(I,pJ,pK,dummy,dummyDIST,*prob_infectCC,*frac_HCV,t,*R0); // *prob_budDIFF,*replTime,
				}
				/* (b) by diffusing virions */
                J=spreadInfection_DIFF(J,grid_id,*prob_infectDIFF,t,dummy,*R0);
				
				
				//* --------------------------------------------------------------------
				// STATISTICS + UPDATE
				
				/* initialize the different ouput-stream lines */
//				outHEPinf<<m<<","<<t;
//				outHEPdeath<<m<<","<<t;
//				outHEPempty<<m<<","<<t;
				
				/* update the infected and infectious cells/lists */
				length_I=I.size();
				/* pointer to list I */
				list<int> *pI;
				pI = &I;
                /* pointer to list K */
                list<int> *pK;
                pK = &K;
				
                J=updateNewlyInfected_J(J,pK,t,&no_clust,*R0);
                K=updateInfected_K(K,pI,t,*prob_dieJ,*prob_infinf,dummy,*R0);
				
                
                /* save the viral load as well */
                Loads=TotalViralLoad(grid_id,m,t+1,*timesteps+1);
                
                /* Before trunc_time no extracellular virus */
                if (t <= *trunc_time) {
                    resetViralGrid(grid_id);
                }
                else if (t == *media_change) {
                    E2=*E2_0;
                    resetViralGrid(grid_id);
                }
                else {
                    /* update the diffusing viral concentration and calculate the total viral load when t>trunc_time */
                    updateViralConcentration(grid_id,*coupling,*rate_dieV,*runsV,dummy,*E2_0,&E2,*degr_e2);
                }
				
				
				//* --------------------------------------------------------------------
				// SAVE THE OUTPUT
				
				/* save the different lists (I,J,D) */
				testS=(t+1)%saveT;                      // Save only each hour (for t=59, one hour, i.e. 60 min, have passed)
				if (testS==0){
//                    Rprintf("time: %d \n",t);
                    /* update vector to do the statistics of the lattice */
                    stat_count=setZerovector(stat_count);
                    stat_count.vec[0]=m;
                    stat_count.vec[1]=t;
                    
                    /* update the type of infection a cell was infected by*/
                    infType=updateInfectType(I,K);
                    /* go through the grid and update cells according to the type */
                    vector<int>::iterator j;
                    j=grid_id.begin();
                    for (unsigned int i=0;i<grid_id.size();i++){
                        switch ((*hepaticlayer[*j]).status) {
                            case 0: stat_count=updateUninfected(hepaticlayer[*j],stat_count,t,*prob_dieUI,dummy); break;
                            case 4: stat_count.vec[3]++; break;
                            case 2: stat_count.vec[4]++; break;
                            case 3: stat_count.vec[5]++;break;
                        }
                        j++;
                    }
                    //Rprintf("stat vectors to be updated");
                    stat_count.vec[6]=Loads.vec[0];
                    stat_count.vec[7]=Loads.vec[2];
                    stat_count.vec[8]=Loads.vec[1];
                    stat_count.vec[9]=infType.vec[0];
                    stat_count.vec[10]=infType.vec[1];
                    stat_count.vec[11]=infType.vec[2];
                    stat_count.vec[12]=infType.vec[3];
                    stat_count.vec[13]=infType.vec[4];
                    stat_count.vec[14]=E2;
                    
                    //saveVIRALGRID(grid_id,t+1,m);
					saveList(K,m,t+1,4);
					saveList(I,m,t+1,2);
                    saveList(E,m,t+1,3);

					//saveHEPINFO(grid_id,t+1,m); // ,*replTime
					saveHepStat(stat_count);
				}
				//Rprintf("here we are 2 %d \n",t);
				
                if (t==(*timesteps-1)) {
                    Rprintf("Last step \n");
                    saveInfType(grid_id);
                }
                				
				
				/* prepare for the next timestep  - clear list D */
				t++;
				/* make time step one hour */
//				t=t+60;
				
				
			}
			
			/* clear the lists (grid_id,I,J,E) for the next simulation run */
            grid_id.clear();
			J.clear();
			I.clear();
            K.clear();
            E.clear();
		}
	}
}
