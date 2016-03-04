/* "etch.c" */
/* This program takes in a nanoparticle configuration and applies MCMC to it
   in order to "etch" it. 
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_rng.h>
#include "etch_simple.h"

//Use GSL to generate random numbers
gsl_rng *rg;

//Control variables - energies are in units of eV
#define T 0.0259
#define esub 3.930 //sublimation energy
#define epsilon 0.3275 //energy per bond
#define npass 12500000
#define d 1e-6 //small factor used when converting to int
#define num_frames 100 //Don't make this too high or you'll crash the computer.

double mu = (-6.5)*epsilon;
double beta = 1.0/T;
int config_freq = npass/num_frames;
int Nsurf = 0, Nvac = 0; //number of surface and vacany atoms
int Nsurf0 = 0, Nvac0 = 0; //To keep track of initial numbers;
int Ntot = 0;

int temp_mc=0;

int nx, ny, nz;

int main(void){

  printf("%d\n", config_freq);

  int Ngrid;
  int i,j,k;
  long l;
  int time = 0;
  long mc_total = 0; //number of total attempts
  long mc_accept = 0;
  long mc_insert = 0; //number of insertion moves
  long mc_delete = 0; //number of deletion moves

  /* Initialize random number generator *******************************************/
  const gsl_rng_type * t;

  t = gsl_rng_default;
  rg = gsl_rng_alloc (t);

  //set random seed
  long unsigned int seed = 2;
  srand((unsigned) seed);
  gsl_rng_env_setup();
  gsl_rng_set(rg, seed);

  /* Need to load initial configuration ******************************************/

  FILE *in;
  char *input = malloc(1024);

  in = fopen("/home/layne/research/nanoparticle-project/init_config.xyz", "r");
  if(in == NULL) printf("problem!\n");
  if(fscanf(in, "%d\n", &Ngrid) != 1) printf("Error: need 1 number.\n");

  if(fscanf(in, "%s\n", input) == 0) printf("Error: need shape name.\n"); 
  printf("Nanoparticle shape: %s\n", input);

  if(strcmp(input, "cylinder")==0){
    nx = 2*nl;
    ny = 2*nl;
    nz = aspect_ratio*2*nl;
  }
  else{
    nx = 2*nl;
    ny = 2*nl;
    nz = 2*nl;
  } 

  //allocate memory (for lattice)
  nn_vec ***lat;
  lat = (nn_vec ***)calloc((size_t)nx, sizeof(nn_vec**));
  for(i = 0; i < nx; i++){
    lat[i] = (nn_vec **) calloc((size_t)ny, sizeof(nn_vec*));
    for(j = 0; j < ny; j++){
      lat[i][j] = (nn_vec*) calloc((size_t)nz, sizeof(nn_vec));
    }
  }



  char atom_curr[2];
  i = 0;
  j = 0;
  k = 0;
  double temp1;
  double temp2;
  double temp3;

  //Now actually read in data
  for(i = 0; i < nx; i++) {
    for(j = 0; j < ny; j++){
      for(k = 0; k < nz; k++){
	      if(fscanf(in, "%s\t" "%lf\t" "%lf\t" "%lf\t" "%d\t"  "%d\n", atom_curr, &temp1, &temp2, &temp3, &lat[i][j][k].n, &lat[i][j][k].occ) != 6) printf("Error: need 6 inputs.\n");
	      if(lat[i][j][k].occ==1){
	        Ntot++;
	        if(lat[i][j][k].n < NN) Nsurf++;
	      }else if(lat[i][j][k].occ==0){
	          if(lat[i][j][k].n > 0) Nvac++;
	        }
      }
    }
  }
  fclose(in);
  in = NULL;

  printf("Total number of atoms at start: %d\n", Ntot);
  printf("number of surface atoms: %d\n", Nsurf);
  printf("number of vacancies: %d\n", Nvac);

  Nsurf0 = Nsurf;
  Nvac0 = Nvac;

  /* Fill arrays with the indices of the locations of surface and vacany atoms *********************************************************************************/

  triple *surf = (triple *) calloc((size_t)(2*Nsurf), sizeof(triple)); //make these twice the size they need to be to accomodate the possibility
  triple *vac = (triple *) calloc((size_t)(2*Nvac), sizeof(triple)); //of getting more than you started with

  int ls=0; 
  int lv=0;
  for(i = 0; i < nx; i++) {
    for(j = 0; j < ny; j++){
      for(k = 0; k < nz; k++){
      	if((i+j+k)%2 == 0){
	        if(lat[i][j][k].occ==1){
	         if(lat[i][j][k].n < NN){
	           surf[ls].x = i;
	            surf[ls].y = j;
	           surf[ls].z = k;
	           ls++;
	          }
	        }else if(lat[i][j][k].occ==0){
	          if(lat[i][j][k].n > 0){
	            vac[lv].x = i;
	            vac[lv].y = j;
	            vac[lv].z = k;
	            lv++;
	         }
      	  }
      	}
      }
    }
  }


  /* Now do Monte Carlo ****************************************************/
  //  test_one_insert(lat, surf, vac);
  //  printf("number of surface atoms: %d\n", Nsurf);
  //  printf("number of vacancies: %d\n", Nvac);
  //  test_two_delete(lat, surf, vac);
  //  test_two_insert(lat, surf, vac);
  //  test_loneatom(lat, surf, vac);

    
  for(l = 0; l < npass; l++){
    if(l%config_freq==0) printf("Pass number %ld\n", l+1);
      int tmp;
      tmp = mc_move(lat, surf, vac);
    if(tmp == 1){
      mc_insert++;
      mc_accept++;
    }
    if(tmp == 2){
      mc_delete++;
      mc_accept++;
    }
    
    if((l+1)%config_freq == 0){
      print_surf(surf, lat, time, Nsurf);
      print_conf(lat, time, Ntot); 
    	time ++;
    }
    mc_total ++;
  }
   
  //Print information//////////////////////////////////////////////////////////
  printf("Total number of atoms at start: %d\n", Ntot);
  printf("final number of surface atoms: %d\n", Nsurf);
  printf("final number of vacancies: %d\n", Nvac);
  printf("Attempted moves: %ld\n", mc_total);
  printf("Accepted moves: %ld\n", mc_delete + mc_insert);
  printf("Deletion moves: %ld\n", mc_delete);
  printf("Insertion moves: %ld\n", mc_insert);
  printf("Accepted ratio: %lf\n", ((double)mc_delete + (double)mc_insert)/((double)mc_total));

///////////////////////////
//Free pointers
  free(surf);
  free(vac);
  for(i = 0; i < nx; i++){
    for(j = 0; j < ny; j++){
      free(lat[i][j]);
    }
    free(lat[i]);
  }
  free(lat);

  return 1; //end of program
}

/* OTHER FUNCTIONS ***********************************************************************/
/*****************************************************************************************/

int mc_move(nn_vec ***r, triple *surf, triple *vac){

  if(Nsurf == 0 || Nvac == 0) return 0;
  if(Nsurf == 1) return 0;
  
  int i;
  double xsi = gsl_rng_uniform(rg);

  /**********************************INSERTION STEP*************************************************/
  if(xsi > 0.5){ //try insertion
    //    printf("Attempting insertion.\n");
    
    int Nstemp = 0;
    int Nvtemp = 0;

    unsigned long int m = gsl_rng_uniform_int(rg, (unsigned long int)Nvac); //randomly select vacancy
    nbs_vec *nbs = (nbs_vec *) calloc((size_t)NN, sizeof(nbs_vec)); //12 nearest neighbors

    neighbors(r, nbs, vac[m].x, vac[m].y, vac[m].z);

    int dM=0; //number of bonds that would be created (lowers energy)
    for(i=0; i<NN; i++){
      if(nbs[i].occ == 1) dM++; //find number of bonds that would be created by inserting an atom
      if(nbs[i].occ == 0 && nbs[i].n == 0) Nvtemp ++; //increase in vacancy atoms if inserted
      if(nbs[i].occ == 1 && nbs[i].n == NN-1) Nstemp ++; //decrease in surface atoms if inserted
    }
    Nvtemp = Nvtemp - 1; //adding an atom in the vacancy would eliminate the vacancy
    Nstemp = Nstemp - 1; //adding a surface atom would counteract the decrease by 1
    double alpha = 1.0/((Nsurf-Nstemp)/(1.0*Nvac));
    double weight = exp(mu*beta + epsilon*dM*beta)*alpha;

    // printf("alpha: %lf\n", alpha);
    //   printf("Acceptance weight: %.15e\n", weight);

    if(weight > 1){ //accept it
      r[vac[m].x][vac[m].y][vac[m].z].occ = 1; //make occupation number = 1
      Ntot++;

      for(i=0; i<NN; i++){
	//increase n for the surrounding sites
	r[nbs[i].x][nbs[i].y][nbs[i].z].n ++;
	nbs[i].n++;

	if(nbs[i].occ == 0 && nbs[i].n == 1){
	  //These sites are now neighbors - add them to vacany list
	  add(vac, &Nvac, nbs[i].x, nbs[i].y, nbs[i].z);
	}
	if(nbs[i].occ == 1 && nbs[i].n == NN){
	  //These atoms are no longer surface atoms - delete them from surface list
	  delete(surf, &Nsurf, nbs[i].x, nbs[i].y, nbs[i].z);
	}
      }

      //finally, add the new atom to the surface list
      add(surf, &Nsurf, vac[m].x, vac[m].y, vac[m].z);
      //and delete the new atom from the vacany list
      delete(vac, &Nvac, vac[m].x, vac[m].y, vac[m].z);

      temp_mc++;
      free(nbs);
      return 1;

    }else{
      double xsi2 = gsl_rng_uniform(rg);
      if(weight > xsi2){
	//do the same stuff as if weight > 1
	r[vac[m].x][vac[m].y][vac[m].z].occ = 1; //make occupation number = 1
  Ntot++;

	for(i=0; i<NN; i++){
	  //increase n for the surrounding sites
	  r[nbs[i].x][nbs[i].y][nbs[i].z].n ++;
	  nbs[i].n++;
	  
	  if(nbs[i].occ == 0 && nbs[i].n == 1){
	    //These sites are now neighbors - add them to vacany list
	    add(vac, &Nvac, nbs[i].x, nbs[i].y, nbs[i].z);
	  }
	  if(nbs[i].occ == 1 && nbs[i].n == NN){
	    //These atoms are no longer surface atoms - delete them from surface list
	    delete(surf, &Nsurf, nbs[i].x, nbs[i].y, nbs[i].z);
	  }
	}
	
	//finally, add the new atom to the surface list
	add(surf, &Nsurf, vac[m].x, vac[m].y, vac[m].z);
	//and delete the new atom from the vacany list
	delete(vac, &Nvac, vac[m].x, vac[m].y, vac[m].z);
	
	temp_mc++;
	free(nbs);
	return 1;

      }else{
	      free(nbs);
	      return 0;
      }
    }    
  }

  /*******************************************DELETION STEP******************************************/
  else{ //try deletion
    //   printf("Attempting deletion.\n");

    if(Nsurf==1) return 0;

    int Nstemp = 0;
    int Nvtemp = 0;

    long unsigned int m = gsl_rng_uniform_int(rg, (long unsigned)Nsurf); //randomly select surface atom
    nbs_vec *nbs = (nbs_vec *) calloc((size_t)NN, sizeof(nbs_vec)); //12 nearest neighbors

    neighbors(r, nbs, surf[m].x, surf[m].y, surf[m].z);
  
    int dM=0; //number of bonds that would be destroyed (increases energy)
    for(i=0; i<NN; i++){
      if(nbs[i].occ == 1) dM++; //find number of bonds that would be destroyed by deleting an atom
      if(nbs[i].occ == 0 && nbs[i].n == 1) Nvtemp ++; //decrease in vacancy sites if deleted
      if(nbs[i].occ == 1 && nbs[i].n == NN) Nstemp ++; //increase in surface atoms if inserted
    }
    Nvtemp = Nvtemp - 1; //deleting the atom would counteract the vacancy decrease by 1
    Nstemp = Nstemp - 1; //deleting the atom would decrease the number of surface atoms by 1
    double alpha = 1.0/((Nvac-Nvtemp)/(1.0*Nsurf));
    double weight = exp(-mu*beta - epsilon*dM*beta)*alpha;

    //    printf("alpha: %.15e\n", alpha);
    //   printf("Acceptance weight: %.15e\n", weight);

   if(weight > 1){ //accept it
      r[surf[m].x][surf[m].y][surf[m].z].occ = 0; //make occupation number = 0
      Ntot--;

      for(i=0; i<NN; i++){
	nbs[i].n --;
	r[nbs[i].x][nbs[i].y][nbs[i].z].n -- ;

	if(nbs[i].occ == 0 && nbs[i].n == 0){
	  //These sites are no longer surface vacancies - delete them from vacancy list
	  delete(vac, &Nvac, nbs[i].x, nbs[i].y, nbs[i].z);
	}
	if(nbs[i].occ == 1 && nbs[i].n == NN-1){
	  //These atoms are now surface atoms - add them to surface list
	  add(surf, &Nsurf, nbs[i].x, nbs[i].y, nbs[i].z);
	}
	
	if(r[nbs[i].x][nbs[i].y][nbs[i].z].n < 0) printf("Warning - nn less than 0!\n");
      }
      //finally, add the new vacancy to the vacany list
      if(r[surf[m].x][surf[m].y][surf[m].z].n > 0) add(vac, &Nvac, surf[m].x, surf[m].y, surf[m].z);
      //and delete the atom from the surface list
      delete(surf, &Nsurf, surf[m].x, surf[m].y, surf[m].z);

      //Now, check for lone atoms and delete them and their vacancies
      for(i=0; i<NN; i++){
	if(nbs[i].occ == 1 && nbs[i].n == 0){
	  //These atoms are now detached from the nanoparticle - delete them from surface list
	  r[nbs[i].x][nbs[i].y][nbs[i].z].occ = 0;
    Ntot--;
	  delete(surf, &Nsurf, nbs[i].x, nbs[i].y, nbs[i].z);
	    //Delete vacancies associated with the atom
	    nbs_vec *nbs_temp = (nbs_vec *) calloc((size_t)NN, sizeof(nbs_vec));
	    neighbors(r, nbs_temp, nbs[i].x, nbs[i].y, nbs[i].z);
	    int q;
	    for(q=0; q<NN; q++){
	      nbs_temp[q].n--;
	      r[nbs_temp[q].x][nbs_temp[q].y][nbs_temp[q].z].n--;
	      if(nbs_temp[q].n==0) delete(vac, &Nvac, nbs_temp[q].x, nbs_temp[q].y, nbs_temp[q].z);
	    }
	}
      }
	free(nbs);
      return 2;
   }else{
     double xsi2 = gsl_rng_uniform(rg);
      if(weight > xsi2){
	//do the same stuff as if weight > 1
	r[surf[m].x][surf[m].y][surf[m].z].occ = 0; //make occupation number = 0
  Ntot--;

	for(i=0; i<NN; i++){
	  nbs[i].n --;
	  r[nbs[i].x][nbs[i].y][nbs[i].z].n -- ;

	  if(nbs[i].occ == 0 && nbs[i].n == 0){
	    //These sites are no longer surface vacancies - delete them from vacancy list
	    delete(vac, &Nvac, nbs[i].x, nbs[i].y, nbs[i].z);
	  }
	  if(nbs[i].occ == 1 && nbs[i].n == NN-1){
	    //These atoms are now surface atoms - add them to surface list
	    add(surf, &Nsurf, nbs[i].x, nbs[i].y, nbs[i].z);
	  }
	
	  if(r[nbs[i].x][nbs[i].y][nbs[i].z].n < 0) printf("Warning - nn less than 0!\n");
	}
	//finally, add the new vacancy to the vacany list
	if(r[surf[m].x][surf[m].y][surf[m].z].n > 0) add(vac, &Nvac, surf[m].x, surf[m].y, surf[m].z);
	//and delete the atom from the surface list
	delete(surf, &Nsurf, surf[m].x, surf[m].y, surf[m].z);

	//Now, check for lone atoms and delete them and their vacancies
	for(i=0; i<NN; i++){
	  if(nbs[i].occ == 1 && nbs[i].n == 0){
	  //These atoms are now detached from the nanoparticle - delete them from surface list
	  r[nbs[i].x][nbs[i].y][nbs[i].z].occ = 0;
    Ntot--;
	  delete(surf, &Nsurf, nbs[i].x, nbs[i].y, nbs[i].z);
	    //Delete vacancies associated with the atom
	    nbs_vec *nbs_temp = (nbs_vec *) calloc((size_t)NN, sizeof(nbs_vec));
	    neighbors(r, nbs_temp, nbs[i].x, nbs[i].y, nbs[i].z);
	    int q;
	    for(q=0; q<NN; q++){
	      nbs_temp[q].n--;
	      r[nbs_temp[q].x][nbs_temp[q].y][nbs_temp[q].z].n--;
	      if(nbs_temp[q].n==0) delete(vac, &Nvac, nbs_temp[q].x, nbs_temp[q].y, nbs_temp[q].z);
	    }
	  }
	}
	free(nbs);
	return 2;
      }else{
	      free(nbs);
	      return 0; 
      }   
   }
  } 
}



void neighbors(nn_vec ***r, nbs_vec *nbs, int i, int j, int k){
  if(i==0 || j==0 || k==0 || i==(nx-1) || j==(ny-1) || k==(nz-1)){
    printf("Error: lattice has grown to edge of bounding cube. Exiting program.\n");
    exit(1);
  }    
  if(i>0 && j>0){
    nbs[0].x = i-1;
    nbs[0].y = j-1;
    nbs[0].z = k;
    nbs[0].n = r[i-1][j-1][k].n;
    nbs[0].occ = r[i-1][j-1][k].occ;
  }if(i>0 && j<(ny-1)){
    nbs[1].x = i-1;
    nbs[1].y = j+1;
    nbs[1].z = k;
    nbs[1].n = r[i-1][j+1][k].n;
    nbs[1].occ = r[i-1][j+1][k].occ;
  }if(i<(nx-1) && j>0){
    nbs[2].x = i+1;
    nbs[2].y = j-1;
    nbs[2].z = k;
    nbs[2].n = r[i+1][j-1][k].n;
    nbs[2].occ = r[i+1][j-1][k].occ;
  }if(i<(nx-1) && j<(ny-1)){
    nbs[3].x = i+1;
    nbs[3].y = j+1;
    nbs[3].z = k;
    nbs[3].n = r[i+1][j+1][k].n;
    nbs[3].occ = r[i+1][j+1][k].occ;
  }if(j>0 && k<(nz-1)){
    nbs[4].x = i;
    nbs[4].y = j-1;
    nbs[4].z = k+1;
    nbs[4].n = r[i][j-1][k+1].n;
    nbs[4].occ = r[i][j-1][k+1].occ;
  }if(i>0 && k<(nz-1)){
    nbs[5].x = i-1;
    nbs[5].y = j;
    nbs[5].z = k+1;
    nbs[5].n = r[i-1][j][k+1].n;
    nbs[5].occ = r[i-1][j][k+1].occ;
  }if(i<(nx-1) && k<(nz-1)){
    nbs[6].x = i+1;
    nbs[6].y = j;
    nbs[6].z = k+1;
    nbs[6].n = r[i+1][j][k+1].n;
    nbs[6].occ = r[i+1][j][k+1].occ;
  }if(j<(ny-1) && k<(nz-1)){
    nbs[7].x = i;
    nbs[7].y = j+1;
    nbs[7].z = k+1;
    nbs[7].n = r[i][j+1][k+1].n;
    nbs[7].occ = r[i][j+1][k+1].occ;
  }if(j>0 && k>0){
    nbs[8].x = i;
    nbs[8].y = j-1;
    nbs[8].z = k-1;
    nbs[8].n = r[i][j-1][k-1].n;
    nbs[8].occ = r[i][j-1][k-1].occ;
  }if(i>0 && k>0){
    nbs[9].x = i-1;
    nbs[9].y = j;
    nbs[9].z = k-1;
    nbs[9].n = r[i-1][j][k-1].n;
    nbs[9].occ = r[i-1][j][k-1].occ;
  }if(i<(nx-1) && k>0){
    nbs[10].x = i+1;
    nbs[10].y = j;
    nbs[10].z = k-1;
    nbs[10].n = r[i+1][j][k-1].n;
    nbs[10].occ = r[i+1][j][k-1].occ;
   }if(j<(ny-1) && k>0){
    nbs[11].x = i;
    nbs[11].y = j+1;
    nbs[11].z = k-1;
    nbs[11].n = r[i][j+1][k-1].n;
    nbs[11].occ = r[i][j+1][k-1].occ;
   }
}


void print_traj(nn_vec ****t){

  long length = npass/config_freq+1;
  int i,j,k;
  long l;
  int N = nx*ny*nz/2;

  FILE *f;
  f = fopen("/home/layne/research/nanoparticle-project/traj.xyz", "w");

  for(l=0; l<length; l++){
    fprintf(f,"%d\n", N); //this time, since we aren't using them again, just record the actual atoms
    fprintf(f, "Step %ld\n", l+1);

    for(i = 0; i < nx; i++) {
      for(j = 0; j < ny; j++){
	for(k = 0; k < nz; k++){
	  if((i+j+k)%2 == 0){
	    if(t[l][i][j][k].occ == 1 && t[l][i][j][k].n < NN)  fprintf(f, "Au\t%d\t%d\t%d\n", i-nx/2, j-ny/2, k-nz/2);
	    //	    if(i==0 || i==2*nl-1 || j==0 || j==2*nl-1 || k==0 || k==2*nl-1) fprintf(f, "H\t%d\t%d\t%d\n", i-nl, j-nl, k-nl);
	    else fprintf(f, "Au\t%d\t%d\t%d\n", 0, 0, 0);
	  }
	}
      }
    }
  }

  fclose(f);
  f = NULL;

}

void print_surf(triple *surf, nn_vec ***r, int t, int N) { 
   int i,j,k;
   FILE *f;
   char filename[80];
   sprintf(filename,"surf_%04d.xyz",t);
   f = fopen(filename, "w");
 
   fprintf(f, "%d\n",N); 
   fprintf(f, "shape\n");
   int q=0;
   while(q<N){
     i = surf[q].x;
     j = surf[q].y;
     k = surf[q].z;
     if(r[i][j][k].occ == 1){
      fprintf(f, "Au\t%d\t%d\t%d\t%d\t%d\n", (i-nx), (j-ny), (k-nz), r[i][j][k].n, r[i][j][k].occ);
     }
     q++;
   }
   fclose(f);
   f = NULL;    
}

void print_conf(nn_vec ***r, int t, int N) { 
   int i,j,k;
   FILE *f;
   char filename[80];
   sprintf(filename,"conf_%04d.xyz",t);
   f = fopen(filename, "w");
 
   fprintf(f, "%d\n",N); 
   fprintf(f, "shape\n");
   for(i = 0; i < nx; i++) {
    for(j = 0; j < ny; j++){
      for(k = 0; k < nz; k++){
        if(r[i][j][k].occ == 1){
          fprintf(f, "Au\t%d\t%d\t%d\t%d\t%d\n", (i-nx/2), (j-ny/2), (k-nz/2), r[i][j][k].n, r[i][j][k].occ);
        }
      }
    }
   }
   fclose(f);
   f = NULL;    
}


void print_final(nn_vec ***r) {

  int i,j,k;
  int N = nx*ny*nz/2;

  FILE *f;
  f = fopen("/home/layne/research/nanoparticle-project/final_config.xyz", "w");

  fprintf(f, "%d\n",2*N); //vmd will think there are twice as many atoms as there actually are, due to unused array elements
  fprintf(f, "shape\n");
  for(i = 0; i < nx; i++) {
    for(j = 0; j < ny; j++){
      for(k = 0; k < nz; k++){
	if((i+j+k)%2 != 0){
	  fprintf(f, "H\t%d\t%d\t%d\t%d\t%d\n", (i-nx/2), (j-ny/2), (k-nz/2), r[i][j][k].n, r[i][j][k].occ);
	}
       	else{
	  if(r[i][j][k].occ == 0) fprintf(f, "H\t%d\t%d\t%d\t%d\t%d\n", (i-nx/2), (j-ny/2), (k-nz/2), r[i][j][k].n, r[i][j][k].occ);
	  else  fprintf(f, "Au\t%d\t%d\t%d\t%d\t%d\n", (i-nx/2), (j-ny/2), (k-nz/2), r[i][j][k].n, r[i][j][k].occ);
	}
      }
    }
  }
  fclose(f);
  f = NULL;    
}

/* Methods for adding or deleting elements from a (dynamically allocated) array *********************************************/

void add(triple *list, int *length, int x, int y, int z){

  //it would be good to check that list isn't just pointing into the nether
  list[*length].x = x;
  list[*length].y = y;
  list[*length].z = z;
  *length = *length + 1;
}

void delete(triple *list, int *length, int x, int y, int z){

  int i,m=0; //setting m=0 might result in bad behavior, but I'm lazy
  for(i = 0; i < *length; i++){
    if((list[i].x == x) && (list[i].y == y) && (list[i].z == z)){
      m = i;
      break;
    }
  }
  for(i = m; i < *length; i++){
    list[i].x = list[i+1].x;
    list[i].y = list[i+1].y;
    list[i].z = list[i+1].z;
  }
  *length = *length - 1;
}

/*********************************************************************
Test functions
*********************************************************************/

/*
int test_one_insert(nn_vec ***r, triple *surf, triple *vac){

	int temp = 0;
	while(temp != 1){
//		printf("Test. Attempting MC move.\n");
	   	temp = mc_move_test(r, surf, vac);
//		printf("%d\n", temp);
	}

	return temp;
}

int test_two_delete(nn_vec ***r, triple *surf, triple *vac){

	int temp = 0;
	while(temp != 2){
		temp = mc_move(r, surf, vac);
	}

	return temp;
}

int test_two_insert(nn_vec ***r, triple *surf, triple *vac){

  int temp =0;
  while(temp != 1) temp = mc_move_test(r, surf,vac);
  return temp;
}

int test_loneatom(nn_vec ***r, triple *surf, triple *vac){

  int i = 0;
  int m=0; //this might result in bad behavior
  for(i=0; i<4; i++){
    if(surf[i].x == nl && surf[i].y == nl && surf[i].z == nl) m = i;
  }

  nbs_vec *nbs = (nbs_vec *) calloc((size_t)NN, sizeof(nbs_vec)); //12 nearest neighbors
  neighbors(r, nbs, surf[m].x, surf[m].y, surf[m].z);

  r[surf[m].x][surf[m].y][surf[m].z].occ = 0; //make occupation number = 0

  for(i=0; i<NN; i++){
    nbs[i].n --;
    r[nbs[i].x][nbs[i].y][nbs[i].z].n -- ;

    if(nbs[i].occ == 0 && nbs[i].n == 0){
      //These sites are no longer surface vacancies - delete them from vacancy list
      delete(vac, &Nvac, nbs[i].x, nbs[i].y, nbs[i].z);
    }

    if(nbs[i].occ == 1 && nbs[i].n == 11){
      //These atoms are now surface atoms - add them to surface list
      add(surf, &Nsurf, nbs[i].x, nbs[i].y, nbs[i].z);
    }
	
    if(r[nbs[i].x][nbs[i].y][nbs[i].z].n < 0) printf("Warning - nn less than 0!\n");
  }

  //finally, add the new vacancy to the vacany list
  if(r[surf[m].x][surf[m].y][surf[m].z].n > 0) add(vac, &Nvac, surf[m].x, surf[m].y, surf[m].z);
  //and delete the atom from the surface list
  delete(surf, &Nsurf, surf[m].x, surf[m].y, surf[m].z);

  //Now, check for lone atoms and delete them and their vacancies
  for(i=0; i<NN; i++){
  if(nbs[i].occ == 1 && nbs[i].n == 0){
	  //These atoms are now detached from the nanoparticle - delete them from surface list
	  r[nbs[i].x][nbs[i].y][nbs[i].z].occ = 0;
	  delete(surf, &Nsurf, nbs[i].x, nbs[i].y, nbs[i].z);
	    //Delete vacancies associated with the atom
	    nbs_vec *nbs_temp = (nbs_vec *) calloc((size_t)NN, sizeof(nbs_vec));
	    neighbors(r, nbs_temp, nbs[i].x, nbs[i].y, nbs[i].z);
	    int q;
	    for(q=0; q<NN; q++){
	      nbs_temp[q].n--;
	      r[nbs_temp[q].x][nbs_temp[q].y][nbs_temp[q].z].n--;
	      if(nbs_temp[q].n==0) delete(vac, &Nvac, nbs_temp[q].x, nbs_temp[q].y, nbs_temp[q].z);
	    }
	}
  }

  return 2;
}


int mc_move_test(nn_vec ***r, triple *surf, triple *vac){

  if(Nsurf == 0 || Nvac == 0) return 0;
  
  int i;
  double xsi = gsl_rng_uniform(rg);

 // **********************************INSERTION STEP*************************************************
  if(xsi > 0.5){ //try insertion
    //    printf("Attempting insertion.\n");
    
    int Nstemp = 0;
    int Nvtemp = 0;

    long unsigned int m = gsl_rng_uniform_int(rg, (long unsigned int)Nvac); //randomly select vacancy
    nbs_vec *nbs = (nbs_vec *) calloc((size_t)12, sizeof(nbs_vec)); //12 nearest neighbors

    neighbors(r, nbs, vac[m].x, vac[m].y, vac[m].z);

    int dM=0; //number of bonds that would be created (lowers energy)
    for(i=0; i<NN; i++){
      if(nbs[i].occ == 1) dM++; //find number of bonds that would be created by inserting an atom
      if(nbs[i].occ == 0 && nbs[i].n == 0) Nvtemp ++; //increase in vacancy atoms if inserted
      if(nbs[i].occ == 1 && nbs[i].n == 11) Nstemp ++; //decrease in surface atoms if inserted
    }
    // ********************This is the testing line*********************************
    if(dM!=3) return 0;

    // ****************************************************************************

    Nvtemp = Nvtemp - 1; //adding an atom in the vacancy would eliminate the vacancy
    Nstemp = Nstemp - 1; //adding a surface atom would counteract the decrease by 1
    double alpha = 1.0/((Nsurf-Nstemp)/(1.0*Nvac));
    double weight = exp(mu*beta + epsilon*dM*beta)*alpha;

    printf("alpha: %lf\n", alpha);
    printf("Acceptance weight: %.15e\n", weight);

    if(weight > 1){ //accept it
      r[vac[m].x][vac[m].y][vac[m].z].occ = 1; //make occupation number = 1

      for(i=0; i<NN; i++){
	//increase n for the surrounding sites
	r[nbs[i].x][nbs[i].y][nbs[i].z].n ++;
	nbs[i].n++;

	if(nbs[i].occ == 0 && nbs[i].n == 1){
	  //These sites are now neighbors - add them to vacany list
	  add(vac, &Nvac, nbs[i].x, nbs[i].y, nbs[i].z);
	}
	if(nbs[i].occ == 1 && nbs[i].n == 12){
	  //These atoms are no longer surface atoms - delete them from surface list
	  delete(surf, &Nsurf, nbs[i].x, nbs[i].y, nbs[i].z);
	}
      }

      //finally, add the new atom to the surface list
      add(surf, &Nsurf, vac[m].x, vac[m].y, vac[m].z);
      //and delete the new atom from the vacany list
      delete(vac, &Nvac, vac[m].x, vac[m].y, vac[m].z);

      temp_mc++;
      return 1;

    }else{
      double xsi2 = gsl_rng_uniform(rg);
      if(weight > xsi2){
	//do the same stuff as if weight > 1
	r[vac[m].x][vac[m].y][vac[m].z].occ = 1; //make occupation number = 1

	for(i=0; i<NN; i++){
	  //increase n for the surrounding sites
	  r[nbs[i].x][nbs[i].y][nbs[i].z].n ++;
	  nbs[i].n++;
	  
	  if(nbs[i].occ == 0 && nbs[i].n == 1){
	    //These sites are now neighbors - add them to vacany list
	    add(vac, &Nvac, nbs[i].x, nbs[i].y, nbs[i].z);
	  }
	  if(nbs[i].occ == 1 && nbs[i].n == 12){
	    //These atoms are no longer surface atoms - delete them from surface list
	    delete(surf, &Nsurf, nbs[i].x, nbs[i].y, nbs[i].z);
	  }
	}
	
	//finally, add the new atom to the surface list
	add(surf, &Nsurf, vac[m].x, vac[m].y, vac[m].z);
	//and delete the new atom from the vacany list
	delete(vac, &Nvac, vac[m].x, vac[m].y, vac[m].z);
	
	temp_mc++;
	return 1;

      }else return 0;
    }       
  }


}

*/
