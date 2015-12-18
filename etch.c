/* "etch.c" */
/* This program takes in a nanoparticle configuration and applies MCMC to it
   in order to "etch" it. 
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_rng.h>
#include "etch.h"

//Use GSL to generate random numbers
gsl_rng *rg;

//Control variables
#define T 1.0
#define mu -10.0
#define epsilon 1.0
#define npass 9000
#define d 1e-6 //small factor used when converting to int
#define config_freq 10

int nside; //to be read in
int Nsurf = 0, Nvac = 0; //number of surface and vacany atoms
int Nsurf0 = 0, Nvac0 = 0; //To keep track of initial numbers;
int Ntot = 0;

int temp_mc=0;

int main(int argc, char *argv[]){

  int Ngrid;
  int i,j,k,l;
  int time = 0;
  int mc_total = 0; //number of total attempts
  int mc_accept = 0;
  int mc_insert = 0; //number of insertion moves
  int mc_delete = 0; //number of deletion moves
  char ignore[1024]; //for skipping lines

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

  vec3 v_curr;
  int n_curr;
  int occ_curr;
  FILE *in;
  in = fopen("/home/layne/research/nanoparticles/init_config.xyz", "r");
  if(in == NULL) printf("problem!\n");
  fscanf(in, "%d\n", &Ngrid);

  nside = (int) (pow(Ngrid, 1.0/3.0) + 1e-5); //number of unit cells per side
  nside = nside/2;

  //allocate memory (for lattice)
  //The lattice is 3d, each side indexes from 0 to 2*nside-1
  nn_vec ***lat;
  lat = (nn_vec ***)calloc((size_t)2*nside, sizeof(nn_vec**));
  for(i = 0; i < 2*nside; i++){
    lat[i] = (nn_vec **) calloc((size_t)2*nside, sizeof(nn_vec*));
    for(j = 0; j < 2*nside; j++){
      lat[i][j] = (nn_vec*) calloc((size_t)2*nside, sizeof(nn_vec));
    }
  }

  fgets(ignore, sizeof(ignore), in); //skip line

  char atom_curr[2];
  i = 0;
  j = 0;
  k = 0;

  //Now actually read in data
  for(i = 0; i < 2*nside; i++) {
    for(j = 0; j < 2*nside; j++){
      for(k = 0; k < 2*nside; k++){
	fscanf(in, "%s\t" "%lf\t" "%lf\t" "%lf\t" "%d\t"  "%d\n", atom_curr, &lat[i][j][k].v.x, &lat[i][j][k].v.y, &lat[i][j][k].v.z, &lat[i][j][k].n, &lat[i][j][k].occ);
	if(lat[i][j][k].occ==1){
	  Ntot++;
	  if(lat[i][j][k].n < 12) Nsurf++;
	}else if(lat[i][j][k].occ==0){
	  if(lat[i][j][k].n > 0) Nvac++;
	}
      }
    }
  }
  fclose(in);
  in = NULL;

  //Create a 'trajectory' array that stores configurations from mc steps
  nn_vec ****traj;
  printf("npass/config_freq + 2: %d\n", npass/config_freq+2);
  traj = (nn_vec ****)calloc((size_t)(npass/config_freq + 2), sizeof(nn_vec***));
  for(l = 0; l < (npass/config_freq+2); l++){
    traj[l] = (nn_vec ***) calloc((size_t)2*nside, sizeof(nn_vec**));
    for(i = 0; i < 2*nside; i++){
      traj[l][i] = (nn_vec **) calloc((size_t)2*nside, sizeof(nn_vec*));
      for(j = 0; j < 2*nside; j++){
	traj[l][i][j] = (nn_vec *) calloc((size_t)2*nside, sizeof(nn_vec));
      }
    }
  }

  for(i = 0; i < 2*nside; i++){
    for(j = 0; j < 2*nside; j++){
      for(k = 0; k < 2*nside; k++){
	traj[time][i][j][k].v.x = lat[i][j][k].v.x;
	traj[time][i][j][k].v.y = lat[i][j][k].v.y;
	traj[time][i][j][k].v.z = lat[i][j][k].v.z;
	traj[time][i][j][k].n = lat[i][j][k].n;
	traj[time][i][j][k].occ = lat[i][j][k].occ;
      }
    }
  }
  time = time + 1;

  printf("Total number of atoms at start: %d\n", Ntot);
  printf("number of surface atoms: %d\n", Nsurf);
  printf("number of vacancies: %d\n", Nvac);

  Nsurf0 = Nsurf;
  Nvac0 = Nvac;

  /* Fill arrays with the indices of the locations of surface and vacany atoms *********************************************************************************/

  triple *surf = (triple *) calloc((size_t)2*Nsurf, sizeof(triple)); //make these twice the size they need to be to accomodate the possibility
  triple *vac = (triple *) calloc((size_t)2*Nvac, sizeof(triple)); //of getting more than you started with

  int ls=0; 
  int lv=0;
  for(i = 0; i < 2*nside; i++) {
    for(j = 0; j < 2*nside; j++){
      for(k = 0; k < 2*nside; k++){
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

  //Create an array that stores the number of surface atoms at every step

  /* Now do Monte Carlo ****************************************************/
  //  test_one_insert(lat, surf, vac);
  //  printf("number of surface atoms: %d\n", Nsurf);
  //  printf("number of vacancies: %d\n", Nvac);
  //  test_two_delete(lat, surf, vac);
  //  test_two_insert(lat, surf, vac);

  
  for(l = 0; l < npass; l++){
    printf("Pass number %d\n", l+1);
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
        for(i = 0; i < 2*nside; i++){
	  for(j = 0; j < 2*nside; j++){
	    for(k = 0; k < 2*nside; k++){
	      traj[time][i][j][k].v.x = lat[i][j][k].v.x;
	      traj[time][i][j][k].v.y = lat[i][j][k].v.y;
	      traj[time][i][j][k].v.z = lat[i][j][k].v.z;
	      traj[time][i][j][k].n = lat[i][j][k].n;
	      traj[time][i][j][k].occ = lat[i][j][k].occ;
	    }
	  }
	  } 
	//	print_config_t(lat, time);
	time ++;
    }
    mc_total ++;
  }
 

  print_traj(traj);
  printf("Total number of atoms at start: %d\n", Ntot);
  printf("final number of surface atoms: %d\n", Nsurf);
  printf("final number of vacancies: %d\n", Nvac);
  printf("Attempted moves: %d\n", mc_total);
  printf("Accepted moves: %d\n", mc_delete + mc_insert);
  printf("Deletion moves: %d\n", mc_delete);
  printf("Insertion moves: %d\n", mc_insert);
  printf("Accepted ratio: %lf\n", (mc_delete + mc_insert)/(1.0*mc_total));

  return 1; //end of program
}

/* OTHER FUNCTIONS ***********************************************************************/
/*****************************************************************************************/

int test_one_insert(nn_vec ***r, triple *surf, triple *vac){

	int temp = 0;
	while(temp != 1){
//		printf("Test. Attempting MC move.\n");
	   	temp = mc_move(r, surf, vac);
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

int mc_move(nn_vec ***r, triple *surf, triple *vac){

  if(Nsurf == 0 || Nvac == 0) return 0;
  if(Nsurf == 1) return 0;
  
  int i,j;
  double xsi = gsl_rng_uniform(rg);

  /**********************************INSERTION STEP*************************************************/
  if(xsi > 0.5){ //try insertion
    //    printf("Attempting insertion.\n");
    
    int Nstemp = 0;
    int Nvtemp = 0;

    int m = gsl_rng_uniform_int(rg, Nvac); //randomly select vacancy
    nn_vec *nbs = (nn_vec *) calloc((size_t)12, sizeof(nn_vec)); //12 nearest neighbors
    nn_vec *nbs_temp = (nn_vec *) calloc((size_t)12, sizeof(nn_vec)); //for debugging

    neighbors(r, nbs, vac[m].x, vac[m].y, vac[m].z);

    int dM=0; //number of bonds that would be created (lowers energy)
    for(i=0; i<NN; i++){
      if(nbs[i].occ == 1) dM++; //find number of bonds that would be created by inserting an atom
      if(nbs[i].occ == 0 && nbs[i].n == 0) Nvtemp ++; //increase in vacancy atoms if inserted
      if(nbs[i].occ == 1 && nbs[i].n == 11) Nstemp ++; //decrease in surface atoms if inserted
    }
    Nvtemp = Nvtemp - 1; //adding an atom in the vacancy would eliminate the vacancy
    Nstemp = Nstemp - 1; //adding a surface atom would counteract the decrease by 1
    double alpha = 1.0/((Nsurf-Nstemp)/(1.0*Nvac));
    double weight = exp(mu/T + epsilon*dM/T)*alpha;

    //    printf("alpha: %lf\n", alpha);
    //   printf("Acceptance weight: %.15e\n", weight);

    if(weight > 1){ //accept it
      r[(int)(vac[m].x+d)][(int)(vac[m].y+d)][(int)(vac[m].z+d)].occ = 1; //make occupation number = 1

      for(i=0; i<NN; i++){
	//increase n for the surrounding sites
	r[(int)(nbs[i].v.x+d)][(int)(nbs[i].v.y+d)][(int)(nbs[i].v.z+d)].n ++;
	nbs[i].n++;

	if(nbs[i].occ == 0 && nbs[i].n == 1){
	  //These sites are now neighbors - add them to vacany list
	  add(vac, &Nvac, nbs[i].v.x, nbs[i].v.y, nbs[i].v.z);
	}
	if(nbs[i].occ == 1 && nbs[i].n == 12){
	  //These atoms are no longer surface atoms - delete them from surface list
	  delete(surf, &Nsurf, nbs[i].v.x, nbs[i].v.y, nbs[i].v.z);
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
	r[(int)(vac[m].x+d)][(int)(vac[m].y+d)][(int)(vac[m].z+d)].occ = 1; //make occupation number = 1

	for(i=0; i<NN; i++){
	  //increase n for the surrounding sites
	  r[(int)(nbs[i].v.x+d)][(int)(nbs[i].v.y+d)][(int)(nbs[i].v.z+d)].n ++;
	  nbs[i].n++;
	  
	  if(nbs[i].occ == 0 && nbs[i].n == 1){
	    //These sites are now neighbors - add them to vacany list
	    add(vac, &Nvac, nbs[i].v.x, nbs[i].v.y, nbs[i].v.z);
	  }
	  if(nbs[i].occ == 1 && nbs[i].n == 12){
	    //These atoms are no longer surface atoms - delete them from surface list
	    delete(surf, &Nsurf, nbs[i].v.x, nbs[i].v.y, nbs[i].v.z);
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

  /*******************************************DELETION STEP******************************************/
  else{ //try deletion
    //   printf("Attempting deletion.\n");

    if(Nsurf==1) return 0;

    int Nstemp = 0;
    int Nvtemp = 0;

    int m = gsl_rng_uniform_int(rg, Nsurf); //randomly select surface atom
    nn_vec *nbs = (nn_vec *) calloc((size_t)NN, sizeof(nn_vec)); //12 nearest neighbors

    neighbors(r, nbs, surf[m].x, surf[m].y, surf[m].z);
  
    int dM=0; //number of bonds that would be destroyed (increases energy)
    for(i=0; i<NN; i++){
      if(nbs[i].occ == 1) dM++; //find number of bonds that would be destroyed by deleting an atom
      if(nbs[i].occ == 0 && nbs[i].n == 1) Nvtemp ++; //decrease in vacancy sites if deleted
      if(nbs[i].occ == 1 && nbs[i].n == 12) Nstemp ++; //increase in surface atoms if inserted
    }
    Nvtemp = Nvtemp - 1; //deleting the atom would counteract the vacancy decrease by 1
    Nstemp = Nstemp - 1; //deleting the atom would decrease the number of surface atoms by 1
    double alpha = 1.0/((Nvac-Nvtemp)/(1.0*Nsurf));
    double weight = exp(-mu/T - epsilon*dM/T)*alpha;

    //    printf("Acceptance weight: %.15e\n", weight);

   if(weight > 1){ //accept it
      r[(int)(surf[m].x+d)][(int)(surf[m].y+d)][(int)(surf[m].z+d)].occ = 0; //make occupation number = 0

      for(i=0; i<NN; i++){
	nbs[i].n --;
	r[(int)(nbs[i].v.x+d)][(int)(nbs[i].v.y+d)][(int)(nbs[i].v.z+d)].n -- ;

	if(nbs[i].occ == 0 && nbs[i].n == 0){
	  //These sites are no longer surface vacancies - delete them from vacancy list
	  delete(vac, &Nvac, nbs[i].v.x, nbs[i].v.y, nbs[i].v.z);
	}
	/*	if(nbs[i].occ == 1 && nbs[i].n == 0){
	  //These atoms are now detached from the nanoparticle - delete them from surface list
	  r[(int)(nbs[i].v.x+d)][(int)(nbs[i].v.y+d)][(int)(nbs[i].v.z+d)].occ = 0;
	  delete(surf, &Nsurf, nbs[i].v.x, nbs[i].v.y, nbs[i].v.z);
	    //Delete vacancies associated with the atom
	    nn_vec *nbs_temp = (nn_vec *) calloc((size_t)NN, sizeof(nn_vec));
	    neighbors(r, nbs_temp, nbs[i].v.x, nbs[i].v.y, nbs[i].v.z);
	    int q;
	    for(q=0; q<NN; q++){
	      nbs_temp[q].n--;
	      r[(int)nbs_temp[q].v.x][(int)nbs_temp[q].v.y][(int)nbs_temp[q].v.z].n--;
	      if(nbs_temp[q].n==0) delete(vac, &Nvac, nbs_temp[q].v.x, nbs_temp[q].v.y, nbs_temp[q].v.z);
	    }
	}
	*/
	if(nbs[i].occ == 1 && nbs[i].n == 11){
	  //These atoms are now surface atoms - add them to surface list
	  add(surf, &Nsurf, nbs[i].v.x, nbs[i].v.y, nbs[i].v.z);
	}
	
	if(r[(int)(nbs[i].v.x+d)][(int)(nbs[i].v.y+d)][(int)(nbs[i].v.z+d)].n < 0) printf("Warning - nn less than 0!\n");
      }
      //finally, add the new vacancy to the vacany list
      if(r[(int)(surf[m].x+d)][(int)(surf[m].y+d)][(int)(surf[m].z+d)].n > 0) add(vac, &Nvac, surf[m].x, surf[m].y, surf[m].z);
      //and delete the atom from the surface list
      delete(surf, &Nsurf, surf[m].x, surf[m].y, surf[m].z);
      return 2;
   }else{
     double xsi2 = gsl_rng_uniform(rg);
      if(weight > xsi2){
	//do the same stuff as if weight > 1
	r[(int)(surf[m].x+d)][(int)(surf[m].y+d)][(int)(surf[m].z+d)].occ = 0; //make occupation number = 0

	for(i=0; i<NN; i++){
	  nbs[i].n --;
	  r[(int)(nbs[i].v.x+d)][(int)(nbs[i].v.y+d)][(int)(nbs[i].v.z+d)].n -- ;

	  if(nbs[i].occ == 0 && nbs[i].n == 0){
	    //These sites are no longer surface vacancies - delete them from vacancy list
	    delete(vac, &Nvac, nbs[i].v.x, nbs[i].v.y, nbs[i].v.z);
	  }
	  /*	  if(nbs[i].occ == 1 && nbs[i].n == 0){
	    //These atoms are now detached from the nanoparticle - delete them from surface list
	    r[(int)nbs[i].v.x][(int)nbs[i].v.y][(int)nbs[i].v.z].occ = 0;
	    delete(surf, &Nsurf, nbs[i].v.x, nbs[i].v.y, nbs[i].v.z);
	    //Delete vacancies associated with the atom
	    nn_vec *nbs_temp = (nn_vec *) calloc((size_t)NN, sizeof(nn_vec));
	    neighbors(r, nbs_temp, nbs[i].v.x, nbs[i].v.y, nbs[i].v.z);
	    int q;
	    for(q=0; q<NN; q++){
	      nbs_temp[q].n--;
	      r[(int)nbs_temp[q].v.x][(int)nbs_temp[q].v.y][(int)nbs_temp[q].v.z].n--;
	      if(nbs_temp[q].n==0) delete(vac, &Nvac, nbs_temp[q].v.x, nbs_temp[q].v.y, nbs_temp[q].v.z);
	    }
	  }
	  */
	  if(nbs[i].occ == 1 && nbs[i].n == 11){
	    //These atoms are now surface atoms - add them to surface list
	    add(surf, &Nsurf, nbs[i].v.x, nbs[i].v.y, nbs[i].v.z);
	  }
	
	  if(r[(int)(nbs[i].v.x+d)][(int)(nbs[i].v.y+d)][(int)(nbs[i].v.z+d)].n < 0) printf("Warning - nn less than 0!\n");
	}
	//finally, add the new vacancy to the vacany list
	if(r[(int)(surf[m].x+d)][(int)(surf[m].y+d)][(int)(surf[m].z+d)].n > 0) add(vac, &Nvac, surf[m].x, surf[m].y, surf[m].z);
	//and delete the atom from the surface list
	delete(surf, &Nsurf, surf[m].x, surf[m].y, surf[m].z);
	return 2;
      }else return 0;    
   }
  } 
}

int mc_move_test(nn_vec ***r, triple *surf, triple *vac){

  if(Nsurf == 0 || Nvac == 0) return 0;
  
  int i,j;
  double xsi = gsl_rng_uniform(rg);

  /**********************************INSERTION STEP*************************************************/
  if(xsi > 0.5){ //try insertion
    //    printf("Attempting insertion.\n");
    
    int Nstemp = 0;
    int Nvtemp = 0;

    int m = gsl_rng_uniform_int(rg, Nvac); //randomly select vacancy
    nn_vec *nbs = (nn_vec *) calloc((size_t)12, sizeof(nn_vec)); //12 nearest neighbors
    nn_vec *nbs_temp = (nn_vec *) calloc((size_t)12, sizeof(nn_vec)); //for debugging

    neighbors(r, nbs, vac[m].x, vac[m].y, vac[m].z);

    int dM=0; //number of bonds that would be created (lowers energy)
    for(i=0; i<NN; i++){
      if(nbs[i].occ == 1) dM++; //find number of bonds that would be created by inserting an atom
      if(nbs[i].occ == 0 && nbs[i].n == 0) Nvtemp ++; //increase in vacancy atoms if inserted
      if(nbs[i].occ == 1 && nbs[i].n == 11) Nstemp ++; //decrease in surface atoms if inserted
    }
    //********************This is the testing line*********************************//
    if(dM!=3) return 0;

    //****************************************************************************//

    Nvtemp = Nvtemp - 1; //adding an atom in the vacancy would eliminate the vacancy
    Nstemp = Nstemp - 1; //adding a surface atom would counteract the decrease by 1
    double alpha = 1.0/((Nsurf-Nstemp)/(1.0*Nvac));
    double weight = exp(mu/T + epsilon*dM/T)*alpha;

    printf("alpha: %lf\n", alpha);
    printf("Acceptance weight: %.15e\n", weight);

    if(weight > 1){ //accept it
      r[(int)(vac[m].x+d)][(int)(vac[m].y+d)][(int)(vac[m].z+d)].occ = 1; //make occupation number = 1

      for(i=0; i<NN; i++){
	//increase n for the surrounding sites
	r[(int)(nbs[i].v.x+d)][(int)(nbs[i].v.y+d)][(int)(nbs[i].v.z+d)].n ++;
	nbs[i].n++;

	if(nbs[i].occ == 0 && nbs[i].n == 1){
	  //These sites are now neighbors - add them to vacany list
	  add(vac, &Nvac, nbs[i].v.x, nbs[i].v.y, nbs[i].v.z);
	}
	if(nbs[i].occ == 1 && nbs[i].n == 12){
	  //These atoms are no longer surface atoms - delete them from surface list
	  delete(surf, &Nsurf, nbs[i].v.x, nbs[i].v.y, nbs[i].v.z);
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
	r[(int)(vac[m].x+d)][(int)(vac[m].y+d)][(int)(vac[m].z+d)].occ = 1; //make occupation number = 1

	for(i=0; i<NN; i++){
	  //increase n for the surrounding sites
	  r[(int)(nbs[i].v.x+d)][(int)(nbs[i].v.y+d)][(int)(nbs[i].v.z+d)].n ++;
	  nbs[i].n++;
	  
	  if(nbs[i].occ == 0 && nbs[i].n == 1){
	    //These sites are now neighbors - add them to vacany list
	    add(vac, &Nvac, nbs[i].v.x, nbs[i].v.y, nbs[i].v.z);
	  }
	  if(nbs[i].occ == 1 && nbs[i].n == 12){
	    //These atoms are no longer surface atoms - delete them from surface list
	    delete(surf, &Nsurf, nbs[i].v.x, nbs[i].v.y, nbs[i].v.z);
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

void neighbors(nn_vec ***r, nn_vec *nbs, int i, int j, int k){
  if(i>0 && j>0){
    nbs[0].v.x = i-1;
    nbs[0].v.y = j-1;
    nbs[0].v.z = k;
    nbs[0].n = r[i-1][j-1][k].n;
    // if(nbs[0].n == 0) printf("invalid\n");
    nbs[0].occ = r[(int)(i-1+d)][(int)(j-1+d)][(int)(k+d)].occ;
  }if(i>0 && j<(2*nside-1)){
    nbs[1].v.x = i-1;
    nbs[1].v.y = j+1;
    nbs[1].v.z = k;
    nbs[1].n = r[(int)(i-1+d)][(int)(j+1+d)][(int)(k+d)].n;
    // if(nbs[1].n == 0) printf("invalid\n");
    nbs[1].occ = r[(int)(i-1+d)][(int)(j+1+d)][(int)(k+d)].occ;
  }if(i<(2*nside-1) && j>0){
    nbs[2].v.x = i+1;
    nbs[2].v.y = j-1;
    nbs[2].v.z = k;
    nbs[2].n = r[(int)(i+1+d)][(int)(j-1+d)][(int)(k+d)].n;
    // if(nbs[2].n == 0) printf("invalid\n");
    nbs[2].occ = r[(int)(i+1+d)][(int)(j-1+d)][(int)(k+d)].occ;
  }if(i<(2*nside-1) && j<(2*nside-1)){
    nbs[3].v.x = i+1;
    nbs[3].v.y = j+1;
    nbs[3].v.z = k;
    nbs[3].n = r[(int)(i+1+d)][(int)(j+1+d)][(int)(k+d)].n;
    // if(nbs[3].n == 0) printf("invalid\n");
    nbs[3].occ = r[(int)(i+1+d)][(int)(j+1+d)][(int)(k+d)].occ;
  }if(j>0 && k<(2*nside-1)){
    nbs[4].v.x = i;
    nbs[4].v.y = j-1;
    nbs[4].v.z = k+1;
    nbs[4].n = r[(int)(i+d)][(int)(j-1+d)][(int)(k+1+d)].n;
    // if(nbs[4].n == 0) printf("invalid\n");
    nbs[4].occ = r[(int)(i+d)][(int)(j-1+d)][(int)(k+1+d)].occ;
  }if(i>0 && k<(2*nside-1)){
    nbs[5].v.x = i-1;
    nbs[5].v.y = j;
    nbs[5].v.z = k+1;
    nbs[5].n = r[(int)(i-1+d)][(int)(j+d)][(int)(k+1+d)].n;
    // if(nbs[5].n == 0) printf("invalid\n");
    nbs[5].occ = r[(int)(i-1+d)][(int)(j+d)][(int)(k+1+d)].occ;
  }if(i<(2*nside-1) && k<(2*nside-1)){
    nbs[6].v.x = i+1;
    nbs[6].v.y = j;
    nbs[6].v.z = k+1;
    nbs[6].n = r[(int)(i+1+d)][(int)(j+d)][(int)(k+1+d)].n;
    // if(nbs[6].n == 0) printf("invalid\n");
    nbs[6].occ = r[(int)(i+1+d)][(int)(j+d)][(int)(k+1+d)].occ;
  }if(j<(2*nside-1) && k<(2*nside-1)){
    nbs[7].v.x = i;
    nbs[7].v.y = j+1;
    nbs[7].v.z = k+1;
    nbs[7].n = r[(int)(i+d)][(int)(j+1+d)][(int)(k+1+d)].n;
    //  if(nbs[7].n == 0) printf("invalid\n");
    nbs[7].occ = r[(int)(i+d)][(int)(j+1+d)][(int)(k+1+d)].occ;
  }if(j>0 && k>0){
    nbs[8].v.x = i;
    nbs[8].v.y = j-1;
    nbs[8].v.z = k-1;
    nbs[8].n = r[(int)(i+d)][(int)(j-1+d)][(int)(k-1+d)].n;
    nbs[8].occ = r[(int)(i+d)][(int)(j-1+d)][(int)(k-1+d)].occ;
  }if(i>0 && k>0){
    nbs[9].v.x = i-1;
    nbs[9].v.y = j;
    nbs[9].v.z = k-1;
    nbs[9].n = r[(int)(i-1+d)][(int)(j+d)][(int)(k-1+d)].n;
    //  if(nbs[9].n == 0) printf("invalid\n");
    nbs[9].occ = r[(int)(i-1+d)][(int)(j+d)][(int)(k-1+d)].occ;
  }if(i<(2*nside-1) && k>0){
    nbs[10].v.x = i+1;
    nbs[10].v.y = j;
    nbs[10].v.z = k-1;
    nbs[10].n = r[(int)(i+1+d)][(int)(j+d)][(int)(k-1+d)].n;
    //  if(nbs[10].n == 0) printf("invalid\n");
    nbs[10].occ = r[(int)(i+1+d)][(int)(j+d)][(int)(k-1+d)].occ;
   }if(j<(2*nside-1) && k>0){
    nbs[11].v.x = i;
    nbs[11].v.y = j+1;
    nbs[11].v.z = k-1;
    nbs[11].n = r[(int)(i+d)][(int)(j+1+d)][(int)(k-1+d)].n;
    //  if(nbs[11].n == 0) printf("invalid\n");
    nbs[11].occ = r[(int)(i+d)][(int)(j+1+d)][(int)(k-1+d)].occ;
   }
}

void print_config(nn_vec ***r, int nside) {

  int i,j,k;
  int N = 4*nside*nside*nside;

  FILE *f;
  f = fopen("/home/layne/research/nanoparticles/config.xyz", "w");

  fprintf(f, "%d\n",2*N); //vmd will think there are twice as many atoms as there actually are, due to unused array elements
  fprintf(f, "initial configuration of gold nanoparticle\n");
  for(i = 0; i < 2*nside; i++) {
    for(j = 0; j < 2*nside; j++){
      for(k = 0; k < 2*nside; k++){
	if((i+j+k)%2 != 0){
	  fprintf(f, "H\t%lf\t%lf\t%lf\t%d\t%d\n", r[i][j][k].v.x, r[i][j][k].v.y, r[i][j][k].v.z, r[i][j][k].n, r[i][j][k].occ);
	}
       	else{
	  if(r[i][j][k].occ <= 0) fprintf(f, "H\t%lf\t%lf\t%lf\t%d\t%d\n", r[i][j][k].v.x, r[i][j][k].v.y, r[i][j][k].v.z, r[i][j][k].n, r[i][j][k].occ);
	  else fprintf(f, "Au\t%lf\t%lf\t%lf\t%d\t%d\n", r[i][j][k].v.x, r[i][j][k].v.y, r[i][j][k].v.z, r[i][j][k].n, r[i][j][k].occ);
	}
      }
    }
  }
  fclose(f);
  f = NULL;    
}

void print_config_t(nn_vec ***r, int time) {

  int length = npass/config_freq+1;
  int i,j,k;
  int N = 4*nside*nside*nside;

  char buffer[200];
  FILE *f;
  sprintf(buffer, "/home/layne/research/nanoparticles/config%d.xyz", time);
  f = fopen(buffer, "w");
  fprintf(f, "%d\n",2*N); //vmd will think there are twice as many atoms as there actually are, due to unused array elements
  fprintf(f, "configuration of gold nanoparticle\n");
  for(i = 0; i < 2*nside; i++) {
    for(j = 0; j < 2*nside; j++){
      for(k = 0; k < 2*nside; k++){
	if((i+j+k)%2 != 0){
	  fprintf(f, "H\t%lf\t%lf\t%lf\t%d\t%d\n", r[i][j][k].v.x, r[i][j][k].v.y, r[i][j][k].v.z, r[i][j][k].n, r[i][j][k].occ);
	}
       	else{
	  if(r[i][j][k].occ <= 0) fprintf(f, "H\t%lf\t%lf\t%lf\t%d\t%d\n", r[i][j][k].v.x, r[i][j][k].v.y, r[i][j][k].v.z, r[i][j][k].n, r[i][j][k].occ);
	  else fprintf(f, "Au\t%lf\t%lf\t%lf\t%d\t%d\n", r[i][j][k].v.x, r[i][j][k].v.y, r[i][j][k].v.z, r[i][j][k].n, r[i][j][k].occ);
	}
      }
    }
  }
  fclose(f);
  f = NULL;    
}

void print_traj(nn_vec ****t){

  int length = npass/config_freq+1;
  int i,j,k,l;
  int N = 4*nside*nside*nside;

  FILE *f;
  f = fopen("/home/layne/research/nanoparticles/traj.xyz", "w");

  for(l=0; l<length; l++){
    fprintf(f,"%d\n", N); //this time, since we aren't using them again, just record the actual atoms
    fprintf(f, "Step %d\n", l+1);

    for(i = 0; i < 2*nside; i++) {
      for(j = 0; j < 2*nside; j++){
	for(k = 0; k < 2*nside; k++){
	  if((i+j+k)%2 == 0){
	    if(t[l][i][j][k].occ == 1 && t[l][i][j][k].n < 12)  fprintf(f, "Au\t%lf\t%lf\t%lf\t%d\t%d\n", t[l][i][j][k].v.x, t[l][i][j][k].v.y, t[l][i][j][k].v.z, t[l][i][j][k].n, t[l][i][j][k].occ);
	    else fprintf(f, "Au\t%lf\t%lf\t%lf\t%d\t%d\n", 0.0, 0.0, 0.0, t[l][i][j][k].n, t[l][i][j][k].occ);
	  }
	}
      }
    }
  }

  fclose(f);
  f = NULL;

}

/* Methods for adding or deleting elements from a (dynamically allocated) array *********************************************/

void add(triple *list, int *length, double x, double y, double z){

  //it would be good to check that list isn't just pointing into the nether
  list[*length].x = (int)x;
  list[*length].y = (int)y;
  list[*length].z = (int)z;
  *length = *length + 1;
}

void delete(triple *list, int *length, double x, double y, double z){

  int i,m;
  for(i = 0; i < *length; i++){
    if((list[i].x == (int)x) && (list[i].y == (int)y) && (list[i].z == (int)z)){
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
