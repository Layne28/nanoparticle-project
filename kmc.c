/* "etch.c" */
/* This program takes in a nanoparticle configuration and applies MCMC to it
   in order to "etch" it. 
*/

#include <stdio.h>
#include <stdlib.h>
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
#define npass 8750000
#define d 1e-6 //small factor used when converting to int
#define num_frames 200 //Don't make this too high or you'll crash the computer.

double mu = (-6.5)*epsilon;
double beta = 1.0/T;
//int config_freq = npass/num_frames;
int config_freq = 1000;
int Nsurf = 0, Nvac = 0; //number of surface and vacany atoms
int Nsurf0 = 0, Nvac0 = 0; //To keep track of initial numbers;
int Ntot = 0;

int temp_mc=0;

int main(void){

  printf("%d\n", config_freq);

  int Ngrid;
  int i,j,k;
  long l;
  int time = 0;
  long mc_total = 0; //number of total attempts
  //long mc_accept = 0;
  long mc_insert = 0; //number of insertion moves
  long mc_delete = 0; //number of deletion moves
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

  FILE *in;
  in = fopen("/home/layne/research/nanoparticle-project/init_config.xyz", "r");
  if(in == NULL) printf("problem!\n");
  if(fscanf(in, "%d\n", &Ngrid) != 1) printf("Error: need 1 number.\n");

  //allocate memory (for lattice)
  //The lattice is 3d, each side indexes from 0 to 2*nl-1
  nn_vec ***lat;
  lat = (nn_vec ***)calloc((size_t)2*nl, sizeof(nn_vec**));
  for(i = 0; i < 2*nl; i++){
    lat[i] = (nn_vec **) calloc((size_t)2*nl, sizeof(nn_vec*));
    for(j = 0; j < 2*nl; j++){
      lat[i][j] = (nn_vec*) calloc((size_t)2*nl, sizeof(nn_vec));
    }
  }

  if(fgets(ignore, sizeof(ignore), in) == NULL) printf("End of file."); //skip line

  char atom_curr[2];
  i = 0;
  j = 0;
  k = 0;
  double temp1;
  double temp2;
  double temp3;

  //Now actually read in data
  for(i = 0; i < 2*nl; i++) {
    for(j = 0; j < 2*nl; j++){
      for(k = 0; k < 2*nl; k++){
	if(fscanf(in, "%s\t" "%lf\t" "%lf\t" "%lf\t" "%d\t"  "%d\n", atom_curr, &temp1, &temp2, &temp3, &lat[i][j][k].n, &lat[i][j][k].occ) != 6) printf("Error: need 6 inputs.\n");
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
    traj[l] = (nn_vec ***) calloc((size_t)2*nl, sizeof(nn_vec**));
    for(i = 0; i < 2*nl; i++){
      traj[l][i] = (nn_vec **) calloc((size_t)2*nl, sizeof(nn_vec*));
      for(j = 0; j < 2*nl; j++){
	traj[l][i][j] = (nn_vec *) calloc((size_t)2*nl, sizeof(nn_vec));
      }
    }
  }

  for(i = 0; i < 2*nl; i++){
    for(j = 0; j < 2*nl; j++){
      for(k = 0; k < 2*nl; k++){
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

  triple *surf = (triple *) calloc((size_t)(2*Nsurf), sizeof(triple)); //make these twice the size they need to be to accomodate the possibility
  triple *vac = (triple *) calloc((size_t)(2*Nvac), sizeof(triple)); //of getting more than you started with
  nbs_vec *queue = (nbs_vec *) calloc((size_t)(2*Nsurf), sizeof(nbs_vec));

  int ls=0; 
  int lv=0;
  int q=0;
  for(i = 0; i < 2*nl; i++) {
    for(j = 0; j < 2*nl; j++){
      for(k = 0; k < 2*nl; k++){
	if((i+j+k)%2 == 0){
	  if(lat[i][j][k].occ==1){
	    if(lat[i][j][k].n < NN){
	      surf[ls].x = i;
	      surf[ls].y = j;
	      surf[ls].z = k;
	      push(queue, &q, i, j, k, lat[i][j][k].n);
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
 
  i = 0;

  //Create an array that stores the number of surface atoms at every step

  /* Now do Monte Carlo ****************************************************/

  kmc(lat, traj, queue, &q);
//  exit(-1);
  printf("%d\n", queue[0].n);


  //Print information//////////////////////////////////////////////////////////
  print_traj(traj); //print whole trajectory
//  print_final(lat); //print final configuration
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
  for(i = 0; i < 2*nl; i++){
    for(j = 0; j < 2*nl; j++){
      free(lat[i][j]);
    }
    free(lat[i]);
  }
  free(lat);

  for(l = 0; l < (npass/config_freq+2); l++){
    for(i = 0; i < 2*nl; i++){
      for(j = 0; j < 2*nl; j++){
	      free(traj[l][i][j]);
      }
      free(traj[l][i]);
    }
    free(traj[l]);
  }
  free(traj);

  return 1; //end of program
}

/* OTHER FUNCTIONS ***********************************************************************/
/*****************************************************************************************/

void kmc(nn_vec ***r, nn_vec ****traj, nbs_vec *queue, int *q){

  int time = 0;

	while(*q > 0){
    int icnt = 0;
    int in = queue[0].n;
    while(in == queue[0].n){
            in = queue[icnt].n;
            icnt++;
    }
    icnt--; //you go ahead by one too many indices in the while loop
    long unsigned int index = gsl_rng_uniform_int(rg, (unsigned long int)icnt);

	  nbs_vec *neighbs = (nbs_vec *) calloc((size_t)NN, sizeof(nbs_vec));
    printf("queue indices: %d %d %d %d\n", queue[index].x, queue[index].y, queue[index].z, queue[index].n);
	  neighbors(r, neighbs, queue[index].x, queue[index].y, queue[index].z);
    r[queue[index].x][queue[index].y][queue[index].z].occ = 0;

    int x0 = queue[0].x;
    int y0 = queue[0].y;
    int z0 = queue[0].z;
    int n0 = queue[0].n;
    queue[0].x = queue[index].x;
    queue[0].y = queue[index].y;
    queue[0].z = queue[index].z;
    queue[0].n = queue[index].n;
    queue[index].x = x0;
    queue[index].y = y0;
    queue[index].z = z0;
    queue[index].n = n0;

	  pop(queue, q);
	  for(int i = 0; i < NN; i++) {
	    int x = neighbs[i].x;
      int y = neighbs[i].y;
      int z = neighbs[i].z;
      int occ = neighbs[i].occ;
      int n = neighbs[i].n;
      r[x][y][z].n --;
	    if(occ == 1 && n == NN){
        push(queue, q, x, y, z, n-1);
      }           
	  }

//   printf("%d\n", queue[0].n);

  int i,j,k;
//  if((time)%config_freq == 0){
  if(time < npass/config_freq){
  for(i = 0; i < 2*nl; i++){
	  for(j = 0; j < 2*nl; j++){
	    for(k = 0; k < 2*nl; k++){
	      traj[time][i][j][k].n = r[i][j][k].n;
	      traj[time][i][j][k].occ = r[i][j][k].occ;
	    }
	  }
	} 
// }
 }
  time ++;
  }
}





void neighbors(nn_vec ***r, nbs_vec *nbs, int i, int j, int k){
  if(i==0 || j==0 || k==0 || i==(2*nl-1) || j==(2*nl-1) || k==(2*nl-1)){
    printf("Error: lattice has grown to edge of bounding cube. Exiting program.\n");
    exit(1);
  }    
  if(i>0 && j>0){
    nbs[0].x = i-1;
    nbs[0].y = j-1;
    nbs[0].z = k;
    nbs[0].n = r[i-1][j-1][k].n;
    nbs[0].occ = r[i-1][j-1][k].occ;
  }if(i>0 && j<(2*nl-1)){
    nbs[1].x = i-1;
    nbs[1].y = j+1;
    nbs[1].z = k;
    nbs[1].n = r[i-1][j+1][k].n;
    nbs[1].occ = r[i-1][j+1][k].occ;
  }if(i<(2*nl-1) && j>0){
    nbs[2].x = i+1;
    nbs[2].y = j-1;
    nbs[2].z = k;
    nbs[2].n = r[i+1][j-1][k].n;
    nbs[2].occ = r[i+1][j-1][k].occ;
  }if(i<(2*nl-1) && j<(2*nl-1)){
    nbs[3].x = i+1;
    nbs[3].y = j+1;
    nbs[3].z = k;
    nbs[3].n = r[i+1][j+1][k].n;
    nbs[3].occ = r[i+1][j+1][k].occ;
  }if(j>0 && k<(2*nl-1)){
    nbs[4].x = i;
    nbs[4].y = j-1;
    nbs[4].z = k+1;
    nbs[4].n = r[i][j-1][k+1].n;
    nbs[4].occ = r[i][j-1][k+1].occ;
  }if(i>0 && k<(2*nl-1)){
    nbs[5].x = i-1;
    nbs[5].y = j;
    nbs[5].z = k+1;
    nbs[5].n = r[i-1][j][k+1].n;
    nbs[5].occ = r[i-1][j][k+1].occ;
  }if(i<(2*nl-1) && k<(2*nl-1)){
    nbs[6].x = i+1;
    nbs[6].y = j;
    nbs[6].z = k+1;
    nbs[6].n = r[i+1][j][k+1].n;
    nbs[6].occ = r[i+1][j][k+1].occ;
  }if(j<(2*nl-1) && k<(2*nl-1)){
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
  }if(i<(2*nl-1) && k>0){
    nbs[10].x = i+1;
    nbs[10].y = j;
    nbs[10].z = k-1;
    nbs[10].n = r[i+1][j][k-1].n;
    nbs[10].occ = r[i+1][j][k-1].occ;
   }if(j<(2*nl-1) && k>0){
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
  int N = 4*nl*nl*nl;

  FILE *f;
  f = fopen("/home/layne/research/nanoparticle-project/traj.xyz", "w");

  for(l=0; l<length; l++){
    fprintf(f,"%d\n", N); //this time, since we aren't using them again, just record the actual atoms
    fprintf(f, "Step %ld\n", l+1);

    for(i = 0; i < 2*nl; i++) {
      for(j = 0; j < 2*nl; j++){
	for(k = 0; k < 2*nl; k++){
	  if((i+j+k)%2 == 0){
	    if(t[l][i][j][k].occ == 1 && t[l][i][j][k].n < 12)  fprintf(f, "Au\t%d\t%d\t%d\n", i-nl, j-nl, k-nl);
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

void print_final(nn_vec ***r) {

  int i,j,k;
  int N = 4*nl*nl*nl;

  FILE *f;
  f = fopen("/home/layne/research/nanoparticle-project/final_config.xyz", "w");

  fprintf(f, "%d\n",2*N); //vmd will think there are twice as many atoms as there actually are, due to unused array elements
  fprintf(f, "initial configuration of gold nanoparticle\n");
  for(i = 0; i < 2*nl; i++) {
    for(j = 0; j < 2*nl; j++){
      for(k = 0; k < 2*nl; k++){
	if((i+j+k)%2 != 0){
	  fprintf(f, "H\t%d\t%d\t%d\t%d\t%d\n", (i-nl), (j-nl), (k-nl), r[i][j][k].n, r[i][j][k].occ);
	}
       	else{
	  if(r[i][j][k].occ == 0) fprintf(f, "H\t%d\t%d\t%d\t%d\t%d\n", (i-nl), (j-nl), (k-nl), r[i][j][k].n, r[i][j][k].occ);
	  else  fprintf(f, "Au\t%d\t%d\t%d\t%d\t%d\n", (i-nl), (j-nl), (k-nl), r[i][j][k].n, r[i][j][k].occ);
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

/* Queue Functions */

void pop(nbs_vec *list, int *length){
  for(int i = 0; i < *length; i++){
    list[i].x = list[i+1].x;
    list[i].y = list[i+1].y;
    list[i].z = list[i+1].z;
    list[i].n = list[i+1].n;
  }
  *length = *length - 1;
}

void push(nbs_vec *list, int *length, int x, int y, int z, int n){
	int n0 = list[0].n;
	int cnt = 0;

  if(*length == 0){
    list[0].x = x;
    list[0].y = y;
    list[0].z = z;
    list[0].n = n;
    *length = *length + 1;
  }
  
  else{

            while(n>=n0 && cnt < *length){
              n0 = list[cnt].n;
              cnt++;
          }
          for(int i = *length-1; i >= cnt; i--){
            list[i+1].x = list[i].x;
            list[i+1].y = list[i].y;
            list[i+1].z = list[i].z;
            list[i+1].n = list[i].n;
          }
          *length = *length + 1;

          list[cnt].x = x;
          list[cnt].y = y;
          list[cnt].z = z;
          list[cnt].n = n;

  }

}


