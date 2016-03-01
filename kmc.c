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
#define T 0.075//0.0259
#define esub 3.930 //sublimation energy
#define epsilon 0.3275 //energy per bond
#define npass 8750000
#define d 1e-6 //small factor used when converting to int
//#define num_frames 500//Don't make this too high or you'll crash the computer.
#define num_frames 100

double mu = (-6.5)*epsilon;
double beta = 1.0/T;
int Nsurf = 0, Nvac = 0; //number of surface and vacany atoms
int Nsurf0 = 0, Nvac0 = 0; //To keep track of initial numbers;
int Ntot = 0;
//int num_frames = 0;
int frame_rate = 0;
int time = 0;

int temp_mc=0;

int main(void){
  int i,j,k;

  int Ngrid;
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
  in = fopen("init_config.xyz", "r");
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
  //num_frames = Ntot/frame_rate +2;
  frame_rate = Ntot/(num_frames-2);
  fprintf(stderr,"Frame rate: %d Number of frames: %d\n", frame_rate, num_frames);

  printf("Total number of atoms at start: %d\n", Ntot);
  printf("number of surface atoms: %d\n", Nsurf);
  printf("number of vacancies: %d\n", Nvac);

  Nsurf0 = Nsurf;
  Nvac0 = Nvac;

  /* Fill arrays with the indices of the locations of surface and vacany atoms *********************************************************************************/

  triple *surf = (triple *) calloc((size_t)(2*Nsurf), sizeof(triple)); //make these twice the size they need to be to accomodate the possibility
  triple *vac = (triple *) calloc((size_t)(2*Nvac), sizeof(triple)); //of getting more than you started with
  nbs_vec *queue = (nbs_vec *) calloc((size_t)(2*Nsurf), sizeof(nbs_vec));
  // checks whether or not a lattice site is in the queue
  int ***queue_query = (int ***) calloc((size_t)(2*nl), sizeof(int**)); //make these twice the size they need to be to accomodate the possibility
  for (i=0;i<2*nl;i++) {
    queue_query[i] = (int **) calloc((size_t)(2*nl), sizeof(int*)); 
    for (j=0; j<2*nl; j++) {
      queue_query[i][j] = (int *) calloc((size_t)(2*nl), sizeof(int)); 
    }
  }
  fprintf(stderr,"Allocated queue\n");

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
	      push(lat, queue, &q, i, j, k, lat[i][j][k].n);
        queue_query[i][j][k]=1;
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
  fprintf(stderr,"Populated queue\n");
  //print_queue(lat,queue,q);
  i = 0;

  //Create an array that stores the number of surface atoms at every step
  /* Now do Monte Carlo ****************************************************/
  fprintf(stderr,"Starting the Kinetic Monte Carlo run.\n");
  kmc(lat, queue, queue_query, &q);
  fprintf(stderr,"Finished Kinetic Monte Carlo run.\n");


  //Print information//////////////////////////////////////////////////////////
  //print_final(lat); //print final configuration
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
      free(queue_query[i][j]);
    }
    free(lat[i]);
    free(queue_query[i]);
  }
  free(lat);

  return 1; //end of program
}

/* OTHER FUNCTIONS ***********************************************************************/
/*****************************************************************************************/

void kmc_peel(nn_vec ***r, nbs_vec *queue, int ***queue_query, int *q){
  nbs_vec *neighbs = (nbs_vec *) calloc((size_t)NN, sizeof(nbs_vec));
	while(*q > 0){
    int icnt=0;
    int n0 = r[queue[0].x][queue[0].y][queue[0].z].n;
    while (r[queue[icnt].x][queue[icnt].y][queue[icnt].z].n==n0 && icnt<*q-1) {
      icnt++;
    }
    long unsigned int index;
    if (icnt>0) {
      index = gsl_rng_uniform_int(rg, (unsigned long int)icnt);
    } else {
      index = 0;
    }
	  neighbors(r, neighbs, queue[index].x, queue[index].y, queue[index].z);
    r[queue[index].x][queue[index].y][queue[index].z].occ = 0;
    queue_query[queue[index].x][queue[index].y][queue[index].z] = 0;
	  pop(queue, q, (int)index);
	  for(int i = 0; i < NN; i++) {
	    int x = neighbs[i].x;
      int y = neighbs[i].y;
      int z = neighbs[i].z;
      int occ = neighbs[i].occ;
      r[x][y][z].n--;
      if (occ==1) {
        // need to change neighbor number on the queue
        if (queue_query[x][y][z]==1) {
          int nindex = 0;
          while (!(queue[nindex].x==x && queue[nindex].y==y && queue[nindex].z==z)) {
            nindex++;
          }
          pop(queue, q, nindex);
          push(r,queue, q, x, y, z, r[x][y][z].n);
        }
        else {
          push(r, queue, q, x, y, z, r[x][y][z].n);
          queue_query[x][y][z] = 1;
        }           
      }
	  }

    // dump the trajectory
    Ntot--;
    if(time%frame_rate == 0) {
      fprintf(stderr,"Printing at time %d, queue length %d, Ntot %d\n",time, *q, Ntot);
      print_conf(r, time, *q);
    }
    time++;
  }
}

void kmc(nn_vec ***r, nbs_vec *queue, int ***queue_query, int *q){
  nbs_vec *neighbs = (nbs_vec *) calloc((size_t)NN, sizeof(nbs_vec));
  double  *event_list = (double *) calloc((size_t)(8*nl*nl*nl), sizeof(double));
  while (*q>0) {
    double R = 0;
    event_list[0] = 0;
    for (int i=0; i<*q; i++) {
      R += exp(-beta*(epsilon*queue[i].n-beta*mu)); // compute the total rate
      event_list[i+1]=R;
    }
    double event = gsl_rng_uniform_pos(rg);
    int index=0;
    while (!(event >= event_list[index]/R && event <= event_list[index+1]/R)) {
      index++;
    }
    //u = gsl_rng_uniform_pos(rg);
    //dt = log(1/u)/R;
    //printf("queue indices: %d %d %d %d\n", queue[index].x, queue[index].y, queue[index].z, queue[index].n);
	  neighbors(r, neighbs, queue[index].x, queue[index].y, queue[index].z);
    r[queue[index].x][queue[index].y][queue[index].z].occ = 0;
	  pop(queue, q, index);
	  for(int i = 0; i < NN; i++) {
	    int x = neighbs[i].x;
      int y = neighbs[i].y;
      int z = neighbs[i].z;
      int occ = neighbs[i].occ;
      r[x][y][z].n--;
      if (occ==1) {
        // need to change neighbor number on the queue
        if (queue_query[x][y][z]==1) {
          int nindex = 0;
          while (!(queue[nindex].x==x && queue[nindex].y==y && queue[nindex].z==z)) {
            nindex++;
          }
          pop(queue, q, nindex);
          push(r,queue, q, x, y, z, r[x][y][z].n);
        }
        else {
          push(r, queue, q, x, y, z, r[x][y][z].n);
          queue_query[x][y][z] = 1;
        }           
      }
	  }

    // dump the trajectory
    Ntot--;
    if(time%frame_rate == 0) {
      print_conf(r, time, *q);
    }
    time++;
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

void print_queue(nn_vec ***lat, nbs_vec *queue, int q) {
  int x,y,z,n;
  int valid = 1, n0=0;
  for (int i=0; i<q; i++) {
    x = queue[i].x;
    y = queue[i].y;
    z = queue[i].z;
    n = lat[x][y][z].n;
    fprintf(stderr,"%d %d %d %d %d\n",i,x,y,z,n);
    if (n<n0) {
      valid = 0;
    }
    n0 = lat[x][y][z].n;
  }
  if (valid==0) exit(-1);
}

void print_conf(nn_vec ***r, int t, int N) {

  int i,j,k;
  FILE *f;
  char filename[80];
  sprintf(filename,"conf_%07d.xyz",t);
  f = fopen(filename, "w");

  fprintf(f, "%d\n",N); //vmd will think there are twice as many atoms as there actually are, due to unused array elements
  fprintf(f, "initial configuration of gold nanoparticle\n");
  for(i = 0; i < 2*nl; i++) {
    for(j = 0; j < 2*nl; j++) {
      for(k = 0; k < 2*nl; k++) {
        if (r[i][j][k].occ == 1 && r[i][j][k].n < NN) {
          fprintf(f, "Au\t%d\t%d\t%d\t%d\t%d\n", (i-nl), (j-nl), (k-nl), r[i][j][k].n, r[i][j][k].occ);
        } //else fprintf(f, "Au\t%d\t%d\t%d\t%d\t%d\n", 0,0,0,0,0);
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

void pop(nbs_vec *list, int *length, int index){
  for(int i = index; i < *length; i++){
    list[i].x = list[i+1].x;
    list[i].y = list[i+1].y;
    list[i].z = list[i+1].z;
    list[i].n = list[i+1].n;
  }
  *length = *length - 1;
}

void push(nn_vec ***lat, nbs_vec *queue, int *length, int x, int y, int z, int n){
	int n0 = queue[0].n;
	int cnt = 0;

  if(*length == 0) {
    queue[0].x = x;
    queue[0].y = y;
    queue[0].z = z;
    queue[0].n = n;
    *length = *length + 1;
  }
  
  else{
    while(n>n0 && cnt < *length){
      n0 = lat[queue[cnt].x][queue[cnt].y][queue[cnt].z].n;
      cnt++;
    }
    // corner case
    if (cnt==*length && n>n0) {
      queue[*length].x = x;
      queue[*length].y = y;
      queue[*length].z = z;
      queue[*length].n = n;
      *length = *length+1;
    }
    else if (cnt==*length && n<n0) {
      queue[*length].x = queue[*length-1].x;
      queue[*length].y = queue[*length-1].y;
      queue[*length].z = queue[*length-1].z;
      queue[*length].n = queue[*length-1].n;
      queue[*length-1].x = x;
      queue[*length-1].y = y;
      queue[*length-1].z = z;
      queue[*length-1].n = n;
      *length = *length+1;
    }
    else {
      if (cnt>0) {
        cnt -= 1;
      }
      for(int i = *length-1; i >= cnt; i--){
        queue[i+1].x = queue[i].x;
        queue[i+1].y = queue[i].y;
        queue[i+1].z = queue[i].z;
        queue[i+1].n = queue[i].n;
      }
      queue[cnt].x = x;
      queue[cnt].y = y;
      queue[cnt].z = z;
      queue[cnt].n = n;
      *length = *length + 1;
    }
    // check insertion
    for (int i=cnt+1;i<*length;i++) {
      if (queue[i].n < queue[cnt].n) {
        fprintf(stderr,"i: %d cnt: %d in: %d cntn: %d length: %d\n",i,cnt,queue[i].n,queue[cnt].n, *length);
        exit(0);
      }
    }
  }
}


