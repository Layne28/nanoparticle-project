/* This program creates a nanoparticle of a particular shape.
   it writes out a configuration in xyz format, to be used
   as an initial configuration for 'etch.c'. 
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "etch.h"

int main(int argc, char *argv[]) {

  //Set  constants
  double d = 4.065; //(Au-Au spacing in Angstroms)
  double a = d*sqrt(2.0); //unit cell length
  int nside = 6; //number of unit cells in each dimension
  double L = a*nside; //length of side of cube
  int N = 4*nside*nside*nside; //total number of atoms (and vacancies)

  int i,j = 0;

  //pick shape
  char *input = "sphere";

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

  //Make fcc lattice
  fcc(lat, nside, a);

  //Create shape on the lattice
  shape(input, lat, a, nside);
  
  //Print out lattice as xyz file
  print_config(lat, nside);
}

/* OTHER FUNCTIONS ***************************************************************************/
/********************************************************************************************/

/* FUNCTIONS FOR MAKING INITIAL CONFIGURATION *************************************************/

void shape(char *name, nn_vec ***r, double a, int nside) {

  int i,j,k;
  double L = a*nside;
  int N = 4*nside*nside*nside;

  if(name == "two_points"){
    r[nside][nside][nside].occ = 1;
    r[nside][nside][nside-2].occ = 1;

    //get nearest neighbors
    for(i = 0; i < 2*nside; i++){
      for(j = 0; j < 2*nside; j++){
	for(k = 0; k < 2*nside; k++){
	  if((i+j+k)%2==0){
	    //Now a bunch of if statements in case you're on the edge of the lattice
	    //Could also do this with mod operator...
	      if(i>0 && j>0){
		if(r[i-1][j-1][k].occ == 1) r[i][j][k].n ++;
	      }if(i>0 && j<(2*nside-1)){
		if(r[i-1][j+1][k].occ == 1) r[i][j][k].n ++;
	      }if(i<(2*nside-1) && j>0){
		if(r[i+1][j-1][k].occ == 1) r[i][j][k].n ++;
	      }if(i<(2*nside-1) && j<(2*nside-1)){
		if(r[i+1][j+1][k].occ == 1) r[i][j][k].n ++;
	      }if(j>0 && k<(2*nside-1)){
		if(r[i][j-1][k+1].occ == 1) r[i][j][k].n ++;
	      }if(i>0 && k<(2*nside-1)){
		if(r[i-1][j][k+1].occ == 1) r[i][j][k].n ++;
	      }if(i<(2*nside-1) && k<(2*nside-1)){
		if(r[i+1][j][k+1].occ == 1) r[i][j][k].n ++;
	      }if(j<(2*nside-1) && k<(2*nside-1)){
		if(r[i][j+1][k+1].occ == 1) r[i][j][k].n ++;
	      }if(j>0 && k>0){
		if(r[i][j-1][k-1].occ == 1) r[i][j][k].n ++;
	      }if(i>0 && k>0){
		if(r[i-1][j][k-1].occ == 1) r[i][j][k].n ++;
	      }if(i<(2*nside-1) && k>0){
		if(r[i+1][j][k-1].occ == 1) r[i][j][k].n ++;
	      }if(j<(2*nside-1) && k>0){
		if(r[i][j+1][k-1].occ == 1) r[i][j][k].n ++;
	      }
	    
	  }
	}
      }
    }
  }

  if(name == "cube"){
    for(i = 0; i < 2*nside; i++){
      for(j = 0; j < 2*nside; j++){
	for(k = 0; k < 2*nside; k++){
	  if((i+j+k)%2==0){
	    r[i][j][k].occ = 1;
	  }
	}
      }
    }
  }

  if(name == "unit_cell"){
    r[nside][nside][nside].occ = 1;
    r[nside+1][nside+1][nside].occ = 1;
    r[nside][nside+1][nside+1].occ = 1;
    r[nside+1][nside][nside+1].occ = 1;

  //calculate nearest neighbors
    double nn_dist = a/sqrt(2.0);
    for(i = 0; i < 2*nside; i++){
      for(j = 0; j < 2*nside; j++){
	for(k = 0; k < 2*nside; k++){
	  if((i+j+k)%2==0){
	    //Now a bunch of if statements in case you're on the edge of the lattice
	      if(i>0 && j>0){
		if(r[i-1][j-1][k].occ == 1) r[i][j][k].n ++;
	      }if(i>0 && j<(2*nside-1)){
		if(r[i-1][j+1][k].occ == 1) r[i][j][k].n ++;
	      }if(i<(2*nside-1) && j>0){
		if(r[i+1][j-1][k].occ == 1) r[i][j][k].n ++;
	      }if(i<(2*nside-1) && j<(2*nside-1)){
		if(r[i+1][j+1][k].occ == 1) r[i][j][k].n ++;
	      }if(j>0 && k<(2*nside-1)){
		if(r[i][j-1][k+1].occ == 1) r[i][j][k].n ++;
	      }if(i>0 && k<(2*nside-1)){
		if(r[i-1][j][k+1].occ == 1) r[i][j][k].n ++;
	      }if(i<(2*nside-1) && k<(2*nside-1)){
		if(r[i+1][j][k+1].occ == 1) r[i][j][k].n ++;
	      }if(j<(2*nside-1) && k<(2*nside-1)){
		if(r[i][j+1][k+1].occ == 1) r[i][j][k].n ++;
	      }if(j>0 && k>0){
		if(r[i][j-1][k-1].occ == 1) r[i][j][k].n ++;
	      }if(i>0 && k>0){
		if(r[i-1][j][k-1].occ == 1) r[i][j][k].n ++;
	      }if(i<(2*nside-1) && k>0){
		if(r[i+1][j][k-1].occ == 1) r[i][j][k].n ++;
	      }if(j<(2*nside-1) && k>0){
		if(r[i][j+1][k-1].occ == 1) r[i][j][k].n ++;
	      }
	    
	  }
	}
      }
    }
  }  




  if(name == "sphere"){
    double R = 0.75*(L/2.0); //change this to adjust size of sphere
    double delta =  1e-5; //small buffer for comparison
    //loop through once to identify which atoms are in the sphere
    for(i = 0; i < 2*nside; i++){
      for(j = 0; j < 2*nside; j++){
	for(k = 0; k < 2*nside; k++){
	  if((i+j+k)%2==0){
	    double rmag = sqrt(dotpdt(&r[i][j][k].v, &r[i][j][k].v));
	    if(rmag < (R+delta)) r[i][j][k].occ = 1;
	  }
	}
      }
    } 
    //calculate nearest neighbors
    double nn_dist = a/sqrt(2.0);
    for(i = 0; i < 2*nside; i++){
      for(j = 0; j < 2*nside; j++){
	for(k = 0; k < 2*nside; k++){
	  if((i+j+k)%2==0){
	    //Now a bunch of if statements in case you're on the edge of the lattice
	      if(i>0 && j>0){
		if(r[i-1][j-1][k].occ == 1) r[i][j][k].n ++;
	      }if(i>0 && j<(2*nside-1)){
		if(r[i-1][j+1][k].occ == 1) r[i][j][k].n ++;
	      }if(i<(2*nside-1) && j>0){
		if(r[i+1][j-1][k].occ == 1) r[i][j][k].n ++;
	      }if(i<(2*nside-1) && j<(2*nside-1)){
		if(r[i+1][j+1][k].occ == 1) r[i][j][k].n ++;
	      }if(j>0 && k<(2*nside-1)){
		if(r[i][j-1][k+1].occ == 1) r[i][j][k].n ++;
	      }if(i>0 && k<(2*nside-1)){
		if(r[i-1][j][k+1].occ == 1) r[i][j][k].n ++;
	      }if(i<(2*nside-1) && k<(2*nside-1)){
		if(r[i+1][j][k+1].occ == 1) r[i][j][k].n ++;
	      }if(j<(2*nside-1) && k<(2*nside-1)){
		if(r[i][j+1][k+1].occ == 1) r[i][j][k].n ++;
	      }if(j>0 && k>0){
		if(r[i][j-1][k-1].occ == 1) r[i][j][k].n ++;
	      }if(i>0 && k>0){
		if(r[i-1][j][k-1].occ == 1) r[i][j][k].n ++;
	      }if(i<(2*nside-1) && k>0){
		if(r[i+1][j][k-1].occ == 1) r[i][j][k].n ++;
	      }if(j<(2*nside-1) && k>0){
		if(r[i][j+1][k-1].occ == 1) r[i][j][k].n ++;
	      }
	    
	  }
	}
      }
    }
  }  
}

void fcc(nn_vec ***r, int nside, double a) {

  //Follows Allen & Tildesley's prescription for creating an FCC lattice
  
  double L = nside*a;
  double d = 0.5*a;
  int i,j,k;
  int n=0;

  //Loop through the fcc lattice points (only points the sum of whose indices is even)
  for(i = 0; i < 2*nside; i++){
    for(j = 0; j < 2*nside; j++){
      if(i%2==0){
	if(j%2==0){
	  for(k = 0; k < 2*nside; k=k+2){
	    r[i][j][k].v.x = i*d - L/2;
	    r[i][j][k].v.y = j*d - L/2;
	    r[i][j][k].v.z = k*d - L/2;
	  }
	}if(j%2!=0){
	  for(k = 1; k < 2*nside; k=k+2){
	    r[i][j][k].v.x = i*d - L/2;
	    r[i][j][k].v.y = j*d - L/2;
	    r[i][j][k].v.z = k*d - L/2;
	  }
	}
      }
      if(i%2!=0){
	if(j%2==0){
	  for(k = 1; k < 2*nside; k=k+2){
	    r[i][j][k].v.x = i*d - L/2;
	    r[i][j][k].v.y = j*d - L/2;
	    r[i][j][k].v.z = k*d - L/2;
	  }
	}if(j%2!=0){
	  for(k = 0; k < 2*nside; k=k+2){
	    r[i][j][k].v.x = i*d - L/2;
	    r[i][j][k].v.y = j*d - L/2;
	    r[i][j][k].v.z = k*d - L/2;
	  }
	}
      }
    }
  }    
  

  /*  //Now enforce "periodic boundaries" (not really, it just makes the lattice into an actual square)
    for(i = 0; i < n; i++){

    if(r[i].v.x > L/2) r[i].v.x -= L;
    if(r[i].v.y > L/2) r[i].v.y -= L;
    if(r[i].v.z > L/2) r[i].v.z -= L;

    if(r[i].v.x < - L/2) r[i].v.x += L;
    if(r[i].v.y < - L/2) r[i].v.y += L;
    if(r[i].v.z < - L/2) r[i].v.z += L;
    } */
}

void print_config(nn_vec ***r, int nside) {

  int i,j,k;
  int N = 4*nside*nside*nside;

  FILE *f;
  f = fopen("/home/layne/research/nanoparticles/init_config.xyz", "w");

  fprintf(f, "%d\n",2*N); //vmd will think there are twice as many atoms as there actually are, due to unused array elements
  fprintf(f, "initial configuration of gold nanoparticle\n");
  for(i = 0; i < 2*nside; i++) {
    for(j = 0; j < 2*nside; j++){
      for(k = 0; k < 2*nside; k++){
	if((i+j+k)%2 != 0){
	  fprintf(f, "H\t%lf\t%lf\t%lf\t%d\t%d\n", r[i][j][k].v.x, r[i][j][k].v.y, r[i][j][k].v.z, r[i][j][k].n, r[i][j][k].occ);
	}
       	else{
	  if(r[i][j][k].occ == 0) fprintf(f, "H\t%lf\t%lf\t%lf\t%d\t%d\n", r[i][j][k].v.x, r[i][j][k].v.y, r[i][j][k].v.z, r[i][j][k].n, r[i][j][k].occ);
	  else  fprintf(f, "Au\t%lf\t%lf\t%lf\t%d\t%d\n", r[i][j][k].v.x, r[i][j][k].v.y, r[i][j][k].v.z, r[i][j][k].n, r[i][j][k].occ);
	}
      }
    }
  }
  fclose(f);
  f = NULL;    
}

/* EXTRA MATH FUNCTIONS *****************************************************************/
/*****************************************************************************************/

double dotpdt(vec3 *a, vec3 *b){

  double value = a->x*b->x + a->y*b->y + a->z*b->z;

  return value;   
}

