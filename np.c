/* This program creates a nanoparticle of a particular shape.
   it writes out a configuration in xyz format, to be used
   as an initial configuration for 'etch.c'. 
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "etch_simple.h"

int nx;
int ny;
int nz;
const char *initial_shape;

int main(void){

  //Set  constants
  //  double d = 4.065; //(Au-Au spacing in Angstroms)
  int i,j = 0;

  //pick shape
  initial_shape = "cube";

  nx = 2*nl;
  ny = 2*nl;
  if(strcmp(initial_shape, "cylinder")==0) nz = aspect_ratio*2*nl;
  else nz = 2*nl;

  //allocate memory (for lattice)
  //The lattice is 3d, each side indexes from 0 to 2*nl-1
  //Except for the spherocylinder configuration, where one side is 4 times the length of the other 
  nn_vec ***lat;

  if(strcmp(initial_shape, "cylinder")==0){
    lat = (nn_vec ***)calloc((size_t)nx, sizeof(nn_vec**));
    for(i = 0; i < nx; i++){
      lat[i] = (nn_vec **) calloc((size_t)ny, sizeof(nn_vec*));
      for(j = 0; j < ny; j++){
        lat[i][j] = (nn_vec*) calloc((size_t)nz, sizeof(nn_vec));
      }
    } 
  }
  else{
    lat = (nn_vec ***)calloc((size_t)nx, sizeof(nn_vec**));
    for(i = 0; i < nx; i++){
     lat[i] = (nn_vec **) calloc((size_t)ny, sizeof(nn_vec*));
      for(j = 0; j < ny; j++){
        lat[i][j] = (nn_vec*) calloc((size_t)nz, sizeof(nn_vec));
      }
    }
  }

  //Create shape on the lattice
  shape(initial_shape, lat);
  
  //Print out lattice as xyz file
  print_shape(lat, initial_shape);
}

/* OTHER FUNCTIONS ***************************************************************************/
/********************************************************************************************/

/* FUNCTIONS FOR MAKING INITIAL CONFIGURATION *************************************************/

void shape(const char *name, nn_vec ***r) {

  int i,j,k;

  if(strcmp(name, "five")==0){
    r[nl][nl][nl].occ = 1;
    r[nl+1][nl+1][nl].occ = 1;
    r[nl+1][nl-1][nl].occ = 1;
    r[nl-1][nl-1][nl].occ = 1;
    r[nl+2][nl][nl].occ = 1;
    record_nn(r);
  }

  if(strcmp(name, "three")==0){
    r[nl+1][nl+1][nl].occ = 1;
    r[nl+1][nl-1][nl].occ = 1;
    r[nl+2][nl][nl].occ = 1;
    record_nn(r);
  }

  if(strcmp(name, "two")==0){
    r[nl][nl][nl].occ = 1;
    r[nl][nl+1][nl-1].occ = 1;
    record_nn(r);
  }

  if(strcmp(name, "one")==0){
    r[nl][nl][nl].occ = 1;
    record_nn(r);
  }

  if(strcmp(name, "cylinder")==0){
    int Ntot = 0;
    
    double d = 4.065; //(Au-Au spacing in Angstroms)
    double a = d*sqrt(2.0); //unit cell length
    double c = 0.5*a;
    double L = a*nl; //length of side of cube
    double R = 0.80*(L/2.0); //change this to adjust size of cube
    double delta =  1e-5; //small buffer for comparison

    for(i = 0; i < nx; i++){
      for(j = 0; j < ny; j++){
      	for(k = 0; k < nz; k++){
	        if((i+j+k)%2==0){
	          vec3 temp;
	          temp.x = i*c - L/2;
	          temp.y = j*c - L/2;
	          temp.z = k*c - aspect_ratio*L/2;
            double r2 = sqrt(temp.x*temp.x + temp.y*temp.y);
            if(fabs(temp.z) < (aspect_ratio*R - R + delta)){
	            if(r2 < R+delta){
	              r[i][j][k].occ = 1;
	              Ntot++;
	           }
            }
            else{
              if(r2 < sqrt(R*R - (fabs(temp.z) - 3*R)*(fabs(temp.z) - 3*R))){
                r[i][j][k].occ = 1;
                Ntot++;
              }
            }
	        }
	      }
      }
    }
    record_nn(r);
    printf("Ntot: %d\n", Ntot);
  }

  if(strcmp(name, "cube")==0){
    int Ntot = 0;

    double d = 4.065; //(Au-Au spacing in Angstroms)
    double a = d*sqrt(2.0); //unit cell length
    double c = 0.5*a;
    double L = a*nl; //length of side of cube
    double R = 0.80*(L/2.0); //change this to adjust size of cube
    double delta =  1e-5; //small buffer for comparison

    for(i = 0; i < nx; i++){
      for(j = 0; j < ny; j++){
      	for(k = 0; k < nz; k++){
	        if((i+j+k)%2==0){
	          vec3 temp;
	          temp.x = i*c - L/2;
	          temp.y = j*c - L/2;
	          temp.z = k*c - L/2;
	          if(fabs(temp.x) < R+delta && fabs(temp.y) < R+delta && fabs(temp.z) < R+delta){
	            r[i][j][k].occ = 1;
	            Ntot++;
	          }
	        }
	      }
      }
    }
    record_nn(r);
    printf("Ntot: %d\n", Ntot);
  }

  if(strcmp(name, "unit_cell")==0){
    r[nl][nl][nl].occ = 1;
    r[nl+1][nl+1][nl].occ = 1;
    r[nl][nl+1][nl+1].occ = 1;
    r[nl+1][nl][nl+1].occ = 1;

    record_nn(r);
  }  

  if(strcmp(name, "sphere")==0){
    int Ntot = 0;

    double d = 4.065; //(Au-Au spacing in Angstroms)
    double a = d*sqrt(2.0); //unit cell length
    double c = 0.5*a;
    double L = a*nl; //length of side of cube
    double R = 0.8*(L/2.0); //change this to adjust size of sphere
    double delta =  1e-5; //small buffer for comparison

    //loop through once to identify which atoms are in the sphere
    for(i = 0; i < nx; i++){
      for(j = 0; j < ny; j++){
      	for(k = 0; k < nz; k++){
	        if((i+j+k)%2==0){
	          vec3 temp;
      	    temp.x = i*c - L/2;
	          temp.y = j*c - L/2;
	          temp.z = k*c - L/2;
	          double rmag = sqrt(dotpdt(&temp, &temp));
      	    if(rmag < (R+delta)){
	           r[i][j][k].occ = 1;
	           Ntot++;
	          }
	        }
	      }
      }
    } 
    record_nn(r);
    printf("Ntot: %d\n", Ntot);
  }

}

void record_nn(nn_vec ***r){
  int i,j,k;
  for(i = 0; i < nx; i++){
    for(j = 0; j < ny; j++){
      for(k = 0; k < nz; k++){
        if((i+j+k)%2==0){
          //Now a bunch of if statements in case you're on the edge of the lattice
            if(i>0 && j>0){
	      if(r[i-1][j-1][k].occ == 1) r[i][j][k].n ++;
	    }if(i>0 && j<(ny-1)){
		if(r[i-1][j+1][k].occ == 1) r[i][j][k].n ++;
	    }if(i<(nx-1) && j>0){
		if(r[i+1][j-1][k].occ == 1) r[i][j][k].n ++;
	    }if(i<(nx-1) && j<(ny-1)){
		if(r[i+1][j+1][k].occ == 1) r[i][j][k].n ++;
            }if(j>0 && k<(nz-1)){
		if(r[i][j-1][k+1].occ == 1) r[i][j][k].n ++;
            }if(i>0 && k<(nz-1)){
		if(r[i-1][j][k+1].occ == 1) r[i][j][k].n ++;
	    }if(i<(nx-1) && k<(nz-1)){
		if(r[i+1][j][k+1].occ == 1) r[i][j][k].n ++;
	    }if(j<(ny-1) && k<(nz-1)){
		if(r[i][j+1][k+1].occ == 1) r[i][j][k].n ++;
            }if(j>0 && k>0){
		if(r[i][j-1][k-1].occ == 1) r[i][j][k].n ++;
            }if(i>0 && k>0){
		if(r[i-1][j][k-1].occ == 1) r[i][j][k].n ++;
            }if(i<(nx-1) && k>0){
		if(r[i+1][j][k-1].occ == 1) r[i][j][k].n ++;
            }if(j<(ny-1) && k>0){
		if(r[i][j+1][k-1].occ == 1) r[i][j][k].n ++;
            }
	    
        }
      }
    }
  }
}

void print_shape(nn_vec ***r, const char *input) {

  int i,j,k;
  int N = nx*ny*nz;

  FILE *f;
  f = fopen("/home/layne/research/nanoparticle-project/init_config.xyz", "w");

  fprintf(f, "%d\n",N); //vmd will think there are twice as many atoms as there actually are, due to unused array elements
  fprintf(f, "%s\n", input);
  for(i = 0; i < nx; i++) {
    for(j = 0; j < ny; j++){
      for(k = 0; k < nz; k++){
      	if((i+j+k)%2 != 0){
	        fprintf(f, "H\t%d\t%d\t%d\t%d\t%d\n", i, j, k, r[i][j][k].n, r[i][j][k].occ);
      	}
       	else{
	        if(r[i][j][k].occ == 0) fprintf(f, "H\t%d\t%d\t%d\t%d\t%d\n", i, j, k, r[i][j][k].n, r[i][j][k].occ);
	        else  fprintf(f, "Au\t%d\t%d\t%d\t%d\t%d\n", i, j, k, r[i][j][k].n, r[i][j][k].occ);
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
