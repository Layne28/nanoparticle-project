/* Analyze configurations generated from etching */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "etch_simple.h"

int nx, ny, nz;
int N;
int confnum = 10;

void fill_layer(nn_vec ***lat, nbs_vec *occup, nbs_vec *layer, int axis, int nlayer);
nbs_vec * get_perimeter(nbs_vec *layer);

int main(void){

  int i,j,k;
  int xmax=nl, ymax=nl, zmax=nl, xmin=nl, ymin=nl, zmin=nl;

  /* Start by loading configuration file & allocating memory */

  FILE *in;
  char *shape = malloc(1024);
  char * filename = malloc(1024);

  sprintf(filename, "conf_%04d.xyz", confnum);
  in = fopen(filename, "r");
  if(in == NULL){
    printf("Couldn't open file; exiting program.\n");
    exit(-1);
  }
  if(fscanf(in, "%d\n", &N) != 1) printf("Error: need 1 number.\n");
  printf("number of atoms: %d\n", N);
  if(fscanf(in, "%s\n", shape) == 0) printf("Error: need shape name.\n"); 
  if(strcmp(shape, "cylinder")==0){
    nx = 2*nl;
    ny = 2*nl;
    nz = aspect_ratio*2*nl;
  }
  else{
    nx = 2*nl;
    ny = 2*nl;
    nz = 2*nl;
  }

  //allocate lattice memory
  nn_vec ***lat;
  lat = (nn_vec ***)calloc((size_t)nx, sizeof(nn_vec**));
  for(i = 0; i < nx; i++){
    lat[i] = (nn_vec **) calloc((size_t)ny, sizeof(nn_vec*));
    for(j = 0; j < ny; j++){
      lat[i][j] = (nn_vec*) calloc((size_t)nz, sizeof(nn_vec));
    }
  }

  nbs_vec *occup;
  occup = (nbs_vec *)calloc((size_t)N, sizeof(nbs_vec));

 //read in data
 //while doing so, get maximum x,z,y indices
 int q = 0;
 char atom_curr[2];
 for(q=0; q<N; q++){
 if(fscanf(in, "%s\t" "%d\t" "%d\t" "%d\t" "%d\t" "%d\n", atom_curr, &occup[q].x, &occup[q].y, &occup[q].z, &occup[q].n, &occup[q].occ) != 6) printf("Error: need 6 inputs.\n");
  occup[q].x += nx/2;
  occup[q].y += ny/2;
  occup[q].z += nz/2;
  i = occup[q].x;
  j = occup[q].y;
  k = occup[q].z;
  lat[i][j][k].n = occup[q].n;
  lat[i][j][k].occ = occup[q].occ;
  if(lat[i][j][k].occ == 1){
    if(i > xmax) xmax = i;
    if(i < xmin) xmin = i;
    if(j > ymax) ymax = j;
    if(j < ymin) ymin = j;
    if(k > zmax) zmax = k;
    if(k < zmax) zmin = k;
  }
 }
 printf("zmax: %d\n", zmax - nz/2);

 fclose(in);
 in = NULL;

 //Now actually do things
 nbs_vec *top = (nbs_vec *)calloc((size_t)nx*ny, sizeof(nbs_vec));
 nbs_vec *next = (nbs_vec *)calloc((size_t)nx*ny, sizeof(nbs_vec));
 fill_layer(lat, occup, top, 3, zmax);
 fill_layer(lat, occup, next, 3, zmax-1);

 print_shape(lat, shape);

 return 1;
 
}

/////////////////////////////////////////////////////

void fill_layer(nn_vec ***lat, nbs_vec *occup, nbs_vec *layer, int axis, int nlayer){
  int q;
  int p=0;
  if(axis==3){ //look at z axis
    for(q=0; q<N; q++){
      if(occup[q].z == nlayer){
        layer[p].x = occup[q].x;
        layer[p].y = occup[q].y;
        layer[p].z = occup[q].z;
        layer[p].n = occup[q].n;
        layer[p].occ = 2;
        lat[layer[p].x][layer[p].y][layer[p].z].occ = 2;
        p++;
      }
    }
  }
}

void print_shape(nn_vec ***r, const char *input) {

  int i,j,k;

  FILE *f;
  char *filename = malloc(1024);
  sprintf(filename, "ana_%04d.xyz", confnum);
  f = fopen(filename, "w");

  fprintf(f, "%d\n",N); //vmd will think there are twice as many atoms as there actually are, due to unused array elements
  fprintf(f, "%s\n", input);
  for(i = 0; i < nx; i++) {
    for(j = 0; j < ny; j++){
      for(k = 0; k < nz; k++){
        if((i+j+k)%2==0){
	        if(r[i][j][k].occ==1){
            fprintf(f, "Au\t%d\t%d\t%d\t%d\t%d\n", i-nx/2, j-ny/2, k-nz/2, r[i][j][k].n, r[i][j][k].occ);
          }
          else if(r[i][j][k].occ==2){
            fprintf(f, "Ag\t%d\t%d\t%d\t%d\t%d\n", i-nx/2, j-ny/2, k-nz/2, r[i][j][k].n, r[i][j][k].occ);          
          }
        }
        }
      }
    }

  fclose(f);
  f = NULL;    
}

