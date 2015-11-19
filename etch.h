/* Header file for nanoparticle etching project */

#ifndef ETCH_H
#define ETCH_H

#define NN 12

//Define a 3-vector
typedef struct {
  double x;
  double y;
  double z;
} vec3;

/*Define a nearest neighbors "vector"
  (real vector + #nearest neighbors + 0 for vacant, 1 for occupied)
*/
typedef struct {
  vec3 v;
  int n;
  int occ;
} nn_vec;

//Define a triple of integers (so that you don't have to make a quadruple array)
typedef struct {
  int x;
  int y;
  int z;
} triple;

// Initialize functions

void shape(char *name, nn_vec ***r, double a, int nside);
void print_config(nn_vec ***r, int nside);
void print_config_t(nn_vec ***r, int time);
void print_traj(nn_vec ****t);
void fcc(nn_vec ***r, int nside, double a);
double dotpdt(vec3 *a, vec3 *b); //Calculates the dot product of two 3D vectors and returns the value calculated
int mc_move(nn_vec ***r, triple *surf, triple *vac);
void neighbors(nn_vec ***r, nn_vec *nbs, int i, int j, int k);
void add(triple *list, int *length, double x, double y, double z);
void delete(triple *list, int *length, double x, double y, double z);

#endif
