/* Header file for nanoparticle etching project */

#ifndef ETCH_H
#define ETCH_H

#define nl 11
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
  int n;
  int occ;
} nn_vec;

//Define a triple of integers (so that you don't have to make a quadruple array)
typedef struct {
  int x;
  int y;
  int z;
} triple;

typedef struct {
	int x;
	int y;
	int z;
	int n;
	int occ;
} nbs_vec;

	

// Initialize functions

void shape(const char *name, nn_vec ***r);
void print_config(nn_vec ***r);
void print_config_t(nn_vec ***r);
void print_traj(nn_vec ****t);
void print_queue(nn_vec ***lat, nbs_vec *queue, int q);
void print_final(nn_vec ***r);
void fcc(nn_vec ***r);
double dotpdt(vec3 *a, vec3 *b); //Calculates the dot product of two 3D vectors and returns the value calculated
int mc_move(nn_vec ***r, triple *surf, triple *vac);
int mc_move_glauber(nn_vec ***r, triple *surf, triple *vac);
void neighbors(nn_vec ***r, nbs_vec *nbs, int i, int j, int k);

//For manipulating linked list
void add(triple *list, int *length, int x, int y, int z);
void delete(triple *list, int *length, int x, int y, int z);

//For manipulating queue
void pop(nbs_vec *list, int *length, int index);
void push(nn_vec ***lat, nbs_vec *queue, int *length, int x, int y, int z, int n);


void kmc_peel(nn_vec ***r, nn_vec ****traj, nbs_vec *queue, int ***queue_query, int *q);
void kmc(nn_vec ***r, nn_vec ****traj, nbs_vec *queue, int *q);

void record_nn(nn_vec ***r);


//Test functions
int test_one_insert(nn_vec ***r, triple *surf, triple *vac);
int test_two_delete(nn_vec ***r, triple *surf, triple *vac);
int test_two_insert(nn_vec ***r, triple *surf, triple *vac);
int test_loneatom(nn_vec ***r, triple *surf, triple *vac);
int mc_move_test(nn_vec ***r, triple *surf, triple *vac);

#endif
