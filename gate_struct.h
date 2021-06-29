
/* gates */
enum gate_list {P, H, CNOT, CZ, SWAP, X, Y, Z};

struct gate_t {
  enum gate_list type;
  long q_i; 
  long q_j;
};
typedef struct gate_t gate;

/* product of gates */
struct gate_prod_t {
  long len;
  gate *g;
};
typedef struct gate_prod_t  gate_prod;

/* PZX form for a P, CZ, CNOT circuit :
   z_v * P_b * Z_B * X_A  */
struct PZX_form_t {
  long n;
  int *v;
  int *b;
  int **B;
  int **A;
};
typedef struct PZX_form_t PZX_form;

/* intermediate form of a stabilizer circuit :
   H_a * P_d * Z_D * H_w * exp(i*k*Pi/4) * x_u * z_v * P_b * Z_B * X_A */

struct normal_form_t {
  long n;
  int *a; 
  int *d;
  int **D;
  int *w;  // Omega in the paper
  int k;   // global phase phi is k*Pi/4
  int *u;
  int *v;
  int *b;
  int **B;
  int **A;
} ;
typedef struct normal_form_t normal_form;

struct graph_state_t {
  long n;
  int **B;
  int **B_red;
  int **A;
  int **A_inv;
  int *v;
};
typedef struct graph_state_t graph_state;
