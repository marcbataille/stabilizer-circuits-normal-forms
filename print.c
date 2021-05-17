#include <stdlib.h>
#include <stdio.h>

#include "gate_struct.h"
#include "constants.h"
#include "aux_ops.h"
#include "print.h"


void print_help(long n) {
  printf("-----------------------------------------------------------\n");
  printf("                           HELP                         \n");
  printf("-----------------------------------------------------------\n");
  printf("**Qubits are numbered from 0 to %ld.\n",n-1);
  printf("**Phase on qubit i --> Pi\n");
  printf("**Hadamard on qubit i --> Hi\n");
  printf("**CNOT with target on i and control on j --> iXj\n");
  printf("**CZ between qubit i and j --> iZj\n");
  printf("**SWAP between qubit i and j --> iSj\n");
  printf("**PauliX on qubit i --> Xi\n");
  printf("**PauliY on qubit i --> Yi\n");
  printf("**PauliZ on qubit i --> Zi\n\n");
 }
void print_error_cmd(void){
  fprintf(stderr, "Invalid arguments. Command usage is : \'stabnf N MODE\', " 
	  "where N is the number of qubits (2 <= N <= %d) and "
	  "MODE is \'man\' (for manual mode) or \'stat\' (for statistics mode).\n", MAX_QUBIT); 
}

void print_error_qb(long n) {
  printf("Invalid value : qubits are numbered from 0 to %ld. Please, try again.\n", n - 1);
}

void print_instructions_man(long n) {
  printf("\n=======================================================================\n");
  printf("             NORMAL FORM FOR A %ld-QUBIT STABILIZER CIRCUIT        \n",n);
  printf("=======================================================================\n");  
  printf("\n==========================   MANUAL MODE   ============================\n\n");    
  printf("\n-----------------------------------------------------------------------\n");
  printf("              ENTER THE CIRCUIT TO WRITE IN NORMAL FORM      \n");
  printf("-----------------------------------------------------------------------\n");
  printf("**Enter the gates of the circuit from left to right.\n");
  printf("**Press ENTER after each gate, including the last gate.\n");
  printf("**Press ctrl + D on Linux or ctrl + Z on Windows to complete the entry.\n");
  printf("**Qubits are numbered from 0 to %ld.\n",n-1);
  printf("**Phase on qubit i --> Pi\n");
  printf("**Hadamard on qubit i --> Hi\n");
  printf("**CNOT with target on i and control on j --> iXj\n");
  printf("**CZ between qubit i and j --> iZj\n");
  printf("**SWAP between qubit i and j --> iSj\n");
  printf("**PauliX on qubit i --> Xi\n");
  printf("**PauliY on qubit i --> Yi\n");
  printf("**PauliZ on qubit i --> Zi\n\n");
 }

void print_information_stat() {
  printf("\n=====================================   STATISTICS MODE   =======================================\n\n");    
  printf("INPUT CIRCUIT is a circuit of CZ and H gates implementing a graph state |G> = Z_B |+...+>  where \n");
  printf("Z_B is a product of CZ gates and |+> = H|0>.\n\n");
  printf("OUTPUT CIRCUIT is a circuit of Pauli-Z, CNOT, CZ and H gates, implementing the same graph state\n");
  printf("in the form |G> = Z_v X_A Z_B' |+...+>  where :\n");  
  printf("\t Z_v is a product of Pauli-Z gate,\n");
  printf("\t X_A is a product of CNOT gate,\n");
  printf("\t Z_B' is a product of CZ gate corresponding to a circuit of depth 1.\n");
  printf("\n================================================================================================\n");    
}

void print_matrix(int **M, long n) {
  long r, c;
  for (r = 0; r < n; ++r) {
    for (c = 0; c < n; ++c) {
      printf("%d ", M[r][c]);
    }
    printf("\n");
  }
}

void print_vector(int *vec, long n) {
  for (long r = 0; r < n; ++r) {
    printf("%d\n", vec[r]);
  }
}

void print_normal_form(normal_form *nf) {
  long n = nf -> n;
  long r, c;
  long offset = (n % 2) ? 0 : 1;
  printf("----------------------------------------------------------------------------------------------------------------------------\n");
  printf("                  NORMAL FORM : exp(i*k*Pi/4) * H_a * Z_u * P_d * Z_D * H_w * Z_v * P_b * Z_B * X_A\n");
  printf("----------------------------------------------------------------------------------------------------------------------------\n\n");
  if (n > MAX_PRINT_DIM) {
    printf("Sorry, matrices are too big to be printed.\n");
  }
  else {
    printf("k");
    printf("   ");
    printf("a");
    printf("   ");
    printf("u");
    printf("   ");
    printf("d");
    printf("   ");
    for (c = 0; c < n/2 - offset; ++c) {
      printf("  ");
    }
    printf("D ");
    for (c = 0; c < n/2; ++c) {
      printf("  ");
    }
    printf("   ");
    printf("w");
    printf("   ");
    printf("v");
    printf("   ");
    printf("b");
    printf("   ");
    for (c = 0; c < n/2 - offset; ++c) {
      printf("  ");
    }
    printf("B ");
    for (c = 0; c < n/2; ++c) {
      printf("  ");
    }
    printf("   ");
    for (c = 0; c < n/2 - offset; ++c) {
      printf("  ");
    }
    printf("A ");
    for (c = 0; c < n/2; ++c) {
      printf("  ");
    }
    printf("\n\n");
    for (r = 0; r < n; ++r) {
      if (r == 0){
	printf("%d",nf -> k);
      } else {
	printf(" ");
      }
      printf("   ");
      printf("%d", nf -> a[r]);
      printf("   ");
      printf("%d", nf -> u[r]);
      printf("   ");
      printf("%d", nf -> d[r]);
      printf("   ");
      for (c = 0; c < n; ++c) {
	printf("%d ", nf -> D[r][c]);
      }
      printf("   ");
      printf("%d", nf -> w[r]);
      printf("   ");
      printf("%d", nf -> v[r]);
      printf("   ");
      printf("%d", nf -> b[r]);
      printf("   ");
      for (c = 0; c < n; ++c) {
	printf("%d ", nf -> B[r][c]);
      }
      printf("   ");
      for (c = 0; c < n; ++c) {
	printf("%d ", nf -> A[r][c]);
      }
      printf("\n");
    }
    printf("\n");
  }
}

 
void print_input(gate_prod  *input, long n) {
  long pos;
  printf("INPUT CIRCUIT     : --/%ld--", n);
  for (pos = 0; pos < input -> len; ++pos) {
    switch (input -> g[pos].type) {
    case P : printf("P%ld--", input -> g[pos].q_i);
      break;
    case H : printf("H%ld--", input -> g[pos].q_i);
      break;
    case X : printf("X%ld--", input -> g[pos].q_i);
      break;
    case Y : printf("Y%ld--", input -> g[pos].q_i);
      break;
    case Z : printf("Z%ld--", input -> g[pos].q_i);
      break;
    case CNOT : printf("CNOT[%ld,%ld]--", input -> g[pos].q_i, input -> g[pos].q_j);
      break;
    case CZ : printf("CZ{%ld,%ld}--", input -> g[pos].q_i, input -> g[pos].q_j);
      break;
    case SWAP : printf("SWAP{%ld,%ld}--", input -> g[pos].q_i, input -> g[pos].q_j);
      break;
    }
  }
  printf("(%ld gates)\n", input -> len);
}

void print_circuit_nf(normal_form *nf, gate_prod *CNOT_prod, int **A_aux) {
  long n = nf -> n;
  long r, c;
  long pos, len = 0;
  int ok;
  int a_aux [n];
  int w_aux [n];
  vector_cp(nf -> a, a_aux, n);
  vector_cp(nf -> w, w_aux, n);
  /* simplify the Hadamard gates for printing */
  for (r = 0; r < n; ++r) {
    if (nf -> a[r] == 1 && nf -> d[r] == 0 && nf -> u[r] == 0) {
      ok = 1;
      for(c = 0; c < n; ++c) {
	if (nf -> D[r][c] != 0) {
	  ok = 0;
	  break;
	}
      }
      if (ok) {
	a_aux[r] = 0;
	w_aux[r] = 0;
      }
    }
  }
  printf("NF CIRCUIT        : --/%ld--", n); 
  matrix_cp(nf -> A, A_aux, n);
  decompose_GL_matrix(A_aux, CNOT_prod, n);
  for (pos = CNOT_prod -> len - 1; pos >= 0; -- pos) {
    printf("CNOT[%ld,%ld]--", CNOT_prod -> g[pos].q_i, CNOT_prod -> g[pos].q_j);
    ++len;
  }
  for (r = 0; r < n; ++r) {
    for (c = r + 1; c < n; ++c) { 
      if (nf -> B[r][c]) {
	printf("CZ{%ld,%ld}--", r, c);
	++len;
      }
    }
  }
  for (r = 0; r < n; ++r) {
    if (nf -> b[r]) {
      printf("P%ld--", r);
      ++len;
    }
  }
  for (r = 0; r < n; ++r) {
    if (nf -> v[r]) {
      printf("Z%ld--", r);
      ++len;
    }
  }
  
  for (r = 0; r < n; ++r) {
    if (w_aux[r]) {
      printf("H%ld--", r);
      ++len;
    }
  }
  for (r = 0; r < n; ++r) {
    for (c = r + 1; c < n; ++c) { 
      if (nf -> D[r][c]) {
	printf("CZ{%ld,%ld}--", r, c);
	++len;
      }
    }
  }
  for (r = 0; r < n; ++r) {
    if (nf -> d[r]) {
      printf("P%ld--", r);
      ++len;
    }
  }
  for (r = 0; r < n; ++r) {
    if (nf -> u[r]) {
      printf("Z%ld--", r);
      ++len;
    }
  }
  for (r = 0; r < n; ++r) {
    if (a_aux[r]) {
      printf("H%ld--", r);
      ++len;
    }
  }
  printf("(%ld gates, phase = %d*Pi/4)\n", len, nf -> k);
}


void print_stabilizer_state(int *u, int *v, int *w, int **B, long n) {
  int is_identity = 1;
  printf("\n---------------------------------------  STABILIZER STATE AND GRAPH STATE  ----------------------------------------\n");
  printf("\nThe stabilizer state |S> resulting from applying your circuit to the state |0...O> is |S> = H Z P |G> ,  where :\n\n");
  printf("H = ");
  for (long i = 0; i < n; ++i) {
    if (u[i]) {
      printf("H%ld ",i);
      is_identity = 0;
    }
  }
  if (is_identity) {
    printf("Id");
  }
  is_identity = 1;
  printf("\n\nZ = ");
  for (long i = 0; i < n; ++i) {
    if (v[i]) {
      printf("Z%ld ",i);
      is_identity = 0;
    }
  }
  if (is_identity) {
    printf("Id");
  }
  is_identity = 1;
  printf("\n\nP = ");
  for (long i = 0; i < n; ++i) {
    if (w[i]) {
      printf("P%ld ",i);
      is_identity = 0;
    }
  }
  if (is_identity) {
    printf("Id");
  }
  
  printf("\n\n|G> is the graph state corresponding to the graph G = { ");
  for (long i = 0; i < n; ++i) {
    for (long j = i + 1; j < n; ++j) {
      if (B[i][j]) {
	printf("{%ld,%ld},",i,j);
      }
    }
  }
  printf("\b }\n");
  if (n <= MAX_PRINT_DIM) {
  printf("\nThe adjacency matrix of G is : \n");
  print_matrix(B,n);
  }
}
  
void print_graph_state(graph_state *gs, gate_prod *CNOT_prod, int **A_aux) {
printf("\n|G> can be written in the form |G> = Z_v X_A Z_B' |+...+>  where :\n");  
  printf("\t Z_v is a product of Pauli-Z gate,\n");
  printf("\t X_A is a product of CNOT gate,\n");
  printf("\t Z_B' is a product of CZ gate corresponding to a circuit of depth 1.\n");
  long n = gs -> n;
  matrix_cp(gs -> A, A_aux, n);
  decompose_GL_matrix(A_aux, CNOT_prod, n);
  printf("\n|G> = ");
  for (long r = 0; r < n; ++r) {
    if (gs -> v[r]) {
      printf("Z%ld ", r);
    }
  }
  for (long pos = 0; pos < CNOT_prod -> len; ++pos){
    printf("CNOT[%ld,%ld] ", CNOT_prod -> g[pos].q_i, CNOT_prod -> g[pos].q_j);
  }
  for (long r = 0; r < n; ++r) {
    for (long c = r + 1; c < n; ++c) {
      if (gs -> B_red[r][c]) {
	printf("CZ{%ld,%ld} ", r, c);
      }
    }
  }
  printf("|+...+>\n");
}
