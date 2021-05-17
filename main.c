#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include <unistd.h>
#include <math.h>

#include "constants.h"
#include "gate_struct.h"
#include "print.h"
#include "aux_ops.h"
#include "cases.h"


int main(int argc, char **argv){

  /* SCANNING COMMAND */
  
  /* argv[0] = n, the number of qubits */
  long n;
  /* argv[2] = "man" for manual mode, "stat" for statistics mode */
  enum {MANUAL, STATISTICS} mode;
  /* to store the first invalid char of a string that represents a long int */
  char *invalid_n;
  if (argc != 3) {
    print_error_cmd();
    exit(EXIT_FAILURE);
  }
  n = strtol(argv[1], &invalid_n, 10);
  if (*invalid_n != '\0') {
    print_error_cmd();
    exit(EXIT_FAILURE);
  }
  if (n < 2 || n > MAX_QUBIT) {
    print_error_cmd();
    exit(EXIT_FAILURE);
  }
  if (strcmp(argv[2],"man")==0) {
    mode = MANUAL;
  } else if (strcmp(argv[2],"stat")==0) {
    mode = STATISTICS;
  } else {
    print_error_cmd();
    exit(EXIT_FAILURE);
  }
  
  long i, j; // qubits
  char *invalid_i, *invalid_j;
  long r; // r for row, c for column
  int sf; // scanf return value
  char buffer[16], buff_i[16], buff_j[16]; // buffers for user entry
  char q_or_c[2]={'c','\0'};  // quit or continue

  /* MEMORY ALLOCATION */

  /* Memory allocation for user's input */
  gate gtab[MAX_INPUT];
  gate_prod *input = (gate_prod *) malloc(sizeof(gate_prod));;
  input -> g = gtab;
  
  /* Memory allocation for a product of P, CZ and CNOT gates.
     This product is an input of the C_to_PZX algorithm in cases 2.2.1 and 3.4.1 (see paper) */
  long max_PZX = n*(n-1)/2 + 3*(n-1); // the length of this circuit is alway less than max_PZX (see paper)
  gate_prod *PZX_prod = (gate_prod *)malloc(sizeof(gate_prod));
  PZX_prod -> g = (gate *) calloc((unsigned long)max_PZX, sizeof(gate));
  
  /* Memory allocation for products of CNOT gates. */
  
  /* CNOT_prod is an array of less than n^2 gates, used as input of the C_to_PZX algorithm 
     in cases 2.2.1 and 3.4.1 (see paper). 
     Also used to store the decomposition of an invertible matrix */
  gate_prod *CNOT_prod = (gate_prod *) malloc(sizeof(gate_prod));
  CNOT_prod -> g = (gate *) calloc((unsigned long)(n * n), sizeof(gate));
  
    
  /* Memory allocation for a GL matrix A_aux */
  int **A_aux = (int**) calloc((unsigned long)n, sizeof(int*));
  for (r = 0; r < n; ++r){
    A_aux[r] = (int *) calloc((unsigned long)n, sizeof(int));
  }

  /* Memory allocation for a symmetric matrix B_aux */
  int **B_aux = (int**) calloc((unsigned long)n, sizeof(int*));
  for (r = 0; r < n; ++r){
    B_aux[r] = (int *) calloc((unsigned long)n, sizeof(int));
  }
  
  /* Memory allocation for the PZX form */
  PZX_form *PZX = (PZX_form *) malloc(sizeof(PZX_form));
  PZX -> n = n;
  PZX -> v = (int *) calloc((unsigned long)n, sizeof(int));
  PZX -> b = (int *) calloc((unsigned long)n, sizeof(int));
  PZX -> B = (int **) calloc((unsigned long)n, sizeof(int*));
  PZX -> A = (int **) calloc((unsigned long)n, sizeof(int*));
  for (r = 0; r < n; ++r){
    PZX -> B[r] = (int *) calloc((unsigned long)n, sizeof(int));
    PZX -> A[r] = (int *) calloc((unsigned long)n, sizeof(int));
  }
  
  /* Memory allocation for the normal form */
  normal_form *nf = (normal_form *) malloc(sizeof(normal_form));
  nf -> n = n;
  nf -> a = (int *) calloc((unsigned long)n, sizeof(int));
  nf -> d = (int *) calloc((unsigned long)n, sizeof(int));  
  nf -> w = (int *) calloc((unsigned long)n, sizeof(int)); 
  nf -> u = (int *) calloc((unsigned long)n, sizeof(int)); 
  nf -> v = (int *) calloc((unsigned long)n, sizeof(int));
  nf -> b = (int *) calloc((unsigned long)n, sizeof(int));
  nf -> D = (int **) calloc((unsigned long)n, sizeof(int*));
  nf -> B = (int **) calloc((unsigned long)n, sizeof(int*));
  nf -> A = (int **) calloc((unsigned long)n, sizeof(int*));
  for (r = 0; r < n; ++r){
    nf -> D[r] = (int *) calloc((unsigned long)n, sizeof(int));
    nf -> B[r] = (int *) calloc((unsigned long)n, sizeof(int));
    nf -> A[r] = (int *) calloc((unsigned long)n, sizeof(int));
  }
  
  /* Memory allocation for the graph state */
  graph_state *gs = (graph_state *) malloc(sizeof(graph_state));
  gs -> n = n;
  gs -> v = (int *) calloc((unsigned long)n, sizeof(int));
  gs -> B = (int **) calloc((unsigned long)n, sizeof(int*));
  gs -> B_red = (int **) calloc((unsigned long)n, sizeof(int*));
  gs -> A = (int **) calloc((unsigned long)n, sizeof(int*));
  gs -> A_inv = (int **) calloc((unsigned long)n, sizeof(int*));
  for (r = 0; r < n; ++r){
    gs -> B[r] = (int *) calloc((unsigned long)n, sizeof(int));
    gs -> B_red[r] = (int *) calloc((unsigned long)n, sizeof(int));
    gs -> A[r] = (int *) calloc((unsigned long)n, sizeof(int));
    gs -> A_inv[r] = (int *) calloc((unsigned long)n, sizeof(int));
  }

  switch (mode) {
  case MANUAL :
    print_instructions_man(n);
    /* reducing a new circuit */
    while (q_or_c[0] == 'c') { 
      input -> len = 0;
      initialize_nf(nf);
      initialize_gs(gs);
      /* starting the induction process :
	 during the while loop, the gate entered by the user is added to the normal form */
      while (1) {
	print_input(input, n);
	print_circuit_nf(nf, CNOT_prod, A_aux);
	printf("gate %ld ? ", input -> len + 1);
	/* scanning user input */
	sf=scanf("%15s",buffer);
	if (sf == EOF) {
	  printf("End of circuit. \n\n");
	  break;
	}
	if (sscanf(buffer,"H%s",buff_i) == 1) {
	  i = strtol(buff_i, &invalid_i, 10);
	  if (*invalid_i != '\0') {
	    print_error_qb(n);
	  } else if (i < 0 || i >= n) {
	    print_error_qb(n);
	  } else { // Hadamard was scanned
	    printf("gate %ld --> H%ld\n", input -> len + 1, i);
	    input -> g[input -> len].type = H;
	    input -> g[input -> len].q_i = i;
	    ++input -> len;
	    merge_Hi_with_nf(i, nf);
	  }
	} else if (sscanf(buffer,"P%s",buff_i) == 1) {
	  i = strtol(buff_i, &invalid_i, 10);
	  if (*invalid_i != '\0') {
	    print_error_qb(n);
	  } else if (i < 0 || i >= n) {
	    print_error_qb(n);
	  } else { // Phase was scanned
	    printf("gate %ld --> P%ld\n", input -> len + 1, i);
	    input -> g[input -> len].type = P;
	    input -> g[input -> len].q_i = i;
	    ++input -> len;
	    merge_Pi_with_nf(i, nf, PZX, PZX_prod, CNOT_prod);
	  }
	} else if (sscanf(buffer,"%[0123456789]X%s", buff_i, buff_j) == 2) {
	  i = strtol(buff_i, &invalid_i, 10);
	  j = strtol(buff_j, &invalid_j, 10);
	  if (*invalid_i != '\0' || *invalid_j !='\0') {
	    print_error_qb(n);
	  } else if (i < 0 || i >= n || j < 0 || j >= n || i==j) {
	    print_error_qb(n);
	  } else { // CNOT was scanned
	    printf("gate %ld --> CNOT[%ld,%ld]\n", input -> len + 1, i, j);
	    input -> g[input -> len].type = CNOT;
	    input -> g[input -> len].q_i = i;
	    input -> g[input -> len].q_j = j;
	    ++input -> len;
	    merge_CNOTij_with_nf(i, j, nf, PZX, PZX_prod, CNOT_prod);
	  }
	} else if (sscanf(buffer,"%[0123456789]Z%s", buff_i, buff_j) == 2) {
	  i = strtol(buff_i, &invalid_i, 10);
	  j = strtol(buff_j, &invalid_j, 10);
	  if (*invalid_i != '\0' || *invalid_j !='\0') {
	    print_error_qb(n);
	  } else if (i < 0 || i >= n || j < 0 || j >= n || i==j) {
	    print_error_qb(n);
	  } else { // CZ was scanned
	    printf("gate %ld --> CZ{%ld,%ld}\n", input -> len + 1, i, j);
	    input -> g[input -> len].type = CZ;
	    input -> g[input -> len].q_i = i;
	    input -> g[input -> len].q_j = j;
	    ++input -> len;
	    merge_CZij_with_nf(i, j, nf, PZX, PZX_prod, CNOT_prod);
	  }
	} else if (sscanf(buffer,"%[0123456789]S%s", buff_i, buff_j) == 2) {
	  i = strtol(buff_i, &invalid_i, 10);
	  j = strtol(buff_j, &invalid_j, 10);
	  if (*invalid_i != '\0' || *invalid_j !='\0') {
	    print_error_qb(n);
	  } else if (i < 0 || i >= n || j < 0 || j >= n || i==j) {
	    print_error_qb(n);
	  } else { // SWAP was scanned
	    printf("gate %ld --> SWAP{%ld,%ld}\n", input -> len + 1, i, j);
	    input -> g[input -> len].type = SWAP;
	    input -> g[input -> len].q_i = i;
	    input -> g[input -> len].q_j = j;
	    ++input -> len;
	    merge_SWAPij_with_nf(i, j, nf);
	  }
	} else if (sscanf(buffer,"X%s",buff_i) == 1) {
	  i = strtol(buff_i, &invalid_i, 10);
	  if (*invalid_i != '\0') {
	    print_error_qb(n);
	  } else if (i < 0 || i >= n) {
	    print_error_qb(n);
	  } else { // Pauli-X was scanned
	    printf("gate %ld --> X%ld\n", input -> len + 1, i);
	    input -> g[input -> len].type = X;
	    input -> g[input -> len].q_i = i;
	    ++input -> len;
	    merge_Xi_with_nf(i, nf, PZX, PZX_prod, CNOT_prod);
	  }
	} else if (sscanf(buffer,"Y%s",buff_i) == 1) {
	  i = strtol(buff_i, &invalid_i, 10);
	  if (*invalid_i != '\0') {
	    print_error_qb(n);
	  } else if (i < 0 || i >= n) {
	    print_error_qb(n);
	  } else { // Pauli-Y was scanned
	    printf("gate %ld --> Y%ld\n", input -> len + 1, i);
	    input -> g[input -> len].type = Y;
	    input -> g[input -> len].q_i = i;
	    ++input -> len;
	    merge_Yi_with_nf(i, nf, PZX, PZX_prod, CNOT_prod);
	  }
	} else if (sscanf(buffer,"Z%s",buff_i) == 1) {
	  i = strtol(buff_i, &invalid_i, 10);
	  if (*invalid_i != '\0') {
	    print_error_qb(n);
	  } else if (i < 0 || i >= n) {
	    print_error_qb(n);
	  } else { // Pauli-Z was scanned
	    printf("gate %ld --> Z%ld\n", input -> len + 1, i);
	    input -> g[input -> len].type = Z;
	    input -> g[input -> len].q_i = i;
	    ++input -> len;
	    merge_Zi_with_nf(i, nf, PZX, PZX_prod, CNOT_prod);
	  }
	} else {
	  printf("Invalid gate. Please check syntax.\n");
	  print_help(n);
	}  
	while (getchar() != '\n');
      }
      /* end of the induction process */
      /* copy vector a of the normal form nf before simplification (for stabilizer state) */
      int a[n];
      vector_cp(nf -> a, a, n);
      /****************************/
      simplify_nf(nf);
      print_normal_form(nf);
      print_stabilizer_state(a, nf -> u, nf -> d, nf -> D, n);
      matrix_cp(nf -> D, B_aux, n);
      compute_gs(gs, B_aux, CNOT_prod);
      print_graph_state(gs, CNOT_prod, A_aux);
      /************************************************************************************************************/
      printf("\nEnter \'c\' to scan another circuit or \'q\' to quit program.\n--> ");
      sf=scanf("%1s",q_or_c);
      while (getchar() != '\n');
      while (q_or_c[0] != 'c' && q_or_c[0] != 'q') {
	printf("Enter \'c\' to scan another circuit or \'q\' to quit program.\n--> ");
	sf=scanf("%1s",q_or_c);
	while (getchar() != '\n');
      }
    }
    break;
  case STATISTICS :
    #define INVALID_MESSAGE "Invalid value."
    print_information_stat();
    long sample_size, input_cz_len, output_2q_len;
    long max_input_len=n*(n-1)/2;
    long output_max_2q_len; 
    long output_min_2q_len; 
    long sum_output_2q_len; 
    int invalid, done; 
    long k, sum_null_val;
    while (q_or_c[0] == 'c') {
      printf("\nEnter the sample size of INPUT circuits (max %d) : ", MAX_SAMPLE);
      if (scanf("%ld", &sample_size) != 1) { 
	invalid = 1;
      } else {
	invalid = sample_size > MAX_SAMPLE || sample_size < 1;
      }
      if (invalid) {
	printf(INVALID_MESSAGE "\n");
	break;
      }
      printf("\n\nFor a %ld-qubit register, max length for CZ circuits of distinct CZ gates is %ld.\n", n, max_input_len);
      printf("Enter the CZ gate count of the INPUT circuits : ");
      if (scanf("%ld", &input_cz_len) != 1) { 
	invalid = 1;
      } else {
	invalid = input_cz_len > max_input_len || input_cz_len < 1;
      }
      if (invalid) {
	printf(INVALID_MESSAGE "\n");
	break;
      }
      output_max_2q_len = 0;
      output_min_2q_len = n*n;
      sum_output_2q_len = 0;
      srandom((unsigned int)time(NULL));
      for (long s = 0; s < sample_size; ++s) {
	for (long i = 0; i < n; ++i) {
	  for (long j = 0; j < n; ++j) {
	    B_aux[i][j]=0;
	  }
	}
	for (long t = 0; t < input_cz_len; ++t) {
	  k = random() % (max_input_len -t) + 1;
	  sum_null_val = 0;
	  done = 0;
	  for (long i = 0; i < n; ++i) {
	    if (done) break;
	    for (long j = i+1; j < n; ++j) {
	      sum_null_val += 1-B_aux[i][j];
	      if (sum_null_val >= k) {
		done = 1;
		B_aux[i][j]=1;
		B_aux[j][i]=1;
		break;
	      }
	    }
	  }
	}
	if (n <= MAX_PRINT_DIM) {
	  printf("\nMatrix of graph %ld : \n",s+1);
	  print_matrix(B_aux,n);
	}
	initialize_gs(gs);
	compute_gs(gs, B_aux, CNOT_prod);
	count_gate_gs(gs, CNOT_prod, A_aux, &output_2q_len);
	sum_output_2q_len += output_2q_len;
 	if (output_2q_len > output_max_2q_len) {
	  output_max_2q_len = output_2q_len;
	}
	if (output_2q_len < output_min_2q_len) {
	  output_min_2q_len = output_2q_len;
	}
      }
      printf("\n********************************************************************************");
      printf("\n*      STATISTICS FOR A SAMPLE OF %3ld GRAPH STATES (%3ld-QUBIT REGISTER)        *",
	     sample_size, n);
      printf("\n********************************************************************************");
      printf("\n*      2-QUBIT COUNT OF INPUT CIRCUITS :             %6ld -->  %4d %%        *", input_cz_len, 100);
      printf("\n********************************************************************************");
      printf("\n*      2-QUBIT COUNT OF OUTPUT CIRCUITS :        MAX %6ld --> %5.1lf %%        *",
	     output_max_2q_len,  (double)output_max_2q_len/(double)input_cz_len*100.);
      printf("\n*                                                AVG %6.0lf --> %5.1lf %%        *",
	     (double)sum_output_2q_len/(double)sample_size, (double)sum_output_2q_len/(double)sample_size/(double)input_cz_len*100.);
      printf("\n*                                                MIN %6ld --> %5.1lf %%        *",
	     output_min_2q_len, (double)output_min_2q_len/(double)input_cz_len*100.);
      printf("\n********************************************************************************");
      printf("\n*      AVERAGE GAIN OF OUTPUT OVER INPUT :           %6.0lf --> %5.1lf %%        *",
	     (double)input_cz_len-(double)sum_output_2q_len/(double)sample_size > 0 ? (double)input_cz_len-(double)sum_output_2q_len/(double)sample_size : 0.,
	     (double)input_cz_len-(double)sum_output_2q_len/(double)sample_size > 0 ? ((double)input_cz_len-(double)sum_output_2q_len/(double)sample_size)/(double)input_cz_len*100 : 0.);
      printf("\n********************************************************************************");
      printf("\nEnter \'c\' for a new sample or \'q\' to quit program.\n--> ");
      sf=scanf("%1s", q_or_c);
      while (getchar() != '\n');
      while (q_or_c[0] != 'c' && q_or_c[0] != 'q') {
	printf("Enter \'c\' to scan another circuit or \'q\' to quit program.\n--> ");
	sf=scanf("%1s",q_or_c);
	while (getchar() != '\n');
      }
    }
  }
  printf("Exiting program.\n");
  
  /* free memory for input circuit */
  free(input);
  
  /* free memory for normal form */
  free(nf -> a);
  free(nf -> d);
  free(nf -> w);
  free(nf -> u);
  free(nf -> v);
  free(nf -> b);
  for (r = 0; r < n; ++r) {
    free(nf -> D[r]);
    free(nf -> B[r]);
    free(nf -> A[r]);
  }
  free(nf -> D);
  free(nf -> B);
  free(nf -> A);
  free(nf);
    
  /* free memory for PZX */
  free(PZX -> v);
  free(PZX -> b);
  for (r = 0; r < n; ++r) {
    free(PZX -> B[r]);
    free(PZX -> A[r]);
  }
  free(PZX -> B);
  free(PZX -> A);
  free(PZX);
  
  /* free memory for gate products */
  free(PZX_prod -> g);
  free(PZX_prod);
  free(CNOT_prod -> g);
  free(CNOT_prod);
    
    
  /* free memory for matrix A_aux */
  for (r = 0; r < n; ++r) {
    free(A_aux[r]);
  }
  free(A_aux);

  /* free memory for matrix B_aux */
  for (r = 0; r < n; ++r) {
    free(B_aux[r]);
  }
  free(B_aux);


  /* free memory for graph state */
  free(gs -> v);
  for (r = 0; r < n; ++r){
    free(gs -> B[r]);
    free(gs -> B_red[r]);
    free(gs -> A[r]);
    free(gs -> A_inv[r]);
  }
  free(gs -> B);
  free(gs -> B_red);
  free(gs -> A);
  free(gs -> A_inv);
  free(gs);

}

    
  
    


