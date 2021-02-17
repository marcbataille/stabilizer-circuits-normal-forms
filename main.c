#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include <unistd.h>

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
  enum mode_ls {manual, statistics};
  enum mode_ls mode;
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
    mode = manual;
  } else if (strcmp(argv[2],"stat")==0) {
    mode = statistics;
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

  /* for stat mode */
  long sample_size;
  long input_len, nf_len, cz_red_len; // length of circuit
  long input_2q_len, nf_2q_len, cz_red_2q_len;  // number of 2 qubit gates
  long sum_input_2q, sum_nf, sum_cz_red, sum_nf_2q, sum_cz_red_2q;
  long s, t; // loop counter  

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
  
  /* This product is to store the matrix A_red_D such that D_red=A_red_D^T*D*A_red_D (see paper) 
     under the form of a transvection product. */
  gate_prod *A_red_D_prod = (gate_prod *) malloc(sizeof(gate_prod));
  A_red_D_prod -> g = (gate *) calloc((unsigned long)(n * n), sizeof(gate));
  
  /* This product is to store the matrix A_red_B such that B_red=A_red_B^T*B*A_red_B (see paper) 
     under the form of a transvection product. */
  gate_prod *A_red_B_prod = (gate_prod *) malloc(sizeof(gate_prod));
  A_red_B_prod -> g = (gate *) calloc((unsigned long)(n * n), sizeof(gate));
  
  /* Memory allocation for the inverse matrices of A_red_D and A_red_B */
  int **A_red_D_inv = (int**) calloc((unsigned long)n, sizeof(int*));
  int **A_red_B_inv = (int**) calloc((unsigned long)n, sizeof(int*));
  for (r = 0; r < n; ++r){
    A_red_D_inv[r] = (int *) calloc((unsigned long)n, sizeof(int));
    A_red_B_inv[r] = (int *) calloc((unsigned long)n, sizeof(int));
  }
  
  /* Memory allocation for a GL matrix A_aux */
  int **A_aux = (int**) calloc((unsigned long)n, sizeof(int*));
  for (r = 0; r < n; ++r){
    A_aux[r] = (int *) calloc((unsigned long)n, sizeof(int));
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
  
  /* Memory allocation for the CZ-reduced normal form */
  CZ_red_normal_form *red_nf = (CZ_red_normal_form *) malloc(sizeof(CZ_red_normal_form));
  red_nf -> n = n;
  red_nf -> a = (int *) calloc((unsigned long)n, sizeof(int));
  red_nf -> d = (int *) calloc((unsigned long)n, sizeof(int));  
  red_nf -> w = (int *) calloc((unsigned long)n, sizeof(int)); 
  red_nf -> u = (int *) calloc((unsigned long)n, sizeof(int)); 
  red_nf -> v = (int *) calloc((unsigned long)n, sizeof(int));
  red_nf -> b = (int *) calloc((unsigned long)n, sizeof(int));
  red_nf -> D_red = (int **) calloc((unsigned long)n, sizeof(int*));
  red_nf -> B_red = (int **) calloc((unsigned long)n, sizeof(int*));
  red_nf -> A1 = (int **) calloc((unsigned long)n, sizeof(int*));
  red_nf -> A2 = (int **) calloc((unsigned long)n, sizeof(int*));
  red_nf -> A3 = (int **) calloc((unsigned long)n, sizeof(int*));
  for (r = 0; r < n; ++r){
    red_nf -> D_red[r] = (int *) calloc((unsigned long)n, sizeof(int));
    red_nf -> B_red[r] = (int *) calloc((unsigned long)n, sizeof(int));
    red_nf -> A1[r] = (int *) calloc((unsigned long)n, sizeof(int));
    red_nf -> A2[r] = (int *) calloc((unsigned long)n, sizeof(int));
    red_nf -> A3[r] = (int *) calloc((unsigned long)n, sizeof(int));
  }
  
  printf("\n=======================================================================\n");
  printf("             NORMAL FORMS FOR A %ld-QUBIT STABILIZER CIRCUIT        \n",n);
  printf("=======================================================================\n");  
 
  switch (mode) {
    
  case manual :
    print_instructions_man(n);
    /* reducing a new circuit */
    while (q_or_c[0] == 'c') { 
      input -> len = 0;
      initialize_nf(nf);
      /* starting the induction process :
	 during the while loop, the gate entered by the user is added to the normal form */
      while (1) {
	initialize_red_nf(red_nf, nf);
	compute_red_nf(red_nf, PZX, A_red_D_prod, A_red_B_prod, CNOT_prod, A_red_D_inv, A_red_B_inv, A_aux);
	print_input(input, n);
	print_circuit_nf(nf, CNOT_prod, A_aux);
	print_circuit_red_nf(red_nf, CNOT_prod, A_aux);  
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
      simplify_nf(nf);
      print_normal_form(nf);
      print_CZ_red_normal_form(red_nf);
      printf("\nEnter \'c\' to scan another circuit or \'q\' to quit program.\n--> ");
      sf=scanf("%1s",q_or_c);
      while (getchar() != '\n');
      while (q_or_c[0] != 'c' && q_or_c[0] != 'q') {
	printf("Enter \'c\' to scan another circuit or \'q\' to quit program.\n--> ");
	sf=scanf("%1s",q_or_c);
	while (getchar() != '\n');
      }
    }
    /* end of case manual */
    break;

  case statistics :
    print_introduction_stat();    
    while (q_or_c[0] == 'c') {
      printf("Enter the number of random stabilizer circuits in your sample and the length of the circuits.\n");
      printf("\nSample size (max %d) ? ", MAX_SAMPLE);
      while (scanf("%ld", &sample_size) != 1) { 
	while (getchar() != '\n'); 
	printf("Please, try again.\n");
      }
      while (getchar() != '\n');
      if (sample_size > MAX_SAMPLE || sample_size < 1) {
	printf("Invalid value.\n");
	break;
      }
      printf("Sample size --> %ld\n", sample_size);
      printf("\nCircuit length (max %d) ? ", MAX_LENGTH);
      while(scanf("%ld", &input_len) != 1) { 
	while (getchar() != '\n'); 
	printf("Please, try again.\n");
      }
      while (getchar() != '\n');
      if (input_len > MAX_LENGTH || input_len < 1) {
	printf("Invalid value.\n");
	break;
      }
      printf("Circuit length --> %lu\n", input_len);
      
      sum_nf = 0;
      sum_nf_2q = 0;
      sum_cz_red = 0;
      sum_cz_red_2q = 0;
      sum_input_2q = 0;
      srandom((unsigned int)time(NULL));
      for (s = 0; s < sample_size; ++s) {
	initialize_nf(nf);
	input_2q_len = 0;
	for (t = 0; t < input_len; ++t) {
	  if (random() % 100 < PROBABILITY) { // CNOT chosen
	    i = random() % n;
	    j = random() % (n-1);
	    if (j >= i) ++j;
	    merge_CNOTij_with_nf(i, j, nf, PZX, PZX_prod, CNOT_prod);
	    input_2q_len += 1;
	  } else if (random() % 100 < 50) { // P chosen
	    i = random() % n;
	    merge_Pi_with_nf(i, nf, PZX, PZX_prod, CNOT_prod);
	  } else { // H chosen
	    i = random() % n;
	    merge_Hi_with_nf(i, nf);
	  }
	}
	initialize_red_nf(red_nf, nf);
	simplify_nf(nf);
	compute_red_nf(red_nf, PZX, A_red_D_prod, A_red_B_prod, CNOT_prod, A_red_D_inv, A_red_B_inv, A_aux);
	simplify_red_nf(red_nf);
	count_gate_nf(nf, CNOT_prod, A_aux, &nf_len, &nf_2q_len );
	count_gate_cz_red(red_nf, CNOT_prod, A_aux, &cz_red_len, &cz_red_2q_len );
	sum_input_2q += input_2q_len;
	sum_nf += nf_len;
	sum_nf_2q += nf_2q_len;
	sum_cz_red += cz_red_len;
	sum_cz_red_2q += cz_red_2q_len;
      }
      printf("\n********************************************************************************");
      printf("\n*    STATISTICS FOR A SAMPLE OF %3ld STABILIZER CIRCUITS (%3ld-QUBIT REGISTER)   *", sample_size, n);
      printf("\n********************************************************************************"); 
      printf("\n*                   *       ALL GATES COUNT     *     2-QUBIT GATES COUNT      *");
      printf("\n********************************************************************************"); 
      printf("\n* INPUT CIRCUIT     *     %6ld -->   100%%     *     %6.0lf -->   100%%        *", input_len, (double)sum_input_2q/(double)sample_size);
      printf("\n********************************************************************************"); 
      printf("\n* NF CIRCUIT        *     %6.0lf --> %5.1lf%%     *     %6.0lf --> %5.1lf%%        *", (double)sum_nf/(double)sample_size, (double)sum_nf/(double)sample_size*100./(double)input_len, (double)sum_nf_2q/(double)sample_size, sum_input_2q ? (double)sum_nf_2q/(double)sum_input_2q * 100. : 100.);
      printf("\n********************************************************************************"); 
      printf("\n* CZ-RED NF CIRCUIT *     %6.0lf --> %5.1lf%%     *     %6.0lf --> %5.1lf%%        *",
	     (double)sum_cz_red/(double)sample_size, (double)sum_cz_red/(double)sample_size*100./(double)input_len, (double)sum_cz_red_2q/(double)sample_size, sum_input_2q ? (double)sum_cz_red_2q/(double)sum_input_2q * 100. : 100.);
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
    /* end of case statistics */
    break;
  } // end of switch
  printf("Exiting program.\n");
  
  /* free memory for input circuit */
  free(input);
  
  /* free memory for nf */
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

  /* free memory for red_nf */
  free(red_nf -> a);
  free(red_nf -> d);
  free(red_nf -> w);
  free(red_nf -> u);
  free(red_nf -> v);
  free(red_nf -> b);
  for (r = 0; r < n; ++r) {
    free(red_nf -> D_red[r]);
    free(red_nf -> B_red[r]);
    free(red_nf -> A1[r]);
    free(red_nf -> A2[r]);
    free(red_nf -> A3[r]);
  }
  free(red_nf -> D_red);
  free(red_nf -> B_red);
  free(red_nf -> A1);
  free(red_nf -> A2);
  free(red_nf -> A3);
  free(red_nf);

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
  free(A_red_D_prod -> g);
  free(A_red_D_prod);
  free(A_red_B_prod -> g);
  free(A_red_B_prod);

  /* free memory for matrices A_red_D_inv and A_red_B_inv */
  for (r = 0; r < n; ++r) {
    free(A_red_D_inv[r]);
    free(A_red_B_inv[r]);
  }
  free(A_red_D_inv);
  free(A_red_B_inv);

  /* free memory for matrix A_aux */
  for (r = 0; r < n; ++r) {
    free(A_aux[r]);
  }
  free(A_aux);
}

    
  
    


