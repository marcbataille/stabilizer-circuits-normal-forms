#include <stdlib.h>
#include <math.h>

#include "gate_struct.h"
#include "aux_ops.h"

/* Multiply the matrix M to the left by the transvection matrix [ij]  */
void left_mult_by_trans(int **M, long i, long j, long n){
  for (long c=0; c < n; ++c) {
    M[i][c] = (M[i][c] + M[j][c]) % 2;
  }
}

/* Multiply the matrix M to the left by a product of CNOT gates  */
void left_mult_by_CNOT_prod(int **M, gate_prod *CNOT_prod, long n) {
  long i, j;
  for (long pos = CNOT_prod -> len - 1; pos >= 0; --pos) {
    i = CNOT_prod -> g[pos].q_i;
    j = CNOT_prod -> g[pos].q_j;
    for (long c = 0; c < n; ++c) {
      M[i][c] = (M[i][c] + M[j][c]) % 2;
    }
  }
}

/* Multiply the matrix M to the right by the transvection matrix [ij]  */
void right_mult_by_trans(int **M, long i, long j, long n){
  for (long r=0; r < n; ++r) {
    M[r][j] = (M[r][j] + M[r][i]) % 2;
  }
}

/* Multiply the matrix M to the right by a product of CNOT gates  */
void right_mult_by_CNOT_prod(int **M, gate_prod *CNOT_prod, long n) {
  long i, j;
  for (long pos = 0; pos < CNOT_prod -> len; ++pos) {
    i = CNOT_prod -> g[pos].q_i;
    j = CNOT_prod -> g[pos].q_j;
    for (long r=0; r < n; ++r) {
      M[r][j] = (M[r][j] + M[r][i]) % 2;
    }
  }
}

/* inverting a CNOT gate product  */
void invert_CNOT_prod(gate_prod *CNOT_prod) {
  long start = 0;
  long end = CNOT_prod -> len -1;
  gate g;
  while (start <= end) {
    if (start == end) {
      break;
    }
    /* swap the gates */
    g = CNOT_prod -> g[end];
    CNOT_prod -> g[end] = CNOT_prod -> g[start];
    CNOT_prod -> g[start] = g;
    ++start;
    --end;
  }
}

/* transpose a CNOT gate product between position start and position end  */
void transpose_CNOT_prod(gate_prod *CNOT_prod, long start, long end) {
  gate g;
  long q;
  while (start <= end) {
    if (start == end) {
      /* swap qubit i and j */
      q = CNOT_prod -> g[end].q_i;
      CNOT_prod -> g[end].q_i = CNOT_prod -> g[end].q_j;
      CNOT_prod -> g[end].q_j = q;
      break;
    }
    /* swap the gates */
    g = CNOT_prod -> g[end];
    CNOT_prod -> g[end] = CNOT_prod -> g[start];
    CNOT_prod -> g[start] = g;
    /* swap qubit i and j in both gates */
    q = CNOT_prod -> g[end].q_i;
    CNOT_prod -> g[end].q_i = CNOT_prod -> g[end].q_j;
    CNOT_prod -> g[end].q_j = q;
    q = CNOT_prod -> g[start].q_i;
    CNOT_prod -> g[start].q_i = CNOT_prod -> g[start].q_j;
    CNOT_prod -> g[start].q_j = q;
    ++start;
    --end;
  }
}
  
/* Multiply the vector vec to the left by the transvection matrix [ij] */
void mult_vec_by_trans(int *vec, long i, long j){
  vec[i] = ( vec[i] + vec[j] ) % 2;
}

/* Multiply the vector vec to the left by a product of CNOT gates */
void mult_vec_by_CNOT_prod(int *vec, gate_prod *CNOT_prod) {
  for (long pos = CNOT_prod -> len - 1; pos >= 0; --pos) {
    vec[CNOT_prod -> g[pos].q_i] = ( vec[CNOT_prod -> g[pos].q_i] + vec[CNOT_prod -> g[pos].q_j] ) % 2;
  }
}

/* Multiply the matrix M to the left by the permutation matrix (ij)  */
void left_mult_by_SWAP(int **M, long i, long j, long n){
  int temp;
  long c;
  for (c=0; c < n; ++c) {
    temp = M[i][c];
    M[i][c] = M[j][c];
    M[j][c] = temp;
  }
}

/* Conjugate the matrix M by the permutation matrix (ij)  */
void conj_by_SWAP(int **M, long i, long j, long n){
  int temp;
  long r, c;
  for (c=0; c < n; ++c) {
    temp = M[i][c];
    M[i][c] = M[j][c];
    M[j][c] = temp;
  }
  for (r=0; r < n; ++r) {
    temp = M[r][i];
    M[r][i] = M[r][j];
    M[r][j] = temp;
  }
}

/* transpose a matrix */
void transpose_matrix(int **M, long n){
  int aux;
  long r, c;
  for (r = 0; r < n; ++r) {
    for (c = r + 1; c < n; ++c) {
      aux = M[r][c];
      M[r][c] = M[c][r];
      M[c][r]= aux;
    }
  }
}

/* Multiply the vector vec to the left by the permutation matrix (ij) */
void mult_vec_by_SWAP(int *vec, long i, long j){
  int temp = vec[i];
  vec[i] = vec[j];
  vec[j] = temp;
}

/* Add the vector source to the vector target */
void vector_add(int *source, int *target, long n) {
  for (long r = 0; r < n; ++r) {
    target[r] = (target[r] + source[r]) %2;
  }
}

/* Scalar product of two vectors modulo 2 */
int scalar_prod(int *u, int *v, long n) {
  int res = 0;
  for (long r = 0; r < n; ++r) {
    res += u[r]*v[r];
  }
  return res % 2;
}

/* Copy the vector source to the vector target*/
void vector_cp(int *source, int *target, long n) {
  for (long r = 0; r < n; ++r) {
    target[r] = source [r];
  }
}

/* Copy the matrix source to the matrix target */ 
void matrix_cp(int **source, int **target, long n) {
  long r, c;
  for (r = 0; r < n; ++r) {
    for (c = 0; c < n; ++c) {
      target[r][c] = source[r][c];
    }
  }
}

/* Conjugate the Pauli form x_u z_v by a CNOT product */
void pauli_conj_CNOT_prod(int *u, int *v, gate_prod *CNOT_prod) {
  gate CNOT;
  long pos;
  for (pos = CNOT_prod -> len - 1; pos >= 0; --pos) {
    CNOT = CNOT_prod -> g[pos];
    mult_vec_by_trans(u, CNOT.q_i, CNOT.q_j);
    mult_vec_by_trans(v, CNOT.q_j, CNOT.q_i);
  }
}

/* Multiply the Pauli form exp(i*k*pi/4) x_u z_v and the Pauli form exp(i*k'*pi/4) x_u' z_v' */
void pauli_mult(int *k_source, int *u_source, int *v_source, int *k_target, int *u_target, int *v_target, long n) {
  *k_target = (*k_target + *k_source + 4 * scalar_prod(u_target, v_source, n)) % 8;
  vector_add(u_source, u_target, n);
  vector_add(v_source, v_target, n);
}

/* Conjugate the Pauli form (exp(i*k*pi/4) x_u z_v by the P_i gate */
void pauli_conj_phase(int *k,int *u, int *v, long i) {
  *k = (*k + 2 * u[i]) % 8;
  v[i] = (v[i] + u[i]) % 2;
}

/* Conjugate the Pauli form (exp(i*k*pi/4) x_u z_v by the Z_B gate */
void pauli_conj_CZ(int *k, int *u, int *v, int **B, long n) {
  long r, c;
  for (r = 0; r < n; ++r) {
    for(c = 0; c < r; ++c) {
      if (B[r][c] == 1) {
	*k = (*k + 4 * u[r] * u[c]) % 8;  
	v[r] = (v[r] + u[c]) % 2;
	v[c] = (v[c] + u[r]) % 2;
      }
    }
  }
}

/* Conjugate the Pauli form x_u z_v by the h gate */
void pauli_conj_h(int *u, int *v, long n) {
  int temp;
  long r;
  for (r = 0; r < n; ++r){
    temp = u[r];
    u[r] = v[r];
    v[r] = temp;
  }
}

void reduce_CZ(int **B, gate_prod *A_red_B_prod, long n) {
  long len = 0;
  long r, c, cc, row, col, pivot;
  for (c = 0; c < n; ++c) {
    pivot = -1;
    for (r = 0; r < n; ++r) {
      if (B[r][c] == 1) {
	pivot = r;
	break;
      }
    }
    if (pivot >= 0) {
      for (r = pivot + 1; r < n; ++r) {
	if (B[r][c] == 1) {
	  for (col = 0; col < n; ++col) {
	    B[r][col] = (B[r][col] + B[pivot][col]) % 2; 
	  }
	  for (row = 0; row < n; ++row) {
	    B[row][r] = (B[row][r] + B[row][pivot]) % 2; 
	  }
	  A_red_B_prod -> g[len].type = CNOT;
	  A_red_B_prod -> g[len].q_i = r;
	  A_red_B_prod -> g[len].q_j = pivot;
	  ++ len;
	}
      }
      for (cc = c + 1; cc < n; ++cc) {
	if (B[pivot][cc] == 1) {
	  for (row = 0; row < n; ++row) {
	    B[row][cc] = (B[row][cc] + B[row][c]) % 2; 
	  }
	  for (col = 0; col < n; ++col) {
	    B[cc][col] = (B[cc][col] + B[c][col]) % 2; 
	  }
	  A_red_B_prod -> g[len].type = CNOT;
	  A_red_B_prod -> g[len].q_i = cc;
	  A_red_B_prod -> g[len].q_j = c;
	  ++ len;
	}
      }
    }
  }
  transpose_CNOT_prod(A_red_B_prod, 0, len - 1); 
  A_red_B_prod -> len = len;
}

/* red_nf is initialized from nf after applying the C-to-NF algorithm (see paper) */
void initialize_red_nf(CZ_red_normal_form *red_nf, normal_form *nf) {
  long r, c;
  long n = red_nf -> n;
  red_nf -> k = nf -> k;
  vector_cp(nf -> a, red_nf -> a, n);
  vector_cp(nf -> w, red_nf -> w, n);
  vector_cp(nf -> d, red_nf -> d, n);
  vector_cp(nf -> b, red_nf -> b, n);
  vector_cp(nf -> u, red_nf -> u, n);
  vector_cp(nf -> v, red_nf -> v, n);
  matrix_cp(nf -> D, red_nf -> D_red, n);
  matrix_cp(nf -> B, red_nf -> B_red, n);
  matrix_cp(nf -> A, red_nf -> A2, n);
  for (r = 0; r < n; ++r) {
    for (c = 0; c < n; ++c) {
      red_nf -> A1[r][c] = 0;
      red_nf -> A3[r][c] = 0;
    }
    red_nf -> A1[r][r] = 1;
    red_nf -> A3[r][r] = 1;
  }
}

void compute_A_inv(int **A_inv, gate_prod *A_prod, long n){
  /* initialize matrix A_inv to identity */
  long r, c;
  long pos;
  for (r = 0; r < n; ++r) {
    for (c = 0; c < n; ++c) {
      A_inv[r][c]=0;
    }
    A_inv[r][r]=1;
  }
  /* multiply A_inv by the inverse of the product of CNOT gates */
  for(pos = 0; pos < A_prod -> len; ++pos) {
    left_mult_by_trans(A_inv, A_prod -> g[pos].q_i, A_prod -> g[pos].q_j, n);
  }
}

void compute_qB_of_A(int **B, int **A, int *vec, long n){
  long c, row, col;
  for (c = 0; c < n; ++c) {
    vec[c] = 0;
    for (row = 0; row < n; ++row) {
      for (col = row + 1; col < n; ++col) {
	vec[c] += B[row][col] * A[row][c] * A[col][c];
      }
    }
    vec[c] = vec[c] % 2;
  }
}

long decompose_GL_lower(int **A, gate_prod *CNOT_prod, long n) {
  long r, c, row, col, width, e, p;
  long m = (long) ceil((0.5 * log((double)n) / log(2)));
  long max = (long) pow(2,(double)m);
  long len = CNOT_prod -> len;
  long pattern[max]; // pattern[p] is the row where the binary pattern equal to p appears for first time.
  for (c = 0; c < n; c = c + m) {
    width = (c + m < n) ? m : n - c;
    for (p = 0; p < max; ++p) {
      pattern[p] = -1;
    }
    for (r = c; r < n; ++r) {
      p = 0;
      e = width;
      for (col= c; col < c + width; ++col) {
	--e;
	p += (long) (pow(2,(double)e) * A[r][col]);  
      }
      if (pattern[p] != -1) {
	if (p != 0) {
	  vector_add(A[pattern[p]], A[r], n);
	  CNOT_prod -> g[len].type = CNOT;
	  CNOT_prod -> g[len].q_i= r;
	  CNOT_prod -> g[len].q_j= pattern[p];
	  ++len;
	}
      } else {
	pattern[p] = r;
      }
    }
    for (col = c; col < c + width; ++ col) { // Gaussian elimination
      if (A[col][col] == 0) { // nead to put 1 on diagonal
	for (row = col + 1 ;row < n; ++row) {
	  if (A[row][col] == 1) { // found a one
	  vector_add(A[row], A[col], n);
	  CNOT_prod -> g[len].type = CNOT;
	  CNOT_prod -> g[len].q_i= col;
	  CNOT_prod -> g[len].q_j= row;
	  ++len;
	  break;
	  }
	}
      }
      for (row = col + 1 ;row < n; ++row) {
	if (A[row][col] == 1) {
	  vector_add(A[col], A[row], n);
	  CNOT_prod -> g[len].type = CNOT;
	  CNOT_prod -> g[len].q_i= row;
	  CNOT_prod -> g[len].q_j= col;
	  ++len;
	}
      }
    }
  }
  CNOT_prod -> len = len;
  return len;
}

void decompose_GL_matrix(int **A, gate_prod *CNOT_prod, long n) {
  CNOT_prod -> len = 0;
  long start, end;
  start = decompose_GL_lower(A, CNOT_prod, n);
  transpose_matrix(A, n);
  end = decompose_GL_lower(A, CNOT_prod, n) - 1;
  transpose_CNOT_prod(CNOT_prod, start, end);
}

void count_gate_nf(normal_form *nf, gate_prod *CNOT_prod, int **A_aux, long *nf_len, long *nf_2q_len ) {
  long n = nf -> n;
  long r, c;
  int ok;
  int a_aux [n];
  int w_aux [n];
  *nf_len = 0;
  *nf_2q_len = 0;
  vector_cp(nf -> a, a_aux, n);
  vector_cp(nf -> w, w_aux, n);
  /* simplify the Hadamard gates for printing */
  for (r = 0; r < n; ++r) {
    if (nf -> a[r] == 1 && nf -> d[r] == 0) {
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
  matrix_cp(nf -> A, A_aux, n);
  decompose_GL_matrix(A_aux, CNOT_prod, n);
  *nf_len += CNOT_prod -> len;
  *nf_2q_len += CNOT_prod -> len;
  for (r = 0; r < n; ++r) {
    for (c = r + 1; c < n; ++c) { 
      if (nf -> B[r][c]) {
	*nf_len += 1;
	*nf_2q_len += 1; 
      }
    }
  }
  for (r = 0; r < n; ++r) {
    if (nf -> b[r]) {
      *nf_len += 1;
    }
  }
  for (r = 0; r < n; ++r) {
    if (nf -> v[r]) {
      *nf_len += 1;
    }
  }
  for (r = 0; r < n; ++r) {
    if (nf -> u[r]) {
      *nf_len += 1;
    }
  }
  for (r = 0; r < n; ++r) {
    if (w_aux[r]) {
      *nf_len += 1;
    }
  }
  for (r = 0; r < n; ++r) {
    for (c = r + 1; c < n; ++c) { 
      if (nf -> D[r][c]) {
	*nf_len += 1;
	*nf_2q_len += 1;
      }
    }
  }
  for (r = 0; r < n; ++r) {
    if (nf -> d[r]) {
      *nf_len += 1;
    }
  }
  for (r = 0; r < n; ++r) {
    if (a_aux[r]) {
      *nf_len += 1;
    }
  }
}


void count_gate_cz_red(CZ_red_normal_form *red_nf, gate_prod *CNOT_prod, int **A_aux, long *cz_red_len, long *cz_red_2q_len ) {
  long r, c;
  long n = red_nf -> n;
  *cz_red_len = 0; 
  *cz_red_2q_len = 0;
  for (r = 0; r < n; ++r) {
    if (red_nf -> b[r]) {
      *cz_red_len += 1;
    }
  }
  matrix_cp(red_nf -> A3, A_aux, n);
  decompose_GL_matrix(A_aux, CNOT_prod, n);
  *cz_red_len += CNOT_prod -> len;
  *cz_red_2q_len += CNOT_prod -> len;
  for (r = 0; r < n; ++r) {
    for (c = r + 1; c < n; ++c) { 
      if (red_nf -> B_red[r][c]) {
	*cz_red_len += 1;
	*cz_red_2q_len += 1;
      }
    }
  }
  matrix_cp(red_nf -> A2, A_aux, n);
  decompose_GL_matrix(A_aux, CNOT_prod, n);
  *cz_red_len += CNOT_prod -> len;
  *cz_red_2q_len += CNOT_prod -> len;
  for (r = 0; r < n; ++r) {
    if (red_nf -> v[r]) {
      *cz_red_len += 1;
    }
  }
  for (r = 0; r < n; ++r) {
    if (red_nf -> u[r]) {
      *cz_red_len += 1;
    }
  }
  for (r = 0; r < n; ++r) {
    if (red_nf -> w[r]) {
      *cz_red_len += 1;
    }
  }
  for (r = 0; r < n; ++r) {
    for (c = r + 1; c < n; ++c) { 
      if (red_nf -> D_red[r][c]) {
	*cz_red_len += 1;
	*cz_red_2q_len += 1;
      }
    }
  }
  matrix_cp(red_nf -> A1, A_aux, n);
  decompose_GL_matrix(A_aux, CNOT_prod, n);
  *cz_red_len += CNOT_prod -> len;
  *cz_red_2q_len += CNOT_prod -> len;
  for (r = 0; r < n; ++r) {
    if (red_nf -> d[r]) {
      *cz_red_len += 1;
    }
  }
  for (r = 0; r < n; ++r) {
    if (red_nf -> a[r]) {
      *cz_red_len += 1;
    }
  }
}
