#include <stdlib.h>

#include "gate_struct.h"
#include "aux_ops.h"
#include "cases.h"

void initialize_nf(normal_form *nf) {
  long r, c;
  long n = nf -> n;
  nf -> k = 0;
  for (r = 0; r < n; ++r) {
    nf -> a[r] = 1;
    nf -> d[r] = 0;
    nf -> w[r] = 1;
    nf -> u[r] = 0;
    nf -> v[r] = 0;
    nf -> b[r] = 0;
    for (c = 0; c < n; ++c) {
      nf -> D[r][c] = 0;
      nf -> B[r][c] = 0;
      nf -> A[r][c] = 0;
    }
    nf -> A[r][r] = 1;
  }
}

void initialize_PZX(PZX_form *PZX) {
  long r, c;
  long n = PZX -> n;
  for (r = 0; r < n; ++r) {
    PZX -> v[r] = 0;
    PZX -> b[r] = 0;
    for (c = 0; c < n; ++c) {
      PZX -> B[r][c] = 0;
      PZX -> A[r][c] = 0;
    }
    PZX -> A[r][r] = 1;
  }
}

void C_to_PZX(gate_prod *PZX_prod, PZX_form *PZX) {
  long n = PZX -> n;
  long i, j, pos;
  gate g;
  for (pos = PZX_prod -> len - 1; pos >= 0; --pos) {
    g = PZX_prod -> g[pos];
    if (g.type == CZ) {
      i = g.q_i;
      j = g.q_j;
      PZX -> B[i][j] = ((PZX -> B[i][j]) + 1) % 2;
      PZX -> B[j][i] = ((PZX -> B[j][i]) + 1) % 2;
    } else if (g.type == P) {
      i = g.q_i;
      PZX -> v[i] = ((PZX -> v[i]) + (PZX -> b[i])) % 2;
      PZX -> b[i] = ((PZX -> b[i]) + 1) % 2;
    } else { // g.type == CNOT
      i = g.q_i;
      j = g.q_j;
      mult_vec_by_trans(PZX -> v, j, i); 
      PZX -> v[j] = ((PZX -> v[j]) + (PZX -> b[i]) * (PZX -> b[j]) + (PZX -> B[i][j])) % 2;
      left_mult_by_trans(PZX -> B, j, i, n);
      right_mult_by_trans(PZX -> B, i, j, n);
      PZX -> B[i][j] = ((PZX -> B[i][j]) + (PZX -> b[i])) % 2;
      PZX -> B[j][i] = ((PZX -> B[j][i]) + (PZX -> b[i])) % 2;
      mult_vec_by_trans(PZX -> b, j, i);
      left_mult_by_trans(PZX -> A, i, j, n);
    }
  }
}

void merge_Hi_with_nf(long i, normal_form *nf) {
  // case 1 (see paper)
  nf -> a[i] = ((nf -> a[i]) + 1) % 2;
}

void case_21(long i, normal_form *nf) {
  nf -> u[i] = ((nf -> u[i]) + (nf -> d[i])) % 2;
  nf -> d[i] = ((nf -> d[i]) + 1) % 2;
}

void case_221(long i, normal_form *nf, PZX_form *PZX, gate_prod *PZX_prod, gate_prod *CNOT_prod) {
  CNOT_prod -> len = 0;
  PZX_prod -> len = 0;
  long r, c; 
  long n = nf -> n;
  /* Constucting the CNOT product, corresponding to 
     the matrix A'^{-T} of case 2.2.1 */
  for (c = 0; c < n; ++c) {
    if (nf -> D[i][c] == 1) {
      CNOT_prod -> g[CNOT_prod -> len].type = CNOT;
      CNOT_prod -> g[CNOT_prod -> len].q_i = c;
      CNOT_prod -> g[CNOT_prod -> len].q_j = i;
      ++CNOT_prod -> len;
    }
  }
  pauli_conj_CNOT_prod(nf -> u, nf -> v, CNOT_prod);
  /* Constructing the input circuit of the PZX algorithm */
  for (r = 0; r < n; ++r) {
    if (nf -> d[r] == 1) {
      PZX_prod -> g[PZX_prod -> len].type = P;
      PZX_prod -> g[PZX_prod -> len].q_i = r;
      ++PZX_prod -> len;
    }
  }
  for (r = 0; r < n; ++r) {
    for (c = 0; c < r; ++c){
      if (nf -> D[r][c] == 1 && r != i && c != i) {
	PZX_prod -> g[PZX_prod -> len].type = CZ;
	PZX_prod -> g[PZX_prod -> len].q_i = r;
	PZX_prod -> g[PZX_prod -> len].q_j = c;
	++PZX_prod -> len;
      }
    }
  }
  for (c = 0; c < n; ++c) {
    if (nf -> D[i][c] == 1) {
      PZX_prod -> g[PZX_prod -> len].type = CZ;
      PZX_prod -> g[PZX_prod -> len].q_i = i;
      PZX_prod -> g[PZX_prod -> len].q_j = c;
      ++PZX_prod -> len;
      PZX_prod -> g[PZX_prod -> len].type = CNOT;
      PZX_prod -> g[PZX_prod -> len].q_i = i;
      PZX_prod -> g[PZX_prod -> len].q_j = c;
      ++PZX_prod -> len;
      PZX_prod -> g[PZX_prod -> len].type = P;
      PZX_prod -> g[PZX_prod -> len].q_i = c;
      ++PZX_prod -> len;
    }
  }
  initialize_PZX(PZX);
  C_to_PZX(PZX_prod, PZX);
  vector_cp(PZX -> b, nf -> d, n);
  matrix_cp(PZX -> B, nf -> D, n);
  vector_add(PZX -> v, nf -> u, n);
  initialize_PZX(PZX);
  vector_cp(nf -> b, PZX -> b, n);
  matrix_cp(nf -> B, PZX -> B, n);
  matrix_cp(nf -> A, PZX -> A, n);
  C_to_PZX(CNOT_prod, PZX);
  vector_add(PZX -> v, nf -> v, n);
  vector_cp(PZX -> b, nf -> b, n);
  matrix_cp(PZX -> B, nf -> B, n);
  matrix_cp(PZX -> A, nf -> A, n);
}

void case_222(long i, normal_form *nf, PZX_form *PZX, gate_prod *PZX_prod, gate_prod *CNOT_prod) {
  long r; 
  long n = nf -> n;
  int u_aux[n];
  int v_aux[n];
  for (r = 0; r < n; ++r) {
    u_aux[r] = 0;
    v_aux[r] = 0;
  }
  int k_aux = 0;
  nf -> k = ((nf -> k) + 1) % 8;
  nf -> a[i] = ((nf -> a[i]) + 1) % 2;
  nf -> d[i] = ((nf -> d[i]) + 1) % 2;
  /* Conjugating the x_i gate by P_i h Z_D and merging the result with the pauli part of nf */
  u_aux[i] = 1;
  pauli_conj_CZ(&k_aux, u_aux, v_aux, nf -> D, n);
  pauli_conj_h(u_aux, v_aux, n);
  pauli_conj_phase(&k_aux, u_aux, v_aux, i);
  pauli_mult(&k_aux, u_aux, v_aux, &(nf -> k), nf -> u, nf -> v, n);
  case_221(i, nf, PZX, PZX_prod, CNOT_prod);
}

void case_22(long i, normal_form *nf, PZX_form *PZX, gate_prod *PZX_prod, gate_prod *CNOT_prod) {
  nf -> v[i] = ((nf -> v[i]) + (nf -> u[i]) + (nf -> b[i])) % 2;
  nf -> b[i] = ((nf -> b[i]) + 1) % 2;
  nf -> k = ((nf -> k) + 2 * (nf -> u[i])) % 8;  
  if ((nf -> d[i]) == 0) { 
    case_221(i, nf, PZX, PZX_prod, CNOT_prod);
  } else { 
    case_222(i, nf, PZX, PZX_prod, CNOT_prod);
  }
}
 
void merge_Pi_with_nf(long i, normal_form *nf, PZX_form *PZX, gate_prod *PZX_prod, gate_prod *CNOT_prod ) {
  // case 2
  if ((nf -> a[i]) == 0) { 
    case_21(i, nf);
  } else { // case 2.2 (see paper)
    case_22(i, nf, PZX, PZX_prod, CNOT_prod);
  }
}
 
void case_31(long i, long j, normal_form *nf) {
  long n = nf -> n;
  /* conjugating x_u*z_v by X[j,i] */
  mult_vec_by_trans(nf -> u, j, i);
  mult_vec_by_trans(nf -> v, i, j);
  /* conjugating P_d*Z_D by X[i,j] */ 
  nf -> u[j] = ((nf -> u[j]) + (nf -> d[i]) * (nf -> d[j]) + (nf -> D[i][j])) % 2;
  left_mult_by_trans(nf -> D, j, i, n);
  right_mult_by_trans(nf -> D, i, j, n);
  nf -> D[i][j] = ((nf -> D[i][j]) + (nf -> d[i])) % 2;
  nf -> D[j][i] = ((nf -> D[j][i]) + (nf -> d[i])) % 2;
  mult_vec_by_trans(nf -> d, j, i);
  /* conjugating P_b*Z_B by X[j,i] */
  nf -> v[i] = ((nf -> v[i]) + (nf -> b[j]) * (nf -> b[i]) + (nf -> B[j][i])) % 2;
  left_mult_by_trans(nf -> B, i, j, n);
  right_mult_by_trans(nf -> B, j, i, n);
  nf -> B[i][j] = ((nf -> B[i][j]) + (nf -> b[j])) % 2;
  nf -> B[j][i] = ((nf -> B[j][i]) + (nf -> b[j])) % 2;
  mult_vec_by_trans(nf -> b, i, j);
  /* merging X_[j,i] with X_A */
  left_mult_by_trans(nf -> A, j, i, n);
}

void case_341(long i, long j, normal_form *nf, PZX_form *PZX, gate_prod *PZX_prod, gate_prod *CNOT_prod) {
  CNOT_prod -> len = 0;
  PZX_prod -> len = 0;
  long r, c; 
  long n = nf -> n;
  int d_aux[n];

  nf -> k =((nf -> k) + 4 * (nf -> u[i]) * (nf -> u[j])) % 8;
  nf -> v[i] = ((nf -> v[i]) + (nf -> u[j])) % 2;
  nf -> v[j] = ((nf -> v[j]) + (nf -> u[i])) % 2;
  nf -> B[i][j] = ((nf -> B[i][j]) + 1) % 2;
  nf -> B[j][i] = ((nf -> B[j][i]) + 1) % 2;

  /* Constucting the CNOT product, corresponding to 
     the matrix A'^{-T} of case 3.4.1 */
  for (c = 0; c < n; ++c) {
    if (nf -> D[i][c] == 1) {
      CNOT_prod -> g[CNOT_prod -> len].type = CNOT;
      CNOT_prod -> g[CNOT_prod -> len].q_i = c;
      CNOT_prod -> g[CNOT_prod -> len].q_j = j;
      ++CNOT_prod -> len;
    }
  }
  for (c = 0; c < n; ++c) {
    if (nf -> D[j][c] == 1) {
      CNOT_prod -> g[CNOT_prod -> len].type = CNOT;
      CNOT_prod -> g[CNOT_prod -> len].q_i = c;
      CNOT_prod -> g[CNOT_prod -> len].q_j = i;
      ++CNOT_prod -> len;
    }
  }
  pauli_conj_CNOT_prod(nf -> u, nf -> v, CNOT_prod);
        
  /* Constructing the input circuit of the PZX algorithm */
  for (r = 0; r < n; ++r) {
    for (c = 0; c < r; ++c){
      if (nf -> D[r][c] == 1 && r != i && c != i && r != j && c != j) {
	PZX_prod -> g[PZX_prod -> len].type = CZ;
	PZX_prod -> g[PZX_prod -> len].q_i = r;
	PZX_prod -> g[PZX_prod -> len].q_j = c;
	++PZX_prod -> len;
      }
    }
  }

  for (c = 0; c < n; ++c) {
    if (nf -> D[i][c] == 1) {
      PZX_prod -> g[PZX_prod -> len].type = CZ;
      PZX_prod -> g[PZX_prod -> len].q_i = i;
      PZX_prod -> g[PZX_prod -> len].q_j = c;
      ++PZX_prod -> len;
      PZX_prod -> g[PZX_prod -> len].type = CNOT;
      PZX_prod -> g[PZX_prod -> len].q_i = j;
      PZX_prod -> g[PZX_prod -> len].q_j = c;
      ++PZX_prod -> len;
    }
  }
  for (c = 0; c < n; ++c) {
    if (nf -> D[j][c] == 1) {
      PZX_prod -> g[PZX_prod -> len].type = CZ;
      PZX_prod -> g[PZX_prod -> len].q_i = j;
      PZX_prod -> g[PZX_prod -> len].q_j = c;
      ++PZX_prod -> len;
      PZX_prod -> g[PZX_prod -> len].type = CNOT;
      PZX_prod -> g[PZX_prod -> len].q_i = i;
      PZX_prod -> g[PZX_prod -> len].q_j = c;
      ++PZX_prod -> len;
    }
  }
  initialize_PZX(PZX);
  C_to_PZX(PZX_prod, PZX);
  matrix_cp(PZX -> B, nf -> D, n);
  vector_add(PZX -> v, nf -> u, n);
  initialize_PZX(PZX);
  vector_cp(nf -> b, PZX -> b, n);
  matrix_cp(nf -> B, PZX -> B, n);
  matrix_cp(nf -> A, PZX -> A, n);
  C_to_PZX(CNOT_prod, PZX);
  vector_add(PZX -> v, nf -> v, n);
  vector_cp(PZX -> b, nf -> b, n);
  matrix_cp(PZX -> B, nf -> B, n);
  matrix_cp(PZX -> A, nf -> A, n);
  if (nf -> d[i] == 0 && nf -> d[j] == 0) { // case 3.4.1.1
    return;
  } else if (nf -> d[i] == 0 && nf -> d[j] == 1) { // case 3.4.1.2
    case_31(i, j, nf);
    case_22(i, nf, PZX, PZX_prod, CNOT_prod);
  } else if (nf -> d[i] == 1 && nf -> d[j] == 0) { // case 3.4.1.3
    case_31(j, i, nf);
    case_22(j, nf, PZX, PZX_prod, CNOT_prod);
  } else { // case 3.4.1.4 : (nf -> d[i] == 1 && nf -> d[j] == 1)
    /* copy nf -> d to d_aux and replace d by the null vector. */
    vector_cp(nf -> d, d_aux, n);
    for (r = 0; r < n; ++r) {
      nf -> d[r] = 0;
    }
    case_31(j, i, nf);
    case_22(j, nf, PZX, PZX_prod, CNOT_prod);
    for (r = 0; r < n; ++r) {
      if (d_aux[r] == 1) {
	case_21(r, nf);
      }
    }
    case_31(i, j, nf);
    case_22(i, nf, PZX, PZX_prod, CNOT_prod);
  }
}

void merge_CNOTij_with_nf(long i, long j, normal_form *nf, PZX_form *PZX, gate_prod *PZX_prod, gate_prod *CNOT_prod) {
  long n = nf -> n;
  if (nf -> a[i] == 0 && nf -> a[j] == 0) { // case 3.1 
    case_31(i, j, nf);
  } else if (nf -> a[i] == 1 && nf -> a[j] == 1) { // case 3.2
    case_31(j, i, nf);
  } else if (nf -> a[i] == 1 && nf -> a[j] == 0) { // case 3.3
    nf -> D[i][j] = ((nf -> D[i][j]) + 1) %2;
    nf -> D[j][i] = ((nf -> D[j][i]) + 1) %2;
  } else { // case 3.4 : (nf -> a[i] == 0 && nf -> a[j] == 1) 
    if (nf -> D[i][j] == 0) { // case 3.4.1
      case_341(i, j, nf, PZX, PZX_prod, CNOT_prod);
    } else { // case 3.4.2 : nf -> D[i][j] == 1)
      nf -> a[i] = ((nf -> a[i]) + 1) % 2;
      nf -> a[j] = ((nf -> a[j]) + 1) % 2;
      nf -> D[i][j] = 0;
      nf -> D[j][i] = 0;
      mult_vec_by_SWAP(nf -> d, i, j);
      conj_by_SWAP(nf -> D, i, j, n);
      mult_vec_by_SWAP(nf -> u, i, j);
      mult_vec_by_SWAP(nf -> v, i, j);
      mult_vec_by_SWAP(nf -> b, i, j);
      conj_by_SWAP(nf -> B, i, j, n);
      left_mult_by_SWAP(nf -> A, i, j, n);
      case_341(i, j, nf, PZX, PZX_prod, CNOT_prod);
    }
  }
}

void merge_SWAPij_with_nf(long i, long j, normal_form *nf) {
  long n = nf -> n;
  mult_vec_by_SWAP(nf -> a, i, j);
  mult_vec_by_SWAP(nf -> d, i, j);
  conj_by_SWAP(nf -> D, i, j, n);
  mult_vec_by_SWAP(nf -> u, i, j);
  mult_vec_by_SWAP(nf -> v, i, j);
  mult_vec_by_SWAP(nf -> b, i, j);
  conj_by_SWAP(nf -> B, i, j, n);
  left_mult_by_SWAP(nf -> A, i, j, n);
}

/* Zij = Hi Xij Hi */
void merge_CZij_with_nf(long i, long j, normal_form *nf, PZX_form *PZX, gate_prod *PZX_prod, gate_prod *CNOT_prod) {
  merge_Hi_with_nf(i, nf);
  merge_CNOTij_with_nf(i, j, nf, PZX, PZX_prod, CNOT_prod);
  merge_Hi_with_nf(i, nf);
}

/* Zi = Pi^2 */
void merge_Zi_with_nf(long i, normal_form *nf, PZX_form *PZX, gate_prod *PZX_prod, gate_prod *CNOT_prod) {
  merge_Pi_with_nf(i, nf, PZX, PZX_prod, CNOT_prod );
  merge_Pi_with_nf(i, nf, PZX, PZX_prod, CNOT_prod );
}

/* X_i = H_i Z_i H_i */
void merge_Xi_with_nf(long i, normal_form *nf, PZX_form *PZX, gate_prod *PZX_prod, gate_prod *CNOT_prod) {
  merge_Hi_with_nf(i, nf);
  merge_Zi_with_nf(i, nf, PZX, PZX_prod, CNOT_prod);
  merge_Hi_with_nf(i, nf);
}

/* Y=iXZ */
void merge_Yi_with_nf(long i, normal_form *nf, PZX_form *PZX, gate_prod *PZX_prod, gate_prod *CNOT_prod) {
  merge_Zi_with_nf(i, nf, PZX, PZX_prod, CNOT_prod); 
  merge_Xi_with_nf(i, nf, PZX, PZX_prod, CNOT_prod); 
  nf -> k = ((nf -> k) + 2) % 8; 
} 

void simplify_nf(normal_form *nf) {
  long n = nf -> n;
  long r, c;
  int ok;
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
	nf -> a[r] = 0;
	nf -> w[r] = 0;
      }
    }
  }
}

void simplify_red_nf(CZ_red_normal_form *red_nf) {
  long n = red_nf -> n;
  long r, c;
  int ok;
  for (r = 0; r < n; ++r) {
    if (red_nf -> a[r] == 1  && red_nf -> d[r] == 0 && red_nf -> A1[r][r] == 1) {
      ok = 1;
      for(c = 0; c < n; ++c) {
	if (c != r && (red_nf -> D_red[r][c] != 0 || red_nf -> A1[r][c] != 0 || red_nf -> A1[c][r] != 0)) {
	  ok = 0;
	  break;
	}
      }
      if (ok) {
	red_nf -> a[r] = 0;
	red_nf -> w[r] = 0;
      }
    }
  }
}

void compute_red_nf(CZ_red_normal_form *red_nf, PZX_form *PZX, gate_prod *A_red_D_prod, gate_prod *A_red_B_prod, gate_prod *CNOT_prod, int **A_red_D_inv, int **A_red_B_inv, int **A_aux) {
  long n = red_nf -> n;
  long pos;
  int vec[n];
  /* computing w', b' and B' (see paper) */
  initialize_PZX(PZX);
  vector_cp(red_nf -> b, PZX -> b, n);
  matrix_cp(red_nf -> B_red, PZX -> B, n); // B_red is matrix B of the normal form
  matrix_cp(red_nf -> A2, PZX -> A, n); // A2 is matrix A of the normal form
  matrix_cp(red_nf -> A2, A_aux, n);
  decompose_GL_matrix(A_aux, CNOT_prod, n); // now CNOT_prod is A
  invert_CNOT_prod(CNOT_prod); // now CNOT_prod is A^{-1}
  C_to_PZX(CNOT_prod, PZX);
  vector_cp(PZX -> b, red_nf -> b, n);
  matrix_cp(PZX -> B, red_nf -> B_red, n);
  /* computing v' (see paper) */
  transpose_CNOT_prod(CNOT_prod, 0, CNOT_prod -> len - 1); // now CNOT_prod is A^{-T}
  vector_cp(PZX -> v, vec, n);
  mult_vec_by_CNOT_prod(vec, CNOT_prod); // vec is w' in the paper
  vector_add(vec, red_nf -> v, n);
  reduce_CZ(red_nf -> B_red, A_red_B_prod, n);
  compute_A_inv(A_red_B_inv, A_red_B_prod, n);    
  compute_qB_of_A(red_nf -> B_red, A_red_B_inv, vec, n);
  mult_vec_by_CNOT_prod(vec, CNOT_prod); // vec is now A^{-T}*q_B_red(A_red_B^{-1})
  vector_add(vec, red_nf -> v, n);
  /* computing v from v' */
  reduce_CZ(red_nf -> D_red, A_red_D_prod, n);
  for (pos = 0; pos < A_red_D_prod -> len; ++pos) {
    mult_vec_by_trans(red_nf -> v, A_red_D_prod -> g[pos].q_i, A_red_D_prod -> g[pos].q_j);
  }
  /* computing u' (see paper) */
  compute_A_inv(A_red_D_inv, A_red_D_prod, n);
  compute_qB_of_A(red_nf -> D_red, A_red_D_inv, vec, n); // vec is now q_D_red(A_red_D^{-1})
  vector_add(vec, red_nf -> u, n);
  /* computing u from u' */
  for (pos = 0; pos < A_red_D_prod -> len; ++pos) {
    mult_vec_by_trans(red_nf -> u, A_red_D_prod -> g[pos].q_j, A_red_D_prod -> g[pos].q_i);
  }
  /* computing A1 */
  for (pos = A_red_D_prod -> len - 1; pos >= 0; -- pos) {
    left_mult_by_trans(red_nf -> A1, A_red_D_prod -> g[pos].q_i, A_red_D_prod -> g[pos].q_j, n);
  }
  /* computing A2 */
  for (pos = 0; pos < A_red_B_prod -> len; ++pos) {
    right_mult_by_trans(red_nf -> A2, A_red_B_prod -> g[pos].q_i, A_red_B_prod -> g[pos].q_j, n);
  }
  for (pos = 0; pos < A_red_D_prod -> len; ++pos) {
    left_mult_by_trans(red_nf -> A2, A_red_D_prod -> g[pos].q_j, A_red_D_prod -> g[pos].q_i, n);
  }
  /* computing A3 */
  for (pos = 0; pos < A_red_B_prod -> len; ++pos) {
    left_mult_by_trans(red_nf -> A3, A_red_B_prod -> g[pos].q_i, A_red_B_prod -> g[pos].q_j, n);
  }
  simplify_red_nf(red_nf);
}

