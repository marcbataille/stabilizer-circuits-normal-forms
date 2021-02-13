

void left_mult_by_trans(int **M, long i, long j, long n);

void left_mult_by_CNOT_prod(int **M, gate_prod *CNOT_prod, long n) ;

void right_mult_by_trans(int **M, long i, long j, long n);

void right_mult_by_CNOT_prod(int **M, gate_prod *CNOT_prod, long n);

void invert_CNOT_prod(gate_prod *CNOT_prod) ;

void transpose_CNOT_prod(gate_prod *CNOT_prod, long start, long end);

void mult_vec_by_trans(int *vec, long i, long j);

void mult_vec_by_CNOT_prod(int *vec, gate_prod *CNOT_prod);

void left_mult_by_SWAP(int **M, long i, long j, long n);

void conj_by_SWAP(int **M, long i, long j, long n);

void transpose_matrix(int **M, long n);

void mult_vec_by_SWAP(int *vec, long i, long j);

void vector_add(int *source, int *target, long n);

int scalar_prod(int *u, int *v, long n);

void vector_cp(int *source, int *target, long n);

void matrix_cp(int **source, int **target, long n);

void pauli_conj_CNOT_prod(int *u, int *v, gate_prod *CNOT_prod);

void pauli_mult(int *k_source, int *u_source, int *v_source, int *k_target, int *u_target, int *v_target, long n);

void pauli_conj_phase(int *k,int *u, int *v, long i);

void pauli_conj_CZ(int *k, int *u, int *v, int **B, long n);

void pauli_conj_h(int *u, int *v, long n);

void count_gate_cz_red(CZ_red_normal_form *red_nf, gate_prod *CNOT_prod, int **A_aux, long *cz_red_len, long *cz_red_2q_len );

void count_gate_nf(normal_form *nf, gate_prod *CNOT_prod, int **A_aux, long *nf_len, long *nf_2q_len );

long decompose_GL_lower(int **A, gate_prod *CNOT_prod, long n);

void decompose_GL_matrix(int **A, gate_prod *CNOT_prod, long n);

void compute_qB_of_A(int **B, int **A, int *vec, long n);

void compute_A_inv(int **A_inv, gate_prod *A_prod, long n);

void reduce_CZ(int **B, gate_prod *A_red_B_prod, long n);
