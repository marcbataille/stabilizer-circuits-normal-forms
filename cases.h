
void initialize_PZX(PZX_form *PZX);

void initialize_nf(normal_form *nf);

void initialize_red_nf(CZ_red_normal_form *red_nf, normal_form *nf);

void C_to_PZX(gate_prod *PZX_prod, PZX_form *PZX);

void merge_Hi_with_nf(long i, normal_form *nf);

void case_21(long i, normal_form *nf);

void case_221(long i, normal_form *nf, PZX_form *PZX, gate_prod *PZX_prod, gate_prod *CNOT_prod);

void case_222(long i, normal_form *nf, PZX_form *PZX, gate_prod *PZX_prod, gate_prod *CNOT_prod);

void case_22(long i, normal_form *nf, PZX_form *PZX, gate_prod *PZX_prod, gate_prod *CNOT_prod);
  
void merge_Pi_with_nf(long i, normal_form *nf, PZX_form *PZX, gate_prod *PZX_prod, gate_prod *CNOT_prod );

void case_31(long i, long j, normal_form *nf);

void case_341(long i, long j, normal_form *nf, PZX_form *PZX, gate_prod *PZX_prod, gate_prod *CNOT_prod);

void merge_CNOTij_with_nf(long i, long j, normal_form *nf, PZX_form *PZX, gate_prod *PZX_prod, gate_prod *CNOT_prod);

void merge_SWAPij_with_nf(long i, long j, normal_form *nf);

void merge_CZij_with_nf(long i, long j, normal_form *nf, PZX_form *PZX, gate_prod *PZX_prod, gate_prod *CNOT_prod);

void merge_Zi_with_nf(long i, normal_form *nf, PZX_form *PZX, gate_prod *PZX_prod, gate_prod *CNOT_prod);

void merge_Xi_with_nf(long i, normal_form *nf, PZX_form *PZX, gate_prod *PZX_prod, gate_prod *CNOT_prod);

void merge_Yi_with_nf(long i, normal_form *nf, PZX_form *PZX, gate_prod *PZX_prod, gate_prod *CNOT_prod);

void simplify_nf(normal_form *nf);

void simplify_red_nf(CZ_red_normal_form *red_nf);

void compute_red_nf(CZ_red_normal_form *red_nf, PZX_form *PZX, gate_prod *A_red_D_prod, gate_prod *A_red_B_prod, gate_prod *CNOT_prod, int **A_red_D_inv, int **A_red_B_inv, int **A_aux);
