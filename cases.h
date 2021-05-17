

void initialize_nf(normal_form *nf);

void initialize_gs(graph_state *gs);

void initialize_PZX(PZX_form *PZX);

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

void compute_gs(graph_state *gs, int **B_aux, gate_prod *CNOT_prod);
