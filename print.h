
void print_help(long n);

void print_error_cmd(void);

void print_error_qb(long n);
      
void print_instructions_man(long n);

void print_introduction_stat(void);

void print_matrix(int **M, long n);
  
void print_vector(int *vec, long n);

void print_normal_form(normal_form *nf);

void print_CZ_red_normal_form(CZ_red_normal_form *red_nf);
    
void print_input(gate_prod  *input, long n);
    
void print_circuit_nf(normal_form *nf, gate_prod *CNOT_prod, int **A_aux);

void print_circuit_red_nf(CZ_red_normal_form *red_nf, gate_prod *CNOT_prod, int **A_aux);

void print_stabilizer_state(int *u, int *v, int *w, int **B, long n);

void print_graph_state(CZ_red_normal_form *red_nf, gate_prod *A_red_D_prod, int **A_red_D_inv);
