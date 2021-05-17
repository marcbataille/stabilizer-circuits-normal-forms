
void print_help(long n);

void print_error_cmd(void);

void print_error_qb(long n);
      
void print_instructions_man(long n);

void print_information_stat();

void print_matrix(int **M, long n);
  
void print_vector(int *vec, long n);

void print_normal_form(normal_form *nf);

void print_input(gate_prod  *input, long n);
    
void print_circuit_nf(normal_form *nf, gate_prod *CNOT_prod, int **A_aux);

void print_stabilizer_state(int *u, int *v, int *w, int **B, long n);

void print_graph_state(graph_state *gs, gate_prod *CNOT_prod, int **A_aux);
