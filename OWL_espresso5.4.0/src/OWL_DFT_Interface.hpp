// Define the external Quantum Espresso F90 subroutines


extern "C" {

  void owl_qe_startup(int* lib_comm, int* nim, int* npl, int* nta, int* nbn, int* ndg, const char* infile);
  void run_pwscf_(int *exit_status);
  void owl_qe_stop(int *exit_status);
  void get_natom_ener(int *natom, double *f_etot);
  void get_pos_array(double *pos_array);
  void get_cell_array(double *cell_array);
  void pass_pos_array(double *pos_array);
  void pass_cell_array(double *cell_array);
  void owl_do_pwscf(int *exit_status);
  void owl_stop_run(int *exit_status);

}

