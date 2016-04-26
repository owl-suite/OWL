// Define the external Quantum Espresso F90 subroutines

extern "C" {
  void wl_qe_startup_(int *my_comm);
  void run_pwscf_(int *exit_status);
  void wl_qe_stop_(int *exit_status);
  void get_natom_ener_(int *natom, double *f_etot);
  void get_pos_array_(double *pos_array);
  void get_cell_array_(double *cell_array);
  void pass_pos_array_(double *pos_array);
  void pass_cell_array_(double *cell_array);
  void wl_do_pwscf_(int *exit_status);
  void wl_stop_run_(int *exit_status);
}

