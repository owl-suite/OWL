// Define the external Quantum Espresso F90 subroutines

extern "C" {
  //void __command_line_options_MOD_get_command_line(const char* command_line, size_t len);
  void command_line_options_mp_get_command_line_(const char* command_line, size_t len);
  void owl_qe_startup_(int *my_comm);
  void run_pwscf_(int *exit_status);
  void owl_qe_stop_(int *exit_status);
  void get_natom_ener_(int *natom, double *f_etot);
  void get_pos_array_(double *pos_array);
  void get_cell_array_(double *cell_array);
  void pass_pos_array_(double *pos_array);
  void pass_cell_array_(double *cell_array);
  void owl_do_pwscf_(int *exit_status);
  void owl_stop_run_(int *exit_status);
}

