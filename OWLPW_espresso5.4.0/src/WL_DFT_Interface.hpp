// Define the external Quantum Espresso F90 subroutines

extern "C" {

#if defined(__ICC) || defined(__INTEL_COMPILER)
  /* Intel ICC/ICPC */
  void command_line_options_mp_get_command_line_(const char* command_line, size_t len);
  //void command_line_options_mp_set_command_line_(int nimage, int npool, int ntg, int nband, int ndiag);

#elif (defined(__GNUC__) || defined(__GNUG__)) //&& !(defined(__clang__) || defined(__INTEL_COMPILER))
  /* GNU GCC/G++ */
  void __command_line_options_MOD_get_command_line(const char* command_line, size_t len);
  //void __command_line_options_MOD_set_command_line(int nimage, int npool, int ntg, int nband, int ndiag);

#endif

  void wl_qe_startup_(int* lib_comm, int* nim, int* npl, int* nta, int* nbn,int* ndg);
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

