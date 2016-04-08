// Define the external Quantum Espresso F90 subroutines

extern "C" {
  void wl_qe_startup_(int *my_comm);
  void run_pwscf_(int *exit_status);
  void wl_qe_stop_(int *exit_status);
  void get_natom_ener_(int *natom, double *f_etot);
  void get_array_(double *pos_array, double *cell_array);
  void wl_pass_array_(double *pos_array);
  void wl_do_pwscf_(int *exit_status);
  void wl_stop_run_(int *exit_status);
}


void writeEnergyFile(char fileName[], double energy)
{
  FILE *energy_file;
  energy_file = fopen(fileName, "a");
  fprintf(energy_file, "%f", energy);
  fclose(energy_file);
}

