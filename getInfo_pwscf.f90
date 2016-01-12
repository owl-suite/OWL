!
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt.
!
!----------------------------------------------------------------------------
SUBROUTINE run_qe_startup ( )
  !----------------------------------------------------------------------------
  !
  ! ... Set up the PWscf calculation for Quantum Espresso
  !
  USE environment,       ONLY : environment_start
  USE mp_global,         ONLY : mp_startup
  USE read_input,        ONLY : read_input_file
  USE command_line_options, ONLY: input_file_
  !
  IMPLICIT NONE
  !
  CALL mp_startup ( )
  CALL environment_start ('PWSCF')
  !
  CALL read_input_file ('PW', input_file_ )
  !
END SUBROUTINE run_qe_startup
  !
!----------------------------------------------------------------------------
SUBROUTINE run_qe_run (exit_status)
  !----------------------------------------------------------------------------
  !
  ! ... Perform the PWscf calculation using the QE codes
  !
  IMPLICIT NONE
  INTEGER :: exit_status
  !
  CALL run_pwscf  ( exit_status )
  !
END SUBROUTINE run_qe_run
  !
!----------------------------------------------------------------------------
SUBROUTINE run_qe_stop (exit_status)
  !----------------------------------------------------------------------------
  !
  ! ... Finish the PWscf calculation for Quantum Espresso
  !
  IMPLICIT NONE
  INTEGER :: exit_status
  !
  CALL stop_run( exit_status )
  CALL do_stop( exit_status )
  !
END SUBROUTINE run_qe_stop
  !
!----------------------------------------------------------------------------
SUBROUTINE get_natom_ener (natom, f_etot)
  !----------------------------------------------------------------------------
  !
  ! ... Extract a total number of atoms and total energy from the QE modules
  !
  USE ener,          ONLY : etot
  USE ions_base,     ONLY : nat
  !
  IMPLICIT NONE
  INTEGER :: natom
  REAL(KIND=8) :: f_etot
  !
  natom = nat     ! ... a total number of atoms in system
  f_etot = etot   ! ... total energy of system (in eV/unit cell)
  !
END SUBROUTINE get_natom_ener
  !
!----------------------------------------------------------------------------
SUBROUTINE get_pos_array (pos_array)
  !----------------------------------------------------------------------------
  !
  ! ... Extract the position array from the Quantum Espresso modules
  !
  USE constants,     ONLY : bohr_radius_angs
  USE cell_base,     ONLY : alat
  USE ions_base,     ONLY : tau, nat
  !
  IMPLICIT NONE
  REAL(KIND=8) :: pos_array(3,nat)
  !
  pos_array(:,:) = tau(:,:) * alat * bohr_radius_angs   ! ... (in angstrom)
  !
END SUBROUTINE get_pos_array
  !
!----------------------------------------------------------------------------
SUBROUTINE get_cell_array (cell_array)
  !----------------------------------------------------------------------------
  !
  ! ... Extract the cell vector array from the Quantum Espresso modules
  !
  USE constants,     ONLY : bohr_radius_angs
  USE cell_base,     ONLY : alat, at
  !
  IMPLICIT NONE
  REAL(KIND=8) :: cell_array(3,3)
  !
  cell_array(:,:) = at(:,:) * alat * bohr_radius_angs   ! ... (in angstrom)
  !
END SUBROUTINE get_cell_array
  !
!----------------------------------------------------------------------------
