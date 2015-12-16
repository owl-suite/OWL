!----------------------------------------------------------------------------
SUBROUTINE getEnergy_startup ( )
  !----------------------------------------------------------------------------
  !
  ! ... Set up the plane wave self-consistent calculations for Quantum Espresso
  !
  USE environment,       ONLY : environment_start
  USE mp_global,         ONLY : mp_startup
  USE read_input,        ONLY : read_input_file
  USE command_line_options, ONLY: input_file_
  !
  IMPLICIT NONE

!  CHARACTER(LEN=256) :: input_file_ = ' '
  !
  CALL mp_startup ( )
  CALL environment_start ('PWSCF')
  !
  CALL read_input_file ('PW', input_file_ )
  !
END SUBROUTINE getEnergy_startup
  !
!----------------------------------------------------------------------------
SUBROUTINE getEnergy_run (exit_status)
  !----------------------------------------------------------------------------
  !
  ! ... Perform the actual calculations using the Quantum Espresso codes
  !
  IMPLICIT NONE

  INTEGER :: exit_status
  !
  CALL run_pwscf  ( exit_status )
  !
  CALL stop_run( exit_status )
  CALL do_stop( exit_status )
  !
END SUBROUTINE getEnergy_run
  !
!----------------------------------------------------------------------------
SUBROUTINE getEnergy_get ( )
  !----------------------------------------------------------------------------
  !
  ! ... Assign the pointers to the energy and position arrays
  !
  USE ener,       ONLY : etot
  !
  IMPLICIT NONE

  real, dimension(:,:), allocatable :: e_array
  !

 allocate (e_array(1,1))
  e_array(1,1) = etot       ! Pass the total energy to the array

END SUBROUTINE getEnergy_get

!----------------------------------------------------------------------------
