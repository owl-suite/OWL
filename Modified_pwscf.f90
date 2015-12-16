!----------------------------------------------------------------------------
PROGRAM pwscf
  !----------------------------------------------------------------------------
  !
  ! ... Main program calling one instance of Plane Wave Self-Consistent Field code
  ! ... Modified to pass the energy and atomic positions to Wang-Landau algorithm
``! ... "Bridge" code between Fortran 90 and C++
  !
  IMPLICIT NONE

  INTEGER :: exit_status
  !
  CALL getEnergy_startup ( )         ! Set up the QE calculations
  CALL getEnergy_run (exit_status)   ! Run the QE calculations
  CALL getEnergy_get ( )             ! Get the energy and atomic positions
  !
  STOP
  !
END PROGRAM pwscf
