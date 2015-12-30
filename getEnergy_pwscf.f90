!
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt.
!
!----------------------------------------------------------------------------
SUBROUTINE getenergy_startup ( )
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

!  CHARACTER(LEN=256) :: input_file_ = ' '
  !
  CALL mp_startup ( )
  CALL environment_start ('PWSCF')
  !
  CALL read_input_file ('PW', input_file_ )
  !
END SUBROUTINE getenergy_startup
  !
!----------------------------------------------------------------------------
SUBROUTINE getenergy_run (exit_status)
  !----------------------------------------------------------------------------
  !
  ! ... Perform the PWscf calculation using the Quantum Espresso codes
  !
  USE ener,          ONLY : etot
  USE ions_base,     ONLY : nat, tau, if_pos, atm, ityp
  USE constants,     ONLY : bohr_radius_angs
  USE cell_base,     ONLY : alat, at
  !
  IMPLICIT NONE

  REAL*8,  ALLOCATABLE :: e_array(:,:)
  REAL*8,  ALLOCATABLE :: pos_array(:,:)
  REAL*8,  ALLOCATABLE :: cell_array(:,:)
  INTEGER             :: na
  INTEGER :: exit_status
  !
  CALL run_pwscf  ( exit_status )
  !
  ALLOCATE (e_array(1,1))         ! Pass the total energy to the array
  e_array(1,1) = etot
  !
  ALLOCATE (pos_array(3,nat))     ! Pass the atomic positions to the array
  pos_array(:,:) = tau(:,:)
  !
  ALLOCATE (cell_array(3,3))      ! Pass the cell vector to the array
  cell_array(:,:) = at(:,:) * alat * bohr_radius_angs
!
!  WRITE(1000000,*) e_array(1,1)    ! Write the energy array
!
!  DO na = 1, nat                   ! Write the position array
!
!    IF ( ANY( if_pos(:,na) == 0 ) ) THEN
!        WRITE(1210000,'(A3,3X,3F14.9,1X,3i4)') &
!                        atm(ityp(na)), pos_array(:,na), if_pos(:,na)
!     ELSE
!        WRITE(1210000,'(A3,3X,3F14.9)') &
!                        atm(ityp(na)), pos_array(:,na)
!     END IF
!
!   END DO
!
! WRITE(1220000,*) cell_array(:,:)    ! Write the cell vector array
!
  CALL stop_run( exit_status )
  CALL do_stop( exit_status )
  !
END SUBROUTINE getenergy_run
  !
!----------------------------------------------------------------------------
