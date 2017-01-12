!
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!----------------------------------------------------------------------------
SUBROUTINE wl_qe_startup (lib_comm,nim,npl,nta,nbn,ndg,infile) BIND(C)
  !----------------------------------------------------------------------------
  !
  ! ... C wrapper for library interface to the Pwscf
  USE ISO_C_BINDING
  ! ... Set up the PWscf calculation for Quantum Espresso
  !
  USE environment,       ONLY : environment_start
  USE mp_global,         ONLY : mp_startup
  USE read_input,        ONLY : read_input_file
  USE command_line_options, ONLY : set_command_line
  !
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: lib_comm, nim, npl, nta, nbn, ndg
  CHARACTER (kind=C_CHAR), INTENT(IN) :: infile(*)
  CHARACTER(LEN=80)  :: input_file
  INTEGER :: i
  !USE command_line_options, ONLY : input_file_
  !INTEGER :: my_comm
  !
  input_file = ' '
  !
  ! ... Copying a string from C to Fortran is a bit ugly.
  DO i=1,80
      IF (infile(i) == C_NULL_CHAR) EXIT
      input_file = TRIM(input_file) // infile(i)
  END DO
  !
  CALL set_command_line (nimage=nim, npool=npl, ntg=nta, &
       nband=nbn, ndiag=ndg )
  CALL mp_startup (lib_comm, diag_in_band_group = .true.)
  CALL environment_start ('PWSCF')
  !
  CALL read_input_file ('PW', input_file )
  !CALL read_input_file ('PW', input_file_ )
  !
END SUBROUTINE wl_qe_startup
  !
!----------------------------------------------------------------------------
SUBROUTINE wl_qe_stop (exit_status) BIND(C)
  !----------------------------------------------------------------------------
  !
  USE ISO_C_BINDING
  !
  ! ... Finish the PWscf calculation for Quantum Espresso
  !
  IMPLICIT NONE
  INTEGER :: exit_status
  !
  CALL stop_run( exit_status )
  CALL do_stop( exit_status )
  !
END SUBROUTINE wl_qe_stop
  !
!----------------------------------------------------------------------------
SUBROUTINE wl_stop_run( exit_status ) BIND(C)
  !----------------------------------------------------------------------------
  !
  USE ISO_C_BINDING
  !
  ! ... Modified subroutine to restart the PWscf calculation after receiving
  ! ... the Wang-Landau inputs
  !
  USE io_global,          ONLY : ionode
  USE mp_global,          ONLY : mp_global_end
  USE environment,        ONLY : environment_end
  USE io_files,           ONLY : iuntmp, seqopn
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(IN) :: exit_status
  LOGICAL             :: exst, opnd, lflag
  !
  lflag = ( exit_status == 0 )
  IF ( lflag ) THEN
     !
     ! ... remove files needed only to restart
     !
     CALL seqopn( iuntmp, 'restart', 'UNFORMATTED', exst )
     CLOSE( UNIT = iuntmp, STATUS = 'DELETE' )
     !
     IF ( ionode ) THEN
        CALL seqopn( iuntmp, 'update', 'FORMATTED', exst )
        CLOSE( UNIT = iuntmp, STATUS = 'DELETE' )
        CALL seqopn( iuntmp, 'para', 'FORMATTED', exst )
        CLOSE( UNIT = iuntmp, STATUS = 'DELETE' )
     END IF
     !
  END IF
  !
  CALL close_files(lflag)
  !
  CALL print_clock_pw()
  !
  CALL clean_pw( .FALSE. )   ! ... Need to restart the subsequent PWscf runs
  !
END SUBROUTINE wl_stop_run
  !
!----------------------------------------------------------------------------
