!
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!----------------------------------------------------------------------------
SUBROUTINE owl_do_pwscf ( exit_status ) 
  !----------------------------------------------------------------------------
  !
  ! ... The modified subroutine to run multiple PWscf calculations, providing 
  ! ... inputs for Wang-Landau algorithm
  !
  USE io_global,        ONLY : stdout, ionode, ionode_id
  USE parameters,       ONLY : ntypx, npk, lmaxx
  USE cell_base,        ONLY : fix_volume, fix_area
  USE control_flags,    ONLY : conv_elec, gamma_only, ethr, lscf, twfcollect
  USE control_flags,    ONLY : conv_ions, istep, nstep, restart, lmd, lbfgs
  USE force_mod,        ONLY : lforce, lstres, sigma, force
  USE check_stop,       ONLY : check_stop_init, check_stop_now
  USE mp_images,        ONLY : intra_image_comm
  USE extrapolation,    ONLY : update_file, update_pot
  USE scf,              ONLY : rho
  USE lsda_mod,         ONLY : nspin
  USE fft_base,         ONLY : dfftp
  USE qmmm,             ONLY : qmmm_initialization, qmmm_shutdown, &
                               qmmm_update_positions, qmmm_update_forces
  !
  IMPLICIT NONE
  INTEGER, INTENT(OUT) :: exit_status
  INTEGER :: idone 
  !
  ! counter of electronic + ionic steps done in this run
  !
  exit_status = 0
  IF ( ionode ) WRITE( unit = stdout, FMT = 9010 ) ntypx, npk, lmaxx
  !
  IF (ionode) CALL plugin_arguments()
  CALL plugin_arguments_bcast( ionode_id, intra_image_comm )
  ! 
  CALL qmmm_initialization()
  !
  ! ... convert to internal variables
  !
  IF ( gamma_only ) WRITE( UNIT = stdout, &
     & FMT = '(/,5X,"gamma-point specific algorithms are used")' )
  !
  ! call to void routine for user defined / plugin patches initializations
  !
  CALL plugin_initialization()
  !
  CALL check_stop_init()
  !
  CALL setup ()
  !
  CALL qmmm_update_positions()
  !
  CALL init_run()
  !
  ! ... dry run: code will stop here if called with exit file present
  ! ... useful for a quick and automated way to check input data
  !
  IF ( check_stop_now() ) THEN
     CALL punch( 'config' )
     exit_status = 255
     RETURN
  ENDIF
  !
  main_loop: DO idone = 1, nstep
     !
     ! ... electronic self-consistency or band structure calculation
     !
     IF ( .NOT. lscf) THEN
        CALL non_scf ()
     ELSE
        CALL electrons()
     END IF
     !
     ! ... code stopped by user or not converged
     !
     IF ( check_stop_now() .OR. .NOT. conv_elec ) THEN
        IF ( check_stop_now() ) exit_status = 255
        IF ( .NOT. conv_elec )  exit_status =  2
        ! workaround for the case of a single k-point
        twfcollect = .FALSE.
        CALL punch( 'config' )
        RETURN
     ENDIF
     !
     ! ... ionic section starts here
     !
     CALL start_clock( 'ions' )
     conv_ions = .TRUE.
     !
     ! ... file in CASINO format written here if required
     !
     IF ( lmd ) THEN
        CALL pw2casino( istep )
     ELSE
        CALL pw2casino( 0 )
     END IF
     !
     ! ... force calculation
     !
     IF ( lforce ) CALL forces()
     !
     ! ... stress calculation
     !
     IF ( lstres ) CALL stress ( sigma )
     !
     ! ... send out forces to MM code in QM/MM run
     !
     CALL qmmm_update_forces( force, rho%of_r, nspin, dfftp)
     !
     IF ( lmd .OR. lbfgs ) THEN
        !
        if (fix_volume) CALL impose_deviatoric_stress(sigma)
        if (fix_area)  CALL  impose_deviatoric_stress_2d(sigma)
        !
        ! ... save data needed for potential and wavefunction extrapolation
        !
        CALL update_file ( )
        !
        ! ... ionic step (for molecular dynamics or optimization)
        !
        CALL move_ions ( idone )
        !
        ! ... then we save restart information for the new configuration
        !
        IF ( idone <= nstep .AND. .NOT. conv_ions ) CALL punch( 'config' )
        !
     END IF
     !
     CALL stop_clock( 'ions' )
     !
     ! ... exit condition (ionic convergence) is checked here
     !
     IF ( conv_ions ) EXIT main_loop
     !
     ! ... receive new positions from MM code in QM/MM run
     !
     CALL qmmm_update_positions()
     !
     ! ... terms of the hamiltonian depending upon nuclear positions
     ! ... are reinitialized here
     !
     IF ( lmd .OR. lbfgs ) THEN
        !
        ! ... update the wavefunctions, charge density, potential
        ! ... update_pot initializes structure factor array as well
        !
        CALL update_pot()
        !
        ! ... re-initialize atomic position-dependent quantities
        !
        CALL hinit1()
        !
     END IF
     ! ... Reset convergence threshold of iterative diagonalization for
     ! ... the first scf iteration of each ionic step (after the first)
     !
     ethr = 1.0D-6
     !
  END DO main_loop
  !
  ! ... save final data file
  !
  CALL punch('all')
  !
  CALL qmmm_shutdown()
  !
  IF ( .NOT. conv_ions )  exit_status =  3
  RETURN
  !
9010 FORMAT( /,5X,'Current dimensions of program PWSCF are:', &
           & /,5X,'Max number of different atomic species (ntypx) = ',I2,&
           & /,5X,'Max number of k-points (npk) = ',I6,&
           & /,5X,'Max angular momentum in pseudopotentials (lmaxx) = ',i2)
  !
END SUBROUTINE owl_do_pwscf
