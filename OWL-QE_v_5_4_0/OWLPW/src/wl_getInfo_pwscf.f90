!
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt.
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
  f_etot = etot   ! ... total energy of system (in Ry)
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
  ! ... Extract the position array from the Quantum Espresso modules
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
SUBROUTINE pass_pos_array (pos_array)
  !----------------------------------------------------------------------------
  !
  ! ... Extract the position array from the Quantum Espresso modules
  ! ... after updating the atomic positions from Wang-Landau codes
  !
  USE constants,     ONLY : bohr_radius_angs
  USE cell_base,     ONLY : alat
  USE ions_base,     ONLY : tau, nat
  !
  IMPLICIT NONE
  REAL(KIND=8) :: pos_array(3,nat)
  !
  ! ... Feed the updated position array from Wang-Landau codes
  !
  tau(:,:) = (pos_array(:,:) / alat) / bohr_radius_angs   ! ... (in Cartesian)
  !
!  tau(1,nat) = 0.7744183 ! ... test; assuming this change is made from WL codes
!  tau(2,nat) = 0.5192893 ! ... test; assuming this change is made from WL codes
!  tau(3,nat) = 0.8074200 ! ... test; assuming this change is made from WL codes
  !
END SUBROUTINE pass_pos_array
  !
!----------------------------------------------------------------------------
SUBROUTINE pass_cell_array (cell_array)
  !----------------------------------------------------------------------------
  !
  ! ... Extract the position array from the Quantum Espresso modules
  ! ... after updating the cell lattice vector from Wang-Landau codes
  !
  USE constants,     ONLY : bohr_radius_angs
  USE cell_base,     ONLY : alat, at, omega, bg
  USE ions_base,     ONLY : tau, nat
  !
  IMPLICIT NONE
  REAL(KIND=8) :: cell_array(3,3)
  !
  ! ... Convert the Cartesian to fractional coordinates for the position array
  ! ... using the old lattice vector
  !
  CALL cryst_to_cart( nat, tau, bg, -1 )
  !
  ! ... Feed the updated cell vector array from Wang-Landau codes
  !
  at(:,:) = (cell_array(:,:) / alat) / bohr_radius_angs   ! ... (in Cartesian)
!  at(3,3) = 1.015   ! ... test; assuming this change is made from WL codes
  !
  ! ... Update the cell volume, position array, and reciprocal lattice vector
  ! ... using the updated lattice vector
  !
  CALL volume( alat, at(1,1), at(1,2), at(1,3), omega )
  !
  CALL cryst_to_cart( nat, tau, at, 1 )
  !
  CALL recips( at(1,1), at(1,2), at(1,3), bg(1,1), bg(1,2), bg(1,3) )
  !
END SUBROUTINE pass_cell_array
  !
!----------------------------------------------------------------------------
