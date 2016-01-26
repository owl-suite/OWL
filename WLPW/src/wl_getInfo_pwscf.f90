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
  f_etot = etot   ! ... total energy of system (in eV/unit cell)
  !
END SUBROUTINE get_natom_ener
  !
!----------------------------------------------------------------------------
SUBROUTINE get_array (pos_array, cell_array)
  !----------------------------------------------------------------------------
  !
  ! ... Extract the position array from the Quantum Espresso modules
  !
  USE constants,     ONLY : bohr_radius_angs
  USE cell_base,     ONLY : alat, at
  USE ions_base,     ONLY : tau, nat
  !
  IMPLICIT NONE
  REAL(KIND=8) :: pos_array(3,nat)
  REAL(KIND=8) :: cell_array(3,3)
  !
  pos_array(:,:) = tau(:,:) * alat * bohr_radius_angs   ! ... (in angstrom)
  cell_array(:,:) = at(:,:) * alat * bohr_radius_angs   ! ... (in angstrom)
  !
END SUBROUTINE get_array
  !
!----------------------------------------------------------------------------
SUBROUTINE wl_pass_array (pos_array)
  !----------------------------------------------------------------------------
  !
  ! ... Extract the position array from the Quantum Espresso modules
  !
  USE constants,     ONLY : bohr_radius_angs
  USE cell_base,     ONLY : alat
  USE ions_base,     ONLY : tau, nat, ityp, if_pos, extfor
  !
  IMPLICIT NONE
  REAL(KIND=8) :: pos_array(3,nat)
  !
  ! ... Feed the position array from Wang-Landau algorithm
  !
  tau(:,:) = (pos_array(:,:) / alat) / bohr_radius_angs   ! ... (in Cartesian)
  tau(3,nat) = 0.55
  !
END SUBROUTINE wl_pass_array
  !
!----------------------------------------------------------------------------
