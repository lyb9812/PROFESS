MODULE TDDFT
!------------------------------------------------------------------------------
! STRUCTURE OF MODULE:
!   MODULE TDDFT
!     |_SUBROUTINE 
!     |_SUBROUTINE 
!
! DESCRIPTION:
!       This module wraps the time-dependent dft calculations.
!
! CONDITIONS AND ASSUMPTIONS:
!
!   02/26/2021 Yuanbo Li
!
!------------------------------------------------------------------------------

                         !>> GLOBAL <<!

  USE Constants, ONLY: DP                ! Double precision


  IMPLICIT NONE
  
  REAL(KIND=DP) :: t_dt ! time interval
  REAL(KIND=DP) :: t_temperature  ! temperature
  INTEGER :: t_Nstep ! number of total steps
  !CHARACTER(LEN=40) :: rho_in ! the name of the input file of the electronic density
  !CHARACTER(LEN=40) :: current_in ! the name of the input file of the current density
  INTEGER :: Hamiltonian_version ! the version of the form of Hamiltonian
  !    0 : H=-\nabla^2/2+v_TF+v_Hartree+v_ext+v_xc
  !    1 : H=-\nabla^2/2+v_TF+v_Hartree+v_ext+v_xc+v_CD
  !    2 : H=v_vW+v_TF+v_Hartree+v_ext+v_xc
  INTEGER :: propagation_method  ! the method to propagate the Madelung wave function
  !    0 : Runge-Kutta method 
  !    1 : Crank-Nicolson method 
  !    2 : Split operator method 
  INTEGER :: ion_method ! the method to move the ions
  !    0 : Erenfest dynamics
  INTEGER :: xc_funct  ! type of exchange-correlation functional
  !    0 : LDA
  !    1 : PBE
  REAL(KIND=DP) :: TDE ! time dependent total energy 
  REAL(KIND=DP) :: TDE0 ! <\psi|H|\psi>+ion
  REAL(KIND=DP) :: IonKinetic ! the total kinetic energy of the ions
  REAL(KIND=DP) :: t_time ! the current time
  REAL(KIND=DP) :: totalcharge ! for testing 
CONTAINS


END MODULE TDDFT
