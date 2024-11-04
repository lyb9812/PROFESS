MODULE Propagation_Split
!------------------------------------------------------------------------------
! STRUCTURE OF MODULE:
!   MODULE Propagation_Split
!     |_SUBROUTINE 
!     |_SUBROUTINE 
!
! DESCRIPTION:
!       This module calculates the propagation of the Madelung wave function by
!       the split operator method.
! REFERENCE:
!       Bandrauk AD, Shen H. 1991. Chern. Phys. Lett. 176:428-31
!
! CONDITIONS AND ASSUMPTIONS:
!
!   04/08/2021 Yuanbo Li
!   not finished yet
!------------------------------------------------------------------------------

                         !>> GLOBAL <<!

  USE Constants, ONLY: DP                ! Double precision
  USE Constants, ONLY: imag 
  USE Hpsi, ONLY: V_psi
  USE TDDFT, ONLY: t_dt  
  USE Sys, ONLY: psi
  USE PlaneWave, ONLY: qTable  
  USE Fourier, ONLY: FFT

  IMPLICIT NONE
  REAL(KIND=DP):: gamma0=1.3512072_DP  ! parameter in the third order split operator method

CONTAINS

SUBROUTINE PropagationSplit()
  IMPLICIT NONE
  
  ! intermediate variables
  COMPLEX(KIND=DP), DIMENSION(SIZE(psi, 1),SIZE(psi, 2), SIZE(psi, 3)):: psi1
  COMPLEX(KIND=DP), DIMENSION(SIZE(psi, 1),SIZE(psi, 2), SIZE(psi, 3)):: psi2
  COMPLEX(KIND=DP), DIMENSION(SIZE(psi, 1),SIZE(psi, 2), SIZE(psi, 3)):: psi3
  COMPLEX(KIND=DP), DIMENSION(SIZE(psi, 1),SIZE(psi, 2), SIZE(psi, 3)):: psi4
  COMPLEX(KIND=DP), DIMENSION(SIZE(psi, 1),SIZE(psi, 2), SIZE(psi, 3)):: psi5
  COMPLEX(KIND=DP), DIMENSION(SIZE(psi, 1),SIZE(psi, 2), SIZE(psi, 3)):: psi6

  CALL StartClock('Propagation_Split')

  ! second order
!  psi1 = FFT(exp(- 0.25_DP * imag*qTable**2 *t_dt)* FFT(REAL(psi)))+ imag* FFT(exp(- 0.25_DP * imag*qTable**2 *t_dt)* FFT(AIMAG(psi)))
!  psi2 = exp(-imag*t_dt* V_psi(psi1))*psi1
!  psi = FFT(exp(- 0.25_DP * imag*qTable**2 *t_dt)* FFT(REAL(psi2)))+ imag* FFT(exp(-0.25_DP * imag*qTable**2 *t_dt)* FFT(AIMAG(psi2)))

! --------------------------------------------------coment for test 2024-06-29
  !  third order
!   psi1 = FFT(exp(-0.25_DP*gamma0*imag*qTable**2*t_dt)*FFT(REAL(psi)))+imag* &
!          FFT(exp(-0.25_DP*gamma0*imag*qTable**2*t_dt)*FFT(AIMAG(psi)))  
!   psi2 = exp(-imag*t_dt*gamma0*V_psi(psi1))*psi1
!   psi3 = FFT(exp(-0.25_DP*(1.0_DP-gamma0)*imag*qTable**2*t_dt)* &
!          FFT(REAL(psi2)))+imag*FFT(exp(-0.25_DP*(1.0_DP-gamma0)* &
!         imag*qTable**2*t_dt)*FFT(AIMAG(psi2)))
!   psi4 = exp(-imag*t_dt*(1.0_DP-2.0_DP*gamma0)*V_psi(psi3))*psi3
!   psi5 = FFT(exp(-0.25_DP*(1.0_DP-gamma0)*imag*qTable**2*t_dt)* & 
!         FFT(REAL(psi4)))+imag*FFT(exp(-0.25_DP*(1.0_DP-gamma0)* &
!          imag*qTable**2*t_dt)*FFT(AIMAG(psi4)))
!   psi6 = exp(-imag*t_dt*gamma0*V_psi(psi5))*psi5
!   psi = FFT(exp(-0.25_DP*gamma0*imag*qTable**2*t_dt)*FFT(REAL(psi6)))+imag* &
!         FFT(exp(-0.25_DP*gamma0*imag*qTable**2*t_dt)*FFT(AIMAG(psi6))) 
! --------------------------------------------------coment for test 2024-06-29

  ! psi1 = FFT(exp(-0.25_DP*gamma0*imag*qTable**2*t_dt)*FFT(psi))
  ! psi2 = exp(-imag*t_dt*gamma0*V_psi(psi1))*psi1
  ! psi3 = FFT(exp(-0.25_DP*(1.0_DP-gamma0)*imag*qTable**2*t_dt)* &
  !        FFT(psi2))
  ! psi4 = exp(-imag*t_dt*(1.0_DP-2.0_DP*gamma0)*V_psi(psi3))*psi3
  ! psi5 = FFT(exp(-0.25_DP*(1.0_DP-gamma0)*imag*qTable**2*t_dt)* &
  !        FFT(psi4))
  ! psi6 = exp(-imag*t_dt*gamma0*V_psi(psi5))*psi5
  ! psi = FFT(exp(-0.25_DP*gamma0*imag*qTable**2*t_dt)*FFT(psi6))

  CALL StopClock('Propagation_Split')
END SUBROUTINE PropagationSplit

END MODULE Propagation_Split
