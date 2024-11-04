MODULE Propagation_RK
!------------------------------------------------------------------------------
! STRUCTURE OF MODULE:
!   MODULE Propagation_RK
!     |_SUBROUTINE 
!     |_SUBROUTINE 
!
! DESCRIPTION:
!       This module calculates the propagation of the Madelung wave function by
!       Runge-Kutta method.
!
! REFERENCE : 
!    1. Rehn, D. A.; Shen, Y.; Buchholz, M. E.; Dubey, M.; Namburu, R.; Reed, E.
!    J. ODE integration schemes for plane-wave real-time time-dependent density
!    functional theory. J. Chem. Phys. 2019, 150, No. 014101.
! CONDITIONS AND ASSUMPTIONS:
!
!   03/03/2021 Yuanbo Li
!
!------------------------------------------------------------------------------

                         !>> GLOBAL <<!

  USE Constants, ONLY: DP                ! Double precision
  USE Constants, ONLY: imag 
  USE Hpsi, ONLY: H_psi
  USE TDDFT, ONLY: t_dt  
  USE Sys, ONLY: psi
  USE Sys, ONLY: psi_1
  USE Sys, ONLY: psi_2

  IMPLICIT NONE
  
CONTAINS

SUBROUTINE PropagationRK()
  IMPLICIT NONE
  COMPLEX(KIND=DP), DIMENSION(SIZE(psi, 1),SIZE(psi, 2), SIZE(psi, 3)):: psi1  
  COMPLEX(KIND=DP), DIMENSION(SIZE(psi, 1),SIZE(psi, 2), SIZE(psi, 3)):: psi2
  COMPLEX(KIND=DP), DIMENSION(SIZE(psi, 1),SIZE(psi, 2), SIZE(psi, 3)):: psi3
  REAL(KIND=DP), DIMENSION(SIZE(psi, 1),SIZE(psi, 2), SIZE(psi, 3)):: rhoR
  COMPLEX(KIND=DP), DIMENSION(SIZE(psi, 1),SIZE(psi, 2), SIZE(psi, 3)):: K1
  COMPLEX(KIND=DP), DIMENSION(SIZE(psi, 1),SIZE(psi, 2), SIZE(psi, 3)):: K2
  COMPLEX(KIND=DP), DIMENSION(SIZE(psi, 1),SIZE(psi, 2), SIZE(psi, 3)):: K3
  COMPLEX(KIND=DP), DIMENSION(SIZE(psi, 1),SIZE(psi, 2), SIZE(psi, 3)):: K4

  CALL StartClock('Propagation_RK')

  !WRITE(*,*) "Run Propagation_RK()"  
  ! RK2
  !psi1 = psi - 0.5_DP * imag * t_dt * H_psi(psi)
  !psi = psi - imag * t_dt * H_psi(psi1)
  
  ! RK4    
  rhoR=ABS(psi_1)**2+ABS(psi_2)**2
  K1 = - imag*t_dt*H_psi(psi_1,rhoR)
  !WRITE(*,*) "K1=",K1(0,0,0)
  psi1 = psi_1 + 0.5_DP * K1
  rhoR=ABS(psi1)**2+ABS(psi_2)**2
  K2 = - imag*t_dt*H_psi(psi1,rhoR)
  psi2 = psi_1 + 0.5_DP * K2
  rhoR=ABS(psi2)**2+ABS(psi_2)**2
  K3 = - imag*t_dt*H_psi(psi2,rhoR)
  psi3 = psi_1 + K3
  rhoR=ABS(psi3)**2+ABS(psi_2)**2
  K4 = - imag*t_dt*H_psi(psi3,rhoR)
  psi_1 = psi_1 + 1.0_DP / 6.0_DP *(K1+2*K2+2*K3+K4)
  !psi_1 = psi_1 + K1

  !rhoR=ABS(psi_1+psi_2)**2
  !K1 = - imag*t_dt*H_psi(psi_2,rhoR)
  !psi1 = psi_2 + 0.5_DP * K1
  !rhoR=ABS(psi1)**2+ABS(psi_1)**2
  !K2 = - imag*t_dt*H_psi(psi1,rhoR)
  !psi2 = psi_2 + 0.5_DP * K2
  !rhoR=ABS(psi2)**2+ABS(psi_1)**2
  !K3 = - imag*t_dt*H_psi(psi2,rhoR)
  !psi3 = psi_2 + K3
  !rhoR=ABS(psi3)**2+ABS(psi_1)**2
  !K4 = - imag*t_dt*H_psi(psi3,rhoR)
  !psi_2 = psi_2 + 1.0_DP / 6.0_DP *(K1+2*K2+2*K3+K4)
  !psi_2 = psi_2 +K1

  ! test taylor 1st
  !psi=psi1

  CALL StopClock('Propagation_RK')
END SUBROUTINE PropagationRK

END MODULE Propagation_RK
