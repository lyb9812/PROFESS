MODULE Propagation_CN
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
! CONDITIONS AND ASSUMPTIONS:
!
!   03/03/2021 Yuanbo Li
!   not finished 
!------------------------------------------------------------------------------

                         !>> GLOBAL <<!

  USE Constants, ONLY: DP                ! Double precision
  USE Constants, ONLY: imag 
  USE Hpsi, ONLY: H_psi
  USE Hpsi, ONLY: V_psi
  USE TDDFT, ONLY: t_dt  
  USE Sys, ONLY: psi

  IMPLICIT NONE

CONTAINS

SUBROUTINE PropagationCN()
  IMPLICIT NONE
  COMPLEX(KIND=DP), DIMENSION(SIZE(psi, 1),SIZE(psi, 2), SIZE(psi,3))::psi1
  COMPLEX(KIND=DP), DIMENSION(SIZE(psi, 1),SIZE(psi, 2), SIZE(psi,3))::psi2
  INTEGER :: ii,jj
  REAL(KIND=DP) :: drho
  !psi1=psi
  !DO ii=1,2
  !  psi2 = (1-imag*H_psi(psi1)*t_dt*0.5)/(1+imag*H_psi(psi1)*t_dt*0.5)*psi1
  !  IF (ABS(ABS(psi2)**2-ABS(psi1)**2)<1.0E-6) THEN
  !    psi1=psi2
  !    exit
  !  ELSE
  !    psi1=0.5*psi1+0.5*psi2
  !  END IF
  !ENDDO

  !psi=psi1

  ! test 2024-6-28
  !psi1=(1-imag*V_psi(psi)*t_dt*0.5)/(1+imag*V_psi(psi)*t_dt*0.5)*psi
  psi=psi

END SUBROUTINE PropagationCN

END MODULE Propagation_CN
