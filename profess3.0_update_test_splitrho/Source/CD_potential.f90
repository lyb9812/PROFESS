MODULE CD_potential
!------------------------------------------------------------------------------
! STRUCTURE OF MODULE:
!   MODULE CD_potential
!     |_SUBROUTINE  CDPotential
!
! DESCRIPTION:
!       This module calculates the current dependent potential .
!   REFERENCE:
!   1. A. J. White, O. Certik, Y. H. Ding, S. X. Hu, and L. A. Collins,
!   Physical Review B 98, 144302 (2018) 
! CONDITIONS AND ASSUMPTIONS:
!
!   02/28/2021 Yuanbo Li
!   there is somenthing wrong when doing FFT 
!
!------------------------------------------------------------------------------

                         !>> GLOBAL <<!

  USE Constants, ONLY: DP                ! Double precision
  USE Constants, ONLY: pi
  USE Constants, ONLY: imag
  USE Constants, ONLY: boltzmann
  USE Fourier, ONLY: FFT
  USE PlaneWave, ONLY: qTable
  USE PlaneWave, ONLY: qVectors
  USE TDDFT, ONLY: t_temperature

  IMPLICIT NONE
 
CONTAINS

FUNCTION CDPotential(psi)
  IMPLICIT NONE

  COMPLEX(KIND=DP), DIMENSION(:,:,:), INTENT(IN) :: psi
  REAL(KIND=DP), DIMENSION(SIZE(psi, 1),SIZE(psi, 2), SIZE(psi, 3)):: rhoR
  REAL(KIND=DP), DIMENSION(SIZE(psi, 1),SIZE(psi, 2), SIZE(psi, 3)):: kF
  REAL(KIND=DP), DIMENSION(SIZE(psi, 1),SIZE(psi, 2), SIZE(psi, 3)):: CT
  REAL(KIND=DP), DIMENSION(SIZE(psi, 1),SIZE(psi, 2), SIZE(psi, 3),3):: current
  COMPLEX(KIND=DP), DIMENSION(SIZE(psi, 1),SIZE(psi, 2), SIZE(psi, 3),3):: current1
  COMPLEX(KIND=DP), DIMENSION(SIZE(qVectors, 1),SIZE(qVectors, 2), SIZE(qVectors, 3),3):: currentQ
  REAL(KIND=DP), DIMENSION(SIZE(psi, 1),SIZE(psi, 2), SIZE(psi, 3)):: CDPotential
  COMPLEX(KIND=DP), DIMENSION(SIZE(qVectors, 1),SIZE(qVectors, 2), SIZE(qVectors, 3)):: CDPotential1
  COMPLEX(KIND=DP), DIMENSION(SIZE(qVectors, 1),SIZE(qVectors, 2), SIZE(qVectors, 3)):: FFTreal
  COMPLEX(KIND=DP), DIMENSION(SIZE(qVectors, 1),SIZE(qVectors, 2), SIZE(qVectors, 3)):: FFTimag
  REAL(KIND=DP)::kF0,CT0  ! average kF  
  INTEGER::ii,jj,kk

  CALL StartClock('CD_Potential')
  
  rhoR = ABS(psi)**2
  FFTreal = FFT(REAL(psi))
  FFTimag = FFT(AIMAG(psi))

  !write(*,*) "size of psi:",SIZE(psi,1),SIZE(psi,2),SIZE(psi,3)
  !write(*,*) "size of kvector:",SIZE(qVectors,1),SIZE(qVectors,2),SIZE(qVectors,3)
  !write(*,*) "max of FFTreal REAL",maxval(REAL(FFTreal))
  !write(*,*) "max of FFTreal IMAG",maxval(AIMAG(FFTreal))
  !write(*,*) "max of FFTimag REAL",maxval(REAL(FFTimag))
  !write(*,*) "max of FFTimag IMAG",maxval(AIMAG(FFTimag))
  !write(*,*) "write FFTImag:  start: ",FFTimag," end"

  !current1(:,:,:,1) = FFT(qVectors(:,:,:,1)*FFTimag) 
  current1(:,:,:,1) = FFT(imag*qVectors(:,:,:,1)*FFTreal) &
+imag*FFT(imag*qVectors(:,:,:,1)*FFTimag) 
  current1(:,:,:,2) = FFT(imag*qVectors(:,:,:,2)*FFTreal) &
+imag*FFT(imag*qVectors(:,:,:,1)*FFTimag)
  current1(:,:,:,3) = FFT(imag*qVectors(:,:,:,3)*FFTreal) &
+imag*FFT(imag*qVectors(:,:,:,3)*FFTimag)

  !write(*,*) "max of current1 x",maxval(ABS(current1(:,:,:,1)))
  !write(*,*) "max of current1 y",maxval(ABS(current1(:,:,:,2)))
  !write(*,*) "max of current1 z",maxval(ABS(current1(:,:,:,3)))

  current(:,:,:,1) = AIMAG(CONJG(psi)*current1(:,:,:,1))  
  current(:,:,:,2) = AIMAG(CONJG(psi)*current1(:,:,:,2))
  current(:,:,:,3) = AIMAG(CONJG(psi)*current1(:,:,:,3))   

  !write(*,*) "max of psi",maxval(abs(psi))
  !write(*,*) "elelment of current",psi(2,2,2),current1(2,2,2,1),current1(2,2,2,2),current1(2,2,2,3)
  !write(*,*) "max of current x",maxval(current(:,:,:,1))
  !write(*,*) "max of current y",maxval(current(:,:,:,2))
  !write(*,*) "max of current z",maxval(current(:,:,:,3))
  !WRITE(*,*) "current"
  !WRITE(*,*) current

  currentQ(:,:,:,1) = FFT(current(:,:,:,1))
  currentQ(:,:,:,2) = FFT(current(:,:,:,2))
  currentQ(:,:,:,3) = FFT(current(:,:,:,3))

  !WRITE(*,*) "currentQ"
  !WRITE(*,*) currentQ
  !write(*,*) "max of currentQ1 REAL",maxval(REAL(currentQ(:,:,:,1)))
  !write(*,*) "max of currentQ1 IMAG",maxval(AIMAG(currentQ(:,:,:,1)))

  kF0=(3.0_DP*pi**2*(SUM(rhoR)/(SIZE(rhoR,1)*SIZE(rhoR,2)*SIZE(rhoR,3))))**(1.0_DP / 3.0_DP)
  kF = (3.0_DP*pi**2*rhoR)**(1.0_DP / 3.0_DP)
  CT = ((1.69271_DP* (2.0_DP*boltzmann*t_temperature/(kF**2))**(0.5_DP))**(3.6_DP) +1)**(1.0_DP / 3.6_DP )
  !CT0 = ((1.69271_DP*(2.0_DP*boltzmann*t_temperature/(kF0**2))**(0.5_DP))**(3.6_DP) +1)**(1.0_DP / 3.6_DP )
  CDPotential1 = currentQ(:,:,:,1)*qVectors(:,:,:,1) & 
                +currentQ(:,:,:,2)*qVectors(:,:,:,2) & 
               +currentQ(:,:,:,3)*qVectors(:,:,:,3)

  !write(*,*) "qvectors: ",qVectors(:,:,:,1)
  !write(*,*) "write each element 2,2,2: ",currentQ(2,2,2,1),qVectors(2,2,2,1),currentQ(2,2,2,2),qVectors(2,2,2,2),currentQ(2,2,2,3),qVectors(2,2,2,3)
  !write(*,*) "max of CDPotential1 REAL",maxval(ABS(REAL(CDPotential1)))
  !write(*,*) "max of CDPotential1 IMAG",maxval(ABS(AIMAG(CDPotential1)))
  !write(*,*) "max of CDPotential1/q REAL",maxval(ABS(REAL(imag*CDPotential1/qTable)))
  !write(*,*) "max of CDPotential1/q IMAG",maxval(ABS(AIMAG(imag*CDPotential1/qTable)))
   !WRITE(*,*) "CT=",CT
  !WRITE(*,*) "kF0=",kF0
  
  !WRITE(*,*) "qTable: ",qTable(1,1,1),qTable
  CDPotential1=imag*CDPotential1/qTable
  !WRITE(*,*) "CDPotential before FFT",CDPotential1
  CDPotential1(1,1,1)=0.0_DP

  CDPotential=FFT(CDPotential1)
  !CDPotential = FFT(imag*CDPotential1/qTable)
  !CDPotential(1,1,1) = 0.0_DP

  !write(*,*) "max of initial CDPotential",maxval(CDPotential)
  DO ii=1,SIZE(CDPotential,1)
    DO jj=1,SIZE(CDPotential,2)
      DO kk=1,SIZE(CDPotential,3)
        IF (ISNAN(CDPotential(ii,jj,kk))) THEN
          CDPotential(ii,jj,kk)=0.00_DP
        END IF
      ENDDO
    ENDDO
  ENDDO 

  !write(*,*) "max of kF",maxval(kF)
  !write(*,*) "max of CT",maxval(CT)
  !write(*,*) "max of CDPotential fixed",maxval(CDPotential)

  !WRITE(*,*) "CDPotential before FFT"
  !WRITE(*,*) CDPotential
  !CDPotential = CT* FFT(CDPotential)*(pi**3)/(2.0_DP * kF0**2 )
  !CDPotential = CT* CDPotential*(pi**3)/(2.0_DP * kF0**2 )
  CDPotential = CT* CDPotential*(pi**3)/(2.0_DP * kF**2 )
  !CDPotential = CT0* CDPotential*(pi**3)/(2.0_DP * kF**2 )
  !CDPotential = CT0* CDPotential*(pi**3)/(2.0_DP * kF0**2 )
  !CDPotential = CT* CDPotential*(pi**3)

  !WRITE(*,*) "kF0",kF0 
  !WRITE(*,*) "qTable"
  !WRITE(*,*) qTable(1,1,1)  

  ! test
  !WRITE(*,*) "CDPotential print start "
  !DO ii=1,SIZE(CDPotential,1)
  !  DO jj=1,SIZE(CDPotential,2)
  !    DO kk=1,SIZE(CDPotential,3)
  !      WRITE(*,*) ii,jj,kk,CDPotential(ii,jj,kk)
  !    ENDDO
  !  ENDDO
  !ENDDO
  !WRITE(*,*) "CDPotential print end "

  !OPEN(unit=401,FILE="current.txt",ACTION='READWRITE')
  !WRITE(401,*) "current print start "
  !DO ii=1,SIZE(CDPotential,1)
  !  DO jj=1,SIZE(CDPotential,2)
  !    DO kk=1,SIZE(CDPotential,3)
  !      WRITE(401,*) ii,jj,kk,current(ii,jj,kk,1),current(ii,jj,kk,2) &
  !         ,current(ii,jj,kk,3)
  !    ENDDO
  !  ENDDO
  !ENDDO
  !WRITE(401,*) "current print end "
  !WRITE(401,*)

  !OPEN(unit=400,FILE="currentQ.txt",ACTION='READWRITE')
  !WRITE(400,*) "currentQ print start "
  !DO ii=1,SIZE(CDPotential,1)
  !  DO jj=1,SIZE(CDPotential,2)
  !    DO kk=1,SIZE(CDPotential,3)
  !      WRITE(400,*) ii,jj,kk,currentQ(ii,jj,kk,1),currentQ(ii,jj,kk,2) &
  !         ,currentQ(ii,jj,kk,3)
  !    ENDDO
  !  ENDDO
  !ENDDO
  !WRITE(400,*) "currentQ print end "
  !WRITE(400,*) 

  CALL StopClock('CD_Potential')
 
  RETURN

END FUNCTION CDPotential

END MODULE CD_potential
