MODULE Hpsi
!------------------------------------------------------------------------------
! STRUCTURE OF MODULE:
!   MODULE Hpsi
!     |_SUBROUTINE  HPsi
!
! DESCRIPTION:
!       This module calculates the result of the Madelung function applied by
!       the Hamiltonian .
!
! CONDITIONS AND ASSUMPTIONS:
!
!   03/01/2021 Yuanbo Li
!
!------------------------------------------------------------------------------

                         !>> GLOBAL <<!

  USE Constants, ONLY: DP                ! Double precision
  USE Constants, ONLY: pi
  USE Constants, ONLY: imag
  USE Fourier, ONLY: FFT
  USE PlaneWave, ONLY: qTable
  USE PlaneWave, ONLY: qVectors
  USE CD_potential, ONLY: CDPotential 
  USE KEDF_TF, ONLY: TFPotential  
  USE KEDF_VW, ONLY: VWPotential
  USE XC_LDA, ONLY: LDAPot
  USE XC_PBE, ONLY: PBEPot 
  USE IonElectron, ONLY: IonElectronPotentialRecip
  USE Hartree, ONLY: JPotentialPlus
  USE CellInfo, ONLY: cell
  USE TDDFT, ONLY: Hamiltonian_version
  USE TDDFT, ONLY: xc_funct
  USE TDDFT, ONLY: TDE0

  IMPLICIT NONE
 
CONTAINS

FUNCTION V_psi(psi,rhoR)
  ! calculate V(psi)
  IMPLICIT NONE

  COMPLEX(KIND=DP), DIMENSION(:,:,:), INTENT(IN) :: psi
  REAL(KIND=DP), DIMENSION(:,:,:), INTENT(IN) :: rhoR  
  !Electron density in real space, spin INDEPENDENT
  REAL(KIND=DP), DIMENSION(SIZE(psi, 1),SIZE(psi, 2), SIZE(psi, 3), 1):: rhoR_S
  !Electron density in real space, number of spin = 1
  REAL(KIND=DP), DIMENSION(SIZE(psi, 1),SIZE(psi, 2), SIZE(psi, 3), 1):: xcPotential
  REAL(KIND=DP), DIMENSION(SIZE(psi, 1),SIZE(psi, 2), SIZE(psi, 3)):: IonElectronPotential
  REAL(KIND=DP), DIMENSION(SIZE(psi, 1),SIZE(psi, 2), SIZE(psi, 3)):: HartreePotential
  COMPLEX(KIND=DP), DIMENSION(SIZE(psi, 1),SIZE(psi, 2), SIZE(psi, 3)):: V_psi

  ! ! calculate the energy for each term or not
  LOGICAL :: calcEnergy=.false.
  REAL(KIND=DP) :: energy

  CALL StartClock('V_psi')  

  !rhoR = ABS(psi)**2
  IonElectronPotential = FFT(IonElectronPotentialRecip(cell%cellReal,cell%ionTable, cell%elementTable))
  SELECT CASE (xc_funct)
    CASE(0)
      rhoR_S(:,:,:,1) = rhoR
      xcPotential = LDAPot(rhoR_S)
    CASE(1)
      CALL PBEPot(rhoR,xcPotential(:,:,:,1),calcEnergy,energy)
  END SELECT
  CALL JPotentialPlus(rhoR, HartreePotential, calcEnergy, energy)
  SELECT CASE (Hamiltonian_version)
    CASE(0)
      V_psi = TFPotential(rhoR)+ xcPotential(:,:,:,1) &
             + IonElectronPotential + HartreePotential
      !V_psi = IonElectronPotential
    CASE(1)
      V_psi = TFPotential(rhoR)- CDPotential(psi) + xcPotential(:,:,:,1) & 
             + IonElectronPotential + HartreePotential
      !write(*,*) "max of V_CD/V_TF",maxval((CDPotential(psi))/TFPotential(rhoR))
        !write(*,*) "max of V_CD",maxval((CDPotential(psi)))
    CASE(2)
      V_psi = TFPotential(rhoR)+VWPotential(rhoR)+ xcPotential(:,:,:,1) &
             + IonElectronPotential + HartreePotential
  END SELECT  

  CALL StopClock('V_psi')
  RETURN

END FUNCTION V_psi

FUNCTION H_psi(psi,rhoR)
  ! calculate H(psi)*psi
  IMPLICIT NONE

  COMPLEX(KIND=DP), DIMENSION(:,:,:), INTENT(IN) :: psi
  REAL(KIND=DP), DIMENSION(:,:,:), INTENT(IN) :: rhoR
  COMPLEX(KIND=DP), DIMENSION(SIZE(psi, 1),SIZE(psi, 2), SIZE(psi, 3)):: KineticPotential
  COMPLEX(KIND=DP), DIMENSION(SIZE(psi, 1),SIZE(psi, 2), SIZE(psi, 3)):: H_psi
  REAL(KIND=DP) :: ii, jj, kk

  CALL StartClock('H_psi')

  IF (Hamiltonian_version==2) THEN
    H_psi = V_psi(psi,rhoR) * psi
  ELSE
    KineticPotential = 0.5_DP * ( FFT(FFT(REAL(psi))*qTable**2) &
      + imag* FFT(FFT(AIMAG(psi))*qTable**2) )
    H_psi = KineticPotential + V_psi(psi,rhoR) * psi 
    !H_psi = KineticPotential
    !H_psi = V_psi(psi) * psi
  END IF

  TDE0=0.0_DP
  DO ii=0,size(psi,1)-1
    DO jj=0,size(psi,2)-1
      DO kk=0,size(psi,3)-1
        TDE0 = TDE0 + REAL(CONJG(psi(ii,jj,kk))*H_psi(ii,jj,kk))
        !WRITE(*,*) "i= ",ii," j=",jj," k=",kk," TDEE = ", TDEE
      ENDDO 
    ENDDO
  ENDDO

  CALL StopClock('H_psi')
  RETURN
END FUNCTION H_psi

END MODULE Hpsi
