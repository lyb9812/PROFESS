MODULE Propagation
!------------------------------------------------------------------------------
! STRUCTURE OF MODULE:
!   MODULE Propagation
!     |_SUBROUTINE 
!     |_SUBROUTINE 
!
! DESCRIPTION:
!       This module warps the calculation of  the movement of the ions and the
!       propagation of the Madelung wavefunction.
!
! CONDITIONS AND ASSUMPTIONS:
!
!   03/05/2021 Yuanbo Li
!
!------------------------------------------------------------------------------

                         !>> GLOBAL <<!

  USE Constants, ONLY: DP                ! Double precision
  USE TDDFT, ONLY: t_Nstep
  USE TDDFT, ONLY: propagation_method
  USE TDDFT, ONLY: TDE
  USE TDDFT, ONLY: TDE0
  USE TDDFT, ONLY: IonKinetic
  USE TDDFT, ONLY: t_dt
  USE TDDFT, ONLY: t_time
  USE TDDFT, ONLY: totalcharge
  USE Propagation_Ion, ONLY: Initialforces
  USE Propagation_Ion, ONLY: PropagationIon
  USE Propagation_Ion, ONLY: Cleanforces
  USE Propagation_CN,ONLY: PropagationCN
  USE Propagation_RK, ONLY: PropagationRK
  USE Propagation_Split, ONLY: PropagationSplit
  USE Output, ONLY : Print_TDDFT
  USE Sys, ONLY: energy
  USE Sys, ONLY: psi
  USE Sys, ONLY: psi_1
  USE Sys, ONLY: psi_2
  USE Sys, ONLY: velocity 
  USE Sys, ONLY: rhoR
  USE CalPotPlus, ONLY: CalculatePotentialPlus
  USE CellInfo, ONLY: cell
  USE CellInfo, ONLY : m1G, m2G, m3G
  USE Hpsi, ONLY: H_psi
  USE Hpsi, ONLY: V_psi
  USE Fourier, ONLY: FFT
  USE KEDF_TF, ONLY: TFEnergy
  USE RefreshIons, ONLY: RefreshIonTerms
  !USE ReadIonFile, ONLY: ReadDensity
  IMPLICIT NONE

CONTAINS

SUBROUTINE Propagate
  IMPLICIT NONE
  INTEGER(KIND=DP) :: ii
 
  CALL StartClock('Propagate')

  CALL Initialforces

  ! test 2024-6-27
  !CALL ReadDensity("rho_in", 1, 1, 1)
  !psi=SQRT(rhoR(:,:,:,1))  

  DO ii=1,t_Nstep
    WRITE(*,*) "step of TDDFT propagation :", ii
    t_time = ii *t_dt

    ! test
    !WRITE(*,*) "psi"
    !WRITE(*,*) psi
 
    SELECT CASE (propagation_method)
      CASE(0)
        CALL PropagationRK
      CASE(1)
        CALL PropagationCN
      CASE(2)
        CALL PropagationSplit
    END SELECT

    ! test 
    !rhoR(:,:,:,1) = ABS(psi)**2    
    rhoR(:,:,:,1) = ABS(psi_1)**2+ABS(psi_2)**2

    CALL PropagationIon
    CALL RefreshIonTerms
    CALL CalEnergy()
    CALL Print_TDDFT
  END DO 
  CALL Cleanforces

  CALL StopClock('Propagate')
  RETURN
END SUBROUTINE Propagate

SUBROUTINE CalEnergy()
  IMPLICIT NONE
  REAL(KIND=DP), DIMENSION(SIZE(psi, 1),SIZE(psi, 2), SIZE(psi, 3),1):: rhoR
  REAL(KIND=DP), DIMENSION(SIZE(psi, 1),SIZE(psi, 2), SIZE(psi, 3)):: pot  
  REAL(KIND=DP), DIMENSION(SIZE(psi, 1),SIZE(psi, 2), SIZE(psi, 3)):: Hpsi
  REAL(KIND=DP), DIMENSION(SIZE(psi, 1),SIZE(psi, 2), SIZE(psi, 3)):: planewave
  REAL(KIND=DP)::TDEE  ! time dependent electronic energy
  REAL(KIND=DP) :: mass
  INTEGER :: ii, jj , kk

  !WRITE(*,*) "SIZE  OF PSI ", SIZE(psi, 1),SIZE(psi, 2), SIZE(psi, 3)
  rhoR(:,:,:,1) = ABS(psi)**2
  CALL CalculatePotentialPlus(rhoR, .TRUE., pot, energy)
  TDEE = 0.0_DP
  !Hpsi = H_psi(psi)
  !Hpsi = V_psi(psi)
  !DO ii=0,size(psi,1)-1
  !  DO jj=0,size(psi,2)-1
  !    DO kk=0,size(psi,3)-1
  !      TDEE = TDEE + REAL(CONJG(psi(ii,jj,kk))*Hpsi(ii,jj,kk))
        !WRITE(*,*) "i= ",ii," j=",jj," k=",kk," TDEE = ", TDEE
  !    ENDDO 
  !  ENDDO
  !ENDDO
  !WRITE(*,*) "psi 1 5 256 : ",psi(1,5,256)
  !WRITE(*,*) "TDEE  ",TDEE
  !TDE = TDEE !- 2.0_DP/3.0_DP*TFEnergy(rhoR(:,:,:,1))! + energy(6)
  !TDE = energy(3) + energy(4) + energy(5)+ energy(2)+ energy(6)
  TDE = energy(1)
  IonKinetic = 0.0_DP
  DO ii=1,size(cell%ionTable,1)
      mass = cell%elementTable(cell%ionTable(ii)%elementID)%mass
      IonKinetic = IonKinetic + 0.5_DP * mass * SUM(velocity(:,ii)**2)
  ENDDO
  TDE = TDE + IonKinetic

  TDE0=TDE0 + IonKinetic + energy(6)
  totalcharge = SUM(rhoR(:,:,:,1))* cell%vol /m1G /m2G /m3G
  !WRITE(*,'(A,3I8)') " m1G,m2G,m3G      : ", m1G, m2G, m3G  

  !planewave = FFT(psi)
  !WRITE(*,*) "Plane wave coefficient (0,0,0) : ", planewave(0,0,0)
  !WRITE(*,*) "Plane wave coefficient (1,1,1) : ", planewave(1,1,1)  
  !WRITE(*,*) "Plane wave coefficient (1,1,0) : ", planewave(1,1,0)  
  !WRITE(*,*) "Plane wave coefficient (1,0,0) : ", planewave(1,0,0)  
  !WRITE(*,*) "Plane wave coefficient (0,1,0) : ", planewave(0,1,0)

  RETURN
END SUBROUTINE CalEnergy

END MODULE Propagation
