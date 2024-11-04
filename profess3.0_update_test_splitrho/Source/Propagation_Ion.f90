MODULE Propagation_Ion
!------------------------------------------------------------------------------
! STRUCTURE OF MODULE:
!   MODULE Propagation_Ion
!     |_SUBROUTINE 
!     |_SUBROUTINE 
!
! DESCRIPTION:
!       This module calculates the movement of the ions.
!
! CONDITIONS AND ASSUMPTIONS:
!
!   03/05/2021 Yuanbo Li
!
!------------------------------------------------------------------------------

                         !>> GLOBAL <<!

  USE Constants, ONLY: DP                ! Double precision
  USE TDDFT, ONLY: t_dt  
  USE Sys, ONLY: psi
  USE Sys, ONLY: psi_1
  USE Sys, ONLY: psi_2
  USE Sys, ONLY: forces 
  USE Sys, ONLY: velocity
  USE Sys, ONLY: frozenIon
  USE CellInfo, ONLY: cell
  USE CellInfo, ONLY: numSpin
  USE CalForces, ONLY : CalculateForces
  USE MathFunctions, ONLY: Vecmul
  USE MathFunctions, ONLY : Inverse
  USE RefreshIons, ONLY: RefreshIonTerms
  USE RefreshIons, ONLY: RescaleDensity

  IMPLICIT NONE
  REAL(KIND=DP), DIMENSION(:,:,:), ALLOCATABLE :: forces1 
  ! forces f(t+\delta t)    forces1  f(t)
  ! First is ion number, second direction (1,2,3 for x,y,z), final index is 1 for total force
 
CONTAINS

SUBROUTINE Initialforces
  IMPLICIT NONE
  REAL(KIND=DP), DIMENSION(SIZE(psi, 1),SIZE(psi, 2), SIZE(psi, 3),numSpin) :: rhoR
  
  ALLOCATE(forces1(size(cell%ionTable,1),3,3))
  rhoR(:,:,:,1) = ABS(psi_1)**2+ABS(psi_2)**2
  CALL CalculateForces(rhoR,forces1)
!  WRITE(*,*) "initialize forces"
  RETURN
END SUBROUTINE Initialforces

SUBROUTINE PropagationIon
  IMPLICIT NONE
  REAL(KIND=DP) :: mass ! mass of ion  
  REAL(KIND=DP), DIMENSION(3) :: fracStep ! Step (due to velocity) in fractional coordinates
  INTEGER(KIND=DP) :: ii,jj
  REAL(KIND=DP), DIMENSION(SIZE(psi, 1),SIZE(psi, 2), SIZE(psi, 3),numSpin) :: rhoR

  CALL StartClock('Propagation_Ion')  
 
!  WRITE(*,*) "Propagation_Ion start"
  rhoR(:,:,:,1) = ABS(psi_1)**2+ABS(psi_2)**2
  CALL RescaleDensity(rhoR)
  CALL CalculateForces(rhoR,forces)
  DO ii = 1, size(cell%ionTable,1)
      mass = cell%elementTable(cell%ionTable(ii)%elementID)%mass
      fracStep = Vecmul(Inverse(cell%cellReal), &
      ! (velocity(:,ii)*t_dt+forces1(ii,:,1)*t_dt**2/(2.0_DP*mass)))
        (velocity(:,ii)*t_dt))
      DO jj = 1, 3
        IF (frozenIon(ii,jj) .EQV. .TRUE. ) THEN
          fracStep(jj) = 0.0_DP
        END IF
      END DO
      cell%ionTable(ii)%coord = MODULO(cell%ionTable(ii)%coord+fracStep, 1._DP)
  END DO
  DO ii=1,size(cell%ionTable,1)
      mass = cell%elementTable(cell%ionTable(ii)%elementID)%mass
      DO jj= 1, 3
        IF (frozenIon(ii,jj) .EQV. .FALSE. ) THEN
          velocity(:,ii) = velocity(:,ii) + & 
          !  (forces(ii,:,1)+forces1(ii,:,1))/mass*t_dt/2._DP
          forces(ii,:,1)/mass*t_dt
        END IF
      END DO
  ENDDO
  forces1 = forces
!  WRITE(*,*) "RefreshIonTerms"
 ! CALL RefreshIonTerms
!  WRITE(*,*) "Propagation_Ion end"
  
  CALL StopClock('Propagation_Ion')
  RETURN
END SUBROUTINE PropagationIon

SUBROUTINE Cleanforces
  DEALLOCATE(forces1)
  RETURN
END SUBROUTINE Cleanforces

END MODULE Propagation_Ion
