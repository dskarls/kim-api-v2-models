!
! CDDL HEADER START
!
! The contents of this file are subject to the terms of the Common Development
! and Distribution License Version 1.0 (the "License").
!
! You can obtain a copy of the license at
! http://www.opensource.org/licenses/CDDL-1.0.  See the License for the
! specific language governing permissions and limitations under the License.
!
! When distributing Covered Code, include this CDDL HEADER in each file and
! include the License file in a prominent location with the name LICENSE.CDDL.
! If applicable, add the following below this CDDL HEADER, with the fields
! enclosed by brackets "[]" replaced with your own identifying information:
!
! Portions Copyright (c) [yyyy] [name of copyright owner]. All rights reserved.
!
! CDDL HEADER END
!

!
! Contributors:
!    Nick Lewkow
!

!****************************************************************************
!**
!**  MODULE model_He_P_AbIn
!**
!**  Ab Initio Ground State He+He Interaction Potential
!**
!**  Reference: Hellmann, R., Bich, E., & Vogel, E. 2007, Mol. Phys., 105, 3013.
!**
!**  Language: Fortran 2003
!**
!**  Release: This file is part of the kim-api.git repository.
!**
!****************************************************************************
module model_He_P_AbIn

use, intrinsic :: iso_c_binding
use kim_model_headers_module
implicit none

save
private
public Compute_Energy_Forces, &
       refresh, &
       destroy, &
       compute_arguments_create, &
       compute_arguments_destroy, &
       model_cutoff, &
       speccode, &
       buffer_type

! Below are the definitions and values of all Model parameters
integer(c_int), parameter :: cd = c_double ! used for literal constants
integer(c_int), parameter :: DIM=3          ! dimensionality of space
integer(c_int), parameter :: speccode = 1   ! internal species code
real(c_double), parameter :: model_cutoff  = 8.0_cd  ! cutoff radius
                                                     ! in angstroms
real(c_double), parameter :: empir_cutoff  = 0.522562533_cd ! cutoff radius for
                                                            ! empirical/ab initio switch
                                                            ! in angstroms
real(c_double), parameter :: model_cutsq   = model_cutoff**2

!-------------------------------------------------------------------------------
! Below are the definitions and values of all additional model parameters
!
! Recall that the Fortran 2003 format for declaring parameters is as follows:
!
! integer(c_int), parameter :: parname = value   ! This defines an integer
!                                                ! parameter called `parname'
!                                                ! with a value equal to
!                                                ! `value' (a number)
!
! real(c_double), parameter :: parname = value   ! This defines a real(c_double)
!                                                ! parameter called `parname'
!                                                ! with a value equal to
!                                                ! `value' (a number)
!-------------------------------------------------------------------------------
real(c_double),PARAMETER :: TOA   =  0.529177249_cd                 ! Bohr to Angstrom
real(c_double),PARAMETER :: A     =  0.307092338615e7_cd            ! Units K
real(c_double),PARAMETER :: a1    = -0.201651289932e1_cd/TOA        ! Units 1/a0
real(c_double),PARAMETER :: a_1   = -0.431646276045e0_cd*TOA        ! Units a0
real(c_double),PARAMETER :: a2    = -0.459521265125e-1_cd/(TOA*TOA) ! Units 1/a0/a0
real(c_double),PARAMETER :: a_2   =  0.138539045980e0_cd*TOA*TOA    ! Units a0*a0
real(c_double),PARAMETER :: d1    =  0.167127323768e-2_cd           ! Dimensionless
real(c_double),PARAMETER :: d2    =  0.178284243205e1_cd/TOA        ! Units 1/a0
real(c_double),PARAMETER :: d3    =  0.176635702255e1_cd            ! Dimensionless
real(c_double),PARAMETER :: b     =  0.203625105759e1_cd/TOA        ! Units 1/a0
real(c_double),PARAMETER :: ep_kb =  10.997898_cd                   ! Units K
real(c_double),PARAMETER :: R_ep  =  5.608068_cd*TOA                ! Units a0
real(c_double),PARAMETER :: sig   =  4.990672_cd*TOA                ! Units a0
real(c_double),PARAMETER :: TOK   = 315776.177845143_cd             ! AU to K
real(c_double),PARAMETER :: TOEV  = 27.21138386_cd                  ! AU to eV
real(c_double),PARAMETER :: Z     =  2.0_cd
real(c_double),PARAMETER :: aa    =  0.88534_cd/SQRT(2.0_cd*Z**0.23_cd)
real(c_double),PARAMETER :: aa1   = 0.1818_cd
real(c_double),PARAMETER :: aa2   = 0.50990_cd
real(c_double),PARAMETER :: aa3   = 0.2802_cd
real(c_double),PARAMETER :: aa4   = 0.02817_cd
real(c_double),PARAMETER :: bb1   = 3.2_cd
real(c_double),PARAMETER :: bb2   = 0.9423_cd
real(c_double),PARAMETER :: bb3   = 0.4029_cd
real(c_double),PARAMETER :: bb4   = 0.2016_cd

real(c_double),PARAMETER,DIMENSION(17) :: C = (/0.0_cd,0.0_cd,0.0_cd,0.0_cd,0.0_cd, &
                                                0.4616213781e6_cd,  &
                                                0.0_cd,             &
                                                0.4460565781e7_cd,  &
                                                0.0_cd,             &
                                                0.5803352873e8_cd,  &
                                                0.0_cd,             &
                                                0.1031677697e10_cd, &
                                                0.0_cd,             &
                                                0.2415716766e11_cd, &
                                                0.0_cd,             &
                                                0.7191492488e12_cd, &
                                                0.0_cd /)

type, bind(c) :: buffer_type
  real(c_double) :: influence_distance
  real(c_double) :: cutoff(1)
  real(c_double) :: cutsq(1)
end type buffer_type

contains

!-------------------------------------------------------------------------------
!
!  Calculate pair potential phi(r)
!
!-------------------------------------------------------------------------------
subroutine calc_phi(R,phi)
  implicit none

  !-- Transferred variables
  real(c_double), intent(in)  :: R
  real(c_double), intent(out) :: phi

  !-- Local variables
  real(c_double) :: fact
  real(c_double) :: ret0
  real(c_double) :: ret1
  real(c_double) :: ret2
  real(c_double) :: atom

  integer(c_int) :: n
  integer(c_int) :: k
  integer(c_int) :: i

  ret1 = 0.0_cd

  IF (R .gt. model_cutoff) THEN
    ! Argument exceeds cutoff radius
    phi = 0.0_cd

  ELSE IF ( (R .gt. empir_cutoff) .AND. (R .lt. model_cutoff) ) THEN
    ! Use ab initio potential
    ret0 = A*EXP(a1*R + a2*(R*R) + a_1/R + a_2/(R*R) + d1*SIN(d2*R+d3))
    DO n=3,8
      ret2  = 0.0_cd
      DO k=0,(2*n)
        fact = 0.0_cd
        IF ( (k == 0) .OR. (k == 1) ) THEN
          fact = 1.0_cd
        ELSE
          fact = 1.0_cd
          DO i=2,k
            fact = fact*i
          END DO ! i
        END IF ! if (k)
        ret2 = ret2 + (b*R)**k/fact
      END DO ! k
      ret1 = ret1 + (C(2*n)/(R**(2*n)))*(1.0_cd-EXP(-b*R)*ret2)
    END DO ! n
    phi = (ret0 - ret1)/TOK ! [au]

  ELSE
    ! Use screened Coulomb potential
    atom = aa1*EXP(-bb1*(R/aa)) + aa2*EXP(-bb2*(R/aa)) &
           + aa3*EXP(-bb3*(R/aa)) + aa4*EXP(-bb4*(R/aa))
    phi = (Z*Z/R)*atom*10_cd**1.37503_cd/TOEV
  END IF

  ! Convert from au to ev
  phi = phi*TOEV

  return
end subroutine calc_phi

!-------------------------------------------------------------------------------
!
!  Calculate pair potential phi(r) and its derivative dphi(r)
!
!-------------------------------------------------------------------------------
subroutine calc_phi_dphi(R,phi,dphi)
  implicit none

  !-- Transferred variables
  real(c_double), intent(in)  :: R
  real(c_double), intent(out) :: phi,dphi

  !-- Local variables
  real(c_double) :: fact
  real(c_double) :: ret0
  real(c_double) :: ret1
  real(c_double) :: ret2
  real(c_double) :: dret0
  real(c_double) :: dret1
  real(c_double) :: atom
  real(c_double) :: datom
  real(c_double) :: ksum1
  real(c_double) :: ksum2

  integer(c_int) :: n
  integer(c_int) :: k
  integer(c_int) :: i

  ret1 = 0.0_cd

  if (R .gt. model_cutoff) then
    ! Argument exceeds cutoff radius
    phi    = 0.0_cd
    dphi   = 0.0_cd

  else if ( (R .gt. empir_cutoff) .and. (R .lt. model_cutoff) ) then
    ! Use ab initio potential
    ret0  = A*EXP(a1*R+ a2*R*R + a_1/R &
            + a_2/(R*R) + d1*SIN(d2*R+d3))
    DO n=3,8
      ret2  = 0.0_cd
      DO k=0,(2*n)
        fact = 0.0_cd
        IF ( (k == 0) .OR. (k == 1) ) THEN
          fact = 1.0_cd
        ELSE
          fact = 1.0_cd
          DO i=2,k
            fact = fact*i
          END DO ! i
        END IF ! if (k)
        ret2 = ret2 + (b*R)**k/fact
      END DO ! k
      ret1 = ret1 + (C(2*n)/(R**(2*n)))*(1.0_cd-EXP(-b*R)*ret2)
    END DO ! n
    phi = (ret0 - ret1)/TOK ! [au]

    ! Calculate derivative of potential
    dret0 = ret0*(a1+2.0_cd*a2*R-a_1/(R*R)-2.d0*a_2/(R*R*R)+ &
            d1*d2*COS(d2*R+d3))
    dret1 = 0.0_cd
    DO n=3,8
      ksum1 = 0.0_cd
      ksum2 = 0.0_cd
      DO k=0,(2*n)
        fact = 0.0_cd
        IF ( (k.eq.0) .or. (k.eq.1) ) THEN
          fact = 1.0_cd
        ELSE
          fact = 1.0_cd
          DO i=2,k
            fact = fact*i
          END DO ! i
        END IF
        ksum1 = ksum1 + (b*R)**k/fact
        ksum2 = ksum2 + (b**k)*real(k,cd)*(R**(k-1))/fact
       END DO ! k
       dret1 = dret1 - 2.0_cd*real(n,cd)*C(2*n)*(R**(-2*n-1))*  &
               (1.0_cd-EXP(-b*R)*ksum1) + C(2*n)/(R**(2*n))*    &
               (b*EXP(-b*R)*ksum1 - EXP(-b*R)*ksum2)
    END DO ! n
    dphi = (dret0 - dret1)/TOK

  else
    ! Use screened Coulomb potential
    atom  = aa1*EXP(-bb1*(R/aa)) + aa2*EXP(-bb2*(R/aa)) &
            + aa3*EXP(-bb3*(R/aa)) + aa4*EXP(-bb4*(R/aa))
    datom = -(aa1*bb1/aa)*EXP(-bb1*(R/aa)) &
            -(aa2*bb2/aa)*EXP(-bb2*(R/aa)) &
            -(aa3*bb3/aa)*EXP(-bb3*(R/aa)) &
            -(aa4*bb4/aa)*EXP(-bb4*(R/aa))

    phi  =  (Z*Z/R)*atom*10_cd**1.37503_cd/TOEV
    dphi = -(Z*Z/(R*R))*atom + (Z*Z/R)*datom
    dphi = dphi*10_cd**1.37503_cd/TOEV
  endif

  ! Convert to eV
  phi  = phi*TOEV
  dphi = dphi*TOEV

  return
end subroutine calc_phi_dphi

!-------------------------------------------------------------------------------
!
! Compute energy and forces on atoms from the positions.
!
!-------------------------------------------------------------------------------
#include "kim_model_compute_log_macros.fd"
subroutine Compute_Energy_Forces(model_compute_handle, &
  model_compute_arguments_handle, ierr) bind(c)
  implicit none

  !-- Transferred variables
  type(kim_model_compute_handle_type), intent(in) :: model_compute_handle
  type(kim_model_compute_arguments_handle_type), intent(in) :: &
    model_compute_arguments_handle
  integer(c_int), intent(out) :: ierr

  !-- Local variables
  real(c_double) :: Rij(DIM)
  real(c_double) :: r,Rsqij,phi,dphi,dEidr = 0.0_cd
  integer(c_int) :: i,j,jj,numnei,atom_ret,comp_energy,comp_force,comp_particleEnergy, &
                    comp_virial
  integer(c_int) :: ierr2

  !-- KIM variables
  integer(c_int), pointer :: N
  real(c_double), pointer :: energy
  real(c_double), pointer :: coor(:,:)
  real(c_double), pointer :: force(:,:)
  real(c_double), pointer :: particleEnergy(:)
  integer(c_int), pointer :: nei1part(:)
  integer(c_int), pointer :: particleSpeciesCodes(:)
  integer(c_int), pointer :: particleContributing(:)
  real(c_double), pointer :: virial(:)

  kim_log_file = __FILE__

  ! Unpack data from KIM object
  !
  ierr = 0
  call kim_model_compute_arguments_get_argument_pointer( &
    model_compute_arguments_handle, &
    kim_compute_argument_name_number_of_particles, N, ierr2)
  ierr = ierr + ierr2
  call kim_model_compute_arguments_get_argument_pointer( &
    model_compute_arguments_handle, &
    kim_compute_argument_name_particle_species_codes, n, particleSpeciesCodes, &
    ierr2)
  ierr = ierr + ierr2
  call kim_model_compute_arguments_get_argument_pointer( &
    model_compute_arguments_handle, &
    kim_compute_argument_name_particle_contributing, n, particleContributing, &
    ierr2)
  ierr = ierr + ierr2
  call kim_model_compute_arguments_get_argument_pointer( &
    model_compute_arguments_handle, &
    kim_compute_argument_name_coordinates, dim, n, coor, ierr2)
  ierr = ierr + ierr2
  call kim_model_compute_arguments_get_argument_pointer( &
    model_compute_arguments_handle, &
    kim_compute_argument_name_partial_energy, energy, ierr2)
  ierr = ierr + ierr2
  call kim_model_compute_arguments_get_argument_pointer( &
    model_compute_arguments_handle, &
    kim_compute_argument_name_partial_forces, dim, n, force, ierr2)
  ierr = ierr + ierr2
  call kim_model_compute_arguments_get_argument_pointer( &
    model_compute_arguments_handle, &
    kim_compute_argument_name_partial_particle_energy, n, particleEnergy, ierr2)
  ierr = ierr + ierr2
  call kim_model_compute_arguments_get_argument_pointer( &
    model_compute_arguments_handle, &
    kim_compute_argument_name_partial_virial, 6, virial, ierr2)
  ierr = ierr + ierr2
  if (ierr /= 0) then
    kim_log_message = "Failed to retrieve data from KIM API compute-arguments object"
    LOG_ERROR()
    return
  endif

  ! Check to see if we have been asked to compute the forces, energyperpart,
  ! energy and virial
  !
  if (associated(energy)) then
    comp_energy = 1
  else
    comp_energy = 0
  end if
  if (associated(force)) then
    comp_force = 1
  else
    comp_force = 0
  end if
  if (associated(particleEnergy)) then
    comp_particleEnergy = 1
  else
    comp_particleEnergy = 0
  end if
  if (associated(virial)) then
    comp_virial = 1
  else
    comp_virial = 0
  end if

  ! Check to be sure that the species are correct
  !
  ierr = 1 ! assume an error
  do i = 1,N
    if (particleSpeciesCodes(i).ne.speccode) then
      kim_log_message = "Unexpected species code detected"
      LOG_ERROR()
      return
    endif
  enddo
  ierr = 0 ! everything is ok

  ! Initialize potential energies, forces, virial term, and electron density
  !
  ! Note: that the variable `particleEnergy' does not need to be initialized
  !       because it's initial value is set during the embedding energy
  !       calculation.
  if (comp_energy.eq.1) energy = 0.0_cd
  if (comp_particleEnergy.eq.1) particleEnergy = 0.0_cd
  if (comp_force.eq.1) force  = 0.0_cd
  if (comp_virial.eq.1) virial = 0.0_cd

  !
  !  Compute energy and forces
  !

  !  Loop over particles and compute energy and forces
  !
  do i = 1,N
    if(particleContributing(i).eq.1) then
      ! Get neighbor list of current atom
      call kim_model_compute_arguments_get_neighbor_list( &
        model_compute_arguments_handle, 1, i, numnei, nei1part, ierr)
      if (ierr /= 0) then
        ! some sort of problem, exit
        kim_log_message = "GetNeighborList failed"
        LOG_ERROR()
        ierr = 1
        return
      endif

      ! Loop over the neighbors of atom i
      !
      do jj = 1, numnei

        j = nei1part(jj)                         ! get neighbor ID

        ! compute relative position vector
        !
        Rij(:) = coor(:,j) - coor(:,i)           ! distance vector between i j

        ! compute energy and forces
        !
        Rsqij = dot_product(Rij,Rij)             ! compute square distance
        if ( Rsqij .lt. model_cutsq ) then       ! particles are interacting?

          r = sqrt(Rsqij)                        ! compute distance
          if (comp_force.eq.1.or.comp_virial.eq.1) then
            call calc_phi_dphi(r,phi,dphi)       ! compute pair potential
                                                 !   and it derivative
            dEidr = 0.5_cd*dphi
          else
            call calc_phi(r,phi)                 ! compute just pair potential
          endif

          ! contribution to energy
          !
          if (comp_particleEnergy.eq.1) then
            particleEnergy(i) = particleEnergy(i) + 0.5_cd*phi   ! accumulate energy
          endif
          if (comp_energy.eq.1) then
            energy = energy + 0.5_cd*phi         !      add half v to total energy
          endif

          ! contribution to virial tensor, virial(i,j)=r(i)*r(j)*(dV/dr)/r
          !
          if (comp_virial.eq.1) then
            virial(1) = virial(1) + Rij(1)*Rij(1)*dEidr/r
            virial(2) = virial(2) + Rij(2)*Rij(2)*dEidr/r
            virial(3) = virial(3) + Rij(3)*Rij(3)*dEidr/r
            virial(4) = virial(4) + Rij(2)*Rij(3)*dEidr/r
            virial(5) = virial(5) + Rij(1)*Rij(3)*dEidr/r
            virial(6) = virial(6) + Rij(1)*Rij(2)*dEidr/r
          endif

          ! contribution to forces
          !
          if (comp_force.eq.1) then
            force(:,i) = force(:,i) + dEidr*Rij/r ! accumulate force on atom i
            force(:,j) = force(:,j) - dEidr*Rij/r ! accumulate force on atom j
          endif
        endif

      enddo  ! Loop on neighbors of atom i
    endif ! Check on whether particle is contributing
  enddo  ! Loop of atom i

  ! Everything is great
  !
  ierr = 0
  return

end subroutine Compute_Energy_Forces

!-------------------------------------------------------------------------------
!
! Model destroy routine
!
!-------------------------------------------------------------------------------
#include "kim_model_destroy_log_macros.fd"
subroutine destroy(model_destroy_handle, ierr) bind(c)
  use, intrinsic :: iso_c_binding
  implicit none

  !-- Transferred variables
  type(kim_model_destroy_handle_type), intent(inout) :: model_destroy_handle
  integer(c_int), intent(out) :: ierr
  type(buffer_type), pointer :: buf; type(c_ptr) :: pbuf

  kim_log_file = __FILE__

  call kim_model_destroy_get_model_buffer_pointer(model_destroy_handle, pbuf)
  call c_f_pointer(pbuf, buf)

  kim_log_message = "deallocating model buffer"
  LOG_INFORMATION()
  deallocate(buf)
  ierr = 0  ! everything is good

  return
end subroutine destroy

!-------------------------------------------------------------------------------
!
! Model refresh routine (REQUIRED)
!
!-------------------------------------------------------------------------------
#include "kim_model_refresh_log_macros.fd"
subroutine refresh(model_refresh_handle, ierr) bind(c)
  use, intrinsic :: iso_c_binding
  implicit none

  !-- Transferred variables
  type(kim_model_refresh_handle_type), intent(inout) :: model_refresh_handle
  integer(c_int), intent(out) :: ierr

  type(buffer_type), pointer :: buf; type(c_ptr) :: pbuf

  kim_log_file = __FILE__

  call kim_model_refresh_get_model_buffer_pointer(model_refresh_handle, pbuf)
  call c_f_pointer(pbuf, buf)

  kim_log_message = "Resetting influence distance and cutoffs"
  LOG_INFORMATION()
  call kim_model_refresh_set_influence_distance_pointer( &
    model_refresh_handle, buf%influence_distance)
  call kim_model_refresh_set_neighbor_list_cutoffs_pointer( &
    model_refresh_handle, 1, buf%cutoff)

  ierr = 0  ! everything is good

  return
end subroutine refresh

!-------------------------------------------------------------------------------
!
! Model compute arguments create routine (REQUIRED)
!
!-------------------------------------------------------------------------------
#include "kim_model_compute_arguments_create_log_macros.fd"
subroutine compute_arguments_create(model_compute_handle, &
  model_compute_arguments_create_handle, ierr) bind(c)
  use, intrinsic :: iso_c_binding
  use kim_model_compute_arguments_create_module, &
    log_entry=>kim_model_compute_arguments_create_log_entry
  implicit none

  !-- Transferred variables
  type(kim_model_compute_handle_type), intent(in) :: model_compute_handle
  type(kim_model_compute_arguments_create_handle_type), intent(inout) :: &
    model_compute_arguments_create_handle
  integer(c_int), intent(out) :: ierr

  integer(c_int) :: ierr2

  ierr = 0
  ierr2 = 0

  ! register arguments
  call kim_model_compute_arguments_create_set_argument_support_status( &
    model_compute_arguments_create_handle, &
    kim_compute_argument_name_partial_energy, &
    kim_support_status_optional, ierr2)
  ierr = ierr + ierr2
  call kim_model_compute_arguments_create_set_argument_support_status( &
    model_compute_arguments_create_handle, &
    kim_compute_argument_name_partial_forces, &
    kim_support_status_optional, ierr2)
  ierr = ierr + ierr2
  call kim_model_compute_arguments_create_set_argument_support_status( &
    model_compute_arguments_create_handle, &
    kim_compute_argument_name_partial_particle_energy, &
    kim_support_status_optional, ierr2)
  ierr = ierr + ierr2
  call kim_model_compute_arguments_create_set_argument_support_status( &
    model_compute_arguments_create_handle, &
    kim_compute_argument_name_partial_virial, &
    kim_support_status_optional, ierr2)
  ierr = ierr + ierr2

  ! register call backs (ProcessDEDrTerm, ProcessD2EDr2Term
  ! NONE

  if (ierr /= 0) then
    ierr = 1
    kim_log_message = "Unable to successfully create compute_arguments object"
    LOG_ERROR()
  endif

  return
end subroutine compute_arguments_create

!-------------------------------------------------------------------------------
!
! Model compute arguments destroy routine (REQUIRED)
!
!-------------------------------------------------------------------------------
#include "kim_model_compute_arguments_destroy_log_macros.fd"
subroutine compute_arguments_destroy(model_compute_handle, &
  model_compute_arguments_destroy_handle, ierr) bind(c)
  use, intrinsic :: iso_c_binding
  use kim_model_compute_arguments_destroy_module, &
    log_entry=>kim_model_compute_arguments_destroy_log_entry
  implicit none

  !-- Transferred variables
  type(kim_model_compute_handle_type), intent(in) :: model_compute_handle
  type(kim_model_compute_arguments_destroy_handle_type), intent(inout) :: &
    model_compute_arguments_destroy_handle
  integer(c_int), intent(out) :: ierr

  integer(c_int) :: ierr2

  ierr = 0
  ierr2 = 0

  ! nothing to do

  return
end subroutine compute_arguments_destroy

end module model_He_P_AbIn

!-------------------------------------------------------------------------------
!
! Model initialization routine (REQUIRED)
!
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!
! Model initialization routine (REQUIRED)
!
!-------------------------------------------------------------------------------
#include "kim_model_create_log_macros.fd"
subroutine create(model_create_handle, requested_length_unit, &
  requested_energy_unit, requested_charge_unit, requested_temperature_unit, &
  requested_time_unit, ierr) bind(c)
  use, intrinsic :: iso_c_binding
  use model_He_P_AbIn
  use kim_model_headers_module
  implicit none

  !-- Transferred variables
  type(kim_model_create_handle_type), intent(inout) :: model_create_handle
  type(kim_length_unit_type), intent(in) :: requested_length_unit
  type(kim_energy_unit_type), intent(in) :: requested_energy_unit
  type(kim_charge_unit_type), intent(in) :: requested_charge_unit
  type(kim_temperature_unit_type), intent(in) :: requested_temperature_unit
  type(kim_time_unit_type), intent(in) :: requested_time_unit
  integer(c_int), intent(out) :: ierr

  !-- Local variables
  integer(c_int), parameter :: one=1

  !-- KIM variables
  integer(c_int) :: ierr2
  type(buffer_type), pointer :: buf

  kim_log_file = __FILE__

  ierr = 0
  ierr2 = 0

  call kim_model_create_set_units(model_create_handle, &
    kim_length_unit_a, &
    kim_energy_unit_ev, &
    kim_charge_unit_unused, &
    kim_temperature_unit_unused, &
    kim_time_unit_unused, &
    ierr2)
  ierr = ierr + ierr2

  ! register species
  call kim_model_create_set_species_code(model_create_handle, &
    kim_species_name_he, speccode, ierr2)
  ierr = ierr + ierr2

  ! register numbering
  call kim_model_create_set_model_numbering(model_create_handle, &
    kim_numbering_one_based, ierr2);
  ierr = ierr + ierr2

  ! register function pointers
  call kim_model_create_set_compute_pointer(model_create_handle, &
    kim_language_name_fortran, c_funloc(Compute_Energy_Forces), ierr2)
  ierr = ierr + ierr2
  call kim_model_create_set_compute_arguments_create_pointer( &
    model_create_handle, kim_language_name_fortran, &
    c_funloc(compute_arguments_create), ierr2)
  ierr = ierr + ierr2
  call kim_model_create_set_compute_arguments_destroy_pointer( &
    model_create_handle, kim_language_name_fortran, &
    c_funloc(compute_arguments_destroy), ierr2)
  ierr = ierr + ierr2
  call kim_model_create_set_destroy_pointer(model_create_handle, &
    kim_language_name_fortran, c_funloc(destroy), ierr2)
  ierr = ierr + ierr2
  call kim_model_create_set_refresh_pointer( &
    model_create_handle, kim_language_name_fortran, &
    c_funloc(refresh), ierr2)
  ierr = ierr + ierr2

  ! allocate buffer
  allocate(buf)

  ! store model buffer in KIM object
  call kim_model_create_set_model_buffer_pointer(model_create_handle, &
    c_loc(buf))

  ! set buffer values
  buf%influence_distance = model_cutoff
  buf%cutoff(1) = model_cutoff

  ! register influence distance
  call kim_model_create_set_influence_distance_pointer( &
    model_create_handle, buf%influence_distance)

  ! register cutoff
  call kim_model_create_set_neighbor_list_cutoffs_pointer(model_create_handle, &
    1, buf%cutoff)

  if (ierr /= 0) then
    ierr = 1
    kim_log_message = "Unable to successfully initialize model"
    LOG_ERROR()
  endif

  return
end subroutine create
