!****************************************************************************
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
! Copyright (c) 2012, Regents of the University of Minnesota.  All rights reserved.
!
!****************************************************************************
!
! MODULE Pair_Exp6_Hogervorst
!
!
! Compute energy and forces on an isolated cluster of atoms using
! the exp-6 potential shifted to have zero energy at cutoff.
! Parameters due to Hogervorst [1] and mixing rule between species
! due to Kong and Chakrabarty [2].
!
! Author : E. B. Tadmor (10-MAR-12), updated to kim api v1.6 4-AUG-14,
!          updated to kim api v2.0 30-JUN-18 by D. S. Karls.
!
! Support atom species:
!
!  1 = "Ar"
!  2 = "Ne"
!
!  References:
!
!  [1] W. Hogervorst, Physica, vol 51, 77-89 (1971). (Like interactions.)
!
!  [2] C. L. Kong and M. R. Chakrabarty, J. Phys. Chem., Vol. 77, 2668-2670
!      (1973). (Mixing rule.)
!
!****************************************************************************

#define THIS_FILE_NAME __FILE__

module Pair_Exp6_Hogervorst

use, intrinsic :: iso_c_binding
use kim_model_headers_module
implicit none

save
private
public Compute_Energy_Forces, &
       compute_arguments_create, &
       compute_arguments_destroy, &
       refresh, &
       destroy, &
       speccodeAr, &
       speccodeNe, &
       cutoff, &
       buffer_type

integer(c_int), parameter :: cd = c_double ! used for literal constants

!
! Define global potential parameters
!
real(c_double), parameter :: kb     = 8.6173e-5_cd ! Boltzmann's constant [eV/K]
real(c_double), parameter :: cutoff = 8.15_cd      ! cutoff radius [A] (arbitrary value)
real(c_double), parameter :: eps11  = 138.0_cd*kb  ! Ar-Ar epsilon parameter [eV]
real(c_double), parameter :: rm11   = 3.77_cd      ! Ar-Ar rm parameter [A]
real(c_double), parameter :: alf11  = 14.8_cd      ! Ar-Ar alpha parameter [-]
real(c_double), parameter :: eps22  = 43.0_cd*kb   ! Ne-Ne epsilon parameter [eV]
real(c_double), parameter :: rm22   = 3.03_cd      ! Ne-Ne rm parameter [A]
real(c_double), parameter :: alf22  = 16.0_cd      ! Ne-Ne alpha parameter [-]
real(c_double), parameter :: eps12  = 68.89_cd*kb  ! Ar-Ne epsilon parameter [eV]
real(c_double), parameter :: rm12   = 3.447_cd     ! Ar-Ne rm parameter [A]
real(c_double), parameter :: alf12  = 15.52_cd     ! Ar-Ne alpha parameter [-]
integer(c_int), parameter :: speccodeAr = 1, speccodeNe = 2

type, bind(c) :: buffer_type
  real(c_double) :: influence_distance
  real(c_double) :: cutoff(1)
  real(c_double) :: cutsq(1)
  integer(c_int) :: padding_neighbor_hints(1)
  integer(c_int) :: half_list_hints(1)
end type buffer_type

contains

!-------------------------------------------------------------------------------
!
!  Calculate exp-6 pair potential phi(r) and its derivative dphi(r)
!  Note: potential is shifted to have zero energy at r=cutoff
!
!-------------------------------------------------------------------------------
subroutine calc_phi_dphi(r,phi,dphi,eps,rm,alf,cutoff)
  implicit none

  !-- Transferred variables
  real(c_double), intent(in)  :: r
  real(c_double), intent(in)  :: eps,rm,alf,cutoff
  real(c_double), intent(out) :: phi,dphi

  !-- Local variables
  real(c_double) soa,amp,ror,cor,phicut

  if (r .gt. cutoff) then
    ! Argument exceeds cutoff radius
    phi = 0.0_cd
    dphi = 0.0_cd
  else
    soa = 6.0_cd/alf
    amp = eps/(1.0_cd-soa)
    cor = cutoff/rm
    ror = r/rm
    phicut = amp*(soa*exp(alf*(1.0_cd-cor))-cor**(-6))
    phi = amp*(soa*exp(alf*(1.0_cd-ror))-ror**(-6)) - phicut
    dphi = -6.0_cd*amp/rm*(exp(alf*(1.0_cd-ror))-ror**(-7))
  endif

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

  !--Local Variables
  integer(c_int), parameter :: DIM=3
  integer(c_int) :: i, j, jj, numnei, ierr2
  integer(c_int) :: comp_energy, comp_force
  real(c_double) Rij(DIM), Rsqij, r, phi, dphi

  !-- KIM variables
  integer(c_int), pointer :: N
  real(c_double), pointer :: energy
  real(c_double), pointer :: coor(:,:)
  real(c_double), pointer :: force(:,:)
  integer(c_int), pointer :: nei1part(:)
  integer(c_int), pointer :: particleSpeciesCodes(:)
  integer(c_int), pointer :: particleContributing(:)
  type(buffer_type), pointer :: buf; type(c_ptr) :: pbuf

  kim_log_file = __FILE__

  ! Get buffer address from KIM API
  call kim_model_compute_get_model_buffer_pointer(model_compute_handle, pbuf)
  call c_f_pointer(pbuf, buf)

  ierr = 0
  ierr2 = 0

  ! Unpack data from KIM object
  !
  ierr = 0
  call kim_model_compute_arguments_get_argument_pointer( &
    model_compute_arguments_handle, &
    kim_compute_argument_name_number_of_particles, N, ierr2)
  ierr = ierr + ierr2
  call kim_model_compute_arguments_get_argument_pointer( &
    model_compute_arguments_handle, &
    kim_compute_argument_name_particle_species_codes, N, particleSpeciesCodes, &
    ierr2)
  ierr = ierr + ierr2
  call kim_model_compute_arguments_get_argument_pointer( &
    model_compute_arguments_handle, &
    kim_compute_argument_name_particle_contributing, N, particleContributing, &
    ierr2)
  ierr = ierr + ierr2
  call kim_model_compute_arguments_get_argument_pointer( &
    model_compute_arguments_handle, &
    kim_compute_argument_name_coordinates, dim, N, coor, ierr2)
  ierr = ierr + ierr2
  call kim_model_compute_arguments_get_argument_pointer( &
    model_compute_arguments_handle, &
    kim_compute_argument_name_partial_energy, energy, ierr2)
  ierr = ierr + ierr2
  call kim_model_compute_arguments_get_argument_pointer( &
    model_compute_arguments_handle, &
    kim_compute_argument_name_partial_forces, dim, N, force, ierr2)
  ierr = ierr + ierr2
  if (ierr /= 0) then
    kim_log_message = "Failed to retrieve data from KIM API compute-arguments object"
    LOG_ERROR()
    return
  endif

  ! Check to see if we have been asked to compute the energy and/or forces
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

  ! Check to be sure that the species are correct
  !
  ierr = 1 ! assume an error
  do i = 1,N
    if (particleSpeciesCodes(i).ne.speccodeAr .and. particleSpeciesCodes(i).ne.speccodeNe) then
      kim_log_message = "Unexpected species code detected"
      LOG_ERROR()
      return
    endif
  enddo
  ierr = 0 ! everything is ok

  ! Compute energy and forces
  if (comp_energy.eq.1) energy = 0.0_cd
  if (comp_force.eq.1)  force  = 0.0_cd

  do i=1,N
    if (particleContributing(i) == 1) then
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

      do jj = 1, numnei

        j = nei1part(jj) ! Get index of neighbor

        if ( i.lt.j ) then
          Rij(:) = coor(:,j) - coor(:,i)
          Rsqij = dot_product(Rij,Rij)

          if ( Rsqij.lt.buf%cutsq(1) ) then
             r = sqrt(Rsqij)
             if (particleSpeciesCodes(i).eq.1 .and. particleSpeciesCodes(j).eq.1) then
                call calc_phi_dphi(r,phi,dphi,eps11,rm11,alf11,buf%cutoff(1))
             elseif (particleSpeciesCodes(i).eq.2 .and. particleSpeciesCodes(j).eq.2) then
                call calc_phi_dphi(r,phi,dphi,eps22,rm22,alf22,buf%cutoff(1))
             else
                call calc_phi_dphi(r,phi,dphi,eps12,rm12,alf12,buf%cutoff(1))
             endif
             if (comp_energy.eq.1) then
               if (particleContributing(j) == 1) then
                 energy = energy + phi
               else
                 energy = energy + 0.5_cd*phi
               endif
             endif

             if (comp_force.eq.1) then
                force(:,i) = force(:,i) + dphi*Rij/r
                force(:,j) = force(:,j) - dphi*Rij/r
             endif

          endif
        endif
     enddo ! End loop over neighbors of atom i
    endif ! Check on whether particle is contributing
enddo ! End primary loop over particles

! No errors
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
  call kim_model_refresh_set_neighbor_list_pointers( &
    model_refresh_handle, 1, buf%cutoff, buf%padding_neighbor_hints, buf%half_list_hints)

  ierr = 0  ! everything is good
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

  ierr = 0

  ! nothing to do

  return
end subroutine compute_arguments_destroy

end module Pair_Exp6_Hogervorst

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
  use Pair_Exp6_Hogervorst
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
    kim_species_name_ar, speccodeAr, ierr2)
  ierr = ierr + ierr2
  call kim_model_create_set_species_code(model_create_handle, &
    kim_species_name_ne, speccodeNe, ierr2)
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
  buf%influence_distance = cutoff
  buf%cutoff(1) = cutoff
  buf%cutsq(1) = cutoff*cutoff
  buf%padding_neighbor_hints(1) = 1
  buf%half_list_hints(1) = 1

  ! register influence distance
  call kim_model_create_set_influence_distance_pointer( &
    model_create_handle, buf%influence_distance)

  ! register cutoff
  call kim_model_create_set_neighbor_list_pointers(model_create_handle, &
    1, buf%cutoff, buf%padding_neighbor_hints, buf%half_list_hints)

  if (ierr /= 0) then
    ierr = 1
    kim_log_message = "Unable to successfully initialize model"
    LOG_ERROR()
  endif
end subroutine create
