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
! Copyright (c) 2013, Regents of the University of Minnesota.
! All rights reserved.
!
! Contributors:
!    Ellad B. Tadmor
!    Ryan S. Elliott
!    Stephen M. Whalen
!

!****************************************************************************
!**
!**  MODULE Pair_Lennard_Jones_Truncated
!**
!**  Lennard-Jones pair potential KIM Model Driver
!**  truncated to have zero energy above the cutoff radius
!**
!**  Note: The energy and its derivatives will be discontinuous at the cutoff.
!**
!**  Language: Fortran 2003
!**
!****************************************************************************
module Pair_Lennard_Jones_Truncated

use, intrinsic :: iso_c_binding
use kim_model_driver_headers_module
implicit none

save
private
public BUFFER_TYPE,               &
       Compute_Energy_Forces,     &
       compute_arguments_create,  &
       compute_arguments_destroy, &
       refresh   ,                &
       destroy,                   &
       calc_phi,                  &
       calc_phi_dphi,             &
       calc_phi_dphi_d2phi,       &
       speccode

! Below are the definitions and values of all Model parameters
integer(c_int), parameter          :: cd = c_double  ! for literal constants
integer(c_int), parameter          :: DIM=3          ! dimensionality of space
integer(c_int), parameter          :: speccode = 1   ! internal species code

!-------------------------------------------------------------------------------
!
!  Definition of Buffer type
!
!-------------------------------------------------------------------------------
type BUFFER_TYPE
  real(c_double) :: influence_distance
  real(c_double) :: cutoff(1)
  real(c_double) :: cutsq(1)
  real(c_double) :: eps(1)
  real(c_double) :: sigma(1)
  integer(c_int) :: padding_neighbor_hints(1)
  integer(c_int) :: half_list_hints(1)
endtype BUFFER_TYPE


contains

!-------------------------------------------------------------------------------
!
!  Calculate pair potential phi(r)
!
!-------------------------------------------------------------------------------
subroutine calc_phi(model_eps,  &
                    model_sigma,    &
                    model_cutoff,r,phi)
  implicit none

  !-- Transferred variables
  real(c_double), intent(in)  :: model_eps
  real(c_double), intent(in)  :: model_sigma
  real(c_double), intent(in)  :: model_cutoff
  real(c_double), intent(in)  :: r
  real(c_double), intent(out) :: phi

  !-- Local variables
  real(c_double) rsq,sor,sor6,sor12

  rsq  = r*r             !  r^2
  sor  = model_sigma/r   !  (sig/r)
  sor6 = sor*sor*sor     !
  sor6 = sor6*sor6       !  (sig/r)^6
  sor12= sor6*sor6       !  (sig/r)^12
  if (r .gt. model_cutoff) then
     ! Argument exceeds cutoff radius
     phi = 0.0_cd
  else
     phi = 4.0_cd*model_eps*(sor12-sor6)
  endif

  return
end subroutine calc_phi

!-------------------------------------------------------------------------------
!
!  Calculate pair potential phi(r) and its derivative dphi(r)
!
!-------------------------------------------------------------------------------
subroutine calc_phi_dphi(model_eps,  &
                         model_sigma,    &
                         model_cutoff,r,phi,dphi)
  implicit none

  !-- Transferred variables
  real(c_double), intent(in)  :: model_eps
  real(c_double), intent(in)  :: model_sigma
  real(c_double), intent(in)  :: model_cutoff
  real(c_double), intent(in)  :: r
  real(c_double), intent(out) :: phi,dphi

  !-- Local variables
  real(c_double) rsq,sor,sor6,sor12

  rsq  = r*r             !  r^2
  sor  = model_sigma/r   !  (sig/r)
  sor6 = sor*sor*sor     !
  sor6 = sor6*sor6       !  (sig/r)^6
  sor12= sor6*sor6       !  (sig/r)^12
  if (r .gt. model_cutoff) then
     ! Argument exceeds cutoff radius
     phi    = 0.0_cd
     dphi   = 0.0_cd
  else
     phi  = 4.0_cd*model_eps*(sor12-sor6)
     dphi = 24.0_cd*model_eps*(-2.0_cd*sor12+sor6)/r
  endif

  return
end subroutine calc_phi_dphi

!-------------------------------------------------------------------------------
!
!  Calculate pair potential phi(r) and its derivatives dphi(r) and d2phi(r)
!
!-------------------------------------------------------------------------------
subroutine calc_phi_dphi_d2phi(model_eps,  &
                               model_sigma,    &
                               model_cutoff,r,phi,dphi,d2phi)
  implicit none

  !-- Transferred variables
  real(c_double), intent(in)  :: model_eps
  real(c_double), intent(in)  :: model_sigma
  real(c_double), intent(in)  :: model_cutoff
  real(c_double), intent(in)  :: r
  real(c_double), intent(out) :: phi,dphi,d2phi

  !-- Local variables
  real(c_double) rsq,sor,sor6,sor12

  rsq  = r*r             !  r^2
  sor  = model_sigma/r   !  (sig/r)
  sor6 = sor*sor*sor     !
  sor6 = sor6*sor6       !  (sig/r)^6
  sor12= sor6*sor6       !  (sig/r)^12
  if (r .gt. model_cutoff) then
     ! Argument exceeds cutoff radius
     phi    = 0.0_cd
     dphi   = 0.0_cd
     d2phi  = 0.0_cd
  else
     phi   = 4.0_cd*model_eps*(sor12-sor6)
     dphi  = 24.0_cd*model_eps*(-2.0_cd*sor12+sor6)/r
     d2phi = 24.0_cd*model_eps*(26.0_cd*sor12-7.0_cd*sor6)/rsq
  endif

  return
end subroutine calc_phi_dphi_d2phi

!-------------------------------------------------------------------------------
!
! Compute energy and forces on atoms from the positions.
!
!-------------------------------------------------------------------------------i
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
  real(c_double) :: r,Rsqij,phi,dphi,d2phi,dEidr,d2Eidr
  integer(c_int) :: i,j,jj,numnei
  integer(c_int) :: ierr2
  integer(c_int) :: comp_force,comp_energy,comp_particleEnergy,comp_process_dEdr, &
                    comp_process_d2Edr2
  type(BUFFER_TYPE), pointer :: buf; type(c_ptr) :: pbuf

  real(c_double), pointer :: Rij(:)
  ! Quantities for computing d2Edr2
  real(c_double), pointer :: Rij_pairs(:,:)
  real(c_double), pointer :: r_pairs(:)
  integer(c_int), pointer :: i_pairs(:), j_pairs(:)

  !-- KIM variables
  real(c_double) :: model_cutoff
  integer(c_int), pointer :: N
  real(c_double), pointer :: energy
  real(c_double), pointer :: coor(:,:)
  real(c_double), pointer :: force(:,:)
  real(c_double), pointer :: particleEnergy(:)
  integer(c_int), pointer :: nei1part(:)
  integer(c_int), pointer :: particleSpeciesCodes(:)
  integer(c_int), pointer :: particleContributing(:)

  kim_log_file = __FILE__

  ! get model buffer from KIM object
  call kim_model_compute_get_model_buffer_pointer(model_compute_handle, pbuf)
  call c_f_pointer(pbuf, buf)

  model_cutoff = buf%influence_distance

  ! Unpack data from KIM object
  !
  ierr = 0
  call kim_model_compute_arguments_get_argument_pointer( &
    model_compute_arguments_handle, &
    kim_compute_argument_name_number_of_particles, N, ierr2)
  ierr = ierr + ierr2

  call kim_model_compute_arguments_get_argument_pointer( &
    model_compute_arguments_handle, &
    kim_compute_argument_name_particle_species_codes, &
    N, particlespeciesCodes, ierr2)
  ierr = ierr + ierr2
  call kim_model_compute_arguments_get_argument_pointer( &
    model_compute_arguments_handle, &
    kim_compute_argument_name_particle_contributing, N, particlecontributing, &
    ierr2)
  ierr = ierr + ierr2
  call kim_model_compute_arguments_get_argument_pointer( &
    model_compute_arguments_handle, &
    kim_compute_argument_name_coordinates, DIM, N, coor, ierr2)
  ierr = ierr + ierr2
  call kim_model_compute_arguments_get_argument_pointer( &
    model_compute_arguments_handle, &
    kim_compute_argument_name_partial_energy, energy, ierr2)
  ierr = ierr + ierr2
  call kim_model_compute_arguments_get_argument_pointer( &
    model_compute_arguments_handle, &
    kim_compute_argument_name_partial_forces, DIM, N, force, ierr2)
  ierr = ierr + ierr2
  call kim_model_compute_arguments_get_argument_pointer( &
    model_compute_arguments_handle, &
    kim_compute_argument_name_partial_particle_energy, N, particleEnergy, ierr2)
  ierr = ierr + ierr2
  if (ierr /= 0) then
    kim_log_message = "get_argument_pointer"
    LOG_ERROR()
    return
  endif

  ! Check to see if we have been asked to compute the energy, forces, energy per particle,
  ! dEdr, and d2Edr2
  !
  if (associated(energy)) then
    comp_energy =  1
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

  ierr = 0
  call kim_model_compute_arguments_is_callback_present( &
    model_compute_arguments_handle, &
    kim_compute_callback_name_process_dedr_term, comp_process_dedr, ierr2)
  ierr = ierr + ierr2
  call kim_model_compute_arguments_is_callback_present( &
    model_compute_arguments_handle, &
    kim_compute_callback_name_process_d2edr2_term, comp_process_d2edr2, ierr2)
  ierr = ierr + ierr2
  if (ierr /= 0) then
     kim_log_message = "get_compute"
     LOG_ERROR()
     return
  endif

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

  allocate( Rij(DIM) )
  if (comp_process_d2Edr2.eq.1) then
     allocate( r_pairs(2)       )
     allocate( Rij_pairs(DIM,2) )
     allocate( i_pairs(2)       )
     allocate( j_pairs(2)       )
  endif

  ! Initialize potential energies, forces
  !
  if (comp_particleEnergy.eq.1) particleEnergy = 0.0_cd
  if (comp_energy.eq.1) energy = 0.0_cd
  if (comp_force.eq.1)  force  = 0.0_cd

  !
  !  Compute energy and forces
  !

  !  Loop over particles and compute energy and forces
  !
  do i = 1,N
    if (particleContributing(i) == 1) then
      ! Set up neighbor list for next particle
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

        j = nei1part(jj)                           ! get neighbor ID

        if ( i.lt. j) then                         ! Effective half list
          ! compute relative position vector
          !
          Rij(:) = coor(:,j) - coor(:,i)           ! distance vector between i j

          ! compute energy and forces
          !
          Rsqij = dot_product(Rij,Rij)             ! compute square distance
          if ( Rsqij .lt. buf%cutsq(1) ) then      ! particles are interacting?

            r = sqrt(Rsqij)                        ! compute distance
            if (comp_process_d2Edr2.eq.1) then
              call calc_phi_dphi_d2phi(buf%eps(1),    &
                                       buf%sigma(1),  &
                                       buf%cutoff(1), &
                                       r,phi,dphi,d2phi) ! compute pair potential
                                                         !   and it derivatives
              if (particleContributing(j).eq.1) then
                dEidr  = dphi
                d2Eidr = d2phi
              else
                dEidr  = 0.5_cd*dphi
                d2Eidr = 0.5_cd*d2phi
              endif
            elseif (comp_force.eq.1.or.comp_process_dEdr.eq.1) then
              call calc_phi_dphi(buf%eps(1),    &
                                 buf%sigma(1),  &
                                 buf%cutoff(1), &
                                 r,phi,dphi)        ! compute pair potential
                                                    !   and it derivative
              if (particleContributing(j).eq.1) then
                dEidr = dphi
              else
                dEidr = 0.5_cd*dphi
              endif
            else
              call calc_phi(buf%eps(1),    &
                            buf%sigma(1),  &
                            buf%cutoff(1), &
                            r,phi)                  ! compute just pair potential
            endif

            ! contribution to energy
            !
            if (comp_particleEnergy.eq.1) then
              particleEnergy(i) = particleEnergy(i) + 0.5_cd*phi     ! accumulate energy
              if (particleContributing(j).eq.1) then
                particleEnergy(j) = particleEnergy(j) + 0.5_cd*phi   ! accumulate energy
              endif
            endif
            if (comp_energy.eq.1) then
              if (particleContributing(j).eq.1) then
                energy = energy + phi              !      add half v to total energy
              else
                energy = energy + 0.5_cd*phi       !      add half v to total energy
              endif
            endif

            ! contribution to process_dEdr
            !
            if (comp_process_dEdr.eq.1) then
              call kim_model_compute_arguments_process_dedr_term( &
                model_compute_arguments_handle, deidr, r, c_loc(Rij(1)), i, j, ierr)
            endif

            ! contribution to process_d2Edr2
            if (comp_process_d2Edr2.eq.1) then
              r_pairs(1) = r
              r_pairs(2) = r
              Rij_pairs(:,1) = Rij
              Rij_pairs(:,2) = Rij
              i_pairs(1) = i
              i_pairs(2) = i
              j_pairs(1) = j
              j_pairs(2) = j

              call kim_model_compute_arguments_process_d2edr2_term( &
                model_compute_arguments_handle, d2eidr, &
                c_loc(r_pairs(1)),     &
                c_loc(Rij_pairs(1,1)), &
                c_loc(i_pairs(1)),     &
                c_loc(j_pairs(1)), ierr)

            endif

            ! contribution to forces
            !
            if (comp_force.eq.1) then
              force(:,i) = force(:,i) + dEidr*Rij/r ! accumulate force on atom i
              force(:,j) = force(:,j) - dEidr*Rij/r ! accumulate force on atom j
            endif

          endif ! Check on whether particle jj is interacting with particle i
        endif ! if ( i.lt.j )
      enddo  ! loop on jj
    endif ! Check on whether particle i is contributing
  enddo  ! infinite do loop (terminated by exit statements above)

  ! Free temporary storage
  !
  deallocate( Rij )
  if (comp_process_d2Edr2.eq.1) then
    deallocate( r_pairs   )
    deallocate( Rij_pairs )
    deallocate( i_pairs   )
    deallocate( j_pairs   )
  endif

  ! No errors
  !
  ierr = 0
  return

end subroutine Compute_Energy_Forces
!-------------------------------------------------------------------------------
!
! Model driver refresh routine
!
!-------------------------------------------------------------------------------
subroutine refresh(model_refresh_handle, ierr) bind(c)
  implicit none

  !-- Transferred variables
  type(kim_model_refresh_handle_type), intent(inout) :: model_refresh_handle
  integer(c_int), intent(out) :: ierr

  !-- Local variables
  type(BUFFER_TYPE), pointer :: buf; type(c_ptr) :: pbuf

  ! Get model buffer from KIM object
  call kim_model_refresh_get_model_buffer_pointer(model_refresh_handle, pbuf)
  call c_f_pointer(pbuf, buf)

  call kim_model_refresh_set_influence_distance_pointer(model_refresh_handle, &
    buf%influence_distance)
  call kim_model_refresh_set_neighbor_list_pointers(model_refresh_handle, &
    1, buf%cutoff, buf%padding_neighbor_hints, buf%half_list_hints)

  ! Set new values in KIM object and buffer
  !
  buf%influence_distance = buf%cutoff(1)
  buf%cutsq(1) = buf%cutoff(1)*buf%cutoff(1)

  ierr = 0
  return

end subroutine refresh

!-------------------------------------------------------------------------------
!
! Model driver destroy routine
!
!-------------------------------------------------------------------------------
subroutine destroy(model_destroy_handle, ierr) bind(c)
  implicit none

  !-- Transferred variables
  type(kim_model_destroy_handle_type), intent(inout) :: model_destroy_handle
  integer(c_int), intent(out) :: ierr

  !-- Local variables
  type(BUFFER_TYPE), pointer :: buf; type(c_ptr) :: pbuf

  ! get model buffer from KIM object
  call kim_model_destroy_get_model_buffer_pointer(model_destroy_handle, pbuf)
  call c_f_pointer(pbuf, buf)

  deallocate( buf )

  ierr = 0
  return

end subroutine destroy

!-------------------------------------------------------------------------------
!
! Model driver compute arguments create routine
!
!-------------------------------------------------------------------------------
#include "kim_model_compute_arguments_create_log_macros.fd"
subroutine compute_arguments_create(model_compute_handle, &
  model_compute_arguments_create_handle, ierr) bind(c)
use kim_model_compute_arguments_create_module, &
    log_entry=>kim_model_compute_arguments_create_log_entry
  implicit none

  !-- Transferred variables
  type(kim_model_compute_handle_type), intent(in) :: model_compute_handle
  type(kim_model_compute_arguments_create_handle_type), intent(inout) :: &
    model_compute_arguments_create_handle
  integer(c_int), intent(out) :: ierr

  integer(c_int) ierr2

  ierr = 0
  ierr2 = 0

  ! register arguments
  call kim_model_compute_arguments_create_set_argument_support_status( &
    model_compute_arguments_create_handle, &
    kim_compute_argument_name_partial_energy, &
    kim_support_status_optional, ierr)
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
  if (ierr /= 0) then
    kim_log_message = "Unable to register arguments support_statuss"
    LOG_ERROR()
    goto 42
  end if

  ! register callbacks
  call kim_model_compute_arguments_create_set_callback_support_status( &
    model_compute_arguments_create_handle, &
    kim_compute_callback_name_process_dedr_term, &
    kim_support_status_optional, ierr)
  call kim_model_compute_arguments_create_set_callback_support_status( &
    model_compute_arguments_create_handle, &
    kim_compute_callback_name_process_d2edr2_term, &
    kim_support_status_optional, ierr2)
  ierr = ierr + ierr2
  if (ierr /= 0) then
    kim_log_message = "Unable to register callbacks support_statuss"
    LOG_ERROR()
    goto 42
  end if

  ierr = 0
  42 continue
  return

end subroutine compute_arguments_create

!-------------------------------------------------------------------------------
!
! Model driver compute arguments destroy routine
!
!-------------------------------------------------------------------------------
#include "kim_model_compute_arguments_destroy_log_macros.fd"
subroutine compute_arguments_destroy(model_compute_handle, &
  model_compute_arguments_destroy_handle, ierr) bind(c)
  use kim_model_compute_arguments_destroy_module, &
      log_entry=>kim_model_compute_arguments_destroy_log_entry
  implicit none

  !-- Transferred variables
  type(kim_model_compute_handle_type), intent(in) :: model_compute_handle
  type(kim_model_compute_arguments_destroy_handle_type), intent(inout) :: &
    model_compute_arguments_destroy_handle
  integer(c_int), intent(out) :: ierr

  ! nothing to be done

  ierr = 0

  return
end subroutine compute_arguments_destroy

end module Pair_Lennard_Jones_Truncated

!-------------------------------------------------------------------------------
!
! Model driver initialization routine (REQUIRED)
!
!-------------------------------------------------------------------------------
#include "kim_model_driver_create_log_macros.fd"
subroutine create(model_driver_create_handle, requested_length_unit, &
  requested_energy_unit, requested_charge_unit, requested_temperature_unit, &
  requested_time_unit,  ierr) bind(c)
  use, intrinsic :: iso_c_binding
  use Pair_Lennard_Jones_Truncated
  use kim_model_driver_headers_module
  implicit none

  !-- Transferred variables
  type(kim_model_driver_create_handle_type), intent(inout) :: model_driver_create_handle
  type(kim_length_unit_type), intent(in), value :: requested_length_unit
  type(kim_energy_unit_type), intent(in), value :: requested_energy_unit
  type(kim_charge_unit_type), intent(in), value :: requested_charge_unit
  type(kim_temperature_unit_type), intent(in), value :: requested_temperature_unit
  type(kim_time_unit_type), intent(in), value :: requested_time_unit
  integer(c_int), intent(out) :: ierr

  !-- Local variables
  integer(c_int), parameter :: cd = c_double ! used for literal constants
  integer(c_int), parameter :: one=1
  integer(c_int) :: number_of_parameter_files
  character(len=1024, kind=c_char) :: parameter_file_name
  integer(c_int) :: ierr2
  character(len=100, kind=c_char) :: in_species
  real(c_double) :: in_cutoff
  real(c_double) :: in_eps
  real(c_double) :: in_sigma
  type(kim_species_name_type) species_name
  type(BUFFER_TYPE), pointer :: buf

  kim_log_file = __FILE__

  ! Register units
  call kim_model_driver_create_set_units( &
    model_driver_create_handle, &
    requested_length_unit, &
    requested_energy_unit, &
    kim_charge_unit_unused, &
    kim_temperature_unit_unused, &
    kim_time_unit_unused, ierr)
  if (ierr /= 0) then
    kim_log_message = "Unable to set units"
    LOG_ERROR()
    goto 42
  end if

  ! Register numbering
  call kim_model_driver_create_set_model_numbering( &
    model_driver_create_handle, kim_numbering_one_based, ierr)
  if (ierr /= 0) then
    kim_log_message = "Unable to set numbering"
    LOG_ERROR()
    goto 42
  end if

  ! Store callback pointers in KIM object
  call kim_model_driver_create_set_compute_pointer( &
    model_driver_create_handle, kim_language_name_fortran, &
    c_funloc(Compute_Energy_Forces), ierr)
  call kim_model_driver_create_set_compute_arguments_create_pointer( &
    model_driver_create_handle, kim_language_name_fortran, &
    c_funloc(compute_arguments_create), ierr2)
  ierr = ierr + ierr2
  call kim_model_driver_create_set_compute_arguments_destroy_pointer( &
    model_driver_create_handle, kim_language_name_fortran, &
    c_funloc(compute_arguments_destroy), ierr2)
  ierr = ierr + ierr2
  call kim_model_driver_create_set_refresh_pointer( &
    model_driver_create_handle, kim_language_name_fortran, &
    c_funloc(refresh), ierr2)
  ierr = ierr + ierr2
  call kim_model_driver_create_set_destroy_pointer( &
    model_driver_create_handle, kim_language_name_fortran, &
    c_funloc(destroy), ierr2)
  ierr = ierr + ierr2
  if (ierr /= 0) then
    kim_log_message = "Unable to store callback pointers"
    LOG_ERROR()
    goto 42
  end if

  ! Process parameter files
  call kim_model_driver_create_get_number_of_parameter_files( &
    model_driver_create_handle, number_of_parameter_files)
  if (number_of_parameter_files .ne. 1) then
    kim_log_message = "Wrong number of parameter files"
    LOG_ERROR()
    ierr = 1
    goto 42
  end if

  ! Read in model parameters from parameter file
  call kim_model_driver_create_get_parameter_file_name( &
    model_driver_create_handle, 1, parameter_file_name, ierr)
  if (ierr /= 0) then
    kim_log_message = "Unable to get parameter file name"
    LOG_ERROR()
    ierr = 1
    goto 42
  end if

  allocate(buf)

  ! Read in model parameters from parameter file
  !
  open(10,file=parameter_file_name,status="old")
    read(10,*,iostat=ierr,err=100) in_species
    read(10,*,iostat=ierr,err=100) in_cutoff
    read(10,*,iostat=ierr,err=100) in_eps
    read(10,*,iostat=ierr,err=100) in_sigma
  close(10)

  goto 200
  100 continue
  ! Reading parameters failed
  ierr = 1
  kim_log_message = "Unable to read model parameters"
  LOG_ERROR()
  goto 42

  200 continue

  buf%influence_distance = in_cutoff
  buf%cutoff(1) = in_cutoff
  buf%cutsq(1) = in_cutoff*in_cutoff
  buf%eps(1) = in_eps
  buf%sigma(1) = in_sigma
  buf%padding_neighbor_hints(1) = 1
  buf%half_list_hints(1) = 1

  ! Store model cutoff in KIM object
  call kim_model_driver_create_set_influence_distance_pointer( &
    model_driver_create_handle, buf%influence_distance)
  call kim_model_driver_create_set_neighbor_list_pointers( &
    model_driver_create_handle, 1, buf%cutoff, &
    buf%padding_neighbor_hints, buf%half_list_hints)

  ! Register buffer in KIM API object
  call kim_model_driver_create_set_model_buffer_pointer( &
    model_driver_create_handle, c_loc(buf))

  ! Register species
  call kim_species_name_from_string(in_species, species_name)
  if (ierr /= 0) then
    kim_log_message = "Unable to set species_name"
    LOG_ERROR()
    goto 42
  end if

  call kim_model_driver_create_set_species_code( &
    model_driver_create_handle, species_name, speccode, ierr)
  if (ierr /= 0) then
    kim_log_message = "Unable to set species code"
    LOG_ERROR()
    goto 42
  end if

  ! Set parameter pointers so they can be changed by the simulator
  call kim_model_driver_create_set_parameter_pointer( &
    model_driver_create_handle, buf%cutoff, "cutoff", ierr)
  if (ierr /= 0) then
    kim_log_message = "set_parameter"
    LOG_ERROR()
    goto 42
  endif

  call kim_model_driver_create_set_parameter_pointer( &
    model_driver_create_handle, buf%eps, "epsilon", ierr)
  if (ierr /= 0) then
    kim_log_message = "set_parameter"
    LOG_ERROR()
    goto 42
  endif

  call kim_model_driver_create_set_parameter_pointer( &
    model_driver_create_handle, buf%sigma, "sigma", ierr)
  if (ierr /= 0) then
    kim_log_message = "set_parameter"
    LOG_ERROR()
    goto 42
  endif

  ierr = 0
  42 continue
  return

end subroutine create
