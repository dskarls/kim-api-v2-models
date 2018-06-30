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
!
! CDDL HEADER END
!

!
! Copyright (c) 2012, Mark R. Gilbert, CCFE Fusion Association.
! All rights reserved.
!
! Contributors:
!    Mark R. Gilbert
!    Ryan S. Elliott
!
! Converted to be compatible with KIM API 2.0 by Daniel S. Karls (UMinn)

!****************************************************************************
!**
!**  MODULE model_driver_PF_cubic_splines
!**
!**  EAM-like potential with cubic splines representing knot functions
!**  magnetic ability also available via B parameter
!**
!**  Language: Fortran 2003
!**
!**
!**
!****************************************************************************

#define THIS_FILE_NAME __FILE__
#define TRUEFALSE(TRUTH) merge(1,0,(TRUTH))

module model_driver_pf_cubic_splines

use, intrinsic :: iso_c_binding
use kim_model_driver_headers_module
implicit none

save
private
public Compute_Energy_Forces, &
       compute_arguments_create, &
       compute_arguments_destroy, &
       refresh, &
       destroy, &
       speccode, &
       buffer_type

type, bind(c) :: buffer_type
  real(c_double) :: influence_distance
  real(c_double) :: cutoff(1)  ! cutoff radius in angstroms
  real(c_double) :: cutsq(1)

  ! new variables for general potential - to be read in - 12/10/06
  real(c_double) :: A, B,Z,r1,r2,a_inter(7),a_rho
  integer(c_int) :: nknotv,nknotp,interpolate_num
  logical(c_bool) :: linear_interpolate
  type(c_ptr) :: vknotcoeff, vknotpoint, pknotcoeff, pknotpoint
end type buffer_type

! Below are the definitions and values of all Model parameters
integer(c_int), parameter          :: cd = c_double  ! used for literal constants
integer(c_int), parameter          :: DIM=3          ! dimensionality of space
integer(c_int), parameter          :: speccode = 1   ! internal species code

real(c_double) , parameter :: pc=1.0_cd, pi=3.141592653589793_cd
! electric constant (in F/m) and electron charge (in J)
real(c_double) , parameter :: E_0=8.854187817e-12_cd, Ec=1.60217653e-19_cd
! Length conversion factor (multiply by this to convert from Bohr to Angstroms)
real(c_double) , parameter :: a_B=0.5291772108_cd

! constants of biersack-ziegler coulomb potential
real(c_double) , parameter :: p1=0.1818_cd, p2=0.5099_cd, p3=0.2802_cd, p4=0.02817_cd
real(c_double) , parameter :: Bt1=-3.2_cd, Bt2=-0.9423_cd, Bt3=-0.4029_cd, Bt4=-0.2016_cd
!

contains

!-------------------------------------------------------------------------------
!
!  Calculate pair potential phi(r)
!
!-------------------------------------------------------------------------------
subroutine calc_phi(r,buf,phi)
  implicit none

  !-- Transferred variables
  real(c_double),    intent(in)  :: r
  type(buffer_type), intent(in)  :: buf
  real(c_double),    intent(out) :: phi

  !-- Local variables
  real(c_double) :: a_s,phi_r,B_Zpre
  integer(c_int) :: i
  real(c_double), pointer :: vknotcoeff(:), vknotpoint(:)

  ! Set local pointers to buffer spline pointers
  call c_f_pointer(buf%vknotcoeff, vknotcoeff, (/buf%nknotv/))
  call c_f_pointer(buf%vknotpoint, vknotpoint, (/buf%nknotv/))

  if (r .gt. buf%cutoff(1)) then
    ! Argument exceeds cutoff radius
    phi = 0.0_cd
  else
    IF (r>=buf%r2) THEN
      phi=0.0_cd
      DO i=1,buf%nknotv
        IF (r<vknotpoint(i)) THEN
          phi=phi+vknotcoeff(i)*(vknotpoint(i)-r)**3
        END IF
      END DO
    ELSE IF (r<=buf%r1) THEN
      a_s=0.88534_cd*a_B/(sqrt(2.0_cd)*(buf%Z**(1.0_cd/3.0_cd)))
      B_Zpre=buf%Z**2*Ec*1.0e10_cd/(4.0_cd*pi*E_0)
      phi_r=p1*exp(Bt1*r/a_s)+p2*exp(Bt2*r/a_s)+p3*exp(Bt3*r/a_s)+ &
            p4*exp(Bt4*r/a_s)
      phi=B_Zpre*phi_r/r
    ELSE
      IF (buf%linear_interpolate) THEN
        phi=buf%a_inter(1)+buf%a_inter(2)*r+buf%a_inter(3)*(r**2)+ &
        buf%a_inter(4)*(r**3)+buf%a_inter(5)*(r**4)+buf%a_inter(6)*(r**5)+buf%a_inter(7)*(r**6)
      ELSE
        phi=exp(buf%a_inter(1)+buf%a_inter(2)*r+buf%a_inter(3)*(r**2)+ &
        buf%a_inter(4)*(r**3)+buf%a_inter(5)*(r**4)+buf%a_inter(6)*(r**5)+buf%a_inter(7)*(r**6))
      END IF
    END IF
  endif

end subroutine calc_phi

!-------------------------------------------------------------------------------
!
!  Calculate pair potential phi(r) and its derivative dphi(r)
!
!-------------------------------------------------------------------------------
subroutine calc_phi_dphi(r,buf,phi,dphi)
  implicit none

  !-- Transferred variables
  real(c_double), intent(in)  :: r
  type(buffer_type), intent(in)  :: buf
  real(c_double), intent(out) :: phi,dphi

  !-- Local variables
  real(c_double) :: a_s,phi_r,B_Zpre,dphi_r
  integer(c_int) :: i
  real(c_double), pointer :: vknotcoeff(:), vknotpoint(:)

  ! Set local pointers to buffer spline pointers
  call c_f_pointer(buf%vknotcoeff, vknotcoeff, (/buf%nknotv/))
  call c_f_pointer(buf%vknotpoint, vknotpoint, (/buf%nknotv/))

  if (r .gt. buf%cutoff(1)) then
    ! Argument exceeds cutoff radius
    phi    = 0.0_cd
    dphi   = 0.0_cd
  else
    IF (r>=buf%r2) THEN
      dphi = 0.0_cd
      phi  = 0.0_cd
      DO i=1,buf%nknotv
        IF (r<vknotpoint(i)) THEN
          dphi=dphi-3*vknotcoeff(i)*(vknotpoint(i)-r)**2
          phi=phi+vknotcoeff(i)*(vknotpoint(i)-r)**3
        END IF
      END DO
    ELSE IF (r<=buf%r1) THEN
      a_s=0.88534_cd*a_B/(sqrt(2.0_cd)*(buf%Z**(1.0_cd/3.0_cd)))
      B_Zpre=buf%Z**2*Ec*1.0e10_cd/(4.0_cd*pi*E_0)
      phi_r=p1*exp(Bt1*r/a_s)+p2*exp(Bt2*r/a_s)+p3*exp(Bt3*r/a_s)+ &
            p4*exp(Bt4*r/a_s)
      dphi_r=(p1*Bt1*exp(Bt1*r/a_s)+p2*Bt2*exp(Bt2*r/a_s)+p3*Bt3*exp(Bt3*r/a_s)+ &
              p4*Bt4*exp(Bt4*r/a_s))/a_s
      dphi=-B_Zpre*phi_r/(r**2)+B_Zpre*dphi_r/r
      phi=B_Zpre*phi_r/r
    ELSE
      IF (buf%linear_interpolate) THEN
        dphi=buf%a_inter(2)+2.0_cd*buf%a_inter(3)*r+3.0_cd*buf%a_inter(4)*(r**2)+ &
             4.0_cd*buf%a_inter(5)*(r**3)+5.0_cd*buf%a_inter(6)*(r**4)+6.0_cd*buf%a_inter(7)*(r**5)
        phi=buf%a_inter(1)+buf%a_inter(2)*r+buf%a_inter(3)*(r**2)+ &
            buf%a_inter(4)*(r**3)+buf%a_inter(5)*(r**4)+buf%a_inter(6)*(r**5)+buf%a_inter(7)*(r**6)
      ELSE
        phi_r=exp(buf%a_inter(1)+buf%a_inter(2)*r+buf%a_inter(3)*(r**2)+ &
              buf%a_inter(4)*(r**3)+buf%a_inter(5)*(r**4)+buf%a_inter(6)*(r**5)+buf%a_inter(7)*(r**6))
        dphi_r=(buf%a_inter(2)+2.0_cd*buf%a_inter(3)*r+3.0_cd*buf%a_inter(4)*(r**2)+ &
                4.0_cd*buf%a_inter(5)*(r**3)+5.0_cd*buf%a_inter(6)*(r**4)+6.0_cd*buf%a_inter(7)*(r**5))
        dphi=phi_r*dphi_r
        phi=exp(buf%a_inter(1)+buf%a_inter(2)*r+buf%a_inter(3)*(r**2)+ &
            buf%a_inter(4)*(r**3)+buf%a_inter(5)*(r**4)+buf%a_inter(6)*(r**5)+buf%a_inter(7)*(r**6))
      END IF
    END IF
  endif

end subroutine calc_phi_dphi

!-------------------------------------------------------------------------------
!
!  Calculate electron density g(r)
!
!-------------------------------------------------------------------------------
subroutine calc_g(r,buf,g)
  implicit none

  !-- Transferred variables
  real(c_double), intent(in)  :: r
  type(buffer_type), intent(in)  :: buf
  real(c_double), intent(out) :: g

  !-- Local variables
  integer(c_int) :: i
  real(c_double), pointer :: pknotcoeff(:), pknotpoint(:)

  ! Set local pointers to buffer spline pointers
  call c_f_pointer(buf%pknotcoeff, pknotcoeff, (/buf%nknotp/))
  call c_f_pointer(buf%pknotpoint, pknotpoint, (/buf%nknotp/))

  if (r .gt. buf%cutoff(1)) then
    ! Argument exceeds cutoff radius
    g = 0.0_cd
  else
    g = 0.0_cd
    DO i=1,buf%nknotp
      IF (r<pknotpoint(i)) THEN
        g=g+pknotcoeff(i)*(pknotpoint(i)-r)**3
      END IF
    END DO
  endif

end subroutine calc_g

!-------------------------------------------------------------------------------
!
!  Calculate electron density derivative dg(r)
!
!-------------------------------------------------------------------------------
subroutine calc_dg(r,buf,dg)
  implicit none

  !-- Transferred variables
  real(c_double), intent(in)  :: r
  type(buffer_type), intent(in)  :: buf
  real(c_double), intent(out) :: dg

  !-- Local variables
  integer(c_int) :: i
  real(c_double), pointer :: pknotcoeff(:), pknotpoint(:)

  ! Set local pointers to buffer spline pointers
  call c_f_pointer(buf%pknotcoeff, pknotcoeff, (/buf%nknotp/))
  call c_f_pointer(buf%pknotpoint, pknotpoint, (/buf%nknotp/))

  if (r .gt. buf%cutoff(1)) then
    ! Argument exceeds cutoff radius
    dg = 0.0_cd
  else
    dg = 0.0_cd
    DO i=1,buf%nknotp
      IF (r<pknotpoint(i)) THEN
        dg=dg-3.0_cd*pknotcoeff(i)*(pknotpoint(i)-r)**2
      END IF
    END DO
  endif

end subroutine calc_dg

!-------------------------------------------------------------------------------
!
!  Calculate embedding function U(rho)
!
!-------------------------------------------------------------------------------
subroutine calc_U(rho,buf,U)
  implicit none

  !-- Transferred variables
  real(c_double),         intent(in)  :: rho
  type(buffer_type), intent(in)  :: buf
  real(c_double),         intent(out) :: U

  IF (rho<pc) THEN
    U=-buf%A*sqrt(rho)-(buf%B/log(2.0_cd))*(1.0_cd-sqrt(rho/pc))*log(2.0_cd-rho/pc)+ &
       buf%a_rho*(rho**2)
  ELSE
    U=-buf%A*sqrt(rho)+buf%a_rho*(rho**2)
  END IF

end subroutine calc_U

!-------------------------------------------------------------------------------
!
!  Calculate embedding function U(rho) and first derivative dU(rho)
!
!-------------------------------------------------------------------------------
subroutine calc_U_dU(rho,buf,U,dU)
  implicit none

  !-- Transferred variables
  real(c_double),         intent(in)  :: rho
  type(buffer_type), intent(in)  :: buf
  real(c_double),         intent(out) :: U,dU

  IF(rho.le.1e-10_cd) THEN
    dU = 0.0_cd
    U  = 0.0_cd
  ELSE
    IF (rho<pc) THEN
      U=-buf%A*sqrt(rho)-(buf%B/log(2.0_cd))*(1.0_cd-sqrt(rho/pc))*log(2.0_cd-rho/pc)+ &
        buf%a_rho*(rho**2)
      dU=-buf%A/(2.0_cd*sqrt(rho))-(buf%B/log(2.0_cd))*( &
        (-1.0_cd/(2.0_cd*sqrt(rho*pc)))*log(2.0_cd-rho/pc)+ &
        (1.0_cd-sqrt(rho/pc))*(-1.0_cd/(2.0_cd*pc-rho)))+2.0_cd*buf%a_rho*rho
    ELSE
      U=-buf%A*sqrt(rho)+buf%a_rho*(rho**2)
      dU=-buf%A/(2.0_cd*sqrt(rho))+2.0_cd*buf%a_rho*rho
    END IF
  END IF

end subroutine calc_U_dU


!-------------------------------------------------------------------------------
!
!  Calculate embedding function U(rho) and first +2nd derivative dU(rho)
!
!-------------------------------------------------------------------------------
subroutine calc_U_dU_ddU(rho,buf,U,dU,ddU)
  implicit none

  !-- Transferred variables
  real(c_double),         intent(in)  :: rho
  type(buffer_type), intent(in)  :: buf
  real(c_double),         intent(out) :: U,dU,ddU

  !-- Local variables
  ! <FILL place any local variable definitions here>
  IF(rho.le.1e-10_cd) THEN
    dU  = 0.0_cd
    U   = 0.0_cd
    ddU = 0.0_cd
  ELSE
    IF (rho<pc) THEN
      U=-buf%A*sqrt(rho)-(buf%B/log(2.0_cd))*(1.0_cd-sqrt(rho/pc))*log(2.0_cd-rho/pc)+ &
        buf%a_rho*(rho**2)
      dU=-buf%A/(2.0_cd*sqrt(rho))-(buf%B/log(2.0_cd))*( &
         (-1.0_cd/(2.0_cd*sqrt(rho*pc)))*log(2.0_cd-rho/pc)+ &
         (1.0_cd-sqrt(rho/pc))*(-1.0_cd/(2.0_cd*pc-rho)))+2.0_cd*buf%a_rho*rho
      ddU=buf%A/(4.0_cd*(rho**(3.0_cd/2.0_cd)))-(buf%B/log(2.0_cd))*( &
          (rho**(-3.0_cd/2.0_cd)/(4.0_cd*sqrt(pc)))*log(2.0_cd-rho/pc)+ &
          (1.0_cd/sqrt(pc*rho))*(1.0_cd/(2.0_cd*pc-rho))- &
          (1.0_cd-sqrt(rho/pc))*(1.0_cd/((2.0_cd*pc-rho)**2)))+2.0_cd*buf%a_rho
    ELSE
      U=-buf%A*sqrt(rho)+buf%a_rho*(rho**2)
      dU=-buf%A/(2.0_cd*sqrt(rho))+2.0_cd*buf%a_rho*rho
      ddU=buf%A/(4.0_cd*(rho**(3.0_cd/2.0_cd)))+2.0_cd*buf%a_rho
    END IF
  END IF

end subroutine calc_U_dU_ddU

!-------------------------------------------------------------------------------
!
! Compute energy and forces on atoms from the positions.
!
!-------------------------------------------------------------------------------
#include "kim_model_compute_log_macros.fd"
subroutine Compute_Energy_Forces(model_compute_handle, &
  compute_arguments_handle, ierr) bind(c)
  implicit none

  !-- Transferred variables
  type(kim_model_compute_handle_type), intent(in) :: model_compute_handle
  type(kim_model_compute_arguments_handle_type), intent(in) :: &
    compute_arguments_handle
  integer(c_int), intent(out) :: ierr

  !-- Local variables
  real(c_double) :: Rij(DIM)
  real(c_double) :: r,Rsqij,phi,dphi,g,dg,dU,U,dphieff = 0.0_cd
  real(c_double) :: dphii,dUi,Ei,dphij,dUj,Ej
  integer(c_int) :: i,j,jj,numnei,comp_energy,comp_force,comp_particleEnergy,comp_virial
  integer(c_int) :: ierr2
  real(c_double), allocatable :: rho(:),derU(:)

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
  type(buffer_type), pointer :: buf; type(c_ptr) :: pbuf

  kim_log_file = __FILE__

  ! Get buffer address from KIM API
  call kim_model_compute_get_model_buffer_pointer(model_compute_handle, pbuf)
  call c_f_pointer(pbuf, buf)

  ! Unpack data from KIM object
  !
  ierr = 0
  call kim_model_compute_arguments_get_argument_pointer( &
    compute_arguments_handle, &
    kim_compute_argument_name_number_of_particles, N, ierr2)
  ierr = ierr + ierr2
  call kim_model_compute_arguments_get_argument_pointer( &
    compute_arguments_handle, &
    kim_compute_argument_name_particle_species_codes, n, particleSpeciesCodes, &
    ierr2)
  ierr = ierr + ierr2
  call kim_model_compute_arguments_get_argument_pointer( &
    compute_arguments_handle, &
    kim_compute_argument_name_particle_contributing, n, particleContributing, &
    ierr2)
  ierr = ierr + ierr2
  call kim_model_compute_arguments_get_argument_pointer( &
    compute_arguments_handle, &
    kim_compute_argument_name_coordinates, dim, n, coor, ierr2)
  ierr = ierr + ierr2
  call kim_model_compute_arguments_get_argument_pointer( &
    compute_arguments_handle, &
    kim_compute_argument_name_partial_energy, energy, ierr2)
  ierr = ierr + ierr2
  call kim_model_compute_arguments_get_argument_pointer( &
    compute_arguments_handle, &
    kim_compute_argument_name_partial_forces, dim, n, force, ierr2)
  ierr = ierr + ierr2
  call kim_model_compute_arguments_get_argument_pointer( &
    compute_arguments_handle, &
    kim_compute_argument_name_partial_particle_energy, n, particleEnergy, ierr2)
  ierr = ierr + ierr2
  call kim_model_compute_arguments_get_argument_pointer( &
    compute_arguments_handle, &
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

  ! Initialize potential energies, forces, virial term, electron density
  !
  ! Note: the variable `particleEnergy' does not need to be initialized
  !       because it's initial value is set during the embedding energy
  !       calculation.
  !
  if (comp_energy.eq.1) energy = 0.0_cd
  if (comp_force.eq.1)  force(1:DIM,1:N) = 0.0_cd
  if (comp_virial.eq.1) virial = 0.0_cd
  allocate( rho(N) )  ! pair functional electron density
  rho(1:N) = 0.0_cd
  if (comp_force.eq.1.or.comp_virial.eq.1) allocate( derU(N) )  ! EAM embedded energy deriv

  !
  !  Compute energy and forces
  !

  !  Loop over particles in the neighbor list a first time,
  !  to compute electron density (=coordination)
  !
  do i = 1,N
    if(particleContributing(i).eq.1) then
      ! Get neighbor list of current atom
      call kim_model_compute_arguments_get_neighbor_list( &
        compute_arguments_handle, 1, i, numnei, nei1part, ierr)
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

        ! compute contribution to electron density
        !
        Rsqij = dot_product(Rij,Rij)             ! compute square distance
        if ( Rsqij .lt. buf%cutsq(1) ) then      ! particles are interacting?
          r = sqrt(Rsqij)                        ! compute distance
          call calc_g(r,buf,g)                   ! compute electron density
          rho(i) = rho(i) + g                    ! accumulate electron density
        endif

      enddo  ! loop on jj
    endif ! Check on whether particle is contributing
  enddo  ! infinite do loop (terminated by exit statements above)

  !  Now that we know the electron densities, calculate embedding part of energy
  !  U and its derivative U' (derU)
  !
  do i = 1,N
    if(particleContributing(i).eq.1) then
      if (comp_force.eq.1.or.comp_virial.eq.1) then
        call calc_U_dU(rho(i),buf,U,dU)          ! compute embedding energy
                                                 !   and its derivative
        derU(i) = dU                             ! store dU for later use
      else
        call calc_U(rho(i),buf,U)                ! compute just embedding energy
      endif

      ! accumulate the embedding energy contribution
      !
      ! Assuming U(rho=0) = 0.0_cd
      !
      if (comp_particleEnergy.eq.1) then         ! accumulate embedding energy contribution
        particleEnergy(i) = U
      endif
      if (comp_energy.eq.1) then
        energy = energy + U
      endif

    endif ! Check on whether particle is contributing
  enddo

  !  Loop over particles in the neighbor list a second time, to compute
  !  the forces and complete energy calculation
  !

  do i = 1,N
    if(particleContributing(i).eq.1) then
      ! Get neighbor list of current atom
      call kim_model_compute_arguments_get_neighbor_list( &
        compute_arguments_handle, 1, i, numnei, nei1part, ierr)
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

        j = nei1part(jj)                       ! get neighbor ID

        ! compute relative position vector
        !
        Rij(:) = coor(:,j) - coor(:,i)         ! distance vector between i j

        ! compute energy and forces
        !
        Rsqij = dot_product(Rij,Rij)           ! compute square distance
        if ( Rsqij .lt. buf%cutsq(1) ) then    ! particles are interacting?

          r = sqrt(Rsqij)                      ! compute distance
          if (comp_force.eq.1.or.comp_virial.eq.1) then
            call calc_phi_dphi(r,buf,phi,dphi) ! compute pair potential
                                               !   and it derivative
            call calc_dg(r,buf,dg)             ! compute elect dens first deriv
            dphii  = 0.5_cd*dphi
            dphij  = 0.0_cd
            dUi    = derU(i)*dg
            dUj    = 0.0_cd
            dphieff = dphii + dphij + dUi + dUj
          else
            call calc_phi(r,buf,phi)           ! compute just pair potential
          endif
          Ei = 0.5_cd*phi
          Ej = 0.0_cd

          ! contribution to energy
          !
          if (comp_particleEnergy.eq.1) then
            particleEnergy(i) = particleEnergy(i) + Ei          ! accumulate energy Ei
            particleEnergy(j) = particleEnergy(j) + Ej          ! accumulate energy Ej
          endif
          if (comp_energy.eq.1) then
            energy = energy + Ei                  ! accumulate energy
            energy = energy + Ej                  ! accumulate energy
          endif

          ! contribution to virial tensor
          !
          if (comp_virial.eq.1) then
            virial(1) = virial(1) + Rij(1)*Rij(1)*dphieff/r
            virial(2) = virial(2) + Rij(2)*Rij(2)*dphieff/r
            virial(3) = virial(3) + Rij(3)*Rij(3)*dphieff/r
            virial(4) = virial(4) + Rij(2)*Rij(3)*dphieff/r
            virial(5) = virial(5) + Rij(1)*Rij(3)*dphieff/r
            virial(6) = virial(6) + Rij(1)*Rij(2)*dphieff/r
          endif

          ! contribution to forces
          !
          if (comp_force.eq.1) then                        ! Ei contribution
            force(:,i) = force(:,i) + dphieff*Rij/r ! accumulate force on atom i
            force(:,j) = force(:,j) - dphieff*Rij/r ! accumulate force on atom j
          endif

        endif

      enddo  ! loop on jj
    endif ! Check on whether particle is contributing
  enddo  ! infinite do loop (terminated by exit statements above)

  ! Free temporary storage
  !
  deallocate( rho )
  if (comp_force.eq.1.or.comp_virial.eq.1) deallocate( derU )

  ! Everything is great
  !
  ierr = 0

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
  type(buffer_type), pointer :: buf; type(c_ptr) :: pbuf

  ! get model buffer from KIM object
  call kim_model_refresh_get_model_buffer_pointer(model_refresh_handle, pbuf)
  call c_f_pointer(pbuf, buf)

  call kim_model_refresh_set_influence_distance_pointer(model_refresh_handle, &
    buf%influence_distance)
  call kim_model_refresh_set_neighbor_list_cutoffs_pointer(model_refresh_handle, &
    1, buf%cutoff)

  ierr = 0
  return
end subroutine refresh

!-------------------------------------------------------------------------------
!
! Model driver destroy routine
!
!-------------------------------------------------------------------------------
subroutine destroy(destroy_handle, ierr) bind(c)
  implicit none

  !-- Transferred variables
  type(kim_model_destroy_handle_type), intent(inout) :: destroy_handle
  integer(c_int), intent(out) :: ierr

  !-- Local variables
  type(buffer_type), pointer :: buf; type(c_ptr) :: pbuf
  real(c_double), pointer :: pknotcoeff(:), pknotpoint(:), vknotcoeff(:), &
    vknotpoint(:)

  ! get model buffer from KIM object
  call kim_model_destroy_get_model_buffer_pointer(destroy_handle, pbuf)
  call c_f_pointer(pbuf, buf)

  ! Set local pointers to buffer spline pointers
  call c_f_pointer(buf%pknotcoeff, pknotcoeff, (/buf%nknotp/))
  call c_f_pointer(buf%pknotpoint, pknotpoint, (/buf%nknotp/))
  call c_f_pointer(buf%vknotcoeff, vknotcoeff, (/buf%nknotv/))
  call c_f_pointer(buf%vknotpoint, vknotpoint, (/buf%nknotv/))

  deallocate( pknotcoeff, pknotpoint )
  deallocate( vknotcoeff, vknotpoint )

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

  ! Register arguments
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
  call kim_model_compute_arguments_create_set_argument_support_status( &
    model_compute_arguments_create_handle, &
    kim_compute_argument_name_partial_virial, &
    kim_support_status_optional, ierr2)
  ierr = ierr + ierr2

  if (ierr /= 0) then
    kim_log_message = "Unable to register arguments support_statuss"
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
  compute_arguments_destroy_handle, ierr) bind(c)
  use kim_model_compute_arguments_destroy_module, &
    log_entry=>kim_model_compute_arguments_destroy_log_entry
  implicit none

  !-- Transferred variables
  type(kim_model_compute_handle_type), intent(in) :: model_compute_handle
  type(kim_model_compute_arguments_destroy_handle_type), intent(inout) :: &
    compute_arguments_destroy_handle
  integer(c_int), intent(out) :: ierr

  ! nothing to be done

  ierr = 0

  return
end subroutine compute_arguments_destroy

end module model_driver_pf_cubic_splines

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
  use model_driver_pf_cubic_splines
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
  integer(c_int), parameter :: cd = c_double  ! used for literal constants
  integer(c_int) :: number_of_parameter_files
  character(len=1024, kind=c_char) :: parameter_file_name
  integer(c_int) :: ierr2
  integer(c_int), parameter :: one=1
  integer(c_int) i, idum
  character (LEN=55) :: text
  character(len=100, kind=c_char) :: in_species
  type(kim_species_name_type) species_name
  logical ldum
  type(buffer_type), pointer :: buf
  real(c_double), pointer :: pknotcoeff(:), pknotpoint(:), vknotcoeff(:), &
    vknotpoint(:)

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

  open(10,file=parameter_file_name,status="old")
    read(10,*,iostat=ierr,err=100) in_species
    read(10,*,iostat=ierr,err=100) buf%cutoff(1)
    read(10,*,iostat=ierr,err=100) buf%Z
    read(10,*,iostat=ierr,err=100) buf%r1, buf%r2, &
      buf%linear_interpolate, buf%interpolate_num, buf%a_rho
    read(10,*,iostat=ierr,err=100) buf%nknotp, buf%nknotv
    !
    ! for A,B, and knot functions there are numbers and logicals to skip
    read(10,*,iostat=ierr,err=100) idum, buf%A, ldum
    !
    !
    read(10,*,iostat=ierr,err=100) idum, buf%B, ldum

    ! Allocate local arrays for spline knots/coeffs
    allocate(pknotcoeff(buf%nknotp),pknotpoint(buf%nknotp))
    allocate(vknotcoeff(buf%nknotv),vknotpoint(buf%nknotv))

    ! Assign c_ptrs in buffer to c-addresses of local pointers
    buf%pknotcoeff = c_loc(pknotcoeff)
    buf%pknotpoint = c_loc(pknotpoint)
    buf%vknotcoeff = c_loc(vknotcoeff)
    buf%vknotpoint = c_loc(vknotpoint)

    i=1
    DO WHILE ((ierr==0).AND.(i.LE.buf%nknotp))
      read(10,*,iostat=ierr,err=100) idum, pknotcoeff(i), ldum, &
        pknotpoint(i)
      i=i+1
    END DO
    i=1
    DO WHILE ((ierr==0).AND.(i.LE.buf%nknotv))
      read(10,*,iostat=ierr,err=100) idum, vknotcoeff(i), ldum, &
        vknotpoint(i)
      i=i+1
    END DO
    i=1
    buf%a_inter=0.0_cd
    DO WHILE ((ierr==0).AND.(i.LE.buf%interpolate_num))
      read(10,*,iostat=ierr,err=100) buf%a_inter(i)
      i=i+1
    END DO
    !PRINT *,Z,r1,r2,A,B,a_rho,nknotp,nknotv,a_inter,model_cutoff,interpolate_num
    !PRINT *,vknotpoint,vknotcoeff,pknotpoint,pknotcoeff
    read(10,'(A55)',iostat=ierr,err=100) text
  close(10)

  write (*,*) 'Potential info: ',TRIM(ADJUSTL(text))
  write (*,*) 'Atomic number: ', nint(buf%Z)

  goto 200
  100 continue
  ! Reading parameters failed
  ierr = 1
  kim_log_message = "Unable to read model parameters"
  LOG_ERROR()
  goto 42
  !PRINT *,Z,r1,r2,A,B,a_rho,nknotp,nknotv,a_inter,model_cutoff,interpolate_num
  !PRINT *,vknotpoint,vknotcoeff,pknotpoint,pknotcoeff

  200 continue

  ! Also set influence distance and cutsq
  buf%influence_distance = buf%cutoff(1)
  buf%cutsq(1) = buf%cutoff(1)*buf%cutoff(1)

  ! Store model cutoff in KIM object
  call kim_model_driver_create_set_influence_distance_pointer( &
    model_driver_create_handle, buf%influence_distance)
  call kim_model_driver_create_set_neighbor_list_cutoffs_pointer( &
    model_driver_create_handle, 1, buf%cutoff)

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

  ierr = 0
  42 continue
  return

end subroutine create
