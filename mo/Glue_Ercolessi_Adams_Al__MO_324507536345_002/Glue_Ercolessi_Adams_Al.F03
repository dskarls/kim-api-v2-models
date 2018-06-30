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
!    Ryan S. Elliott
!    Ellad B. Tadmor
!    Valeriu Smirichinski
!    Stephen M. Whalen
!

!****************************************************************************
!**
!**  MODULE Glue_Ercolessi_Adams_Al
!**
!**  Ercolessi-Adams pair functional model for Al
!**
!**  Reference: F. Ercolessi and J. B. Adams, Europhys. Lett. 26, 583 (1994).
!**             http://www.sissa.it/furio/potentials/Al/
!**
!**  Language: Fortran 2003
!**
!**
!****************************************************************************


#define THIS_FILE_NAME __FILE__
#define TRUEFALSE(TRUTH) merge(1,0,(TRUTH))

module Glue_Ercolessi_Adams_Al

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
integer(c_int), parameter :: cd = c_double  ! used for literal constants
integer(c_int), parameter :: DIM = 3        ! dimensionality of space
integer(c_int), parameter :: speccode = 1   ! internal species code
real(c_double), parameter :: model_cutoff  = 5.55805441821810_cd ! cutoff radius
                                                                 ! in angstroms
real(c_double), parameter :: model_cutsq   = model_cutoff**2

!-------------------------------------------------------------------------------
! Below are the definitions and values of all additional model parameters
!
! Recall that the Fortran 2003 format for declaring parameters is as follows:
!
integer(c_int), parameter :: npt= 17
integer(c_int), parameter :: nuu= 13
real(c_double), parameter :: rhomin = 0.0_cd

! pair potential data
!
real(c_double), parameter :: x(npt) =   (/    .202111069753385E+01_cd, &
                                              .227374953472558E+01_cd, &
                                              .252638837191732E+01_cd, &
                                              .277902720910905E+01_cd, &
                                              .303166604630078E+01_cd, &
                                              .328430488349251E+01_cd, &
                                              .353694372068424E+01_cd, &
                                              .378958255787597E+01_cd, &
                                              .404222139506771E+01_cd, &
                                              .429486023225944E+01_cd, &
                                              .454749906945117E+01_cd, &
                                              .480013790664290E+01_cd, &
                                              .505277674383463E+01_cd, &
                                              .530541558102636E+01_cd, &
                                              .555805441821810E+01_cd, &
                                              .555807968210182E+01_cd, &
                                              .555810494598553E+01_cd /)

real(c_double), parameter :: yv2(npt) = (/    .196016472197158E+01_cd, &
                                              .682724240745344E+00_cd, &
                                              .147370824539188E+00_cd, &
                                             -.188188235860390E-01_cd, &
                                             -.576011902692490E-01_cd, &
                                             -.519846499644276E-01_cd, &
                                             -.376352484845919E-01_cd, &
                                             -.373737879689433E-01_cd, &
                                             -.531351030124350E-01_cd, &
                                             -.632864983555742E-01_cd, &
                                             -.548103623840369E-01_cd, &
                                             -.372889232343935E-01_cd, &
                                             -.188876517630154E-01_cd, &
                                             -.585239362533525E-02_cd, &
                                              .000000000000000E+00_cd, &
                                              .000000000000000E+00_cd, &
                                              .000000000000000E+00_cd /)

real(c_double), parameter :: Bv2(npt) = (/   -.702739315585347E+01_cd, &
                                             -.333140549270729E+01_cd, &
                                             -.117329394261502E+01_cd, &
                                             -.306003283486901E+00_cd, &
                                             -.366656699104026E-01_cd, &
                                              .588330899204400E-01_cd, &
                                              .384220572312032E-01_cd, &
                                             -.390223173707191E-01_cd, &
                                             -.663882722510521E-01_cd, &
                                             -.312918894386669E-02_cd, &
                                              .590118945294245E-01_cd, &
                                              .757939459148246E-01_cd, &
                                              .643822548468606E-01_cd, &
                                              .399750987463792E-01_cd, &
                                              .177103852679117E-05_cd, &
                                             -.590423369301474E-06_cd, &
                                              .590654950414731E-06_cd /)

real(c_double), parameter :: Cv2(npt) = (/    .877545959718548E+01_cd, &
                                              .585407125495837E+01_cd, &
                                              .268820820643116E+01_cd, &
                                              .744718689404422E+00_cd, &
                                              .321378734769888E+00_cd, &
                                              .566263292669091E-01_cd, &
                                             -.137417679148505E+00_cd, &
                                             -.169124163201523E+00_cd, &
                                              .608037039066423E-01_cd, &
                                              .189589640245655E+00_cd, &
                                              .563784150384640E-01_cd, &
                                              .100486298765028E-01_cd, &
                                             -.552186092621482E-01_cd, &
                                             -.413902746758285E-01_cd, &
                                             -.116832934994489E+00_cd, &
                                              .233610871054729E-01_cd, &
                                              .233885865725971E-01_cd /)

real(c_double), parameter :: Dv2(npt) = (/   -.385449887634130E+01_cd, &
                                             -.417706040200591E+01_cd, &
                                             -.256425277368288E+01_cd, &
                                             -.558557503589276E+00_cd, &
                                             -.349316054551627E+00_cd, &
                                             -.256022933201611E+00_cd, &
                                             -.418337423301704E-01_cd, &
                                              .303368330939646E+00_cd, &
                                              .169921006301015E+00_cd, &
                                             -.175759761362548E+00_cd, &
                                             -.611278214082881E-01_cd, &
                                             -.861140219824535E-01_cd, &
                                              .182451950513387E-01_cd, &
                                             -.995395392057973E-01_cd, &
                                              .184972909229936E+04_cd, &
                                              .362829766922787E+00_cd, &
                                              .362829766922787E+00_cd /)

! electron density  data
!
real(c_double), parameter :: yrh(npt) = (/    .865674623712589E-01_cd, &
                                              .925214702944478E-01_cd, &
                                              .862003123832002E-01_cd, &
                                              .762736292751052E-01_cd, &
                                              .606481841271735E-01_cd, &
                                              .466030959588197E-01_cd, &
                                              .338740138848363E-01_cd, &
                                              .232572661705343E-01_cd, &
                                              .109046405489829E-01_cd, &
                                              .524910605677597E-02_cd, &
                                              .391702419142291E-02_cd, &
                                              .308277776293383E-02_cd, &
                                              .250214745349505E-02_cd, &
                                              .147220513798186E-02_cd, &
                                              .000000000000000E+00_cd, &
                                              .000000000000000E+00_cd, &
                                              .000000000000000E+00_cd /)

real(c_double), parameter :: Brh(npt) = (/    .608555214104682E-01_cd, &
                                             -.800158928716306E-02_cd, &
                                             -.332089451111092E-01_cd, &
                                             -.521001991705069E-01_cd, &
                                             -.618130637429111E-01_cd, &
                                             -.529750064268036E-01_cd, &
                                             -.442210477548108E-01_cd, &
                                             -.473645664984640E-01_cd, &
                                             -.390741582571631E-01_cd, &
                                             -.101795580610560E-01_cd, &
                                             -.318316981110289E-02_cd, &
                                             -.281217210746153E-02_cd, &
                                             -.236932031483360E-02_cd, &
                                             -.683554708271547E-02_cd, &
                                             -.638718204858808E-06_cd, &
                                              .212925486831149E-06_cd, &
                                             -.212983742465787E-06_cd /)

real(c_double), parameter :: Crh(npt) = (/   -.170233687052940E+00_cd, &
                                             -.102317878901959E+00_cd, &
                                              .254162872544396E-02_cd, &
                                             -.773173610292656E-01_cd, &
                                              .388717099948882E-01_cd, &
                                             -.388873819867093E-02_cd, &
                                              .385388290924526E-01_cd, &
                                             -.509815666327127E-01_cd, &
                                              .837968231208082E-01_cd, &
                                              .305743500420042E-01_cd, &
                                             -.288110886134041E-02_cd, &
                                              .434959924771674E-02_cd, &
                                             -.259669459714693E-02_cd, &
                                             -.150816117849093E-01_cd, &
                                              .421356801161513E-01_cd, &
                                             -.842575249165724E-02_cd, &
                                             -.843267014952237E-02_cd /)

real(c_double), parameter :: Drh(npt) = (/    .896085612514625E-01_cd, &
                                              .138352319847830E+00_cd, &
                                             -.105366473134009E+00_cd, &
                                              .153300619856764E+00_cd, &
                                             -.564184148788224E-01_cd, &
                                              .559792096400504E-01_cd, &
                                             -.118113795329664E+00_cd, &
                                              .177827488509794E+00_cd, &
                                             -.702220789044304E-01_cd, &
                                             -.441413511810337E-01_cd, &
                                              .954024354744484E-02_cd, &
                                             -.916498550800407E-02_cd, &
                                             -.164726813535368E-01_cd, &
                                              .754928689733184E-01_cd, &
                                             -.667110847110954E+03_cd, &
                                             -.912720300911022E-01_cd, &
                                             -.912720300911022E-01_cd /)

! Embedding function data
!
real(c_double), parameter :: xuu(nuu) = (/    .000000000000000E+00_cd, &
                                              .100000000000000E+00_cd, &
                                              .200000000000000E+00_cd, &
                                              .300000000000000E+00_cd, &
                                              .400000000000000E+00_cd, &
                                              .500000000000000E+00_cd, &
                                              .600000000000000E+00_cd, &
                                              .700000000000000E+00_cd, &
                                              .800000000000000E+00_cd, &
                                              .900000000000000E+00_cd, &
                                              .100000000000000E+01_cd, &
                                              .110000000000000E+01_cd, &
                                              .120000000000000E+01_cd /)

real(c_double), parameter :: yuu(nuu) = (/    .000000000000000E+00_cd, &
                                             -.113953324143752E+01_cd, &
                                             -.145709859805864E+01_cd, &
                                             -.174913308002738E+01_cd, &
                                             -.202960322136630E+01_cd, &
                                             -.225202324967546E+01_cd, &
                                             -.242723053979436E+01_cd, &
                                             -.255171976467357E+01_cd, &
                                             -.260521638832322E+01_cd, &
                                             -.264397894381693E+01_cd, &
                                             -.265707884842034E+01_cd, &
                                             -.264564149400021E+01_cd, &
                                             -.260870604452106E+01_cd /)

real(c_double), parameter :: Buu(nuu) = (/   -.183757286015853E+02_cd, &
                                             -.574233124410516E+01_cd, &
                                             -.236790436375322E+01_cd, &
                                             -.307404645857774E+01_cd, &
                                             -.251104850116555E+01_cd, &
                                             -.196846462620234E+01_cd, &
                                             -.154391254686695E+01_cd, &
                                             -.846780636273251E+00_cd, &
                                             -.408540363905760E+00_cd, &
                                             -.286833282404628E+00_cd, &
                                             -.309389414590161E-06_cd, &
                                              .236958014464143E+00_cd, &
                                              .503352368511243E+00_cd /)

real(c_double), parameter :: Cuu(nuu) = (/    .830779120415016E+02_cd, &
                                              .432560615333001E+02_cd, &
                                             -.951179272978074E+01_cd, &
                                              .245037178153561E+01_cd, &
                                              .317960779258630E+01_cd, &
                                              .224623095704576E+01_cd, &
                                              .199928983630817E+01_cd, &
                                              .497202926962879E+01_cd, &
                                             -.589626545953876E+00_cd, &
                                              .180669736096520E+01_cd, &
                                              .106163236918694E+01_cd, &
                                              .130795086934864E+01_cd, &
                                              .135599267112235E+01_cd /)

real(c_double), parameter :: Duu(nuu) = (/   -.132739501694005E+03_cd, &
                                             -.175892847543603E+03_cd, &
                                              .398738817043878E+02_cd, &
                                              .243078670350231E+01_cd, &
                                             -.311125611846847E+01_cd, &
                                             -.823137069125319E+00_cd, &
                                              .990913144440207E+01_cd, &
                                             -.185388527186089E+02_cd, &
                                              .798774635639692E+01_cd, &
                                             -.248354997259420E+01_cd, &
                                              .821061667205675E+00_cd, &
                                              .160139339245701E+00_cd, &
                                              .160139339245701E+00_cd /)

! Buffer
!
! (irlast is the last entry into the pair potential and electron density
! spline tables. It is stored to speed calculations since the next value
! of r is likely close to the last one.)
!
! (ielast is the last entry into the embedding function spline table.
! It is stored to speed calculations since the next value of rho is
! likely close to the last one.)

type, bind(c) :: buffer_type
  real(c_double) :: influence_distance
  real(c_double) :: cutoff(1)
  integer(c_int) :: irlast
  integer(c_int) :: ielast
end type buffer_type

contains

!-------------------------------------------------------------------------------
!
!  Calculate pair potential phi(r)
!
!-------------------------------------------------------------------------------
subroutine calc_phi(r,phi,irlast)
implicit none

!-- Transferred variables
real(c_double), intent(in)    :: r
real(c_double), intent(out)   :: phi
integer(c_int), intent(inout) :: irlast

!-- Local variables
integer(c_int) i
real(c_double) dx

if (r .gt. model_cutoff) then
   ! Argument exceeds cutoff radius
   phi = 0.0_cd
else
   i = seval_i(npt,r,x,irlast)
   dx = r - x(i)
   phi = yv2(i) + dx*(Bv2(i) + dx*(Cv2(i) + dx*Dv2(i)))
endif

end subroutine calc_phi

!-------------------------------------------------------------------------------
!
!  Calculate pair potential phi(r) and its derivative dphi(r)
!
!-------------------------------------------------------------------------------
subroutine calc_phi_dphi(r,phi,dphi,irlast)
implicit none

!-- Transferred variables
real(c_double), intent(in)    :: r
real(c_double), intent(out)   :: phi,dphi
integer(c_int), intent(inout) :: irlast

!-- Local variables
integer(c_int) i
real(c_double) dx

if (r .gt. model_cutoff) then
   ! Argument exceeds cutoff radius
   phi    = 0.0_cd
   dphi   = 0.0_cd
else
   i = seval_i(npt,r,x,irlast)
   dx = r - x(i)
   phi  = yv2(i) + dx*(Bv2(i) + dx*(Cv2(i) + dx*Dv2(i)))
   dphi = Bv2(i) + dx*(2.0_cd*Cv2(i) + 3.0_cd*dx*Dv2(i))
endif

end subroutine calc_phi_dphi

!-------------------------------------------------------------------------------
!
!  Calculate electron density g(r)
!
!-------------------------------------------------------------------------------
subroutine calc_g(r,g,irlast)
implicit none

!-- Transferred variables
real(c_double), intent(in)    :: r
real(c_double), intent(out)   :: g
integer(c_int), intent(inout) :: irlast

!-- Local variables
integer(c_int) i
real(c_double) dx

if (r .gt. model_cutoff) then
   ! Argument exceeds cutoff radius
   g = 0.0_cd
else
   i = seval_i(npt,r,x,irlast)
   dx = r - x(i)
   g = yrh(i) + dx*(Brh(i) + dx*(Crh(i) + dx*Drh(i)))
endif

end subroutine calc_g

!-------------------------------------------------------------------------------
!
!  Calculate electron density derivative dg(r)
!
!-------------------------------------------------------------------------------
subroutine calc_dg(r,dg,irlast)
implicit none

!-- Transferred variables
real(c_double), intent(in)    :: r
real(c_double), intent(out)   :: dg
integer(c_int), intent(inout) :: irlast

!-- Local variables
integer(c_int) i
real(c_double) dx

if (r .gt. model_cutoff) then
   ! Argument exceeds cutoff radius
   dg = 0.0_cd
else
   i = seval_i(npt,r,x,irlast)
   dx = r - x(i)
   dg = Brh(i) + dx*(2.0_cd*Crh(i) + 3.0_cd*dx*Drh(i))
endif

end subroutine calc_dg

!-------------------------------------------------------------------------------
!
!  Calculate embedding function U(rho)
!
!-------------------------------------------------------------------------------
subroutine calc_U(rho,U,ielast)
implicit none

!-- Transferred variables
real(c_double), intent(in)    :: rho
real(c_double), intent(out)   :: U
integer(c_int), intent(inout) :: ielast

!-- Local variables
integer(c_int) i
real(c_double) dx

if (rho .le. rhomin) then
   ! Argument less than the minimum stored value
   U = 0.0_cd
else
   i = seval_i(nuu,rho,xuu,ielast)
   dx = rho - xuu(i)
   U = yuu(i) + dx*(Buu(i) + dx*(Cuu(i) + dx*Duu(i)))
endif

end subroutine calc_U

!-------------------------------------------------------------------------------
!
!  Calculate embedding function U(rho) and first derivative dU(rho)
!
!-------------------------------------------------------------------------------
subroutine calc_U_dU(rho,U,dU,ielast)
implicit none

!-- Transferred variables
real(c_double), intent(in)    :: rho
real(c_double), intent(out)   :: U,dU
integer(c_int), intent(inout) :: ielast

!-- Local variables
integer(c_int) i
real(c_double) dx

if (rho .le. rhomin) then
   ! Argument less than the minimum stored value
   U  = 0.0_cd
   dU = 0.0_cd
else
   i = seval_i(nuu,rho,xuu,ielast)
   dx = rho - xuu(i)
   U  = yuu(i) + dx*(Buu(i) + dx*(Cuu(i) + dx*Duu(i)))
   dU = Buu(i) + dx*(2.0_cd*Cuu(i) + 3.0_cd*dx*Duu(i))
endif

end subroutine calc_U_dU

!-------------------------------------------------------------------------------
!
!  This function performs a binary search to find the index i
!  for evaluating the cubic spline function
!
!    seval = y(i) + B(i)*(u-x(i)) + C(i)*(u-x(i))**2 + D(i)*(u-x(i))**3
!
!    where  x(i) .lt. u .lt. x(i+1), using horner's rule
!
!  if  u .lt. x(1) then  i = 1  is used.
!  if  u .ge. x(n) then  i = n  is used.
!
!  input..
!
!    n = the number of data points
!    u = the abscissa at which the spline is to be evaluated
!    x = the array of data abscissas
!    i = current value of i
!
!  if  u  is not in the same interval as the previous call, then a
!  binary search is performed to determine the proper interval.
!
!-------------------------------------------------------------------------------
integer(c_int) function seval_i(n, u, x, i)
implicit none

!--Transferred variables
integer(c_int), intent(in)    ::  n
real(c_double), intent(in)    ::  u, x(n)
integer(c_int), intent(inout) ::  i

!--Local variables
integer(c_int) j, k

if ( i .ge. n ) i = 1
if ( u .lt. x(i) ) go to 10
if ( u .le. x(i+1) ) go to 30

!  binary search
!
10 i = 1
   j = n+1
20 k = (i+j)/2
   if ( u .lt. x(k) ) j = k
   if ( u .ge. x(k) ) i = k
   if ( j .gt. i+1 ) go to 20

!  got i, return
!
30 seval_i = i
   return

end function seval_i

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
real(c_double) :: r,Rsqij,phi,dphi,g,dg,dU,U,dphieff = 0.0_cd
real(c_double) :: dphii,dUi,Ei,dphij,dUj,Ej
integer(c_int) :: i,j,jj,numnei,comp_force,comp_particleEnergy,comp_virial,comp_energy
integer(c_int) :: ierr2
real(c_double), allocatable :: rho(:),derU(:)
integer(c_int) :: irlast, ielast

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

! Initialize interpolation interval indices
irlast = 1
ielast = 1

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
if (comp_force.eq.1) force  = 0.0_cd
if (comp_virial.eq.1) virial = 0.0_cd

allocate( rho(N) )  ! pair functional electron density
rho = 0.0_cd

! EAM embedded energy deriv
if (comp_force.eq.1.or.comp_virial.eq.1) allocate( derU(N) )

!  Loop over particles in the neighbor list a first time,
!  to compute electron density (=coordination)
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

         j = nei1part(jj)                            ! get neighbor ID

         ! compute relative position vector
         !
         Rij(:) = coor(:,j) - coor(:,i)              ! distance vector between i j

         ! compute contribution to electron density
         !
         Rsqij = dot_product(Rij,Rij)                ! compute square distance
         if ( Rsqij .lt. model_cutsq ) then          ! particles are interacting?
            r = sqrt(Rsqij)                          ! compute distance
            call calc_g(r,g,irlast)                  ! compute electron density
            rho(i) = rho(i) + g                      ! accumulate electron density
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
       call calc_U_dU(rho(i),U,dU,ielast)      ! compute embedding energy
                                               !   and its derivative
       derU(i) = dU                            ! store dU for later use
     else
       call calc_U(rho(i),U,ielast)            ! compute just embedding energy
     endif

     ! accumulate the embedding energy contribution
     !
     ! Assuming U(rho=0) = 0.0_cd
     !
     if (comp_particleEnergy.eq.1) then        ! accumulate embedding energy contribution
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

         j = nei1part(jj)                            ! get neighbor ID

         ! compute relative position vector
         !
         Rij(:) = coor(:,j) - coor(:,i)              ! distance vector between i j

         ! compute energy and forces
         !
         Rsqij = dot_product(Rij,Rij)                ! compute square distance
         if ( Rsqij .lt. model_cutsq ) then          ! particles are interacting?

            r = sqrt(Rsqij)                          ! compute distance
            if (comp_force.eq.1.or.comp_virial.eq.1) then
               call calc_phi_dphi(r,phi,dphi,irlast) ! compute pair potential
                                                     !   and it derivative
               call calc_dg(r,dg,irlast)             ! compute elect dens first deriv
               dphii  = 0.5_cd*dphi
               dphij  = 0.0_cd
               dUi    = derU(i)*dg
               dUj    = 0.0_cd
               dphieff = dphii + dphij + dUi + dUj
            else
               call calc_phi(r,phi,irlast)           ! compute just pair potential
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
            if (comp_force.eq.1) then                  ! Ei contribution
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
ierr = 0

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
  call kim_model_refresh_set_neighbor_list_cutoffs_pointer( &
    model_refresh_handle, 1, buf%cutoff)

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

end module Glue_Ercolessi_Adams_Al

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
use Glue_Ercolessi_Adams_Al
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
  kim_species_name_al, speccode, ierr2)
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
end subroutine create

