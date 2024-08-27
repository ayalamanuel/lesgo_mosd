!!
!!  Copyright (C) 2009-2017  Johns Hopkins University
!!
!!  This file is part of lesgo.
!!
!!  lesgo is free software: you can redistribute it and/or modify
!!  it under the terms of the GNU General Public License as published by
!!  the Free Software Foundation, either version 3 of the License, or
!!  (at your option) any later version.
!!
!!  lesgo is distributed in the hope that it will be useful,
!!  but WITHOUT ANY WARRANTY; without even the implied warranty of
!!  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!!  GNU General Public License for more details.
!!
!!  You should have received a copy of the GNU General Public License
!!  along with lesgo.  If not, see <http://www.gnu.org/licenses/>.
!!

!*******************************************************************************
subroutine wallstress
!*******************************************************************************
!
! This subroutine calculates the wall stress txz, tyz (w-nodes) and dudz,
! dvdz (w-nodes) at the first z-location k = 1. The wall stress is calculated
! depending on lower boundary condition lbc_mom. This subroutine should only
! be called after ensuring coord==0
!
! Options for lbc_mom:
!   0 - stress free
!       txz, tyz, dudz, and dvdz are all 0
!
!   1 - DNS wall boundary conditions
!       calculates wall stress values from the first grid point
!
!   2 - Equilibirum wall model
!       See John D. Albertson's dissertation, eqns (2.46)-(2.52)
!       Also see E. Bou-Zeid, C. Meneveau & M.B. Parlange, "A scale-dependent
!           Lagrangian dynamic model for large eddy simulation of complex
!           turbulent flows" (2005) -- Appendix
!
!   3 - Integral wall model
!       See X.I.A. Yang, J. Sadique, R. Mittal & C. Meneveau, "Integral wall
!           model for large eddy simulations of wall-bounded turbulent flows." (2015)
!
!   4 - MOving Surface Gradient (MOSD) model 
!       Combination of stress from the equilibrium wall model and from
!	the potential flow-based model. See (future Ayala paper) 
!

use types, only : rprec
use param, only : lbc_mom
use param, only : ubc_mom, coord, nproc, nz ! these necessary only for upper bc
use messages, only : error
use iwmles, only : iwm_wallstress
use sim_param, only : txz, tyz, dudz, dvdz
implicit none
character(*), parameter :: sub_name = 'wallstress'

! Lower boundary condition
if (coord == 0) then
    select case (lbc_mom)
        ! Stress free
        case (0)
            call ws_free_lbc

        ! DNS wall
        case (1)
            call ws_dns_lbc

        ! Equilibrium wall model
        case (2)
            call ws_equilibrium_lbc

        ! Integral wall model (not implemented for top wall)
        case (3)
            call iwm_wallstress()

        ! Equilibrium Moving Surface Gradient (EMSG)
        case(4)
            call ws_mosd_lbc

        ! Otherwise, invalid
        case default
            call error (sub_name, 'invalid lbc_mom')
    end select
end if

if (coord == nproc-1) then
    select case (ubc_mom)
        ! Stress free
        case (0)
            call ws_free_ubc

        ! DNS wall
        case (1)
            call ws_dns_ubc

        ! Equilibrium wall model
        case (2)
            call ws_equilibrium_ubc

        ! Integral wall model (not implemented for top wall)
        case (3)
            call error(sub_name, 'invalid ubc_mom')

        ! Otherwise, invalid
        case default
            call error(sub_name, 'invalid ubc_mom')
    end select
end if

contains

!*******************************************************************************
subroutine ws_free_lbc
!*******************************************************************************
implicit none

txz(:, :, 1) = 0._rprec
tyz(:, :, 1) = 0._rprec
dudz(:, :, 1) = 0._rprec
dvdz(:, :, 1) = 0._rprec

end subroutine ws_free_lbc

!*******************************************************************************
subroutine ws_free_ubc
!*******************************************************************************
implicit none

txz(:, :,nz) = 0._rprec
tyz(:, :,nz) = 0._rprec
dudz(:,:,nz) = 0._rprec
dvdz(:,:,nz) = 0._rprec

end subroutine ws_free_ubc

!*******************************************************************************
subroutine ws_dns_lbc
!*******************************************************************************
use param, only : nx, ny, nu_molec, z_i, u_star, dz
use param, only : ubot
use sim_param , only : u, v
implicit none
integer :: i, j

do j = 1, ny
    do i = 1, nx
        dudz(i,j,1) = ( u(i,j,1) - ubot ) / (0.5_rprec*dz)
        dvdz(i,j,1) = v(i,j,1) / (0.5_rprec*dz)
        txz(i,j,1) = -nu_molec/(z_i*u_star)*dudz(i,j,1)
        tyz(i,j,1) = -nu_molec/(z_i*u_star)*dvdz(i,j,1)
    end do
end do

end subroutine ws_dns_lbc

!*******************************************************************************
subroutine ws_dns_ubc
!*******************************************************************************
use param, only : nx, ny, nu_molec, z_i, u_star, dz
use param, only : utop
use sim_param , only : u, v
implicit none
integer :: i, j

do j = 1, ny
    do i = 1, nx
        dudz(i,j,nz) = ( utop - u(i,j,nz-1) ) / (0.5_rprec*dz)
        dvdz(i,j,nz) = -v(i,j,nz-1) / (0.5_rprec*dz)
        txz(i,j,nz) = -nu_molec/(z_i*u_star)*dudz(i,j,nz)
        tyz(i,j,nz) = -nu_molec/(z_i*u_star)*dvdz(i,j,nz)
    end do
end do

end subroutine ws_dns_ubc

!*******************************************************************************
subroutine ws_equilibrium_lbc
!*******************************************************************************
use param, only : dz, ld, nx, ny, vonk, zo, amp, wave_freq, wave_n, total_time, dx, worb
use sim_param, only : u, v, ustar_lbc, w_orb
use test_filtermodule
#ifdef PPSCALARS
use scalars, only : obukhov, phi_m
#endif

implicit none

integer :: i, j
real(rprec), dimension(nx, ny) :: denom, u_avg
real(rprec), dimension(ld, ny) :: u1, v1, x_grid
real(rprec) :: const

u1 = u(:,:,1)
v1 = v(:,:,1)
call test_filter(u1)
call test_filter(v1)
denom = log(0.5_rprec*dz/zo)
u_avg = sqrt(u1(1:nx,1:ny)**2+v1(1:nx,1:ny)**2)
#ifdef PPSCALARS
call obukhov(u_avg)
#else
ustar_lbc = u_avg*vonk/denom
#endif

if (worb) then
do j = 1,ny
        do i = 1,ld  !changed this from ld to nx just for trial of 4dz grid point velocity (10/12/22) 12pm MA
        x_grid(i,j) = (i-1)*dx
        w_orb(i,j) = amp*wave_freq*sin(wave_n*x_grid(i,j) - wave_freq*total_time)
        end do
end do
end if 


do j = 1, ny
    do i = 1, nx
        const = -(ustar_lbc(i,j)**2)/u_avg(i,j)
        txz(i,j,1) = const*u1(i,j)
        tyz(i,j,1) = const*v1(i,j)
        !this is as in Moeng 84
#ifdef PPSCALARS
        dudz(i,j,1) = ustar_lbc(i,j)/(0.5_rprec*dz*vonK)*u(i,j,1)/u_avg(i,j)   &
            * phi_m(i,j)
        dvdz(i,j,1) = ustar_lbc(i,j)/(0.5_rprec*dz*vonK)*v(i,j,1)/u_avg(i,j)   &
            * phi_m(i,j)
#else
        dudz(i,j,1) = ustar_lbc(i,j)/(0.5_rprec*dz*vonK)*u(i,j,1)/u_avg(i,j)
        dvdz(i,j,1) = ustar_lbc(i,j)/(0.5_rprec*dz*vonK)*v(i,j,1)/u_avg(i,j)
#endif
        dudz(i,j,1) = merge(0._rprec,dudz(i,j,1),u(i,j,1).eq.0._rprec)
        dvdz(i,j,1) = merge(0._rprec,dvdz(i,j,1),v(i,j,1).eq.0._rprec)
    end do
end do

end subroutine ws_equilibrium_lbc


!*******************************************************************************
subroutine ws_mosd_lbc
!*******************************************************************************
use param, only : dx, dz, ld, nx, ny, vonk, zo, dx, total_time, amp, wave_n,    &
        wave_freq, cd, pi, z_i, u_star, nu_molec, wave_spec, dt, total_time_dim,&
        smooth_eqm, zgrid_match, dy, dt, amp2, wave_n2, wave_freq2, L_x,L_y,    &
        kp_spec
use sim_param, only : u, v, eqmxz, eqmyz, wpmxz, wpmyz, eta, detadx, detady,    &
        detadt, u_orb, w_orb, ddtw_orb, ur_eqm, vr_eqm, Cx_wave, Cy_wave,       &
        ur_wpm, vr_wpm, u_LES, eta_spectrum, omega_wave, eta_hat_o,             &
        grad_eta_mag, dgrad_etadt, CdotGradGradeta, unsxz, unsyz, ur_mag_wpm,   &
        H_wpm, nx_wpm, ny_wpm, uns_convec_x1, uns_convec_x2, uns_convec_y1,     &
        uns_convec_y2, corr_1, corr_2, corr_3, corr_4, corr_5, corr_6,          &
        eta_hat_o_filter, eta_filter

use test_filtermodule
use fft
use derivatives, only : ddx_fd, ddy_fd
use grid_m
use functions, only : cell_indx
use wave_spectrum, only : frequency_shift
implicit none

integer :: jx, jy, dz_match_ind
real(rprec) :: dz_match, const_time, z_dx_low, z_dx_up, z_diff, threshold_wave, &
               k_min, L_x_wave, L_y_wave
real(rprec), pointer, dimension(:) :: z
real(rprec), dimension(ld, ny) :: u1, v1, u_dx, v_dx, x_grid, u_eqm, v_eqm
real(rprec), dimension(nx, ny) :: u_avg, beta_1, beta_2, beta_11, beta_21, &
                                  Re_delta, Re_fit, Re_delta1, Re_fit1, A1, A, u_orb_1,  &
                                  eta_new, eta_1, eta_filter_1                      &
                                  C_mag, ur_mag_eqm, H_arg,              &
                                  deta2dx2, deta2dy2, detadxdy, detadydx,      &
                                  grad_eta_mag_new, uns_mosd
complex(rprec), dimension(nx,ny) :: eta_hat, eta_shift, u_orb_hat, u_orb_shift, &
                                    eta_hat_filter, eta_filter_shift
complex(rprec), dimension(nx/2+1,ny) :: eta_hat_2, u_orb_hat_2, eta_hat_filter_2
complex(kind=8) :: ii
ii = cmplx(0.0d0,1.0d0)
nullify(z)
z => grid % z

do jy = 1,ny
do jx= 1,ld
x_grid(jx,jy) = (jx-1)*dx
end do
end do 

! Calculating the velocity input for EQM, this is set up to manage any matching location 
! specified in lesgo.conf. (e.g. zgrid_match=1.5 is at 2nd grid point and 2.5 is at 3rd)
dz_match = zgrid_match*dz
dz_match_ind = ceiling(zgrid_match)
u_eqm = u(1:ld,1:ny,dz_match_ind)
v_eqm = v(1:ld,1:ny,dz_match_ind)

! Calculating the velocity input for WPM at a height dx. Not using the linear_interp 
! because it does not work. The velocity is also filtered at delta scale.
! NOTE: ensure that height dx lies in the coord = 0.
z_dx_low = cell_indx('k',dz,dx)
z_dx_up = z_dx_low + 1
z_diff = dx - z(z_dx_low)
u_dx = u(1:ld,1:ny,z_dx_low) + (u(1:ld,1:ny,z_dx_up) - u(1:ld,1:ny,z_dx_low))*z_diff/dz
v_dx = v(1:ld,1:ny,z_dx_low) + (v(1:ld,1:ny,z_dx_up) - v(1:ld,1:ny,z_dx_low))*z_diff/dz
call test_filter(u_dx)
call test_filter(v_dx)

! Calculating the velocity input for dudz, dvdz at 1st grid point. Velocity is filtered
! and velocity magnitude is calculated.
u1 = u(1:ld,1:ny,1)
v1 = v(1:ld,1:ny,1)
call test_filter(u1)
call test_filter(v1)
u_avg = sqrt(u1(1:nx,1:ny)**2+v1(1:nx,1:ny)**2)

! Calculating the wave surface, spatial and time derivatives of surface, the 
! orbital velocities, and phase velocity of wave.
if (wave_spec) then
        do jy = 1,ny
                do jx = 1,nx
        eta_hat(jx,jy) = eta_hat_o(jx,jy)*exp(-ii*omega_wave(jx,jy)*total_time_dim)
        eta_hat_filter(jx,jy) = eta_hat_o_filter(jx,jy)*exp(-ii*omega_wave(jx,jy)*total_time_dim)
        u_orb_hat(jx,jy) = eta_hat(jx,jy)*omega_wave(jx,jy)
                end do
        end do
        call frequency_shift(eta_hat,eta_shift)
        call frequency_shift(eta_hat_filter,eta_filter_shift)
        call frequency_shift(u_orb_hat,u_orb_shift)

        eta_hat_2(:,:) = eta_shift(:nx/2+1,1:ny)
        eta_hat_filter_2(:,:) = eta_filter_shift(:nx/2+1,1:ny)
        u_orb_hat_2(:,:) = u_orb_shift(:nx/2+1,1:ny)

        call dfftw_execute_dft_c2r(back_wave, eta_hat_2, eta_1)
        call dfftw_execute_dft_c2r(back_wave, eta_hat_filter_2, eta_filter_1)
        call dfftw_execute_dft_c2r(back_wave, u_orb_hat_2, u_orb_1)

        ! Since the surface spectrum and orbital velocities are calculated using dimensions,
        ! a normalization must be applied. Because of the FFT, the surface spectrum has to be 
        ! devided by 2
        eta_new(1:nx,1:ny) = eta_1(1:nx,1:ny)/2.0_rprec/z_i
        eta_filter(1:nx,1:ny) = eta_filter_1(1:nx,1:ny)/2.0_rprec
        u_orb(1:nx,1:ny) = u_orb_1(1:nx,1:ny)/2.0_rprec/u_star
        
        ! Calculating the time derivative and spatial derivatives  of the surface
        const_time = 1.0_rprec/dt
        detadt(1:nx,1:ny) = const_time*(eta_new(1:nx,1:ny) - eta(1:nx,1:ny))
        eta(1:nx,1:ny) = eta_new(1:nx,1:ny)
        call ddx_fd(eta, detadx)
        call ddy_fd(eta, detady)
        detadx(nx,:) = -eta(nx,:)/dx
        detady(:,ny) = -eta(:,ny)/dy

        call ddx_fd(detadx,deta2dx2)
        call ddy_fd(detady,deta2dy2)
        deta2dx2(nx,:) = -detadx(nx,:)/dx
        deta2dy2(:,ny) = -detady(:,ny)/dy

        call ddy_fd(detadx,detadxdy)
        call ddx_fd(detady,detadydx)
        detadxdy(:,ny) = -detadx(:,ny)/dy
        detadydx(nx,:) = -detady(nx,:)/dx
else
        eta(1:nx,1:ny) = amp*cos(wave_n*x_grid(1:nx,1:ny) - wave_freq*total_time) + &
                   amp2*cos(wave_n2*x_grid(1:nx,1:ny) - wave_freq2*total_time)
        detadx(1:nx,1:ny) = amp*wave_n*sin(wave_freq*total_time - wave_n*x_grid(1:nx,1:ny))&
                      + amp2*wave_n2*sin(wave_freq2*total_time - wave_n2*x_grid(1:nx,1:ny))
        detady(1:nx,1:ny) = 0.0_rprec
        detadt(1:nx,1:ny) = amp*wave_freq*sin(wave_n*x_grid(1:nx,1:ny)   - wave_freq*total_time)&
                      + amp2*wave_freq2*sin(wave_n2*x_grid(1:nx,1:ny)  - wave_freq2*total_time)
        u_orb(1:nx,1:ny) = amp*wave_freq*cos(wave_n*x_grid(1:nx,1:ny)  - wave_freq*total_time)&
                        + amp2*wave_freq2*cos(wave_n2*x_grid(1:nx,1:ny) - wave_freq2*total_time)
        w_orb(1:nx,1:ny) = 0.0_rprec
        deta2dx2(1:nx,1:ny) = -amp*wave_n**2*cos(-wave_freq*total_time + wave_n*x_grid(1:nx,1:ny))&
                      - amp2*wave_n2**2*cos(-wave_freq2*total_time + wave_n2*x_grid(1:nx,1:ny))
        deta2dy2(1:nx,1:ny) = 0.0_rprec
        detadxdy(1:nx,1:ny) = 0.0_rprec
        ddtw_orb(1:nx,1:ny) = -amp*wave_freq*wave_freq*cos(wave_n*x_grid(1:nx,1:ny) - wave_freq*total_time)
end if


! Calculating the general phase velocity and the magnitude of surface gradient 
! for WPM model
grad_eta_mag_new(1:nx,1:ny) = sqrt((detadx(1:nx,1:ny))**2 + (detady(1:nx,1:ny))**2)
dgrad_etadt(1:nx,1:ny) = (1.0_rprec/dt)*(grad_eta_mag_new(1:nx,1:ny) - grad_eta_mag(1:nx,1:ny))
grad_eta_mag(1:nx,1:ny) = grad_eta_mag_new(1:nx,1:ny)

!threshold_wave = 343.0_rprec/u_star
!L_x_wave = L_x*z_i
!L_y_wave = L_y*z_i
!k_min = sqrt((2*pi/L_x_wave)**2 + (2*pi/L_y_wave)**2)
k_min = 0.25_rprec*(kp_spec/z_i)
threshold_wave = sqrt(9.81_rprec/k_min)/u_star

Cx_wave(1:nx,1:ny) = -detadt(1:nx,1:ny)*detadx(1:nx,1:ny)*(1/grad_eta_mag(1:nx,1:ny)**2)
Cx_wave(1:nx,1:ny) = min(max(Cx_wave(1:nx,1:ny),-threshold_wave),threshold_wave)
Cy_wave(1:nx,1:ny) = -detadt(1:nx,1:ny)*detady(1:nx,1:ny)*(1/grad_eta_mag(1:nx,1:ny)**2)
Cy_wave(1:nx,1:ny) = min(max(Cx_wave(1:nx,1:ny),-threshold_wave),threshold_wave)
C_mag(1:nx,1:ny) = sqrt(Cx_wave(1:nx,1:ny)**2 + Cy_wave(1:nx,1:ny)**2)

! Calcualting the realtive velocities used in EQM and WPM models.
! NOTE: only the relative velocity for EQM is test filtered. u_dx and v_dx
! are already filtered and we CANNOT filter Cx_wave, Cy_wave because they
! won't cancel out with the surface gradients specially if wave_spec = true
ur_wpm(1:nx,1:ny) = u_dx(1:nx,1:ny) - Cx_wave(1:nx,1:ny)
vr_wpm(1:nx,1:ny) = v_dx(1:nx,1:ny) - Cy_wave(1:nx,1:ny)

ur_eqm(1:nx,1:ny) = u_eqm(1:nx,1:ny) - u_orb(1:nx,1:ny)
vr_eqm(1:nx,1:ny) = v_eqm(1:nx,1:ny)
call test_filter(ur_eqm)
call test_filter(vr_eqm)

! Calculating the parameters for the EQM model (see Meneveau JoT 2020 or Ayala BLM 2024)
! We calcualte 2 sets because the second is used for the calculation of dudz,dvdz
u_LES(1:nx,1:ny)  = sqrt((ur_eqm(1:nx,1:ny))**2 + (vr_eqm(1:nx,1:ny))**2)
!Re_delta(1:nx,1:ny) = ((u_LES(1:nx,1:ny)*u_star)*((dz_match-eta_spectrum(1:nx,1:ny))*z_i))/nu_molec
Re_delta(1:nx,1:ny) = ((u_LES(1:nx,1:ny)*u_star)*((dz_match)*z_i))/nu_molec
beta_1(1:nx,1:ny) = (1_rprec + 0.155_rprec*Re_delta(1:nx,1:ny)**(-0.03_rprec))**(-1_rprec)
beta_2(1:nx,1:ny) = 1.7_rprec - (1_rprec + 36_rprec*Re_delta(1:nx,1:ny)**(-0.75_rprec))**(-1_rprec)
Re_fit(1:nx,1:ny) = (0.005_rprec)**(beta_1(1:nx,1:ny)-0.5_rprec)*Re_delta(1:nx,1:ny)**(beta_1(1:nx,1:ny))* &
                    (1_rprec + (0.005_rprec*Re_delta(1:nx,1:ny))**(-beta_2(1:nx,1:ny)))**((beta_1(1:nx,1:ny) &
                    -0.5_rprec)/beta_2(1:nx,1:ny))

Re_delta1(1:nx,1:ny) = ((u_avg(1:nx,1:ny)*u_star)*(0.5_rprec*dz*z_i))/nu_molec
beta_11(1:nx,1:ny) = (1_rprec + 0.155_rprec*Re_delta1(1:nx,1:ny)**(-0.03_rprec))**(-1_rprec)
beta_21(1:nx,1:ny) = 1.7_rprec - (1_rprec + 36_rprec*Re_delta1(1:nx,1:ny)**(-0.75_rprec))**(-1_rprec)
Re_fit1(1:nx,1:ny) = (0.005_rprec)**(beta_11(1:nx,1:ny)-0.5_rprec)*Re_delta1(1:nx,1:ny)**(beta_11(1:nx,1:ny))* &
                    (1_rprec + (0.005_rprec*Re_delta1(1:nx,1:ny))**(-beta_21(1:nx,1:ny)))**((beta_11(1:nx,1:ny) &
                    -0.5_rprec)/beta_21(1:nx,1:ny))

! Calculating the WPM component
nx_wpm(1:nx,1:ny) = detadx(1:nx,1:ny)/grad_eta_mag(1:nx,1:ny)
ny_wpm(1:nx,1:ny) = detady(1:nx,1:ny)/grad_eta_mag(1:nx,1:ny)
ur_mag_wpm(1:nx,1:ny) = sqrt((ur_wpm(1:nx,1:ny)*nx_wpm(1:nx,1:ny))**2 + (vr_wpm(1:nx,1:ny)*ny_wpm(1:nx,1:ny))**2)
H_arg(1:nx,1:ny) = (ur_wpm(1:nx,1:ny)*detadx(1:nx,1:ny) + vr_wpm(1:nx,1:ny)*detady(1:nx,1:ny))
H_wpm(1:nx,1:ny) = (H_arg(1:nx,1:ny) + abs(H_arg(1:nx,1:ny)))/(2*H_arg(1:nx,1:ny))


corr_1(1:nx,1:ny) = ur_mag_wpm(1:nx,1:ny)**2*ny_wpm(1:nx,1:ny)
corr_2(1:nx,1:ny) = grad_eta_mag(1:nx,1:ny)**2*ny_wpm(1:nx,1:ny)
corr_3(1:nx,1:ny) = H_wpm(1:nx,1:ny)*ny_wpm(1:nx,1:ny)
corr_4(1:nx,1:ny) = ur_mag_wpm(1:nx,1:ny)**2*H_wpm(1:nx,1:ny)*ny_wpm(1:nx,1:ny)
corr_5(1:nx,1:ny) = grad_eta_mag(1:nx,1:ny)**2*H_wpm(1:nx,1:ny)*ny_wpm(1:nx,1:ny)
corr_6(1:nx,1:ny) = (ur_mag_wpm(1:nx,1:ny)**2)*(grad_eta_mag(1:nx,1:ny)**2)*ny_wpm(1:nx,1:ny)


wpmxz(1:nx,1:ny)  = -(1/pi)*(ur_mag_wpm(1:nx,1:ny)**2)*nx_wpm(1:nx,1:ny)*(grad_eta_mag(1:nx,1:ny)**2)*H_wpm(1:nx,1:ny)
wpmyz(1:nx,1:ny)  = -(1/pi)*(ur_mag_wpm(1:nx,1:ny)**2)*ny_wpm(1:nx,1:ny)*(grad_eta_mag(1:nx,1:ny)**2)*H_wpm(1:nx,1:ny)

! Calculating the EQM component
ur_mag_eqm(1:nx,1:ny) = sqrt((ur_eqm(1:nx,1:ny))**2 + (vr_eqm(1:nx,1:ny))**2)
if (smooth_eqm) then 
        eqmxz(1:nx,1:ny) = -(Re_fit(1:nx,1:ny)/Re_delta(1:nx,1:ny))**2 * u_eqm(1:nx,1:ny)*u_LES(1:nx,1:ny) 
        eqmyz(1:nx,1:ny) = -(Re_fit(1:nx,1:ny)/Re_delta(1:nx,1:ny))**2 * v_eqm(1:nx,1:ny)*u_LES(1:nx,1:ny)
        dudz(1:nx,1:ny,1) = 1/(0.5_rprec*dz*vonk)*(Re_fit1(1:nx,1:ny)/Re_delta1(1:nx,1:ny))*u(1:nx,1:ny,1)
        dvdz(1:nx,1:ny,1) = 1/(0.5_rprec*dz*vonk)*(Re_fit1(1:nx,1:ny)/Re_delta1(1:nx,1:ny))*v(1:nx,1:ny,1)
else
        !A(1:nx,1:ny) = ((Re_fit(1:nx,1:ny)/Re_delta(1:nx,1:ny))**(6.0_rprec) +   &
        !(2.5_rprec*log((dz_match-eta_spectrum(1:nx,1:ny))/zo))**(-6.0_rprec))**(0.33333_rprec) 
        !eqmxz(1:nx,1:ny) = - A(1:nx,1:ny) * u_eqm(1:nx,1:ny)*u_LES(1:nx,1:ny)
        !eqmyz(1:nx,1:ny) = - A(1:nx,1:ny) * v_eqm(1:nx,1:ny)*u_LES(1:nx,1:ny)

        A(1:nx,1:ny) = ((Re_fit(1:nx,1:ny)/Re_delta(1:nx,1:ny))**(6.0_rprec) +   &
        (2.5_rprec*log((dz_match)/zo))**(-6.0_rprec))**(0.33333_rprec)
        eqmxz(1:nx,1:ny) = - A(1:nx,1:ny) * u_eqm(1:nx,1:ny)*u_LES(1:nx,1:ny)
        eqmyz(1:nx,1:ny) = - A(1:nx,1:ny) * v_eqm(1:nx,1:ny)*u_LES(1:nx,1:ny)

        A1(1:nx,1:ny) = ((Re_fit1(1:nx,1:ny)/Re_delta1(1:nx,1:ny))**(6.0_rprec) + & 
        (2.5_rprec*log((0.5_rprec*dz)/zo))**(-6.0_rprec))**(0.16666_rprec)
        dudz(1:nx,1:ny,1) = A1(1:nx,1:ny)/(vonk*0.5_rprec*dz)*u(1:nx,1:ny,1)
        dvdz(1:nx,1:ny,1) = A1(1:nx,1:ny)/(vonk*0.5_rprec*dz)*v(1:nx,1:ny,1)
end if

! Calculating the unsteady term for MOSD
CdotGradGradeta(1:nx,1:ny) = Cx_wave(1:nx,1:ny)*(detadx(1:nx,1:ny)*deta2dx2(1:nx,1:ny) + &
detady(1:nx,1:ny)*detadxdy(1:nx,1:ny))/grad_eta_mag(1:nx,1:ny) + Cy_wave(1:nx,1:ny)*(detadx(1:nx,1:ny)*detadydx(1:nx,1:ny) + &
detady(1:nx,1:ny)*deta2dy2(1:nx,1:ny))/grad_eta_mag(1:nx,1:ny)

!corr_1(1:nx,1:ny) = Cx_wave(1:nx,1:ny)/grad_eta_mag(1:nx,1:ny)
!corr_2(1:nx,1:ny) = Cx_wave(1:nx,1:ny)*deta2dx2(1:nx,1:ny)
!corr_3(1:nx,1:ny) = Cx_wave(1:nx,1:ny)*detadx(1:nx,1:ny)
!corr_4(1:nx,1:ny) = deta2dx2(1:nx,1:ny)/grad_eta_mag(1:nx,1:ny)
!corr_5(1:nx,1:ny) = detadx(1:nx,1:ny)/grad_eta_mag(1:nx,1:ny)
!corr_6(1:nx,1:ny) = deta2dx2(1:nx,1:ny)*detadx(1:nx,1:ny)


uns_mosd(1:nx,1:ny) = 1/(4*pi)*ur_mag_wpm(1:nx,1:ny)*sqrt(dx**2 + dy**2)*grad_eta_mag(1:nx,1:ny)*(dgrad_etadt(1:nx,1:ny) +  &
CdotGradGradeta(1:nx,1:ny))

unsxz(1:nx,1:ny) =  cd*uns_mosd(1:nx,1:ny)*nx_wpm(1:nx,1:ny)*H_wpm(1:nx,1:ny)
unsyz(1:nx,1:ny) =  cd*uns_mosd(1:nx,1:ny)*ny_wpm(1:nx,1:ny)*H_wpm(1:nx,1:ny)

txz(1:nx,1:ny,1) = eqmxz(1:nx,1:ny) +  wpmxz(1:nx,1:ny) + unsxz(1:nx,1:ny)
tyz(1:nx,1:ny,1) = eqmyz(1:nx,1:ny) +  wpmyz(1:nx,1:ny) + unsyz(1:nx,1:ny)

dudz(1:nx,1:ny,1) = merge(0._rprec,dudz(1:nx,1:ny,1),u(1:nx,1:ny,1).eq.0._rprec)
dvdz(1:nx,1:ny,1) = merge(0._rprec,dvdz(1:nx,1:ny,1),v(1:nx,1:ny,1).eq.0._rprec)


end subroutine ws_mosd_lbc


!*******************************************************************************
subroutine ws_equilibrium_ubc
!*******************************************************************************
use param, only : dz, ld, nx, ny, vonk, zo
use sim_param, only : u, v
use test_filtermodule
implicit none
integer :: i, j
real(rprec), dimension(nx, ny) :: denom, u_avg, ustar
real(rprec), dimension(ld, ny) :: u1, v1
real(rprec) :: const


u1 = u(:,:,nz-1)
v1 = v(:,:,nz-1)
call test_filter(u1)
call test_filter(v1)
denom = log(0.5_rprec*dz/zo)
u_avg = sqrt(u1(1:nx,1:ny)**2+v1(1:nx,1:ny)**2)
ustar = u_avg*vonk/denom

do j = 1, ny
    do i = 1, nx
        const = (ustar(i,j)**2)/u_avg(i,j) ! diff sign for upper b.c.
        txz(i,j,nz) = const*u1(i,j)
        tyz(i,j,nz) = const*v1(i,j)
        !this is as in Moeng 84
        dudz(i,j,nz) = -ustar(i,j)/(0.5_rprec*dz*vonK)*u(i,j,nz-1)/u_avg(i,j)
        dvdz(i,j,nz) = -ustar(i,j)/(0.5_rprec*dz*vonK)*v(i,j,nz-1)/u_avg(i,j)
        dudz(i,j,nz) = merge(0._rprec,dudz(i,j,nz),u(i,j,nz-1).eq.0._rprec)
        dvdz(i,j,nz) = merge(0._rprec,dvdz(i,j,nz),v(i,j,nz-1).eq.0._rprec)
    end do
end do

end subroutine ws_equilibrium_ubc

end subroutine wallstress
