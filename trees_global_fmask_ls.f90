!!
!!  Copyright (C) 2009-2013  Johns Hopkins University
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

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module trees_global_fmask_ls
use param, only : ld, nx, ny, nz, dx, dy, dz
use trees_base_ls
implicit none

save
private

public :: global_fmask
public :: global_fmask_init
public :: calc_global_fmask_ta
public :: read_global_fmask

character (*), parameter :: mod_name = 'trees_global_fmask_ls'
character (*), parameter :: raw_suffix = '.out'
character (*), parameter :: txt_suffix = '.dat'
character (*), parameter :: MPI_suffix = '.c'
character (*), parameter :: gfmask_base = 'global_fmask'

integer, parameter :: fratiomax = 2
    !--truncation of Gaussian filter at fratiomax*delta
integer, parameter :: idelta = 4
    !--filter width delta divided by dx: delta/dx

logical, parameter :: do_filter_global_fmask = .true.

!--these are input to calc_global_fmask_ta
integer :: np
logical :: MPI_split

!real (rp) :: global_fmask(ld, ny, nz)  !--experimental
real(rp), allocatable, dimension(:,:,:) :: global_fmask
    !--nonzero where unres force is to be applied and contents may be filtered
    !--could dimension in a smarter way to save space

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
contains
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine global_fmask_init()
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
implicit none

allocate( global_fmask(ld,ny,nz) )

return
end subroutine global_fmask_init


!--this is usually called from trees_pre_ls
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine calc_global_fmask_ta ( nnp, MMPI_split )
implicit none

integer, intent (in) :: nnp
logical, intent (in) :: MMPI_split

character (*), parameter :: sub_name = mod_name // '.calc_global_fmask_ta'

integer :: i

real (rp) :: r

type (branch_type) :: br  !--to simplify argument passing

!---------------------------------------------------------------------

!--set module copies
np = nnp
MPI_split = MMPI_split

global_fmask = 0.0_rp

do i = 1, n_tree

    r = tree_array(i) % ratio

    br % l = tree_array(i) % l
    br % d = tree_array(i) % d

    !--could make this more detailed: add_base, etc.
    if ( add_cap ) br % l = br % l + 0.5_rp * br % d

    br % gen = tree_array(i) % trunk % gen  !--0 with current convention

    br % x0 = tree_array(i) % trunk % x0
    br % abs_dir = tree_array(i) % trunk % abs_dir
    br % rel_dir = tree_array(i) % trunk % rel_dir
    br % x_hat = tree_array(i) % trunk % x_hat
    br % y_hat = tree_array(i) % trunk % y_hat
    br % z_hat = tree_array(i) % trunk % z_hat

    br % n_sub_branch = tree_array(i) % trunk % n_sub_branch

    call calc_global_fmask_br ( r, br )

end do

if ( do_filter_global_fmask ) call filter_global_fmask ()

!--module copies of MPI_suffix must be set before calling these routines
call write_global_fmask ()
call write_fmt_global_fmask ()


end subroutine calc_global_fmask_ta

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!--this will be used in simulation run-time code, therefore needs to
!  know some MPI stuff from param
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine read_global_fmask ()
use param, only : coord, nproc
implicit none

integer, parameter :: lun = 1

character (128) :: fname

!---------------------------------------------------------------------

#ifdef PPMPI

    write ( fname, '(a,a,a,i0)' ) gfmask_base, raw_suffix, MPI_suffix, coord

#else

    write ( fname, '(a,a)' ) gfmask_base, raw_suffix

#endif

open ( lun, file=fname, action='read', position='rewind', &
       form='unformatted' )
read (lun) global_fmask
close (lun)

end subroutine read_global_fmask

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine write_global_fmask ()
implicit none

character (*), parameter :: sub_name = mod_name // '.write_global_fmask'

integer, parameter :: lun = 1

character (128) :: fname

!logical :: do_fmt

integer :: ip
integer :: lbz, ubz

!---------------------------------------------------------------------

if ( MPI_split ) then

    do ip = 0, np-1

        write ( fname, '(a,a,a,i0)' ) gfmask_base, raw_suffix, MPI_suffix, ip

        open ( lun, file=fname, action='write', position='rewind',  &
               form='unformatted' )

        lbz = ip * (nz - 1) / np + 1  !--1 level (local)
        ubz = lbz + (nz - 1) / np  !--nz level (local)

        call mesg ( sub_name, '(ip,lbz,ubz)=', (/ ip, lbz, ubz /) )

        write ( lun ) global_fmask(:, :, lbz:ubz)

        close ( lun )

    end do

else

    write ( fname, '(a,a)' ) gfmask_base, raw_suffix
    open ( lun, file=fname, action='write', position='rewind',  &
           form='unformatted' )
    write ( lun ) global_fmask
    close ( lun )

end if

end subroutine write_global_fmask

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine write_fmt_global_fmask ()
implicit none

integer, parameter :: lun = 1

character (32) :: fmt
!character (64) :: fname

integer :: i, j, k

real (rp) :: x(nd)

!---------------------------------------------------------------------

open ( lun, file=gfmask_base // txt_suffix, action='write',  &
       position='rewind' )
write ( lun, '(a)' ) 'variables = "x" "y" "z" "global_fmask"'
write ( lun, '(3(a,i0))' ) 'zone, f=point, i=', nx, ', j=', ny, ', k=', nz-1

fmt = '(4(es13.6,1x))'

do k = 1, nz - 1
    do j = 1, ny
        do i = 1, nx

            x(1) = pt_of_grid ( i, 1, 1 )  !--u-nodes
            x(2) = pt_of_grid ( j, 2, 1 )
            x(3) = pt_of_grid ( k, 3, 1 )

            write ( lun, fmt ) x(1), x(2), x(3), global_fmask( i, j, k )

        end do
    end do
end do

close ( lun )

end subroutine write_fmt_global_fmask

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
recursive subroutine calc_global_fmask_br ( r, br )
implicit none

!--perhaps define a type to hold all this, heck, even define a branch_type
real (rp), intent (in) :: r
type (branch_type), intent (in) :: br

character (*), parameter :: sub_name = mod_name // '.calc_global_fmask_br'

integer, parameter :: gen_max = 6  !--depth of recursion

integer :: i

real (rp) :: twist
real (rp) :: x_tmp(nd), y_tmp(nd)

type ( branch_type ) :: sbr


!---------------------------------------------------------------------

if ( br % gen <= gen_max ) call calc_global_fmask ( br )

if ( br % gen < gen_max ) then  !--recursion

    do i = 1, br % n_sub_branch

        !--prepare sub-branch sbr
        sbr = br

        sbr % gen = br % gen + 1

        sbr % l = r * (br % l)
        sbr % d = r * (br % d)
        sbr % abs_dir = tree_array(n_tree) % rel_dir(1, i) * br % x_hat +  &
                        tree_array(n_tree) % rel_dir(2, i) * br % y_hat +  &
                        tree_array(n_tree) % rel_dir(3, i) * br % z_hat

        sbr % z_hat = sbr % abs_dir
        sbr % x_hat = cross_product ( sbr % z_hat, br % abs_dir )
        sbr % y_hat = cross_product ( sbr % z_hat, sbr % x_hat )

        !--if this coordinate system is degenerate (sbr % [xy]_hat = 0)
        !  then just use parents local coordinate system
        if ( maxval (abs ( sbr % y_hat ) ) < epsilon ( 1.0_rp ) ) then

            sbr % x_hat = br % x_hat
            sbr % y_hat = br % y_hat
            sbr % z_hat = br % z_hat

        else  !--normalize

            sbr % x_hat = (sbr % x_hat) / mag ( sbr % x_hat )
            sbr % y_hat = (sbr % y_hat) / mag ( sbr % y_hat )
            sbr % z_hat = (sbr % z_hat) / mag ( sbr % z_hat )

        end if

        twist = tree_array(n_tree) % twist(i)  !--careful if trees different

        !--apply twist
        x_tmp = cos (twist) * (sbr % x_hat) + sin (twist) * (sbr % y_hat)
        y_tmp = cos (twist) * (sbr % y_hat) - sin (twist) * (sbr % x_hat)
        sbr % x_hat = x_tmp
        sbr % y_hat = y_tmp

        !--may need to check to see if other options apply
        if ( add_cap .and. sub_branches_outside ) then

            sbr % x0 = (br % x0) +  (0.5_rp * br % d) * (sbr % abs_dir) +  &
                       ( tree_array(n_tree) % root_height(i) *  &
                         ( br % l - 0.5_rp * br % d) * br % abs_dir )
                           !--cap is already added to l here
        else

          call error (sub_name, 'current x0_sub only valid for' // n_l //  &
                                'add_cap .and. sub_branches_outside')

        end if

        call calc_global_fmask_br ( r, sbr)

    end do

end if


end subroutine calc_global_fmask_br

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!--see calc_mask. unfortunately, had to re-invent wheel here
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine calc_global_fmask ( br )
implicit none

type (branch_type), intent (in) :: br

character (*), parameter :: sub_name = mod_name // '.calc_global_fmask'

logical, parameter :: use_loose_mask = .true.  !--expands mask slightly

integer :: i, j, k

real (rp) :: d_para
real (rp) :: eps
real (rp) :: x(nd)
real (rp) :: x_perp(nd)

!---------------------------------------------------------------------

eps = epsilon ( 1._rp )  !--this is the fudge factor in our mask

if ( use_loose_mask ) eps = eps + dx

!--unfortunately, loop over all points in the domain
do k = 1, nz - 1
    do j = 1, ny
        do i = 1, nx

            x(1) = pt_of_grid ( i, 1, 1 )  !--u-nodes
            x(2) = pt_of_grid ( j, 2, 1 )
            x(3) = pt_of_grid ( k, 3, 1 )

            x = x - br % x0

            d_para = dot_product (x, br % abs_dir)

            if ( ( 0 <= d_para ) .and. ( d_para <= br % l + eps ) ) then

                x_perp = x - d_para * (br % abs_dir)

                !--square cross-section test
                if ( maxval ( abs ( x_perp ) ) <= 0.5_rp * (br % d) + eps )  &
                    global_fmask( i, j, k ) = 1.0_rp

            end if

        end do
    end do
end do

end subroutine calc_global_fmask

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!--uses slow filtering algorithm, since only done once
!--this treats everything outside the domain as 0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine filter_global_fmask ()
implicit none

character (*), parameter :: sub_name = mod_name // '.filter_global_fmask'

integer :: i, j, k
integer :: ii, jj, kk
integer :: iimin, jjmin, kkmin
integer :: iimax, jjmax, kkmax

real (rp) :: total_in, total
real (rp) :: filtval
!real (rp) :: x(nd), xp(nd)
real (rp), allocatable :: wksp( :, :, : )

!---------------------------------------------------------------------

allocate ( wksp(ld, ny, nz) )

wksp = global_fmask

total_in = sum ( wksp )  !--to calculate normalization factor

do k = 1, nz - 1

  !x(3) = pt_of_grid ( k, 3, 1 )

  do j = 1, ny

    !x(2) = pt_of_grid ( j, 2, 1 )

    do i = 1, nx

      !x(1) = pt_of_grid ( i, 1, 1 )

      filtval = 0.0_rp

      !--perhaps a factor of 0.5 here?
      kkmin = max ( 1, k - fratiomax * idelta )
      kkmax = min ( nz - 1, k + fratiomax * idelta )

      jjmin = max ( 1, j - fratiomax * idelta )
      jjmax = min ( ny, j + fratiomax * idelta )

      iimin = max ( 1, i - fratiomax * idelta )
      iimax = min ( nx, i + fratiomax * idelta )

      do kk = kkmin, kkmax

        !xp(3) = pt_of_grid ( kk, 3, 1 )

        do jj = jjmin, jjmax

          !xp(2) = pt_of_grid ( jj, 2, 1 )

          do ii = iimin, iimax

            !xp(1) = pt_of_grid ( ii, 1, 1 )

            !--the normalization will be slightly off since discrete
            filtval = filtval +  &
                      truncgauss_kernel_3d ( i-ii, j-jj, k-kk ) *  &
                      wksp(ii, jj, kk)

          end do
        end do
      end do

      global_fmask( i, j, k ) = filtval

    end do
  end do
end do

where ( global_fmask < epsilon ( 1.0_rp ) ) global_fmask = 0.0_rp

!--watch accumulated error in doing these sums
!--no dV required since we are doing the
total = sum ( global_fmask )
global_fmask = global_fmask * total_in / total

deallocate ( wksp )


end subroutine filter_global_fmask

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!--requires external normalization
!--assumes filter width equal to grid spacing, and that grid spacing
!  equal in all directions
!--(|x|/delta)^2 = (i^2 + j^2 + k^2); dx and delta cancel
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
function truncgauss_kernel_3d ( i, j, k )
implicit none

real (rp) :: truncgauss_kernel_3d

integer, intent (in) :: i, j, k

integer, parameter :: m_max = ( fratiomax * idelta )**2

integer :: m

logical, save :: do_init = .true.

!real (rp) :: arg
real (rp), save :: table(0:m_max)

!---------------------------------------------------------------------

if ( do_init ) then
    do m = 0, m_max
        table(m) = exp ( -6.0_rp * m / idelta**2 )
    end do
    do_init = .false.
end if

m = i**2 + j**2 + k**2

if ( m <= m_max ) then
    truncgauss_kernel_3d = table(m)
else
    truncgauss_kernel_3d = 0.0_rp
end if

end function truncgauss_kernel_3d
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
end module trees_global_fmask_ls
