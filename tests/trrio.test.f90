! This file is part of atomlib
! https://github.com/JureCerar/atomlib
!
! Copyright (C) 2022 Jure Cerar
!
! This program is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with this program.  If not, see <https://www.gnu.org/licenses/>.

program main
  use iso_fortran_env, only: INT64, REAL64
  use xdrfor
  use atomlib_trrio
  implicit none
  type(xdrfile), pointer :: xd => null()
  type(trnheader) :: sh
  real, parameter :: delta = 0.001
  character(128) :: arg, message
  integer(INT64), allocatable :: offset(:)
  real, allocatable :: coor(:,:)
  real :: box(3,3)
  integer :: i, nframes, status

  do i = 1, command_argument_count()
    call get_command_argument (i, arg, STATUS=status)
    if (status /= 0) error stop

    ! Read file
    call trr_open (arg, xd, nframes, offset)
    call trr_do_header (xd, .True., sh)
    call trr_do_skip (xd, sh)
    call trr_do_header (xd, .True., sh)
    if (allocated(coor)) deallocate (coor)
    allocate (coor(3,sh%natoms))
    call trr_do_data (xd, sh, coor, box)
    call trr_close (xd)

    ! Check if data was read correctly
    if (nframes /= 3) error stop
    if (sh%natoms /= 44) error stop
    if (sh%step /= 2) error stop
    if (abs(sh%time - 2.000) > delta) error stop
    if (abs(sh%lambda - 0.000) > delta) error stop
    if (abs(coor(1,1) - 0.6300) > delta) error stop
 
    ! Read file w/ error handling
    call trr_open (arg, xd, nframes, offset, STAT=status, ERRMSG=message)
    if (status /= 0) error stop message
    call trr_do_header (xd, .True., sh, STAT=status, ERRMSG=message)
    if (status /= 0) error stop message
    call trr_do_skip (xd, sh, STAT=status, ERRMSG=message)
    if (status /= 0) error stop message
    call trr_do_header (xd, .True., sh, STAT=status, ERRMSG=message)
    if (status /= 0) error stop message
    if (allocated(coor)) deallocate (coor)
    allocate (coor(3,sh%natoms))
    call trr_do_data (xd, sh, coor, box, STAT=status, ERRMSG=message)
    if (status /= 0) error stop message
    call trr_close (xd, STAT=status, ERRMSG=message)
    if (status /= 0) error stop message
    
    ! Write file
    call xdr_open (xd, "out.trr", .False.)
    call trr_do_header (xd, .False., sh)
    call trr_do_data (xd, sh, coor, box)
    call xdr_close (xd)
    
    ! Write file w/ error handling 
    call xdr_open (xd, "out.trr", .False., STAT=status, ERRMSG=message)
    if (status /= 0) error stop message
    call trr_do_header (xd, .False., sh, STAT=status, ERRMSG=message)
    if (status /= 0) error stop message
    call trr_do_data (xd, sh, coor, box, STAT=status, ERRMSG=message)
    if (status /= 0) error stop message
    call xdr_close (xd, STAT=status, ERRMSG=message)
    if (status /= 0) error stop message

  end do

end program main