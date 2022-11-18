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
  use atomlib_xtcio
  implicit none
  real, parameter :: delta = 0.001
  type(xdrfile), pointer :: xd => null()
  character(128) :: arg, message
  integer(INT64), allocatable :: offset(:)
  real, allocatable :: coor(:,:)
  real :: time, prec, box(3,3)
  integer :: i, nframes, natoms, step, status

  do i = 1, command_argument_count()
    call get_command_argument (i, arg, STATUS=status)
    if (status /= 0) error stop

    ! Read file
    call xtc_open (arg, xd, nframes, offset)
    call xtc_do_header (xd, natoms, step, time)
    call xtc_do_skip (xd, natoms)
    call xtc_do_header (xd, natoms, step, time)
    if (allocated(coor)) deallocate (coor)
    allocate (coor(3,natoms))
    call xtc_do_data (xd, .True., natoms, prec, coor, box)
    call xtc_close (xd)

    ! Check if data was read correctly
    if (nframes /= 3) error stop
    if (natoms /= 44) error stop
    if (step /= 2) error stop
    if (abs(prec - 1E4) > delta) error stop
    if (abs(time - 2.000) > delta) error stop
    if (abs(coor(1,1) - 0.6300) > delta) error stop
    if (abs(box(1,1) - 2.000) > delta) error stop
 
    ! Read file w/ error handling
    call xtc_open (arg, xd, nframes, offset, STAT=status, ERRMSG=message)
    if (status /= 0) error stop message
    call xtc_do_header (xd, natoms, step, time, STAT=status, ERRMSG=message)
    if (status /= 0) error stop message
    call xtc_do_skip (xd, natoms, STAT=status, ERRMSG=message)
    if (status /= 0) error stop message
    call xtc_do_header (xd, natoms, step, time, STAT=status, ERRMSG=message)
    if (status /= 0) error stop message
    if (allocated(coor)) deallocate (coor)
    allocate (coor(3,natoms))
    call xtc_do_data (xd, .True., natoms, prec, coor, box, STAT=status, ERRMSG=message)
    if (status /= 0) error stop message
    call xtc_close (xd, STAT=status, ERRMSG=message)
    if (status /= 0) error stop message
    
    
    ! Write file
    call xdr_open (xd, "out.xtc", .False.)
    call xtc_do_header (xd, natoms, step, time)
    call xtc_do_data (xd, .False., natoms, prec, coor, box)
    call xdr_close (xd)
    
    ! Write file w/ error handling 
    call xdr_open (xd, "out.xtc", .False., STAT=status, ERRMSG=message)
    if (status /= 0) error stop message
    call xtc_do_header (xd, natoms, step, time, STAT=status, ERRMSG=message)
    if (status /= 0) error stop message
    call xtc_do_data (xd, .False., natoms, prec, coor, box, STAT=status, ERRMSG=message)
    if (status /= 0) error stop message
    call xdr_close (xd, STAT=status, ERRMSG=message)
    if (status /= 0) error stop message

  end do

end program main