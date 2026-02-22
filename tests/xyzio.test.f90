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
  use iso_fortran_env, only: INT32, INT64, REAL32, REAL64
  use atomlib_xyzio
  implicit none
  real, parameter :: delta = 0.001
  character(128) :: arg, comment, message
  integer(INT64), allocatable :: offset(:)
  real, allocatable :: coor(:,:)
  character(8), allocatable :: name(:)
  integer :: i, unit, nframes, natoms, status

  do i = 1, command_argument_count()
    call get_command_argument (i, arg, STATUS=status)
    if (status /= 0) error stop

    ! Read file
    call xyz_open (arg, unit, nframes, offset)
    call xyz_read_header (unit, comment, natoms)
    if (allocated(coor)) deallocate (coor, name)
    allocate (coor(3,natoms), name(natoms))
    call xyz_read_data (unit, natoms, name, coor)
    call xyz_read_header (unit, comment, natoms)
    call xyz_read_coor (unit, natoms, coor)
    call xyz_read_header (unit, comment, natoms)
    call xyz_skip_data (unit, natoms)
    call xyz_close (unit)

    ! Check if data was read correctly
    if (nframes /= 3) error stop
    if (natoms /= 44) error stop
    if (name(1) /= "N") error stop
    if (abs(coor(1,1) - 6.304) > delta) error stop
 
    ! Read file w/ error handling
    call xyz_open (arg, unit, nframes, offset, STAT=status, ERRMSG=message)
    if (status /= 0) error stop message
    call xyz_read_header (unit, comment, natoms, STAT=status, ERRMSG=message)
    if (status /= 0) error stop message
    if (allocated(coor)) deallocate (coor, name)
    allocate (coor(3,natoms), name(natoms))
    call xyz_read_data (unit, natoms, name, coor, STAT=status, ERRMSG=message)
    if (status /= 0) error stop message
    call xyz_read_header (unit, comment, natoms, STAT=status, ERRMSG=message)
    if (status /= 0) error stop message
    call xyz_read_coor (unit, natoms, coor, STAT=status, ERRMSG=message)
    if (status /= 0) error stop message
    call xyz_read_header (unit, comment, natoms, STAT=status, ERRMSG=message)
    if (status /= 0) error stop message
    call xyz_skip_data (unit, natoms, STAT=status, ERRMSG=message)
    if (status /= 0) error stop message
    call xyz_close (unit, STAT=status, ERRMSG=message)
    if (status /= 0) error stop message
    
    ! Write file
    open (NEWUNIT=unit, FILE="out.xyz", ACTION="write", STATUS="unknown")
    call xyz_write (unit, comment, natoms, name, coor)
    close (unit)
    
    ! Write file w/ error handling 
    open (NEWUNIT=unit, FILE="out.xyz", ACTION="write", STATUS="unknown", IOSTAT=status, IOMSG=message)
    if (status /= 0) error stop message
    call xyz_write (unit, comment, natoms, name, coor, STAT=status, ERRMSG=message)
    if (status /= 0) error stop message
    close (unit, IOSTAT=status, IOMSG=message)
    if (status /= 0) error stop message

  end do

end program main