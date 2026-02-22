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
  use atomlib_dcdio
  implicit none
  real, parameter :: delta = 0.001
  character(128) :: arg, message
  integer(INT64), allocatable :: offset(:)
  character(128), allocatable :: remark(:)
  real, allocatable :: coor(:,:)
  real :: timestep, box(3,3)
  integer :: start_time, every_time, end_time
  integer :: i, natoms, nframes, unit, status

  do i = 1, command_argument_count()
    call get_command_argument (i, arg, STATUS=status)
    if (status /= 0) error stop

    ! Read file
    call dcd_open (arg, unit, nframes, offset)
    call dcd_read_header (unit, nframes, remark, start_time, every_time, end_time, timestep, natoms)
    call dcd_skip_data (unit, natoms)
    if (allocated(coor)) deallocate (coor)
    allocate (coor(3,natoms))
    call dcd_read_data (unit, natoms, coor, box)
    call dcd_close (unit)

    ! Check if data was read correctly
    if (nframes /= 3) error stop
    if (natoms /= 44) error stop
    if (start_time /= 0) error stop
    if (every_time /= 1) error stop
    if (end_time /= 3) error stop
    if (abs(timestep - 1.000) > delta) error stop
    if (abs(coor(1,1) - 6.300) > delta) error stop
    if (abs(box(1,1) - 20.00) > delta) error stop
 
    ! Read file w/ error handling
    call dcd_open (arg, unit, nframes, offset, STAT=status, ERRMSG=message)
    if (status /= 0) error stop message
    call dcd_read_header (unit, nframes, remark, start_time, every_time, end_time, &
    & timestep, natoms, STAT=status, ERRMSG=message)
    if (status /= 0) error stop message
    call dcd_skip_data (unit, natoms, STAT=status, ERRMSG=message)
    if (status /= 0) error stop message
    if (allocated(coor)) deallocate (coor, STAT=status, ERRMSG=message)
    allocate (coor(3,natoms), STAT=status, ERRMSG=message)
    if (status /= 0) error stop message
    call dcd_read_data (unit, natoms, coor, box, STAT=status, ERRMSG=message)
    if (status /= 0) error stop message
    call dcd_close (unit, STAT=status, ERRMSG=message)
    if (status /= 0) error stop message
    
    ! Write file
    open (NEWUNIT=unit, FILE="out.dcd", FORM="unformatted", ACCESS="stream", ACTION="readwrite", STATUS="unknown")
    call dcd_write_header (unit, remark, start_time, every_time, end_time, timestep, natoms)
    call dcd_write_data (unit, natoms, coor, box)
    close (unit)
    
    ! Write file w/ error handling 
    open (NEWUNIT=unit, FILE="out.dcd", FORM="unformatted", ACCESS="stream", & 
    & ACTION="readwrite", STATUS="unknown", IOSTAT=status, IOMSG=message)
    if (status /= 0) error stop message
    call dcd_write_header (unit, remark, start_time, every_time, end_time, & 
    & timestep, natoms, STAT=status, ERRMSG=message)
    if (status /= 0) error stop message
    call dcd_write_data (unit, natoms, coor, box, STAT=status, ERRMSG=message)
    if (status /= 0) error stop message
    close (unit, IOSTAT=status, IOMSG=message)
    if (status /= 0) error stop message

  end do

end program main