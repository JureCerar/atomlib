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
  use atomlib_crdio
  implicit none
  real, parameter :: delta = 0.001
  character(128) :: arg, message
  integer(INT64), allocatable :: offset(:)
  integer, allocatable :: atomi(:), resi(:)
  character(8), allocatable :: atomn(:), resn(:), resic(:), segid(:)
  real, allocatable :: wfact(:), coor(:,:)
  logical :: extended
  integer :: i, unit, nframes, natoms, status

  do i = 1, command_argument_count()
    call get_command_argument (i, arg, STATUS=status)
    if (status /= 0) error stop

    ! Read file
    call crd_open (arg, unit, nframes, offset)
    call crd_read_header (unit, extended, natoms)
    if (allocated(coor)) deallocate (atomi, atomn, resi, resn, segid, resic, wfact, coor)
    allocate (atomi(natoms), atomn(natoms), resi(natoms), resn(natoms), & 
    & segid(natoms), resic(natoms), wfact(natoms), coor(3,natoms))
    call crd_read_data (unit, extended, natoms, atomi, atomn, resi, resn, & 
    & segid, resic, wfact, coor)
    call crd_read_header (unit, extended, natoms)
    call crd_read_coor (unit, extended, natoms, coor)
    call crd_read_header (unit, extended, natoms)
    call crd_skip_data (unit, natoms)
    call crd_close (unit)

    ! Check if data was read correctly
    if (nframes /= 3) error stop
    if (natoms /= 44) error stop
    if (atomn(1) /= "N") error stop
    if (atomi(1) /= 1) error stop
    if (resn(1) /= "ALA") error stop
    if (resi(1) /= 1) error stop
    if (segid(1) /= "PROT") error stop
    if (resic(1) /= "A") error stop
    if (abs(wfact(1) - 1.000) > delta) error stop
    if (abs(coor(1,1) - 6.300) > delta) error stop
 
    ! Read file w/ error handling
    call crd_open (arg, unit, nframes, offset, STAT=status, ERRMSG=message)
    if (status /= 0) error stop message
    call crd_read_header (unit, extended, natoms, STAT=status, ERRMSG=message)
    if (allocated(coor)) deallocate (atomi, atomn, resi, resn, segid, resic, wfact, coor)
    allocate (atomi(natoms), atomn(natoms), resi(natoms), resn(natoms), & 
    & segid(natoms), resic(natoms), wfact(natoms), coor(3,natoms))
    call crd_read_data (unit, extended, natoms, atomi, atomn, resi, resn, & 
    & segid, resic, wfact, coor, STAT=status, ERRMSG=message)
    if (status /= 0) error stop message
    call crd_read_header (unit, extended, natoms, STAT=status, ERRMSG=message)
    if (status /= 0) error stop message
    call crd_read_coor (unit, extended, natoms, coor, STAT=status, ERRMSG=message)
    if (status /= 0) error stop message
    call crd_read_header (unit, extended, natoms, STAT=status, ERRMSG=message)
    if (status /= 0) error stop message
    call crd_skip_data (unit, natoms, STAT=status, ERRMSG=message)
    if (status /= 0) error stop message
    call crd_close (unit, STAT=status, ERRMSG=message)
    if (status /= 0) error stop message
    
    ! Write file
    open (NEWUNIT=unit, FILE="out.crd", ACTION="write", STATUS="unknown")
    call crd_write (unit, natoms, atomi, atomn, resi, resn, segid, resic, wfact, coor)
    close (unit)
    
    ! Write file w/ error handling 
    open (NEWUNIT=unit, FILE="out.crd", ACTION="write", STATUS="unknown", IOSTAT=status, IOMSG=message)
    if (status /= 0) error stop message
    call crd_write (unit, natoms, atomi, atomn, resi, resn, segid, resic, wfact, coor, STAT=status, ERRMSG=message)
    if (status /= 0) error stop message
    close (unit, IOSTAT=status, IOMSG=message)
    if (status /= 0) error stop message

  end do

end program main