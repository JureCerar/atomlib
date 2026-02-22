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
  use atomlib_pdbio
  implicit none
  real, parameter :: delta = 0.001
  character(128) :: arg, message
  integer(INT64), allocatable :: offset(:)
  integer, allocatable :: atomi(:), resi(:), charge(:)
  character(8), allocatable :: atomn(:), resn(:), chain(:), resic(:), altloc(:), element(:)
  real, allocatable :: qfact(:), bfact(:), coor(:,:)
  real :: box(3,3)
  integer :: i, unit, nframes, natoms, status

  do i = 1, command_argument_count()
    call get_command_argument (i, arg, STATUS=status)
    if (status /= 0) error stop

    ! Read file
    call pdb_open (arg, unit, nframes, offset)
    call pdb_get_natoms (unit, natoms)
    if (allocated(coor)) deallocate (atomi, atomn, resi, resn, chain, resic, &
    & altloc, qfact, bfact, element, charge, coor)
    allocate (atomi(natoms), atomn(natoms), resi(natoms), resn(natoms), chain(natoms), resic(natoms), &
    & altloc(natoms), qfact(natoms), bfact(natoms), element(natoms), charge(natoms), coor(3,natoms))
    call pdb_read_data (unit, natoms, atomi, atomn, resi, resn, chain, resic, &
    & altloc, qfact, bfact, element, charge, coor, box)
    call pdb_get_natoms (unit, natoms)
    call pdb_read_coor (unit, natoms, coor, box)
    call pdb_get_natoms (unit, natoms)
    call pdb_skip_data (unit, natoms)
    call pdb_close (unit)

    ! Check if data was read correctly
    if (nframes /= 3) error stop
    if (natoms /= 44) error stop
    if (atomn(1) /= "N") error stop
    if (atomi(1) /= 1) error stop
    if (resn(1) /= "ALA") error stop
    if (resi(1) /= 1) error stop
    if (chain(1) /= "A") error stop
    if (resic(1) /= "") error stop
    if (altloc(1) /= "") error stop
    if (abs(qfact(1) - 1.000) > delta) error stop
    if (abs(bfact(1) - 0.000) > delta) error stop
    if (element(1) /= "N") error stop
    if (charge(1) /= -1) error stop
    if (abs(coor(1,1) - 6.304) > delta) error stop
    if (abs(box(1,1) - 20.000) > delta) error stop
    
    ! Read file w/ error handling
    call pdb_open (arg, unit, nframes, offset, STAT=status, ERRMSG=message)
    if (status /= 0) error stop message
    call pdb_get_natoms (unit, natoms, STAT=status, ERRMSG=message)
    if (status /= 0) error stop message
    if (allocated(coor)) deallocate (atomi, atomn, resi, resn, chain, resic, &
    & altloc, qfact, bfact, element, charge, coor, STAT=status, ERRMSG=message)
    allocate (atomi(natoms), atomn(natoms), resi(natoms), resn(natoms), chain(natoms), resic(natoms), &
    & altloc(natoms), qfact(natoms), bfact(natoms), element(natoms), charge(natoms), coor(3,natoms))
    call pdb_read_data (unit, natoms, atomi, atomn, resi, resn, chain, resic, &
    & altloc, qfact, bfact, element, charge, coor, box, STAT=status, ERRMSG=message)
    if (status /= 0) error stop message
    call pdb_get_natoms (unit, natoms, STAT=status, ERRMSG=message)
    if (status /= 0) error stop message
    call pdb_read_coor (unit, natoms, coor, box, STAT=status, ERRMSG=message)
    if (status /= 0) error stop message
    call pdb_get_natoms (unit, natoms, STAT=status, ERRMSG=message)
    if (status /= 0) error stop message
    call pdb_skip_data (unit, natoms, STAT=status, ERRMSG=message)
    if (status /= 0) error stop message
    call pdb_close (unit, STAT=status, ERRMSG=message)
    if (status /= 0) error stop message
    
    ! Write file
    open (NEWUNIT=unit, FILE="out.pdb", ACTION="write", STATUS="unknown")
    call pdb_write (unit, natoms, atomi, atomn, resi, resn, chain, resic, &
    & altloc, qfact, bfact, element, charge, coor, box)
    close (unit)
    
    ! Write file w/ error handling 
    open (NEWUNIT=unit, FILE="out.pdb", ACTION="write", STATUS="unknown", IOSTAT=status, IOMSG=message)
    if (status /= 0) error stop message
    call pdb_write (unit, natoms, atomi, atomn, resi, resn, chain, resic, &
    & altloc, qfact, bfact, element, charge, coor, box, STAT=status, ERRMSG=message)
    if (status /= 0) error stop message
    close (unit, IOSTAT=status, IOMSG=message) 
    if (status /= 0) error stop message

  end do

end program main