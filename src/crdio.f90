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

module atomlib_crdio
  implicit none
  private
  public :: crd_open, crd_close
  public :: crd_read_header, crd_read_data, crd_read_coor, crd_skip_data, crd_write

  ! TODO:
  ! - [ ] Add read for binary files (I can't find any format documentation).
  ! - [ ] Add comments to read header. 
 
  ! Global constants
  integer, parameter :: DIM = 3 ! Default dimension
  integer, parameter :: ERRLEN = 128 ! Length of error string
  integer, parameter :: CRD_MAX_NUM = 100000 ! Max. values of residue and atom index
  integer, parameter :: CRD_MAX_LEN = 140 ! Max. length of comment string

  ! Structure of CRD/COR file for "normal" format i.e less than 100000 atoms and with 
  ! residue and with residue or atom names less than five characters long:
  !
  ! COMMENTS (multiple lines starting with "*")
  ! NATOMS (I5)
  ! ATOMI RESI   RESN   ATOMN  X  Y  Z   SEGID   RESID Weighting
  !    I5   I5 X   A4 X    A4   3F10.5 X    A4 X    A4     F10.5
  ! ...
  ! 
  ! And for expanded format for more than 100000 atoms (upto 10**10) and with residue or
  ! atom names upto 8 character long (supported by CHARMM version C31A1 and later):
  !
  ! COMMENTS (multiple lines starting with "*")
  ! NATOMS (I10 followed by 'EXT')
  ! ATOMI RESI    RESN    ATOMN  X  Y  Z    SEGID    RESID   Weighting
  !   I10  I10 2X   A8 2X    A8  3F20.10 2X    A8 2X    A8      F20.10
  ! ...
  !
  ! See: https://www.ks.uiuc.edu/Research/vmd/plugins/molfile/corplugin.html
  !      http://thegrantlab.org/bio3d/reference/write.crd.html
  !      https://www.ks.uiuc.edu/Research/vmd/plugins/doxygen/corplugin_8c-source.html

contains

! Check CRD/COR file and construct frame offset table.
subroutine crd_open (file, unit, nframes, offset, stat, errmsg)
  use iso_fortran_env, only: INT64, IOSTAT_END
  implicit none
  character(*), intent(in) :: file ! File name.
  integer, intent(out) :: unit ! Fortran file handle.
  integer, intent(out) :: nframes ! Num. of frames in file
  integer(INT64), allocatable, intent(out) :: offset(:) ! Stream position offsets.
  integer, intent(out), optional :: stat ! Error status code. Returns zero if no error.
  character(*), intent(out), optional :: errmsg ! Error message.
  character(ERRLEN) :: message
  integer(INT64) :: pos
  integer(INT64), allocatable :: temp(:)
  integer :: status, natoms
  logical :: extended

  catch: block

    ! Open file as stream.
    open (FILE=file, NEWUNIT=unit,  ACCESS="stream", FORM="formatted", &
    &     STATUS="old", ACTION="read", IOSTAT=status, IOMSG=message)
    if (status /= 0) exit catch

    ! Start at offset table of size 128. If needed, we will expand it as we go along.
    if (allocated(offset)) deallocate (offset, STAT=status, ERRMSG=message)
    allocate (offset(128), SOURCE=0_INT64, STAT=status, ERRMSG=message)
    if (status /= 0) exit catch

    nframes = 0

    while: do while (.true.)
      ! Store current position
      inquire (UNIT=unit, POS=pos, IOSTAT=status, IOMSG=message)
      if (status /= 0) exit catch

      ! Read header and if EOF then there are (probably) no more frames in trajectory.
      call crd_read_header (unit, extended, natoms, STAT=status, ERRMSG=message)
      if (status == IOSTAT_END) exit while
      if (status /= 0) exit catch

      ! Expand offset table if needed
      if (nframes > size(offset)) then
        allocate (temp(2 * size(offset)), STAT=status, ERRMSG=message)
        if (status /= 0) exit catch
        temp(1:nframes) = offset(1:nframes)
        call move_alloc (temp, offset)
      end if
      
      nframes = nframes + 1
      offset(nframes) = pos

      ! Skip remaining frame and only check for errors.
      call crd_skip_data (unit, natoms, STAT=status, ERRMSG=message)
      if (status /= 0) exit catch

    end do while

    ! Rewind to start of the file for actual processing.
    rewind (unit, IOSTAT=status, IOMSG=message)
    if (status /= 0) exit catch

    ! Reallocate offset table to actual size.
    allocate (temp(nframes), STAT=status, ERRMSG=message)
    if (status /= 0) exit catch
    temp = offset(1:nframes)
    call move_alloc (temp, offset)

  end block catch

  ! Error handling
  if (present(stat)) then
    stat = status
  else if (status /= 0) then
    error stop message
  end if
  if (present(errmsg)) errmsg = trim(message)

end subroutine crd_open

! Close CRD/COR file.
subroutine crd_close (unit, stat, errmsg)
  implicit none
  integer, intent(in) :: unit ! Fortran file handle.
  integer, intent(out), optional :: stat ! Error status code. Returns zero if no error.
  character(*), intent(out), optional :: errmsg ! Error message.
  character(ERRLEN) :: message
  integer :: status

  close (UNIT=unit, IOSTAT=status, IOMSG=message)
  if (present(stat)) then
    stat = status
  else if (status /= 0) then
    error stop message
  end if
  if (present(errmsg)) errmsg = trim(message) 

end subroutine crd_close

! Open and check CRD/COR file. Construct frame offset table. 
subroutine crd_read_header (unit, extended, natoms, stat, errmsg)
  implicit none
  integer, intent(in) :: unit ! Fortran file handle.
  logical, intent(out) :: extended ! Is file in normal or extended format.
  integer, intent(out) :: natoms ! Number of atoms
  integer, intent(out), optional :: stat ! Error status code. Returns zero if no error.
  character(*), intent(out), optional :: errmsg ! Error message.
  character(ERRLEN) :: message
  character(CRD_MAX_LEN+2) :: buffer
  character(9) :: action, access
  logical :: opened
  integer :: status

  catch: block

    ! Check if file unit is opened for reading and NOT in binary.
    inquire (UNIT=unit, OPENED=opened, ACTION=action, ACCESS=access, IOSTAT=status, IOMSG=message)
    if (.not. opened .or. index(action, "READ") == 0) then
      message = "File unit not opened for reading"
      status = 1
    else if (index(access, "DIRECT") /= 0) then
      message = "Binary CRD/COR files are not supported"
      status = 1
    end if
    if (status /= 0) exit catch

    ! Skip comments in header. Comments always start with an asterisk ("*").
    do while (.true.)
      read (unit, "(a)", IOSTAT=status, IOMSG=message) buffer
      if (status /= 0)  exit catch

      buffer = trim(adjustl(buffer))
      if (index(buffer, "*") /= 1) exit ! while

    end do

    ! Ups, we read to far. No problem; last read line should contain number of atoms.
    read (buffer, *, IOSTAT=status, IOMSG=message) natoms
    if (status /= 0) exit catch

    ! Check if number of atoms makes sense
    if (natoms < 0) then
      message = "Invalid number of atoms"
      status = 1
      exit catch
    end if

    ! Check if file is in extended IO format
    extended = .false.
    extended = extended .or. (index(buffer, "EXT") /= 0)
    extended = extended .or. (natoms >= CRD_MAX_NUM)

  end block catch

  if (present(stat)) then
    stat = status
  else if (status /= 0) then
    error stop message
  end if
  if (present(errmsg)) errmsg = trim(message)

end subroutine crd_read_header

! Read data in CRD/COR format.
subroutine crd_read_data (unit, extended, natoms, atomi, atomn, resi, resn, & 
& segid, resic, wfact, coor, stat, errmsg)
  implicit none
  integer, intent(in) :: unit ! Fortran file handle.
  logical, intent(in) :: extended ! Is file in normal or extended format.
  integer, intent(in) :: natoms ! Number of atoms.
  integer, intent(out) :: atomi(natoms) ! Atom index
  character(*), intent(out) :: atomn(natoms) ! Atom name
  integer, intent(out) :: resi(natoms) ! Residue index
  character(*), intent(out) :: resn(natoms) ! Residue name
  character(*), intent(out) :: segid(natoms) ! Segment identifier
  character(*), intent(out) :: resic(natoms) ! Alternative residue identifier.
  real, intent(out) :: wfact(natoms) ! Weighting array value
  real, intent(out) :: coor(DIM,natoms) ! Particle's cartesian coordinates in [A] units.
  integer, intent(out), optional :: stat ! Error status code. Returns zero if no error.
  character(*), intent(out), optional :: errmsg ! Error message.
  character(ERRLEN) :: message
  character(64) :: fmt
  character(9) :: action, access
  logical :: opened
  integer :: i, status

  catch: block

    ! Check if file unit is opened for reading and NOT in binary.
    inquire (UNIT=unit, OPENED=opened, ACTION=action, ACCESS=access, IOSTAT=status, IOMSG=message)
    if (.not. opened .or. index(action, "READ") == 0) then
      message = "File unit not opened for reading"
      status = 1
    else if (index(access, "DIRECT") /= 0) then
      message = "Binary CRD/COR files are not supported"
      status = 1
    end if
    if (status /= 0) exit catch

    ! Select format based if normal or extended format.
    if (extended) then
      fmt = "(i10,i10,2x,a8,2x,a8,3(f20.10),2x,a8,2x,a8,f20.10)"
    else
      fmt = "(i5,i5,x,a4,x,a4,3(f10.5),x,a4,x,a4,f10.5)"
    end if

    do i = 1, natoms
      read (unit, fmt, IOSTAT=status, IOMSG=message) atomi(i), resi(i), resn(i), atomn(i), &
      & coor(:,i), segid(i), resic(i), wfact(i)
      if (status /= 0) exit catch

    end do

  end block catch

  if (present(stat)) then
    stat = status
  else if (status /= 0) then
    error stop message
  end if
  if (present(errmsg)) errmsg = trim(message)

end subroutine crd_read_data

! Read data in CRD/COR format and only store coordinates and box.
subroutine crd_read_coor (unit, extended, natoms, coor, stat, errmsg)
  implicit none
  integer, intent(in) :: unit ! Fortran file handle.
  logical, intent(in) :: extended ! Is file in normal or extended format.
  integer, intent(in) :: natoms ! Number of atoms
  real, intent(out) :: coor(DIM,natoms) ! Particle's cartesian coordinates in [A] units.
  integer, intent(out), optional :: stat ! Error status code. Returns zero if no error.
  character(*), intent(out), optional :: errmsg ! Error message.
  character(ERRLEN) :: message
  character(64) :: fmt
  character(9) :: action, access
  logical :: opened
  integer :: i, status

  catch: block

    ! Check if file unit is opened for reading and NOT in binary.
    inquire (UNIT=unit, OPENED=opened, ACTION=action, ACCESS=access, IOSTAT=status, IOMSG=message)
    if (.not. opened .or. index(action, "READ") == 0) then
      message = "File unit not opened for reading"
      status = 1
    else if (index(access, "DIRECT") /= 0) then
      message = "Binary CRD files are not supported"
      status = 1
    end if
    if (status /= 0) exit catch

    ! Select format based if normal or extended format.
    if (extended) then
      fmt = "(40x,3(f20.10))"
    else
      fmt = "(20x,3(f10.5))"
    end if

    do i = 1, natoms
      read (unit, fmt, IOSTAT=status, IOMSG=message) coor(:,i)
      if (status /= 0) exit catch
    end do

  end block catch

  if (present(stat)) then
    stat = status
  else if (status /= 0) then
    error stop message
  end if
  if (present(errmsg)) errmsg = trim(message)

end subroutine crd_read_coor

! Read through CRD/COR file but do not store any data.
subroutine crd_skip_data (unit, natoms, stat, errmsg)
  implicit none
  integer, intent(in) :: unit ! Fortran file handle.
  integer, intent(in) :: natoms ! Number of atoms
  integer, intent(out), optional :: stat ! Error status code. Returns zero if no error.
  character(*), intent(out), optional :: errmsg ! Error message.
  character(ERRLEN) :: message
  character(9) :: action, access
  logical :: opened
  integer :: i, status

  catch: block

    ! Check if file unit is opened for reading and NOT in binary.
    inquire (UNIT=unit, OPENED=opened, ACTION=action, ACCESS=access, IOSTAT=status, IOMSG=message)
    if (.not. opened .or. index(action, "READ") == 0) then
      message = "File unit not opened for reading"
      status = 1
    else if (index(access, "DIRECT") /= 0) then
      message = "Binary CRD files are not supported"
      status = 1
    end if
    if (status /= 0) exit catch

    do i = 1, natoms
      read (unit, *, IOSTAT=status, IOMSG=message)
      if (status /= 0) exit catch
    end do

  end block catch

  if (present(stat)) then
    stat = status
  else if (status /= 0) then
    error stop message
  end if
  if (present(errmsg)) errmsg = trim(message)

end subroutine crd_skip_data

! Write data in CRD/COR format. Will automaticaly write in expanded format when needed.
subroutine crd_write (unit, natoms, atomi, atomn, resi, resn, segid, resic, wfact, coor, stat, errmsg)
  implicit none
  integer, intent(in) :: unit ! Fortran file handle.
  integer, intent(in) :: natoms ! Number of atoms 
  integer, intent(in) :: atomi(natoms) ! Atom index
  character(*), intent(in) :: atomn(natoms) ! Atom name
  integer, intent(in) :: resi(natoms) ! Residue index
  character(*), intent(in) :: resn(natoms) ! Residue name
  character(*), intent(in) :: segid(natoms) ! Segment identifier
  character(*), intent(in) :: resic(natoms) ! Residue identifier
  real, intent(in) :: wfact(natoms) ! Weighting array value
  real, intent(in) :: coor(DIM,natoms) ! Particle coordinates
  integer, intent(out), optional :: stat ! Error status code. Returns zero if no error.
  character(*), intent(out), optional :: errmsg ! Error message.
  character(ERRLEN) :: message
  character(9) :: action, access 
  character(64) :: fmt
  logical :: opened, extended
  integer :: i, status

  catch: block

    ! Check if file unit is opened for writing and NOT in binary.
    inquire (UNIT=unit, OPENED=opened, ACTION=action, ACCESS=access, IOSTAT=status, IOMSG=message)
    if (.not. opened .or. index(action, "WRITE") == 0) then
      message = "File unit not opened for writing"
      status = 1
    else if (index(access, "DIRECT") /= 0) then
      message = "Binary CRD/COR files are not supported"
      status = 1
    end if
    if (status /= 0) exit catch

    ! Determine if normal or extended format will be used.
    extended = .false.
    extended = extended .or. (natoms >= CRD_MAX_NUM)
    extended = extended .or. any(len_trim(atomn) > 5)
    extended = extended .or. any(len_trim(resn) > 5)

    ! Add empty comment
    write (unit, "('*',x,a140)", IOSTAT=status, IOMSG=message) ""
    if (status /= 0) exit catch

    if (extended) then
      write (unit, "(i10,2x,'EXT')", IOSTAT=status, IOMSG=message) natoms
    else 
      write (unit, "(i5)", IOSTAT=status, IOMSG=message) natoms
    end if
    if (status /= 0) exit catch

    if (extended) then
      fmt = "(i10,i10,2x,a8,2x,a8,3(f20.10),2x,a8,2x,a8,f20.10)"
    else
      fmt = "(i5,i5,x,a4,x,a4,3(f10.5),x,a4,x,a4,f10.5)"
    end if

    ! All strings are right adjusted.
    do i = 1, natoms
      write (unit, fmt, IOSTAT=status, IOMSG=message) atomi(i), resi(i), adjustl(resn(i)), &
      & adjustl(atomn(i)), coor(:,i), adjustl(segid(i)), adjustl(resic(i)), wfact(i)
      if (status /= 0) exit catch
    end do
 
  end block catch

  if (present(stat)) then
    stat = status
  else if (status /= 0) then
    error stop message
  end if
  if (present(errmsg)) errmsg = trim(message)

end subroutine crd_write

end module atomlib_crdio