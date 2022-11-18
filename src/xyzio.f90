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

module atomlib_xyzio
  implicit none
  private
  public :: xyz_open, xyz_close, xyz_read_header, xyz_read_data, xyz_read_coor, xyz_skip_data, xyz_write

  ! Global constants
  integer, parameter :: DIM = 3 ! Default dimension
  integer, parameter :: ERRLEN = 256 ! Length of error message

  ! Structure of XYZ file:
  !
  ! NATOMS
  ! COMMENT
  ! NAME    X  Y  Z (name is optional, format is not defined)
  ! ...
  ! 
  ! See: <reference>

contains

! Util: Count number of tokens in a string seperated by delimiter.
integer function cnttok (string, del)
  implicit none
  character(*), intent(in) :: string
  character(*), intent(in) :: del
  integer :: i

  cnttok = 0
  i = len_trim(string)
  do while (i > 0)
    cnttok = cnttok + 1
    i = index(trim(string(:i)), del, BACK=.true.)
  end do

end function cnttok

! Open and check XYZ file. Construct frame offset table. 
subroutine xyz_open (file, unit, nframes, offset, stat, errmsg)
  use iso_fortran_env, only: INT64, IOSTAT_END
  implicit none
  character(*), intent(in) :: file ! Name of file.
  integer, intent(out) :: unit ! Fortran file handle.
  integer, intent(out) :: nframes ! Num. of frames in file.
  integer(INT64), allocatable, intent(out) :: offset(:) ! Stream position offsets.
  integer, intent(out), optional :: stat ! Error status code. Returns zero if no error.
  character(*), intent(out), optional :: errmsg ! Error message.
  character(ERRLEN) :: message, comment
  integer(INT64) :: pos
  integer(INT64), allocatable :: temp(:)
  integer :: status, natoms

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
      ! Store current position.
      inquire (UNIT=unit, POS=pos, IOSTAT=status, IOMSG=message)
      if (status /= 0) exit catch

      ! Read header and if EOF then there are (probably) no more frames in trajectory.
      call xyz_read_header (unit, comment, natoms, STAT=status, ERRMSG=message)
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
      call xyz_skip_data (unit, natoms, STAT=status, ERRMSG=message)
      if (status /= 0) exit catch

    end do while    

    ! Rewind to start of the file for actual processing.
    rewind (unit, IOSTAT=status, IOMSG=message)
    if (status /= 0) exit catch

    ! Allocate offset table to actual size.
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

end subroutine xyz_open

! Close XYZ file.
subroutine xyz_close (unit, stat, errmsg)
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

end subroutine xyz_close

! Read header in XYZ format.
subroutine xyz_read_header (unit, comment, natoms, stat, errmsg)
  implicit none
  integer, intent(in) :: unit ! Fortran file handle.
  character(*), intent(out) :: comment ! Comment line
  integer, intent(out) :: natoms ! Number of atoms
  integer, intent(out), optional :: stat ! Error status code. Returns zero if no error.
  character(*), intent(out), optional :: errmsg ! Error message.
  character(ERRLEN) :: message
  character(16) :: action
  logical :: opened
  integer :: status

  catch: block

    ! Check if file unit is opened for reading
    inquire (UNIT=unit, OPENED=opened, ACTION=action, IOSTAT=status, IOMSG=message)
    if (.not. opened .or. index(action, "READ") == 0) then
      message = "File unit not opened for reading"
      status = 1
    end if
    if (status /= 0) exit catch 

    read (unit, *, IOSTAT=status, IOMSG=message) natoms
    if (status /= 0) exit catch 

    if (natoms < 0) then
      message = "Invalid number of atoms"
      status = 1
      exit catch
    end if

    read (unit, "(a)", IOSTAT=status, IOMSG=message) comment
    if (status /= 0) exit catch 

  end block catch

  if (present(stat)) then
    stat = status
  else if (status /= 0) then
    error stop message
  end if
  if (present(errmsg)) errmsg = trim(message)

end subroutine xyz_read_header

! Read data in XYZ format.
subroutine xyz_read_data (unit, natoms, name, coor, stat, errmsg)
  implicit none
  integer, intent(in) :: unit ! Fortran file handle.
  integer, intent(in) :: natoms ! Number of atoms
  character(*), intent(out) :: name(natoms) ! Atom names usually element symbol. 
  real, intent(out) :: coor(DIM,natoms) ! Atom cartesian coordinates in [A] units.
  integer, intent(out), optional :: stat ! Error status code. Returns zero if no error.
  character(*), intent(out), optional :: errmsg ! Error message.
  character(ERRLEN) :: message
  character(128) :: buffer
  character(16) :: action
  logical :: opened
  integer :: i, ncol, status

  catch: block

    ! Check if file unit is opened for reading
    inquire (UNIT=unit, OPENED=opened, ACTION=action, IOSTAT=status, IOMSG=message)
    if (.not. opened .or. index(action, "READ") == 0) then
      message = "File unit not opened for reading"
      status = 1
    end if
    if (status /= 0) exit catch 
  
    ! Read 1st line from a file and count number of columns on the line. First column
    ! * that contains atom names is optional in XYZ file. After that rewind the line.
    read (unit, "(a)", IOSTAT=status, IOMSG=message) buffer
    if (status /= 0) exit catch
    backspace (unit, IOSTAT=status, IOMSG=message)
    if (status /= 0) exit catch

    ! Count number of columns and check if everithing is OK. We expect 3 or more columns.
    ! "Unexpected number of columns: got '", ncol, "' expected '", DIM, "+'"
    ncol = cnttok(buffer, del=" ")
    if (ncol < DIM) then
      message = "Unexpected number of columns" 
      status = 1
      exit catch
    end if

    ! Read file data according to number of columns in the file.
    if (ncol > 3) then
      do i = 1, natoms
        read (unit, *, IOSTAT=status, IOMSG=message) name(i), coor(:,i)
        if (status /= 0) exit catch
      end do
    else
      do i = 1, natoms
        read (unit, *, IOSTAT=status, IOMSG=message) coor(:,i)
        if (status /= 0) exit catch
      end do
      name(:) = ""
    end if

  end block catch

  if (present(stat)) then
    stat = status
  else if (status /= 0) then
    error stop message
  end if
  if (present(errmsg)) errmsg = trim(message)

end subroutine xyz_read_data

! Read data in XYZ format and only store coordinates.
subroutine xyz_read_coor (unit, natoms, coor, stat, errmsg)
  implicit none
  integer, intent(in) :: unit ! Fortran file handle.
  integer, intent(in) :: natoms ! Number of atoms
  real, intent(out) :: coor(DIM,natoms) ! Atom cartesian coordinates in [A] units.
  integer, intent(out), optional :: stat ! Error status code. Returns zero if no error.
  character(*), intent(out), optional :: errmsg ! Error message.
  character(ERRLEN) :: message
  character(128) :: buffer
  character(16) :: action, dummy
  logical :: opened
  integer :: i, ncol, status

  catch: block

    ! Check if file unit is opened for reading
    inquire (UNIT=unit, OPENED=opened, ACTION=action, IOSTAT=status, IOMSG=message)
    if (.not. opened .or. index(action, "READ") == 0) then
      message = "File unit not opened for reading"
      status = 1
    end if
    if (status /= 0) exit catch 

    ! Read 1st line from a file and count number of columns on the line. First column
    ! * that contains atom names is optional in XYZ file. After that rewind the line.
    read (unit, "(a)", IOSTAT=status, IOMSG=message) buffer
    if (status /= 0) exit catch
    backspace (unit, IOSTAT=status, IOMSG=message)
    if (status /= 0) exit catch

    ! Count number of columns and check if everithing is OK. We expect 3 or more columns.
    ! "Unexpected number of columns: got '", ncol, "' expected '", DIM, "+'"
    ncol = cnttok(buffer, del=" ")
    if (ncol < DIM) then
      message = "Unexpected number of columns" 
      status = 1
      exit catch
    end if

    ! Read file data according to number of columns in the file.
    if (ncol > 3) then
      do i = 1, natoms
        read (unit, *, IOSTAT=status, IOMSG=message) dummy, coor(:,i)
        if (status /= 0) exit catch
      end do
    else
      do i = 1, natoms
        read (unit, *, IOSTAT=status, IOMSG=message) coor(:,i)
        if (status /= 0) exit catch
      end do
    end if

  end block catch

  if (present(stat)) then
    stat = status
  else if (status /= 0) then
    error stop message
  end if
  if (present(errmsg)) errmsg = trim(message)

end subroutine xyz_read_coor

! Read through XYZ file but not storing any data.
subroutine xyz_skip_data (unit, natoms, stat, errmsg)
  implicit none
  integer, intent(in) :: unit ! Fortran file handle.
  integer, intent(in) :: natoms ! Number of atoms
  integer, intent(out), optional :: stat ! Error status code. Returns zero if no error.
  character(*), intent(out), optional :: errmsg ! Error message.
  character(ERRLEN) :: message
  character(16) :: action
  logical :: opened
  integer :: i, status

  catch: block

    ! Check if file unit is opened for reading
    inquire (UNIT=unit, OPENED=opened, ACTION=action, IOSTAT=status, IOMSG=message)
    if (.not. opened .or. index(action, "READ") == 0) then
      message = "File unit not opened for reading"
      status = 1
    end if
    if (status /= 0) exit catch 

    if (natoms < 0) then
      message = "Invalid number of atoms"
      status = 1
      exit catch
    end if

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

end subroutine xyz_skip_data

! Write data in XYZ format.
subroutine xyz_write (unit, comment, natoms, name, coor, stat, errmsg)
  implicit none
  integer, intent(in) :: unit ! Fortran file handle.
  character(*), intent(in) :: comment ! Comment line
  integer, intent(in) :: natoms ! Number of atoms
  character(*), intent(in) :: name(natoms) ! Atom names usually element symbol. 
  real, intent(in) :: coor(DIM,natoms) ! Atom cartesian coordinates in [A] units.
  integer, intent(out), optional :: stat ! Error status code. Returns zero if no error.
  character(*), intent(out), optional :: errmsg ! Error message.
  character(ERRLEN) :: message
  character(9) :: action
  logical :: opened
  integer :: i, status

  catch: block

    ! Check if file unit is opened for writing.
    inquire (UNIT=unit, OPENED=opened, ACTION=action, IOSTAT=status, IOMSG=message)
    if (.not. opened .or. index(action, "WRITE") == 0) then
      message = "File unit not opened for writing"
      status = 1
    end if
    if (status /= 0) exit catch 

    write (unit, "(i0)", IOSTAT=status, IOMSG=message) natoms
    if (status /= 0) exit catch

    write (unit, "(a)", IOSTAT=status, IOMSG=message) trim(comment)
    if (status /= 0) exit catch

    do i = 1, natoms
      write (unit, "(a5,3(5x,f10.5))", IOSTAT=status, IOMSG=message) adjustl(name(i)), coor(:,i)
      if (status /= 0) exit catch
    end do

  end block catch

  if (present(stat)) then
    stat = status
  else if (status /= 0) then
    error stop message
  end if
  if (present(errmsg)) errmsg = trim(message)

end subroutine xyz_write

end module atomlib_xyzio
